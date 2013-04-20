#-------------------------------------------------------------------------------
# Script producing the figures/tables for the student paper competition entry.
#
# Adam Loy
# December 2012
#-------------------------------------------------------------------------------

### Preliminaries
library(lme4)
library(arm)       # for nicer print outs from lmer
library(ggplot2)
library(nortest)   # for tests of normality
library(boot)      # for simulation envelopes

### Reading in the radon data provided by Gelman and Hill (2007)
srrs2 <- read.table ("/Users/adam/Documents/Thesis/Gelman and Hill 2007/ARM_Data/radon/srrs2.dat", header=T, sep=",")

## Restricting attention to Minnesota
mn <- srrs2$state=="MN"
radon <- srrs2$activity[mn]
log.radon <- log (ifelse (radon==0, .1, radon))
floor <- srrs2$floor[mn]       # 0 for basement, 1 for first floor
n <- length(radon)
y <- log.radon
basement <- floor

## Getting county index variable
county.name <- as.vector(srrs2$county[mn])
uniq <- unique(county.name)
J <- length(uniq)
county <- rep (NA, J)
for (i in 1:J){
  county[county.name==uniq[i]] <- i
}

## Reading the county-level data
srrs2.fips <- srrs2$stfips*1000 + srrs2$cntyfips
cty <- read.table ("/Users/adam/Documents/Thesis/Gelman and Hill 2007/Book_Codes/Ch.12/cty.dat", header=T, sep=",")
usa.fips <- 1000*cty[,"stfips"] + cty[,"ctfips"]
usa.rows <- match (unique(srrs2.fips[mn]), usa.fips)
uranium <- cty[usa.rows,"Uppm"]
u <- log (uranium)
uranium <- u[county]


### Fitting a model with random intercept for county and 
### random slope for basement.
fm <- lmer(log.radon ~ basement + uranium + (basement | county), REML = FALSE)
lev1 <- resid(fm)
rint <- ranef(fm)[[1]][,1]
rslope <- ranef(fm)[[1]][,2]

bmat <- matrix(t(ranef(fm)[[1]]), ncol=1)

#-------------------------------------------------------------------------------
# Normal quantile plots for the sample residuals
#-------------------------------------------------------------------------------

sim_env <- function(x, conf = .95){
  n <- length(x)
  P <- ppoints(x)
  z <- qnorm(P)
  a <- as.numeric(HLMdiag:::qqlineInfo(x)[1])
  b <- as.numeric(HLMdiag:::qqlineInfo(x)[2])
  zz <- qnorm(1 - (1 - conf)/2)
  SE <- (b/dnorm(z)) * sqrt(P * (1 - P)/n)
  fit.value <- a + b * z
  upper <- fit.value + zz * SE
  lower <- fit.value - zz * SE
  return(data.frame(lower, upper))
}


naive1.env <- sim_env(lev1)
naive1.qq <- qqnorm(lev1, plot.it = FALSE)
naive1.qq <- lapply(naive1.qq, sort)

ggplot(data.frame(lev1), aes(sample = lev1)) + 
	geom_abline(intercept = HLMdiag:::qqlineInfo(lev1)[1], 
              slope = HLMdiag:::qqlineInfo(lev1)[2], colour = I("grey60")) + 
	stat_qq() + 
	geom_ribbon(aes(x = naive1.qq$x, ymin = naive1.env$lower, ymax = naive1.env$upper), alpha = .2) + 
	xlab("Normal Quantiles") + ylab("Sample Quantiles") + 
	ylim(-4,4) + 
	theme_bw()
ggsave(filename = "/Users/adam/Documents/Thesis/Dissertation/paper competition 2012/figures/raw-lev1-qq.pdf", width = 3, height = 3)

naive.rint.env <- sim_env(rint)
naive.rint.qq <- qqnorm(rint, plot.it = FALSE)
naive.rint.qq <- lapply(naive.rint.qq, sort)

ggplot(data.frame(rint), aes(sample = rint)) + 
	geom_abline(intercept = HLMdiag:::qqlineInfo(rint)[1], 
              slope = HLMdiag:::qqlineInfo(rint)[2], colour = I("grey60")) + 
	stat_qq() + 
	geom_ribbon(aes(x = naive.rint.qq$x, ymin = naive.rint.env$lower, ymax = naive.rint.env$upper), alpha = .2) + 
	xlab("Normal Quantiles") + ylab("Sample Quantiles") + 
	ylim(-.2, .2) + 
	theme_bw()
ggsave(filename = "/Users/adam/Documents/Thesis/Dissertation/paper competition 2012/figures/raw-intercept-qq.pdf", width = 3, height = 3)

ggplot(data.frame(rslope), aes(sample = rslope)) + 
	geom_abline(intercept = HLMdiag:::qqlineInfo(rslope)[1], 
              slope = HLMdiag:::qqlineInfo(rslope)[2], colour = I("grey60")) + 
	stat_qq() + 
	geom_ribbon(aes(x = naive.rint.qq$x, ymin = naive.rint.env$lower, ymax = naive.rint.env$upper), alpha = .2) +
	xlab("Normal Quantiles") + ylab("Sample Quantiles") + 
	ylim(-.6, .6) + 
	theme_bw()
ggsave(filename = "/Users/adam/Documents/Thesis/Dissertation/paper competition 2012/figures/raw-slope-qq.pdf", width = 3, height = 3)

#-------------------------------------------------------------------------------
# Using lineups
#-------------------------------------------------------------------------------

library(nullabor)
library(HLMdiag)
library(plyr)

set.seed(987654321)
y.sim <- simulate(fm, nsim = 19)				# simulate responses
mod.sim <- apply(y.sim, 2, refit, object = fm)	# refit model 

## Level-1
r1 <- data.frame(resid = resid(fm))
sim.lev1.resid <- ldply(mod.sim, function(x) data.frame(resid = resid(x)))
sim.lev1.resid$.n <- as.numeric(sapply(strsplit(sim.lev1.resid$.id, "\\_"), "[", 2))

qplot(sample = resid, data = r1, stat = "qq") %+%
	lineup(true = r1, sample = sim.lev1.resid) + 
	facet_wrap(~ .sample, ncol = 4) + 
	xlab("Normal Quantiles") + ylab("Order Statistics")
	
## Level-2
r2 <- data.frame(X.Intercept. = rint, basement = rslope)
r2 <- transform(r2, band = sim_env(basement), x = sort(qqnorm(basement, plot.it=FALSE)$x))


sim.lev2.resid <- ldply(mod.sim, function(x) data.frame(ranef(x)[[1]]))
sim.lev2.resid$.n <- as.numeric(sapply(strsplit(sim.lev2.resid$.id, "\\_"), "[", 2))
sim.lev2.resid <- ddply(sim.lev2.resid, .(.id), transform, band = sim_env(basement), x = sort(qqnorm(basement, plot.it=FALSE)$x))

# Intercept
qplot(sample = X.Intercept., data = r2, stat = "qq") %+%
	lineup(true = r2, sample = sim.lev2.resid) + 
	facet_wrap(~ .sample, ncol = 4) + 
	xlab("Normal Quantiles") + ylab("Sample Quantiles")

# Slope	
library(grid)
set.seed(123456)
qplot(sample = basement, data = r2, stat = "qq") %+%
	lineup(true = r2, sample = sim.lev2.resid) + 
	facet_wrap(~ .sample, ncol = 5) + 
	geom_ribbon(aes(x = x, ymin = band.lower, ymax = band.upper), alpha = .25) + 
	xlab("Normal Quantiles") + ylab("Sample Quantiles") +  
	theme_bw() + 
	theme(panel.margin = unit(0, "lines"))
ggsave("test.jpeg")

#-------------------------------------------------------------------------------
# Least confounded residuals
#-------------------------------------------------------------------------------

BlockZ <- function(object) {
  Z <- getME(object, "Z")
  
  grp.size <- table(object@flist)
  ngrps <- length(grp.size)
  nranef <- dim(ranef(object)[[1]])[2]
  
  base.ord <- seq(from = 1, by = ngrps, length.out = nranef)
  ord <- base.ord + rep(0:(ngrps - 1), each = nranef)
  
  perm.mat <- t(as(ord, "pMatrix"))
  
  return(Z %*% perm.mat)
}

mcrotate <- function(A, B) {
  r <- rankMatrix(B)
  
  B.svd <- svd(B)
  Cr.diag <- B.svd$d[1:r]
  Tr <- B.svd$u[, 1:r]
  
  A.star <- Diagonal( x = 1 / sqrt(Cr.diag) ) %*% t(Tr) %*% 
    A %*% 
    Tr %*% Diagonal( x = 1 / sqrt(Cr.diag) )
  
  A.star.svd <- svd( A.star )
  
  W <- Tr %*% Diagonal( x = 1 / sqrt( Cr.diag ) ) %*% A.star.svd$u
  
  return(W)
}

X <- getME(fm, "X")
Z <- BlockZ(fm) 

n <- nrow(X)
p <- ncol(X)
ngrps <- unname( summary(fm)@ngrps )
	           
vc <- VarCorr(fm)
D <- kronecker( Diagonal(ngrps), as.matrix(vc[[1]]))
R <- Diagonal( n = n, x = attr(vc, "sc")^2 )

V <- Z %*% D %*% t(Z) + R
V.chol <- chol( V )
Vinv  <- chol2inv( V.chol )

xvix <- t(X) %*% Vinv %*% X
xvi <- t(X) %*% Vinv
rhs <- crossprod(xvi, solve(xvix, xvi))
P <- Vinv - rhs

### Level-1
A1 <- R %*% P %*% Z %*% D %*% t(Z) %*% P %*% R
B1 <- R %*% P %*% R

M1 <- mcrotate(A = A1, B = B1)
lcr1 <- as.matrix( t(M1) %*% resid(fm) )

qplot(sample = lcr1, stat = "qq", ylab = "Order Statistics", 
      xlab = "Normal Quantiles") + theme_bw() + 
  geom_abline(intercept = 0, slope = 1)



### Random intercept
L.int <- kronecker(Diagonal(ngrps), c(1, 0))

A.int <- t(L.int) %*% D %*% t(Z) %*% P %*% R %*% P %*% Z %*% D %*% L.int
B.int <- t(L.int) %*% D %*% t(Z) %*% P %*% Z %*% D %*% L.int

M.int <- mcrotate(A = A.int, B = B.int)
lcr2.int <- as.matrix( t(M.int) %*% rint )


lcrint.env <- sim_env(lcr2.int)
lcrint.qq <- qqnorm(lcr2.int, plot.it = FALSE)
lcrint.qq <- lapply(lcrint.qq, sort)

qplot(sample = lcr2.int, stat = "qq", ylab = "Sample Quantiles", 
      xlab = "Normal Quantiles") + theme_bw() + 
  geom_abline(intercept = HLMdiag:::qqlineInfo(lcr2.int)[1], 
              slope = HLMdiag:::qqlineInfo(lcr2.int)[2])

rsv <- ranef(fm, postVar=TRUE)
rslope.var <- apply(attr(rsv[[1]], "postVar"), 3, function(x) x[2,2])
qplot(x = reorder(county, var, identity), y = rslope, data = data.frame(county = str_trim(unique(county.name)), rslope = rslope, var = rslope.var)) + 
geom_linerange(aes(ymin= rslope-2*var, ymax=rslope+2*var)) + 
coord_flip() + 
theme(axis.text.y = element_text(size=6, colour="black"))


ggplot(data.frame(lcr2.int), aes(sample = lcr2.int)) + 
	geom_abline(intercept = HLMdiag:::qqlineInfo(lcr2.int)[1], 
              slope = HLMdiag:::qqlineInfo(lcr2.int)[2], colour = I("grey60")) + 
	stat_qq() + 
	geom_ribbon(aes(x = lcrint.qq$x, ymin = lcrint.env$lower, ymax = lcrint.env$upper), alpha = .2, fullrange=TRUE) +
	xlab("Normal Quantiles") + ylab("Sample Quantiles") +  
	ylim(c(-3.6, 3.6)) + 
	theme_bw()
ggsave(filename = "/Users/adam/Documents/Thesis/Dissertation/paper competition 2012/figures/lcr-intercept-qq.pdf", width = 3, height = 3)

image(t(M.int)[c(19, 29, 59),])
ggfluctuation(as.table(as.matrix(t(M.int)[c(19, 29, 59),])), type = "colour")
out.df <- data.frame( row = c(19, 29, 59), as.matrix(zapsmall(t(M.int)[c(19, 29, 59),])))
out.melt <- melt(out.df, id.vars="row")
out.melt$col <- as.numeric(out.melt$variable)
out.melt$max <- pmax(out.melt$value, 0)
out.melt$min <- pmin(out.melt$value, 0)

qplot(x = col, y = value, data=out.melt, facets = ~ row, geom = c("point", "line")) + coord_flip() + xlab("Column") + ylab("Weight")

ggplot(out.melt, aes(xmin = col-.5, xmax = col + .5, ymin = min, ymax = max)) + geom_rect() + facet_grid(facets = row~.) + geom_hline(yintercept=0) + ylab("Weights") + xlab("Column") + theme_bw()

ggfluctuation(as.table(as.matrix(t(M.int)[c(19, 29, 59),])))

### Random slope
L.slope <- kronecker(Diagonal(ngrps), c(0, 1))

# A.slope <- as.matrix(t(L.slope) %*% D %*% t(Z) %*% P %*% R %*% P %*% Z %*% D %*% L.slope)
A.slope <-  as.matrix( attr(vc, "sc")^2 * crossprod(P %*% Z %*% D %*% L.slope) )
B.slope <- t(L.slope) %*% D %*% t(Z) %*% P %*% Z %*% D %*% L.slope

M.slope <- mcrotate(A = A.slope, B = B.slope)
lcr2.slope <- as.matrix( t(M.slope) %*% rslope )

# QQ plot of rotate residuals
qplot(sample = lcr2.slope, stat = "qq", ylab = "Order Statistics", 
      xlab = "Normal Quantiles") + theme_bw() + 
  geom_abline(intercept = HLMdiag:::qqlineInfo(lcr2.slope)[1], 
              slope = HLMdiag:::qqlineInfo(lcr2.slope)[2])
 
lcrslope.env <- sim_env(lcr2.slope)
lcrslope.qq <- qqnorm(lcr2.slope, plot.it = FALSE)
lcrslope.qq <- lapply(lcrslope.qq, sort)

ggplot(data.frame(lcr2.slope), aes(sample = lcr2.slope)) + 
	geom_abline(intercept = 0, 
              slope = 1, colour = I("grey60")) + 
	stat_qq() + 
	geom_ribbon(aes(x = lcrslope.qq$x, ymin = lcrslope.env$lower, ymax = lcrslope.env$upper), alpha = .2) + 
	xlab("Normal Quantiles") + ylab("Sample Quantiles") + 
	ylim(-3.6, 3.6) + 
	theme_bw()
ggsave(filename = "/Users/adam/Documents/Thesis/Dissertation/paper competition 2012/figures/lcr-slope-qq.pdf", width = 3, height = 3)

#-------------------------------------------------------------------------------
# Exploring weight matrix
#-------------------------------------------------------------------------------

qplot(sample = rslope, stat = "qq")

wts.means <- apply(as.matrix(t(M.slope)), 2, mean)
wts.var <- apply(as.matrix(t(M.slope)), 2, var)

image(t(M.slope))

wts.df <- data.frame(wts = wts.means, raw = rslope, 
	quants = qqnorm(rslope, plot.it=F)$x)
wts.df$color <- "Mid"
wts.df$color[c(50, 32, 73, 51, 11, 6, 85, 60, 52, 82, 16, 12, 68, 42)] <- "Tail"
wts.df$color[c(73, 51, 11, 6, 85, 60, 52, 82, 16, 12)] <- "Tail" 
 
wts.df$color <- factor(wts.df$color)
wts.df$color <- factor(wts.df$color, levels(wts.df$color)[2:1])

### QQplot and histogram

qplot(x = wts, data = wts.df, geom = "histogram", fill = factor(color),
  binwidth=0.1, colour = I("grey90"), size=I(0.25)) + 
  xlab("weights") +
  scale_fill_manual(values=c("#E41A1C", "#377EB8")) + 
  theme_bw() + 
  theme(legend.position = "none") 
ggsave(filename = "hist-wts-tail.pdf", width=3, height=3)

ggplot(wts.df, aes(x = quants, y = raw, sample = raw)) + 
	geom_abline(intercept = HLMdiag:::qqlineInfo(rslope)[1], 
              slope = HLMdiag:::qqlineInfo(rslope)[2], colour = I("grey60")) + 
	geom_point(aes(colour = factor(color))) + 
	geom_point(data = subset(wts.df, color == "Tail"), colour = "#E41A1C") +
	xlab("Normal Quantiles") + ylab("Sample Quantiles") + 
	ylim(-.6, .6) + 
	theme_bw() + 
	scale_colour_manual(values=c("#E41A1C", "#377EB8")) +
	theme(legend.position = "none")
ggsave(filename = "qq-wts-tail.pdf", width=3, height=3)


wts.df$color <- "Large"
wts.df$color[c(84,83,81,80,78,67,64,63,56,54,47,39,34,28,27,24,18,17,15,14,9,7,4,2,1,79,74,70,66,62,61,59,58,57,49,45,44,37,36,35,26,25,21,19,10)] <- "Small"

qplot(x = wts, data = wts.df, geom = "histogram", fill = abs(wts) > .1,
  binwidth=0.1, colour = I("black")) + 
  xlab("weights") +
  scale_fill_manual(values=c("#E41A1C", "#377EB8")) + 
  theme_bw() + 
  theme(legend.position = "none") 
ggsave(filename = "hist-wts-mid.pdf", width=3, height=3)

ggplot(wts.df, aes(x = quants, y = raw, sample = raw)) + 
	geom_abline(intercept = HLMdiag:::qqlineInfo(rslope)[1], 
              slope = HLMdiag:::qqlineInfo(rslope)[2], colour = I("grey60")) + 
	geom_point(aes(colour = abs(wts) < .1)) + 
	geom_point(data = subset(wts.df, abs(wts) < .1), colour = "#E41A1C") +
	xlab("Normal Quantiles") + ylab("Sample Quantiles") + 
	ylim(-.6, .6) + 
	theme_bw() + 
	scale_colour_manual(values=c("#377EB8", "#E41A1C")) +
	theme(legend.position = "none")
ggsave(filename = "qq-wts-mid.pdf", width=3, height=3)
   

iplot(x = wts.df$quants, y = wts.df$raw, xlab = "Theoretical Quantiles", 
		ylab = "Sample Quantiles")
ihist(wts.df$wts)
iplot(x = wts.df$quants, y = wts.df$wts)
iplot(x = as.numeric(table(county)), y = wts.df$wts)


wts.means <- apply(as.matrix(t(M.int)), 2, mean)
wts.var <- apply(as.matrix(t(M.int)), 2, var)

wts.df <- data.frame(wts = wts.means, raw = rint, 
	quants = qqnorm(rint, plot.it=F)$x, low = wts.means < -.5, 
	high = wts.means > .5)


qplot(x = wts, y = raw, data = wts.df)
qplot(sample = raw, stat = "qq", data = wts.df) +  
	geom_point(aes(x = qqnorm(rslope, plot.it=FALSE)$x, y = qqnorm(rslope, plot.it=FALSE)$y), colour = wts.df$high)


tmslope <- as.matrix(zapsmall(t(M.slope)))
colnames(tmslope) <- unique(str_trim(county.name))

### Clustered heat map for M.slope
heatmap(tmslope)
hc.rows <- hclust(dist(tmslope))
plot(hc.rows)

hc.cols <- hclust(dist(as.matrix(t(tmslope))))
plot(hc.cols)

# heatmap for first cluster
heatmap(tmslope[cutree(hc.rows,k=4)==4,], Colv=as.dendrogram(hc.cols), scale='none')

# heatmap for second cluster
palette <- colorRampPalette(c('#f0f3ff','#0033BB'))(256)
heatmap(tmslope[cutree(hc.rows,k=4)==4,], Colv=as.dendrogram(hc.cols), scale='none', col=palette)

ggfluctuation <- function (table, type = "size", floor = 0, ceiling = max(table$freq, 
    na.rm = TRUE)) 
{
    .Deprecated()
    if (is.table(table)) 
        table <- as.data.frame(t(table))
    oldnames <- names(table)
    names(table) <- c("x", "y", "result")
    table <- transform(table, x = as.factor(x), y = as.factor(y), 
        freq = result)
    if (type == "size") {
        table <- transform(table, freq = sqrt(pmin(freq, ceiling)/ceiling), 
            border = ifelse(is.na(freq), "grey90", ifelse(freq > 
                ceiling, "grey30", "grey50")))
        table[is.na(table$freq), "freq"] <- 1
        table <- subset(table, freq * ceiling >= floor)
    }
    if (type == "size") {
        nx <- length(levels(table$x))
        ny <- length(levels(table$y))
        p <- ggplot(table, aes_string(x = "x", y = "y", height = "freq", 
            width = "freq", fill = "border")) + geom_tile(colour = "white") + 
            scale_fill_identity() + theme(aspect.ratio = ny/nx)
    }
    else {
        p <- ggplot(table, aes_string(x = "x", y = "y", fill = "freq")) + 
            geom_tile(colour = "white") + scale_fill_gradient2(expression(sign(w[ij])*sqrt(abs(w[ij])))
, low = "red", mid="white", high = "blue")
    }
    p$xlabel <- oldnames[1]
    p$ylabel <- oldnames[2]
    p
}

library(stringr)
tmslope <- as.matrix(zapsmall(t(M.slope)))

colnames(tmslope) <- str_trim(unique(county.name))
rownames(tmslope) <- as.character(1:84)
heatmap.info <- heatmap(tmslope)

# reordered.tmslope <- tmslope[heatmap.info$rowInd, heatmap.info$colInd]
# reordered.tmslope2 <- reordered.tmslope[,c(4:84, 3:1, 85)]
reordered.tmslope <- tmslope[, order(rslope.var)]
reordered.tmslope2.tform <- sign(reordered.tmslope) * sqrt(abs(reordered.tmslope))
# reordered.tmslope2.tform <- (reordered.tmslope2)^(1/3)


ggfluctuation(as.table(reordered.tmslope2.tform), type = "colour") + 
  theme(axis.text.x = theme_text(angle = 90, hjust=1, vjust=.5, size=6, colour="black"),
    axis.text.y = theme_text(size=6, colour="black")) + 
  xlab("County") + ylab("Index of Rotated Residuals") + 
  coord_equal() + 
  theme(panel.margin = unit(0, "lines"))

qplot(x = reorder(county, rslope, identity), y = rslope, data = data.frame(county = str_trim(unique(county.name)), rslope = rslope)) + coord_flip() + 
theme(axis.text.y = element_text(size=6, colour="black"))


### clustered heat map of Lambda_r^-1/2 t(T_r)
B.slope.svd <- svd(B.slope)

sqrt.cov <- as.matrix(zapsmall(Diagonal( x = 1 / sqrt(B.slope.svd$d[1:84]) ) %*% t(B.slope.svd$u[, 1:84])))
colnames(sqrt.cov) <- str_trim(unique(county.name))
rownames(sqrt.cov) <- as.character(1:84)

tform.sqrt.cov <- sign(sqrt.cov) * sqrt(abs(sqrt.cov))

ggfluctuation(as.table(tform.sqrt.cov[, order(rslope.var)]), type = "colour") +
  theme(axis.text.x = theme_text(angle = 90, hjust=1, vjust=.5, size=6, colour="black"), axis.text.y = theme_text(size=6, colour="black")) +
  xlab("County") + ylab("Index of Rel. Precision Factor") + 
  coord_equal() + 
  theme(panel.margin = unit(0, "lines"))
  
# confounding
confound.df <- data.frame(confound = diag(A.slope) / diag(B.slope), 
  county = unique(str_trim(county.name)))

  
### heat map of U
numerator.svd <- svd( sqrt.cov %*% A.slope %*% t(sqrt.cov) )
Ut <- t(numerator.svd$u)

ggfluctuation(as.table(as.matrix(Ut)), type = "colour")

#-------------------------------------------------------------------------------
# Simulated LC errors 
#-------------------------------------------------------------------------------

set.seed(987654321)

sim <- simulate(fm)
fm.sim <- refit(fm, sim)

X.sim <- getME(fm.sim, "X")
Z.sim <- BlockZ(fm.sim) 

n.sim <- nrow(X.sim)
p.sim <- ncol(X.sim)
ngrps.sim <- unname( summary(fm.sim)@ngrps )
	           
vc.sim <- VarCorr(fm.sim)
D.sim <- kronecker( Diagonal(ngrps.sim), as.matrix(vc.sim[[1]]))
R.sim <- Diagonal( n = n.sim, x = attr(vc.sim, "sc")^2 )

V.sim <- Z.sim %*% D.sim %*% t(Z.sim) + R.sim
V.chol.sim <- chol( V.sim )
Vinv.sim  <- chol2inv( V.chol.sim )

xvix.sim <- t(X.sim) %*% Vinv.sim %*% X.sim
xvi.sim <- t(X.sim) %*% Vinv.sim
rhs.sim <- crossprod(xvi.sim, solve(xvix.sim, xvi.sim))
P.sim <- Vinv.sim - rhs.sim

### Level-1
A1.sim <- R.sim %*% P.sim %*% Z.sim %*% D.sim %*% t(Z.sim) %*% P.sim %*% R.sim
B1.sim <- R.sim %*% P.sim %*% R.sim

M1.sim <- mcrotate(A = A1.sim, B = B1.sim)
lcr1.sim <- as.matrix( t(M1.sim) %*% resid(fm.sim) )

### Random intercept
L.int.sim <- kronecker(Diagonal(ngrps.sim), c(1, 0))

A.int.sim <- t(L.int.sim) %*% D.sim %*% t(Z.sim) %*% P.sim %*% R.sim %*% P.sim %*% Z.sim %*% D.sim %*% L.int.sim
B.int.sim <- t(L.int.sim) %*% D.sim %*% t(Z.sim) %*% P.sim %*% Z.sim %*% D.sim %*% L.int.sim

M.int.sim <- mcrotate(A = A.int.sim, B = B.int.sim)
lcr2.int.sim <- as.matrix( t(M.int.sim) %*% ranef(fm.sim)[[1]][,1] )

image(t(M.int.sim))

j <- apply(t(M.int.sim), 1, function(x) sum(abs(x)))



### Random slope
L.slope.sim <- kronecker(Diagonal(ngrps.sim), c(0, 1))

# A.slope <- as.matrix(t(L.slope) %*% D %*% t(Z) %*% P %*% R %*% P %*% Z %*% D %*% L.slope)
A.slope.sim <-  as.matrix( attr(vc.sim, "sc")^2 * crossprod(P.sim %*% Z.sim %*% D.sim %*% L.slope.sim) )
B.slope.sim <- t(L.slope.sim) %*% D.sim %*% t(Z.sim) %*% P.sim %*% Z.sim %*% D.sim %*% L.slope.sim

M.slope.sim <- mcrotate(A = A.slope.sim, B = B.slope.sim)
lcr2.slope.sim <- as.matrix( t(M.slope.sim) %*% ranef(fm.sim)[[1]][,2] )



#-------------------------------------------------------------------------------
# Recoverred errors 
#-------------------------------------------------------------------------------

### Level-1
A1 <- R %*% P %*% R
A1.svd <- svd(A1)

U <- A1.svd$u[,1:(n-p)]
C <- Diagonal(n = n - p, x = 1/ sqrt( A1.svd$d[1:(n-p)] ) )
Q <- U %*% C

lre1 <- as.matrix( t(Q) %*% resid(fm) )

B <- varimax(as.matrix(Q))
lre1.vmx <- t(B$loadings) %*% resid(fm)

### Level-2 -- Intercept
L.int <- kronecker(Diagonal(ngrps), c(1, 0))
A2i <-  t(L.int) %*% D %*% t(Z) %*% P %*% Z %*% D %*% L.int
A2i.svd <- svd(A2i)

U <- A2i.svd$u[,1:84]
C <- Diagonal(n = 84, x = 1/ sqrt( A2i.svd$d[1:84] ) )
Q <- U %*% C

lre2i <- as.matrix( t(Q) %*% rint )

B <- varimax(as.matrix(Q))
lre2i.vmx <- t(B$loadings) %*% rint

### Level-2 -- slope
L.slope <- kronecker(Diagonal(ngrps), c(0, 1))
A2s <-  t(L.slope) %*% D %*% t(Z) %*% P %*% Z %*% D %*% L.slope
A2s.svd <- svd(A2s)

U <- A2s.svd$u[,1:84]
C <- Diagonal(n = 84, x = 1/ sqrt( A2s.svd$d[1:84] ) )
Q <- U %*% C

lre2s <- as.matrix( t(Q) %*% rslope )

B <- varimax(as.matrix(Q))
lre2s.vmx <- t(B$loadings) %*% rslope



###########################
# maps 

srrs2 <- read.table ("data/srrs2.dat", header=T, sep=",")
## Restricting attention to Minnesota
mn <- subset(srrs2, state=="MN")

## clean names of counties: lower case and remove all white space at the end of names
mn$county <- tolower(as.character(mn$county))
mn$county <- gsub(" +$", "", as.character(mn$county))

library(ggplot2)
library(maps)
counties <- map_data("county")

# restrict to minnesota
counties <- subset(counties, region=="minnesota")
counties$region <- counties$subregion # no point in keeping the regions
qplot(data=counties, x=long, y=lat, group=group, geom="polygon")

# county level summary of radon data
library(plyr)
mn <- cbind(mn, uranium)
counties.stat <- ddply(mn, .(county), summarize, 
                     radon=mean(activity, na.rm=T),
                    sdradon=sd(activity, na.rm=T),
                    uranium=unique(uranium))

ggplot(counties.stat, aes(map_id = county)) + geom_map(aes(fill=uranium), map = counties) + 
  expand_limits(x = counties$long, y = counties$lat) + xlab("longitude") + ylab("latitude") + coord_map() + theme_bw()

# county level sample size
ggplot(counties.stat, aes(map_id = county)) + geom_map(aes(fill=radon), map = counties) + 
  expand_limits(x = counties$long, y = counties$lat) + 
  xlab("longitude") + ylab("latitude") + coord_map() + 
  theme(
  	axis.text.x = element_blank(), 
  	axis.text.y = element_blank(),
  	axis.title.x = element_blank(), 
  	axis.title.y = element_blank(),
  	axis.ticks = element_line(colour=rgb(0,0,0,alpha=0)),
  	panel.background =  element_blank(),
  	panel.grid.major =  element_blank(),
  	panel.grid.minor =  element_blank(),
  	legend.title = element_blank())
ggsave(file="map.png", width=2, height=3)

twocounties <- subset(mn, county %in% c("hennepin", "winona"))
twocounties$basement <- twocounties$floor==0
tc.stat <- ddply(twocounties, .(county, basement), summarise, activity=mean(activity), sd=sd(activity)/length(activity))

qplot(factor(basement), activity, data=twocounties, facets=~county) + 
  geom_errorbar(aes(ymin=activity-1.96*sd, ymax=activity+1.96*sd), data=tc.stat, width=0.5, size=0.75, colour="steelblue")+
  geom_point(data=tc.stat, size=3.5, colour="steelblue") + theme_bw()+theme(aspect.ratio=1) + ylab("radon") + xlab("basement")
ggsave(file="figures/radon-twocounties.pdf", width=6, height=3)
