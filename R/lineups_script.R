#-------------------------------------------------------------------------------
# Script fitting the model(s) to the radon data; and obtaining residuals
# for lineups. 
#
# Adam Loy
# April 2013
#-------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------
### Preliminaries
# -----------------------------------------------------------------------------------------
library(lme4)
library(plyr)
library(reshape2)
library(stringr)
library(nullabor)
library(grid)
library(ggplot2)

radon <- read.csv(file.choose()) ## find radon_for_sims.csv

### Function to reorganize Z like it should be in 511
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

### Function to calculate the marginal variance of random effects
lev2.marginal.var <- function(.model) {
  y <- .model@y
  X <- getME(.model, "X")
  Z <- BlockZ(.model)
  n <- nrow(X)
  ngrps <- unname(sapply(.model@flist, function(x) length(levels(x))))
  
  # Constructing V = Cov(Y)
  sig0 <- attr(VarCorr(.model), "sc") # sigma(.model)
  
  ZDZt <- sig0^2 * crossprod( .model@A )
  R    <- Diagonal( n = n, x = sig0^2 )
  D    <- kronecker( Diagonal(ngrps), bdiag(VarCorr(.model)) )
  V    <- Diagonal(n) + ZDZt
  
  # Inverting V
  V.chol <- chol( V )
  Vinv   <- chol2inv( V.chol )

  bse <- crossprod( chol(Vinv) %*% Z %*% D ) # Marginal COV. used by Lange and Ryan
  bse.diag <- diag(bse)

  semat <- matrix(sqrt(bse.diag), ncol = 2, byrow = TRUE)

  return(semat)
}

### Function for weighted empirical CDF
wecdf <- function(x, weights) {
    stopifnot(length(x) == length(weights))
    sw <- sum(weights)
    if (length(x) < 1) 
        stop("'x' must have 1 or more non-missing values")
    stopifnot(all(weights >= 0))
    ox <- order(x)
    x  <- x[ox]
    w  <- weights[ox]
    vals <- sort(unique(x))
    xmatch <- factor(match(x, vals), levels = seq_along(vals))
    wmatch <- tapply(w, xmatch, sum)
    wmatch[is.na(wmatch)] <- 0
    rval <- approxfun(vals, cumsum(wmatch) / sw, method = "constant", 
        yleft = 0, yright = 1, f = 0, ties = "ordered")
    class(rval) <- c("ecdf", "stepfun", class(rval))
    attr(rval, "call") <- sys.call()
    return(rval)
}   

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


weighted_sim_env <- function(quants, conf = .95){
  n <- length(quants)
  P <- pnorm(quants)
  z <- quants
  a <- 0
  b <- 1
  zz <- qnorm(1 - (1 - conf)/2)
  SE <- (b/dnorm(z)) * sqrt(P * (1 - P)/n)
  fit.value <- a + b * z
  upper <- fit.value + zz * SE
  lower <- fit.value - zz * SE
  return(data.frame(lower, upper))
}

# -----------------------------------------------------------------------------------------
### Fitting the models
# -----------------------------------------------------------------------------------------

# Random intercept, random slope model
fm <- lmer(log.radon ~ basement + uranium + (basement | county), data = radon)

# -----------------------------------------------------------------------------------------
### Calculating the original residuals
# -----------------------------------------------------------------------------------------

# Level-1 residuals (i.e., the error terms)
e <- resid(fm)

# Random effects (i.e., the level-2 residuals)
b <- ranef(fm)[[1]] # notice that this is actually a matrix

# -----------------------------------------------------------------------------------------
### Simulating null residuals via the parametric bootstrap
# -----------------------------------------------------------------------------------------
set.seed(987654321)
sim.y   <- simulate(fm, nsim = 19)                        ## A 919 x 19 matrix of responses
sim.mod <- apply(sim.y, 2, refit, object = fm)            ## a list of models

### We are really interested in a lineup for the random slopes for our paper
# Simulated random slopes
sim.b1 <- llply(sim.mod, function(x) ranef(x)[[1]][,2])   ## a list of random slopes
sim.b1 <- melt( do.call("rbind", sim.b1) )[,-2]           ## changing to a data frame
names(sim.b1) <- c("sample", "basement")                  ## setting colnames for faceting
sim.b1        <- arrange(sim.b1, sample)                  ## ordering by simulation
# sim.b1.df <- rbind(sim.b1,
# 			 cbind(sample = "true", slope = b[,2]))       ## adding the true r. slopes
# sim.b1.df$slope <- as.numeric(sim.b1.df$slope)            ## making residuals numeric


### Creating a the lineup
b1 <- transform(b, band = sim_env(basement), x = sort(qqnorm(basement, plot.it=FALSE)$x))

# some tweaks for nullabor

sim.b1$.n <- as.numeric( str_extract(sim.b1$sample, "\\d+") )
sim.b1 <- ddply(sim.b1, .(.n), transform, band = sim_env(basement), x = sort(qqnorm(basement, plot.it=FALSE)$x))

qplot(sample = basement, data = b1, stat = "qq") %+%
	lineup(true = b1, sample = sim.b1) + 
	facet_wrap(~ .sample, ncol = 5) + 
	geom_ribbon(aes(x = x, ymin = band.lower, ymax = band.upper), alpha = .25) + 
	xlab("Normal Quantiles") + ylab("Sample Quantiles") +  
	theme_bw() + 
	theme(panel.margin = unit(0, "lines"))

# That looks right in nullabor, so let's make a data frame for facetting
b1.facet <- b1[,-1]
b1.facet$sample <- "true"
sim.b1.facet <- sim.b1[,-3]

sim.b1.df.facetting <- rbind(b1.facet, sim.b1.facet)

# -----------------------------------------------------------------------------------------
### Calculations for a weighted Q-Q plot
# -----------------------------------------------------------------------------------------

# Marginal standard devation for original model
msd <- lev2.marginal.var(fm) ## ignore the notes

Fn.b1 <- wecdf(b[,2] / msd[,2], weights = msd[,1]^2)  ## create the weighted CDF function
p.b1  <- Fn.b1(b[,2] / msd[,2])                       ## calc. probs
qplot(x = qnorm(p.b1), y = b[,2] / msd[,2])

### Incorporating these into a simulation

# standardized random effects
b1.std <- b[,2] / msd[,2]

# Simulation for weighted Q-Q plots for random slopes
sim2.b1 <- llply(sim.mod, function(x) {
	se <- lev2.marginal.var(x)[,1]
	std.b1 <- ranef(x)[[1]][,1] / se
	Fn <- wecdf(std.b1, weights = se^2)
	p  <- Fn(std.b1)
	qs <- qnorm(p)
	RVAL <- data.frame(slope = std.b1, quants = qs)  ## list of data frames for plotting
	return(RVAL)
} )   

sim2.b1 <- do.call("rbind", sim2.b1)                             ## changing to a single data frame
sim2.b1$.n <- as.numeric(str_extract(rownames(sim2.b1), "\\d+")) ## adding col. for sim2. number


### Creating a the lineup
b1.std.df <- data.frame(slope = b1.std, quants = qnorm(p.b1))
b1.wtd    <- transform(b1.std.df, band = weighted_sim_env(quants))
b1.wtd    <- b1.wtd[complete.cases(b1.wtd),]

# some tweaks for nullabor

sim2.b1$.n <- as.numeric( str_extract(sim.b1$sample, "\\d+") )
sim2.b1 <- ddply(sim2.b1, .(.n), transform, band = weighted_sim_env(quants))
sim2.b1 <- sim2.b1[complete.cases(sim2.b1),]

qplot(x = quants, y = slope, data = b1.wtd, geom = "point") %+%
	lineup(true = b1.wtd, sample = sim2.b1) + 
	facet_wrap(~ .sample, ncol = 5) + 
	geom_ribbon(aes(x = quants, ymin = band.lower, ymax = band.upper), alpha = .25) + 
	xlab("Normal Quantiles") + ylab("Sample Quantiles") +  
	theme_bw() + 
	theme(panel.margin = unit(0, "lines"))

# That looks right in nullabor, so let's make a data frame for facetting
b1.wtd.facet <- b1.wtd
b1.wtd.facet$sample <- "true"
sim2.b1.facet <- sim2.b1
names(sim2.b1.facet)[3] <- "sample"

sim2.b1.df.facetting <- rbind(b1.wtd.facet, sim2.b1.facet)