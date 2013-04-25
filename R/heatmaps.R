#-------------------------------------------------------------------------------
# Script producing heatmaps for the rotated residuals paper.
#
# Adam Loy
# April 2013
#-------------------------------------------------------------------------------

### Preliminaries
library(lme4)
library(arm)       # for nicer print outs from lmer
library(ggplot2)
library(nortest)   # for tests of normality
library(RcppEigen)
library(inline)
library(grid)

setwd("~/Documents/Thesis/Dissertation/eresids-chapter/minconfounded_chapter/")

source("../simulations/functions/cpp_functions.R")

# Data for the continuous random slope example
radon <- read.csv("data/radon_for_sims.csv")

# Simulated models under normality situation
normsims <- readRDS("../simulations/continuous slope radon sims/sim_models/sige2_sigb1/norm_norm_REMLmodels.RDS")

# the fitted model
fm <- normsims[[1]]

rm(normsims)

#-------------------------------------------------------------------------------
# Calculating the rotated residuals
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

mcrotate <- function(A, B, s) {
  r <- rankMatrix(B)
  
  B.svd <- svd(B)
  Cr.diag <- B.svd$d[1:r]
  Tr <- B.svd$u[, 1:r]
  
  A.star <- Diagonal( x = 1 / sqrt(Cr.diag) ) %*% t(Tr) %*% 
    A %*% 
    Tr %*% Diagonal( x = 1 / sqrt(Cr.diag) )
  
  A.star.svd <- svd( A.star )
  
  index <- seq(r, length.out = s, by = -1)
  index <- sort(index)
  W <- Tr %*% Diagonal( x = 1 / sqrt( Cr.diag ) ) %*% A.star.svd$u[,index]
  
  return(W)
}

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

### Extract necessary components of the fitted model
y <- fm@y
X <- getME(fm, "X")
Z <- BlockZ(fm) 

n <- nrow(X)
p <- ncol(X)
ngrps <- unname( summary(fm)@ngrps )
	           
vc <- VarCorr(fm)
Di <- Diagonal(x = c(vc[[1]], vc[[2]])) / (unname(attr(vc, "sc")))^2
D  <- kronecker( Diagonal(ngrps), Di )

Aslot <- fm@A # ZDZ'
zdzt <- crossprod( fm@A )
V  <- Diagonal( n ) + zdzt
V.chol <- chol( V )
Vinv  <- chol2inv( V.chol ) 

XVXinv <- solve( t(X) %*% Vinv %*% X )
VinvX  <- Vinv %*% X
M      <- VinvX %*% XVXinv %*% t(VinvX)
P      <- cxxmatsub(as.matrix(Vinv), as.matrix(M))

### Rotating the random slope
L.int <- kronecker(Diagonal(ngrps), c(1, 0))
L.slope <- kronecker(Diagonal(ngrps), c(0, 1))

pzdl <- P %*% Z %*% D %*% L.slope
A.slope <-  crossprod( pzdl )
B.slope <- t(L.slope) %*% D %*% t(Z) %*% P %*% Z %*% D %*% L.slope

pzdl <- P %*% Z %*% D %*% L.int
A.int <-  crossprod( pzdl )
B.int <- t(L.int) %*% D %*% t(Z) %*% P %*% Z %*% D %*% L.int


fc.int   <- diag(ginv(as.matrix(B.int)) %*% A.int)
fc.slope <- diag(ginv(as.matrix(B.slope)) %*% A.slope)

s <- 70 # length of rotated residual vector
### For the random intercept
W.int <- as.matrix( mcrotate(A = A.int, B = B.int, s = s) )
W.vmx.int <- varimax( W.int, normalize = FALSE )$loadings

Wt.int <- as.matrix( zapsmall( t(W.int) ) )
Wt.vmx.int <- as.matrix( zapsmall( t(W.vmx.int) ) )

# colnames(Wt.int) <- colnames(Wt.vmx.int) <- unique( radon$county.name )
colnames(Wt.int) <- colnames(Wt.vmx.int) <- unique( radon$county.name )
rownames(Wt.int) <- rownames(Wt.vmx.int) <- seq_len(nrow(Wt.int))

# Pulling off variances of random ints
rint.var <- diag(B.int) # we could order by this, but fc makes more sense

reordered.Wt.int <- Wt.int[, order(fc.int, decreasing = TRUE)]
reordered.Wt.vmx.int <- Wt.vmx.int[, order(fc.int, decreasing = TRUE)]

reordered.Wt.int.tform <- sign(reordered.Wt.int) * sqrt(abs(reordered.Wt.int))
reordered.Wt.vmx.int.tform <- sign(reordered.Wt.vmx.int) * sqrt(abs(reordered.Wt.vmx.int))

ggfluctuation(as.table(reordered.Wt.int.tform), type = "colour") + 
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=.5, size=6, colour="black"),
    axis.text.y = element_text(size=6, colour="black")) + 
  xlab("County") + ylab("Index of Rotated Residuals") + 
  coord_equal() + 
  theme(panel.margin = unit(0, "lines")) + 
  ggtitle(paste("s = ", s))

ggsave(filename = paste("RandomIntercept_s", s, ".pdf", sep=""), width = 8,  height = 8, units = "in")

ggfluctuation(as.table(reordered.Wt.vmx.int.tform), type = "colour") + 
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=.5, size=6, colour="black"),
    axis.text.y = element_text(size=6, colour="black")) + 
  xlab("County") + ylab("Index of Rotated Residuals") + 
  coord_equal() + 
  theme(panel.margin = unit(0, "lines")) + 
  ggtitle(paste("s = ", s, "(varimax rotation)"))

ggsave(filename = paste("RandomIntercept_s", s, "_varimax.pdf", sep=""), width = 8,  height = 8, units = "in")



s <- 50 # length of rotated residual vector
### For the random slope
W.slope <- as.matrix( mcrotate(A = A.slope, B = B.slope, s = s) )
W.vmx.slope <- varimax( W.slope, normalize = FALSE )$loadings

Wt.slope <- as.matrix( zapsmall( t(W.slope) ) )
Wt.vmx.slope <- as.matrix( zapsmall( t(W.vmx.slope) ) )

# colnames(Wt.slope) <- colnames(Wt.vmx.slope) <- unique( radon$county.name )
colnames(Wt.slope) <- colnames(Wt.vmx.slope) <- unique( radon$county.name )
rownames(Wt.slope) <- rownames(Wt.vmx.slope) <- seq_len(nrow(Wt.int))

# Pulling off variances of random slopes
rslope.var <- diag(B.slope) # we could order by this, but fc makes more sense

reordered.Wt.slope <- Wt.slope[, order(fc.slope, decreasing = TRUE)]
reordered.Wt.vmx.slope <- Wt.vmx.slope[, order(fc.slope, decreasing = TRUE)]

reordered.Wt.slope.tform <- sign(reordered.Wt.slope) * sqrt(abs(reordered.Wt.slope))
reordered.Wt.vmx.slope.tform <- sign(reordered.Wt.vmx.slope) * sqrt(abs(reordered.Wt.vmx.slope))

ggfluctuation(as.table(reordered.Wt.slope.tform), type = "colour") + 
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=.5, size=6, colour="black"),
    axis.text.y = element_text(size=6, colour="black")) + 
  xlab("County") + ylab("Index of Rotated Residuals") + 
  coord_equal() + 
  theme(panel.margin = unit(0, "lines")) + 
  ggtitle(paste("s = ", s))
  
 ggsave(filename = paste("RandomSlope_s", s, ".pdf", sep=""), width = 8,  height = 8, units = "in")

  
ggfluctuation(as.table(reordered.Wt.vmx.slope.tform), type = "colour") + 
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=.5, size=6, colour="black"),
    axis.text.y = element_text(size=6, colour="black")) + 
  xlab("County") + ylab("Index of Rotated Residuals") + 
  coord_equal() + 
  theme(panel.margin = unit(0, "lines")) + 
  ggtitle(paste("s = ", s, "(varimax rotation)"))
  
ggsave(filename = paste("RandomSlope_s", s, "_varimax.pdf", sep=""), width = 8,  height = 8, units = "in")

