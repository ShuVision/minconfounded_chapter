#-------------------------------------------------------------------------------
# Script to render Figure 6 heatmaps
#-------------------------------------------------------------------------------

library(lme4)
library(ggplot2)
library(inline)
library(RcppEigen)
library(MASS)

source('cpp_functions.R')
source('utility_functions.R')


nclass <- 20

set.seed(123)
nstud <- c(rpois(15, 15), rpois(5, 5))

class.effect <- rnorm(nclass)*2

student.effect <- rnorm(sum(nstud))*3

# Data we need in order to generate our analysis.
data.set = data.frame(class.id = rep(1:nclass, time=nstud),
                  class.effect = rep(class.effect, times=nstud),
                  student.id = 1:(sum(nstud)),
                  student.effect = student.effect)
data.set$outcomes = data.set$class.effect + data.set$student.effect

head(data.set)
fm <- lmer(outcomes ~ 1 + (1|class.id), data=data.set)


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
  index <- sort(index[index >= 0])
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
            geom_tile(colour = "white") + scale_fill_gradient2(expression(w[ij])
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
# Di <- Diagonal(x = c(vc[[1]], vc[[2]])) / (unname(attr(vc, "sc")))^2
Di <- bdiag(vc) / sigma(fm)^2
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

### Rotating the random intercept
A <-  crossprod( P %*% Z %*% D )
B <- D %*% t(Z) %*% P %*% Z %*% D
s <- rankMatrix(B)

fc.int   <- diag(ginv(as.matrix(B)) %*% A)

W.int <- as.matrix( mcrotate(A = A, B = B, s = s) )
W.vmx.int <- varimax( W.int, normalize = FALSE )$loadings

Wt.int <- as.matrix( zapsmall( t(W.int) ) )
Wt.vmx.int <- as.matrix( zapsmall( t(W.vmx.int) ) )

# colnames(Wt.int) <- colnames(Wt.vmx.int) <- unique( radon$county.name )
colnames(Wt.int) <- colnames(Wt.vmx.int) <- unique( fm@frame$class.id )
rownames(Wt.int) <- rownames(Wt.vmx.int) <- seq_len(nrow(Wt.int))

# Pulling off variances of random ints
rint.var <- diag(B) # we could order by this, but fc makes more sense

reordered.Wt.int <- Wt.int[, order(fc.int, decreasing = TRUE)]
reordered.Wt.vmx.int <- Wt.vmx.int[, order(fc.int, decreasing = TRUE)]

reordered.Wt.int.tform <- sign(reordered.Wt.int) * sqrt(abs(reordered.Wt.int))
reordered.Wt.vmx.int.tform <- sign(reordered.Wt.vmx.int) * sqrt(abs(reordered.Wt.vmx.int))


ggfluctuation(as.table(reordered.Wt.int), type = "colour") + 
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=.5, size=6, colour="black"),
    axis.text.y = element_text(size=6, colour="black")) + 
  xlab("Group") + ylab("Index of Rotated Residuals") + 
  coord_equal() + 
  theme(panel.margin = unit(0, "lines"))

ggfluctuation(as.table(reordered.Wt.vmx.int), type = "colour") + 
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=.5, size=6, colour="black"),
    axis.text.y = element_text(size=6, colour="black")) + 
  xlab("Group") + ylab("Index of Rotated Residuals") + 
  coord_equal() + 
  theme(panel.margin = unit(0, "lines"))