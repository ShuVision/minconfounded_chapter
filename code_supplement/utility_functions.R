#-------------------------------------------------------------------------------
# Common functions useful for any simulation study with HLM residuals.
# 
# March 2013
# Adam Loy
#-------------------------------------------------------------------------------

## Function to reorganize Z 
## (This is not needed with the development version of lme4)
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

## Function calculating the marginal variance of the random effects
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