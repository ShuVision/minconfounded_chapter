#-------------------------------------------------------------------------------
# Common functions useful for any simulation study with HLM residuals.
# 
# March 2013
# Adam Loy
#-------------------------------------------------------------------------------

## function to calculate ("by hand") the level-1 residuals
lev1.resid <- function(.model) {
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
  
  betahat <- solve(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv %*% y
  mr <- y - X %*% betahat
  
  bvec <- D %*% t(Z) %*% Vinv %*% mr
  
  return(as.vector(y - X %*% betahat - Z %*% bvec))
}

## function to calculate ("by hand") the standardized level-1 residuals
lev1.std.resid <- function(.model) {
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
  
  # Calculating P
  XVXinv <- solve( t(X) %*% Vinv %*% X )
  VinvX  <- Vinv %*% X
  M      <- VinvX %*% XVXinv %*% t(VinvX)
  P      <- cxxmatsub(as.matrix(Vinv), as.matrix(M))

  r <- R %*% P %*% y
  
  return(as.vector(r / sqrt(diag(R %*% P %*% R))))
}

## function to calculate ("by hand") the standardized level-1 residuals
## discussed by Pinheiro and Bates.
lev1.pinheiro.resid <- function(.model) {
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
  
  betahat <- solve(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv %*% y
  mr <- y - X %*% betahat
  
  bvec <- D %*% t(Z) %*% Vinv %*% mr
  
  return(as.vector(y - X %*% betahat - Z %*% bvec) / sig0)
}

## function to calculate ("by hand") the marginal residuals
marginal.resid <- function(.model) {
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

  betahat <- solve(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv %*% y

  return(as.vector(y - X %*% betahat))
}

# function to calculate ("by hand") the cholesky residuals
cholesky.resid <- function(.model) {
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
  
  Lt <- chol(Vinv)

  betahat <- solve(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv %*% y
  mr <- y - X %*% betahat

  return(as.vector(Lt %*% mr))
}

## function to calculate ("by hand") the level-2 residuals, i.e. EBLUPS
lev2.resid <- function(.model) {
  y <- .model@y
  X <- getME(.model, "X")
  Z <- BlockZ(.model)
  n <- nrow(X)
  nranef <- ncol(ranef(.model)[[1]])
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
  
  betahat <- solve(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv %*% y
  mr <- y - X %*% betahat
  
  bvec <- D %*% t(Z) %*% Vinv %*% mr
  if(nranef > 1) {
  	bmat <- matrix(bvec, ncol = nranef, byrow = TRUE)
    return(bmat)
  } else {
  	return(bvec)
  }
  
}

## Standardized EBLUPS for use with Lange and Ryan's approach.
## The marginal covariance is used to standardized rather than
## the usual matrix.
lev2.langeryan.resid <- function(.model) {
  y <- .model@y
  X <- getME(.model, "X")
  Z <- BlockZ(.model)
  n <- nrow(X)
  nranef <- ncol(ranef(.model)[[1]])
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
  
  betahat <- solve(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv %*% y
  mr <- y - X %*% betahat
  
  bvec <- D %*% t(Z) %*% Vinv %*% mr
  #bmat <- matrix(bvec, ncol = 2, byrow = TRUE)

  bse <- crossprod( chol(Vinv) %*% Z %*% D ) # Marginal COV. used by Lange and Ryan
  
  bse.diag <- diag(bse)
  
  rvec <- bvec / sqrt(bse.diag)
  
  if(nranef > 1) {
  	rmat <- matrix(rvec, ncol = nranef, byrow = TRUE)
    return(rmat)
  } else {
  	return(rvec)
  }

}

## EBLUPs standardized by Goldstein's method
lev2.goldstein.resid <- function(.model) {
  y <- .model@y
  X <- getME(.model, "X")
  Z <- BlockZ(.model)
  n <- nrow(X)
  nranef <- ncol(ranef(.model)[[1]])
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
  
  betahat <- solve(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv %*% y
  mr <- y - X %*% betahat
  bvec <- D %*% t(Z) %*% Vinv %*% mr

  XVXinv <- solve( t(X) %*% Vinv %*% X )
  VinvX  <- Vinv %*% X
  M      <- VinvX %*% XVXinv %*% t(VinvX)
  P      <- cxxmatsub(as.matrix(Vinv), as.matrix(M))
  
  bse <- D %*% t(Z) %*% P %*% Z %*% D # Diagnostic COV used by Goldstein
  bse.diag <- diag(bse)
  
  rvec <- bvec / sqrt(bse.diag)

  if(nranef > 1) {
  	rmat <- matrix(rvec, ncol = nranef, byrow = TRUE)
    return(rmat)
  } else {
  	return(rvec)
  }
}