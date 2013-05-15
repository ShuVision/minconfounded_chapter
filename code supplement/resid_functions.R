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

## function to calculate ("by hand") the level-2 residuals
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

## Lange and Ryan stuff -- these are standardized level-2 residuals,
## but the marginal covariance is used to standardized rather than
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

## Level-2 residuals standardized by Goldstein's method
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


## Level-1 minimally confounded residuals
lev1.mcresid <- function(.model) {
  y <- .model@y
  X <- getME(.model, "X")
  Z <- BlockZ(.model)
  n <- nrow(X)
  p <- ncol(X)
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
  
  r1 <- as.vector(y - X %*% betahat - Z %*% bvec)
  
  XVXinv <- solve( t(X) %*% Vinv %*% X )
  VinvX  <- Vinv %*% X
  M      <- VinvX %*% XVXinv %*% t(VinvX)
  P      <- cxxmatsub(as.matrix(Vinv), as.matrix(M))

  P.svd <- svd(P)
  lambda <- P.svd$d
  U <- P.svd$u

  RVAL <- rep(NA, n - p)
  for(i in 1:(n - p)) {
  	RVAL[i] <- lambda[i] * t(U[,i]) %*% y
  }
  
  return(RVAL / sig0^2)
}


### Function carrying out the Reduction Method of McDonald et al. (1979)
mcrotate <- function(A, B, s) {
  r <- rankMatrix(B)
  if(is.null(s)) s <- r
  
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

mcrotate.alt <- function(A, B, s = NULL) {
  r <- rankMatrix(B)
  if(is.null(s)) s <- r
  
  B.svd <- svd(B)
  Cr.diag <- B.svd$d[1:r]
  Tr <- B.svd$u[, 1:r]
  
  A.star <- Diagonal( x = 1 / sqrt(Cr.diag) ) %*% t(Tr) %*% A %*% Tr %*% Diagonal( x = 1 / sqrt(Cr.diag) )

#  A.star <- Diagonal(n = length(Cr.diag), x = 1 / mean(sqrt(Cr.diag)) ) %*% t(Tr) %*% A %*% Tr %*% Diagonal(n = length(Cr.diag), x = 1 / mean(sqrt(Cr.diag)) )
  
  A.star.svd <- svd( A.star )
  
  index <- seq(r, length.out = s, by = -1)
  index <- sort(index[index >= 0])
  W <- Tr %*% Diagonal(n = length(Cr.diag), x = 1 / mean(sqrt(Cr.diag)) ) %*% A.star.svd$u[,index]
  
  return(W)
}

mcrotate.alt2 <- function(A, B, s = NULL) {
  r <- rankMatrix(B)
  if(is.null(s)) s <- r
  
  B.svd <- svd(B)
  Cr.diag <- B.svd$d[1:r]
  Tr <- B.svd$u[, 1:r]
  
  A.star <- Diagonal( x = 1 / sqrt(Cr.diag) ) %*% t(Tr) %*% A %*% Tr %*% Diagonal( x = 1 / sqrt(Cr.diag) )

#  A.star <- t(Tr) %*% A %*% Tr
  
  A.star.svd <- svd( A.star )

  index <- seq(r, length.out = s, by = -1)
  index <- sort(index[index >= 0])  
  W <- Tr %*% A.star.svd$u[,index]
  
  return(W)
}


mcresid2.alt <- function(.mod, .L, s = NULL, .varimax = FALSE, .qr = FALSE){
	       y <- .mod@y
           X <- getME(.mod, "X")
           Z <- BlockZ(.mod)
           
           n <- nrow(X)
           p <- ncol(X)
           ngrps <- unname( summary(.mod)@ngrps )
           
           vc <- VarCorr(.mod)
           Di <- bdiag( VarCorr(.mod) ) / (unname(attr(vc, "sc")))^2
           D  <- kronecker( Diagonal(ngrps), Di )
           
           Aslot <- .mod@A # ZDZ'
           zdzt <- crossprod( .mod@A )
           V  <- Diagonal( n ) + zdzt
           V.chol <- chol( V )
           Vinv  <- chol2inv( V.chol ) 
           
           XVXinv <- solve( t(X) %*% Vinv %*% X )
           VinvX  <- Vinv %*% X
           M      <- VinvX %*% XVXinv %*% t(VinvX)
           P      <- cxxmatsub(as.matrix(Vinv), as.matrix(M))
           
           betahat <- solve(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv %*% y
           mr <- y - X %*% betahat
  
           bvec <- D %*% t(Z) %*% Vinv %*% mr
                      
           pzdl <- P %*% Z %*% D %*% .L
           A <- crossprod( pzdl )
           B <- t(.L) %*% D %*% t(Z) %*% P %*% Z %*% D %*% .L ## diagnostic se
           W <- try( mcrotate.alt(A, B, s) )
           if( class(W) == "try-error") {W <- NA} else {W <- as.matrix(W)}

           if( .varimax == TRUE) {
           	 W <- try( varimax(W, normalize = FALSE)$loadings )
             if( class(W) == "try-error" ) W <- NA 
           }
           
		   if( .qr == TRUE ) {
             W <- try( qr.Q( qr(W) ) )
             if( class(W) == "try-error" ) W <- NA 
           }
             
           return( as.numeric( t(W) %*% as.numeric( t(.L) %*% bvec ) ) )
         }

mcresid2.alt2 <- function(.mod, .L, s = NULL, .varimax = FALSE){
	       y <- .mod@y
           X <- getME(.mod, "X")
           Z <- BlockZ(.mod)
           
           n <- nrow(X)
           p <- ncol(X)
           ngrps <- unname( summary(.mod)@ngrps )
           
           vc <- VarCorr(.mod)
           Di <- bdiag( VarCorr(.mod) ) / (unname(attr(vc, "sc")))^2
           D  <- kronecker( Diagonal(ngrps), Di )
           
           Aslot <- .mod@A # ZDZ'
           zdzt <- crossprod( .mod@A )
           V  <- Diagonal( n ) + zdzt
           V.chol <- chol( V )
           Vinv  <- chol2inv( V.chol ) 
           
           XVXinv <- solve( t(X) %*% Vinv %*% X )
           VinvX  <- Vinv %*% X
           M      <- VinvX %*% XVXinv %*% t(VinvX)
           P      <- cxxmatsub(as.matrix(Vinv), as.matrix(M))
           
           betahat <- solve(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv %*% y
           mr <- y - X %*% betahat
  
           bvec <- D %*% t(Z) %*% Vinv %*% mr
                      
           pzdl <- P %*% Z %*% D %*% .L
           A <- crossprod( pzdl )
           B <- t(.L) %*% D %*% t(Z) %*% P %*% Z %*% D %*% .L ## diagnostic se
           W <- try( mcrotate.alt2(A, B, s) )
           if( class(W) == "try-error") {W <- NA} else {W <- as.matrix(W)}

           if( .varimax == TRUE) {
           	 W <- try( varimax(W, normalize = FALSE)$loadings )
             if( class(W) == "try-error" ) W <- NA 
           }
             
           return( as.numeric( t(W) %*% as.numeric( t(.L) %*% bvec ) ) )
         }


mcresid2 <- function(.mod, .L, s = NULL, .varimax = FALSE, .standardize = FALSE, .quartimax = FALSE, .qr = FALSE){
	       y <- .mod@y
           X <- getME(.mod, "X")
           Z <- BlockZ(.mod)
           
           n <- nrow(X)
           p <- ncol(X)
           ngrps <- unname( summary(.mod)@ngrps )
           
           vc <- VarCorr(.mod)
           Di <- bdiag( VarCorr(.mod) ) / (unname(attr(vc, "sc")))^2
           D  <- kronecker( Diagonal(ngrps), Di )
           
           Aslot <- .mod@A # ZDZ'
           zdzt <- crossprod( .mod@A )
           V  <- Diagonal( n ) + zdzt
           V.chol <- chol( V )
           Vinv  <- chol2inv( V.chol ) 
           
           XVXinv <- solve( t(X) %*% Vinv %*% X )
           VinvX  <- Vinv %*% X
           M      <- VinvX %*% XVXinv %*% t(VinvX)
           P      <- cxxmatsub(as.matrix(Vinv), as.matrix(M))
           
           betahat <- solve(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv %*% y
           mr <- y - X %*% betahat
  
           bvec <- D %*% t(Z) %*% Vinv %*% mr
           
           if(.standardize == TRUE) {
           		bse <- crossprod( chol(Vinv) %*% Z %*% D ) # Marginal COV. used by Lange and Ryan
 		   		bse.diag <- diag(bse)
                bvec <- bvec / sqrt(bse.diag)
           }
                      
           pzdl <- P %*% Z %*% D %*% .L
           A <- crossprod( pzdl )
           B <- t(.L) %*% D %*% t(Z) %*% P %*% Z %*% D %*% .L ## diagnostic se
           W <- try( mcrotate(A, B, s) )
           if( class(W) == "try-error") {W <- NA} else {W <- as.matrix(W)}

		   if( .qr == TRUE ) {
             W <- try( qr.Q( qr(W) ) )
             if( class(W) == "try-error" ) W <- NA 
           }
               
           if( .varimax == TRUE) {
           	 W <- try( varimax(W, normalize = FALSE)$loadings )
             if( class(W) == "try-error" ) W <- NA 
           }
           
           if( .quartimax == TRUE ){
           	 W <- try( quartimax(W, normalize = FALSE, maxit = 10000)$loadings )
           	 if( class(W) == "try-error" ) W <- NA 
           }
                        
           return( as.numeric( t(W) %*% as.numeric( t(.L) %*% bvec ) ) )
         }


FC <- function(.mod, .L) {
	       y <- .mod@y
           X <- getME(.mod, "X")
           Z <- BlockZ(.mod)
           
           n <- nrow(X)
           p <- ncol(X)
           ngrps <- unname( summary(.mod)@ngrps )
           
           vc <- VarCorr(.mod)
           Di <- bdiag( VarCorr(.mod) ) / (unname(attr(vc, "sc")))^2
           D  <- kronecker( Diagonal(ngrps), Di )
           
           Aslot <- .mod@A # ZDZ'
           zdzt <- crossprod( .mod@A )
           V  <- Diagonal( n ) + zdzt
           V.chol <- chol( V )
           Vinv  <- chol2inv( V.chol ) 
           
           XVXinv <- solve( t(X) %*% Vinv %*% X )
           VinvX  <- Vinv %*% X
           M      <- VinvX %*% XVXinv %*% t(VinvX)
           P      <- cxxmatsub(as.matrix(Vinv), as.matrix(M))
                                 
           pzdl <- P %*% Z %*% D %*% .L
           A <- as.matrix( crossprod( pzdl ) )
           B <- as.matrix( t(.L) %*% D %*% t(Z) %*% P %*% Z %*% D %*% .L ) ## diagnostic se
           
           RVAL <- try( mean( diag( ginv(B) %*% (A) ) ) )
           if( class(RVAL) == "try-error" ) RVAL <- NA
           return( RVAL )
}

FC2 <- function(.mod, .L, s = NULL, .varimax = FALSE, .standardize = FALSE, .quartimax = FALSE, .qr = FALSE){
	       y <- .mod@y
           X <- getME(.mod, "X")
           Z <- BlockZ(.mod)
           
           n <- nrow(X)
           p <- ncol(X)
           ngrps <- unname( summary(.mod)@ngrps )
           
           vc <- VarCorr(.mod)
           Di <- bdiag( VarCorr(.mod) ) / (unname(attr(vc, "sc")))^2
           D  <- kronecker( Diagonal(ngrps), Di )
           
           Aslot <- .mod@A # ZDZ'
           zdzt <- crossprod( .mod@A )
           V  <- Diagonal( n ) + zdzt
           V.chol <- chol( V )
           Vinv  <- chol2inv( V.chol ) 
           
           XVXinv <- solve( t(X) %*% Vinv %*% X )
           VinvX  <- Vinv %*% X
           M      <- VinvX %*% XVXinv %*% t(VinvX)
           P      <- cxxmatsub(as.matrix(Vinv), as.matrix(M))
           
           if(.standardize == TRUE) {
           		bse <- crossprod( chol(Vinv) %*% Z %*% D ) # Marginal COV. used by Lange and Ryan
 		   		bse.diag <- diag(bse)
                bvec <- bvec / sqrt(bse.diag)
           }
                      
           pzdl <- P %*% Z %*% D %*% .L
           A <- crossprod( pzdl )
           B <- t(.L) %*% D %*% t(Z) %*% P %*% Z %*% D %*% .L ## diagnostic se
           W <- try( mcrotate(A, B, s) )
           if( class(W) == "try-error") {W <- NA} else {W <- as.matrix(W)}

		   if( .qr == TRUE ) {
             W <- try( qr.Q( qr(W) ) )
             if( class(W) == "try-error" ) W <- NA 
           }
               
           if( .varimax == TRUE) {
           	 W <- try( varimax(W, normalize = FALSE)$loadings )
             if( class(W) == "try-error" ) W <- NA 
           }
           
           if( .quartimax == TRUE ){
           	 W <- try( quartimax(W, normalize = FALSE, maxit = 10000)$loadings )
           	 if( class(W) == "try-error" ) W <- NA 
           }
           
           if(is.matrix(W)){
           		A <- as.matrix(A)
           		B <- as.matrix(B)
           		return( mean( diag( solve( t(W) %*% B %*% W ) %*% ( t(W) %*% A %*% W ) ) ) )
           } else{
           	 return(NA)
           }           
           
         }