#-------------------------------------------------------------------------------
# This script calculates FC for minimally confounded residuals for models fit to 
# a subset of the radon data.
#
# Adam Loy
# April 2013
#-------------------------------------------------------------------------------

library(lme4)
library(HLMdiag)
library(plyr)
library(inline)
library(RcppEigen)
library(MASS)

REMLmodels <- c("exp_exp_REMLmodels.RDS", "exp_norm_REMLmodels.RDS",  "exp_t_REMLmodels.RDS", 
"norm_exp_REMLmodels.RDS", "norm_norm_REMLmodels.RDS", "norm_t_REMLmodels.RDS",
"t_exp_REMLmodels.RDS", "t_norm_REMLmodels.RDS", "t_t_REMLmodels.RDS")

source("cpp_functions.R")
source("resid_functions.R")
source("utility_functions.R")


FC2 <- function(.mod, .L, s = NULL, .varimax = FALSE){
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
           A <- crossprod( pzdl )
           B <- t(.L) %*% D %*% t(Z) %*% P %*% Z %*% D %*% .L ## diagnostic se
           W <- try( HLMdiag:::mcrotate(A, B, s) )
           if( class(W) == "try-error") {W <- NA} else {W <- as.matrix(W)}

           if( .varimax == TRUE) {
           	 W <- try( varimax(W, normalize = FALSE)$loadings )
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

#-------------------------------------------------------------------------------
# Calculating simulated residuals
#-------------------------------------------------------------------------------

### CHANGE THIS INPUTS TO TEST DIFFERENT SIMULATION SETTINGS
s <- NULL ### EITHER NULL, 55, 50, 45, 40, 35, or 30

var.settings <- c("sige2_sigb1", "sige1_sigb1", "sige1_sigb2")

ngrps <- 60
L.b0  <- kronecker(Diagonal(ngrps), c(1, 0))
L.b1  <- kronecker(Diagonal(ngrps), c(0, 1))

for(v in seq(3L)) {

cat(var.settings[v], "\n")

b0.fc <- b0.fc.vmx <- vector("list", length = 9L)
b1.fc <- b1.fc.vmx <- vector("list", length = 9L)

for(i in seq(9L)) {
	cat(REMLmodels[i], "\n")
	models <- readRDS(paste("sim_models", var.settings[v], REMLmodels[i], sep = "/"))

	b0.fc[[i]] <- sapply(models, FC2, s = s, .L = L.b0)	
	b0.fc.vmx[[i]] <- sapply(models, FC2, s = s, .L = L.b0, .varimax = TRUE)

	b1.fc[[i]] <- sapply(models, FC2, s = s, .L = L.b1)
	b1.fc.vmx[[i]] <- sapply(models, FC2, s = s, .L = L.b1, .varimax = TRUE)
	
	rm(models)
}

if(is.null(s)) {
save(b0.fc, b0.fc.vmx, b1.fc, b1.fc.vmx, file = paste("fc/srankB/minconf_level2_fc_", var.settings[v], ".RData", sep = ""))
	
} else{
save(b0.fc, b0.fc.vmx, b1.fc, b1.fc.vmx, file = paste("fc/s", s,"/minconf_level2_fc_", var.settings[v], ".RData", sep = ""))
}

}