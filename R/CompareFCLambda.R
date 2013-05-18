### Comparing Lambda^-1/2 to FC

### Preliminaries
library(lme4)
library(MASS)
library(RcppEigen)
library(inline)
library(nortest)
library(ggplot2)

setwd("~/Documents/Thesis/Dissertation/eresids-chapter/minconfounded_chapter")

# Reading in the data
radon <- read.csv("data/original_radon.csv")

### Loading self-written functions
source('~/Documents/Thesis/Dissertation/eresids-chapter/simulations/functions/cpp_functions.R')
source('~/Documents/Thesis/Dissertation/eresids-chapter/simulations/functions/resid_functions.R')
source('~/Documents/Thesis/Dissertation/eresids-chapter/simulations/functions/utility_functions.R')

### Fitted model
fm <- lmer(log.radon ~ basement + uranium + (basement | county), data = radon)

ngrps <- summary(fm)@ngrps
L.b0  <- kronecker(Diagonal(ngrps), c(1, 0))
L.b1  <- kronecker(Diagonal(ngrps), c(0, 1))

### Initial fraction of confounding
FC(fm, L.b0)  # 0.72
FC(fm, L.b1)  # 0.70

### Initial trace values
tr(fm, L.b0)  # 61.53067
tr(fm, L.b1)  # 59.5131

### Reducing the trace
reduced.tr <- data.frame(s = seq(1, 85, by = 1), b0 = NA, b1 = NA)

for(i in seq( nrow(reduced.tr) )){
	reduced.tr[i, "b0"] <- tr2(.mod = fm, .L = L.b0, s = reduced.tr[i, "s"])
	reduced.tr[i, "b1"] <- tr2(.mod = fm, .L = L.b1, s = reduced.tr[i, "s"])
}

### Reducing the fraction of confounding
reduced.fc <- reduced.tr
reduced.fc[,2:3] <- reduced.fc[,2:3] / 85

### Taking elements from Lambda
Lambda <- function(.mod, .L) {
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
	
	B.svd <- svd(B)
	
	return(B.svd$d)
}


lambda.b0 <- Lambda(fm, L.b0)
lambda.b1 <- Lambda(fm, L.b1)

### Comparing the two
qplot(x = 1:85, y = )