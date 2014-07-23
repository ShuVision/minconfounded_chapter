### Simulation for plot to select the s-dimensional space

### Preliminaries
library(lme4.0)
library(ggplot2)
library(inline)
library(RcppEigen)
library(MASS)
library(HLMdiag)

source('cpp_functions.R')
source('utility_functions.R')

tr2 <- function(.mod, .L, s = NULL){
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
           
           if(is.matrix(W)){
           		A <- as.matrix(A)
           		B <- as.matrix(B)
           		return( sum( diag( solve( t(W) %*% B %*% W ) %*% ( t(W) %*% A %*% W ) ) ) )
           } else{
           	 return(NA)
           }           
           
         }


### Simulating data
set.seed(6132013)
m  <- 50       ## no. of groups
ni <- 30       ## no. obs. per group
N  <- ni * m   ## total no. of obs.

sig.e <- 3
sig.b <- 1

b <- rnorm(m, 0, sig.b)
e <- rnorm(N, 0, sig.e)

grp.var <- rnorm(m)

sim.df <- data.frame(group = rep(1:m, times = ni), grp.var = rep(grp.var, times = ni), b = rep(b, times = ni), e)

# Introducing some small groups
n.small.grps <- 10
which.groups <- sample.int(m, size = n.small.grps)
size <- 5

for(i in which.groups) {
	index <- which(sim.df$group == i)
	#eliminate <- sample(index, size = length(index) - size)
	eliminate <- sample(index, size = length(index) - rpois(1, size))
	sim.df <- sim.df[-eliminate,]
}

sim.df$y <- with(sim.df, 1 + .3 * grp.var + b + e)

mod <- lmer(y ~ grp.var + (1 | group), data = sim.df)

### plotting
### Reducing the trace
reduced.tr <- data.frame(s = seq(1, m, by = 1), b0 = NA)

for(i in seq( nrow(reduced.tr) )){
	reduced.tr[i, "b0"] <- tr2(.mod = mod, .L = Diagonal(m), s = reduced.tr[i, "s"])
}

qplot(x = m-s, y = b0/m, data = reduced.tr[-m,], geom = c("point", "line")) + 
	xlab("dimension reduction") + 
	ylab("fraction of confounding") + 
	theme_bw() + 
	ylim(c(0, .31))

### Simulating an example -- a less obvious elbow
m  <- 50       ## no. of groups
ni <- c(rpois(40, 30), rpois(10, 5) + 1)        ## no. obs. per group
N  <- sum(ni)   ## total no. of obs.

sig.e <- 3
sig.b <- 1

b <- rnorm(m, 0, sig.b)
e <- rnorm(N, 0, sig.e)

#e[1:40] <- rnorm(40, 0, 3*sig.e)

grp.var <- rnorm(m)

sim.df <- data.frame(group = rep(1:m, times = ni), grp.var = rep(grp.var, times = ni), b = rep(b, times = ni), e)

sim.df$y <- with(sim.df, 1 + .3 * grp.var + b + e)

mod <- lmer(y ~ grp.var + (1 | group), data = sim.df)

### plotting
### Reducing the trace
reduced.tr2 <- data.frame(s = seq(1, m, by = 1), b0 = NA)

for(i in seq( nrow(reduced.tr) )){
	reduced.tr2[i, "b0"] <- tr2(.mod = mod, .L = Diagonal(m), s = reduced.tr[i, "s"])
}

qplot(x = m - s, y = b0 / m , data = reduced.tr2, geom = c("point", "line")) + 
	xlab("dimension reduction") + 
	ylab("fraction of confounding") + 
	theme_bw() + 
	ylim(c(0, .31))
