#-------------------------------------------------------------------------------
# Script for simulating the null distribution for using the Lange-Ryan
# residuals to test for normality.
#
# Adam Loy
# June/July 2014
#-------------------------------------------------------------------------------

library(lme4.0)
library(lme4)
library(plyr)
library(MASS)
library(nortest)
library(mnormt)

### Funtion for the marginal variance
lev2.marginal.var <- function(.model) {
  y <- getME(.model, "y")
  X <- getME(.model, "X")
  Z <- BlockZ(.model)
  n <- nrow(X)
  ngrps <- unname(sapply(.model@flist, function(x) length(levels(x))))
  
  # Constructing V = Cov(Y)
  sig0 <- attr(VarCorr(.model), "sc") # sigma(.model)
  
  ZDZt <- sig0^2 * crossprod( getME(.model, "A") )
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

### Function for proper arrangement of the Z matrix
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


### Function for Lange and Ryan's residuals
lev2.langeryan.resid <- function(.model) {
  y <- getME(.model, "y")
  X <- getME(.model, "X")
  Z <- BlockZ(.model)
  n <- nrow(X)
  nranef <- ncol(ranef(.model)[[1]])
  ngrps <- unname(sapply(.model@flist, function(x) length(levels(x))))
  
  # Constructing V = Cov(Y)
  sig0 <- attr(VarCorr(.model), "sc") # sigma(.model)
  
  ZDZt <- sig0^2 * crossprod( getME(.model, "A") )
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

# Function to simulate from the HLM
sim.hlm <- function(.mod, nsim, e.dsn) {
  vc <- VarCorr( .mod )
  D  <- as.matrix( bdiag(vc) )
  
  sigma.err <- attr(vc, "sc")
	
	dims <- .mod@devcomp$dims
	n <- dims[["n"]]
	m <- dims[["q"]] / nrow(D)
	
	## Simulating error terms
	if(e.dsn == "norm") {
		e  <- rnorm(n = nsim * n, mean = 0, sd = sigma.err)
	} 
	if(e.dsn == "t") {
		e  <- (sigma.err / sqrt(3)) * rt(n = nsim * n, df = 3)
	}
	if(e.dsn == "exp") {
		e  <- sigma.err * ( rexp(n = nsim * n) - 1 )
	}
	e <- matrix(e, nc = nsim)
	
	## Simulating random intercept
	b.sim <- mvrnorm(n = nsim * m, mu = c(0, 0), Sigma = D)
	
	b0 <- matrix(b.sim[,1], nc = nsim)
	b1 <- matrix(b.sim[,2], nc = nsim)
	
	## Generating y
	b <- rbind(b0, b1)
	y <- getME(.mod, "X") %*% fixef(.mod) + getME(.mod, "Z") %*% b + e
	
	y.df <- as.data.frame( as.matrix( y) )
	colnames(y.df) <- paste("sim_", 1:ncol(y.df), sep = "")
	
	return( y.df )
}

source("normality_functions.R")



### Fitting a model with random intercept for county and 
### random slope for basement.

radon <- read.csv("original_radon.csv") # in code supplement
fm <- lmer(log.radon ~ basement + uranium + (basement | county), data = radon, REML = FALSE)


### Simulating Lange and Ryan's residuals under a properly specified model.
### This 1000 simulations is to establish the null distribution for a 
### bootstrap test.

null_sims <- simulate(fm, nsim = 1000, seed = 9291632)
null_refit <- lapply(null_sims, refit, object = fm)
null_lr_resids <- lapply(null_refit, lev2.langeryan.resid)

b0.null.test.stats <- sapply(null_refit, test.statistic.langeryan, column = 1)
b1.null.test.stats <- sapply(null_refit, test.statistic.langeryan, column = 2)



### Simulations to check properties of Lange and Ryan's residuals

# Normal error terms
# set.seed(7231985)
# normal_sims <- sim.hlm(fm, nsim = 1000, e.dsn = "norm")
# normal_refit <- lapply(normal_sims, refit, object = fm)
normal_refit <- readRDS("~/Documents/eresids-chapter/simulations/paper_sim_study/small_sim/norm_models.RDS")
normal_lr_resids <- lapply(normal_refit, lev2.langeryan.resid)
b0_normal_lr_stats  <- sapply(normal_refit, test.statistic.langeryan, column = 1)
b1_normal_lr_stats  <- sapply(normal_refit, test.statistic.langeryan, column = 1)

# Heavy-tailed error terms
# set.seed(7231985)
# heavy_sims <- sim.hlm(fm, nsim = 1000, e.dsn = "t")
# heavy_refit <- lapply(heavy_sims, refit, object = fm)
heavy_refit <- readRDS("~/Documents/eresids-chapter/simulations/paper_sim_study/small_sim/t_models.RDS")
heavy_lr_resids <- lapply(heavy_refit, lev2.langeryan.resid)
b0_heavy_lr_stats  <- sapply(heavy_refit, test.statistic.langeryan, column = 1)
b1_heavy_lr_stats  <- sapply(heavy_refit, test.statistic.langeryan, column = 1)


# Skewed error terms
# set.seed(7231985)
# skewed_sims <- sim.hlm(fm, nsim = 1000, e.dsn = "exp")
# skewed_refit <- lapply(skewed_sims, refit, object = fm)
skewed_refit <- readRDS("~/Documents/eresids-chapter/simulations/paper_sim_study/small_sim/exp_models.RDS")
skewed_lr_resids <- lapply(skewed_refit, lev2.langeryan.resid)
b0_skewed_lr_stats  <- sapply(skewed_refit, test.statistic.langeryan, column = 1)
b1_skewed_lr_stats  <- sapply(skewed_refit, test.statistic.langeryan, column = 1)


### Running the tests for normality
b0_lr_stats <- list(b0_normal_lr_stats, b0_heavy_lr_stats, b0_skewed_lr_stats)
dist.combos <- matrix(c("normal", "t", "exp", rep("normal", 3)), ncol = 2, byrow = F)
colnames(dist.combos) <- c("error", "ranef")
langeryan.test(test.stats = b0_lr_stats, settings = dist.combos, 
               null.dsn = b0.null.test.stats)
