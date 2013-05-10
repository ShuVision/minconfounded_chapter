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

# Simulated level-1 residuals
sim.e <- llply(sim.mod, resid)                            ## a list of residuals
sim.e        <- melt( do.call("cbind", sim.e) )[,-1]      ## changing to a data frame
names(sim.e) <- c("sample", "residual")                   ## setting colnames for faceting
sim.e.df <- rbind(sim.e, 
			cbind(sample = "true", residual = e))         ## adding the true residuals

# Simulated random intercepts
sim.b0 <- llply(sim.mod, function(x) ranef(x)[[1]][,1])   ## a list of random intercepts
sim.b0 <- melt( do.call("rbind", sim.b0) )[,-2]           ## changing to a data frame
names(sim.b0) <- c("sample", "intercept")                 ## setting colnames for faceting
sim.b0        <- arrange(sim.b0, sample)                  ## ordering by simulation
sim.b0.df <- rbind(sim.b0,
			 cbind(sample = "true", intercept = b[,1]))   ## adding the true r. intercepts
sim.b0.df$intercept <- as.numeric(sim.b0.df$intercept)    ## making residuals numeric

### We are really interested in a lineup for the random slopes for our paper
# Simulated random slopes
sim.b1 <- llply(sim.mod, function(x) ranef(x)[[1]][,2])   ## a list of random slopes
sim.b1 <- melt( do.call("rbind", sim.b1) )[,-2]           ## changing to a data frame
names(sim.b1) <- c("sample", "slope")                     ## setting colnames for faceting
sim.b1        <- arrange(sim.b1, sample)                  ## ordering by simulation
sim.b1.df <- rbind(sim.b1,
			 cbind(sample = "true", slope = b[,2]))       ## adding the true r. slopes
sim.b1.df$slope <- as.numeric(sim.b1.df$slope)            ## making residuals numeric

# -----------------------------------------------------------------------------------------
### Calculations for a weighted Q-Q plot
# -----------------------------------------------------------------------------------------

library(ggplot2)

# Marginal standard devation for original model
msd <- lev2.marginal.var(fm) ## ignore the notes

Fn.b0 <- wecdf(b[,1] / msd[,1], weights = msd[,1]^2)  ## create the weighted CDF function
p.b0  <- Fn.b0(b[,1] / msd[,1])                       ## calc. probs
qplot(x = qnorm(p.b0), y = b[,1] / msd[,1])

Fn.b1 <- wecdf(b[,2] / msd[,2], weights = msd[,1]^2)  ## create the weighted CDF function
p.b1  <- Fn.b0(b[,2] / msd[,2])                       ## calc. probs
qplot(x = qnorm(p.b1), y = b[,2] / msd[,2])


### Incorporating these into a simulation

# standardized random effects
b0.std <- b[,1] / msd[,1]
b1.std <- b[,2] / msd[,2]

# Simulation for weighted Q-Q plots for random intercepts
sim.b0 <- llply(sim.mod, function(x) {
	se <- lev2.marginal.var(x)[,1]
	std.b0 <- ranef(x)[[1]][,1] / se
	Fn <- wecdf(std.b0, weights = se^2)
	p  <- Fn(std.b0)
	qs <- qnorm(p)
	RVAL <- data.frame(intercept = std.b0, quants = qs)  ## list of data frames for plotting
	return(RVAL)
} )   

sim.b0 <- do.call("rbind", sim.b0)                        ## changing to a single data frame
sim.b0$sample <- str_extract(rownames(sim.b0), "\\d+")    ## adding col. for sim. number
sim.b0.df <- rbind(sim.b0,
			 data.frame(sample = "true", 
			            intercept = b0.std,
			            quants = qnorm(p.b0)))            ## adding the true r. intercepts

# Simulation for weighted Q-Q plots for random slopes
sim.b0 <- llply(sim.mod, function(x) {
	se <- lev2.marginal.var(x)[,1]
	std.b0 <- ranef(x)[[1]][,1] / se
	Fn <- wecdf(std.b0, weights = se^2)
	p  <- Fn(std.b0)
	qs <- qnorm(p)
	RVAL <- data.frame(intercept = std.b0, quants = qs)  ## list of data frames for plotting
	return(RVAL)
} )   

sim.b0 <- do.call("rbind", sim.b0)                        ## changing to a single data frame
sim.b0$sample <- str_extract(rownames(sim.b0), "\\d+")    ## adding col. for sim. number
sim.b0.df <- rbind(sim.b0,
			 data.frame(sample = "true", 
			            intercept = b0.std,
			            quants = qnorm(p.b0)))            ## adding the true r. intercepts


# Simulation for weighted Q-Q plots for random slopes
sim.b1 <- llply(sim.mod, function(x) {
	se <- lev2.marginal.var(x)[,1]
	std.b1 <- ranef(x)[[1]][,1] / se
	Fn <- wecdf(std.b1, weights = se^2)
	p  <- Fn(std.b1)
	qs <- qnorm(p)
	RVAL <- data.frame(intercept = std.b1, quants = qs)  ## list of data frames for plotting
	return(RVAL)
} )   

sim.b1 <- do.call("rbind", sim.b1)                        ## changing to a single data frame
sim.b1$sample <- str_extract(rownames(sim.b1), "\\d+")    ## adding col. for sim. number
sim.b1.df <- rbind(sim.b1,
			 data.frame(sample = "true", 
			            intercept = b1.std,
			            quants = qnorm(p.b1)))            ## adding the true r. intercepts
