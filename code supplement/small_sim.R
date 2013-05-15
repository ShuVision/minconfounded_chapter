#-------------------------------------------------------------------------------
# Script simulating from the model fit to the radon data
#-------------------------------------------------------------------------------

### Packages
library(lme4)
library(MASS)

### Loading data
radon <- read.csv("original_radon.csv")

### Function used to simulate models
sim.hlm <- function(.mod, nsim, e.dsn) {
  vc <- VarCorr( .mod )
  D  <- as.matrix( bdiag(vc) )
  
  sigma.err <- attr(vc, "sc")
	
	dims <- .mod@dims
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

#-------------------------------------------------------------------------------
# Set up for simulation study
#-------------------------------------------------------------------------------

## fitted model
fm <- lmer(log.radon ~ basement + uranium + (basement | county), data = radon)

## Distributional settings
e.dsn <- c("norm", "t", "exp")

#-------------------------------------------------------------------------------
# Simulate data and refit the models - you need a sim_models folder
#-------------------------------------------------------------------------------

## Set seed and warning handler
set.seed(7231985)
options(warn = 2) # turns warnings into errors

## Simulating
for(i in seq_along(e.dsn)) {
	# Simulate
	sim.fm <- sim.hlm( fm, nsim = 1000, e.dsn = e.dsn[i])
                   	   
	# Excluding those simulations with convergence issues 
	# and replacing them with new simulations

	repeat{
  		sim.models <- lapply(sim.fm, function(x) try(refit(object = fm, newresp = x)))
  
  		conv.issues <- which(sapply(sim.models, function(x) class(x) == "try-error"))
  
  		if(length(conv.issues) == 0) break
  
  		# remove simulations with convergence issues
  		sim.fm <- sim.fm[,-conv.issues]
  
  		if(length(sim.fm) == 1000) break
  
  		# Additional simulations to achieve correct number of simulations
  		sim.fm <- cbind( sim.fm, sim.hlm(fm, nsim = length(conv.issues), 
  		                         e.dsn = e.dsn[i])
	}
	
	saveRDS(sim.models, paste(e.dsn[i], "_models.RDS", sep = ""))
}



