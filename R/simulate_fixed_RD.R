#-------------------------------------------------------------------------------
# Script simulating from the model fit to the radon data
#-------------------------------------------------------------------------------

### Packages
library(lme4)
library(HLMdiag)

### Loading data
radon <- read.csv("full_rank_radon.csv")

### Function used to simulate models
sim.hlm <- function(.mod, nsim, e.dsn, b0.dsn, b1.dsn, sigma.err, sigma.b0, sigma.b1) {
  vc <- VarCorr( .mod )
	
	dims <- .mod@dims
	n <- dims[["n"]]
	m <- dims[["q"]] / dims[["nt"]]
	
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
	if(b0.dsn == "norm") {
		b0  <- rnorm(n = nsim * m, mean = 0, sd = sigma.b0)
	} 
	if(b0.dsn == "t") {
		b0  <- (sigma.b0 / sqrt(3)) * rt(n = nsim * m, df = 3)
	}
	if(b0.dsn == "exp") {
		b0  <- sigma.b0 * ( rexp(n = nsim * m) - 1 )
	}
	b0 <- matrix(b0, nc = nsim)

	## Simulating random slope
	if(b1.dsn == "norm") {
		b1  <- rnorm(n = nsim * m, mean = 0, sd = sigma.b1)
	} 
	if(b1.dsn == "t") {
		b1  <- (sigma.b1 / sqrt(3)) * rt(n = nsim * m, df = 3)
	}
	if(b1.dsn == "exp") {
		b1  <- sigma.b1 * ( rexp(n = nsim * m) - 1 )
	}
	b1 <- matrix(b1, nc = nsim)
	
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

## Variance structures
sigma.err <- c(1, 2, 1)
sigma.b   <- c(1, 1, 2)

## fitted model using ML
(fm.ml <- lmer(log.radon ~ basement + house.uranium + uranium + (1 | county) + (0 + house.uranium | county), data = radon, REML = FALSE))

## Distributional settings
e.dsn <- c("norm", "t", "exp")
b.dsn <- c("norm", "t", "exp")
dsn.combos <- expand.grid(e.dsn, b.dsn, KEEP.OUT.ATTRS = FALSE)
colnames(dsn.combos) <- c("e.dsn", "b.dsn")

#-------------------------------------------------------------------------------
# Simulate data and refit the models - you need a sim_models folder
#-------------------------------------------------------------------------------

## fitted model using REML
(fm.reml <- lmer(log.radon ~ basement + house.uranium + uranium + (1 | county) + (0 + house.uranium | county), data = radon, REML = TRUE))

## Set seed and warning handler
set.seed(7231985)
options(warn = 2) # turns warnings into errors

## Simulating
for(v in 1:3) {
for(i in 1:nrow(dsn.combos)) {
	# Simulate
	sim.fm <- sim.hlm( fm.ml, nsim = 1000, e.dsn = dsn.combos[,"e.dsn"][i], 
	                   b0.dsn = dsn.combos[,"b.dsn"][i], 
	                   b1.dsn = dsn.combos[,"b.dsn"][i],
	                   sigma.err = sigma.err[v], sigma.b0 = sigma.b[v], 
	                   sigma.b1 = sigma.b[v] )
                   	   
	# Excluding those simulations with convergence issues 
	# and replacing them with new simulations

	repeat{
  		sim.models <- lapply(sim.fm, 
                       function(x) try(refit(object = fm.reml, newresp = x)))
  
  		conv.issues <- which(sapply(sim.models, function(x) class(x) == "try-error"))
  
  		if(length(conv.issues) == 0) break
  
  		# remove simulations with convergence issues
  		sim.fm <- sim.fm[,-conv.issues]
  
  		if(length(sim.fm) == 1000) break
  
  		# Additional simulations to achieve correct number of simulations
  		sim.fm <- cbind( sim.fm, sim.hlm(fm.ml, nsim = length(conv.issues), 
  		                         e.dsn = dsn.combos[,"e.dsn"][i], 
	                             b0.dsn = dsn.combos[,"b.dsn"][i], 
	                             b1.dsn = dsn.combos[,"b.dsn"][i],
	                             sigma.err = sigma.err[v], sigma.b0 = sigma.b[v], 
	                   			 sigma.b1 = sigma.b[v] ))
	}
	
	saveRDS(sim.models, paste("sim_models/sige", sigma.err[v], "_sigb", sigma.b[v],"/", paste(dsn.combos[,"e.dsn"][i], dsn.combos[,"b.dsn"][i], "REMLmodels.RDS", sep = "_"), sep = ""))
}
}

#-------------------------------------------------------------------------------
# Calculating simulated rotated random effects - you need a residuals folder
# and a folder for each s
#-------------------------------------------------------------------------------

# Vector of files in sim_models folder
REMLmodels <- c("exp_exp_REMLmodels.RDS", "exp_norm_REMLmodels.RDS",  	"exp_t_REMLmodels.RDS", "norm_exp_REMLmodels.RDS", "norm_norm_REMLmodels.RDS", 
	"norm_t_REMLmodels.RDS", "t_exp_REMLmodels.RDS", "t_norm_REMLmodels.RDS",
	 "t_t_REMLmodels.RDS")

s <- NULL ### CHANGE THIS FOR EACH S

var.settings <- c("sige2_sigb1", "sige1_sigb1", "sige1_sigb2")

ngrps <- 60
L.b0  <- kronecker(Diagonal(ngrps), c(1, 0))
L.b1  <- kronecker(Diagonal(ngrps), c(0, 1))

for(v in seq(3L)) {

b0.mcresid <- b0.mcresid.vmx <- vector("list", length = 9L)
b1.mcresid <- b1.mcresid.vmx <- vector("list", length = 9L)

for(i in seq(9L)) {
	b0.mcresid[[i]] <- lapply(models, rotate_ranef, s = s, .L = L.b0)
	b0.mcresid.vmx[[i]] <- lapply(models, rotate_ranef, s = s, .L = L.b0, .varimax = TRUE)

	b1.mcresid[[i]] <- lapply(models, rotate_ranef, s = s, .L = L.b1)
	b1.mcresid.vmx[[i]] <- lapply(models, rotate_ranef, s = s, .L = L.b1, .varimax = TRUE)
	
	rm(models)
}

if(is.null(s)) {
save(b0.mcresid, b0.mcresid.vmx, b1.mcresid, b1.mcresid.vmx, file = paste("residuals/srankB/minconf_level2_residuals_", var.settings[v], ".RData", sep = ""))
	
} else{
save(b0.mcresid, b0.mcresid.vmx, b1.mcresid, b1.mcresid.vmx, file = paste("residuals/s", s,"/minconf_level2_residuals_", var.settings[v], ".RData", sep = ""))
}

}

#-------------------------------------------------------------------------------
# Testing normality of the simulated residuals
#-------------------------------------------------------------------------------

# Vector of files in sim_models folder
REMLmodels <- c("exp_exp_REMLmodels.RDS", "exp_norm_REMLmodels.RDS",  	"exp_t_REMLmodels.RDS", "norm_exp_REMLmodels.RDS", "norm_norm_REMLmodels.RDS", 
	"norm_t_REMLmodels.RDS", "t_exp_REMLmodels.RDS", "t_norm_REMLmodels.RDS",
	 "t_t_REMLmodels.RDS")

# Distribution settings
dist.combos <- t(sapply(strsplit(REMLmodels, "_"), 
				FUN = function(x) return(c(error = x[1], `ranef` = x[2]))))


#-------------------------------------------------------------------------------
# Testing normality of the simulated residuals
#-------------------------------------------------------------------------------

source('normality_functions.R')

### CHANGE THESE INPUTS TO TEST DIFFERENT SIMULATION SETTINGS
s <- 30                      ## EITHER NULL, 55, 50, 45, 40, 35, or 30
var.setting <- "sige1_sigb2" ## EITHER "sige2_sigb1", "sige1_sigb1", or "sige1_sigb2"

if(is.null(s)) {
load(paste("residuals/srankB/minconf_level2_residuals_", var.setting, ".RData", sep = ""))
} else{
load(paste("residuals/s", s,"/minconf_level2_residuals_", var.setting, ".RData", sep = ""))	
}


### Random intercept
b0.mcr          <- test.simulation.results(sims = b0.mcresid, settings = dist.combos)
b0.mcr.vmx      <- test.simulation.results(sims = b0.mcresid.vmx, settings = dist.combos)

b0 <- merge(x = b0.mcr, y = b0.mcr.vmx, by = c("error", "ranef", "alpha"))
print.format(b0, ranef=TRUE)

### Random slope

### CHANGE THESE INPUTS TO TEST DIFFERENT SIMULATION SETTINGS
s <- 30                      ## EITHER NULL, 55, 50, 45, 40, 35, or 30
var.setting <- "sige1_sigb2" ## EITHER "sige2_sigb1", "sige1_sigb1", or "sige1_sigb2"

if(is.null(s)) {
load(paste("residuals/srankB/minconf_level2_residuals_", var.setting, ".RData", sep = ""))
} else{
load(paste("residuals/s", s,"/minconf_level2_residuals_", var.setting, ".RData", sep = ""))	
}

b1.mcr          <- test.simulation.results(sims = b1.mcresid, settings = dist.combos)
b1.mcr.vmx      <- test.simulation.results(sims = b1.mcresid.vmx, settings = dist.combos)

b1 <- merge(x = b1.mcr, y = b1.mcr.vmx, by = c("error", "ranef", "alpha"))
print.format(b1, ranef=TRUE)