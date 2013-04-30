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