################################################################################
# Simulation script exploring the impact of misspecifying the distribution
# of the random effects in a linear mixed-effects model.
#
# Adam Loy
# June 2014
################################################################################

### Prelimiaries
library(lme4)
library(plyr)
library(MASS)

#-------------------------------------------------------------------------------
# Simulation functions
#-------------------------------------------------------------------------------

# Function to extract pieces of interest
extract <- function(model) {
  fixef <- fixef(model) # getME(model, "beta")
  varcomp <- as.data.frame(VarCorr(model))$sdcor
  blups <- ranef(model)[[1]]
  return(list(fixef = fixef, varcomp = varcomp, blups = blups))
}

# squared deviation
se <- function(actual, predicted) {
  (actual - predicted)^2
}


# Summarizing a simulation where the above pieces were extracted
sim_summary <- function(results, actual) {
  actual_fixef <- actual$fixef
  actual_varcomp <- actual$varcomp
#   actual_blups <- actual$blups
  
  # mean square error of fixef estimates
  mse_fixef <- ldply(results, .fun = function(x) se(actual_fixef, x$fixef))
  mse_fixef <- apply(mse_fixef[,-1], 2, mean)
  
  # bias of fixef estimates
  bias_fixef <- ldply(results, .fun = function(x) x$fixef - actual_fixef)
  bias_fixef <- apply(bias_fixef[,-1], 2, mean)
  
  # sd of fixef estimates
  sd_fixef <- ldply(results, .fun = function(x) x$fixef)
  sd_fixef <- apply(sd_fixef[,-1], 2, sd)
  
  fixef_results <- rbind(mse_fixef, bias_fixef, sd_fixef)
  rownames(fixef_results) <- c("MSE", "Bias", "SD")
  
  # mean square error of variance component estimates
  mse_vc <- ldply(results, .fun = function(x) se(actual_varcomp, x$varcomp))
  mse_vc <- apply(mse_vc[,-1], 2, mean)
  
  # bias square error of variance component estimates
  bias_vc <- ldply(results, .fun = function(x) x$varcomp - actual_varcomp)
  bias_vc <- apply(bias_vc[,-1], 2, mean)
  
  # sd square error of variance component estimates
  sd_vc <- ldply(results, .fun = function(x) x$varcomp)
  sd_vc <- apply(sd_vc[,-1], 2, sd)
  
  vc_results <- rbind(mse_vc, bias_vc, sd_vc)
  rownames(vc_results) <- c("MSE", "Bias", "SD")
  colnames(vc_results) <- c("sigma_b0", "sigma_b1", "sigma_e")
  
  
  cat("Fixed effects:\n")
  print(fixef_results)
  cat("\nVariance components:\n")
  print(vc_results)
  
}

#-------------------------------------------------------------------------------
# Fit the model that is tha basis for simulation
#-------------------------------------------------------------------------------

radon <- read.csv(file.choose()) # go find "original_rank_radon.csv"

## fitted model using ML
fm <- lmer(log.radon ~ basement + uranium + (basement | county), data = radon)


#-------------------------------------------------------------------------------
# Assuming normal random effects when random effects are normal
#-------------------------------------------------------------------------------

sim_y <- simulate(fm, nsim = 10)
sim_fits <- apply(sim_y, 2, refit, object = fm)
sim_parts <- llply(sim_fits, extract)

#-------------------------------------------------------------------------------
# JCGS AE's suggestion: set the variance components for REs to zero.
#-------------------------------------------------------------------------------

### Function used to simulate models
sim.hlm <- function(.mod, nsim, e.dsn) {  
  sigma.err <- sigma(.mod)
  
  dims <- .mod@devcomp$dims
  n <- dims[["n"]]
  
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
  
  
  ## Generating y with random effects set to zero
  y <- getME(.mod, "X") %*% fixef(.mod)
  y <- aaply(y, 2, function(x) x + e)
  
  y.df <- as.data.frame( as.matrix( y) )
  colnames(y.df) <- paste("sim_", 1:ncol(y.df), sep = "")
  
  return( y.df )
}


## Set seed and warning handler
set.seed(7231985)
options(warn = 2) # turns warnings into errors


sim.fm <- sim.hlm( fm, nsim = 1000, e.dsn = "norm")

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
  sim.fm <- cbind( sim.fm, sim.hlm(fm, nsim = length(conv.issues), e.dsn = "norm")
}

# confints
cis.wald <- llply(sim.models, confint, method = "Wald")

ci.int <- ldply(cis.wald, function(x) x[1,])
sum(ci.int[,1] < fixef(fm)[1] | ci.int[,2] > fixef(fm)[1]) / 1000

ci.basement <- ldply(cis.wald, function(x) x[2,])
sum(ci.basement[,1] < fixef(fm)[1] | ci.basement[,2] > fixef(fm)[1]) / 1000

ci.uranium <- ldply(cis.wald, function(x) x[3,])
sum(ci.uranium[,1] < fixef(fm)[1] | ci.uranium[,2] > fixef(fm)[1]) / 1000


# other summaries
sim2_actual <- extract(fm)
sim2_actual$varcomp <- c(0, 0, sigma(fm))
sim2_parts <- llply(sim.models, extract)
sim_summary(results = sim2_parts, actual = sim2_actual)



#-------------------------------------------------------------------------------
# Assuming normal random effects when random effects are t
#-------------------------------------------------------------------------------

