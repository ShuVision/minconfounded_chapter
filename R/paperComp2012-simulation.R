#-------------------------------------------------------------------------------
# Script simulating from the model fit to the radon data.
#
# Adam Loy
# November 2012
#-------------------------------------------------------------------------------

### Preliminaries
library(lme4)
library(arm)
library(nortest)
library(multicore)

# Function to reorganize Z
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

### Reading in the radon data provided by Gelman and Hill (2007)
srrs2 <- read.table ("srrs2.dat", header=T, sep=",")

## Restricting attention to Minnesota
mn <- srrs2$state=="MN"
radon <- srrs2$activity[mn]
log.radon <- log (ifelse (radon==0, .1, radon))
floor <- srrs2$floor[mn]       # 0 for basement, 1 for first floor
n <- length(radon)
y <- log.radon
basement <- floor

## Getting county index variable
county.name <- as.vector(srrs2$county[mn])
uniq <- unique(county.name)
J <- length(uniq)
county <- rep (NA, J)
for (i in 1:J){
  county[county.name==uniq[i]] <- i
}

## Reading the county-level data
srrs2.fips <- srrs2$stfips*1000 + srrs2$cntyfips
cty <- read.table ("cty.dat", header=T, sep=",")
usa.fips <- 1000*cty[,"stfips"] + cty[,"ctfips"]
usa.rows <- match (unique(srrs2.fips[mn]), usa.fips)
uranium <- cty[usa.rows,"Uppm"]
u <- log (uranium)
uranium <- u[county]


### Fitting a model with random intercept for county and 
### random slope for basement.
fm <- lmer(log.radon ~ basement + uranium + (basement | county), REML = FALSE)
display(fm)

### Simulating new responses and refitting the model
set.seed(7231985)

sim.fm <- simulate(fm, nsim = 1000)

## excluding those simulations with convergence issues
options(warn = 2) # turns warnings into errors

repeat{
  sim.models <- lapply(sim.fm, 
                       function(x) try(refit(object = fm, newresp = x)))
  
  conv.issues <- which(sapply(sim.models, function(x) class(x) == "try-error"))
  
  if(length(conv.issues) == 0) break
  
  # remove simulations with convergence issues
  sim.fm <- sim.fm[,-conv.issues]
  
  if(length(sim.fm) == 1000) break
  
  # additional simulations to achieve correct number of simulations
  sim.fm <- cbind(sim.fm, simulate(fm, nsim = length(conv.issues)))
}

### Extracting residuals from the simulations
lev1 <- lapply(sim.models, resid)
rint <- lapply(sim.models, function(x) ranef(x)[[1]][,1])
rslope <- lapply(sim.models, function(x) ranef(x)[[1]][,2])

save(lev1, rint, rslope, file = "rawresids-13Dec12.RData")

### Function to calculate the minimally confounded residuals
### for a specified level
mcrotate <- function(A, B) {
  r <- rankMatrix(B)
  
  B.svd <- svd(B)
  Cr.diag <- B.svd$d[1:r]
  Tr <- B.svd$u[, 1:r]
  
  A.star <- Diagonal( x = 1 / sqrt(Cr.diag) ) %*% t(Tr) %*% 
    A %*% 
    Tr %*% Diagonal( x = 1 / sqrt(Cr.diag) )
  
  A.star.svd <- svd( A.star )
  
  W <- Tr %*% Diagonal( x = 1 / sqrt( Cr.diag ) ) %*% A.star.svd$u
  
  return(W)
}

### Extracting A and B for minimally confounded level-1 residuals
mcresiduals <- 
  mclapply(sim.models,
         function(.mod){
           X <- getME(.mod, "X")
           Z <- BlockZ(.mod)
           
           n <- nrow(X)
           p <- ncol(X)
           ngrps <- unname( summary(.mod)@ngrps )
           
           vc <- VarCorr(.mod)
		   D <- kronecker( Diagonal(ngrps), as.matrix(vc[[1]]))
		   R <- Diagonal( n = n, x = attr(vc, "sc")^2 )
			
		   V <- Z %*% D %*% t(Z) + R
		   V.chol <- chol( V )
		   Vinv  <- chol2inv( V.chol )
			
		   xvix <- t(X) %*% Vinv %*% X
		   xvi <- t(X) %*% Vinv
		   rhs <- crossprod(xvi, solve(xvix, xvi))
		   P <- Vinv - rhs
			           
           # For level-1 min. confounded resids
           A1 <- R %*% P %*% Z %*% D %*% t(Z) %*% P %*% R
		   B1 <- R %*% P %*% R

           W1 <- try( mcrotate(A1, B1) )
           if( class(W1) == "try-error")  W1 <- NA else W1 <- as.matrix(W1)
           
           # For level-2 min. confounded resids for intercept
           L.int <- kronecker(Diagonal(ngrps), c(1, 0))
           pzdl.int <- P %*% Z %*% D %*% L.int
             
           A2.int <- t(L.int) %*% D %*% t(Z) %*% P %*% R %*% P %*% Z %*% D %*% L.int
           B2.int <- t(L.int) %*% D %*% t(Z) %*% P %*% Z %*% D %*% L.int
           W2.int <- try( mcrotate(A2.int, B2.int) )
           if( class(W2.int) == "try-error")  W2.int <- NA else W2.int <- as.matrix(W2.int)
           
           # For level-2 min. confounded resids for slope
           L.slope <- kronecker(Diagonal(ngrps), c(0, 1))
           
           A2.slope <- as.matrix( attr(vc, "sc")^2 * crossprod(P %*% Z %*% D %*% L.slope) )
           B2.slope <- t(L.slope) %*% D %*% t(Z) %*% P %*% Z %*% D %*% L.slope
           W2.slope <- try( mcrotate(A2.slope, B2.slope) )
           if( class(W2.slope) == "try-error")  W2.slope <- NA else W2.slope <- as.matrix(W2.slope)
           
           return(list( mclev1 = t( W1 ) %*% resid(.mod), 
                        mcrint = t( W2.int ) %*% ranef(.mod)[[1]][,1],
                        mcrslope = t( W2.slope ) %*% ranef(.mod)[[1]][,2]
                        ))
         })


### Calculating tests of normality from the simulations
norTests <- function(.list) {
  lapply(.list, 
         function(r) {
           sw  <- try(shapiro.test(r))
           ad  <- try(ad.test(r))
           cvm <- try(cvm.test(r))
           ks  <- try(lillie.test(r))
           
           pv <- data.frame(SW = ifelse(class(sw) == "try-error", NA, sw$p.value),
                            AD = ifelse(class(ad) == "try-error", NA, ad$p.value), 
                            CVM = ifelse(class(cvm) == "try-error", NA, cvm$p.value), 
                            KS = ifelse(class(ks) == "try-error", NA, ks$p.value))
           
           return(pv)  
         })
}

### Summarizing tests
summarizeTests <- function(pvals) {
  apply(pvals, 2, 
        FUN = function(pvs) {sum(pvs < .05, na.rm = TRUE) / sum(!is.na(pvs))})
}


### Calculate level-2 minimally confounded residuals for intercept
mclev2.int <- lapply(mcresiduals, function(x) x$mcrint)
test.mc.rint <- do.call("rbind", norTests(mclev2.int))
round(summarizeTests(test.mc.rint), 2)

saveRDS(mclev2.int, "radon-mcrint-13Dec12.RDS")

### Calculate level-1 minimally confounded residuals
mclev2.slope <- lapply(mcresiduals, function(x) x$mcrslope)
test.mc.rslope <- do.call("rbind", norTests(mclev2.slope))
round(summarizeTests(test.mc.rslope), 2)

saveRDS(mclev2.slope, "radon-mcrslope-13Dec12.RDS")

### Calculate level-1 minimally confounded residuals
mclev1 <- lapply(mcresiduals, function(x) x$mclev1)
test.mc.lev1 <- do.call("rbind", norTests(mclev1))
round(summarizeTests(test.mc.lev1), 2)

saveRDS(mclev1, "radon-mclev1-13Dec12.RDS")
