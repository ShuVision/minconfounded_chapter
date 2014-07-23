#-------------------------------------------------------------------------------
# Script the graphics revisiting the radon data
#
# Adam Loy
# April 2013
#-------------------------------------------------------------------------------

### Preliminaries
library(lme4.0)
library(MASS)
library(RcppEigen)
library(inline)
library(nortest)

# Readingin the data
radon <- read.csv("original_radon.csv")

### Loading self-written functions
source('cpp_functions.R')
source('utility_functions.R')


tr2 <- function(.mod, .L, s = NULL){
	       y <- getME(.mod, "y")
           X <- getME(.mod, "X")
           Z <- BlockZ(.mod)
           
           n <- nrow(X)
           p <- ncol(X)
           ngrps <- unname( summary(.mod)@ngrps )
           
           vc <- VarCorr(.mod)
           Di <- bdiag( VarCorr(.mod) ) / (unname(attr(vc, "sc")))^2
           D  <- kronecker( Diagonal(ngrps), Di )
           
           Aslot <- getME(.mod, "A") # ZDZ'
           zdzt <- crossprod( getME(.mod, "A") )
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



### Fitted model
fm <- lmer(log.radon ~ basement + uranium + (basement | county), data = radon)

ngrps <- summary(fm)$ngrps
L.b0  <- kronecker(Diagonal(ngrps), c(1, 0))
L.b1  <- kronecker(Diagonal(ngrps), c(0, 1))


### Initial fraction of confounding
tr2(fm, L.b0) / 85 # 0.72
tr2(fm, L.b1) / 85  # 0.70


### Reducing the trace
reduced.tr <- data.frame(s = seq(30, 85, by = 1), b0 = NA, b1 = NA)

for(i in seq( nrow(reduced.tr) )){
	reduced.tr[i, "b0"] <- tr2(.mod = fm, .L = L.b0, s = reduced.tr[i, "s"])
	reduced.tr[i, "b1"] <- tr2(.mod = fm, .L = L.b1, s = reduced.tr[i, "s"])
}

# adding results from the AD test of the rotated residuals
reduced.tr.melt <- melt(reduced.tr, id.vars=1, variable.name="ranef")
reduced.tr.melt$p.value <- NA
for(i in 1:nrow(reduced.tr.melt)) {
	if(reduced.tr.melt[i,"ranef"] == "b0"){
		rot <- rotate_ranef(.mod = fm, .L = L.b0, s = reduced.tr.melt[i,"s"], .varimax=TRUE)
	} else{
		rot <- rotate_ranef(.mod = fm, .L = L.b1, s = reduced.tr.melt[i,"s"], .varimax=TRUE)
	}
	
	ad.result <- ad.test(rot)
	reduced.tr.melt$p.value[i] <- ad.result$p.value
}

reduced.tr.melt$ad.reject <- reduced.tr.melt$p.value < .05

### Figure 10
levels(reduced.tr.melt$ranef) <- c("Random intercept", "Random slope")
qplot(x = 85 - s, y = value / 85, data = reduced.tr.melt, geom = "line", facets = ~ ranef) + 
	geom_point(aes(shape = ad.reject, colour = ad.reject)) + 
	scale_colour_manual("AD test", labels = c("Fail to reject", "Reject"), values = c("black", "red")) + 
	scale_shape_manual("AD test", labels = c("Fail to reject", "Reject"), values = c(1, 16)) + 
	xlab("dimension reduction (s - 85)") + 
	ylab("fraction of confounding") +
	theme_bw() + 
	theme(legend.position = "bottom")


### Obtaining the rotated random effects for Q-Q plots
b0.rot <- rotate_ranef(.mod = fm, .L = L.b0, s = 65, .varimax = TRUE)
b1.rot <- rotate_ranef(.mod = fm, .L = L.b1, s = 65, .varimax = TRUE)

### constructing Q-Q plots
sim_env <- function(x, conf = .95){
  n <- length(x)
  P <- ppoints(x)
  z <- qnorm(P)
  a <- as.numeric(HLMdiag:::qqlineInfo(x)[1])
  b <- as.numeric(HLMdiag:::qqlineInfo(x)[2])
  zz <- qnorm(1 - (1 - conf)/2)
  SE <- (b/dnorm(z)) * sqrt(P * (1 - P)/n)
  fit.value <- a + b * z
  upper <- fit.value + zz * SE
  lower <- fit.value - zz * SE
  return(data.frame(lower, upper))
}


env.b0 <- sim_env(b0.rot)
b0.qq <- qqnorm(b0.rot, plot.it = FALSE)
b0.qq <- lapply(b0.qq, sort)

env.b1 <- sim_env(b1.rot)
b1.qq <- qqnorm(b1.rot, plot.it = FALSE)
b1.qq <- lapply(b1.qq, sort)


ggplot(data.frame(b0.rot), aes(sample = b0.rot)) + 
	geom_abline(intercept = HLMdiag:::qqlineInfo(b0.rot)[1], 
              slope = HLMdiag:::qqlineInfo(b0.rot)[2], colour = I("grey60")) + 
	stat_qq() + 
	geom_ribbon(aes(x = b0.qq$x, ymin = env.b0$lower, ymax = env.b0$upper), alpha = .2) + 
	xlab("Normal Quantiles") + ylab("Sample Quantiles") + 
	theme_bw()


ggplot(data.frame(b1.rot), aes(sample = b1.rot)) + 
	geom_abline(intercept = HLMdiag:::qqlineInfo(b1.rot)[1], 
              slope = HLMdiag:::qqlineInfo(b1.rot)[2], colour = I("grey60")) + 
	stat_qq() + 
	geom_ribbon(aes(x = b1.qq$x, ymin = env.b1$lower, ymax = env.b1$upper), alpha = .2) + 
	xlab("Normal Quantiles") + ylab("Sample Quantiles") + 
	theme_bw()