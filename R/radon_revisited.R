#-------------------------------------------------------------------------------
# Script the graphics revisiting the radon data
#
# Adam Loy
# April 2013
#-------------------------------------------------------------------------------

### Preliminaries
library(lme4)
library(MASS)
library(RcppEigen)
library(inline)

setwd("~/Documents/Thesis/Dissertation/eresids-chapter/minconfounded_chapter")

# Readingin the data
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


### Reducing the fraction of confounding
reduced.fc <- data.frame(s = seq(30, 80, by = 10), b0 = NA, b1 = NA)

for(i in seq( nrow(reduced.fc) )){
	reduced.fc[i, "b0"] <- FC2(.mod = fm, .L = L.b0, s = reduced.fc[i, "s"])
	reduced.fc[i, "b1"] <- FC2(.mod = fm, .L = L.b1, s = reduced.fc[i, "s"])
}

### Obtaining the rotated random effects for Q-Q plots
b0.rot <- mcresid2(.mod = fm, .L = L.b0, s = 30)
b1.rot <- mcresid2(.mod = fm, .L = L.b1, s = 30)

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
ggsave(filename = "/Users/adam/Documents/Thesis/Dissertation/eresids-chapter/minconfounded_chapter/figures/rotatedQQ-intercept.pdf", width = 3, height = 3)


ggplot(data.frame(b1.rot), aes(sample = b1.rot)) + 
	geom_abline(intercept = HLMdiag:::qqlineInfo(b1.rot)[1], 
              slope = HLMdiag:::qqlineInfo(b1.rot)[2], colour = I("grey60")) + 
	stat_qq() + 
	geom_ribbon(aes(x = b1.qq$x, ymin = env.b1$lower, ymax = env.b1$upper), alpha = .2) + 
	xlab("Normal Quantiles") + ylab("Sample Quantiles") + 
	theme_bw()
ggsave(filename = "/Users/adam/Documents/Thesis/Dissertation/eresids-chapter/minconfounded_chapter/figures/rotatedQQ-slope.pdf", width = 3, height = 3)
