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
library(nortest)
library(ggplot2)

setwd("~/Documents/Thesis/Dissertation/eresids-chapter/minconfounded_chapter")

# Reading in the data
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

### Initial trace values
tr(fm, L.b0)  # 61.53067
tr(fm, L.b1)  # 59.5131

### Reducing the trace
reduced.tr <- data.frame(s = seq(1, 85, by = 1), b0 = NA, b1 = NA)

for(i in seq( nrow(reduced.tr) )){
	reduced.tr[i, "b0"] <- tr2(.mod = fm, .L = L.b0, s = reduced.tr[i, "s"])
	reduced.tr[i, "b1"] <- tr2(.mod = fm, .L = L.b1, s = reduced.tr[i, "s"])
}

qplot(x = 85 - s, y = b0 / 85, data = reduced.tr, geom = c("point", "line")) + 
	xlab("dimension reduction") + 
	ylab("fraction of confounding") +
	theme_bw()
	
qplot(x = 85 - s, y = b1 / 85, data = reduced.tr, geom = c("point", "line")) +
	xlab("dimension reduction") + 
	ylab("fraction of confounding") + 
	theme_bw()
	
qplot(x = 85 - s, y = b0 / tr(fm, L.b0), data = reduced.tr, geom = c("point", "line"))
qplot(x = 85 - s, y = b1 / tr(fm, L.b1), data = reduced.tr, geom = c("point", "line"))


### Reducing the fraction of confounding
reduced.fc <- reduced.tr
reduced.fc[,2:3] <- reduced.fc[,2:3] / 85

qplot(x = 85 - s, y = b0, data = reduced.fc, geom = c("point", "line")) + ylim(c(0,1))
qplot(x = 85 - s, y = b1, data = reduced.fc, geom = c("point", "line"))


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

b0.rot <- mcresid2(.mod = fm, .L = L.b0, s = 80)
b1.rot <- mcresid2(.mod = fm, .L = L.b1, s = 85)

b0.rot <- vector("list", length = 56)
b0.rot <- lapply(30:85, mcresid2, .mod = fm, .L = L.b0)
sapply(b0.rot, FUN = function(x) {
	j <- ad.test(x)
	return(j$p.value)
})

ad.test(b0.rot)
ad.test(b1.rot)

cvm.test(b0.rot)
cvm.test(b1.rot)

lillie.test(b0.rot)
lillie.test(b1.rot)