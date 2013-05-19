### Simulation for plot to select the s-dimensional space

### Preliminaries
library(lme4)
library(ggplot2)
library(inline)
library(RcppEigen)
library(MASS)

source('~/Documents/Thesis/Dissertation/eresids-chapter/simulations/functions/cpp_functions.R')
source('~/Documents/Thesis/Dissertation/eresids-chapter/simulations/functions/resid_functions.R')
source('~/Documents/Thesis/Dissertation/eresids-chapter/simulations/functions/utility_functions.R')



### Simulating data
m  <- 50       ## no. of groups
ni <- 30       ## no. obs. per group
N  <- ni * m   ## total no. of obs.

sig.e <- 3
sig.b <- 1

b <- rnorm(m, 0, sig.b)
e <- rnorm(N, 0, sig.e)

grp.var <- rnorm(m)

sim.df <- data.frame(group = rep(1:m, times = ni), grp.var = rep(grp.var, times = ni), b = rep(b, times = ni), e)

# Introducing some small groups
n.small.grps <- 10
which.groups <- sample.int(m, size = n.small.grps)
size <- 5

for(i in which.groups) {
	index <- which(sim.df$group == i)
	#eliminate <- sample(index, size = length(index) - size)
	eliminate <- sample(index, size = length(index) - rpois(1, size))
	sim.df <- sim.df[-eliminate,]
}

# Introducing variability
n.var.grps <- 5
which.groups2 <- sample.int(m, size = n.var.grps)
mult.factor <- 3

for(i in which.groups2) {
	index <- which(sim.df$group == i)
	sim.df[index, "e"] <- rnorm(length(index), 0, mult.factor * sig.e)	
}

sim.df$y <- with(sim.df, 1 + .3 * grp.var + b + e)

mod <- lmer(y ~ grp.var + (1 | group), data = sim.df)

### plotting
### Reducing the trace
reduced.tr <- data.frame(s = seq(1, m, by = 1), b0 = NA)

for(i in seq( nrow(reduced.tr) )){
	reduced.tr[i, "b0"] <- tr2(.mod = mod, .L = Diagonal(m), s = reduced.tr[i, "s"])
}

qplot(x = m-s, y = b0/m, data = reduced.tr[-m,], geom = c("point", "line")) + 
	xlab("dimension reduction") + ylab("fraction of confounding") + theme_bw()


### Simulating an example -- a less obvious elbow
m  <- 85       ## no. of groups
ni <- c(rpois(60, 30), rpois(25, 5) + 1)        ## no. obs. per group
N  <- sum(ni)   ## total no. of obs.

sig.e <- 3
sig.b <- 1

b <- rnorm(m, 0, sig.b)
e <- rnorm(N, 0, sig.e)

#e[1:40] <- rnorm(40, 0, 3*sig.e)

grp.var <- rnorm(m)

sim.df <- data.frame(group = rep(1:m, times = ni), grp.var = rep(grp.var, times = ni), b = rep(b, times = ni), e)

sim.df$y <- with(sim.df, 1 + .3 * grp.var + b + e)

mod <- lmer(y ~ grp.var + (1 | group), data = sim.df)

### plotting
### Reducing the trace
reduced.tr <- data.frame(s = seq(1, m, by = 1), b0 = NA)

for(i in seq( nrow(reduced.tr) )){
	reduced.tr[i, "b0"] <- tr2(.mod = mod, .L = Diagonal(m), s = reduced.tr[i, "s"])
}

qplot(x = m - s, y = b0 / m , data = reduced.tr, geom = c("point", "line"))
