#-------------------------------------------------------------------------------
# Script summarizes the tests of normality of rotated random effects calculated 
# from residuals simulated from a model fit to the radon data
#
# Adam Loy
# May 2013
#-------------------------------------------------------------------------------

### Preliminaries
library(ggplot2)
library(plyr)
library(nortest)
library(reshape2)
library(stringr)

source("normality_functions.R")
source("utility_functions.R")

setwd("residuals/") # where the simulation results live

#-------------------------------------------------------------------------------
# Code to summarize simlation results
#-------------------------------------------------------------------------------

### Initializing
REMLmodels <- c("exp_exp_REMLmodels.RDS", "exp_norm_REMLmodels.RDS",  "exp_t_REMLmodels.RDS", 
"norm_exp_REMLmodels.RDS", "norm_norm_REMLmodels.RDS", "norm_t_REMLmodels.RDS",
"t_exp_REMLmodels.RDS", "t_norm_REMLmodels.RDS", "t_t_REMLmodels.RDS")
dist.combos <- t(sapply(strsplit(REMLmodels, "_"), FUN = function(x) return(c(error = x[1], `ranef` = x[2]))))

var.settings <- c("sige2_sigb1", "sige1_sigb1", "sige1_sigb2")

s <- paste("s", c(seq(30, 55, by = 5), "rankB"), sep = "")

### Summarizing type I error and power
power.summary <- NULL

cnt <- 0
for(i in seq(s)) {
	for(j in seq(var.settings)){
		load(paste(s[i], "/minconf_level2_residuals_", var.settings[j], ".RData", sep = ""))
		
		b0.mcr     <- test.simulation.results(sims = b0.mcresid, settings = dist.combos)
		b0.mcr.vmx <- test.simulation.results(sims = b0.mcresid.vmx, settings = dist.combos)

		b1.mcr     <- test.simulation.results(sims = b1.mcresid, settings = dist.combos)
		b1.mcr.vmx <- test.simulation.results(sims = b1.mcresid.vmx, settings = dist.combos)

		b0 <- merge(x = b0.mcr, y = b0.mcr.vmx, by = c("error", "ranef", "alpha"))
		colnames(b0)[4:9] <- c("AD.b0", "CVM.b0", "KS.b0", "AD.b0.vmx", "CVM.b0.vmx", "KS.b0.vmx")
		b1 <- merge(x = b1.mcr, y = b1.mcr.vmx, by = c("error", "ranef", "alpha"))
		colnames(b1)[4:9] <- c("AD.b1", "CVM.b1", "KS.b1", "AD.b1.vmx", "CVM.b1.vmx", "KS.b1.vmx")
		
		res <- merge(x = b0, y = b1, by = c("error", "ranef", "alpha"))
		
		res$s <- s[i]
		res$var.settings <- var.settings[j]
		
		power.summary <- rbind(power.summary, res)
	}
}



# converting s column to numeric
power.summary$s_int <- as.character(power.summary$s)
power.summary$s_slope <- as.character(power.summary$s)

power.summary$s_int[power.summary$s_int == "srankB"] <- "s58"
power.summary$s_slope[power.summary$s_slope == "srankB"] <- "s59"

power.summary$s_int <- as.numeric(str_extract(power.summary$s_int, "\\d{1,2}"))
power.summary$s_slope <- as.numeric(str_extract(power.summary$s_slope, "\\d{1,2}"))

# dropping original s
power.summary <- subset(power.summary, select = -s) 

# creating columns for normality test, random effect, and rotation
power.melted <- melt(power.summary, id.vars = c("error", "ranef", "alpha", "var.settings", "s_int", "s_slope"), value.name = "prop.reject")
power.melted$nortest <- str_extract(power.melted$variable, "[A-Z]+")
power.melted$random.effect <- str_extract(power.melted$variable, "[b]\\d")
power.melted$rotation <- !is.na(str_extract(power.melted$variable, "[vmx]{1,3}"))

# creating factors where necessary
power.melted$error <- factor(power.melted$error)
power.melted$ranef <- factor(power.melted$ranef)
power.melted$nortest <- factor(power.melted$nortest)
power.melted$random.effect <- factor(power.melted$random.effect)
power.melted$rotation <- factor(power.melted$rotation, labels = c("Rotated", "Varimax Rotated"))
power.melted$s <- NA
power.melted$s[power.melted$random.effect == "b0"] <- power.melted$s_int[power.melted$random.effect == "b0"]
power.melted$s[power.melted$random.effect == "b1"] <- power.melted$s_slope[power.melted$random.effect == "b1"]

power.melted <- subset(power.melted, select = -c(s_int, s_slope))

# now the same for the eblups
power.melted.original <-  melt(power.summary.original, id.vars = c("error", "ranef", "alpha", "var.settings", "s"), value.name = "prop.reject")
power.melted.original$nortest <- str_extract(power.melted.original$variable, "[A-Z]+")
power.melted.original$random.effect <- str_extract(power.melted.original$variable, "[b]\\d")
power.melted.original$rotation <- "EBLUP"

power.df <- rbind(power.melted, power.melted.original)

#-------------------------------------------------------------------------------
# Code for figures
#-------------------------------------------------------------------------------

propreject <- power.df

### Figure 8
ggplot(data = subset(propreject, subset = ranef == "norm" & nortest == "AD" &  alpha == 0.05 & random.effect == "b1" & var.settings == "sige2_sigb1" & rotation != "EBLUP"),
aes(x=s, y=prop.reject, colour=error, group=error, linetype=error, shape=error)) + 
geom_point(size=3) + geom_smooth(size=1, se=FALSE, method="lm") + 
	facet_grid(. ~ rotation, scales = "free") + 
	xlab("subspace dimension s") +
	ylab("proportion of tests rejected") + theme_bw() +
	scale_shape_manual("Error distribution", labels = c("Skewed", "Normal", "Heavy tailed"), values=c(1, 17, 15)) +
	scale_linetype_discrete("Error distribution", labels = c("Skewed", "Normal", "Heavy tailed")) + 	scale_color_brewer("Error distribution", labels = c("Skewed", "Normal", "Heavy tailed"), palette="Set2") + theme(legend.position="bottom", legend.key.width = unit(3, "line"))  + 
	ylim(0, .25) 


### Figure 9
ggplot(data = subset(propreject, subset = ranef != "norm" & random.effect == "b1" & nortest == "AD" & var.settings == "sige2_sigb1" & alpha == 0.05 & rotation != "EBLUP"),
aes(x=s,y = prop.reject, colour=error, linetype=ranef, shape=error, group= interaction(error, ranef))) +
geom_point(size=3) +geom_smooth(method="lm", size=1, se=F) +
	facet_wrap(~rotation) +
	xlab("subspace dimension s") + 
	ylab("proportion of tests rejected") + 
	theme_bw() + 
	scale_color_brewer("Error distribution", labels = c("Skewed", "Normal", "Heavy tailed"), palette="Set2") + 
	scale_linetype_discrete("Random effects\ndistribution", labels = c("Skewed", "Heavy tailed")) + theme_bw() +
	scale_shape_manual("Error distribution", labels = c("Skewed", "Normal", "Heavy tailed"), values=c(1, 17, 15)) + theme(legend.position="bottom", legend.key.width = unit(3, "line"))   + ylim(c(0,0.4))