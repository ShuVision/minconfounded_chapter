#-------------------------------------------------------------------------------
# Script summarizes the tests of normality of rotated random effects calculated 
# from residuals simulated from a model fit to the radon data discussed by 
# Gelman and Hill.
#
# Adam Loy
# May 2013
#-------------------------------------------------------------------------------

### Preliminaries
library(ggplot2)
library(plyr)

### Reading in the data
propreject <- read.csv(file.choose())   # choose PropRejectSummary.csv
load(file.choose()) # choose GoldStandardPower.RData

goldstandard <- ldply(results.05, identity)
goldstandard$n <- rep(30:59, each = 3) 

#-------------------------------------------------------------------------------
# Code for the current figures
#-------------------------------------------------------------------------------

# Summarizing results of AD tests for random intercept
qplot(x = s, y = prop.reject, data = subset(propreject, subset = ranef == "norm" & nortest == "AD" &  alpha == 0.05 & random.effect == "b0" & var.settings == "sige2_sigb1" & rotation != "EBLUP"), geom = c("point", "smooth"), colour = error, group = error, linetype = error, se = F, method = "lm") + 
	facet_grid(. ~ rotation, scales = "free") + 
	xlab("s") +
	ylab("proportion of tests rejected") + 
	scale_linetype_discrete("Error distribution", labels = c("Skewed", "Normal", "Heavy tailed")) + 	scale_color_discrete("Error distribution", labels = c("Skewed", "Normal", "Heavy tailed"))

# Summarizing results of AD tests for random slope
ggplot(data = subset(propreject, subset = ranef == "norm" & nortest == "AD" &  alpha == 0.05 & random.effect == "b1" & var.settings == "sige2_sigb1" & rotation != "EBLUP"),
aes(x=s, y=prop.reject, colour=error, group=error, linetype=error, shape=error)) + 
geom_point(size=3) + geom_smooth(size=1, se=FALSE, method="lm") + 
	facet_grid(. ~ rotation, scales = "free") + 
	xlab("subspace dimension s") +
	ylab("proportion of tests rejected") + theme_bw() +
	scale_shape_manual("Error distribution", labels = c("Skewed", "Normal", "Heavy tailed"), values=c(1, 17, 15)) +
	scale_linetype_discrete("Error distribution", labels = c("Skewed", "Normal", "Heavy tailed")) + 	scale_color_brewer("Error distribution", labels = c("Skewed", "Normal", "Heavy tailed"), palette="Set2") + theme(legend.position="bottom", legend.key.width = unit(3, "line"))  + 
	ylim(0, .25) 

ggsave("ad_slope_results.pdf", width=8, height=5)

### Plotting the power
# Summarizing results of AD tests for random intercept
qplot(x = s, y = prop.reject, data = subset(propreject, subset = ranef != "norm" & random.effect == "b0" & nortest == "AD" & var.settings == "sige2_sigb1" & alpha == 0.05 & rotation != "EBLUP"), geom = c("point", "smooth"), facets = ~ rotation, linetype = ranef, colour = error, group = interaction(error, ranef), se = F, method = "lm") + 
	xlab("s") + 
	ylab("proportion of tests rejected") + 
	theme_bw() + 
	scale_color_discrete("Error distribution", labels = c("Skewed", "Normal", "Heavy tailed")) + 
	scale_linetype_discrete("Random effects\ndistribution", labels = c("Skewed", "Heavy tailed"))

# Summarizing results of AD tests for random slopes
ggplot(data = subset(propreject, ranef != "norm" & random.effect == "b1" & 
								 nortest == "AD" & var.settings == "sige2_sigb1" & 
								 alpha == 0.05 & rotation != "EBLUP"),
		aes(x=s,y = prop.reject, colour=error)) +
	geom_point(aes(shape=error), size=3) +
	geom_smooth(aes(linetype=ranef,  group= interaction(error, ranef)), method="lm", size=1, se=F) +
	facet_wrap(~rotation) +
	xlab("subspace dimension s") + 
	ylab("proportion of tests rejected") + 
	theme_bw() +  ylim(c(0,0.4)) + 
	scale_color_brewer("Error distribution", labels = c("Skewed", "Normal", "Heavy tailed"), palette="Set2", guide=FALSE) + 
	scale_linetype_discrete("Random effects\ndistribution", labels = c("Skewed", "Heavy tailed")) + 
	theme_bw() +
	theme(legend.position="bottom", legend.key.width = unit(3, "line"))   +
	scale_shape_manual(values=c(1, 17, 15), guide=FALSE)  +
	geom_text(aes(x=x,y=y,label=label), size=4, data=
		data.frame(x=c(55,55,52.5,52.5,55,52.5), 
				   y=c(0.09,0.04,0.22, 0.22, 0.06, 0.32), 
				   label = c("Skewed", "Normal error", "Heavy tailed error",
				   			 "Skewed error", "Normal error", "Heavy tailed"), 
				   rotation=c("Rotated", "Rotated", "Rotated", 
				   			  "Varimax Rotated", "Varimax Rotated", "Varimax Rotated"),
				   error=c("exp","norm","t", "exp","norm","t")))
	
ggsave("ad_slope_power.pdf", width=8, height=5.25)
#-------------------------------------------------------------------------------
# Adding power for the AD test
#-------------------------------------------------------------------------------

# Summarizing results of AD tests for random intercept
qplot(x = s, y = prop.reject, data = subset(propreject, subset = ranef != "norm" & random.effect == "b0" & nortest == "AD" & var.settings == "sige2_sigb1" & alpha == 0.05 & rotation != "EBLUP"), geom = c("point", "smooth"), facets = ~ rotation, linetype = ranef, colour = error, group = interaction(error, ranef), se = F, method = "lm") + 
	geom_line(aes(x = n, y = AD, group = .id), data = subset(goldstandard, .id != "normal"), inherit.aes = FALSE) + 
	xlab("s") + 
	ylab("proportion of tests rejected") + 
	theme_bw() + 
	scale_color_discrete("Error distribution", labels = c("Skewed", "Normal", "Heavy tailed")) + 
	scale_linetype_discrete("Random effects\ndistribution", labels = c("Skewed", "Heavy tailed"))

# Summarizing results of AD tests for random slopes
qplot(x = s, y = prop.reject, data = subset(propreject, subset = ranef != "norm" & random.effect == "b1" & nortest == "AD" & var.settings == "sige2_sigb1" & alpha == 0.05 & rotation != "EBLUP"), geom = c( "smooth"), facets = ~ rotation, colour = error, linetype = ranef, group = interaction(error, ranef), se = F, method = "lm") + 
	geom_line(aes(x = n, y = AD, group = .id), data = subset(goldstandard, .id != "normal"), inherit.aes = FALSE) + 
	xlab("s") + 
	ylab("proportion of tests rejected") + 
	theme_bw() + 
	scale_color_discrete("Error distribution", labels = c("Skewed", "Normal", "Heavy tailed")) + 
	scale_linetype_discrete("Random effects\ndistribution", labels = c("Skewed", "Heavy tailed"))



#-------------------------------------------------------------------------------
# Trying to add in a third panel for the EBLUPs
#-------------------------------------------------------------------------------

propreject$rotation <- factor(propreject$rotation, levels = levels(propreject$rotation)[c(2,3,1)])
propreject$s <- factor(propreject$s)

### Plotting the type I error
# Summarizing results of AD tests for random intercept
qplot(x = s, y = prop.reject, data = subset(propreject, subset = ranef == "norm" & nortest == "AD" &  alpha == 0.05 & random.effect == "b0" & var.settings == "sige2_sigb1"), geom = c("point", "smooth"), colour = error, group = error, linetype = error, se = F, method = "lm") + 
	facet_grid(. ~ rotation, scales = "free") + 
	xlab("s") +
	ylab("proportion of tests rejected") + 
	scale_linetype_discrete("Error distribution", labels = c("Skewed", "Normal", "Heavy tailed")) + 	scale_color_discrete("Error distribution", labels = c("Skewed", "Normal", "Heavy tailed")) + 
	ylim(0, .2)

# Summarizing results of AD tests for random slope
qplot(x = s, y = prop.reject, data = subset(propreject, subset = ranef == "norm" & nortest == "AD" &  alpha == 0.05 & random.effect == "b1" & var.settings == "sige2_sigb1"), geom = c("point", "smooth"), colour = error, group = error, linetype = error, se = F, method = "lm") + 
	facet_grid(. ~ rotation, scales = "free") + 
	xlab("s") +
	ylab("proportion of tests rejected") + 
	scale_linetype_discrete("Error distribution", labels = c("Skewed", "Normal", "Heavy tailed")) + 	scale_color_discrete("Error distribution", labels = c("Skewed", "Normal", "Heavy tailed")) + 
	ylim(0, .4)

### Plotting the power
# Summarizing results of AD tests for random intercept
qplot(x = 2, y = prop.reject, data = subset(propreject, subset = ranef != "norm" & random.effect == "b0" & nortest == "AD" & var.settings == "sige2_sigb1" & alpha == 0.05), geom = c("point", "smooth"), facets = ~ rotation, colour = error, group = interaction(error, ranef), se = F, method = "lm") + xlab("s") + ylab("proportion of tests rejected") + theme_bw() + scale_color_discrete("Error distribution", labels = c("Skewed", "Normal", "Heavy tailed")) + scale_linetype_discrete("Random effects\ndistribution", labels = c("Skewed", "Heavy tailed"))

# Summarizing results of AD tests for random slopes
qplot(x = 2, y = prop.reject, data = subset(power.melted, subset = ranef != "norm" & random.effect == "b1" & nortest == "AD" & var.settings == "sige2_sigb1" & alpha == 0.05), geom = c("smooth"), facets = ~ varimax, colour = error, group = interaction(error, ranef), se = F, method = "lm", linetype = ranef) + xlab("s") + ylab("proportion of tests rejected") + theme_bw() + scale_color_discrete("Error distribution", labels = c("Skewed", "Normal", "Heavy tailed")) + scale_linetype_discrete("Random effects\ndistribution", labels = c("Skewed", "Heavy tailed"))

