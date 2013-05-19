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
library(grid)

### Reading in the data
propreject <- read.csv(file.choose())   # choose PropRejectSummary.csv
load(file.choose()) # choose GoldStandardPower.RData

goldstandard <- ldply(results.05, identity)
goldstandard$n <- rep(30:59, each = 3) 
goldstandard$ranef <- factor(goldstandard$.id)
levels(goldstandard$ranef)[c(1,3)] <- c("exp", "t")

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
qplot(x = s, y = prop.reject, data = subset(propreject, subset = ranef == "norm" & nortest == "AD" &  alpha == 0.05 & random.effect == "b1" & var.settings == "sige2_sigb1" & rotation != "EBLUP"), geom = c("point", "smooth"), colour = error, group = error, linetype = error, se = F, method = "lm") + 
	facet_grid(. ~ rotation, scales = "free") + 
	xlab("s") +
	ylab("proportion of tests rejected") + 
	scale_linetype_discrete("Error distribution", labels = c("Skewed", "Normal", "Heavy tailed")) + 	scale_color_discrete("Error distribution", labels = c("Skewed", "Normal", "Heavy tailed")) + 
	ylim(0, .4)

### Plotting the power
# Summarizing results of AD tests for random intercept
qplot(x = s, y = prop.reject, data = subset(propreject, subset = ranef != "norm" & random.effect == "b0" & nortest == "AD" & var.settings == "sige2_sigb1" & alpha == 0.05 & rotation != "EBLUP"), geom = c("point", "smooth"), facets = ~ rotation, linetype = ranef, colour = error, group = interaction(error, ranef), se = F, method = "lm") + 
	xlab("s") + 
	ylab("proportion of tests rejected") + 
	theme_bw() + 
	scale_color_discrete("Error distribution", labels = c("Skewed", "Normal", "Heavy tailed")) + 
	scale_linetype_discrete("Random effects\ndistribution", labels = c("Skewed", "Heavy tailed"))

# Summarizing results of AD tests for random slopes
qplot(x = s, y = prop.reject, data = subset(propreject, subset = ranef != "norm" & random.effect == "b1" & nortest == "AD" & var.settings == "sige2_sigb1" & alpha == 0.05 & rotation != "EBLUP"), geom = c("point", "smooth"), facets = ~ rotation, colour = error, linetype = ranef, group = interaction(error, ranef), se = F, method = "lm") + 
	xlab("s") + 
	ylab("proportion of tests rejected") + 
	theme_bw() + 
	scale_color_discrete("Error distribution", labels = c("Skewed", "Normal", "Heavy tailed")) + 
	scale_linetype_discrete("Random effects\ndistribution", labels = c("Skewed", "Heavy tailed"))


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
qplot(x = s, y = prop.reject, data = subset(propreject, subset = ranef != "norm" & random.effect == "b1" & nortest == "AD" & var.settings == "sige2_sigb1" & alpha == 0.05 & rotation != "EBLUP"), geom = "point", facets = ~ rotation, colour = error, linetype = ranef, shape = ranef, group = interaction(error, ranef), size = I(2.5)) + 
	geom_smooth(se = F, method = "lm",  size = I(.75)) + 
	scale_color_brewer( "Error distribution", palette = "Set2", labels = c("Exponential", "Normal", expression(t[3]))) + 
	geom_line(aes(x = n, y = AD, group = .id, linetype = ranef), data = subset(goldstandard, .id != "normal"), inherit.aes = FALSE, size = I(.75)) + 
#	geom_line(aes(x = n, y = AD, group = .id), data = subset(goldstandard, .id != "normal"), inherit.aes = FALSE, colour = "#FFD92F") + 
	xlab("s") + 
	ylab("proportion of tests rejected") + 
	theme_bw() +  
	scale_linetype_discrete("Random effects\ndistribution", labels = c("Exponential",  expression(t[3]))) + 
	scale_shape_discrete("Random effects\ndistribution", labels = c("Exponential",  expression(t[3]))) + 
	theme(legend.position = "bottom", legend.key.width = unit(3, "line"))



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

