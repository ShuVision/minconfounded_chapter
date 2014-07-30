#-------------------------------------------------------------------------------
# This script summarizes FC for minimally confounded residuals for models fit to 
# a subset of the radon data.
#
# Adam Loy
# May 2013
#-------------------------------------------------------------------------------

### Preliminaries
library(ggplot2)
library(grid)
library(reshape2)
library(plyr)

### Reading in the data
fc.summary <- read.csv(file.choose()) # read in FCsummary.csv

# correcting fc
fc.summary$fc.b0 <- fc.summary$fc.b0 * fc.summary$s_int
fc.summary$fc.b1 <- fc.summary$fc.b1 * fc.summary$s_slope

fc.summary$fc.b0 <- fc.summary$fc.b0 / 60
fc.summary$fc.b1 <- fc.summary$fc.b1 / 60

### formatting df for ggplot
fc.melted <- melt(fc.summary, id.vars=c("sigma_e", "sigma_b", "e_dsn", "b_dsn", "s_int", "s_slope"), value.name = "fc")

# Some labels for the plots
fc.melted$`Var. structure` <- with(fc.melted, interaction(sigma_e, sigma_b), drop = T)
fc.melted$`Var. structure` <- factor(fc.melted$`Var. structure`, levels = levels(fc.melted$`Var. structure`)[c(2,1,3,4)])
levels(fc.melted$variable) <- c("intercept", "slope")

# Making an 's' column compatible with facetting
fc.melted$s <- NA
fc.melted$s[fc.melted$variable == "intercept"] <- fc.melted$s_int[fc.melted$variable == "intercept"]
fc.melted$s[fc.melted$variable == "slope"] <- fc.melted$s_slope[fc.melted$variable == "slope"]
names(fc.melted)[9] <- "var.struct"

ggplot(aes(x = s, y = fc, colour=var.struct), data=fc.melted) +
facet_grid(.~ variable) + xlim(30,60) +
xlab("s") + ylab("fraction of confounding") + 
geom_line(aes(group=interaction(e_dsn, b_dsn, var.struct), linetype=b_dsn)) + 
scale_linetype_discrete("Random Effects", labels = c("Exponential", "Normal", expression(t[3]))) + 
geom_point(aes(shape=b_dsn), size=2.5) + 
scale_shape_manual("Random Effects", values=c(1, 17, 15), labels = c("Exponential", "Normal", expression(t[3]))) +
scale_colour_brewer("Variance structure", palette="Set2", guide=FALSE) +
theme_bw() + 
theme(legend.position="bottom", legend.key.width = unit(3, "line")) +
geom_text(aes(x=x, y=y, label=label), data=data.frame(
x=c(52,52,52,
54,54,54), 
y=c(0.42,0.2,0,
0.48,0.28,0.10), 
var.struct=c("2.1", "1.1", "1.2",
			 "2.1", "1.1", "1.2"), 
variable=c("intercept","intercept", "intercept",
"slope","slope","slope"), 
label=c("list(sigma[epsilon]^2==4,sigma[b]^2==1)",
		"list(sigma[epsilon]^2==1,sigma[b]^2==1)",
		"list(sigma[epsilon]^2==1,sigma[b]^2==4)",
		"list(sigma[epsilon]^2==4,sigma[b]^2==1)",
		"list(sigma[epsilon]^2==1,sigma[b]^2==1)",
		"list(sigma[epsilon]^2==1,sigma[b]^2==4)")), 
parse=TRUE)

ggsave("fc_by_s.pdf", width = 6.5, height=3.75, units="in")
