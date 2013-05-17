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


qplot(x = s, y = fc, data = fc.melted, geom = c("point"), group = `Var. structure`,  colour = `Var. structure`, facets = ~ variable,  linetype = `Var. structure`, size=I(2.5), shape=b_dsn) + 
xlim(30, 60) + #geom_smooth(se = FALSE, size=1.5) +
scale_color_brewer("Variance structure", palette="Set2", labels = c(expression(paste(sigma[epsilon]^2==4, ", ", sigma[b]^2==1)), expression(paste(sigma[epsilon]^2==1, ", ", sigma[b]^2==1)), expression(paste(sigma[epsilon]^2==1, ", ", sigma[b]^2==4)))) + 
scale_linetype_discrete("Variance structure", labels = c(expression(paste(sigma[epsilon]^2==4, ", ", sigma[b]^2==1)), expression(paste(sigma[epsilon]^2==1, ", ", sigma[b]^2==1)), expression(paste(sigma[epsilon]^2==1, ", ", sigma[b]^2==4)))) + 
xlab("s") + ylab("fraction of confounding") + theme_bw() + 
theme(legend.position="bottom", legend.key.width = unit(3, "line")) + geom_line(aes(group=interaction(e_dsn, b_dsn, `Var. structure`))) + scale_shape_manual("Random Effects", values=c(1, 17, 15), labels = c("Exponential", "Normal", expression(t[3])))

ggsave("fc_by_s.pdf", width = 6.5, height=3.75, units="in")
