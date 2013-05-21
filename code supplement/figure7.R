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

### SET WORKING DIRECTORY TO WHERE RESULTS OF calFC LIVE

### Initializing
REMLmodels <- c("exp_exp_REMLmodels.RDS", "exp_norm_REMLmodels.RDS",  "exp_t_REMLmodels.RDS", 
"norm_exp_REMLmodels.RDS", "norm_norm_REMLmodels.RDS", "norm_t_REMLmodels.RDS",
"t_exp_REMLmodels.RDS", "t_norm_REMLmodels.RDS", "t_t_REMLmodels.RDS")
dist.combos <- t(sapply(strsplit(REMLmodels, "_"), FUN = function(x) return(c(error = x[1], `ranef` = x[2]))))

var.settings <- c("sige2_sigb1", "sige1_sigb1", "sige1_sigb2")

s <- paste("s", c(seq(30, 55, by = 5), "rankB"), sep = "")

### Summarizing FC
fc.summary <- data.frame(s = rep(s, each = 9*3), 
                         sigma_e = rep(rep(c(2, 1, 1), each = 9), times = 7),
                         sigma_b = rep(rep(c(1, 1, 2), each = 9), times = 7),
                         e_dsn = rep(dist.combos[,1], length.out = 189),
                         b_dsn = rep(dist.combos[,2], length.out = 189),
                         fc.b0 = NA, fc.b0.vmx = NA, fc.b1 = NA, fc.b1.vmx = NA)
cnt <- 0
for(i in seq(s)) {
	for(j in seq(var.settings)){
		load(paste(s[i], "/minconf_level2_fc_", var.settings[j], ".RData", sep = ""))
		
		index <- seq( (cnt * 9) + 1, length = 9)
		
		fc.summary[index, "fc.b0"]     <- sapply(b0.fc, mean, na.rm = T)
		fc.summary[index, "fc.b0.vmx"] <- sapply(b0.fc.vmx, mean, na.rm = T)
		fc.summary[index, "fc.b1"]     <- sapply(b1.fc, mean, na.rm = T)
		fc.summary[index, "fc.b1.vmx"] <- sapply(b1.fc.vmx, mean, na.rm = T)
		
		cnt <- cnt + 1
		
		rm(b0.fc, b0.fc.vmx, b1.fc, b1.fc.vmx)
	}
}
# Notice that, as expected, the FC for the varimax version is (roughly) the same 
# as the original rotation. Because of this, we drop the varimax versions from 
# the data frame.
fc.summary <- subset(fc.summary, select = -c(fc.b0.vmx, fc.b1.vmx))


### Plotting the FC
# converting s column to numeric
library(stringr)
fc.summary$s_int <- as.character(fc.summary$s)
fc.summary$s_slope <- as.character(fc.summary$s)

fc.summary$s_int[fc.summary$s_int == "srankB"] <- "s58"
fc.summary$s_slope[fc.summary$s_slope == "srankB"] <- "s59"

fc.summary$s_int <- as.numeric(str_extract(fc.summary$s_int, "\\d{1,2}"))
fc.summary$s_slope <- as.numeric(str_extract(fc.summary$s_slope, "\\d{1,2}"))

# dropping original s
fc.summary <- subset(fc.summary, select = -s) 

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

### Figure 7
qplot(x = s, y = fc, data = fc.melted, geom = c("point"), group = `Var. structure`,  colour = `Var. structure`, facets = ~ variable,  linetype = `Var. structure`, size=I(2.5), shape=b_dsn) + 
xlim(30, 60) + #geom_smooth(se = FALSE, size=1.5) +
scale_color_brewer("Variance structure", palette="Set2", labels = c(expression(paste(sigma[epsilon]^2==4, ", ", sigma[b]^2==1)), expression(paste(sigma[epsilon]^2==1, ", ", sigma[b]^2==1)), expression(paste(sigma[epsilon]^2==1, ", ", sigma[b]^2==4)))) + 
scale_linetype_discrete("Variance structure", labels = c(expression(paste(sigma[epsilon]^2==4, ", ", sigma[b]^2==1)), expression(paste(sigma[epsilon]^2==1, ", ", sigma[b]^2==1)), expression(paste(sigma[epsilon]^2==1, ", ", sigma[b]^2==4)))) + 
xlab("s") + ylab("fraction of confounding") + theme_bw() + 
theme(legend.position="bottom", legend.key.width = unit(3, "line")) + geom_line(aes(group=interaction(e_dsn, b_dsn, `Var. structure`))) + scale_shape_manual("Random Effects", values=c(1, 17, 15), labels = c("Exponential", "Normal", expression(t[3])))