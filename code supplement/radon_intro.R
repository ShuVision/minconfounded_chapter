#-------------------------------------------------------------------------------
# Script the graphics introducing the radon data
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
### Preliminaries
#-------------------------------------------------------------------------------

## Load packages
library(ggplot2)   # for plotting
library(grid)      # for plotting
library(maps)      # for MN map
library(nullabor)  # for lineup
library(HLMdiag)   # for Q-Q plots
library(plyr)      # for data manipulation
library(reshape2)  # for data manipulation
library(stringr)   # for string manipulation

## Read in the data
radon.data <- read.csv("original_radon.csv")
srrs2 <- read.table ("srrs2.dat", header=T, sep=",")

## Restrict attention to Minnesota
mn <- subset(srrs2, state=="MN")

#-------------------------------------------------------------------------------
### Figures 1 and 2
#-------------------------------------------------------------------------------
## clean names of counties: 
## lower case and remove all white space at the end of names
mn$county <- tolower(as.character(mn$county))
mn$county <- gsub(" +$", "", as.character(mn$county))

counties <- map_data("county")

## restrict to minnesota
counties <- subset(counties, region=="minnesota")
counties$region <- counties$subregion # no point in keeping the regions

## county level summary of radon data
library(plyr)
mn <- cbind(mn, uranium = radon.data$uranium)
counties.stat <- ddply(mn, .(county), summarize, 
					 n=sum(complete.cases(activity)),
                     radon=mean(activity, na.rm=T),
                    sdradon=sd(activity, na.rm=T),
                    uranium=unique(uranium))

### Figure 1
ggplot(counties.stat, aes(map_id = county)) + geom_map(aes(fill=radon), map = counties) + 
  scale_fill_gradient("radon activity\n(log pCi/L)") + 
  expand_limits(x = counties$long, y = counties$lat) + 
  xlab("longitude") + ylab("latitude") + coord_map() + 
  theme(
  	axis.text.x = element_blank(), 
  	axis.text.y = element_blank(),
  	axis.title.x = element_blank(), 
  	axis.title.y = element_blank(),
  	axis.ticks = element_line(colour=rgb(0,0,0,alpha=0)),
  	panel.background =  element_blank(),
  	panel.grid.major =  element_blank(),
  	panel.grid.minor =  element_blank())


## Restrict attention to two counties
twocounties <- subset(mn, county %in% c("hennepin", "winona"))
twocounties$basement <- twocounties$floor==0
tc.stat <- ddply(twocounties, .(county, basement), summarise, n=length(activity), sd=sd(activity)/length(activity), activity=mean(activity))

## Figure 2
qplot(factor(basement), activity, data=twocounties, facets=~county) + 
  geom_errorbar(aes(ymin=activity-qt(p=.975, df = n-1, lower.tail = F)*sd, ymax=activity+qt(p=.975, df = n-1, lower.tail = F)*sd), data=tc.stat, width=0.5, size=0.75, colour="steelblue")+
  geom_point(data=tc.stat, size=3.5, colour="steelblue") + theme_bw()+theme(aspect.ratio=1) + ylab("radon") + xlab("basement")


#-------------------------------------------------------------------------------
### Figure 3 and 4
#-------------------------------------------------------------------------------
### Fitting a model with random intercept for county and 
### random slope for basement.
fm <- lmer(log.radon ~ basement + uranium + (basement | county), data = radon.data, REML = FALSE)

# Level-1 residuals (i.e., the error terms)
e <- resid(fm)

# Random effects (i.e., the level-2 residuals)
b <- ranef(fm)[[1]]

# Envelope for normal deviates
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

# Q-Q plot for the error terms
e.env <- sim_env(e)
e.qq <- qqnorm(e, plot.it = FALSE)
e.qq <- lapply(e.qq, sort)

ggplot(data.frame(e), aes(sample = e)) + 
	geom_abline(intercept = HLMdiag:::qqlineInfo(e)[1], 
              slope = HLMdiag:::qqlineInfo(e)[2], colour = I("grey60")) + 
	stat_qq() + 
	geom_ribbon(aes(x = e.qq$x, ymin = e.env$lower, ymax = e.env$upper), alpha = .2) + 
	xlab("Normal Quantiles") + ylab("Sample Quantiles") + 
	ylim(-4,4) + 
	theme_bw()

# Q-Q plot for the random intercept
b0.env <- sim_env(b[,1])
b0.qq <- qqnorm(b[,1], plot.it = FALSE)
b0.qq <- lapply(b0.qq, sort)

ggplot(data.frame(b[,1]), aes(sample = b[,1])) + 
	geom_abline(intercept = HLMdiag:::qqlineInfo(b[,1])[1], 
              slope = HLMdiag:::qqlineInfo(b[,1])[2], colour = I("grey60")) + 
	stat_qq() + 
	geom_ribbon(aes(x = b0.qq$x, ymin = b0.env$lower, ymax = b0.env$upper), alpha = .2) + 
	xlab("Normal Quantiles") + ylab("Sample Quantiles") + 
	ylim(-.2, .2) + 
	theme_bw()
	

#-------------------------------------------------------------------------------
### Figure 4
#-------------------------------------------------------------------------------
set.seed(987654321)
sim.y   <- simulate(fm, nsim = 19)                        ## A 919 x 19 matrix of responses
sim.mod <- apply(sim.y, 2, refit, object = fm)            ## a list of models

### We are really interested in a lineup for the random slopes for our paper
# Simulated random slopes
sim.b1 <- llply(sim.mod, function(x) ranef(x)[[1]][,2])   ## a list of random slopes
sim.b1 <- melt( do.call("rbind", sim.b1) )[,-2]           ## changing to a data frame
names(sim.b1) <- c("sample", "basement")                  ## setting colnames for faceting
sim.b1        <- arrange(sim.b1, sample)                  ## ordering by simulation

### Creating a the lineup
b1 <- transform(b, band = sim_env(basement), x = sort(qqnorm(basement, plot.it=FALSE)$x))

# some tweaks for nullabor
sim.b1$.n <- as.numeric( str_extract(sim.b1$sample, "\\d+") )
sim.b1 <- ddply(sim.b1, .(.n), transform, band = sim_env(basement), x = sort(qqnorm(basement, plot.it=FALSE)$x))

### Figure 4
set.seed(123456)
qplot(sample = basement, data = b1, stat = "qq") %+%
	lineup(true = b1, sample = sim.b1) + 
	facet_wrap(~ .sample, ncol = 5) + 
	geom_ribbon(aes(x = x, ymin = band.lower, ymax = band.upper), alpha = .25) + 
	xlab("Normal Quantiles") + ylab("Sample Quantiles") +  
	theme_bw() + 
	theme(panel.margin = unit(0, "lines"))
