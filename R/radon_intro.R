#-------------------------------------------------------------------------------
# Script the graphics introducing the radon data
#
# Adam Loy
# April 2013
#-------------------------------------------------------------------------------

setwd("~/Documents/Thesis/Dissertation/eresids-chapter/minconfounded_chapter")

radon.data <- read.csv("data/radon_for_sims.csv")

srrs2 <- read.table ("data/srrs2.dat", header=T, sep=",")

## Restricting attention to Minnesota
mn <- subset(srrs2, state=="MN")

## clean names of counties: lower case and remove all white space at the end of names
mn$county <- tolower(as.character(mn$county))
mn$county <- gsub(" +$", "", as.character(mn$county))

library(ggplot2)
library(maps)
counties <- map_data("county")

# restrict to minnesota
counties <- subset(counties, region=="minnesota")
counties$region <- counties$subregion # no point in keeping the regions
qplot(data=counties, x=long, y=lat, group=group, geom="polygon")

# county level summary of radon data
library(plyr)
mn <- cbind(mn, uranium = radon.data$uranium)
counties.stat <- ddply(mn, .(county), summarize, 
                     radon=mean(activity, na.rm=T),
                    sdradon=sd(activity, na.rm=T),
                    uranium=unique(uranium))

ggplot(counties.stat, aes(map_id = county)) + geom_map(aes(fill=uranium), map = counties) + 
  expand_limits(x = counties$long, y = counties$lat) + xlab("longitude") + ylab("latitude") + coord_map() + theme_bw()

# county level sample size
ggplot(counties.stat, aes(map_id = county)) + geom_map(aes(fill=radon), map = counties) + 
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

ggsave(filename=file.choose(new=T), width = 4, height = 4, units = "in")

twocounties <- subset(mn, county %in% c("hennepin", "winona"))
twocounties$basement <- twocounties$floor==0
tc.stat <- ddply(twocounties, .(county, basement), summarise, sd=sd(activity)/length(activity), activity=mean(activity))

qplot(factor(basement), activity, data=twocounties, facets=~county) + 
  geom_errorbar(aes(ymin=activity-1.96*sd, ymax=activity+1.96*sd), data=tc.stat, width=0.5, size=0.75, colour="steelblue")+
  geom_point(data=tc.stat, size=3.5, colour="steelblue") + theme_bw()+theme(aspect.ratio=1) + ylab("radon") + xlab("basement")