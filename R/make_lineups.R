
std_lineup <- function(dframe) {
  require(ggplot2)
  print(ggplot(aes(x=naive1.qq.x, y=naive1.qq.y), data = dframe) + 
          geom_smooth(aes(naive1.qq.x, naive1.env.fit.value), colour="grey50", se=FALSE, method="loess") +
          facet_wrap(~.sample, ncol=5) +
          geom_point() + 
          geom_ribbon(aes(x = naive1.qq.x, ymin = naive1.env.lower, ymax = naive1.env.upper),alpha = .2)+
          labs(x = "", y = "") +
          theme_bw() + 
          xlab("Normal Quantiles") + 
          ylab("Sample Quantiles") +
          theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
                axis.ticks = element_blank())
  )	
}

rot_lineup <- function(dframe) {
  require(ggplot2)
  print(  ggplot(aes(x=naive1.qq.x, y=naive1.qq.y-naive1.env.fit.value), data = dframe) + 
            facet_wrap(~.sample, ncol=5) + 
            geom_hline(yintercept=0, colour="grey30")+
            geom_point() + 
            geom_ribbon(aes(x = naive1.qq.x, 
                            ymin = naive1.env.lower-naive1.env.fit.value, ymax = naive1.env.upper-naive1.env.fit.value),alpha = .2)+
            labs(x = "", y = "") +
            theme_bw() + 
            xlab("Normal Quantiles") + 
            ylab("Sample Quantiles") +
            theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
                  axis.ticks = element_blank()))
}

ctrl_lineup <- function(dframe) {
  require(ggplot2)
  print(ggplot(aes(x=naive1.qq.x, y=naive1.qq.y), data = dframe) + 
          geom_smooth(aes(naive1.qq.x, naive1.env.fit.value), colour="grey50", se=FALSE, method="loess") +
          geom_point() + 
          facet_wrap(~.sample, ncol=5) +
          labs(x = "", y = "") +
          theme_bw() + 
          theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
                axis.ticks = element_blank()) 
  )
}

  
lp <- read.csv(file.choose()) # read data/radondata-1-1-0-85-0-14-10.csv
 
idx <- c(1:16, 20, 17:19)
lp$.sample <- c(1:15,17:20, 16)[lp$.n]
std_lineup(lp)
ggsave(filename="qq-lineup-std.pdf")
rot_lineup(lp)
ggsave(filename="qq-lineup-rot.pdf")
 
   
  

