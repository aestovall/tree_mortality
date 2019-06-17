#Start with trees_all DF
ALLtrees<-read.csv("output/ALLtrees_v2.csv")[,-1]
ALLtrees$z_bin<-5*floor(ALLtrees$zmax/5)

#We aggregate to calculate total mortality over the time period
mort<-aggregate(count~mort_year+z_bin,FUN = "sum", data = ALLtrees[ALLtrees$dead==1,],drop = FALSE)
mort$dead<-mort$count
mort$total<-merge.data.frame(mort,aggregate(count~z_bin,FUN = "sum", data = ALLtrees,drop = FALSE), by = "z_bin")$count.y

# The number of dead trees / total = the mortality over the observed period
mort$rt<-mort$dead/mort$total

#order the mortality years for the figure
mort$mort_year<-factor(mort$mort_year)
mort$mort_year<- factor(mort$mort_year, levels = rev(levels(mort$mort_year)))

#classify the trees into small, medium, and large
mort$size_class<-"size"
mort$size_class[mort$z_bin<=15]<-"small"
mort$size_class[mort$z_bin>15&mort$z_bin<30]<-"medium"
mort$size_class[mort$z_bin>=30]<-"large"
mort$size_class<-factor(mort$size_class)

#create a "count" column
mort$count<-1

#calculate how long the time between image aquisitions to correct the mortality rate
mort$tsd<-1
mort$tsd[mort$mort_year=="2016"]<-2016-2014
mort$tsd[mort$mort_year=="2014"]<-2014-2012
mort$tsd[mort$mort_year=="2012"]<-2012-2009
mort$tsd[mort$mort_year=="2010"]<-2010-2009
mort$tsd[mort$mort_year=="2009"]<-0

#mortality divided by time period = mortality rate (%/year)
mort$rt[is.na(mort$rt)]<-0
mort$rt_per<-mort$rt
mort$rt_per[mort$mort_year!="2009"] <- (mort$rt/mort$tsd)[mort$mort_year!="2009"]

#Adds the period of observation and makes this a factor four our figure
mort$year[mort$mort_year=="2016"]<-"2014-2016"
mort$year[mort$mort_year=="2014"]<-"2012-2014"
mort$year[mort$mort_year=="2012"]<-"2010-2012"
mort$year[mort$mort_year=="2010"]<-"2009-2010"
mort$year[mort$mort_year=="2009"]<-0


#create Figure 1A
p_rt_cum<-ggplot(mort[mort$z_bin>=5&mort$z_bin<=65&mort$mort_year!="2009",], 
                 aes(x = z_bin, y = rt_per*100*2, group = rev(as.factor(year)))) +
  theme_bw() +
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_area(aes(fill = year))+
  geom_area(aes(fill = year))+
  scale_x_continuous(limits = c(5,65), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,70), expand = c(0, 0)) +
  geom_vline(xintercept = 15, linetype = 2, size = 0.5, color = "black") +
  geom_vline(xintercept = 30, linetype = 2, size = 0.5, color = "black") +
  annotate("text", label = "Small", x = 10, y = 19, size = 4, colour = "black", angle = 90, hjust = 0)+
  annotate("text", label = "Medium", x = 22.5, y = 19, size = 4, colour = "black", angle = 90, hjust = 0)+
  annotate("text", label = "Large", x = 47.5, y = 19, size = 4, colour = "black", angle = 90, hjust = 0)+
  # xlim(5,60)+
  scale_fill_viridis(discrete = TRUE, direction = 1, name = "Mortality Period") +
  theme(legend.position = c(0.2, 0.77)) +
  labs(x=bquote('Tree Height (m)'), y=bquote('Cumulative Mortality (%)')) +
  theme(axis.title.x = element_text(vjust=-0.5)) +
  theme(axis.title.y = element_text(vjust=1.5)) +
  theme(plot.title = element_text(vjust=2.5))

# Now we aggregate our results to visualize the mortality rate over time
# for the three size classes of trees. 
death_rate<-aggregate(rt_per~mort_year+size_class,FUN = "mean",mort[mort$z_bin>=5&mort$z_bin<=65,])
death_rate$mort_year <- factor(death_rate$mort_year, levels = rev(levels(death_rate$mort_year)))
death_rate$error<-aggregate(rt_per~mort_year+size_class,FUN = "sd",mort[mort$z_bin>=5&mort$z_bin<=65,])[,3]
death_rate$n<-aggregate(count~mort_year+size_class,FUN = "sum",mort[mort$z_bin>=5&mort$z_bin<=65,])[,3]

#get the error bars -> SE*1.96 = CI
death_rate$sderr<-death_rate$error/death_rate$n

death_rate$error_min<-death_rate$rt_per - death_rate$sderr*1.96
death_rate$error_max<-death_rate$rt_per + death_rate$sderr*1.96
death_rate$error_min[death_rate$error_min<0]<-0

#Convert to Date, so our observations line up with the summer aquisitions
death_rate$year<-death_rate$mort_year
death_rate$mort_year<-as.Date(as.character(paste(death_rate$year,"0705",sep = "")),"%Y%m%d")


#Make Figure 1B
options(warn=-1) # suppress warnings from loess fit
p_rt<-ggplot(death_rate[death_rate$year!="2009",], aes(x = mort_year, y = rt_per*100, 
                                                       group = size_class,
                                                       color = size_class, fill = size_class)) +
  theme_bw() +
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
 geom_errorbar(aes(x = mort_year, ymin = error_min*100, ymax = error_max*100), width = 100) +
  stat_smooth(se = FALSE)+
  geom_point(aes(shape = size_class), size = 2) +
  ylim(0,22) +
  xlim(as.Date("20100101","%Y%m$d"),as.Date("20161231","%Y%m%d"))+
  theme(legend.position = c(0.2, 0.8)) +
  scale_color_manual(name = "Tree Size",
                     labels = c("Large", "Medium","Small"),
                     values = c("tomato2","royalblue3", "springgreen3")) +
  scale_fill_manual(name = "Tree Size",
                    labels = c("Large", "Medium","Small"),
                    values = c("tomato2","royalblue3", "springgreen3")) +
  scale_shape_manual(name = "Tree Size",
                     labels = c("Large", "Medium","Small"),
                     values = c(21,22,24)) +
  labs(x=bquote('Year'), y=bquote('Mortality Rate (% '*yr^{-1}*')')) +
  theme(axis.title.x = element_text(vjust=-0.5)) +
  theme(axis.title.y = element_text(vjust=1.5)) +
  theme(plot.title = element_text(vjust=2.5))

#Make the final plot
library(cowplot)
plot_grid(p_rt_cum,p_rt, labels = c("a", "b"))
options(warn=0)
ggsave("figures/Figure1.pdf", width = 3.75*2, height = 3.75, units = "in")

