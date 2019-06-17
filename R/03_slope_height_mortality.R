library(ggplot2)
library(cowplot)
library(viridis)

ALLtrees<-read.csv("output/ALLtrees_v2.csv")[,-1]
colnames(ALLtrees)[2]<-"zmax"
# ALLtrees$z_bin<-5*floor(ALLtrees$max/5)

ALLtrees<-ALLtrees[ALLtrees$mort_year=="2016"|ALLtrees$mort_year=="Live",]

ALLtrees_bin<-data.frame(dead = ALLtrees$dead,
                         count = ALLtrees$count,
                         zmax = 5*floor(ALLtrees$zmax/5),
                         cover = 0.1*floor(ALLtrees$cover/0.1),
                         tmaxn = 1*floor(ALLtrees$tmaxn/1),
                         tmax = 1*floor(ALLtrees$tmax/1),
                         pptn = 10*floor(ALLtrees$pptn/10),
                         ppt = 10*floor(ALLtrees$ppt/10),
                         slope = 5*floor(ALLtrees$slope/5),
                         aws = 1*floor(ALLtrees$aws/1),
                         vpdmax = 1*floor(ALLtrees$vpdmax/1)
)

####Mortality height relationships######

stats_ls<-list()

### AWS    ####
mort<-aggregate(count ~ zmax + aws, FUN = "sum", ALLtrees_bin, drop = FALSE, na.action = NULL)
mort$dead2016<-aggregate(dead ~ zmax + aws, FUN = "sum", ALLtrees_bin, drop = FALSE, na.action = NULL)[,3]
mort$rt<-mort$dead2016/mort$count/2
mort$rt[is.nan(mort$rt)]<-NA
mort<-mort[mort$count>10,]

v_ls<-list()
range<-seq(round(min(ALLtrees_bin$aws,na.rm = TRUE),0),round(max(ALLtrees_bin$aws,na.rm = TRUE),0),by=1)
for(i in 1:length(range)){
  if(nrow(na.exclude(mort[mort$aws>range[i]&mort$aws<=range[i+1],]))==0) next
  fit_tree<-lm(rt ~ zmax, weights = count, data = na.exclude(mort[mort$aws>range[i]&mort$aws<=range[i+1],]))
  v_ls[[i]] <- data.frame(temp = range[i], zrt=coef(fit_tree)[2], ct = sum(mort$count[mort$aws>range[i]&mort$aws<=range[i+1]], na.rm = TRUE),
                          se = sigma(fit_tree), minci = confint(fit_tree)[2,1], maxci = confint(fit_tree)[2,2], p = summary(fit_tree)$coefficients[2,4])
}

v_all<-do.call(rbind,v_ls)
v_all$var<-expression(SOIL[AWC])
v_all<-v_all[v_all$p<0.05,]

ggplot(mort, aes(y = rt*100, x = zmax, group = factor(aws), color = factor(aws))) + 
  theme_bw() + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  stat_smooth(aes(group = factor(0.5*floor(aws/0.5)), weight=count), 
              method = "lm", se = FALSE,
              size = 0.5) +
  geom_jitter(aes(size = count), width = 5, alpha = 1, shape = 21) + 
  scale_color_viridis(name = "Available Water Capacity",
                      discrete = TRUE, direction = -1) +
  labs(x=bquote('Tree Height (m)'), y=bquote('Mortality Rate (% yr '^{-1}*')')) +
  theme(axis.title.x = element_text(vjust=-0.5)) +
  theme(axis.title.y = element_text(vjust=1.5)) +
  theme(plot.title = element_text(vjust=2.5)) +
  theme(legend.position='none')


gg<-ggplot(v_all, aes(y = zrt*100, x = temp, size = ct)) + 
  theme_bw() + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    stat_smooth(aes(weight=1/(se^2)), 
              method = "lm", se = TRUE,
              color = "black",
              size = 1) +
  geom_point(alpha = 1, shape = 21, fill = "red") + 
    scale_fill_viridis(name = "Available Water Capacity",
                     discrete = TRUE, direction = -1) +
  labs(x=bquote('Available Water Capacity (mm)'), y=bquote('Slope of mortality-height relationship')) +
  theme(axis.title.x = element_text(vjust=-0.5)) +
  theme(axis.title.y = element_text(vjust=1.5)) +
  theme(plot.title = element_text(vjust=2.5)) +
  theme(legend.position='none')

gg1<-ggplot_build(gg)
str(gg1)
gg2<-gg1$data[[1]]

# rsqr<-summary(lm(v_all[,2]~v_all[,1], weights = v_all[,3]))$r.squared
rsqr<-summary(lm(v_all[,2]~v_all[,1], weights = 1/v_all$se^2))$adj.r.squared
p_value<-summary(lm(v_all[,2]~v_all[,1], weights =  1/v_all$se^2))$coefficients[2,4]
beta1<-summary(lm(v_all[,2]~v_all[,1], weights =  1/v_all$se^2))$coefficients[2,1]*100
beta_sigma<-summary(lm(v_all[,2]~v_all[,1], weights =  1/v_all$se^2))$coefficients[2,2]*100

if(p_value<0.001) p_value <- " < 0.001" else p_value <- paste(' ==',round(p_value,3))

stats_ls[[1]]<-data.frame(var = "aws", beta = beta1, sigma = beta_sigma, r2 = rsqr, p = p_value)

aws<-ggplot(v_all, aes(y = zrt*100, x = temp, size = ct)) + 
  theme_bw() + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    stat_smooth(aes(weight=1/(se^2)), 
              method = "lm", se = FALSE,
              color = "black",
              size = 1) +
  coord_cartesian(ylim=c(0,1.5))+
  geom_line(data = gg2, aes(x = x, y = ymin), size = 0.5, linetype = 2) +
  geom_line(data = gg2, aes(x = x, y = ymax), size = 0.5, linetype = 2) +
  geom_errorbar(aes(ymin = minci*100, ymax = maxci*100), size = 0.5, width = 0)+
  geom_point(alpha = 1, shape = 21, fill = "red") + 
  scale_size(name = "", breaks = c(10000,50000,100000,250000,500000), labels = c("<10,000","50,000","100,000","250,000",">500,000"))+
    scale_fill_viridis(name = "Mean Temperature (degrees C)",
                     discrete = TRUE, direction = -1) +
  labs(x=bquote('Available Water Storage (mm)'), y=bquote('Slope of mortality-height relationship')) +
  theme(axis.title.x = element_text(vjust=-0.5)) +
  theme(axis.title.y = element_text(vjust=1.5)) +
  theme(plot.title = element_text(vjust=2.5),axis.text=element_text(color = "black", size=8.75)) +
  annotate("text",-Inf, Inf, label = paste('italic(R)^2==',round(rsqr,2), sep = ""), parse = TRUE, vjust = 1.5, hjust = -0.2) +
  annotate("text",-Inf, Inf, label = paste("italic(P)",p_value), parse = TRUE, vjust = 4, hjust = -0.2) +
  theme(legend.position='none')


### TMAX   ####
mort<-aggregate(count ~ zmax + tmax, FUN = "sum", ALLtrees_bin, drop = FALSE, na.action = NULL)
mort$dead2016<-aggregate(dead ~ zmax + tmax, FUN = "sum", ALLtrees_bin, drop = FALSE, na.action = NULL)[,3]
mort$rt<-mort$dead2016/mort$count/2
mort$rt[is.nan(mort$rt)]<-NA
mort<-mort[mort$count>10,]

v_ls<-list()
range<-seq(round(min(ALLtrees_bin$tmax,na.rm = TRUE),0),round(max(ALLtrees_bin$tmax,na.rm = TRUE),0),by=1)
for(i in 1:length(range)){
  if(nrow(na.exclude(mort[mort$tmax>range[i]&mort$tmax<=range[i+1],]))==0) next
  fit_tree<-lm(rt ~ zmax, weights = count, data = na.exclude(mort[mort$tmax>range[i]&mort$tmax<=range[i+1],]))
  v_ls[[i]] <- data.frame(temp = range[i], zrt=coef(fit_tree)[2], ct = sum(mort$count[mort$tmax>range[i]&mort$tmax<=range[i+1]], na.rm = TRUE),
                          se = sigma(fit_tree), minci = confint(fit_tree)[2,1], maxci = confint(fit_tree)[2,2], p = summary(fit_tree)$coefficients[2,4])
}

v_all<-do.call(rbind,v_ls)
v_all<-v_all[v_all$p<0.05,]

ggplot(mort, aes(y = rt*100, x = zmax, group = factor(tmax), color = factor(tmax))) + 
  theme_bw() + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    stat_smooth(aes(group = factor(0.5*floor(tmax/0.5)), weight=count), 
              method = "lm", se = FALSE,
              size = 0.5) +
  geom_jitter(aes(size = count), width = 5, alpha = 1, shape = 21) + 
    scale_color_viridis(name = "",
                      discrete = TRUE, direction = -1) +
  labs(x=bquote('Tree Height (m)'), y=bquote('Mortality Rate (% yr '^{-1}*')')) +
  theme(axis.title.x = element_text(vjust=-0.5)) +
  theme(axis.title.y = element_text(vjust=1.5)) +
  theme(plot.title = element_text(vjust=2.5)) +
  theme(legend.position='none')

gg<-ggplot(v_all, aes(y = zrt*100, x = temp, size = ct)) + 
  theme_bw() + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    stat_smooth(aes(weight=1/(se^2)), 
              method = "lm", se = TRUE,
              color = "black",
              size = 1) +
  geom_point(alpha = 1, shape = 21, fill = "red") + 
    scale_fill_viridis(name = "Max Temperature (degrees C)",
                     discrete = TRUE, direction = -1) +
  labs(x=bquote('Max Temperature (degrees C)'), y=bquote('Slope of mortality-height relationship')) +
  theme(axis.title.x = element_text(vjust=-0.5)) +
  theme(axis.title.y = element_text(vjust=1.5)) +
  theme(plot.title = element_text(vjust=2.5)) +
  theme(legend.position='none')

gg1<-ggplot_build(gg)
str(gg1)
gg2<-gg1$data[[1]]

# rsqr<-summary(lm(v_all[,2]~v_all[,1], weights = v_all[,3]))$r.squared
rsqr<-summary(lm(v_all[,2]~v_all[,1], weights = 1/v_all$se^2))$adj.r.squared
p_value<-summary(lm(v_all[,2]~v_all[,1], weights = 1/v_all$se^2))$coefficients[2,4]
beta1<-summary(lm(v_all[,2]~v_all[,1], weights = 1/v_all$se^2))$coefficients[2,1]*100
beta_sigma<-summary(lm(v_all[,2]~v_all[,1], weights =  1/v_all$se^2))$coefficients[2,2]*100

if(p_value<0.001) p_value <- " < 0.001" else p_value <- paste(' ==',round(p_value,3))

stats_ls[[2]]<-data.frame(var = "tmax", beta = beta1, sigma = beta_sigma, r2 = rsqr, p = p_value)

tmax<-ggplot(v_all, aes(y = zrt*100, x = temp, size = ct)) + 
  theme_bw() + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  stat_smooth(aes(weight=1/(se^2)), 
              method = "lm", se = FALSE,
              color = "black",
              size = 1) +
  coord_cartesian(ylim=c(0,1.5))+
  geom_line(data = gg2, aes(x = x, y = ymin), size = 0.5, linetype = 2) +
  geom_line(data = gg2, aes(x = x, y = ymax), size = 0.5, linetype = 2) +
  geom_errorbar(aes(ymin = minci*100, ymax = maxci*100), size = 0.5, width = 0)+
  geom_point(alpha = 1, shape = 21, fill = "red") + 
  scale_size(name = "", breaks = c(10000,50000,100000,250000,500000), labels = c("<10,000","50,000","100,000","250,000",">500,000"))+
    scale_fill_viridis(name = "Mean Temperature (degrees C)",
                     discrete = TRUE, direction = -1) +
  labs(x=bquote('Drought Temperature ('^o*' C)'), y=bquote('Slope of mortality-height relationship')) +
  theme(axis.title.x = element_text(vjust=-0.5)) +
  theme(axis.title.y = element_text(vjust=1.5)) +
  theme(plot.title = element_text(vjust=2.5),axis.text=element_text(color = "black", size=8.75)) +
  annotate("text",Inf, Inf, label = paste('italic(R)^2==',round(rsqr,2), sep = ""), parse = TRUE, vjust = 1.5, hjust = 1.25) +
  annotate("text",Inf, Inf, label = paste("italic(P)",p_value), parse = TRUE, vjust = 4, hjust = 1.28) +
    theme(legend.position='none')

### tmaxm  ####
mort<-aggregate(count ~ zmax + tmaxm, FUN = "sum", ALLtrees_bin, drop = FALSE, na.action = NULL)
mort$dead2016<-aggregate(dead ~ zmax + tmaxm, FUN = "sum", ALLtrees_bin, drop = FALSE, na.action = NULL)[,3]
mort$rt<-mort$dead2016/mort$count/2
mort$rt[is.nan(mort$rt)]<-NA
mort<-mort[mort$count>10,]

v_ls<-list()
range<-seq(round(min(ALLtrees_bin$tmaxm,na.rm = TRUE),0),round(max(ALLtrees_bin$tmaxm,na.rm = TRUE),0),by=1)
for(i in 1:length(range)){
  if(nrow(na.exclude(mort[mort$tmaxm>range[i]&mort$tmaxm<=range[i+1],]))==0) next
  fit_tree<-lm(rt ~ zmax, weights = count, data = na.exclude(mort[mort$tmaxm>range[i]&mort$tmaxm<=range[i+1],]))
  v_ls[[i]] <- data.frame(temp = range[i], zrt=coef(fit_tree)[2], ct = sum(mort$count[mort$tmaxm>range[i]&mort$tmaxm<=range[i+1]], na.rm = TRUE),
                          se = sigma(fit_tree), minci = confint(fit_tree)[2,1], maxci = confint(fit_tree)[2,2], p = summary(fit_tree)$coefficients[2,4])
}

v_all<-do.call(rbind,v_ls)
v_all<-v_all[v_all$p<0.05,]


ggplot(mort, aes(y = rt*100, x = zmax, group = factor(tmaxm), color = factor(tmaxm))) + 
  theme_bw() + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  stat_smooth(aes(group = factor(0.5*floor(tmaxm/0.5)), weight=count), 
              method = "lm", se = FALSE,
              # formula = y ~ poly(x, 2), 
              # color = "black", 
              size = 0.5) +
  geom_jitter(aes(size = count), width = 5, alpha = 1, shape = 21) + 
  scale_color_viridis(name = "",
                      discrete = TRUE, direction = -1) +
  labs(x=bquote('Tree Height (m)'), y=bquote('Mortality Rate (% yr '^{-1}*')')) +
  theme(axis.title.x = element_text(vjust=-0.5)) +
  theme(axis.title.y = element_text(vjust=1.5)) +
  theme(plot.title = element_text(vjust=2.5)) +
  theme(legend.position='none')


gg<-ggplot(v_all, aes(y = zrt*100, x = temp, size = ct)) + 
  theme_bw() + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  stat_smooth(aes(weight=1/(se^2)), 
              method = "lm", se = TRUE,
              # formula = y ~ poly(x, 2),
              color = "black",
              size = 1) +
  geom_point(alpha = 1, shape = 21, fill = "red") + 
  scale_fill_viridis(name = "",
                     discrete = TRUE, direction = -1) +
  labs(x=bquote('Drought Temperature ('^o*' C)'), y=bquote('Slope of mortality-height relationship')) +
  theme(axis.title.x = element_text(vjust=-0.5)) +
  theme(axis.title.y = element_text(vjust=1.5)) +
  theme(plot.title = element_text(vjust=2.5)) +
  theme(legend.position='none')

gg1<-ggplot_build(gg)
str(gg1)
gg2<-gg1$data[[1]]

# rsqr<-summary(lm(v_all[,2]~v_all[,1], weights = v_all[,3]))$r.squared
rsqr<-summary(lm(v_all[,2]~v_all[,1], weights = 1/v_all$se^2))$adj.r.squared
p_value<-summary(lm(v_all[,2]~v_all[,1], weights = 1/v_all$se^2))$coefficients[2,4]
beta1<-summary(lm(v_all[,2]~v_all[,1], weights = 1/v_all$se^2))$coefficients[2,1]*100
beta_sigma<-summary(lm(v_all[,2]~v_all[,1], weights =  1/v_all$se^2))$coefficients[2,2]*100

if(p_value<0.001) p_value <- " < 0.001" else p_value <- paste(' ==',round(p_value,3))

stats_ls[[4]]<-data.frame(var = "tmaxm", beta = beta1, sigma = beta_sigma, r2 = rsqr, p = p_value)

tmaxm<-ggplot(v_all, aes(y = zrt*100, x = temp, size = ct)) + 
  theme_bw() + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  stat_smooth(aes(weight=1/(se^2)), 
              method = "lm", se = FALSE,

              color = "black",
              size = 1) +
  coord_cartesian(ylim=c(0,1.5))+
  geom_line(data = gg2, aes(x = x, y = ymin), size = 0.5, linetype = 2) +
  geom_line(data = gg2, aes(x = x, y = ymax), size = 0.5, linetype = 2) +
  geom_errorbar(aes(ymin = minci*100, ymax = maxci*100), size = 0.5, width = 0)+
  geom_point(alpha = 1, shape = 21, fill = "red") + 
  scale_size(name = "", breaks = c(10000,50000,100000,250000,500000), labels = c("<10,000","50,000","100,000","250,000",">500,000"))+
  scale_fill_viridis(name = "Mean Temperature (degrees C)",
                     discrete = TRUE, direction = -1) +
  labs(x=bquote('Drought Temperature ('^o*' C)'), y=bquote('Slope of mortality-height relationship')) +
  theme(axis.title.x = element_text(vjust=-0.5)) +
  theme(axis.title.y = element_text(vjust=1.5)) +
  theme(plot.title = element_text(vjust=2.5),axis.text=element_text(color = "black", size=8.75)) +
  annotate("text",Inf, Inf, label = paste('italic(R)^2==',round(rsqr,2), sep = ""), parse = TRUE, vjust = 1.5, hjust = 1.25) +
  annotate("text",Inf, Inf, label = paste("italic(P)",p_value), parse = TRUE, vjust = 4, hjust = 1.28) +
  theme(legend.position='none')

### tmaxn  ####
mort<-aggregate(count ~ zmax + tmaxn, FUN = "sum", ALLtrees_bin, drop = FALSE, na.action = NULL)
mort$dead2016<-aggregate(dead ~ zmax + tmaxn, FUN = "sum", ALLtrees_bin, drop = FALSE, na.action = NULL)[,3]
mort$rt<-mort$dead2016/mort$count/2
mort$rt[is.nan(mort$rt)]<-NA
mort<-mort[mort$count>10,]

v_ls<-list()
range<-seq(round(min(ALLtrees_bin$tmaxn,na.rm = TRUE),0),round(max(ALLtrees_bin$tmaxn,na.rm = TRUE),0),by=1)
for(i in 1:length(range)){
  if(nrow(na.exclude(mort[mort$tmaxn>range[i]&mort$tmaxn<=range[i+1],]))==0) next
  fit_tree<-lm(rt ~ zmax, weights = count, data = na.exclude(mort[mort$tmaxn>range[i]&mort$tmaxn<=range[i+1],]))
  v_ls[[i]] <- data.frame(temp = range[i], zrt=coef(fit_tree)[2], ct = sum(mort$count[mort$tmaxn>range[i]&mort$tmaxn<=range[i+1]], na.rm = TRUE),
                          se = sigma(fit_tree), minci = confint(fit_tree)[2,1], maxci = confint(fit_tree)[2,2], p = summary(fit_tree)$coefficients[2,4])
}

v_all<-do.call(rbind,v_ls)
v_all<-v_all[v_all$p<0.05,]


ggplot(mort, aes(y = rt*100, x = zmax, group = factor(tmaxn), color = factor(tmaxn))) + 
  theme_bw() + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    stat_smooth(aes(group = factor(0.5*floor(tmaxn/0.5)), weight=count), 
              method = "lm", se = FALSE,
              size = 0.5) +
  geom_jitter(aes(size = count), width = 5, alpha = 1, shape = 21) + 
    scale_color_viridis(name = "",
                      discrete = TRUE, direction = -1) +
  labs(x=bquote('Tree Height (m)'), y=bquote('Mortality Rate (% yr '^{-1}*')')) +
  theme(axis.title.x = element_text(vjust=-0.5)) +
  theme(axis.title.y = element_text(vjust=1.5)) +
  theme(plot.title = element_text(vjust=2.5)) +
  theme(legend.position='none')

gg<-ggplot(v_all, aes(y = zrt*100, x = temp, size = ct)) + 
  theme_bw() + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    stat_smooth(aes(weight=1/(se^2)), 
              method = "lm", se = TRUE,
              # formula = y ~ poly(x, 2),
              color = "black",
              size = 1) +
  geom_point(alpha = 1, shape = 21, fill = "red") + 
  scale_fill_viridis(name = "",
                     discrete = TRUE, direction = -1) +
  labs(x=bquote('Historical Temperature ('^o*' C)'), y=bquote('Slope of mortality-height relationship')) +
  theme(axis.title.x = element_text(vjust=-0.5)) +
  theme(axis.title.y = element_text(vjust=1.5)) +
  theme(plot.title = element_text(vjust=2.5)) +
  theme(legend.position='none')

gg1<-ggplot_build(gg)
str(gg1)
gg2<-gg1$data[[1]]

# rsqr<-summary(lm(v_all[,2]~v_all[,1], weights = v_all[,3]))$r.squared
rsqr<-summary(lm(v_all[,2]~v_all[,1], weights = 1/v_all$se^2))$adj.r.squared
p_value<-summary(lm(v_all[,2]~v_all[,1], weights = 1/v_all$se^2))$coefficients[2,4]
beta1<-summary(lm(v_all[,2]~v_all[,1], weights = 1/v_all$se^2))$coefficients[2,1]*100
beta_sigma<-summary(lm(v_all[,2]~v_all[,1], weights =  1/v_all$se^2))$coefficients[2,2]*100

if(p_value<0.001) p_value <- " < 0.001" else p_value <- paste(' ==',round(p_value,3))

stats_ls[[4]]<-data.frame(var = "tmaxn", beta = beta1, sigma = beta_sigma, r2 = rsqr, p = p_value)

tmaxn<-ggplot(v_all, aes(y = zrt*100, x = temp, size = ct)) + 
  theme_bw() + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  stat_smooth(aes(weight=1/(se^2)), 
              method = "lm", se = FALSE,
              color = "black",
              size = 1) +
  coord_cartesian(ylim=c(0,1.5))+
  geom_line(data = gg2, aes(x = x, y = ymin), size = 0.5, linetype = 2) +
  geom_line(data = gg2, aes(x = x, y = ymax), size = 0.5, linetype = 2) +
  geom_errorbar(aes(ymin = minci*100, ymax = maxci*100), size = 0.5, width = 0)+
  geom_point(alpha = 1, shape = 21, fill = "red") + 
  scale_size(name = "", breaks = c(10000,50000,100000,250000,500000), labels = c("<10,000","50,000","100,000","250,000",">500,000"))+
  scale_fill_viridis(name = "Mean Temperature (degrees C)",
                     discrete = TRUE, direction = -1) +
  labs(x=bquote('Historical Temperature ('^o*' C)'), y=bquote('Slope of mortality-height relationship')) +
  theme(axis.title.x = element_text(vjust=-0.5)) +
  theme(axis.title.y = element_text(vjust=1.5)) +
  theme(plot.title = element_text(vjust=2.5),axis.text=element_text(color = "black", size=8.75)) +
  annotate("text",Inf, Inf, label = paste('italic(R)^2==',round(rsqr,2), sep = ""), parse = TRUE, vjust = 1.5, hjust = 1.25) +
  annotate("text",Inf, Inf, label = paste("italic(P)",p_value), parse = TRUE, vjust = 4, hjust = 1.28) +
  theme(legend.position='none')


### pptn   ####
mort<-aggregate(count ~ zmax + pptn, FUN = "sum", ALLtrees_bin, drop = FALSE, na.action = NULL)
mort$dead2016<-aggregate(dead ~ zmax + pptn, FUN = "sum", ALLtrees_bin, drop = FALSE, na.action = NULL)[,3]
mort$rt<-mort$dead2016/mort$count/2
mort$rt[is.nan(mort$rt)]<-NA
mort<-mort[mort$count>10,]

v_ls<-list()
range<-seq(round(min(ALLtrees_bin$pptn,na.rm = TRUE),0),round(max(ALLtrees_bin$pptn,na.rm = TRUE),0),by=1)
for(i in 1:length(range)){
  if(nrow(na.exclude(mort[mort$pptn>range[i]&mort$pptn<=range[i+1],]))==0) next
  fit_tree<-lm(rt ~ zmax, weights = count, data = na.exclude(mort[mort$pptn>range[i]&mort$pptn<=range[i+1],]))
  v_ls[[i]] <- data.frame(temp = range[i], zrt=coef(fit_tree)[2], ct = sum(mort$count[mort$pptn>range[i]&mort$pptn<=range[i+1]], na.rm = TRUE),
                          se = sigma(fit_tree), minci = confint(fit_tree)[2,1], maxci = confint(fit_tree)[2,2], p = summary(fit_tree)$coefficients[2,4])
}

v_all<-do.call(rbind,v_ls)
v_all<-v_all[v_all$p<0.05,]


ggplot(mort, aes(y = rt*100, x = zmax, group = factor(pptn), color = factor(pptn))) + 
  theme_bw() + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  stat_smooth(aes(group = factor(0.5*floor(pptn/0.5)), weight=count), 
              method = "lm", se = FALSE,
              size = 0.5) +
  geom_jitter(aes(size = count), width = 5, alpha = 1, shape = 21) + 
  scale_color_viridis(name = "Mean Temperature (degrees C)",
                      discrete = TRUE, direction = -1) +
  labs(x=bquote('Tree Height (m)'), y=bquote('Mortality Rate (% yr '^{-1}*')')) +
  theme(axis.title.x = element_text(vjust=-0.5)) +
  theme(axis.title.y = element_text(vjust=1.5)) +
  theme(plot.title = element_text(vjust=2.5)) +
  theme(legend.position='none')
ggsave("Mort_rt_z_soil_all.pdf", width = 3.75, height = 3.75, units = "in")

gg<-ggplot(v_all, aes(y = zrt*100, x = temp, size = ct)) + 
  theme_bw() + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  stat_smooth(aes(weight=1/(se^2)), 
              method = "lm", se = TRUE,
              # formula = y ~ poly(x, 2),
              color = "black",
              size = 1) +
  geom_point(alpha = 1, shape = 21, fill = "red") + 
  scale_fill_viridis(name = "Mean Temperature (degrees C)",
                     discrete = TRUE, direction = -1) +
  labs(x=bquote('Available Water Capacity (mm)'), y=bquote('Slope of mortality-height relationship')) +
  theme(axis.title.x = element_text(vjust=-0.5)) +
  theme(axis.title.y = element_text(vjust=1.5)) +
  theme(plot.title = element_text(vjust=2.5)) +
  theme(legend.position='none')

gg1<-ggplot_build(gg)
str(gg1)
gg2<-gg1$data[[1]]

# rsqr<-summary(lm(v_all[,2]~v_all[,1], weights = v_all[,3]))$r.squared
rsqr<-summary(lm(v_all[,2]~v_all[,1], weights = 1/v_all$se^2))$adj.r.squared
p_value<-summary(lm(v_all[,2]~v_all[,1], weights = 1/v_all$se^2))$coefficients[2,4]
beta1<-summary(lm(v_all[,2]~v_all[,1], weights = 1/v_all$se^2))$coefficients[2,1]*100
beta_sigma<-summary(lm(v_all[,2]~v_all[,1], weights =  1/v_all$se^2))$coefficients[2,2]*100

if(p_value<0.001) p_value <- " < 0.001" else p_value <- paste(' ==',round(p_value,3))

stats_ls[[5]]<-data.frame(var = "pptn", beta = beta1, sigma = beta_sigma, r2 = rsqr, p = p_value)

pptn<-ggplot(v_all, aes(y = zrt*100, x = temp, size = ct)) + 
  theme_bw() + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  stat_smooth(aes(weight=1/(se^2)), 
              method = "lm", se = FALSE,
              # formula = y ~ poly(x, 2),
              color = "black",
              size = 1) +
  coord_cartesian(ylim=c(0,1.5))+
  geom_line(data = gg2, aes(x = x, y = ymin), size = 0.5, linetype = 2) +
  geom_line(data = gg2, aes(x = x, y = ymax), size = 0.5, linetype = 2) +
  geom_errorbar(aes(ymin = minci*100, ymax = maxci*100), size = 0.5, width = 0)+
  geom_point(alpha = 1, shape = 21, fill = "red") + 
  scale_size(name = "", breaks = c(10000,50000,100000,250000,500000), labels = c("<10,000","50,000","100,000","250,000",">500,000"))+
  scale_fill_viridis(name = "Mean Temperature (degrees C)",
                     discrete = TRUE, direction = -1) +
  labs(x=bquote('Historical Precipitation (mm)'), y=bquote('Slope of mortality-height relationship')) +
  theme(axis.title.x = element_text(vjust=-0.5)) +
  theme(axis.title.y = element_text(vjust=1.5)) +
  theme(plot.title = element_text(vjust=2.5),axis.text=element_text(color = "black", size=8.75)) +
  annotate("text",Inf, Inf, label = paste('italic(R)^2==',round(rsqr,2), sep = ""), parse = TRUE, vjust = 1.5, hjust = 1.25) +
  annotate("text",Inf, Inf, label = paste("italic(P)",p_value), parse = TRUE, vjust = 4, hjust = 1.28) +
  theme(legend.position='none')


### ppt    ####
mort<-aggregate(count ~ zmax + ppt, FUN = "sum", ALLtrees_bin, drop = FALSE, na.action = NULL)
mort$dead2016<-aggregate(dead ~ zmax + ppt, FUN = "sum", ALLtrees_bin, drop = FALSE, na.action = NULL)[,3]
mort$rt<-mort$dead2016/mort$count/2
mort$rt[is.nan(mort$rt)]<-NA
mort<-mort[mort$count>10,]

v_ls<-list()
range<-seq(round(min(ALLtrees_bin$ppt,na.rm = TRUE),0),round(max(ALLtrees_bin$ppt,na.rm = TRUE),0),by=1)
range<-seq(round(min(ALLtrees_bin$ppt,na.rm = TRUE),2),round(max(ALLtrees_bin$ppt,na.rm = TRUE),2),by=0.01)
for(i in 1:length(range)){
  if(nrow(na.exclude(mort[mort$ppt>range[i]&mort$ppt<=range[i+1],]))==0) next
  fit_tree<-lm(rt ~ zmax, weights = count, data = na.exclude(mort[mort$ppt>range[i]&mort$ppt<=range[i+1],]))
  v_ls[[i]] <- data.frame(temp = range[i], zrt=coef(fit_tree)[2], ct = sum(mort$count[mort$ppt>range[i]&mort$ppt<=range[i+1]], na.rm = TRUE),
                          se = sigma(fit_tree), minci = confint(fit_tree)[2,1], maxci = confint(fit_tree)[2,2], p = summary(fit_tree)$coefficients[2,4])
}

v_all<-do.call(rbind,v_ls)
v_all<-v_all[v_all$p<0.05,]


ggplot(mort, aes(y = rt*100, x = zmax, group = factor(ppt), color = factor(ppt))) + 
  theme_bw() + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  stat_smooth(aes(group = factor(0.5*floor(ppt/0.5)), weight=count), 
              method = "lm", se = FALSE,
              size = 0.5) +
  geom_jitter(aes(size = count), width = 5, alpha = 1, shape = 21) +
  scale_color_viridis(name = "",
                      discrete = TRUE, direction = -1) +
  labs(x=bquote('Tree Height (m)'), y=bquote('Mortality Rate (% yr '^{-1}*')')) +
  theme(axis.title.x = element_text(vjust=-0.5)) +
  theme(axis.title.y = element_text(vjust=1.5)) +
  theme(plot.title = element_text(vjust=2.5)) +
  theme(legend.position='none')

gg<-ggplot(v_all, aes(y = zrt*100, x = temp, size = ct)) + 
  theme_bw() + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  stat_smooth(aes(weight=1/(se^2)), 
              method = "lm", se = TRUE,
              color = "black",
              size = 1) +
  geom_point(alpha = 1, shape = 21, fill = "red") + 
  scale_fill_viridis(name = "Mean Temperature (degrees C)",
                     discrete = TRUE, direction = -1) +
  labs(x=bquote('Precipitation (mm)'), y=bquote('Slope of mortality-height relationship')) +
  theme(axis.title.x = element_text(vjust=-0.5)) +
  theme(axis.title.y = element_text(vjust=1.5)) +
  theme(plot.title = element_text(vjust=2.5)) +
  theme(legend.position='none')

gg1<-ggplot_build(gg)
str(gg1)
gg2<-gg1$data[[1]]

# rsqr<-summary(lm(v_all[,2]~v_all[,1], weights = v_all[,3]))$r.squared
rsqr<-summary(lm(v_all[,2]~v_all[,1], weights = 1/v_all$se^2))$adj.r.squared
p_value<-summary(lm(v_all[,2]~v_all[,1], weights = 1/v_all$se^2))$coefficients[2,4]
beta1<-summary(lm(v_all[,2]~v_all[,1], weights = 1/v_all$se^2))$coefficients[2,1]*100
beta_sigma<-summary(lm(v_all[,2]~v_all[,1], weights =  1/v_all$se^2))$coefficients[2,2]*100

if(p_value<0.001) p_value <- " < 0.001" else p_value <- paste(' ==',round(p_value,3))

stats_ls[[6]]<-data.frame(var = "ppt", beta = beta1, sigma = beta_sigma, r2 = rsqr, p = p_value)

ppt<-ggplot(v_all, aes(y = zrt*100, x = temp, size = ct)) + 
  theme_bw() + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  stat_smooth(aes(weight=1/(se^2)), 
              method = "lm", se = FALSE,
              color = "black",
              size = 1) +
  coord_cartesian(ylim=c(0,1.5))+
  geom_line(data = gg2, aes(x = x, y = ymin), size = 0.5, linetype = 2) +
  geom_line(data = gg2, aes(x = x, y = ymax), size = 0.5, linetype = 2) +
  geom_errorbar(aes(ymin = minci*100, ymax = maxci*100), size = 0.5, width = 0)+
  geom_point(alpha = 1, shape = 21, fill = "red") + 
  scale_size(name = "", breaks = c(10000,50000,100000,250000,500000), labels = c("<10,000","50,000","100,000","250,000",">500,000"))+
  scale_fill_viridis(name = "",
                     discrete = TRUE, direction = -1) +
  labs(x=bquote('Drought Precipitation (mm)'), y=bquote('Slope of mortality-height relationship')) +
  theme(axis.title.x = element_text(vjust=-0.5)) +
  theme(axis.title.y = element_text(vjust=1.5)) +
  theme(plot.title = element_text(vjust=2.5),axis.text=element_text(color = "black", size=8.75)) +
  annotate("text",Inf, Inf, label = paste('italic(R)^2==',round(rsqr,2), sep = ""), parse = TRUE, vjust = 1.5, hjust = 1.25) +
  annotate("text",Inf, Inf, label = paste("italic(P)",p_value), parse = TRUE, vjust = 4, hjust = 1.28) +
  theme(legend.position='none')

### Aspect ####
mort<-aggregate(count ~ zmax + aspect, FUN = "sum", ALLtrees_bin, drop = FALSE, na.action = NULL)
mort$dead2016<-aggregate(dead ~ zmax + aspect, FUN = "sum", ALLtrees_bin, drop = FALSE, na.action = NULL)[,3]
mort$rt<-mort$dead2016/mort$count/2
mort$rt[is.nan(mort$rt)]<-NA
mort<-mort[mort$count>10,]

v_ls<-list()
range<-seq(round(min(ALLtrees_bin$aspect,na.rm = TRUE),0),round(max(ALLtrees_bin$aspect,na.rm = TRUE),0),by=1)
for(i in 1:length(range)){
  if(nrow(na.exclude(mort[mort$aspect>range[i]&mort$aspect<=range[i+1],]))==0) next
  fit_tree<-lm(rt ~ zmax, weights = count, data = na.exclude(mort[mort$aspect>range[i]&mort$aspect<=range[i+1],]))
  v_ls[[i]] <- data.frame(temp = range[i], zrt=coef(fit_tree)[2], ct = sum(mort$count[mort$aspect>range[i]&mort$aspect<=range[i+1]], na.rm = TRUE),
                          se = sigma(fit_tree), minci = confint(fit_tree)[2,1], maxci = confint(fit_tree)[2,2], p = summary(fit_tree)$coefficients[2,4])
}

v_all<-do.call(rbind,v_ls)
v_all<-v_all[v_all$p<0.05,]

ggplot(mort, aes(y = rt*100, x = zmax, group = factor(aspect), color = factor(aspect))) + 
  theme_bw() + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  stat_smooth(aes(group = factor(0.5*floor(aspect/0.5)), weight=count), 
              method = "lm", se = FALSE,
              size = 0.5) +
  geom_jitter(aes(size = count), width = 5, alpha = 1, shape = 21) + 
  scale_color_viridis(name = "Aspect (o)",
                      discrete = TRUE, direction = -1) +
  labs(x=bquote('Tree Height (m)'), y=bquote('Mortality Rate (% yr '^{-1}*')')) +
  theme(axis.title.x = element_text(vjust=-0.5)) +
  theme(axis.title.y = element_text(vjust=1.5)) +
  theme(plot.title = element_text(vjust=2.5)) +
  theme(legend.position='none')


gg<-ggplot(v_all, aes(y = zrt*100, x = temp, size = ct)) + 
  theme_bw() + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  stat_smooth(aes(weight=1/(se^2)), 
              method = "lm", se = TRUE,
              # formula = y ~ poly(x, 2),
              color = "black",
              size = 1) +
  geom_point(alpha = 1, shape = 21, fill = "red") + 
  scale_fill_viridis(name = "Mean Temperature (degrees C)",
                     discrete = TRUE, direction = -1) +
  labs(x=bquote('Aspect (o)'), y=bquote('Slope of mortality-height relationship')) +
  theme(axis.title.x = element_text(vjust=-0.5)) +
  theme(axis.title.y = element_text(vjust=1.5)) +
  theme(plot.title = element_text(vjust=2.5)) +
  theme(legend.position='none')

gg1<-ggplot_build(gg)
str(gg1)
gg2<-gg1$data[[1]]

# rsqr<-summary(lm(v_all[,2]~v_all[,1], weights = v_all[,3]))$r.squared
rsqr<-summary(lm(v_all[,2]~v_all[,1], weights = 1/v_all$se^2))$adj.r.squared
p_value<-summary(lm(v_all[,2]~v_all[,1], weights =  1/v_all$se^2))$coefficients[2,4]
beta1<-summary(lm(v_all[,2]~v_all[,1], weights =  1/v_all$se^2))$coefficients[2,1]*100
beta_sigma<-summary(lm(v_all[,2]~v_all[,1], weights =  1/v_all$se^2))$coefficients[2,2]*100

if(p_value<0.001) p_value <- " < 0.001" else p_value <- paste(' ==',round(p_value,3))

stats_ls[[7]]<-data.frame(var = "aspect", beta = beta1, sigma = beta_sigma, r2 = rsqr, p = p_value)

aspect<-ggplot(v_all, aes(y = zrt*100, x = temp, size = ct)) + 
  theme_bw() + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  stat_smooth(aes(weight=1/(se^2)), 
              method = "lm", se = FALSE,
              # formula = y ~ poly(x, 2),
              color = "black",
              size = 1) +
  coord_cartesian(ylim=c(0,1.5))+
  geom_line(data = gg2, aes(x = x, y = ymin), size = 0.5, linetype = 2) +
  geom_line(data = gg2, aes(x = x, y = ymax), size = 0.5, linetype = 2) +
  geom_errorbar(aes(ymin = minci*100, ymax = maxci*100), size = 0.5, width = 0)+
  geom_point(alpha = 1, shape = 21, fill = "red") + 
  scale_size(name = "", breaks = c(10000,50000,100000,250000,500000), labels = c("<10,000","50,000","100,000","250,000",">500,000"))+
  scale_fill_viridis(name = "",
                     discrete = TRUE, direction = -1) +
  labs(x=bquote('Aspect (o)'), y=bquote('Slope of mortality-height relationship')) +
  theme(axis.title.x = element_text(vjust=-0.5)) +
  theme(axis.title.y = element_text(vjust=1.5)) +
  theme(plot.title = element_text(vjust=2.5),axis.text=element_text(color = "black", size=8.75)) +
  annotate("text",Inf, Inf, label = paste('italic(R)^2==',round(rsqr,2), sep = ""), parse = TRUE, vjust = 1.5, hjust = 1.25) +
  annotate("text",Inf, Inf, label = paste("italic(P)",p_value), parse = TRUE, vjust = 4, hjust = 1.28) +
  theme(legend.position='none')

### slope  ####
mort<-aggregate(count ~ zmax + slope, FUN = "sum", ALLtrees_bin, drop = FALSE, na.action = NULL)
mort$dead2016<-aggregate(dead ~ zmax + slope, FUN = "sum", ALLtrees_bin, drop = FALSE, na.action = NULL)[,3]
mort$rt<-mort$dead2016/mort$count/2
mort$rt[is.nan(mort$rt)]<-NA
mort<-mort[mort$count>10,]

v_ls<-list()
range<-seq(round(min(ALLtrees_bin$slope,na.rm = TRUE),0),round(max(ALLtrees_bin$slope,na.rm = TRUE),0),by=1)
for(i in 1:length(range)){
  if(nrow(na.exclude(mort[mort$slope>range[i]&mort$slope<=range[i+1],]))==0) next
  fit_tree<-lm(rt ~ zmax, weights = count, data = na.exclude(mort[mort$slope>range[i]&mort$slope<=range[i+1],]))
  v_ls[[i]] <- data.frame(temp = range[i], zrt=coef(fit_tree)[2], ct = sum(mort$count[mort$slope>range[i]&mort$slope<=range[i+1]], na.rm = TRUE),
                          se = sigma(fit_tree), minci = confint(fit_tree)[2,1], maxci = confint(fit_tree)[2,2], p = summary(fit_tree)$coefficients[2,4])
}

v_all<-do.call(rbind,v_ls)
v_all<-v_all[v_all$p<0.05,]


ggplot(mort, aes(y = rt*100, x = zmax, group = factor(slope), color = factor(slope))) + 
  theme_bw() + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  stat_smooth(aes(group = factor(0.5*floor(slope/0.5)), weight=count), 
              method = "lm", se = FALSE,
              size = 0.5) +
  geom_jitter(aes(size = count), width = 5, alpha = 1, shape = 21) +
  scale_color_viridis(name = "Slope (degrees)",
                      discrete = TRUE, direction = -1) +
  labs(x=bquote('Tree Height (m)'), y=bquote('Mortality Rate (% yr '^{-1}*')')) +
  theme(axis.title.x = element_text(vjust=-0.5)) +
  theme(axis.title.y = element_text(vjust=1.5)) +
  theme(plot.title = element_text(vjust=2.5)) +
  theme(legend.position='none')


gg<-ggplot(v_all, aes(y = zrt*100, x = temp, size = ct)) + 
  theme_bw() + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  stat_smooth(aes(weight=1/(se^2)), 
              method = "lm", se = TRUE,
              color = "black",
              size = 1) +
  geom_point(alpha = 1, shape = 21, fill = "red") + 
  scale_fill_viridis(name = "Mean Temperature (degrees C)",
                     discrete = TRUE, direction = -1) +
  labs(x=bquote('Slope ('^o*')'), y=bquote('Slope of mortality-height relationship')) +
  theme(axis.title.x = element_text(vjust=-0.5)) +
  theme(axis.title.y = element_text(vjust=1.5)) +
  theme(plot.title = element_text(vjust=2.5)) +
  theme(legend.position='none')

gg1<-ggplot_build(gg)
str(gg1)
gg2<-gg1$data[[1]]

# rsqr<-summary(lm(v_all[,2]~v_all[,1], weights = v_all[,3]))$r.squared
rsqr<-summary(lm(v_all[,2]~v_all[,1], weights = 1/v_all$se^2))$adj.r.squared
p_value<-summary(lm(v_all[,2]~v_all[,1], weights = 1/v_all$se^2))$coefficients[2,4]
beta1<-summary(lm(v_all[,2]~v_all[,1], weights = 1/v_all$se^2))$coefficients[2,1]*100
beta_sigma<-summary(lm(v_all[,2]~v_all[,1], weights =  1/v_all$se^2))$coefficients[2,2]*100

if(p_value<0.001) p_value <- " < 0.001" else p_value <- paste(' ==',round(p_value,3))

stats_ls[[8]]<-data.frame(var = "slope", beta = beta1, sigma = beta_sigma, r2 = rsqr, p = p_value)

slope<-ggplot(v_all, aes(y = zrt*100, x = temp, size = ct)) + 
  theme_bw() + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  stat_smooth(aes(weight=1/(se^2)), 
              method = "lm", se = FALSE,
              color = "black",
              size = 1) +
  coord_cartesian(ylim=c(0,1.5))+
  geom_line(data = gg2, aes(x = x, y = ymin), size = 0.5, linetype = 2) +
  geom_line(data = gg2, aes(x = x, y = ymax), size = 0.5, linetype = 2) +
  geom_errorbar(aes(ymin = minci*100, ymax = maxci*100), size = 0.5, width = 0)+
  geom_point(alpha = 1, shape = 21, fill = "red") + 
  scale_size(name = "", breaks = c(10000,50000,100000,250000,500000), labels = c("<10,000","50,000","100,000","250,000",">500,000"))+
  scale_fill_viridis(name = "Mean Temperature (degrees C)",
                     discrete = TRUE, direction = -1) +
  labs(x=bquote('Slope ('^o*')'), y=bquote('Slope of mortality-height relationship')) +
  theme(axis.title.x = element_text(vjust=-0.5)) +
  theme(axis.title.y = element_text(vjust=1.5)) +
  theme(plot.title = element_text(vjust=2.5),axis.text=element_text(color = "black", size=8.75)) +
  annotate("text",Inf, Inf, label = paste('italic(R)^2==',round(rsqr,2), sep = ""), parse = TRUE, vjust = 1.5, hjust = 1.25) +
  annotate("text",Inf, Inf, label = paste("italic(P)",p_value), parse = TRUE, vjust = 4, hjust = 1.28) +
  theme(legend.position='none')



### cover  ####
mort<-aggregate(count ~ zmax + cover, FUN = "sum", ALLtrees_bin, drop = FALSE, na.action = NULL)
mort$dead2016<-aggregate(dead ~ zmax + cover, FUN = "sum", ALLtrees_bin, drop = FALSE, na.action = NULL)[,3]
mort$rt<-mort$dead2016/mort$count/2
mort$rt[is.nan(mort$rt)]<-NA
mort<-mort[mort$count>10,]

v_ls<-list()
range<-seq(round(min(ALLtrees_bin$cover,na.rm = TRUE),2),round(max(ALLtrees_bin$cover,na.rm = TRUE),2),by=0.01)
for(i in 1:length(range)){
  if(nrow(na.exclude(mort[mort$cover>range[i]&mort$cover<=range[i+1],]))==0) next
  fit_tree<-lm(rt ~ zmax, weights = count, data = na.exclude(mort[mort$cover>range[i]&mort$cover<=range[i+1],]))
  v_ls[[i]] <- data.frame(temp = range[i], zrt=coef(fit_tree)[2], ct = sum(mort$count[mort$cover>range[i]&mort$cover<=range[i+1]], na.rm = TRUE),
                          se = sigma(fit_tree), minci = confint(fit_tree)[2,1], maxci = confint(fit_tree)[2,2], p = summary(fit_tree)$coefficients[2,4])
}

v_all<-do.call(rbind,v_ls)
v_all<-v_all[v_all$p<0.05,]


ggplot(na.omit(mort), aes(y = rt*100, x = zmax, group = factor(cover), color = factor(cover))) + 
  theme_bw() + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  stat_smooth(aes(group = factor(0.1*floor(cover/0.1)), weight=count), 
              method = "lm", se = FALSE,
              size = 0.5) +
  geom_jitter(aes(size = count), width = 5, alpha = 1, shape = 21) +
  scale_color_viridis(name = "Cover (%)",
                      discrete = TRUE, direction = 1) +
  labs(x=bquote('Tree Height (m)'), y=bquote('Mortality Rate (% yr '^{-1}*')')) +
  theme(axis.title.x = element_text(vjust=-0.5)) +
  theme(axis.title.y = element_text(vjust=1.5)) +
  theme(plot.title = element_text(vjust=2.5)) +
  theme(legend.position='none')
    

gg<-ggplot(v_all, aes(y = zrt*100, x = temp*100, size = ct)) + 
  theme_bw() + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  
  #we also tried a polynomial fit, but felt two linear fits were more representitive of the trend
  
  # stat_smooth(aes(weight=1/(se^2)),
  #             method = "lm",
  #             formula=y ~ poly(x, 2, raw=TRUE),
  #             se = TRUE,
  #             color = "black",
  #             size = 1) +
  # stat_smooth(aes(weight=1/(se^2)),
  #             method = "lm", se = TRUE,
  #             color = "black",
  #             size = 1) +

  stat_smooth(data = v_all[v_all$temp<=0.65,], aes(weight=1/(se^2)),
            method = "lm",
            formula=y ~ x,
            se = TRUE,
            color = "black",
            size = 1) +
  stat_smooth(data = v_all[v_all$temp>=0.5,], aes(weight=1/(se^2)),
              method = "lm",
              formula=y ~ x,
              se = TRUE,
              color = "black",
              size = 1) +
  geom_point(alpha = 1, shape = 21, fill = "red") + 
  scale_fill_viridis(name = "Mean Temperature (degrees C)",
                     discrete = TRUE, direction = -1) +
  labs(x=bquote('Tree Cover(%)'), y=bquote('Slope of mortality-height relationship')) +
  theme(axis.title.x = element_text(vjust=-0.5)) +
  theme(axis.title.y = element_text(vjust=1.5)) +
  theme(plot.title = element_text(vjust=2.5)) +
  theme(legend.position='none')

gg1<-ggplot_build(gg)
str(gg1)
gg2<-gg1$data[[1]]
gg2<-rbind(gg2,gg1$data[[2]])




v_all1<-v_all

v_all<-v_all1[v_all1$temp<=0.65,]

rsqr<-summary(lm(v_all[,2]~v_all[,1], weights = 1/v_all$se^2))$adj.r.squared
p_value<-summary(lm(v_all[,2]~v_all[,1], weights = 1/v_all$se^2))$coefficients[2,4]
beta1<-summary(lm(v_all[,2]~v_all[,1], weights = 1/v_all$se^2))$coefficients[2,1]*100
beta_sigma<-summary(lm(v_all[,2]~v_all[,1], weights =  1/v_all$se^2))$coefficients[2,2]*100

if(p_value<0.001) p_value <- " < 0.001" else p_value <- paste(' ==',round(p_value,3))

v_all<-v_all1[v_all1$temp>=0.5,]

rsqr_2<-summary(lm(v_all[,2]~v_all[,1], weights = 1/v_all$se^2))$adj.r.squared
p_value_2<-summary(lm(v_all[,2]~v_all[,1], weights = 1/v_all$se^2))$coefficients[2,4]
beta1_2<-summary(lm(v_all[,2]~v_all[,1], weights = 1/v_all$se^2))$coefficients[2,1]*100
beta_sigma_2<-summary(lm(v_all[,2]~v_all[,1], weights =  1/v_all$se^2))$coefficients[2,2]*100

if(p_value_2<0.001) p_value_2 <- " < 0.001" else p_value_2 <- paste(' ==',round(p_value_2,3))

stats_ls[[9]]<-rbind(data.frame(var = "cover0_50", beta = beta1, sigma = beta_sigma, r2 = rsqr, p = p_value),
                     data.frame(var = "cover50_100", beta = beta1_2, sigma = beta_sigma_2, r2 = rsqr_2, p = p_value_2)
                     )

v_all<-v_all1

cover<-ggplot(v_all, aes(y = zrt*100, x = temp*100, size = ct)) + 
  theme_bw() + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  
  #we also tried a polynomial fit, but felt two linear fits were more representitive of the trend
  
  # stat_smooth(aes(weight=1/(se^2)),
  #             method = "lm", se = FALSE,
  #             color = "black",
  #             size = 1) +
  # stat_smooth(aes(weight=1/(se^2)),
  #             method = "lm",
  #             formula=y ~ poly(x, 2, raw=TRUE),
  #             se = FALSE,
  #             color = "black",
  #             size = 1) +

  stat_smooth(data = v_all[v_all$temp<=0.65,], aes(weight=1/(se^2)),
              method = "lm",
              formula=y ~ x,
              se = FALSE,
              color = "black",
              size = 1) +
  stat_smooth(data = v_all[v_all$temp>=0.5,], aes(weight=1/(se^2)),
              method = "lm",
              formula=y ~ x,
              se = FALSE,
              color = "black",
              size = 1) +
  coord_cartesian(ylim=c(0,1.5))+
  geom_line(data = gg2, aes(x = x, y = ymin), size = 0.5, linetype = 2) +
  geom_line(data = gg2, aes(x = x, y = ymax), size = 0.5, linetype = 2) +
  geom_errorbar(aes(ymin = minci*100, ymax = maxci*100), size = 0.5, width = 0)+
  geom_point(alpha = 1, shape = 21, fill = "red") + 
  scale_size(name = "", breaks = c(10000,50000,100000,250000,500000), labels = c("<10,000","50,000","100,000","250,000",">500,000"))+
  scale_fill_viridis(name = "Mean Temperature (degrees C)",
                     discrete = TRUE, direction = -1) +
  labs(x=bquote('Tree Cover (%)'), y=bquote('Slope of mortality-height relationship')) +
  theme(axis.title.x = element_text(vjust=-0.5)) +
  theme(axis.title.y = element_text(vjust=1.5)) +
  theme(plot.title = element_text(vjust=2.5),axis.text=element_text(color = "black", size=8.75)) +
  annotate("text",-Inf, Inf, label = paste('italic(R)^2==',round(rsqr,2), sep = ""), parse = TRUE, vjust = 1.5, hjust = -0.2) +
  annotate("text",-Inf, Inf, label = paste("italic(P)",p_value), parse = TRUE, vjust = 4, hjust = -0.17) +
  
  annotate("text",Inf, Inf, label = paste('italic(R)^2==',round(rsqr_2,2), sep = ""), parse = TRUE, vjust = 1.5, hjust = 1.25) +
  annotate("text",Inf, Inf, label = paste("italic(P)",p_value_2), parse = TRUE, vjust = 4, hjust = 1.28) +
  theme(legend.position='none')


### VPDMAX ####
mort<-aggregate(count ~ zmax + vpdmax, FUN = "sum", ALLtrees_bin, drop = FALSE, na.action = NULL)
mort$dead2016<-aggregate(dead ~ zmax + vpdmax, FUN = "sum", ALLtrees_bin, drop = FALSE, na.action = NULL)[,3]
mort$rt<-mort$dead2016/mort$count/2
mort$rt[is.nan(mort$rt)]<-NA
mort<-mort[mort$count>10,]

v_ls<-list()
range<-seq(round(min(ALLtrees_bin$vpdmax,na.rm = TRUE),0),round(max(ALLtrees_bin$vpdmax,na.rm = TRUE),0),by=1)
for(i in 1:length(range)){
  if(nrow(na.exclude(mort[mort$vpdmax>range[i]&mort$vpdmax<=range[i+1],]))==0) next
  fit_tree<-lm(rt ~ zmax, weights = count, data = na.exclude(mort[mort$vpdmax>range[i]&mort$vpdmax<=range[i+1],]))
  v_ls[[i]] <- data.frame(temp = range[i], zrt=coef(fit_tree)[2], ct = sum(mort$count[mort$vpdmax>range[i]&mort$vpdmax<=range[i+1]], na.rm = TRUE),
                          se = sigma(fit_tree), minci = confint(fit_tree)[2,1], maxci = confint(fit_tree)[2,2], p = summary(fit_tree)$coefficients[2,4])
}

v_all<-do.call(rbind,v_ls)
v_all<-v_all[v_all$p<0.05,]

ggplot(mort, aes(y = rt*100, x = zmax, group = factor(vpdmax), color = factor(vpdmax))) + 
  theme_bw() + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  stat_smooth(aes(group = factor(0.5*floor(vpdmax/0.5)), weight=count), 
              method = "lm", se = FALSE,
              # formula = y ~ poly(x, 2), 
              # color = "black", 
              size = 0.5) +
  geom_jitter(aes(size = count), width = 5, alpha = 1, shape = 21) + 
  scale_color_viridis(name = "Mean Temperature (degrees C)",
                      discrete = TRUE, direction = -1) +
  labs(x=bquote('Tree Height (m)'), y=bquote('Mortality Rate (% yr '^{-1}*')')) +
  theme(axis.title.x = element_text(vjust=-0.5)) +
  theme(axis.title.y = element_text(vjust=1.5)) +
  theme(plot.title = element_text(vjust=2.5)) +
  theme(legend.position='none')

gg<-ggplot(v_all, aes(y = zrt*100, x = temp*0.10, size = ct)) + 
  theme_bw() + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  stat_smooth(aes(weight=1/(se^2)), 
              method = "lm", se = TRUE,
              # formula = y ~ poly(x, 2),
              color = "black",
              size = 1) +
  geom_point(alpha = 1, shape = 21, fill = "red") + 
  scale_fill_viridis(name = "Mean Temperature (degrees C)",
                     discrete = TRUE, direction = -1) +
  labs(x=bquote('VPD_Max (kPa)'), y=bquote('Slope of mortality-height relationship')) +
  theme(axis.title.x = element_text(vjust=-0.5)) +
  theme(axis.title.y = element_text(vjust=1.5)) +
  theme(plot.title = element_text(vjust=2.5)) +
  theme(legend.position='none')

gg1<-ggplot_build(gg)
str(gg1)
gg2<-gg1$data[[1]]

# rsqr<-summary(lm(v_all[,2]~v_all[,1], weights = v_all[,3]))$r.squared
rsqr<-summary(lm(v_all[,2]~v_all[,1], weights = 1/v_all$se^2))$adj.r.squared
p_value<-summary(lm(v_all[,2]~v_all[,1], weights = 1/v_all$se^2))$coefficients[2,4]
beta1<-summary(lm(v_all[,2]~v_all[,1], weights = 1/v_all$se^2))$coefficients[2,1]*100
beta_sigma<-summary(lm(v_all[,2]~v_all[,1], weights =  1/v_all$se^2))$coefficients[2,2]*100

if(p_value<0.001) p_value <- " < 0.001" else p_value <- paste(' ==',round(p_value,3))

stats_ls[[10]]<-data.frame(var = "vpdmax", beta = beta1, sigma = beta_sigma, r2 = rsqr, p = p_value)

vpdmax<-ggplot(v_all, aes(y = zrt*100, x = temp*0.10, size = ct)) + 
  theme_bw() + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  stat_smooth(aes(weight=1/(se^2)), 
              method = "lm", se = FALSE,
              color = "black",
              size = 1) +
  geom_line(data = gg2, aes(x = x, y = ymin), size = 0.5, linetype = 2) +
  geom_line(data = gg2, aes(x = x, y = ymax), size = 0.5, linetype = 2) +
  geom_errorbar(aes(ymin = minci*100, ymax = maxci*100), size = 0.5, width = 0)+
  geom_point(alpha = 1, shape = 21, fill = "red") + 
  scale_size(name = "", breaks = c(10000,50000,100000,250000,500000), labels = c("<10,000","50,000","100,000","250,000",">500,000"))+
  scale_fill_viridis(name = "Mean Temperature (degrees C)",
                     discrete = TRUE, direction = -1) +
  labs(x=bquote('Vapor Pressure Deficit (kPa)'), y=bquote('Slope of mortality-height relationship')) +
  theme(axis.title.x = element_text(vjust=-0.5)) +
  theme(axis.title.y = element_text(vjust=1.5)) +
  theme(plot.title = element_text(vjust=2.5),axis.text=element_text(color = "black", size=8.75)) +
  annotate("text",-Inf, Inf, label = paste('italic(R)^2==',round(rsqr,2), sep = ""), parse = TRUE, vjust = 1.5, hjust = -0.2) +
  annotate("text",-Inf, Inf, label = paste("italic(P)",p_value), parse = TRUE, vjust = 4, hjust = -0.2) +
  theme(legend.position='none')

#grab the legend for the final figure
leg<-ggplot(v_all, aes(y = zrt*100, x = temp, size = ct)) + 
  geom_point(alpha = 1, shape = 21, fill = "red") + 
  scale_size(name = "", breaks = c(10000,50000,100000,250000), labels = c("<10,000","50,000","100,000",">250,000"))+
  scale_fill_continuous(guide = guide_legend()) +
  theme(legend.position="bottom", legend.justification = c(0.5, 0.5))


# write the stats to a table ####
write.csv(do.call(rbind, stats_ls),"output/table_v2.csv")


#   Make the final plot ####

pg<-plot_grid(vpdmax+labs(y=bquote(beta['MORTALITY-HEIGHT']) ),tmax+labs(y=""), ppt+labs(y=""),
              aws+labs(y=bquote(beta['MORTALITY-HEIGHT']) ), cover+labs(y=""), slope+labs(y=""),
              labels = c("a", "b", "c", "d", "e", "f"),
              nrow  = 2)

plot_grid(get_legend(leg), pg, nrow = 2, rel_heights = c(0.05,1))

ggsave("figures/Figure3.pdf", width = 3.75*1.5/2*3, height = 3.75*1.5, units = "in")

