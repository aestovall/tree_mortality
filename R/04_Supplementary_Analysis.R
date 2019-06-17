library(scales)
library(ggplot2)


#read in data
ALLtrees<-read.csv("output/ALLtrees_v2.csv")[,-1]
ALLtrees<-ALLtrees[ALLtrees$mort_year=="2016"|ALLtrees$mort_year=="Live",]


####INTERACTION FIGURES#######
#TEMPERATURE
mort<-aggregate(count~zmax + tmaxm, FUN = "sum", ALLtrees_bin)
dead<-aggregate(dead ~zmax + tmaxm, FUN = "sum", ALLtrees_bin)
mort$rt<-dead$dead/mort$count/2
remove(dead)
mort$rt[is.nan(mort$rt)]<-NA
mort<-mort[mort$count>10,]
colnames(mort)[2]<-"var"

tree_fit_rt1<-lm(rt~zmax+zmax:var, weights = count, data = mort)
summary(tree_fit_rt1)
coef(tree_fit_rt1)

df_ls<-list()
for (i in 1:length(seq(1,65, by = 0.5))) df_ls[[i]]<-data.frame(var = seq(min(mort$var),max(mort$var),by = 0.001), zmax = seq(1,65, by = 0.5)[i])
var_int<-data.frame(do.call(rbind, df_ls),rt = predict(tree_fit_rt1, newdata = do.call(rbind, df_ls)))

tmax_p<-ggplot(var_int, aes(x = zmax, y = var*100, fill = rt*100)) + 
  geom_raster() + geom_contour(aes(z = rt*100), color = "black") + 
  scale_fill_viridis(name = "Mortality Rate",
                     option = "magma",
                     # limits=c(5, 40),
                     rescaler = function(x, to = c(0, 1), from = NULL) {
                       ifelse(x<40, 
                              scales::rescale(x,
                                              to = to,
                                              from = c(min(x, na.rm = TRUE), 40)),
                              1)}) +
  theme_bw() + geom_hline(yintercept = 0, linetype = 2)+
  labs(x=bquote('Tree Height (m)'), y=bquote(Delta*'T%')) +
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + 
  theme(legend.position = "none")

#PRECIPITATION
mort<-aggregate(count~zmax + ppt, FUN = "sum", ALLtrees_bin)
dead<-aggregate(dead ~zmax + ppt, FUN = "sum", ALLtrees_bin)
mort$rt<-dead$dead/mort$count/2
remove(dead)
mort$rt[is.nan(mort$rt)]<-NA
mort<-mort[mort$count>10,]
colnames(mort)[2]<-"var"

tree_fit_rt1<-lm(rt~zmax+zmax:var, weights = count, data = mort)
summary(tree_fit_rt1)

# tree_fit_rt1<-lm(rt~zmax + var, weights = count, data = mort)
# summary(tree_fit_rt1)
# coef(tree_fit_rt1)

df_ls<-list()
for (i in 1:length(seq(1,65, by = 0.5))) df_ls[[i]]<-data.frame(var = seq(min(mort$var),max(mort$var), by = 0.001), zmax = seq(1,65, by = 0.5)[i])
var_int<-data.frame(do.call(rbind, df_ls),rt = predict(tree_fit_rt1, newdata = do.call(rbind, df_ls)))

ppt_p<-ggplot(var_int, aes(x = zmax, y = var*100, fill = rt*100)) + 
  geom_raster() + geom_contour(aes(z = rt*100), color = "black") + 
  scale_fill_viridis(name = "Mortality Rate",
                     option = "magma",
                     # limits=c(5, 40),
                     rescaler = function(x, to = c(0, 1), from = NULL) {
                       ifelse(x<40, 
                              scales::rescale(x,
                                              to = to,
                                              from = c(min(x, na.rm = TRUE), 40)),
                              1)}) +
  theme_bw() + geom_hline(yintercept = 0, linetype = 2)+
  labs(x=bquote('Tree Height (m)'), y=bquote(Delta*'PPT%')) + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + 
  theme(legend.position = "none")

#AWS
mort<-aggregate(count~zmax + aws, FUN = "sum", ALLtrees_bin)
dead<-aggregate(dead ~zmax + aws, FUN = "sum", ALLtrees_bin)
mort$rt<-dead$dead/mort$count/2
remove(dead)
mort$rt[is.nan(mort$rt)]<-NA
mort<-mort[mort$count>10,]
colnames(mort)[2]<-"var"

tree_fit_rt1<-lm(rt~zmax+zmax:var, weights = count, data = mort)
summary(tree_fit_rt1)
coef(tree_fit_rt1)

df_ls<-list()
for (i in 1:length(seq(1,65, by = 0.5))) df_ls[[i]]<-data.frame(var = seq(min(mort$var),max(mort$var), by = 0.01), zmax = seq(1,65, by = 0.5)[i])
var_int<-data.frame(do.call(rbind, df_ls),rt = predict(tree_fit_rt1, newdata = do.call(rbind, df_ls)))

aws_p<-ggplot(var_int, aes(x = zmax, y = var, fill = rt*100)) + 
  geom_raster() + geom_contour(aes(z = rt*100), color = "black") + 
  scale_fill_viridis(name = "Mortality Rate",
                     option = "magma",
                     # limits=c(5, 40),
                     rescaler = function(x, to = c(0, 1), from = NULL) {
                       ifelse(x<40, 
                              scales::rescale(x,
                                              to = to,
                                              from = c(min(x, na.rm = TRUE), 40)),
                              1)}) +
  theme_bw() + geom_hline(yintercept = 0, linetype = 2)+
  labs(x=bquote('Tree Height (m)'), y=bquote('SOIL'['AWS'])) + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + 
  theme(legend.position = "none")

#COVER
mort<-aggregate(count~zmax + cover, FUN = "sum", ALLtrees_bin)
dead<-aggregate(dead ~zmax + cover, FUN = "sum", ALLtrees_bin)
mort$rt<-dead$dead/mort$count/2
remove(dead)
mort$rt[is.nan(mort$rt)]<-NA
mort<-mort[mort$count>10,]
colnames(mort)[2]<-"var"

tree_fit_rt1<-lm(rt~zmax+zmax:I(var^2), weights = count, data = mort)
summary(tree_fit_rt1)
coef(tree_fit_rt1)

df_ls<-list()
for (i in 1:length(seq(1,65, by = 0.5))) df_ls[[i]]<-data.frame(var = seq(min(mort$var),max(mort$var), by = 0.01), zmax = seq(1,65, by = 0.5)[i])
var_int<-data.frame(do.call(rbind, df_ls),rt = predict(tree_fit_rt1, newdata = do.call(rbind, df_ls)))

leg<-ggplot(var_int, aes(x = zmax, y = var*100, fill = rt*100)) + 
  geom_raster() + geom_contour(aes(z = rt*100), color = "black") + 
  scale_fill_viridis(name = bquote('Mortality Rate (% yr'^{-1}*')'),
                     option = "magma",
                     limits=c(5, 40)) +
  theme_bw() + geom_hline(yintercept = 0, linetype = 2)+
  labs(x=bquote('Tree Height (m)'), y=bquote('COVER%')) + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + 
  theme(legend.position = "bottom") + guides(fill = guide_colorbar(title.position = "top"))

cover_p<-ggplot(var_int, aes(x = zmax, y = var*100, fill = rt*100)) + 
  geom_raster() + geom_contour(aes(z = rt*100), color = "black") + 
  scale_fill_viridis(name = "Mortality Rate",
                     option = "magma",
                     # limits=c(5, 40),
                     rescaler = function(x, to = c(0, 1), from = NULL) {
                       ifelse(x<40, 
                              scales::rescale(x,
                                              to = to,
                                              from = c(min(x, na.rm = TRUE), 40)),
                              1)}) +
  theme_bw() + 
  labs(x=bquote('Tree Height (m)'), y=bquote('COVER%')) + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + 
  theme(legend.position = "none")

library(cowplot)
plot_grid(plot_grid(tmax_p,ppt_p,aws_p,cover_p, labels = c("a","b", "c","d")),get_legend(leg), nrow = 2, rel_heights = c(1,0.11))
ggsave("FINAL_figures/ED_ENV_interation.pdf", width = 3.75*1.8, height = 3.75*2, units = "in")

