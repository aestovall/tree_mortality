library(scales)
library(ggplot2)


#read in data
ALLtrees<-read.csv("output/ALLtrees_v2.csv")[,-1]
ALLtrees<-ALLtrees[ALLtrees$mort_year=="2016"|ALLtrees$mort_year=="Live",]


######### MORTALITY PROBABILITY ####
te<-aggregate(elevation ~ tmaxn, FUN = "mean", ALLtrees)
tmax_m<-lm(tmaxn~elevation, te)
ALLtrees$tmaxm<-predict(tmax_m, newdata = data.frame(elevation = ALLtrees$elevation))
remove(te,tmax_m,ppt_m)

ppt_m<-lm(ppt~elevation + I(cos(ALLtrees$aspect/180*pi)) + I(sin(ALLtrees$aspect/180*pi)) + tmaxn, ALLtrees)
ALLtrees$pptm<-predict(ppt_m, newdata = data.frame(elevation = ALLtrees$elevation,
                                                   `I(cos(ALLtrees$aspect/180*pi))` = cos(ALLtrees$aspect/180*pi),
                                                   `I(sin(ALLtrees$aspect/180*pi))` = sin(ALLtrees$aspect/180*pi),
                                                   tmaxn = ALLtrees$tmaxn))
remove(te,tmax_m,ppt_m)

#first the reduced VPD model
tree_fit3<-glm(dead~scale(log(zmax)) + scale(cover) + scale(vpdmax) + scale(aws) + scale(log(slope)),
               family = "binomial",
               data = ALLtrees)

# rms::vif(tree_fit3)

prob1<-exp(cbind(P = coef(tree_fit3), confint(tree_fit3)))
prob1.df<-as.data.frame(prob1)
prob1.df$var<-rownames(prob1.df)
prob1.df<-prob1.df[-1,]
prob1.df_all<-prob1.df
remove(prob1.df)
prob1.df_all$var<-factor(prob1.df_all$var, levels = prob1.df_all$var[order(abs(prob1.df_all$P))])

prob1.df_all_VPD<-prob1.df_all

#now the full model with temp and precip
tree_fit3<-glm(dead~scale(log(zmax)) + scale(cover) + scale((tmaxm-tmaxn)/tmaxn) + scale((pptm-pptn)/pptn) + scale(aws) + scale(log(slope)),
               family = "binomial",
               data = ALLtrees)

# rms::vif(tree_fit3)

prob1<-exp(cbind(P = coef(tree_fit3), confint(tree_fit3)))
prob1.df<-as.data.frame(prob1)
prob1.df$var<-rownames(prob1.df)
prob1.df<-prob1.df[-1,]
prob1.df_all<-prob1.df
remove(prob1.df)
prob1.df_all$var<-factor(prob1.df_all$var, levels = prob1.df_all$var[order(abs(prob1.df_all$P))])



# combine the two models to show relative importance of VPD
prob1.df_all
prob1.df_all1<-rbind(prob1.df_all[1,],prob1.df_all_VPD[3,],prob1.df_all[2:6,])
prob1.df_all1$var<-factor(prob1.df_all1$var, levels = prob1.df_all1$var[order(abs(prob1.df_all1$P))])


# MORTALITY RATE ####
ALLtrees_bin<-data.frame(dead = ALLtrees$dead,
                         count = ALLtrees$count,
                         zmax = 5*floor(ALLtrees$zmax/5),
                         cover = 0.1*floor(ALLtrees$cover/0.1),
                         tmaxm =  0.02*floor((ALLtrees$tmaxm-ALLtrees$tmaxn)/ALLtrees$tmaxn/0.02),
                         ppt = 0.02*floor((ALLtrees$pptm-ALLtrees$pptn)/ALLtrees$pptn/0.02),
                         slope = 5*floor(ALLtrees$slope/5),
                         aws = 2*floor(ALLtrees$aws/2),
                         vpdmax = 1*floor(ALLtrees$vpdmax/1)
)

#first the VPD model
mort<-aggregate(count~zmax + cover + vpdmax + slope + aws, FUN = "sum", ALLtrees_bin)
dead<-aggregate(dead ~zmax + cover + vpdmax + slope + aws, FUN = "sum", ALLtrees_bin)
mort$rt<-dead$dead/mort$count/2
remove(dead)
mort$rt[is.nan(mort$rt)]<-NA
mort<-mort[mort$count>10,]

tree_fit_rt<-lm(rt~scale(I(log(zmax))) + scale(cover) + scale(vpdmax) + scale(aws) + scale(slope), weights = count, data = mort)

# summary(tree_fit_rt)
# car::vif(tree_fit_rt)

rt.df<-as.data.frame(cbind(coef(tree_fit_rt),confint(tree_fit_rt)))
rt.df$var<-rownames(rt.df)
rt.df$t_value<-coef(summary(tree_fit_rt))[,"t value"]
rt.df<-rt.df[-1,]
rt.df$var<-factor(rt.df$var, levels = rt.df$var[order(prob1.df_all$P)])
colnames(rt.df)[1]<-"rt"

rt.df.vpd<-rt.df

#now the full model with temp and precip
mort<-aggregate(count~zmax + cover + tmaxm + ppt + slope + aws, FUN = "sum", ALLtrees_bin)
dead<-aggregate(dead ~zmax + cover + tmaxm + ppt + slope + aws, FUN = "sum", ALLtrees_bin)
mort$rt<-dead$dead/mort$count/2
remove(dead)
mort$rt[is.nan(mort$rt)]<-NA
mort<-mort[mort$count>10,]

tree_fit_rt<-lm(rt~scale(I(log(zmax))) + scale(cover) + scale(tmaxm) + scale(ppt) + scale(aws) + scale(slope), weights = count, data = mort)

# summary(tree_fit_rt)
# car::vif(tree_fit_rt)

rt.df<-as.data.frame(cbind(coef(tree_fit_rt),confint(tree_fit_rt)))
rt.df$var<-rownames(rt.df)
rt.df$t_value<-coef(summary(tree_fit_rt))[,"t value"]
rt.df<-rt.df[-1,]
rt.df$var<-factor(rt.df$var, levels = rt.df$var[order(prob1.df_all$P)])
colnames(rt.df)[1]<-"rt"

rt.df.all1<-rbind(rt.df[1,],rt.df.vpd[3,],rt.df[2:6,])
rt.df.all1$var <- rownames(rt.df.all1)
rt.df.all1$var<-factor(rt.df.all1$var, levels = rt.df.all1$var[order(prob1.df_all1$P)])

rt_v2<-rbind(data.frame(rt.df.all1, model = "M1"),data.frame(rt.df.vpd, model = "M2"))

###MAKE FIGURE####
prob_p<-ggplot(prob1.df_all1, aes(y = P, x = var, fill = P)) +
  scale_fill_gradient2(midpoint = 1,low = "darkblue", mid = "white", high = "red") +
  geom_hline(yintercept = 1, linetype = 2,size = 0.3)+
  scale_y_log10(limits = c(0.9, 1.4), breaks = seq(0.9,1.3, 0.1)) +
  labs(x=bquote(NULL), y=bquote('Odds Ratio (log scale)')) + 
  geom_col(width = 0.5, size = 0.1, color = "black") +
  geom_errorbar(aes(x = var, ymax = `97.5 %`, ymin = `2.5 %`), color = "black", width = 0.3, size = 0.3) +
  coord_trans(y = "log10") +
  # geom_point(data=prob1.df_all_VPD, aes(y=P, x=var), color = "black", fill = "black", shape = 4) +
  scale_x_discrete(NULL, labels = rev(c(expression(TREE[Z]),
                                        expression(Delta*'T%'),
                                        expression('*'*VPD[MAX]),
                                        expression(SOIL[AWS]),
                                        expression(COVER*'%'),
                                        expression(SLOPE),
                                        expression(Delta*'PPT%')))) +
  coord_flip() + theme_bw() + theme(panel.grid.minor = element_blank()) +
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75,color = "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(), panel.border = element_blank(),
        axis.line.y = element_line(colour = 'black', size=0.3, linetype='solid'),
        axis.line.x = element_line(colour = 'black', size=0.3, linetype='solid'),
        legend.position = "none", axis.text.y = element_text(hjust = 0))

rt_p<-ggplot(rt.df.all1, aes(y = rt*100, x = var, fill = rt)) +
  scale_fill_gradient2(midpoint = 0,low = "darkblue", mid = "white",
                       high = "red") +
  geom_hline(yintercept = 0, linetype = 2, size = 0.3) +
  labs(x=bquote(NULL), y=bquote(Delta*' Mortality Rate (% yr'^{-1}*')')) +
  geom_bar(stat = 'identity', width = 0.5, size = 0.1, color = "black") +
  geom_errorbar(aes(x = var, ymax = `97.5 %`*100, ymin = `2.5 %`*100), color = "black", width = 0.3, size = 0.3) + 
  # geom_point(data=rt.df.vpd, aes(y=rt*100, x=var), color = "black", fill = "black", shape = 4) +
  coord_flip()+
  theme_bw() + theme(axis.text.y = element_blank()) + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(color = "black", size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.line.y = element_line(colour = 'black', size=0.3, linetype='solid'),
        axis.line.x = element_blank()) + theme_bw() + theme(panel.grid.minor = element_blank()) +
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(), panel.border = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.3, linetype='solid'),
        axis.line.y = element_blank(),legend.position = "none", axis.text.y = element_blank(), axis.ticks.y=element_blank())

library(cowplot)
plot_grid(prob_p,rt_p, rel_widths = c(1.2,0.8))
ggsave("figures/Figure2.pdf", width = 3.75*1.5, height = 3.75*.6, units = "in")




#SUPPLEMENTARY ANALYSIS#########################################

prob_p<-ggplot(prob1.df_all_VPD, aes(y = P, x = var, fill = P)) +
  scale_fill_gradient2(midpoint = 1,low = "darkblue", mid = "white", high = "red") +
  geom_hline(yintercept = 1, linetype = 2,size = 0.3)+
  scale_y_log10(limits = c(0.9, 1.4), breaks = seq(0.9,1.3, 0.1)) +
  labs(x=bquote(NULL), y=bquote('Odds Ratio (log scale)')) + 
  geom_col(width = 0.5, size = 0.1, color = "black") +
  geom_errorbar(aes(x = var, ymax = `97.5 %`, ymin = `2.5 %`), color = "black", width = 0.3, size = 0.3) +
  coord_trans(y = "log10") +
  # geom_point(data=prob1.df_all_VPD, aes(y=P, x=var), color = "black", fill = "black", shape = 4) +
  scale_x_discrete(NULL, labels = rev(c(expression(TREE[Z]),
                                        expression(VPD[MAX]),
                                        expression(COVER*'%'),
                                        expression(SOIL[AWS]),
                                        expression(SLOPE)))) +
  coord_flip() + theme_bw() + theme(panel.grid.minor = element_blank()) +
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75,color = "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(), panel.border = element_blank(),
        axis.line.y = element_line(colour = 'black', size=0.3, linetype='solid'),
        axis.line.x = element_line(colour = 'black', size=0.3, linetype='solid'),
        legend.position = "none", axis.text.y = element_text(hjust = 0))

rt_p<-ggplot(rt.df.vpd, aes(y = rt*100, x = var, fill = rt)) +
  scale_fill_gradient2(midpoint = 0,low = "darkblue", mid = "white",
                       high = "red") +
  geom_hline(yintercept = 0, linetype = 2, size = 0.3) +
  labs(x=bquote(NULL), y=bquote(Delta*' Mortality Rate (% yr'^{-1}*')')) +
  geom_bar(stat = 'identity', width = 0.5, size = 0.1, color = "black") +
  geom_errorbar(aes(x = var, ymax = `97.5 %`*100, ymin = `2.5 %`*100), color = "black", width = 0.3, size = 0.3) + 
  coord_flip()+
  theme_bw() + theme(axis.text.y = element_blank()) + 
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(color = "black", size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.line.y = element_line(colour = 'black', size=0.3, linetype='solid'),
        axis.line.x = element_blank()) + theme_bw() + theme(panel.grid.minor = element_blank()) +
  theme(text = element_text(size = 8.75, colour = "black"), axis.text=element_text(size=8.75)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(), panel.border = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.3, linetype='solid'),
        axis.line.y = element_blank(),legend.position = "none", axis.text.y = element_blank(), axis.ticks.y=element_blank())


plot_grid(prob_p,rt_p, rel_widths = c(1.2,0.8))
ggsave("figures/Figure2_supplement.pdf", width = 3.75*1.5, height = 3.75*.6, units = "in")

