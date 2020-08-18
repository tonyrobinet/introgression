library(ggplot2)
library(labeling)
library(tidyr)
library(dplyr)
library(grid)
library(gridExtra)
library(ggrepel)
library(cowplot)
mytheme <- theme_classic()
mytheme$axis.line.x <- mytheme$axis.line.y <- mytheme$axis.line


"~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/admixture/noms_837DLAB_ATLMED_ICESNAME_sreg_xy.csv" %>%
  read.csv() %>% as_tibble() -> noms
"~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/admixture/dist_SINE_long_lat_23reg_VIGO_CORO.csv" %>%
  read.csv() %>% as_tibble() -> dist_SINE_lat
"~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/admixture/modele_admixture_PAG/ATL-MED_admixture_100iterations/Merged_827_Chips_10_Genomes_1012_SNPs_NoRep.2.Q" %>%
  read.csv(sep=" ",header=T) %>% as_tibble() -> data
"~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/admixture/modele_admixture_PAG/ATL-MED_admixture_100iterations/Merged_827_Chips_10_Genomes_1012_SNPs_NoRep.2.Q_se" %>%
  read.csv(sep=" ",header=F) %>% as_tibble() %>% rename(se1=V1,se2=V2) -> data_se

merge(noms,dist_SINE_lat,all=T) %>% as_tibble() %>% arrange(rank) -> noms # 837 rows
data %>% mutate(id_specimen=noms$id_specimen) %>% mutate(CLST=noms$CLST) %>% mutate(dist_SINE=noms$dist_SINE) %>% mutate(group=noms$group) %>%
  mutate(se=data_se$se2) %>% mutate(MIN=Q2-se) %>% mutate(MAX=Q2+se) %>% mutate(MIN=ifelse(MIN<0,0,MIN)) %>% mutate(MAX=ifelse(MAX>1,1,MAX)) -> data
data %>% arrange(dist_SINE,CLST) %>% mutate(POS=c(1:837))-> data

data %>% group_by(CLST) %>% summarise_all(funs(mean)) %>% as_tibble() -> data_mean
merge(dist_SINE_lat, data_mean, by="CLST") %>% as_tibble() -> data_mean
data_mean %>% filter(dist_SINE.x>=0) %>% mutate(CLST=CLST, dist_SINE=dist_SINE.x, group=group.x, long=long, lat=lat, Q1=Q1, Q2=Q2, se=se, MIN=MIN, MAX=MAX, POS=POS) -> data_mean_ATL


# noms %>% arrange(CLST) %>% print(n = Inf)
# data %>% arrange(CLST) %>% print(n = Inf)
# data_mean %>% print(n = Inf)
# data_mean_ATL %>% print(n = Inf)
# dist_SINE_lat %>% print(n = Inf)


# Make the plot
h <- ggplot(data, aes(x=POS,y=Q2)) + geom_point(aes(colour=CLST))
MED_ancestry_individual <- h + geom_errorbar(aes(ymin=MIN,ymax=MAX,colour=CLST)) +
  labs(x="Individuals within localities ranked by latitude", y = "Mediterranean ancestry") +
  mytheme + theme(legend.justification=c(0.8,1), legend.position=c(0.8,1)) +
  theme(legend.title = element_blank(),legend.background=element_rect(fill="grey95"),
        legend.text=element_text(size=11), axis.text=element_text(size=12), axis.title=element_text(size=14)) + 
  guides(colour = guide_legend(nrow = 3)) +  scale_color_discrete(breaks=c("EMED","WMED","SINE","PENI","PORT","VIGO","CORU","ASTU","CANT","AQUI","CHAR","LOIR","PBRE","GONB","CRNW","ISWI","WALE","BSEI","SGEO","DOVE","ANGL","NEUK")) 

a <- ggplot(data_mean_ATL, aes(x=dist_SINE,y=Q2, group=CLST)) + geom_point(size=4, shape=21, fill="red") + geom_text(aes(label=CLST), hjust=-1, vjust=0) + 
  geom_errorbar(aes(ymin=MIN,ymax=MAX), color="red") +
  labs(x="distance to SINE by the plateau (km)", y = "Mean MED ancestry") + mytheme + 
  theme(legend.title = element_blank(), legend.position="none", plot.background = element_rect(fill = "transparent"),
        axis.title.x = element_text(size=14), axis.title.y = element_text(size=14),
        axis.text.x  = element_text(size=14), axis.text.y  = element_text(size=14)) + expand_limits(x=3050)


# + geom_vline(xintercept=c(10.5,154.5,509.5),lty=2) +
#  annotate("text",x=-20,y=0.58,label='atop(bold("MED"))',parse=T) +
#  annotate("text",x=85,y=0.58,label='atop(bold("Portugal"))',parse=T) +
#  annotate("text",x=320,y=0.58,label='atop(bold("Biscay"))',parse=T) +
#  annotate("text",x=675,y=0.58,label='atop(bold("Celtic Sea - Channel - North Sea"))',parse=T) +
#  annotate("text",x=85,y=0.56,label=mean_POR) +
#  annotate("text",x=320,y=0.56,label=mean_BIS) +
#  annotate("text",x=675,y=0.56,label=mean_NOR) +
#  annotate("text",x=85,y=0.50,label="[0.037,0.076]") +
#  annotate("text",x=320,y=0.50,label="[0.013,0.041]") +
#  annotate("text",x=675,y=0.50,label="[0.004,0.022]")

##################################
## Simple regression of admixture vs. distance to SINE
######################################################
# avec droite ajustée
data_mean_ATL_sansSINE <- data_mean_ATL %>% filter(CLST!="SINE")
e1 <- ggplot(data_mean_ATL_sansSINE, aes(x=dist_SINE, y=Q2)) + geom_point(size=4)
e2 <- ggplot(data_mean_ATL, aes(x=dist_SINE, y=Q2)) + geom_point(size=4)
regression_sans_barrieres <- e2 + geom_errorbar(aes(ymin=MIN,ymax=MAX)) + scale_x_continuous(limits = c(-50, 3050)) + geom_smooth(method = "lm", se = TRUE, colour="grey10") + 
  labs(x="Distance to SINE", y = "Mean MED ancestry", size=20) + mytheme + 
  theme(legend.position="none",axis.text=element_text(size=16), axis.title=element_text(size=16))
attach(data_mean_ATL)
summary(lm(Q1~dist_SINE))

attach(data_mean_ATL_sansSINE)
summary(lm(Q1~dist_SINE))

# sans droite ajustée
e2 + geom_errorbar(aes(ymin=MIN,ymax=MAX)) + scale_x_continuous(limits = c(-50, 3050)) +
  labs(x="Distance to SINE", y = "Mean MED ancestry", size=20) + mytheme + 
  theme(legend.position="none",axis.text=element_text(size=16), axis.title=element_text(size=16))


# Test of the linear regression model
lin.mod <- lm(data_mean_ATL_sansSINE$Q2~data_mean_ATL_sansSINE$dist_SINE)
summary(lin.mod)

## Adjusted R-squared:  0.7405
## 18 degrees of freedom
## p=6.89e-07 ***

lin.mod <- lm(data_mean_ATL$Q2~data_mean_ATL$dist_SINE)
summary(lin.mod)
## Multiple R-squared:  0.6251, Adjusted R-squared:  0.6026
## F-statistic: 31.33 on 1 and 19 DF,  p-value: 2.132e-05
# The model is significant but quite poorly fits the data…

## Piecewise regression of admixture vs. distance to SINE
#########################################################

setwd("~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/admixture/modele_admixture_PAG/")
read.csv("mean_ancestrality_distSINE.csv")->mean_anc_dist_SINE

#############
# Fist break position
#####################

# First, we exclude the SINE location which has an increased mean Mediterranean ancestry due the presence of hybrid genotypes
mean_anc_dist_SINE_sansSINE<-mean_anc_dist_SINE[-22,] %>% mutate(POS=c(1:(length(mean_anc_dist_SINE$LOC)-1)))

x = mean_anc_dist_SINE_sansSINE$dist_SINE
y = mean_anc_dist_SINE_sansSINE$moy_anc_MED

# Now we want to test if a break in the distribution can improve the model fit
Break<-x[1:(length(x))]

# We score and plot the residual standard error obtained using different breakpoint values
d1=NULL
for (i in 1:length(Break)) 
{
  model<-lm(y~(x<Break[i])*x + (x>=Break[i])*x)
  d1[i]<-summary(model)[[6]]
}

as.data.frame(cbind(Break,d1))->data_break1

# plot(Break,d,xlab="Break position",ylab="Residual standard error",pch=16, cex=1.2)
plot(Break,d1,xlab="barrier position",ylab="Linear models residues", col="darkblue", pch=19, cex=1.2, cex.axis=1.2, cex.lab=1.5) +
  abline(v=Break[order(d1)][1],lty=2)

# The most likely break occurs at 
Break[order(d1)][1] # [1] 2078: CRNW
# And here is the model summary
piecewise_model1<-lm(y~(x<Break[order(d1)][1])*x + (x>=Break[order(d1)][1])*x)
summary(piecewise_model1)

# Multiple R-squared:  0.8843,	Adjusted R-squared:  0.8639 
# F-statistic: 43.31 on 3 and 17 DF,  p-value: 3.549e-08

#############
# Second break position
#######################
# on enlève toutes les localités au nord de GONB
mean_anc_dist_SINE_sansSINE %>% as_tibble() %>% filter(dist_SINE<2078) -> mean_anc_dist_SINE_sansSINE

x = mean_anc_dist_SINE_sansSINE$dist_SINE
y = mean_anc_dist_SINE_sansSINE$moy_anc_MED
length(x) # 12 localités

Break2<-x[1:(length(x))]
d2=NULL
for (i in 1:length(Break2)) 
{
  model<-lm(y~(x<Break2[i])*x + (x>=Break2[i])*x)
  d2[i]<-summary(model)[[6]] # residus
}

as.data.frame(cbind(Break2,d2))->data_break2

plot(Break2,d2,xlab="barrier position",ylab="Linear models residues", col="darkblue", pch=19, cex=1.2, cex.axis=1.2, cex.lab=1.5) +
  abline(v=Break2[order(d2)][1],lty=2)

# The most likely break occurs at 
Break2[order(d2)][1] # [1] 954 ASTU
# And here is the model summary
piecewise_model2<-lm(y~(x<Break2[order(d2)][1])*x + (x>=Break2[order(d2)][1])*x)
summary(piecewise_model2)

# Adjusted R-squared:  0.8164 
# F-statistic:  17.3 on 3 and 8 DF,  p-value: 0.0007397
# This time the model fits much better the data (R²=0.93)

# plot(Break,d,xlab="Break position",ylab="Residual standard error",pch=16, cex=1.2)

break_graph1 <- ggplot(data=data_break1) + geom_point(aes(x=Break,y=d1)) + mytheme + geom_vline(xintercept=2070, colour="Red") + scale_x_continuous(limits = c(0,3000)) +
  mytheme + theme(axis.text=element_text(size=16), axis.title=element_text(size=16)) + labs(y="Residual standard error", x="break position") +
  annotate("text",x=1900,y=0.00555,label='atop(bold("GONB"))',parse=T)
break_graph2 <- ggplot(data=data_break2) + geom_point(aes(x=Break2,y=d2)) + mytheme + geom_vline(xintercept=900, colour="Red") + scale_x_continuous(limits = c(0,3000)) +
  mytheme + theme(axis.text=element_text(size=16), axis.title=element_text(size=16)) + labs(y="Residual standard error", x="break position") +
  annotate("text",x=750,y=0.00645,label='atop(bold("CORU"))',parse=T)


# The most likely break occurs at 
Break[order(d2)][1] # [1] 1290: CORU
# And here is the model summary
piecewise_model2<-lm(y~(x<Break[order(d2)][1])*x + (x>=Break[order(d2)][1])*x)
summary(piecewise_model2)
# Multiple R-squared:  0.8665,	Adjusted R-squared:  0.8164 
# F-statistic:  17.3 on 3 and 8 DF,  p-value: 0.0007397

# But is it significantly better than a simple regression model excluding SINE?
# The simple regression model excluding SINE
simple_model<-lm(y~x)
summary(simple_model)

# Multiple R-squared:  0.6983,	Adjusted R-squared:  0.6682 
# F-statistic: 23.15 on 1 and 10 DF,  p-value: 0.0007115

# Compare the two models using ANOVA 
anova(piecewise_model1,simple_model,test="Chisq")
## Analysis of Variance Table
## 
## Model 1: y ~ (x < Break[order(d)][1]) * x + (x >= Break[order(d)][1]) * 
##     x
## Model 2: y ~ x
##   Res.Df        RSS Df   Sum of Sq  Pr(>Chi)    
## 1     15 0.00014299                             
## 2     17 0.00051656 -2 -0.00037357 3.094e-09 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# The reduction in the residual sum of squares in the piecewise model is statistically significant despite its two additional parameters

# So now let’s plot the model

setwd("~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/admixture/modele_admixture_PAG/")
read.csv("mean_ancestrality_distSINE.csv") %>% as_tibble() -> mean_anc_dist_SINE
mean_anc_dist_SINE_sansSINE<-mean_anc_dist_SINE[-22,]

data_mean_ATL %>% filter(CLST!="SINE") -> data_mean_ATL_sansSINE
data_mean_ATL_sansSINE %>% arrange(dist_SINE)

regression_avec_barrieres <- ggplot(data_mean_ATL, aes(x=dist_SINE, y=Q2, colour=group)) + 
  geom_errorbar(aes(ymin=MIN,ymax=MAX, color=group)) +
  geom_point(data=subset(data_mean_ATL, group="Portugal"), size=3) +
  geom_point(data=subset(data_mean_ATL, group="Biscay"), size=3) + 
  geom_point(data=subset(data_mean_ATL, group="North"), size=3) +
  geom_smooth(method="lm", se=F) + labs(x="distance to SINE by the plateau (km)", y = "Mean MED ancestry", size=20) + mytheme + 
  theme(legend.title = element_blank(),legend.background=element_rect(fill="grey95"),legend.position=c(1,0.85),legend.justification=c(1,0),
        legend.text=element_text(size=14), axis.text=element_text(size=16), axis.title=element_text(size=16)) + guides(colour = guide_legend(nrow = 3)) +
  scale_color_discrete(breaks=c("Portugal","Biscay","North")) + scale_x_continuous(limits = c(0,3000))


######

b <- arrangeGrob(break_graph1, top = textGrob("b", x = unit(0, "npc")
                                              , y   = unit(1, "npc"), just=c("left","top"),
                                              gp=gpar(col="black", fontsize=18, fontfamily="Times Roman")))
c <- arrangeGrob(break_graph2, top = textGrob("c", x = unit(0, "npc")
                                              , y   = unit(1, "npc"), just=c("left","top"),
                                              gp=gpar(col="black", fontsize=18, fontfamily="Times Roman")))
x <- arrangeGrob(regression_sans_barrieres, top = textGrob("c", x = unit(0, "npc")
                                                           , y   = unit(1, "npc"), just=c("left","top"),
                                                           gp=gpar(col="black", fontsize=18, fontfamily="Times Roman")))
d <- arrangeGrob(regression_avec_barrieres, top = textGrob("d", x = unit(0, "npc")
                                                           , y   = unit(1, "npc"), just=c("left","top"),
                                                           gp=gpar(col="black", fontsize=18, fontfamily="Times Roman")))
a <- arrangeGrob(a, top = textGrob("a", x = unit(0, "npc")
                                   , y   = unit(1, "npc"), just=c("left","top"),
                                   gp=gpar(col="black", fontsize=18, fontfamily="Times Roman")))

grid.arrange(a,b,c,d)
