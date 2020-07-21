library(ggplot2)
library(segmented)
library(labeling)
library(tidyr)
library(dplyr)
library(grid)
library(gridExtra)
mytheme <- theme_classic()
mytheme$axis.line.x <- mytheme$axis.line.y <- mytheme$axis.line

setwd("~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/admixture/modele_admixture_PAG")




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
merge(dist_SINE_lat,data_mean, by="CLST") %>% as_tibble() -> data_mean
data_mean %>% filter(dist_SINE.x>=0) %>% mutate(CLST=CLST, dist_SINE=dist_SINE.x, group=group.x, long=long, lat=lat, Q1=Q1,Q2=Q2,se=se,MIN=MIN,MAX=MAX,POS=POS) %>%
  select(CLST,dist_SINE,group,long,lat,Q1,Q2,se,MIN,MAX) -> data_mean_ATL

data_mean_ATL_sansSINE <- data_mean_ATL %>% filter(CLST!="SINE")
e1 <- ggplot(data_mean_ATL_sansSINE, aes(x=dist_SINE, y=Q2)) + geom_point(size=4)
e2 <- ggplot(data_mean_ATL, aes(x=dist_SINE, y=Q2)) + geom_point(size=4)
regression_sans_barrieres <- e2 + geom_errorbar(aes(ymin=MIN,ymax=MAX)) + scale_x_continuous(limits = c(-50, 3050)) + geom_smooth(method = "lm", se = TRUE, colour="grey10") + 
  labs(x="Distance to SINE", y = "Mean MED ancestrality", size=20) + mytheme + 
  theme(legend.position="none",axis.text=element_text(size=16), axis.title=element_text(size=16))

# sans droite ajustée
e1 + geom_errorbar(aes(ymin=MIN,ymax=MAX)) + scale_x_continuous(limits = c(-50, 3050)) +
  labs(x="Distance to SINE", y = "Mean MED ancestrality", size=20) + mytheme + 
  theme(legend.position="none",axis.text=element_text(size=16), axis.title=element_text(size=16))


# Test of the linear regression model
x <- data_mean_ATL_sansSINE$Q2
y <- data_mean_ATL_sansSINE$dist_SINE
my.lm <- lm(x~y)
summary(my.lm) # Adjusted R-squared:  0.7405 p=6.89e-07

my.seg <- segmented(my.lm, seg.Z=~y, psi= 500) # 
summary(my.seg) # Adjusted R-squared: 0.8121 p=2.39e-07
davies.test(my.lm,seg.Z=~y,k=20) # 'best' at = 2323.8, n.points = 20, p-value = 0.04357


library(strucchange)
bp <- breakpoints(x ~ y, h = 3)
breakpoints(x ~ y, data = data.frame(x, y))
plot(x ~ y, pch = 19)
lines(fitted(bp, breaks = 1) ~ x, col = 4, lwd = 1.5)
lines(fitted(bp, breaks = 2) ~ x, col = 2, lwd = 1.5)

data("down")
fit.glm<-glm(cases/births~age, weight=
               births, family=binomial, data=down)
fit.seg<-segmented(fit.glm, seg.Z=~age,
                   psi=25)
davies.test(fit.glm,"age",k=5)
plot(down$age,down$births)
