library(MASS)
library(polynom)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(tidyr)
library(dplyr)

setwd("/Volumes/data/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/polymorphisme")
"polymorph_latitude_pops_827ATL1012loci.csv" %>% read.csv() %>% as_tibble() -> data_polym

# %polymoprhisme en fonction de N

data_polym %>% mutate_if(is.factor,as.numeric) -> data_polym
dim(data_polym) # [1] 20  5
attach(data_polym)


# regression logarithmique
f <- function(x,a,b) {a * log(x) + b}
fit_ln <- nls(polymorphism~f(N,a,b),start=c(a=1,b=1), model=T)
# 0.06048 * log(x) + 0.72597
summary(fit_ln)
co_ln <- coef(fit_ln)
a<-co_ln[1]
b<-co_ln[2]

# regressions linéaires
summary(lm(polymorphism~latitude)) # # adj r²=0.3911, p<0.01
summary(lm(Ec_polym_N~latitude)) # adj r²=0.7515, p<0.001
summary(lm(Hobs~latitude)) # adj r²=0.6815, p<0.001

haut_gauche <- ggplot(data=data_polym, aes(x=N,y=polymorphism)) + geom_point(colour="black", size=1.5) +
  geom_smooth(method="nls", colour="orange", formula=y~a*log(x)+b, se=F, method.args = list(start=c(a,b))) + theme_classic() + 
  labs(x="N (individuals)", y="Polymorphism 1012 SNP") + theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.position="none")

bas_gauche <- ggplot(data=data_polym, aes(x=latitude,y=Ec_polym_N)) + geom_point(colour="black", size=1.5) + 
  geom_smooth(aes(y=Ec_polym_N), colour="dodgerblue", method=lm) + theme_classic() + 
  labs(x="Latitude (°N)", y="Residues (Polymorphism ~ N)") + theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.position="none")

haut_droite <- ggplot(data=data_polym, aes(x=latitude,y=polymorphism)) + geom_point(colour="black", size=1.5) + 
  geom_smooth(aes(y=polymorphism), colour="dodgerblue", method=lm) + theme_classic() + 
  labs(x="Latitude (°N)", y="Polymorphism 1012 SNP") + theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.position="none")

bas_droite <- ggplot(data=data_polym, aes(x=latitude,y=Hobs)) + geom_point(colour="black", size=1.5) + 
  geom_smooth(aes(y=Hobs), colour="forestgreen", method=lm) +theme_classic() + 
  labs(x="Latitude (°N)", y="H obs") + theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.position="none")

a <- arrangeGrob(haut_gauche, top = textGrob("a", x = unit(0, "npc")
                                             , y   = unit(1, "npc"), just=c("left","top"),
                                             gp=gpar(col="black", fontsize=18, fontfamily="Times Roman")))
b <- arrangeGrob(haut_droite, top = textGrob("b", x = unit(0, "npc")
                                             , y   = unit(1, "npc"), just=c("left","top"),
                                             gp=gpar(col="black", fontsize=18, fontfamily="Times Roman")))
c <- arrangeGrob(bas_gauche, top = textGrob("c", x = unit(0, "npc")
                                            , y   = unit(1, "npc"), just=c("left","top"),
                                            gp=gpar(col="black", fontsize=18, fontfamily="Times Roman")))
d <- arrangeGrob(bas_droite, top = textGrob("d", x = unit(0, "npc")
                                            , y   = unit(1, "npc"), just=c("left","top"),
                                            gp=gpar(col="black", fontsize=18, fontfamily="Times Roman")))

grid.arrange(a,b,c,d,nrow=2, ncol=2)