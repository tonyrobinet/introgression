library(MASS)
library(polynom) #
library(ggplot2)
library(cowplot) #
library(grid)
library(gridExtra) #
library(tidyr)
library(dplyr)
library(wesanderson)#
library(stringr)
library(adegenet)
library(vegan)
# mytheme <- theme_classic()
theme_set(theme_minimal())

######################
    Hobs Vs Latitude
######################

setwd("/Volumes/data/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/polymorphisme")
"polymorph_latitude_pops_827ATL1012loci.csv" %>% read.csv() %>% as_tibble() -> data_polym
data_polym %>% mutate_if(is.factor,as.numeric) -> data_polym
dim(data_polym) # [1] 20  5
head(data_polym)
attach(data_polym)
# regression logarithmique
f <- function(x,m,n) {m * log(x) + n}
fit_ln <- nls(polymorphism~f(N,m,n),start=c(m=1,n=1), model=T)
# 0.06048 * log(x) + 0.72597
summary(fit_ln)
co_ln <- coef(fit_ln)
m<-co_ln[1]
n<-co_ln[2]
# regressions linéaires
summary(lm(polymorphism~N))
summary(lm(Ec_polym_N~latitude))
summary(lm(Hobs~latitude))

a <- ggplot(data=data_polym, aes(x=latitude,y=Ec_polym_N)) + geom_point(colour="black", size=1.5) + geom_smooth(method="lm", size=.7) + theme_classic() + 
  labs(x="Latitude", y="Standardized polymorphism") + theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.position="none")

a1 <- ggplot(data=data_polym, aes(x=latitude,y=Hobs)) + geom_point(colour="black", size=1.5) + geom_smooth(method="lm", size=.7) +
  theme_classic() + labs(x="Latitude", y="Hobs") + theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.position="none")

a2 <- ggplot(data=data_polym, aes(x=N,y=polymorphism)) + geom_point(colour="black", size=1.5) + 
  geom_smooth(method="nls", colour="orange", formula=y~m*log(x)+n, se=F, method.args = list(start=c(m,n))) +
  theme_classic() + labs(x="N", y="Polymorphism") + theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.position="none")

a1 <- plot_grid(a1 , labels = "a")
a2 <- plot_grid(a2 , labels = "b")
grid.arrange(a1, a2, nrow=2)

########
    IBD
########

setwd("~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/IBD/21localities/")
dgeo_fst <- "IBD_1012loci_21sreg.csv" %>% read.csv() %>% as_tibble
dgeo_fst_2datasets <- "IBD_neutres_outliers_2datasets.csv" %>% read.csv() %>% as_tibble
dgeo_fst_sud_GONB <- "IBD_sud_GONB.csv" %>% read.csv() %>% as_tibble
dgeo_fst_nord_GONB <- "IBD_nord_GONB.csv" %>% read.csv() %>% as_tibble
combi_nord_sud_GONB <- rbind(dgeo_fst_nord_GONB,dgeo_fst_sud_GONB) %>% as_tibble()

attach(dgeo_fst)
summary(lm(FST_Rousset_1012loci ~ dist))

c <- ggplot(dgeo_fst, aes(x=dist, y=FST_Rousset_1012loci)) + geom_jitter(size=.7) + theme_classic() +
  scale_color_manual(values=wes_palette(n=2, name="Darjeeling1")) + geom_smooth(method="lm", size=.7) + expand_limits(x=c(0,3050)) + 
                  theme(axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), 
                  legend.title = element_blank(),legend.background=element_rect(fill="grey95"), 
                  legend.position="none", legend.text=element_text(size=11), axis.text.x=element_text(size=13), axis.text.y=element_text(size=13)) +
  labs(x="distance by the plateau", y="Fst / (1  - Fst)")



####################
    RDA1 Vs Latitude
####################
genind_data <- "/Volumes/data/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/fichiers_ped_gtx/21localities/FID_IID_1012loci_827ATL_21reg.gtx" %>% read.genetix()
setwd("/Volumes/data/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/RDA/")

## replace missing genotypes with mean allele freqeuncy ##
data_Nona = tab(genind_data, NA.method="mean")
write.table(data_Nona,"827_Chips_1012_Nona.txt",sep="\t",quote=F)
facteurs=read.csv("827_RDA_factors.csv",header=TRUE)
## load genotype file previiously prepared ##
data_dep=read.table("827_Chips_1012_Nona.txt",header=TRUE, sep="\t")
dep = data_dep[,2:1013]

rda_lat <- rda(dep ~ Latitude, scale=F, data=facteurs)
head(rda_lat$CCA$wa) # RDA1
cbind(rda_lat$CCA$wa,facteurs) %>% as_tibble -> y
y %>% mutate(Region_code = str_replace(Region_code, "BYNORTH", "North Biscay")) -> y
y %>% mutate(Region_code = str_replace(Region_code, "BYSOUTH", "South Biscay")) -> y
y %>% mutate(Region_code = str_replace(Region_code, "CELTIC", "Celtic sea")) -> y
y %>% mutate(Region_code = str_replace(Region_code, "CHANNEL", "Channel")) -> y
y %>% mutate(Region_code = str_replace(Region_code, "NORTH", "North sea")) -> y
y %>% mutate(Region_code = str_replace(Region_code, "PORTUGAL", "Portugal")) -> y
y %>% mutate(Region_code = str_replace(Region_code, "IRISH", "Irish sea")) -> y

levels(y$Region_code)

rda1_Vs_latitude <- ggplot(y, aes(y=RDA1,x=Latitude)) 
# rda1_Vs_latitude + geom_point(aes(colour=Region_code),size=1) + mytheme
# b <- rda1_Vs_latitude + geom_boxplot(aes(group=y$Region_code, colour=y$Region_code),size=0.5)  + theme_classic() + theme(legend.title = element_blank(),legend.background=element_rect(fill="grey95"), legend.position="none", axis.text=element_text(size=13)) + scale_fill_discrete(breaks=c("Portugal","South Biscay","North Biscay","Celtic sea","Channel","Irish sea","North sea"))

b <- rda1_Vs_latitude + geom_boxplot(aes(group=y$Region_code, colour=y$Region_code),size=0.5)  + theme_classic() + 
   theme(axis.text=element_text(size=13), legend.position = "none")


ab <- plot_grid(
    a , b
    , ncol = 2
    , align = "h", axis = c("b","t")
    , labels = "auto")

a <- plot_grid(a , labels = "a") 
b <- plot_grid(b , labels = "b") 
c <- plot_grid(c , labels = "c") 

grid.arrange(ab, c, nrow=2)

plot_grid(ab, c, ncol=1, align="v")
