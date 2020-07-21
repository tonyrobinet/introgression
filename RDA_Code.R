
###############
##### RDA #####
###############

##### In plink, generate a .raw file #####
##### plink --file 827_Chips_1012_SNPs_NoRep --recodeA --out 827_Chips_1012_SNPs_NoRep #####
library(adegenet)
library(ggplot2)
library(gridExtra)
library(stringr)
mytheme <- theme_classic()
mytheme$axis.line.x <- mytheme$axis.line.y <- mytheme$axis.line

genind_data <- "~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/fichiers_ped_gtx/21localities/FID_IID_1012loci_827ATL_21reg.gtx" %>% read.genetix()


setwd("~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/RDA/")

## Replace missing genotypes in adegenet ##
## From the RAW plink file, generate a .txt for adegenet just for the 827 chip samples and name it 827_Chips_1012_genind.txt ##
# genotypes = read.table("827_Chips_1012_SNPs_NoRep_genind.txt", header=TRUE)
## Extract genotype data ##
# data = genotypes[1:827,2:1013]
## data to genind object ##
# genind_data = df2genind(data,sep="\t",ncode=2,ind.names=genotypes[1:827,1],ploidy=2,type="codom")

## replace missing genotypes with mean allele freqeuncy ##
data_Nona = tab(genind_data, NA.method="mean")
write.table(data_Nona,"827_Chips_1012_Nona.txt",sep="\t",quote=F)

## RDA analysis to detect plate effects ##
library(vegan)
## load factor file ##
facteurs=read.csv("827_RDA_factors.csv",header=TRUE)
## load genotype file previiously prepared ##
data_dep=read.table("827_Chips_1012_Nona.txt",header=TRUE, sep="\t")
dep = data_dep[,2:1013]


rda_all <- rda(dep ~ Latitude+Longitude, scale=F, data=facteurs)
head(summary(rda_all))
#                         RDA1   RDA2
# Proportion Explained  0.5435 0.4565
# Cumulative Proportion 0.5435 1.0000
head(rda_all$CCA$wa) # RDA1 RDA2
cbind(rda_all$CCA$wa,facteurs) %>% as_tibble() -> x
x %>% mutate(Region_code = str_replace(Region_code, "BYNORTH", "North Biscay")) -> x
x %>% mutate(Region_code = str_replace(Region_code, "BYSOUTH", "South Biscay")) -> x
x %>% mutate(Region_code = str_replace(Region_code, "CELTIC", "Celtic sea")) -> x
x %>% mutate(Region_code = str_replace(Region_code, "CHANNEL", "Channel")) -> x
x %>% mutate(Region_code = str_replace(Region_code, "NORTH", "North sea")) -> x
x %>% mutate(Region_code = str_replace(Region_code, "PORTUGAL", "Portugal")) -> x
x %>% mutate(Region_code = str_replace(Region_code, "IRISH", "Irish sea")) -> x

rda1_Vs_rda2 <- ggplot(x, aes(x=RDA1,y=RDA2))
ggplot1 <- rda1_Vs_rda2 + geom_point(aes(colour=Region_code),size=1) + 
  mytheme + theme(legend.title = element_blank(),legend.background=element_rect(fill="grey95"))

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
ggplot2 <- rda1_Vs_latitude + geom_boxplot(aes(group=y$Region_code, colour=y$Region_code),size=0.5) + 
  mytheme + theme(legend.title = element_blank(),legend.background=element_rect(fill="grey95")) + 
  scale_fill_discrete(breaks=c("Portugal","South Biscay","North Biscay","Celtic sea","Channel","Irish sea","North sea"))

grid.arrange(ggplot1,ggplot2)


# Test par permutations sur chaque variable explicative (tests marginaux)
anova(rda1,by="margin",step=100)


## Récupérer les contributions des SNPs à chacun des axes contraints après avois retiré l'effet des deux autres facteurs 
SNP_plate_RDA1 = scores(rda_plate,display="sp")[,1]
SNP_plate_RDA2 = scores(rda_plate,display="sp")[,2]

dat <- data.frame(xx = c(as.numeric(SNP_plate_RDA1),as.numeric(SNP_plate_RDA2)),yy = c(rep("RDA1",length(as.numeric(SNP_plate_RDA1))),rep("RDA2",length(as.numeric(SNP_plate_RDA2)))))


