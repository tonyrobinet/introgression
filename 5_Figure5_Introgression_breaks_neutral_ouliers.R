# DLAB_distribution_Fst_MED-ATL_outliers_neutres.R

library(dplyr)
library(tidyr)
library(adegenet)
library(hierfstat)
library(ggplot2)
library(grid)
library(gridExtra)
mytheme <- theme_classic()
mytheme$axis.line.x <- mytheme$axis.line.y <- mytheme$axis.line

setwd("/Volumes/data/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/fichiers_ped_gtx/23localities/")
"noms_837DLAB_ATLMED_ICESNAME_sreg_xy.csv" %>% read.csv(.,header=T) %>% as_tibble() -> noms_827ATL_10MED_reg
geno_827ATL_10MED_52outliers_ATLMED<-read.genetix("FID_IID_52outliers_837ALTMED_23reg.gtx")
geno_827ATL_10MED_960neutres_ATLMED<-read.genetix("FID_IID_960neutral_837ALTMED_23reg.gtx")
geno_827ATL_10MED_1012loci_ATLMED<-read.genetix("FID_IID_1012loci_837ALTMED_23reg.gtx")
as.factor(noms_827ATL_10MED_reg$group) -> geno_827ATL_10MED_52outliers_ATLMED$pop -> geno_827ATL_10MED_960neutres_ATLMED$pop -> geno_827ATL_10MED_1012loci_ATLMED$pop

setwd("/Volumes/data/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/fichiers_ped_gtx/21localities/")
geno_827ATL_1012loci_ATL<-read.genetix("FID_IID_1012loci_827ATL_21reg.gtx")
geno_827ATL_32outliers_ATL<-read.genetix("FID_IID_32outliers_827ATL_21reg.gtx")
geno_827ATL_980neutral_ATL<-read.genetix("FID_IID_980neutral_827ATL_21reg.gtx")
"noms_827DLAB_ICESNAME_sreg_xy.csv" %>% read.csv(.,header=T) %>% as_tibble() -> noms_827ATL_sreg
as.factor(noms_827ATL_sreg$CLST) -> geno_827ATL_1012loci_ATL$pop -> geno_827ATL_32outliers_ATL$pop -> geno_827ATL_980neutral_ATL$pop

# FST MED-ATL by locus on neutrals and outliers with hierfstat basic.stats()
x <- genind2hierfstat(geno_827ATL_10MED_1012loci_ATLMED)
MEDATL_1012loci_stats <- basic.stats(x)
MEDATL_1012loci_stats$perloc %>% as_tibble() -> MEDATL_1012loci_tbl
x <- genind2hierfstat(geno_827ATL_1012loci_ATL)
ATL_1012loci_stats <- basic.stats(x)
ATL_1012loci_stats$perloc %>% as_tibble() -> ATL_1012loci_tbl

y <- genind2hierfstat(geno_827ATL_10MED_52outliers_ATLMED)
MEDATL_outliers_stats <- basic.stats(y)
MEDATL_outliers_stats$perloc %>% as_tibble() %>% mutate(SNP_type="outliers") -> MEDATL_outliers_tbl
z <- genind2hierfstat(geno_827ATL_10MED_960neutres_ATLMED)
MEDATL_neutres_stats <- basic.stats(z)
MEDATL_neutres_stats$perloc %>% as_tibble() %>% mutate(SNP_type="neutrals") -> MEDATL_neutres_tbl
bind_rows(MEDATL_outliers_tbl,MEDATL_neutres_tbl) -> data_bylocus_MEDATL
write.csv(data_bylocus_MEDATL,"data_bylocus_MEDATL.csv")

y <- genind2hierfstat(geno_827ATL_32outliers_ATL)
intraATL_outliers_stats <- basic.stats(y)
intraATL_outliers_stats$perloc %>% as_tibble() %>% mutate(SNP_type="outliers") -> intraATL_outliers_tbl
z <- genind2hierfstat(geno_827ATL_980neutral_ATL)
intraATL_neutres_stats <- basic.stats(z)
intraATL_neutres_stats$perloc %>% as_tibble() %>% mutate(SNP_type="neutrals") -> intraATL_neutres_tbl

# densities
distri_SNPs_Fst_MEDATL <- ggplot(data_bylocus_MEDATL, aes(x=Fst, colour=SNP_type)) + geom_density(aes(x=Fst,y=..density.., fill=SNP_type), alpha=.2) + 
  mytheme + expand_limits(x=1.00) + labs(x="Fst MED-ATL", y="scaled density") +
  theme(legend.title = element_blank(),legend.background=element_rect(fill="grey95"),
        legend.position="none", legend.text = element_text(size=27),
        axis.title.x = element_text(size=27), axis.title.y = element_text(size=27),
        axis.text.x  = element_text(size=27), axis.text.y  = element_text(size=27))

# lines
distri_SNPs_Fst_MEDATL <- ggplot(data_bylocus_MEDATL, aes(x=Fst, colour=SNP_type)) + geom_line(stat="density") + 
  mytheme + expand_limits(x=1.00) + labs(x="Fst MED-ATL", y="scaled density") +
  theme(legend.title = element_blank(),legend.background=element_rect(fill="grey95"),
        legend.position="none", legend.text = element_text(size=27),
        axis.title.x = element_text(size=27), axis.title.y = element_text(size=27),
        axis.text.x  = element_text(size=27), axis.text.y  = element_text(size=27))

# histograms
distri_SNPs_Fst_MEDATL <- ggplot(data_bylocus_MEDATL, aes(x=Fst, colour=SNP_type)) + geom_histogram(binwidth=.005, position="dodge") +
  mytheme + expand_limits(x=1.00) + labs(x="Fst MED-ATL", y="no. SNPs") +
  theme(legend.title = element_blank(),legend.background=element_rect(fill="grey95"),
        legend.position="none", legend.text = element_text(size=27),
        axis.title.x = element_text(size=27), axis.title.y = element_text(size=27),
        axis.text.x  = element_text(size=27), axis.text.y  = element_text(size=27))


# mean MAF par localité MED-ATL outliers / MED-ATL neutres
setwd("/Volumes/data/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/freq_alleliques/")
"./Freq_alleliq_827ATL_1012loci_21reg.csv" %>% read.csv(sep=",") %>% as_tibble() -> freq_827ATL_1012loci_21reg
"./dist_SINE_long_lat_21reg.csv" %>% read.csv(sep=",") %>% as_tibble() -> distSINE_21reg

freq_827ATL_1012loci_21reg <- merge(freq_827ATL_1012loci_21reg, distSINE_21reg, all=T) %>% as_tibble()
freq_827ATL_1012loci_21reg %>% group_by(CLST) %>% summarise(avg_p=mean(MAF), sd_p=sd(MAF), dist_SINE=mean(dist_SINE)) -> mean_freq_827ATL_1012loci


"/volumes/data/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/freq_alleliques/lists_neutres_outliers/list_52outliers_ATLMED_qval0.01_bayescan.txt" %>%
  read.csv(sep=" ",header=F) %>% as_tibble() -> list_52outliersATLMED
freq_827ATL_1012loci_21reg %>% filter(SNP %in% list_52outliersATLMED$V1) -> freq_827ATL_52outliersATLMED_21reg
freq_827ATL_1012loci_21reg %>% filter(!SNP %in% list_52outliersATLMED$V1) -> freq_827ATL_960neutresATLMED_21reg

meanfreq_827ATL_52outliersATLMED_21reg <- freq_827ATL_52outliersATLMED_21reg %>% group_by(CLST) %>% 
  summarise(avg_p=mean(MAF), sd_p=sd(MAF), dist_SINE=mean(dist_SINE))
meanfreq_827ATL_960neutresATLMED_21reg <- freq_827ATL_960neutresATLMED_21reg %>% group_by(CLST) %>% 
  summarise(avg_p=mean(MAF), sd_p=sd(MAF), dist_SINE=mean(dist_SINE))

meanfreq_827ATL_52outliersATLMED_21reg <- mutate(meanfreq_827ATL_52outliersATLMED_21reg,SNP_type="outliers")
meanfreq_827ATL_960neutresATLMED_21reg <- mutate(meanfreq_827ATL_960neutresATLMED_21reg,SNP_type="neutrals")

meanfreq_960neutres_52outliers_ATLMED <- bind_rows(meanfreq_827ATL_52outliersATLMED_21reg,meanfreq_827ATL_960neutresATLMED_21reg)

d52 <- meanfreq_960neutres_52outliers_ATLMED

meanMAF_MEDATLoutliers_localities <- ggplot(data=d52,aes(x=dist_SINE, y=avg_p, group=SNP_type)) + geom_point(aes(colour=d52$SNP_type, shape=d52$SNP_type), size=5) + 
  mytheme + labs(x="Distance to SINE by the plateau", y="mean MAF") + # les p chutent très bas (<0.05)
  theme(legend.title = element_blank(),legend.background=element_rect(fill="grey95"),
        legend.position="none", legend.text = element_text(size=27),
        axis.title.x = element_text(size=27), axis.title.y = element_text(size=27),
        axis.text.x  = element_text(size=27), axis.text.y  = element_text(size=27)) + 
  scale_x_continuous(limits=c(0,3050)) + scale_y_continuous(limits=c(0,0.22))


a <- arrangeGrob(distri_SNPs_Fst_MEDATL, top = textGrob("a", x = unit(0, "npc")
                                                        , y   = unit(1, "npc"), just=c("left","top"),
                                                        gp=gpar(col="black", fontsize=27, fontfamily="Times Roman")))

b <- arrangeGrob(meanMAF_MEDATLoutliers_localities, top = textGrob("b", x = unit(0, "npc")
                                                                   , y   = unit(1, "npc"), just=c("left","top"),
                                                                   gp=gpar(col="black", fontsize=27, fontfamily="Times Roman")))
grid.arrange(a,b,ncol=2)