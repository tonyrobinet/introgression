## 10_FigureS7_Introgression_breaks_and_genomic_divergence
##########################################################

library(adegenet)
library(poppr)
library(vegan)
library(dplyr)
library(tidyr)
library(gridExtra)
library(ggplot2)
library(ggrepel)
set.seed(99)
mytheme <- theme_classic()
mytheme$axis.line.x <- mytheme$axis.line.y <- mytheme$axis.line

# il faut extraire les loci portant l'introgression
# soit par leur Fst ATL-MED Vs. Fst ATL, soit par leur delta-p MED-ATL
# il faut regarder quel P (P1 ou P2) de la sortie admixture .2P correspondt aux freq alléliques de la pop MED
# il s'agit du P1
# les loci dont delta-p > 0.3 en faveur de la pop MED sont les loci introgressés (29 loci)
# on extrait leurs freq alleliques ATL dans les localités pour visualiser la queue d'introgression
setwd("~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/freq_alleliques/")

"./Freq_alleliq_827ATL_1012loci_21reg.csv" %>% read.csv(sep=",") %>% as_tibble() -> freq_827ATL_1012loci_21reg
"./dist_SINE_long_lat_21reg.csv" %>% read.csv(sep=",") %>% as_tibble() -> distSINE_21reg

"~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/admixture/modele_admixture_PAG/mean_ancestrality_distSINE.csv" %>% read.csv(sep=",") %>%
  as_tibble() -> mean_ancestrality_distSINE


# adjoindre les dist_SINE à chaque CLST: facile avec merge()
freq_827ATL_1012loci_21reg <- merge(freq_827ATL_1012loci_21reg, distSINE_21reg, all=T) %>% as_tibble()

freq_827ATL_1012loci_21reg %>% group_by(CLST) %>% summarise(avg_p=mean(MAF), sd_p=sd(MAF), dist_SINE=mean(dist_SINE)) -> mean_freq_827ATL_1012loci


# filtrer les locus introgressés uniquement
"~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/freq_alleliques/lists_neutres_outliers/list_7outliers_ATL_qval0.05_bayescan.txt" %>%
  read.csv(sep=" ",header=F) %>% as_tibble() -> list_7outliersATL
freq_827ATL_1012loci_21reg %>% filter(SNP %in% list_7outliersATL$V1) -> freq_827ATL_7outliersATL_21reg
freq_827ATL_1012loci_21reg %>% filter(!SNP %in% list_7outliersATL$V1) -> freq_827ATL_1005neutralATL_21reg

"~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/freq_alleliques/lists_neutres_outliers/list_52outliers_ATLMED_qval0.01_bayescan.txt" %>%
  read.csv(sep=" ",header=F) %>% as_tibble() -> list_52outliersATLMED
freq_827ATL_1012loci_21reg %>% filter(SNP %in% list_52outliersATLMED$V1) -> freq_827ATL_52outliersATLMED_21reg
freq_827ATL_1012loci_21reg %>% filter(!SNP %in% list_52outliersATLMED$V1) -> freq_827ATL_960neutresATLMED_21reg

"Freq_alleliq_761ATL_10MED_1012loci_2CLST.csv" %>% read.csv(sep=",") %>% as_tibble() -> freq_827ATL_10MED_1012loci_2CLST
freq_827ATL_10MED_1012loci_2CLST %>% filter(SNP %in% list_52outliersATLMED$V1) %>% mutate(SNP_type="outliers") -> freq_827ATL_10MED_52outliers_2CLST
freq_827ATL_10MED_1012loci_2CLST %>% filter(!SNP %in% list_52outliersATLMED$V1) %>% mutate(SNP_type="neutres") -> freq_827ATL_10MED_960neutres_2CLST

bind_rows(freq_827ATL_10MED_52outliers_2CLST,freq_827ATL_10MED_960neutres_2CLST) %>% arrange(desc(MAF))-> freq_827ATL_10MED_2CLST
freq_827ATL_10MED_2CLST %>% group_by(SNP) %>% summarise(var=var(MAF),moy=mean(MAF), Fst= var*(1-moy)) -> Fst_ATLMED_1012loci
Fst_ATLMED_1012loci %>% filter(SNP %in% list_52outliersATLMED$V1) %>% arrange(desc(Fst)) %>% mutate(SNP_type="outliers") -> Fst_ATLMED_52outliers
Fst_ATLMED_1012loci %>% filter(!SNP %in% list_52outliersATLMED$V1) %>% arrange(desc(Fst)) %>% mutate(SNP_type="neutres") -> Fst_ATLMED_960neutres
bind_rows(Fst_ATLMED_52outliers,Fst_ATLMED_960neutres) -> Fst_ATLMED_52outliers_960neutres

cut(Fst_ATLMED_52outliers_960neutres$Fst, seq(from = 0, to = 1, by = 0.05)) -> Fst_ATLMED_52outliers_960neutres$cut
with(Fst_ATLMED_52outliers_960neutres, table(cut)) %>% as_tibble() -> freq_occ_Fst_ATLMED_52outliers_960neutres

"./Merged_827_Chips_14_Genomes_1012_SNPs_NoRep.2.P.csv" %>% read.csv() %>% as_tibble -> res_admix_MED_ATL_2P

cut(res_admix_MED_ATL_2P$Delta.p, seq(from = 0, to = 1, by = 0.05)) -> res_admix_MED_ATL_2P$value.cut
with(res_admix_MED_ATL_2P, table(value.cut, useNA='ifany')) %>% as_tibble() -> freq_occ_delta_p_admix_MED_ATL
## Plot the distribution of 1012 delta-p MED-ATL after 100 Admixture longruns
ggplot(freq_occ_delta_p_admix_MED_ATL,aes(x=value.cut,y=n)) + geom_histogram(stat = "identity", binwidth = 2) + mytheme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Distribution of 1012 delta-p MED-ATL 
        after 100 Admixture longruns")

library(data.table)
res_admix_MED_ATL_2P %>%
  filter(Delta.p<0.1) -> MED_ATL_0.0_0.1 # 233 locus
res_admix_MED_ATL_2P %>%
  filter(Delta.p>0.2) %>% filter(Delta.p<0.3) -> MED_ATL_0.2_0.3 # 150 locus
res_admix_MED_ATL_2P %>%
  filter(Delta.p>0.4) %>% filter(Delta.p<0.5) -> MED_ATL_0.4_0.5 # 77 locus
res_admix_MED_ATL_2P %>%
  filter(Delta.p>0.6) %>% filter(Delta.p<0.7) -> MED_ATL_0.6_0.7 # 58 locus
res_admix_MED_ATL_2P %>%
  filter(Delta.p>0.8) %>% filter(Delta.p<0.9) -> MED_ATL_0.8_0.9 # 35 locus
res_admix_MED_ATL_2P %>%
  filter(Delta.p>0.9) -> MED_ATL_0.9_1.0 # 42 locus

# attention à prendre les bonnes données (20 ou 21 reg)
freq_827ATL_1012loci_21reg -> freq_827ATL_1012loci

freq_827ATL_1012loci %>% filter(SNP %in% MED_ATL_0.0_0.1$POS) %>% group_by(CLST) %>% summarise(avg_p=mean(MAF), sd_p=sd(MAF), dist_SINE=mean(dist_SINE)) %>% mutate(range_delta_p="0.0-0.1") -> mean_freq_827ATL_0.0_0.1
freq_827ATL_1012loci %>% filter(SNP %in% MED_ATL_0.2_0.3$POS) %>% group_by(CLST) %>% summarise(avg_p=mean(MAF), sd_p=sd(MAF), dist_SINE=mean(dist_SINE)) %>% mutate(range_delta_p="0.2-0.3") -> mean_freq_827ATL_0.2_0.3
freq_827ATL_1012loci %>% filter(SNP %in% MED_ATL_0.4_0.5$POS) %>% group_by(CLST) %>% summarise(avg_p=mean(MAF), sd_p=sd(MAF), dist_SINE=mean(dist_SINE)) %>% mutate(range_delta_p="0.4-0.5") -> mean_freq_827ATL_0.4_0.5
freq_827ATL_1012loci %>% filter(SNP %in% MED_ATL_0.6_0.7$POS) %>% group_by(CLST) %>% summarise(avg_p=mean(MAF), sd_p=sd(MAF), dist_SINE=mean(dist_SINE)) %>% mutate(range_delta_p="0.6-0.7") -> mean_freq_827ATL_0.6_0.7
freq_827ATL_1012loci %>% filter(SNP %in% MED_ATL_0.8_0.9$POS) %>% group_by(CLST) %>% summarise(avg_p=mean(MAF), sd_p=sd(MAF), dist_SINE=mean(dist_SINE)) %>% mutate(range_delta_p="0.8-0.9") -> mean_freq_827ATL_0.8_0.9
freq_827ATL_1012loci %>% filter(SNP %in% MED_ATL_0.9_1.0$POS) %>% group_by(CLST) %>% summarise(avg_p=mean(MAF), sd_p=sd(MAF), dist_SINE=mean(dist_SINE)) %>% mutate(range_delta_p="0.9-1.0") -> mean_freq_827ATL_0.9_1.0
data <- bind_rows(mean_freq_827ATL_0.0_0.1,mean_freq_827ATL_0.2_0.3,mean_freq_827ATL_0.4_0.5,mean_freq_827ATL_0.6_0.7,mean_freq_827ATL_0.8_0.9,mean_freq_827ATL_0.9_1.0)

# reprendre les 5 SNPs outliers bayescan intra-ATL qui sont aussi outliers MED-ATL, et représenter l'évolution de leur fréquence allélique respective p/p à la distance à SINE
ggplot(data=freq_827ATL_5outliersATLMED_21reg, aes(x=dist_SINE, y=avg_p)) + geom_line()


###################################
# introgression profiles Vs delta-p
###################################

ggplot(data, aes(x=dist_SINE,y=avg_p, group=range_delta_p)) + geom_line(aes(colour=data$range_delta_p)) + mytheme +
  theme(legend.position="right", legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=16),
        axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),
        axis.text.x  = element_text(size=16), axis.text.y  = element_text(size=16)) +
  scale_y_continuous(limits = c(0, 0.25)) + 
  labs(x="distance to SINE", y="mean MAF") + scale_color_discrete(expression(paste("Range of 
ATL-MED ",Delta,"-p")), labels=c("0.0-0.1 (233 loci)", "0.2-0.3 (150 loci)", "0.4-0.5 (77 loci)", 
                                 "0.6-0.7 (58 loci)", "0.8-0.9 (35 loci)", "0.9-1.0 (42 loci)"))
