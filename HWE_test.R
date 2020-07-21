library(adegenet)
library(hierfstat)
library(pegas)
library(parallel)
no_cores <- detectCores()
cl <- makeCluster(no_cores)
library(dplyr)
library(stringr)
# install the latest version of diveRsity
# install_github("diveRsity", "kkeenan02")

## preparation des donnees en format HIERFSTAT
## 827 ATL sur 1012 loci / 1005 neutral / 7 outliers
# read.genetix

"/Volumes/data/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/Fst/21localities/noms_827DLAB_ICESNAME_sreg_xy.csv" %>%
  read.csv(.,header=T) %>% as_tibble -> noms_827DLAB
"/Volumes/data/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/Fst/23localities/noms_837DLAB_ATLMED_ICESNAME_sreg_xy.csv" %>%
  read.csv(.,header=T) %>% as_tibble -> noms_827ATL_10MED
"/Volumes/data/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/Fst/21localities/noms_837DLAB_ICESNAME_sreg_xy.csv" %>%
  read.csv(.,header=T) %>% as_tibble -> noms_837DLAB
noms_827ATL_10MED <- right_join(noms_827DLAB,noms_827ATL_10MED, by="id_specimen")
noms_827ATL_10MED <- noms_827ATL_10MED %>% mutate(reg = if_else(CLST.y=="EMED","EMED", if_else(CLST.y=="WMED","WMED", as.character(reg)))) %>%
  mutate(CLST.x=CLST.y) %>% mutate(reg==as.factor(reg))
# noms_827ATL_10MED %>% print(n=Inf)

genind_827DLAB_1012loci <- read.genetix("/Volumes/data/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/fichiers_ped_gtx/21localities/FID_IID_1012loci_827ATL_21reg.gtx")
genind_827DLAB_980neutres <- read.genetix("/Volumes/data/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/fichiers_ped_gtx/21localities/FID_IID_980neutral_827ATL_21reg.gtx")
genind_837DLAB_1012loci <- read.genetix("/Volumes/data/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/fichiers_ped_gtx/23localities/FID_IID_1012loci_837ALTMED_23reg.gtx")

genind_827DLAB_1012loci$pop <- noms_827DLAB$reg
genind_837DLAB_1012loci$pop <- noms_827ATL_10MED$reg
genind_837DLAB_1012loci$pop <- noms_837DLAB$reg

##############################
test.between(genind_827DLAB_1012loci, test.lev = genind_827DLAB_1012loci$pop, rand.unit = , nperm = 100)
##############################
as.data.frame(genind_827DLAB_1012loci$tab) -> dat_genind
dim(dat_genind) #[1]  827 2024
dat_genind[1:20,1:10]
as.factor(levels$ssreg) -> dat_genind$Pop
boot.ppfis(dat=dat_genind,nboot=100,quant=c(0.025,0.975),diploid=TRUE,dig=4)
wc(dat_genind, diploid=T)

##############################
# sur tous loci
##############################
# HWE ATL-MED
################
hw.test(genind_837DLAB_1012loci, B=0) %>% as_tibble() -> hw_test # test HW par locus uniquement
p.adjust(hw_test$`Pr(chi^2 >)`, method = "fdr", n = length(hw_test$`Pr(chi^2 >)`)) %>% tibble::enframe(name = NULL) -> hw_test_ATLMED_p_fdr
hw_test_ATLMED_p_fdr$value.cut = cut(hw_test_ATLMED_p_fdr$value, breaks=c(20))
with(hw_test_ATLMED_p_fdr, table(value.cut, useNA='ifany')) %>% as_tibble() -> freq_occ_pval_chi2_HWE_ATLMED_p_fdr # 189 loci en HWE departure (18.68%)

seppop(genind_837DLAB_1012loci) -> temp
hw.test(temp$MED, B=0) %>% tibble::enframe(name = NULL) -> hw_test # test HW par locus uniquement
p.adjust(hw_test$`Pr(chi^2 >)`, method = "fdr", n = length(hw_test$`Pr(chi^2 >)`)) %>% as_tibble() -> hw_test_MED_p_fdr
hw_test_MED_p_fdr$value.cut = cut(hw_test_MED_p_fdr$value, breaks=c(20))
with(hw_test_MED_p_fdr, table(value.cut, useNA='ifany')) %>% as_tibble() -> freq_occ_pval_chi2_HWE_MED_p_fdr # 13 loci en HWE departure (1.28%)


# HWE ATL entier
################
hw.test(genind_827DLAB_1012loci, B=0) %>% as_tibble() -> hw_test

hw_test$value.cut = cut(hw_test$`Pr(chi^2 >)`, breaks=c(2000)) # Correction de Bonferroni tests multiples: p<0.0005 seuil pour 1000 loci, donc 2000 classes entre 0 et 1.
with(hw_test, table(value.cut, useNA='ifany')) %>% as_tibble() -> freq_occ_pval_chi2_HWE_all # 47 loci en HW departure (4.7%)
# mais la correction de Bonferroni est hyper stringente.....
# utilisons la correction FDR False Discovery Rate (cf.Waples 2014)
p.adjust(hw_test$`Pr(chi^2 >)`, method = "fdr", n = length(hw_test$`Pr(chi^2 >)`)) %>% as_tibble() -> hw_test_all_p_fdr
hw_test_all_p_fdr$value.cut = cut(hw_test_all_p_fdr$value, breaks=c(20))
with(hw_test_all_p_fdr, table(value.cut, useNA='ifany')) %>% as_tibble() -> freq_occ_pval_chi2_HWE_pop_pres_p_fdr # FDR procedure less stringent
# donc 60 loci sur 1012 sont en HW departure (en-dessous de p<0.05)
# soit une proportion de 0.0593 = 5.93% (au lieu de < 5%)
# donc à l'échelle globale, HWE n'est pas respecté (plus de 5% de loci en departure)
# le cas est similaire lorsqu'on enlève le Portugal introgressé
seppop(genind_827DLAB_1012loci) -> temp
repool(temp$Biscay,temp$Celtic_Irish,temp$Channel_North) -> genind_sansPortugal_1012loci
hw.test(genind_sansPortugal_1012loci, B=0) %>% as_tibble() -> hw_test_sansPortugal
p.adjust(hw_test_sansPortugal$`Pr(chi^2 >)`, method = "fdr", n = length(hw_test_sansPortugal$`Pr(chi^2 >)`)) %>% as_tibble() -> hw_test_sansPortugal_p_fdr
hw_test_sansPortugal_p_fdr$value.cut = cut(hw_test_sansPortugal_p_fdr$value, breaks=c(20))
with(hw_test_sansPortugal_p_fdr, table(value.cut, useNA='ifany')) %>% as_tibble() -> freq_occ_pval_chi2_HWE_sans_Portugal # 63 loci pas HWE (0.0623 = 6.23%)
# l'équilbre revient en enlevant aussi le Biscay
repool(temp$Celtic_Irish,temp$Channel_North) -> genind_sansPortugal_sansBiscay_1012loci
hw.test(genind_sansPortugal_sansBiscay_1012loci, B=0) %>% as_tibble() -> hw_test_sansPortugal_sansBiscay
p.adjust(hw_test_sansPortugal_sansBiscay$`Pr(chi^2 >)`, method = "fdr", n = length(hw_test_sansPortugal_sansBiscay$`Pr(chi^2 >)`)) %>% as_tibble() -> hw_test_sansPortugal_sansBiscay_p_fdr
hw_test_sansPortugal_sansBiscay_p_fdr$value.cut = cut(hw_test_sansPortugal_sansBiscay_p_fdr$value, breaks=c(20))
with(hw_test_sansPortugal_sansBiscay_p_fdr, table(value.cut, useNA='ifany')) %>% as_tibble() -> freq_occ_pval_chi2_HWE_sans_Portugal_sansBiscay # 42 loci pas HWE (0.0415 = 4.15%)
# North + Portugal
repool(temp$Channel_North,temp$Portugal) -> genind_North_Portugal_1012loci
hw.test(genind_North_Portugal_1012loci, B=0) %>% as_tibble() -> hw_test_North_Portugal
p.adjust(hw_test_North_Portugal$`Pr(chi^2 >)`, method = "fdr", n = length(hw_test_North_Portugal$`Pr(chi^2 >)`)) %>% as_tibble() -> hw_test_North_Portugal_p_fdr
hw_test_North_Portugal_p_fdr$value.cut = cut(hw_test_North_Portugal_p_fdr$value, breaks=c(20))
with(hw_test_North_Portugal_p_fdr, table(value.cut, useNA='ifany')) %>% as_tibble() -> freq_occ_pval_chi2_HWE_North_Portugal # 34 loci pas HWE (3.36%)


###########
# HWE par pop présumée
######################
genind_827DLAB_1012loci$pop <- noms_827DLAB$pop_pres
seppop(genind_827DLAB_1012loci) -> temp
hw.test(temp$Biscay, B=0) %>% as_tibble() -> hw_test_Biscay
hw.test(temp$Celtic_Irish, B=0) %>% as_tibble() -> hw_test_Celtic_Irish
hw.test(temp$Channel_North, B=0) %>% as_tibble() -> hw_test_Channel_North
hw.test(temp$Portugal, B=0) %>% as_tibble() -> hw_test_Portugal


p.adjust(hw_test_Biscay$`Pr(chi^2 >)`, method = "fdr", n = length(hw_test_Biscay$`Pr(chi^2 >)`)) %>% as_tibble() -> hw_test_Biscay_p_fdr
p.adjust(hw_test_Celtic_Irish$`Pr(chi^2 >)`, method = "fdr", n = length(hw_test_Celtic_Irish$`Pr(chi^2 >)`)) %>% as_tibble() -> hw_test_Celtic_Irish_p_fdr
p.adjust(hw_test_Channel_North$`Pr(chi^2 >)`, method = "fdr", n = length(hw_test_Channel_North$`Pr(chi^2 >)`)) %>% as_tibble() -> hw_test_Channel_North_p_fdr
p.adjust(hw_test_Portugal$`Pr(chi^2 >)`, method = "fdr", n = length(hw_test_Portugal$`Pr(chi^2 >)`)) %>% as_tibble() -> hw_test_Portugal_p_fdr

hw_test_Biscay_p_fdr$value.cut = cut(hw_test_Biscay_p_fdr$value, breaks=c(20)) #
hw_test_Celtic_Irish_p_fdr$value.cut = cut(hw_test_Celtic_Irish_p_fdr$value, breaks=c(20)) #
hw_test_Channel_North_p_fdr$value.cut = cut(hw_test_Channel_North_p_fdr$value, breaks=c(20))
hw_test_Portugal_p_fdr$value.cut = cut(hw_test_Portugal_p_fdr$value, breaks=c(20))

with(hw_test_Biscay_p_fdr, table(value.cut, useNA='ifany')) %>% as_tibble() -> freq_occ_pval_chi2_HWE_pop_pres_Biscay # 30 (2.96%)
with(hw_test_Celtic_Irish_p_fdr, table(value.cut, useNA='ifany')) %>% as_tibble() -> freq_occ_pval_chi2_HWE_pop_pres_Celtic_Irish # 27 (2.67%)
with(hw_test_Channel_North_p_fdr, table(value.cut, useNA='ifany')) %>% as_tibble() -> freq_occ_pval_chi2_HWE_pop_pres_Channel_North # 28 (2.77%)
with(hw_test_Portugal_p_fdr, table(value.cut, useNA='ifany')) %>% as_tibble() -> freq_occ_pval_chi2_HWE_pop_pres_Portugal # 13 (1.28%)

# à l'échelle de chaque pop présumée, HWE est respecté (moins de 5% de loci de departure)
# donc les populations sont à l'équilibre de HW, mais l'ATL global n'est pas à l'équilibre. C'est un effet Walhund.


##############################
# sur loci neutres
##############################
# HWE ATL entier
################
hw.test(genind_827DLAB_980neutres, B=0) %>% as_tibble() -> hw_test # test HW par locus

hw_test$value.cut = cut(hw_test$`Pr(chi^2 >)`, breaks=c(2000)) # Correction de Bonferroni tests multiples: p<0.0005 seuil pour 1000 loci, donc 2000 classes entre 0 et 1.
with(hw_test, table(value.cut, useNA='ifany')) %>% as_tibble() -> freq_occ_pval_chi2_HWE_all # 45 loci en HW departure

p.adjust(hw_test$`Pr(chi^2 >)`, method = "fdr", n = length(hw_test$`Pr(chi^2 >)`)) %>% as_tibble() -> hw_test_all_p_fdr
hw_test_all_p_fdr$value.cut = cut(hw_test_all_p_fdr$value, breaks=c(20))
with(hw_test_all_p_fdr, table(value.cut, useNA='ifany')) %>% as_tibble() -> freq_occ_pval_chi2_HWE_all_p_fdr # 58 (5.92%, même proportion que tous loci)


###########
# HWE par pop présumée
######################
genind_827DLAB_980neutres$pop <- noms_827DLAB$pop_pres
seppop(genind_827DLAB_980neutres) -> temp
hw.test(temp$Biscay, B=0) %>% as_tibble() -> hw_test_Biscay
hw.test(temp$Celtic_Irish, B=0) %>% as_tibble() -> hw_test_Celtic_Irish
hw.test(temp$Channel_North, B=0) %>% as_tibble() -> hw_test_Channel_North
hw.test(temp$Portugal, B=0) %>% as_tibble() -> hw_test_Portugal

p.adjust(hw_test_Biscay$`Pr(chi^2 >)`, method = "fdr", n = length(hw_test_Biscay$`Pr(chi^2 >)`)) %>% as_tibble() -> hw_test_Biscay_p_fdr
p.adjust(hw_test_Celtic_Irish$`Pr(chi^2 >)`, method = "fdr", n = length(hw_test_Celtic_Irish$`Pr(chi^2 >)`)) %>% as_tibble() -> hw_test_Celtic_Irish_p_fdr
p.adjust(hw_test_Channel_North$`Pr(chi^2 >)`, method = "fdr", n = length(hw_test_Channel_North$`Pr(chi^2 >)`)) %>% as_tibble() -> hw_test_Channel_North_p_fdr
p.adjust(hw_test_Portugal$`Pr(chi^2 >)`, method = "fdr", n = length(hw_test_Portugal$`Pr(chi^2 >)`)) %>% as_tibble() -> hw_test_Portugal_p_fdr

hw_test_Biscay_p_fdr$value.cut = cut(hw_test_Biscay_p_fdr$value, breaks=c(20)) #
hw_test_Celtic_Irish_p_fdr$value.cut = cut(hw_test_Celtic_Irish_p_fdr$value, breaks=c(20)) #
hw_test_Channel_North_p_fdr$value.cut = cut(hw_test_Channel_North_p_fdr$value, breaks=c(20))
hw_test_Portugal_p_fdr$value.cut = cut(hw_test_Portugal_p_fdr$value, breaks=c(20))
with(hw_test_Biscay_p_fdr, table(value.cut, useNA='ifany')) %>% as_tibble() -> freq_occ_pval_chi2_HWE_pop_pres_Biscay # 29
with(hw_test_Celtic_Irish_p_fdr, table(value.cut, useNA='ifany')) %>% as_tibble() -> freq_occ_pval_chi2_HWE_pop_pres_Celtic_Irish # 25
with(hw_test_Channel_North_p_fdr, table(value.cut, useNA='ifany')) %>% as_tibble() -> freq_occ_pval_chi2_HWE_pop_pres_Channel_North # 25
with(hw_test_Portugal_p_fdr, table(value.cut, useNA='ifany')) %>% as_tibble() -> freq_occ_pval_chi2_HWE_pop_pres_Portugal # 12


###########
# HWE par paire de pops     script utilisé pour calculer les HWE sur 7 régions    
#######################
genind_827DLAB_1012loci$pop <- noms_827DLAB$reg
seppop(genind_827DLAB_1012loci) -> temp
summary(temp)
repool(temp$IRISH, temp$NORTH, temp$CHANNEL, temp$CELTIC, temp$BYNORTH, temp$BYSOUTH, temp$PORTUGAL) -> x
seppop(x) -> temp
hw.test(x, B=0) %>% as_tibble() -> hw_test_all
hw.test(temp$NORTH, B=0) %>% as_tibble() -> hw_test_NORTH
hw.test(temp$IRISH, B=0) %>% as_tibble() -> hw_test_IRISH
hw.test(temp$CHANNEL, B=0) %>% as_tibble() -> hw_test_CHANNEL
hw.test(temp$CELTIC, B=0) %>% as_tibble() -> hw_test_CELTIC
hw.test(temp$BYNORTH, B=0) %>% as_tibble() -> hw_test_BYNORTH
hw.test(temp$BYSOUTH, B=0) %>% as_tibble() -> hw_test_BYSOUTH
hw.test(temp$PORTUGAL, B=0) %>% as_tibble() -> hw_test_PORTUGAL

p.adjust(hw_test_all$`Pr(chi^2 >)`, method = "fdr", n = length(hw_test_all$`Pr(chi^2 >)`)) %>% as_tibble() -> hw_test_all_p_fdr
p.adjust(hw_test_NORTH$`Pr(chi^2 >)`, method = "fdr", n = length(hw_test_NORTH$`Pr(chi^2 >)`)) %>% as_tibble() -> hw_test_NORTH_p_fdr
p.adjust(hw_test_IRISH$`Pr(chi^2 >)`, method = "fdr", n = length(hw_test_IRISH$`Pr(chi^2 >)`)) %>% as_tibble() -> hw_test_IRISH_p_fdr
p.adjust(hw_test_CHANNEL$`Pr(chi^2 >)`, method = "fdr", n = length(hw_test_CHANNEL$`Pr(chi^2 >)`)) %>% as_tibble() -> hw_test_CHANNEL_p_fdr
p.adjust(hw_test_CELTIC$`Pr(chi^2 >)`, method = "fdr", n = length(hw_test_CELTIC$`Pr(chi^2 >)`)) %>% as_tibble() -> hw_test_CELTIC_p_fdr
p.adjust(hw_test_BYNORTH$`Pr(chi^2 >)`, method = "fdr", n = length(hw_test_BYNORTH$`Pr(chi^2 >)`)) %>% as_tibble() -> hw_test_BYNORTH_p_fdr
p.adjust(hw_test_BYSOUTH$`Pr(chi^2 >)`, method = "fdr", n = length(hw_test_BYSOUTH$`Pr(chi^2 >)`)) %>% as_tibble() -> hw_test_BYSOUTH_p_fdr
p.adjust(hw_test_PORTUGAL$`Pr(chi^2 >)`, method = "fdr", n = length(hw_test_PORTUGAL$`Pr(chi^2 >)`)) %>% as_tibble() -> hw_test_PORTUGAL_p_fdr

hw_test_all_p_fdr$value.cut = cut(hw_test_all_p_fdr$value, breaks=c(20))
hw_test_NORTH_p_fdr$value.cut = cut(hw_test_NORTH_p_fdr$value, breaks=c(20))
hw_test_IRISH_p_fdr$value.cut = cut(hw_test_IRISH_p_fdr$value, breaks=c(20))
hw_test_CHANNEL_p_fdr$value.cut = cut(hw_test_CHANNEL_p_fdr$value, breaks=c(20))
hw_test_CELTIC_p_fdr$value.cut = cut(hw_test_CELTIC_p_fdr$value, breaks=c(20))
hw_test_BYNORTH_p_fdr$value.cut = cut(hw_test_BYNORTH_p_fdr$value, breaks=c(20))
hw_test_BYSOUTH_p_fdr$value.cut = cut(hw_test_BYSOUTH_p_fdr$value, breaks=c(20))
hw_test_PORTUGAL_p_fdr$value.cut = cut(hw_test_PORTUGAL_p_fdr$value, breaks=c(20))

with(hw_test_all_p_fdr, table(value.cut, useNA='ifany')) %>% as_tibble() -> freq_occ_pval_chi2_HWE_all # 60
with(hw_test_NORTH_p_fdr, table(value.cut, useNA='ifany')) %>% as_tibble() -> freq_occ_pval_chi2_HWE_NORTH # 9
with(hw_test_IRISH_p_fdr, table(value.cut, useNA='ifany')) %>% as_tibble() -> freq_occ_pval_chi2_HWE_IRISH # 9
with(hw_test_CHANNEL_p_fdr, table(value.cut, useNA='ifany')) %>% as_tibble() -> freq_occ_pval_chi2_HWE_CHANNEL # 26
with(hw_test_CELTIC_p_fdr, table(value.cut, useNA='ifany')) %>% as_tibble() -> freq_occ_pval_chi2_HWE_CELTIC # 26
with(hw_test_BYNORTH_p_fdr, table(value.cut, useNA='ifany')) %>% as_tibble() -> freq_occ_pval_chi2_HWE_BYNORTH # 47
with(hw_test_BYSOUTH_p_fdr, table(value.cut, useNA='ifany')) %>% as_tibble() -> freq_occ_pval_chi2_HWE_BYSOUTH # 18
with(hw_test_PORTUGAL_p_fdr, table(value.cut, useNA='ifany')) %>% as_tibble() -> freq_occ_pval_chi2_HWE_PORTUGAL # 27



####################################
dat1<-sim.genot(nbpop=4,nbloc=20,nbal=10,f=c(0,0.2,0.4,0.6))
boot.ppfis(dat1, diploid=F)
dat1$Pop
####################################