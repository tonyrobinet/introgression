library(dplyr)
library(tidyr)
library(stringr)

c("Ecosse-Mer du Nord","Mer celtique","Mer d'Irlande_Pays de Galles","Manche Ouest","Gascogne France","Gascogne Espagne")-> PPOL_pool_names
c("Ecosse_mer_du_Nord","Sud Irlande","Manche ouest","Gascogne France","Gascogne Espagne","Portugal")-> MSUR_pool_names

##############
### FILTRATION mars 2018 sorties discoSnpRad sur R1
##############
# PPOL SNPs
AD_Disco_PPOL_P3_cauto_SNPs <- "/Volumes/data/sync/poissons/genetique_poissons/RAD-seq_PPOL/bio-info/discoSNPRAD/sorties_discoSnpRad/PPOL_disco_P3_c-auto_k_31.csv" %>% read.csv() %>% as_tibble() %>% filter(str_detect(ID, "SNP"))
AD_Disco_MSUR_P3_cauto_SNPs <- "/Volumes/data/sync/poissons/genetique_poissons/RAD-seq_MSUR/bio-info/discoSNPRAD/raw_MSUR_disco_P3_c-auto_k_31_c_auto_D_100_P_3_b_1_coherent_sorted_with_clusters.csv"  %>% read.csv() %>% as_tibble() %>% filter(str_detect(ID, "SNP"))
# au moins un score phred > 30
AD_Disco_PPOL_P3_cauto_SNPs_phr30 <- AD_Disco_PPOL_P3_cauto_SNPs %>% filter(.,max(ph1,ph2,ph3)>30,max(ph1.1,ph2.1,ph3.1)>30,max(ph1.2,ph2.2,ph3.2)>30,max(ph1.3,ph2.3,ph3.3)>30,max(ph1.4,ph2.4,ph3.4)>30,max(ph1.5,ph2.5,ph3.5)>30)
# 50x min dans chaque pool (juste un contrôle, car déjà paramétrisé dans stacks et Disco)
AD_Disco_PPOL_P3_cauto_SNPs_50x <- filter(AD_Disco_PPOL_P3_cauto_SNPs_phr30,G1>=50,G2>=50,G3>=50,G4>=50,G5>=50,G6>=50)
# les loci sont polymorphiques dans tous les pools
AD_Disco_PPOL_P3_cauto_SNPs_50x_polyALL <- filter(AD_Disco_PPOL_P3_cauto_SNPs_50x,G1_1>0,G1_2>0,G2_1>0,G2_2>0,G3_1>0,G3_2>0,G4_1>0,G4_2>0,G5_1>0,G5_2>0,G6_1>0,G6_2>0)
# on calcule les FA
AD_Disco_PPOL_P3_cauto_SNPs_50x_polyALL %>% transmute(.,X.CHROM=X.CHROM,G1_1p=G1_1/G1,G1_2p=G1_2/G1,G2_1p=G2_1/G2,G2_2p=G2_2/G2,G3_1p=G3_1/G3,
                                                      G3_2p=G3_2/G3,G4_1p=G4_1/G4,G4_2p=G4_2/G4,G5_1p=G5_1/G5,G5_2p=G5_2/G5,
                                                      G6_1p=G6_1/G6,G6_2p=G6_2/G6) -> PPOL_SNPs_FA
# PPOL INDELs
# au moins un score phred > 30
AD_Disco_PPOL_P3_cauto_INDELs_phr30 <- AD_Disco_PPOL_P3_cauto_INDELs %>% filter(.,max(ph1,ph2,ph3)>30,max(ph1.1,ph2.1,ph3.1)>30,max(ph1.2,ph2.2,ph3.2)>30,max(ph1.3,ph2.3,ph3.3)>30,max(ph1.4,ph2.4,ph3.4)>30,max(ph1.5,ph2.5,ph3.5)>30)
# 50x min dans chaque pool (juste un contrôle, car déjà paramétrisé dans stacks et Disco)
AD_Disco_PPOL_P3_cauto_INDELs_50x <- filter(AD_Disco_PPOL_P3_cauto_INDELs_phr30,G1>=50,G2>=50,G3>=50,G4>=50,G5>=50,G6>=50)
# les loci sont polymorphiques dans tous les pools
AD_Disco_PPOL_P3_cauto_INDELs_50x_polyALL <- filter(AD_Disco_PPOL_P3_cauto_INDELs_50x,G1_1>0,G1_2>0,G2_1>0,G2_2>0,G3_1>0,G3_2>0,G4_1>0,G4_2>0,G5_1>0,G5_2>0,G6_1>0,G6_2>0)
# on calcule les FA
AD_Disco_PPOL_P3_cauto_INDELs_50x_polyALL %>% transmute(.,X.CHROM=X.CHROM,G1_1p=G1_1/G1,G1_2p=G1_2/G1,G2_1p=G2_1/G2,G2_2p=G2_2/G2,G3_1p=G3_1/G3,
                                                      G3_2p=G3_2/G3,G4_1p=G4_1/G4,G4_2p=G4_2/G4,G5_1p=G5_1/G5,G5_2p=G5_2/G5,
                                                      G6_1p=G6_1/G6,G6_2p=G6_2/G6) -> PPOL_INDELs_FA


# MSUR SNPs
AD_Disco_PPOL_P3_cauto_INDELs <- "/Volumes/data/sync/poissons/genetique_poissons/RAD-seq_PPOL/bio-info/discoSNPRAD/sorties_discoSnpRad/PPOL_disco_P3_c-auto_k_31.csv" %>% read.csv() %>% as_tibble() %>% filter(str_detect(ID, "INDEL"))
AD_Disco_MSUR_P3_cauto_INDELs <- "/Volumes/data/sync/poissons/genetique_poissons/RAD-seq_MSUR/bio-info/discoSNPRAD/raw_MSUR_disco_P3_c-auto_k_31_c_auto_D_100_P_3_b_1_coherent_sorted_with_clusters.csv"  %>% read.csv() %>% as_tibble() %>% filter(str_detect(ID, "INDEL"))

AD_Disco_MSUR_P3_cauto_SNPs_phr30 <- AD_Disco_MSUR_P3_cauto_SNPs %>% filter(.,max(ph1,ph2,ph3)>30,max(ph1.1,ph2.1,ph3.1)>30,max(ph1.2,ph2.2,ph3.2)>30,max(ph1.3,ph2.3,ph3.3)>30,max(ph1.4,ph2.4,ph3.4)>30,max(ph1.5,ph2.5,ph3.5)>30)
AD_Disco_MSUR_P3_cauto_SNPs <- select(AD_Disco_MSUR_P3_cauto_SNPs_phr30,X.CHROM, POS, G1,G1_1,G1_2,G2,G2_1,G2_2,G3,G3_1,G3_2,G4,G4_1,G4_2,G5,G5_1,G5_2,G6,G6_1,G6_2)
AD_Disco_MSUR_P3_cauto_SNPs_50x <- filter(AD_Disco_MSUR_P3_cauto_SNPs,G1>=50,G2>=50,G3>=50,G4>=50,G5>=50,G6>=50)
AD_Disco_MSUR_P3_cauto_SNPs_50x_polyALL <- filter(AD_Disco_MSUR_P3_cauto_SNPs_50x,G1_1>0,G1_2>0,G2_1>0,G2_2>0,G3_1>0,G3_2>0,G4_1>0,G4_2>0,G5_1>0,G5_2>0,G6_1>0,G6_2>0)
AD_Disco_MSUR_P3_cauto_SNPs_50x_polyALL %>% transmute(.,X.CHROM=X.CHROM,G1_1p=G1_1/G1,G1_2p=G1_2/G1,G2_1p=G2_1/G2,G2_2p=G2_2/G2,G3_1p=G3_1/G3,
                                                      G3_2p=G3_2/G3,G4_1p=G4_1/G4,G4_2p=G4_2/G4,G5_1p=G5_1/G5,G5_2p=G5_2/G5,
                                                      G6_1p=G6_1/G6,G6_2p=G6_2/G6) -> MSUR_SNPs_FA

# MSUR INDELs
AD_Disco_MSUR_P3_cauto_INDELs_phr30 <- AD_Disco_MSUR_P3_cauto_INDELs %>% filter(.,max(ph1,ph2,ph3)>30,max(ph1.1,ph2.1,ph3.1)>30,max(ph1.2,ph2.2,ph3.2)>30,max(ph1.3,ph2.3,ph3.3)>30,max(ph1.4,ph2.4,ph3.4)>30,max(ph1.5,ph2.5,ph3.5)>30)
AD_Disco_MSUR_P3_cauto_INDELs <- select(AD_Disco_MSUR_P3_cauto_INDELs_phr30,X.CHROM, POS, G1,G1_1,G1_2,G2,G2_1,G2_2,G3,G3_1,G3_2,G4,G4_1,G4_2,G5,G5_1,G5_2,G6,G6_1,G6_2)
AD_Disco_MSUR_P3_cauto_INDELs_50x <- filter(AD_Disco_MSUR_P3_cauto_INDELs,G1>=50,G2>=50,G3>=50,G4>=50,G5>=50,G6>=50)
AD_Disco_MSUR_P3_cauto_INDELs_50x_polyALL <- filter(AD_Disco_MSUR_P3_cauto_INDELs_50x,G1_1>0,G1_2>0,G2_1>0,G2_2>0,G3_1>0,G3_2>0,G4_1>0,G4_2>0,G5_1>0,G5_2>0,G6_1>0,G6_2>0)
AD_Disco_MSUR_P3_cauto_INDELs_50x_polyALL %>% transmute(.,X.CHROM=X.CHROM,G1_1p=G1_1/G1,G1_2p=G1_2/G1,G2_1p=G2_1/G2,G2_2p=G2_2/G2,G3_1p=G3_1/G3,
                                                      G3_2p=G3_2/G3,G4_1p=G4_1/G4,G4_2p=G4_2/G4,G5_1p=G5_1/G5,G5_2p=G5_2/G5,
                                                      G6_1p=G6_1/G6,G6_2p=G6_2/G6) -> MSUR_INDELs_FA


'PPOL: nombre total de SNPs'
tally(PPOL_SNPs_FA)
'MSUR: nombre total de SNPs'
tally(MSUR_SNPs_FA)
'PPOL: nombre total de INDELs'
tally(PPOL_INDELs_FA)
'MSUR: nombre total de INDELs'
tally(MSUR_INDELs_FA)

# manque filtre: 0.1>fréquences alléliques>0.9


##############
### FILTRATION mai 2020 sorties discoSnpRad sur R1 + R2
##############
# PPOL SNPs
AD_Disco_PPOL_P1_cauto_SNPs <- "/Volumes/data/sync/poissons/genetique_poissons/RAD-seq_PPOL/bio_info_PPOL/DiscoSnpRad/2020/Galaxy76-[VCF_with_DiscoSnpRAD_on_R1R2_clonefiltered_IDuniques].vcf.csv" %>% read.csv(sep="\t") %>% as_tibble() %>% filter(str_detect(ID, "SNP"))
AD_Disco_MSUR_P1_cauto_SNPs <- "/Volumes/data/sync/poissons/genetique_poissons/RAD-seq_MSUR/bio-info_MSUR/discoSNPRAD/2018/raw_MSUR_disco_P3_c-auto_k_31_c_auto_D_100_P_3_b_1_coherent_sorted_with_clusters.csv"  %>% read.csv() %>% as_tibble() %>% filter(str_detect(ID, "SNP"))
# au moins un score phred > 30
AD_Disco_PPOL_P3_cauto_SNPs_phr30 <- AD_Disco_PPOL_P3_cauto_SNPs %>% filter(.,max(ph1,ph2,ph3)>30,max(ph1.1,ph2.1,ph3.1)>30,max(ph1.2,ph2.2,ph3.2)>30,max(ph1.3,ph2.3,ph3.3)>30,max(ph1.4,ph2.4,ph3.4)>30,max(ph1.5,ph2.5,ph3.5)>30)
# 50x min dans chaque pool (juste un contrôle, car déjà paramétrisé dans stacks et Disco)
AD_Disco_PPOL_P3_cauto_SNPs_50x <- filter(AD_Disco_PPOL_P3_cauto_SNPs_phr30,G1>=50,G2>=50,G3>=50,G4>=50,G5>=50,G6>=50)
# les loci sont polymorphiques dans tous les pools
AD_Disco_PPOL_P3_cauto_SNPs_50x_polyALL <- filter(AD_Disco_PPOL_P3_cauto_SNPs_50x,G1_1>0,G1_2>0,G2_1>0,G2_2>0,G3_1>0,G3_2>0,G4_1>0,G4_2>0,G5_1>0,G5_2>0,G6_1>0,G6_2>0)
# on calcule les FA
AD_Disco_PPOL_P3_cauto_SNPs_50x_polyALL %>% transmute(.,X.CHROM=X.CHROM,G1_1p=G1_1/G1,G1_2p=G1_2/G1,G2_1p=G2_1/G2,G2_2p=G2_2/G2,G3_1p=G3_1/G3,
                                                      G3_2p=G3_2/G3,G4_1p=G4_1/G4,G4_2p=G4_2/G4,G5_1p=G5_1/G5,G5_2p=G5_2/G5,
                                                      G6_1p=G6_1/G6,G6_2p=G6_2/G6) -> PPOL_SNPs_FA
# PPOL INDELs
# au moins un score phred > 30
AD_Disco_PPOL_P3_cauto_INDELs_phr30 <- AD_Disco_PPOL_P3_cauto_INDELs %>% filter(.,max(ph1,ph2,ph3)>30,max(ph1.1,ph2.1,ph3.1)>30,max(ph1.2,ph2.2,ph3.2)>30,max(ph1.3,ph2.3,ph3.3)>30,max(ph1.4,ph2.4,ph3.4)>30,max(ph1.5,ph2.5,ph3.5)>30)
# 50x min dans chaque pool (juste un contrôle, car déjà paramétrisé dans stacks et Disco)
AD_Disco_PPOL_P3_cauto_INDELs_50x <- filter(AD_Disco_PPOL_P3_cauto_INDELs_phr30,G1>=50,G2>=50,G3>=50,G4>=50,G5>=50,G6>=50)
# les loci sont polymorphiques dans tous les pools
AD_Disco_PPOL_P3_cauto_INDELs_50x_polyALL <- filter(AD_Disco_PPOL_P3_cauto_INDELs_50x,G1_1>0,G1_2>0,G2_1>0,G2_2>0,G3_1>0,G3_2>0,G4_1>0,G4_2>0,G5_1>0,G5_2>0,G6_1>0,G6_2>0)
# on calcule les FA
AD_Disco_PPOL_P3_cauto_INDELs_50x_polyALL %>% transmute(.,X.CHROM=X.CHROM,G1_1p=G1_1/G1,G1_2p=G1_2/G1,G2_1p=G2_1/G2,G2_2p=G2_2/G2,G3_1p=G3_1/G3,
                                                        G3_2p=G3_2/G3,G4_1p=G4_1/G4,G4_2p=G4_2/G4,G5_1p=G5_1/G5,G5_2p=G5_2/G5,
                                                        G6_1p=G6_1/G6,G6_2p=G6_2/G6) -> PPOL_INDELs_FA


# MSUR SNPs
AD_Disco_PPOL_P3_cauto_INDELs <- "/Volumes/data/sync/poissons/genetique_poissons/RAD-seq_PPOL/bio-info/discoSNPRAD/sorties_discoSnpRad/PPOL_disco_P3_c-auto_k_31.csv" %>% read.csv() %>% as_tibble() %>% filter(str_detect(ID, "INDEL"))
AD_Disco_MSUR_P3_cauto_INDELs <- "/Volumes/data/sync/poissons/genetique_poissons/RAD-seq_MSUR/bio-info/discoSNPRAD/raw_MSUR_disco_P3_c-auto_k_31_c_auto_D_100_P_3_b_1_coherent_sorted_with_clusters.csv"  %>% read.csv() %>% as_tibble() %>% filter(str_detect(ID, "INDEL"))

AD_Disco_MSUR_P3_cauto_SNPs_phr30 <- AD_Disco_MSUR_P3_cauto_SNPs %>% filter(.,max(ph1,ph2,ph3)>30,max(ph1.1,ph2.1,ph3.1)>30,max(ph1.2,ph2.2,ph3.2)>30,max(ph1.3,ph2.3,ph3.3)>30,max(ph1.4,ph2.4,ph3.4)>30,max(ph1.5,ph2.5,ph3.5)>30)
AD_Disco_MSUR_P3_cauto_SNPs <- select(AD_Disco_MSUR_P3_cauto_SNPs_phr30,X.CHROM, POS, G1,G1_1,G1_2,G2,G2_1,G2_2,G3,G3_1,G3_2,G4,G4_1,G4_2,G5,G5_1,G5_2,G6,G6_1,G6_2)
AD_Disco_MSUR_P3_cauto_SNPs_50x <- filter(AD_Disco_MSUR_P3_cauto_SNPs,G1>=50,G2>=50,G3>=50,G4>=50,G5>=50,G6>=50)
AD_Disco_MSUR_P3_cauto_SNPs_50x_polyALL <- filter(AD_Disco_MSUR_P3_cauto_SNPs_50x,G1_1>0,G1_2>0,G2_1>0,G2_2>0,G3_1>0,G3_2>0,G4_1>0,G4_2>0,G5_1>0,G5_2>0,G6_1>0,G6_2>0)
AD_Disco_MSUR_P3_cauto_SNPs_50x_polyALL %>% transmute(.,X.CHROM=X.CHROM,G1_1p=G1_1/G1,G1_2p=G1_2/G1,G2_1p=G2_1/G2,G2_2p=G2_2/G2,G3_1p=G3_1/G3,
                                                      G3_2p=G3_2/G3,G4_1p=G4_1/G4,G4_2p=G4_2/G4,G5_1p=G5_1/G5,G5_2p=G5_2/G5,
                                                      G6_1p=G6_1/G6,G6_2p=G6_2/G6) -> MSUR_SNPs_FA

# MSUR INDELs
AD_Disco_MSUR_P3_cauto_INDELs_phr30 <- AD_Disco_MSUR_P3_cauto_INDELs %>% filter(.,max(ph1,ph2,ph3)>30,max(ph1.1,ph2.1,ph3.1)>30,max(ph1.2,ph2.2,ph3.2)>30,max(ph1.3,ph2.3,ph3.3)>30,max(ph1.4,ph2.4,ph3.4)>30,max(ph1.5,ph2.5,ph3.5)>30)
AD_Disco_MSUR_P3_cauto_INDELs <- select(AD_Disco_MSUR_P3_cauto_INDELs_phr30,X.CHROM, POS, G1,G1_1,G1_2,G2,G2_1,G2_2,G3,G3_1,G3_2,G4,G4_1,G4_2,G5,G5_1,G5_2,G6,G6_1,G6_2)
AD_Disco_MSUR_P3_cauto_INDELs_50x <- filter(AD_Disco_MSUR_P3_cauto_INDELs,G1>=50,G2>=50,G3>=50,G4>=50,G5>=50,G6>=50)
AD_Disco_MSUR_P3_cauto_INDELs_50x_polyALL <- filter(AD_Disco_MSUR_P3_cauto_INDELs_50x,G1_1>0,G1_2>0,G2_1>0,G2_2>0,G3_1>0,G3_2>0,G4_1>0,G4_2>0,G5_1>0,G5_2>0,G6_1>0,G6_2>0)
AD_Disco_MSUR_P3_cauto_INDELs_50x_polyALL %>% transmute(.,X.CHROM=X.CHROM,G1_1p=G1_1/G1,G1_2p=G1_2/G1,G2_1p=G2_1/G2,G2_2p=G2_2/G2,G3_1p=G3_1/G3,
                                                        G3_2p=G3_2/G3,G4_1p=G4_1/G4,G4_2p=G4_2/G4,G5_1p=G5_1/G5,G5_2p=G5_2/G5,
                                                        G6_1p=G6_1/G6,G6_2p=G6_2/G6) -> MSUR_INDELs_FA


'PPOL: nombre total de SNPs'
tally(PPOL_SNPs_FA)
'MSUR: nombre total de SNPs'
tally(MSUR_SNPs_FA)
'PPOL: nombre total de INDELs'
tally(PPOL_INDELs_FA)
'MSUR: nombre total de INDELs'
tally(MSUR_INDELs_FA)

# manque filtre: 0.1>fréquences alléliques>0.9

devtools::install_github(repo="knausb/vcfR")
library(vcfR)
library(dplyr)

pkg <- "pinfsc50"
vcf_file <- system.file("extdata", "pinf_sc50.vcf.gz", package = pkg)
dna_file <- system.file("extdata", "pinf_sc50.fasta", package = pkg)
gff_file <- system.file("extdata", "pinf_sc50.gff", package = pkg)
vcf <- read.vcfR( vcf_file, verbose = FALSE )
vcfR2tidy(vcf)
dna <- ape::read.dna(dna_file, format = "fasta")
gff <- read.table(gff_file, sep="\t", quote="")
chrom <- create.chromR(name='Supercontig', vcf=vcf, seq=dna, ann=gff)

vcf_PPOL_R1R2 <- read.vcfR("/Volumes/data/sync/poissons/genetique_poissons/RAD-seq_PPOL/bio_info_PPOL/DiscoSnpRad/2020/Galaxy76-[VCF_with_DiscoSnpRAD_on_R1R2_clonefiltered_IDuniques].vcf.csv")
x <- vcfR2genind(vcf_PPOL_R1R2)
x$pop <- as.factor(c("G1","G2","G3","G4","G5","G6"))
vcf_PPOL_R1R2 <- addID(vcf_PPOL_R1R2)

y <- vcfR2genlight(vcf_PPOL_R1R2)
z <- vcfR2tidy(vcf_PPOL_R1R2)
fix <- z$fix %>% select(.,ID,UL,UR,REF,ALT,Ty,Rk)

# utilise script "fasta2csv.R" dans le dossier workflows pour passer le fasta en csv, puis :
fasta_PPOL_R1R2 <- read.csv("/Volumes/data/sync/poissons/genetique_poissons/RAD-seq_PPOL/bio_info_PPOL/DiscoSnpRad/2020/Galaxy77-[Fasta_with_DiscoSnpRAD_on_R1R2_clonefiltered].csv",
                            sep="\t", header=F) %>% as_tibble()
fasta_PPOL_R1R2 %>% select(.,V1,V26) %>% write.csv(., quote=F,row.names = F, 
                                                   file="/Volumes/data/sync/poissons/genetique_poissons/RAD-seq_PPOL/bio_info_PPOL/DiscoSnpRad/2020/Galaxy77-[Fasta_with_DiscoSnpRAD_on_R1R2_clonefiltered_SIMPLE].fasta.csv")
fasta_PPOL_R1R2 <- fasta_PPOL_R1R2 %>% select(.,V1,V26) %>% transmute(V1=sub('.', '', V1),V26=V26) # enlève > au début du nom des séquences du fasta ('.' signifie le 1er char)    

# pour écrire un fasta à partir d'un csv, voir https://www.biostars.org/p/11615/
# ruby -ne 'puts "" + $_.split(",").first(2).join("\n")' seqs.csv > seqs.fasta

AD <- extract.gt(vcf_PPOL_R1R2, element = "AD", as.numeric=TRUE) %>% as_tibble() %>% 
  transmute(G1_AD=G1,G2_AD=G2,G3_AD=G3,G4_AD=G4,G5_AD=G5,G6_AD=G6)
AD <- bind_cols(z$fix,AD) %>% select(.,ID,G1_AD,G2_AD,G3_AD,G4_AD,G5_AD,G6_AD)
DP <- extract.gt(vcf_PPOL_R1R2, element = "DP", as.numeric=TRUE) %>% as_tibble() %>% transmute(G1_DP=G1,G2_DP=G2,G3_DP=G3,G4_DP=G4,G5_DP=G5,G6_DP=G6)
DP <- bind_cols(z$fix,DP) %>% select(.,ID,G1_DP,G2_DP,G3_DP,G4_DP,G5_DP,G6_DP)
AD_DP <- left_join(AD,DP, by="ID")
AD_DP_fix <- left_join(AD_DP,fix, by="ID")

pop <- as.factor(c("G1", "G2", "G3", "G4", "G5", "G6"))
gendiff_PPOL_R1R2 <- genetic_diff(vcf_PPOL_R1R2, pops = pop, method = 'nei') %>% as_tibble()
gendiff_PPOL_R1R2 <- gendiff_PPOL_R1R2 %>% select(.,Ht,Gst,Htmax,Gstmax,Gprimest)
AD_DP_gendiff <- bind_cols(gendiff_PPOL_R1R2,AD_DP_fix)
# référence pour expliquer Gst (Nei) : Verity and Nichols 2014, https://pubmed.ncbi.nlm.nih.gov/25039308/
gt.to.popsum(vcf_PPOL_R1R2)
pairwise_genetic_diff(vcf_PPOL_R1R2, pop, method = "nei")

MAF_DP50x <- AD_DP_gendiff %>% filter(.,G1_DP>50,G2_DP>50,G3_DP>50,G4_DP>50,G5_DP>50,G6_DP>50) %>%
  transmute(ID=ID, Ty=Ty, Ht=Ht, Gprimest=Gprimest, G1_MAF=1-(G1_AD/G1_DP),G2_MAF=1-(G2_AD/G2_DP),G3_MAF=1-(G3_AD/G3_DP),G4_MAF=1-(G4_AD/G4_DP),G5_MAF=1-(G5_AD/G5_DP),G6_MAF=1-(G6_AD/G6_DP),
            UL=UL,UR=UR,REF=REF,ALT=ALT)

MAF0.1_DP50x <- MAF_DP50x %>% filter(.,G1_MAF>0.01,G2_MAF>0.01,G3_MAF>0.01,G4_MAF>0.01,G5_MAF>0.01,G6_MAF>0.01)
MAF0.1_DP50x_flanq30 <- MAF0.1_DP50x %>% mutate(ID= substr(ID,1,nchar(ID)-7)) %>% filter(.,UL>=30,UR>=30)
MAF0.1_DP50x_flanq50 <- MAF0.1_DP50x %>% mutate(ID= substr(ID,1,nchar(ID)-7)) %>% filter(.,UL>=50,UR>=50)

INDEL_MAF0.1_DP50x_flanq50 <- MAF0.1_DP50x_flanq50 %>% filter(.,Ty=="INS")
write.csv2(INDEL_MAF0.1_DP50x_flanq50, quote=F, 
           "/Volumes/data/sync/poissons/genetique_poissons/RAD-seq_PPOL/bio_info_PPOL/DiscoSnpRad/2020/INDEL_MAF0.1_DP50x_flanq50.csv")
write.csv2(INDEL_MAF0.1_DP50x_flanq50$ID, quote=F, row.names=F,
           "/Volumes/data/sync/poissons/genetique_poissons/RAD-seq_PPOL/bio_info_PPOL/DiscoSnpRad/2020/INDEL_MAF0.1_DP50x_flanq50_justIDs.csv")
SNP_MAF0.1_DP50x_flanq50 <- MAF0.1_DP50x_flanq50 %>% filter(.,Ty=="SNP")
write.csv2(SNP_MAF0.1_DP50x_flanq50, quote=F, 
           "/Volumes/data/sync/poissons/genetique_poissons/RAD-seq_PPOL/bio_info_PPOL/DiscoSnpRad/2020/SNP_MAF0.1_DP50x_flanq50.csv")
write.csv2(SNP_MAF0.1_DP50x_flanq50$ID, quote=F, row.names=F,
           "/Volumes/data/sync/poissons/genetique_poissons/RAD-seq_PPOL/bio_info_PPOL/DiscoSnpRad/2020/SNP_MAF0.1_DP50x_flanq50_justIDs.csv")

# il faudrait ajouter les indices Ho, He, Gst, et peut-être tester les HWE par locus...

# sur le choix de G'st: https://pubmed.ncbi.nlm.nih.gov/25039308/
# library(adegenet)
# pools <- genind2genpop(x)
# library(hierfstat)
# pools_hier <- genind2hierfstat(x, x$pop)

# fasta_PPOL_R1R2 <- ape::read.dna("/Volumes/data/sync/poissons/genetique_poissons/RAD-seq_PPOL/bio_info_PPOL/DiscoSnpRad/2020/Galaxy77-[Fasta_with_DiscoSnpRAD_on_R1R2_clonefiltered].fasta", format = "fasta")
# class(fasta_PPOL_R1R2) # [1] "DNAbin"

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")
BiocManager::install("Rsamtools")
library("Biostrings")
library("Rsamtools")
library("tidyverse")
# fasta_PPOL_R1R2 <- readDNAStringSet("/Volumes/data/sync/poissons/genetique_poissons/RAD-seq_PPOL/bio_info_PPOL/DiscoSnpRad/2020/Galaxy77-[Fasta_with_DiscoSnpRAD_on_R1R2_clonefiltered].fasta", format="fasta",
#                  nrec=-1L, skip=0L, seek.first.rec=FALSE,
#                  use.names=TRUE, with.qualities=FALSE)




