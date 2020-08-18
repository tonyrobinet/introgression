## Script in R and bash (PLINK) to extract and filter SNPs from "Final_report.txt" (Illumina Bead Studio ouput) to Plink files.

## Generate a Marker_map file in which A1 and A2 alleles correspond to the ref/alt alleles indicated in the VCF ##
## The ref/alt allele are extracted from the HQ_Phased_14G_Bass_Chip_3K.vcf using the 2722 SNP list in VCFtools ##
## vcftools --vcf HQ_Phased_14G_Bass_Chip_3K.vcf --positions 2722_SNPs_list --out 14_Genomes_2722_SNPs --recode ##
## The A1 and A2 alleles columns: grep -v '#' 14_Genomes_2722_SNPs.recode.vcf | cut -f 4,5 > 2722_ref_alt_list  ##

## Morgan, A. P. (2016). argyle: an R package for analysis of Illumina genotyping arrays. G3: Genes, Genomes, Genetics, 6(2), 281-286.
## Purcell, S., Neale, B., Todd-Brown, K., Thomas, L., Ferreira, M. A., Bender, D., ... & Sham, P. C. (2007). PLINK: a tool set for whole-genome association and population-based linkage analyses. The American journal of human genetics, 81(3), 559-575.

setwd("~/Recherche/Loup/Genotypage_SNP/3K_Results")
library('argyle')
snps<-read.table('Marker_map.txt',header=TRUE)
attach(snps)
head(snps)

## Import of genotypes and hybridization intensities from BeadStudio reports ##
geno <- read.beadstudio(prefix ="", snps, in.path = "./")
summary(geno)

## Quality control ##
run.sample.qc(geno)
## Use qcplot() to generate the default QC plots ##
qcplot(geno) ## Lot of N due to uncalled loci ## 

## Identify low-performing markers based on ##
## (1) high rate of missingness                                                ##            
## (2) low or zero minor-allele frequency in the population of interest        ##
## (3) higher-than-expected heterozygosity given its minor-allele frequencies  ##

## Distribution of calls per marker ##
calls <- summarize.calls(geno, "marker")

## Distribution of missingness per marker ##
hist(calls$N,col="grey",border=NA,xlab="count of no-calls",ylab="frequency",main=NULL)
## Display 20% of max.N per locus = 846 * 0.2 = 169 ##
abline(v=169)

## Distribution of heterozygous calls per marker ##
hist(calls$H,col="grey",border=NA,xlab="count of heterozygous calls",ylab="frequency",main=NULL)
## Display 55% of max.H per locus = 846 * 0.55 = 465 ##
abline(v=465)

## Distribution of homozygous calls per marker ##
hist(with(calls,A+B),col="grey",border=NA,xlab="count of homozygous calls",ylab="frequency",main=NULL)
## Display 45% of min.hom per locus = 846 * 0.45 = 380 ##
abline(v=380)

## Set these thresholds to marker-level checks ##
geno <- run.marker.qc(geno, max.N=169, max.H=465, min.hom=380) ## Flagging 398 markers (14.6%) ##
## Apply filters to all flagged markers to exclude low-performing markers ##
geno <- apply.filters(geno)
summary(geno) ## 2324 markers remain ##

## plot the characterisics of the remaining low-performing markers with a maf < 0.01 ##
freqplot(subset(geno, chr != ""), max.N = 0.2, max.H = 0.55, min.maf = 0.01)

## plot samples QC ##
qcplot(geno)


## Set 'N' filter to all samples to flag samples with excess of no-calls ##
## 20% of max.N per individual = 2324 * 0.2 = 464 ##
geno <- run.sample.qc(geno, max.N = 465) ## Dropping 18 samples ##
## apply filter to all flagged samples to exclude low-quality samples ##
geno <- apply.filters(geno,"samples") 
summary(geno) ## 2324 markers remain in 828 samples ##

## plot remaining markers and samples to see remaining low-performing ones ##
freqplot(geno, max.N = 0.2, max.H = 0.55, min.maf = 0.01) 
## 1174 markers with a maf < 0.01 remain to be excluded in plink ##
qcplot(geno)

## Summarize genotype calls by sample ##
sum_samples <- summarize.calls(geno, by = c("samples"))
## Summarize genotype calls by locus ##
sum_markers <- summarize.calls(geno, by = c("markers"))

## Visualize cluster plots of intensity for 69 markers with plate effects (detected by Tony) ##
plot.clusters(geno,"LG1A_2658882")
plot.clusters(geno,"LG1A_11183632")
plot.clusters(geno,"LG1A_22962296")
plot.clusters(geno,"LG2_7173119")
plot.clusters(geno,"LG2_22770874")
plot.clusters(geno,"LG2_25111252")
plot.clusters(geno,"LG3_5918314")
##...
## Obviously there are miscalled loci that we want to remove from the dataset ##

## Export to plink format ##
write.plink(geno,prefix="./2324_filtered_SNPs")

##### IN PLINK APPLY MINOR ALLLELE FREQUECY THRESHOLD AND GENOTYPE CALL RATE #####
##### plink --bfile 2324_filtered_SNPs --maf 0.01 --geno 0.1 --make-bed --out maf_0.01_miss_0.1 #####
## >> maf_0.01_miss_0.1 contains 1069 SNPs ##
##### plink --bfile 2324_filtered_SNPs --maf 0.01 --geno 0.2 --make-bed --out maf_0.01_miss_0.2 #####
## >> maf_0.01_miss_0.2 contains 1114 SNPs ##
##### WE WILL USE THE MOST STRINGENT/COMPLETE DATASET WITH MAF 0.01 & MISS 0.1 #####
## Generate plink .ped .map files ##
## plink --bfile maf_0.01_miss_0.1 --recode --out 828_Chips_1069_SNPs ##

## In bash generate a list of 1069 retained SNPs to be extracted from the VCF ##
## cut -f 2 maf_0.01_miss_0.1.bim | sed 's/_/\t/g' > 1069_SNPs_list ##
## Extract from the VCF ##
## vcftools --vcf HQ_Phased_14G_Bass_Chip_3K.vcf --positions 1069_SNPs_list --out 14_Genomes_1069_SNPs --recode ##
## And Generate plink .ped .map files ##  
## vcftools --vcf 14_Genomes_1069_SNPs.recode.vcf --plink --out 14_Genomes_1069_SNPs ##

## Edit 14_Genomes_1069_SNPs.ped ##
## Replace 0 0	0	0	by 0 0	0	-9	##
## Change samples names in both ped files and concatenate them into Merged_828_Chips_14_Genomes_1069_SNPs.ped ##
## Use the 828_Chips_1069_SNPs.map file to make a Merged_828_Chips_14_Genomes_1069_SNPs.map file ##
## Finally, generate .bed .fam .bim files with plink ##
## plink --file Merged_828_Chips_14_Genomes_1069_SNPs --make-bed --out Merged_828_Chips_14_Genomes_1069_SNPs ##
## Remove duplicate individual ZDLAB0012 and generate new plink files Merged_827_Chips_14_Genomes_1069_SNPs_NoRep ##
## Generate plink .raw file named Merged_827_Chips_14_Genomes_1069_SNPs_NoRep.raw for entering adegenet in R ##
## plink --file Merged_827_Chips_14_Genomes_1069_SNPs_NoRep --recodeA --out Merged_827_Chips_14_Genomes_1069_SNPs_NoRep ##

##### DETECT LOCI WITH PLATE EFFECTS IN PLINK #####
##### A folder is created with the new plink file Merged_827_Chips_1069_SNPs_NoRep (SNP chip data only) #####
##### A pheno file with 9 columns encoding the plate location of each sample: 2 affected 1 non affected #####
##### plink --file Merged_827_Chips_1069_SNPs_NoRep --pheno pheno.txt --mpheno 1 --assoc --allow-no-sex #####
##### plink --file Merged_827_Chips_1069_SNPs_NoRep --pheno pheno.txt --mpheno 2 --assoc --allow-no-sex #####
##### plink --file Merged_827_Chips_1069_SNPs_NoRep --pheno pheno.txt --mpheno 3 --assoc --allow-no-sex #####
##### plink --file Merged_827_Chips_1069_SNPs_NoRep --pheno pheno.txt --mpheno 4 --assoc --allow-no-sex #####
##### plink --file Merged_827_Chips_1069_SNPs_NoRep --pheno pheno.txt --mpheno 5 --assoc --allow-no-sex #####
##### plink --file Merged_827_Chips_1069_SNPs_NoRep --pheno pheno.txt --mpheno 6 --assoc --allow-no-sex #####
##### plink --file Merged_827_Chips_1069_SNPs_NoRep --pheno pheno.txt --mpheno 7 --assoc --allow-no-sex #####
##### plink --file Merged_827_Chips_1069_SNPs_NoRep --pheno pheno.txt --mpheno 8 --assoc --allow-no-sex #####
##### plink --file Merged_827_Chips_1069_SNPs_NoRep --pheno pheno.txt --mpheno 9 --assoc --allow-no-sex #####
##### OUtPUT FILES ARE NAMED Plate_01_plink.assoc, Plate_02_plink.assoc, Plate_03_plink.assoc...        #####

## R ##
setwd("~/Recherche/Loup/Genotypage_SNP/3K_Results/Plink_CT_test_for_plate_effects/")
library('qqman')
assoc_1 <- read.table("Plate_01_plink.assoc",header=T)
manhattan(assoc_1)
head(assoc_1[order(assoc_1$P),],10)

assoc_2 <- read.table("Plate_02_plink.assoc",header=T)
manhattan(assoc_2)
head(assoc_2[order(assoc_2$P),],10)

assoc_3 <- read.table("Plate_03_plink.assoc",header=T)
manhattan(assoc_3)
head(assoc_3[order(assoc_3$P),],10)

assoc_4 <- read.table("Plate_04_plink.assoc",header=T)
manhattan(assoc_4)
head(assoc_4[order(assoc_4$P),],10)

assoc_5 <- read.table("Plate_05_plink.assoc",header=T)
manhattan(assoc_5)
head(assoc_5[order(assoc_5$P),],30)

assoc_6 <- read.table("Plate_06_plink.assoc",header=T)
manhattan(assoc_6)
head(assoc_6[order(assoc_6$P),],10)

assoc_7 <- read.table("Plate_07_plink.assoc",header=T)
manhattan(assoc_7)
head(assoc_7[order(assoc_7$P),],10)

assoc_8 <- read.table("Plate_08_plink.assoc",header=T)
manhattan(assoc_8)
head(assoc_8[order(assoc_8$P),],10)

assoc_9 <- read.table("Plate_09_plink.assoc",header=T)
manhattan(assoc_9)
head(assoc_9[order(assoc_9$P),],20)

assoc <- rbind(assoc_1,assoc_2,assoc_3,assoc_4,assoc_5,assoc_6,assoc_7,assoc_8,assoc_9)
SIGNIF <- assoc[assoc$P<1e-4,]
write.table(SIGNIF,"Plate_effect_loci.txt",quote=F,row.names=F,sep="\t")
write.table(sort(unique(SIGNIF$SNP)),"Loci_to_exclude.txt",quote=F,row.names=F,sep="\t")

##### IN PLINK EXCLUDE THE LIST OF 57 ARTEFACTUAL LCOI FROM THE DATASET ######
## plink --file Merged_827_Chips_14_Genomes_1069_SNPs_NoRep --exclude Loci_to_exclude --out Merged_827_Chips_14_Genomes_1012_SNPs_NoRep --make-bed ##
## plink --bfile Merged_827_Chips_14_Genomes_1012_SNPs_NoRep --recode --out Merged_827_Chips_14_Genomes_1012_SNPs_NoRep ##


## Back to Argyle ##
## Analyze the new filtered dataset ##
setwd("~/Recherche/Loup/Genotypage_SNP/3K_Results")
library(argyle)
geno.final <- read.plink("./Merged_827_Chips_14_Genomes_1012_SNPs_NoRep.bed")
summary(geno.final)

## plot remaining markers and samples to see remaining low-performing ones ##
freqplot(geno.final)

## Rough PCA analysis ##
pops=read.table("Pops_14G.txt")
PCA = pca(geno.final,K=5)
PCA$fid=pops[,3]

library(ggplot2)
## Colored PCA ##
ggplot(PCA, aes(x=PC1,y=PC2)) +
  geom_point(aes(colour = fid),size = 3)

ggplot(PCA, aes(x=PC1,y=PC3)) +
  geom_point(aes(colour = fid),size = 3)

ggplot(PCA, aes(x=PC1,y=PC4)) +
  geom_point(aes(colour = fid),size = 3)

ggplot(PCA, aes(x=PC1,y=PC5)) +
  geom_point(aes(colour = fid),size = 3)

ggplot(PCA, aes(x=PC2,y=PC3)) +
  geom_point(aes(colour = fid),size = 3)

ggplot(PCA, aes(x=PC3,y=PC4)) +
  geom_point(aes(colour = fid),size = 3)

DAT <- data.frame(PCA$fid,PCA$PC1)
ggplot(DAT, aes(DAT$PCA.fid,DAT$PCA.PC1)) + geom_boxplot()


