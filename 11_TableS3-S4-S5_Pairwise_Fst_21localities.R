library(adegenet)
library(hierfstat)
library(pegas)
library(parallel)
library(dplyr)
no_cores <- detectCores()
cl <- makeCluster(no_cores)

## preparation des donnees en format HIERFSTAT
## 827 ATL sur 1012 loci / 1005 neutral / 7 outliers
# read.genetix


"~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/Fst/21localities/noms_827DLAB_ICESNAME_sreg_xy.csv" %>% read.csv(.,header=T) %>% as_tibble -> noms_827DLAB
"~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/fichiers_ped_gtx/23localities/noms_837DLAB_ATLMED_ICESNAME_sreg_xy.csv" %>% 
  read.csv(sep=",",header=T) %>% as_tibble -> noms_827ATL_10MED

setwd("~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/fichiers_ped_gtx/21localities/")
geno_827ATL_1012loci_21sreg <- read.genetix("FID_IID_1012loci_827ATL_21reg.gtx")
as.factor(noms_827DLAB$CLST) -> geno_827ATL_1012loci_21sreg$pop
geno_827ATL_32outliers_21sreg <- read.genetix("FID_IID_32outliers_827ATL_21reg.gtx")
as.factor(noms_827DLAB$CLST) -> geno_827ATL_32outliers_21sreg$pop
geno_827ATL_980neutral_21sreg <- read.genetix("FID_IID_980neutral_827ATL_21reg.gtx")
as.factor(noms_827DLAB$CLST) -> geno_827ATL_980neutral_21sreg$pop

geno_827ATL_10MED_52outliers_9reg <- read.genetix("~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/fichiers_ped_gtx/23localities/FID_IID_52outliers_837ALTMED_23reg.gtx")
as.factor(noms_827ATL_10MED$reg) -> geno_827ATL_10MED_52outliers_9reg$pop

geno_827ATL_32outliers_21sreg <- read.genetix("~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/fichiers_ped_gtx/21localities/FID_IID_32outliers_827ATL_21reg.gtx")
as.factor(noms_827DLAB$reg) -> geno_827ATL_32outliers_21sreg$pop


# Global Fst on all loci
x=pvalues=mat.obs=mat.perm=NPERM=allTests=NULL
x <- geno_827ATL_1012loci_21sreg
mat.obs <- pairwise.fst(x, res.type="matrix") # compute original Fst matrix
NBPERM <- 5000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM,function(i) pairwise.fst(x, pop=sample(pop(x)), res.type="matrix"), mc.cores=no_cores)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) mat.perm[[k]][i,j])), mat.obs[i,j], alter="greater")$pvalue
  }
}
data.frame(p=unlist(allTests))->pvalues
data.frame(p=unlist(allTests_827ATL_1012loci_21sreg))->pvalues_827ATL_1012loci_21sreg
write.table(pvalues,row.names=T,col.names=T,sep="\t",quote=F,"~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/Fst/21localities/pvalues_827ATL_1012loci_21sreg.txt")
as.dist(mat.obs, upper=T)
write.csv(mat.obs, "Fst_827ATL_1012loci_21sreg.csv")

# Fst on outliers 
x <- geno_827ATL_10MED_52outliers_9reg
mat.obs <- wc(x, diploid = TRUE)
NBPERM <- 5000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) wc(x,  diploid = TRUE), mc.cores=no_cores)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) mat.perm[[k]][i,j])), mat.obs[i,j], alter="greater")$pvalue
  }
}
allTests->allTests_827ATL_10MED_52outliers_9reg
data.frame(p=unlist(allTests_827ATL_1012loci_21sreg))->pvalues_827ATL_1012loci_21sreg
setwd("~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/Fst/21localities/")
write.table(pvalues_827ATL_1012loci_21sreg,row.names=T,col.names=T,sep="\t",quote=F,"pvalues_827ATL_1012loci_21sreg.txt")
as.dist(mat.obs)->Fst_827ATL_1012loci_21sreg
write.csv(mat.obs, "Fst_827ATL_1012loci_21sreg.csv")

gstat.randtest(geno_827ATL_1012loci_21sreg,pop=geno_827ATL_1012loci_21sreg$pop,method = "global", nsim=500)


setwd("~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/fichiers_ped_gtx/23localities/")
"noms_837DLAB_ATLMED_ICESNAME_sreg_xy.csv" %>% read.csv(.,header=T) %>% as_tibble() -> noms_827ATL_10MED_reg
geno_827ATL_10MED_183outliers_23sreg<-read.genetix("FID_IID_183outliers_837ALTMED_23reg.gtx")
noms_827ATL_10MED_reg$CLST->geno_827ATL_10MED_183outliers_23sreg$pop
geno_827ATL_10MED_829neutral_23sreg<-read.genetix("FID_IID_829neutral_837ALTMED_23reg.gtx")
noms_827ATL_10MED_reg$CLST->geno_827ATL_10MED_829neutral_23sreg$pop
geno_827ATL_10MED_52outliers_ATLMED<-read.genetix("FID_IID_52outliers_837ALTMED_23reg.gtx")
geno_827ATL_10MED_960neutres_ATLMED<-read.genetix("FID_IID_960neutral_837ALTMED_23reg.gtx")
noms_827ATL_10MED_reg$group -> geno_827ATL_10MED_52outliers_ATLMED$pop -> geno_827ATL_10MED_960neutres_ATLMED$pop

# FST MED-ATL by locus on neutrals and outliers
y <- genind2hierfstat(geno_827ATL_10MED_52outliers_ATLMED)
MEDATL_outliers_stats <- basic.stats(y)
MEDATL_outliers_stats$perloc %>% as_tibble() %>% mutate(SNP_type="outliers") -> MEDATL_outliers_tbl
z <- genind2hierfstat(geno_827ATL_10MED_960neutres_ATLMED)
MEDATL_neutres_stats <- basic.stats(z)
MEDATL_neutres_stats$perloc %>% as_tibble() %>% mutate(SNP_type="neutrals") -> MEDATL_neutres_tbl
bind_rows(MEDATL_outliers_tbl,MEDATL_neutres_tbl) -> data_bylocus_MEDATL
cut(data_bylocus_MEDATL$Fst, seq(from = 0, to = 1, by = 0.05)) -> data_bylocus_MEDATL$cut
with(data_bylocus_MEDATL, table(cut)) %>% as_tibble() -> freq_occ_Fst_ATLMED_52outliers_960neutres

## on the 1012 loci within 21 sreg
x <- geno_827ATL_1012loci_21sreg
(mat.obs <- pairwise.fst(x, res.type="matrix")) # compute original Fst matrix
NBPERM <- 5000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) pairwise.fst(x, pop=sample(pop(x)), res.type="matrix"), mc.cores=no_cores)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) mat.perm[[k]][i,j])), mat.obs[i,j], alter="greater")$pvalue
  }
}
allTests->allTests_827ATL_1012loci_21sreg
data.frame(p=unlist(allTests_827ATL_1012loci_21sreg))->pvalues_827ATL_1012loci_21sreg
setwd("~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/Fst/21localities/")
write.table(pvalues_827ATL_1012loci_21sreg,row.names=T,col.names=T,sep="\t",quote=F,"pvalues_827ATL_1012loci_21sreg.txt")
as.dist(mat.obs)->Fst_827ATL_1012loci_21sreg
write.csv(mat.obs, "Fst_827ATL_1012loci_21sreg.csv")

## on 32 outliers within 21sreg
x<-geno_827ATL_32outliers_21sreg
(mat.obs <- pairwise.fst(x, res.type="matrix")) # compute original Fst matrix
NBPERM <- 5000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) pairwise.fst(x, pop=sample(pop(x)), res.type="matrix"), mc.cores=no_cores)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) mat.perm[[k]][i,j])), mat.obs[i,j], alter="greater")$pvalue
  }
}
allTests->allTests_827ATL_32outliers_21sreg
data.frame(p=unlist(allTests_827ATL_32outliers_21sreg))->pvalues_827ATL_32outliers_21sreg
setwd("~/sync/genetique_poissons/genotypage_bar/resultats/analyses_finales/Fst/21localities/")
write.table(pvalues_827ATL_32outliers_21sreg,row.names=T,col.names=T,sep="\t",quote=F,"pvalues_827ATL_32outliers_21sreg.txt")
as.dist(mat.obs)->Fst_827ATL_32outliers_21sreg
write.csv(mat.obs, "Fst_827ATL_32outliers_21sreg.csv")

## on 980 neutral within 21sreg
x<-geno_827ATL_980neutral_21sreg
(mat.obs <- pairwise.fst(x, res.type="matrix")) # compute original Fst matrix
NBPERM <- 5000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) pairwise.fst(x, pop=sample(pop(x)), res.type="matrix"), mc.cores=no_cores)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) mat.perm[[k]][i,j])), mat.obs[i,j], alter="greater")$pvalue
  }
}
allTests->allTests_827ATL_980neutral_21sreg
data.frame(p=unlist(allTests_827ATL_980neutral_21sreg))->pvalues_827ATL_980neutral_21sreg
setwd("~/sync/genetique_poissons/genotypage_bar/resultats/analyses_finales/Fst/21localities/")
write.table(pvalues_827ATL_980neutral_21sreg,row.names=T,col.names=T,sep="\t",quote=F,"pvalues_827ATL_980neutral_21sreg.txt")
as.dist(mat.obs)->Fst_827ATL_980neutral_21sreg
write.csv(mat.obs, "Fst_827ATL_980neutral_21sreg.csv")

## on 827ATL 10 MED 183 outliers within 23sreg
x<-geno_827ATL_10MED_183outliers_23sreg
mat.obs <- pairwise.fst(x, res.type="matrix") # compute original Fst matrix
write.csv(mat.obs, "Fst_827ATL_10MED_183outliers.csv")
NBPERM <- 5000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) pairwise.fst(x, pop=sample(pop(x)), res.type="matrix"), mc.cores=no_cores)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) mat.perm[[k]][i,j])), mat.obs[i,j], alter="greater")$pvalue
  }
}
allTests->allTests_827ATL_10MED_183outliers
data.frame(p=unlist(allTests_827ATL_10MED_183outliers))->pvalues_827ATL_10MED_183outliers
setwd("~/sync/genetique_poissons/genotypage_bar/resultats/analyses_finales/Fst/23localities/")
write.table(pvalues_827ATL_10MED_183outliers,row.names=T,col.names=T,sep="\t",quote=F,"pvalues_827ATL_10MED_183outliers.txt")
as.dist(mat.obs)->Fst_827ATL_10MED_183outliers_23reg

## on 827ATL 10 MED 829 neutral within 23sreg
x<-geno_827ATL_10MED_829neutral_23sreg
mat.obs <- pairwise.fst(x, res.type="matrix")
NBPERM <- 5000
mat.perm <- mclapply(1:NBPERM, function(i) pairwise.fst(x, pop=sample(pop(x)), res.type="matrix"), mc.cores=no_cores)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) mat.perm[[k]][i,j])), mat.obs[i,j], alter="greater")$pvalue
  }
}
allTests->allTests_827ATL_10MED_829neutral_23reg
data.frame(p=unlist(allTests_827ATL_10MED_829neutral_23reg))->pvalues_827ATL_10MED_829neutral_23reg
setwd("~/sync/genetique_poissons/genotypage_bar/resultats/analyses_finales/Fst/23localities/")
write.csv(mat.obs, "Fst_827ATL_10MED_829neutral_23reg.csv")
write.table(pvalues_827ATL_10MED_829neutral_23reg,row.names=T,col.names=T,sep="\t",quote=F,"pvalues_827ATL_10MED_829neutral_23reg.txt")
as.dist(mat.obs)->Fst_827ATL_10MED_829neutral_23reg
