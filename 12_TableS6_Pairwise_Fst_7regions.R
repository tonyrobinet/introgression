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

geno_827ATL_10MED_1012loci_9reg <- read.genetix("~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/fichiers_ped_gtx/23localities/FID_IID_1012loci_837ALTMED_23reg.gtx")
as.factor(noms_827ATL_10MED$reg) -> geno_827ATL_10MED_1012loci_9reg$pop

geno_827ATL_10MED_52outliers_9reg <- read.genetix("~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/fichiers_ped_gtx/23localities/FID_IID_52outliers_837ALTMED_23reg.gtx")
as.factor(noms_827ATL_10MED$reg) -> geno_827ATL_10MED_52outliers_9reg$pop
geno_827ATL_10MED_960neutral_9reg <- read.genetix("~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/fichiers_ped_gtx/23localities/FID_IID_960neutral_837ALTMED_23reg.gtx")
as.factor(noms_827ATL_10MED$reg) -> geno_827ATL_10MED_960neutral_9reg$pop

geno_827ATL_32outliers_7reg <- read.genetix("~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/fichiers_ped_gtx/21localities/FID_IID_32outliers_827ATL_21reg.gtx")
as.factor(noms_827DLAB$reg) -> geno_827ATL_32outliers_7reg$pop
geno_827ATL_980neutral_7reg <- read.genetix("~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/fichiers_ped_gtx/21localities/FID_IID_980neutral_827ATL_21reg.gtx")
as.factor(noms_827DLAB$reg) -> geno_827ATL_980neutral_7reg$pop

geno_827ATL_1012loci_7reg <- read.genetix("/Volumes/data/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/fichiers_ped_gtx/21localities/FID_IID_1012loci_827ATL_21reg.gtx")
as.factor(noms_827DLAB$reg) -> geno_827ATL_1012loci_7reg$pop

#    over 9 regions (7 ATL and 2 MED) : on the 1012 loci
x=pvalues=mat.obs=mat.perm=NPERM=allTests=NULL
x <- geno_827ATL_10MED_1012loci_9reg
(mat.obs <- pairwise.fst(x, res.type="matrix")) # compute original Fst matrix
write.csv(as.matrix(mat.obs),"~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/Fst/9regions/matobs_827ATL_10MED_1012loci_9reg.csv")
NBPERM <- 5000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM,function(i) pairwise.fst(x, pop=sample(pop(x)), res.type="matrix"), mc.cores=no_cores)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) mat.perm[[k]][i,j])), mat.obs[i,j], alter="greater")$pvalue
  }
}
data.frame(p=unlist(allTests))->pvalues
write.table(pvalues,row.names=T,col.names=T,sep="\t",quote=F,"~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/Fst/9regions/pvalues_827ATL_10MED_1012loci_9reg_N5000.txt")
as.dist(mat.obs, upper=T)


#    over 9 regions (7 ATL and 2 MED): on 52 outliers ATL-MED
x=pvalues=mat.obs=mat.perm=NPERM=allTests=NULL
x <- geno_827ATL_10MED_52outliers_9reg
(mat.obs <- pairwise.fst(x, res.type="matrix")) # compute original Fst matrix
write.csv(as.matrix(mat.obs),"~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/Fst/9regions/matobs_827ATL_10MED_52outliers_9reg.csv")
NBPERM <- 5000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM,function(i) pairwise.fst(x, pop=sample(pop(x)), res.type="matrix"), mc.cores=no_cores)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) mat.perm[[k]][i,j])), mat.obs[i,j], alter="greater")$pvalue
  }
}
data.frame(p=unlist(allTests))->pvalues
write.table(pvalues,row.names=T,col.names=T,sep="\t",quote=F,"~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/Fst/9regions/pvalues_827ATL_10MED_52outliers_9reg_N5000.txt")
as.dist(mat.obs, upper=T)


#    over 9 regions (7 ATL et 2 MED): on 960 neutral ATL-MED
x=pvalues=mat.obs=mat.perm=NPERM=allTests=NULL
x <- geno_827ATL_10MED_960neutral_9reg
(mat.obs <- pairwise.fst(x, res.type="matrix")) # compute original Fst matrix
as.dist(mat.obs, upper=T)
write.csv(mat.obs,"~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/Fst/9regions/matobs_827ATL_10MED_960neutral_9reg.csv")
NBPERM <- 5000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM,function(i) pairwise.fst(x, pop=sample(pop(x)), res.type="matrix"), mc.cores=no_cores)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) mat.perm[[k]][i,j])), mat.obs[i,j], alter="greater")$pvalue
  }
}
data.frame(p=unlist(allTests))->pvalues
write.table(pvalues,row.names=T,col.names=T,sep="\t",quote=F,"~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/Fst/9regions/pvalues_827ATL_10MED_960neutral_9reg_N5000.txt")

#    over 7 regions ATL: on the 1012 loci
x=pvalues=mat.obs=mat.perm=NPERM=allTests=NULL
x <- geno_827ATL_1012loci_7reg
(mat.obs <- pairwise.fst(x, res.type="matrix")) # compute original Fst matrix
as.dist(mat.obs, upper=T)
write.csv(mat.obs,"/Volumes/data/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/Fst/7regions/matobs_827ATL_1012loci_7reg.csv")
NBPERM <- 5000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM,function(i) pairwise.fst(x, pop=sample(pop(x)), res.type="matrix"), mc.cores=no_cores)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) mat.perm[[k]][i,j])), mat.obs[i,j], alter="greater")$pvalue
  }
}
data.frame(p=unlist(allTests))->pvalues
write.table(pvalues,row.names=T,col.names=T,sep="\t",quote=F,"/volumes/data/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/Fst/7regions/pvalues_827ATL_1012loci_7reg_N5000.txt")
as.dist(pvalues, upper=T)


#    over 7 regions ATL: on 32 outliers intra-ATL
x=pvalues=mat.obs=mat.perm=NPERM=allTests=NULL
x <- geno_827ATL_32outliers_7reg
(mat.obs <- pairwise.fst(x, res.type="matrix")) # compute original Fst matrix
as.dist(mat.obs, upper=T)
write.csv(mat.obs,"~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/Fst/7regions/matobs_827ATL_32outliers_7reg.csv")
NBPERM <- 5000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM,function(i) pairwise.fst(x, pop=sample(pop(x)), res.type="matrix"), mc.cores=no_cores)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) mat.perm[[k]][i,j])), mat.obs[i,j], alter="greater")$pvalue
  }
}
data.frame(p=unlist(allTests))->pvalues
write.table(pvalues,row.names=T,col.names=T,sep="\t",quote=F,"~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/Fst/7regions/pvalues_827ATL_32outliers_7reg_N5000.txt")
as.dist(pvalues, upper=T)


#    over 7 regions ATL: on 980 neutres intra-ATL
x=pvalues=mat.obs=mat.perm=NPERM=allTests=NULL
x <- geno_827ATL_980neutral_7reg
(mat.obs <- pairwise.fst(x, res.type="matrix")) # compute original Fst matrix
as.dist(mat.obs, upper=T)
write.csv(mat.obs,"~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/Fst/7regions/matobs_827ATL_980neutral_7reg.csv")
NBPERM <- 5000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM,function(i) pairwise.fst(x, pop=sample(pop(x)), res.type="matrix"), mc.cores=no_cores)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) mat.perm[[k]][i,j])), mat.obs[i,j], alter="greater")$pvalue
  }
}
data.frame(p=unlist(allTests))->pvalues
write.table(pvalues,row.names=T,col.names=T,sep="\t",quote=F,"~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/Fst/7regions/pvalues_827ATL_980neutral_7reg_N5000.txt")
as.dist(pvalues, upper=T)

# end
print("happy end")
