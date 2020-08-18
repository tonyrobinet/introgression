#################
## test de Mantel
#################
library(adegenet)
library(tidyr)
library(dplyr)

setwd("/Volumes/data/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/Fst/21localities/")

read.csv("noms_827DLAB_ICESNAME_sreg_xy.csv",header=T)->noms_827DLAB
read.csv("noms_837DLAB_ICESNAME_sreg_xy.csv",header=T)-> noms_827ATL_10MED_reg

setwd("/Volumes/data/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/fichiers_ped_gtx/21localities/")
geno_827ATL_32outliers_21sreg<-read.genetix("FID_IID_32outliers_827ATL_21reg.gtx")
as.factor(noms_827DLAB$CLST)->geno_827ATL_32outliers_21sreg$pop
geno_827ATL_980neutral_21sreg<-read.genetix("FID_IID_980neutral_827ATL_21reg.gtx")
as.factor(noms_827DLAB$CLST)->geno_827ATL_980neutral_21sreg$pop

seppop(geno_827ATL_980neutral_21sreg)->x
seppop(geno_827ATL_32outliers_21sreg)->y
read.csv("/Volumes/data/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/IBD/21localities/9nordGONB_geoDist.csv") -> Dgeo_nordGONB
dist(Dgeo_nordGONB) -> Dgeo_nordGONB
read.csv("/Volumes/data/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/IBD/21localities/12sudGONB_geoDist.csv") -> Dgeo_sudGONB
dist(Dgeo_sudGONB) -> Dgeo_sudGONB


# nord GONB neutres
repool(x$ANGL,x$BSEI,x$CRNW,x$DOVE,x$IRIS,x$ISWI,x$NEUK,x$SGEO,x$WALE)->geno_827ATL_980neutral_9pops_nordGONB
popNames(geno_827ATL_980neutral_9pops_nordGONB) # [1] "ANGL" "BSEI" "CRNW" "DOVE" "IRIS" "ISWI" "NEUK" "SGEO" "WALE"
genpop_827ATL_980neutral_nordGONB <- genind2genpop(geno_827ATL_980neutral_9pops_nordGONB)
Dgen_neutres_nordGONB <- dist.genpop(genpop_827ATL_980neutral_nordGONB,method=1)
mantel_neutres_nordGONB <- mantel.randtest(Dgen_neutres_nordGONB,Dgeo_nordGONB,nrepet=10000)
# Simulated p-value: 0.01509849 / 0.00659934

# sud GONB neutres
repool(x$AQUI,x$ASTU,x$CANT,x$CHAR,x$CORU,x$GONB,x$LOIR,x$PBRE,x$PENI,x$PORT,x$SINE,x$VIGO)->geno_827ATL_980neutral_12pops_sudGONB
popNames(geno_827ATL_980neutral_12pops_sudGONB) #  [1] "AQUI" "ASTU" "CANT" "CHAR" "GALI" "GONB" "LOIR" "PBRE" "PENI" "PORT" "SINE"
genpop_827ATL_980neutral_sudGONB <- genind2genpop(geno_827ATL_980neutral_12pops_sudGONB)
Dgen_neutres_sudGONB <- dist.genpop(genpop_827ATL_980neutral_sudGONB,method=1)
mantel_neutres_sudGONB <- mantel.randtest(Dgen_neutres_sudGONB,Dgeo_sudGONB,nrepet=10000)
# Simulated p-value: 0.5317468

# nord GONB outliers
repool(y$ANGL,y$BSEI,y$CRNW,y$DOVE,y$IRIS,y$ISWI,y$NEUK,y$SGEO,y$WALE)->geno_827ATL_32outliers_9pops_nordGONB
popNames(geno_827ATL_32outliers_9pops_nordGONB) # [1] "ANGL" "BSEI" "CRNW" "DOVE" "IRIS" "ISWI" "NEUK" "SGEO" "WALE"
genpop_827ATL_32outliers_nordGONB <- genind2genpop(geno_827ATL_32outliers_9pops_nordGONB)
Dgen_outliers_nordGONB <- dist.genpop(genpop_827ATL_32outliers_nordGONB,method=1)
mantel_outliers_nordGONB <- mantel.randtest(Dgen_outliers_nordGONB,Dgeo_nordGONB,nrepet=10000)
# Simulated p-value: 0.01109889 / 0.06539346

# sud GONB outliers
repool(y$AQUI,y$ASTU,y$CANT,y$CHAR,y$CORU,y$GONB,y$LOIR,y$PBRE,y$PENI,y$PORT,y$SINE,y$VIGO)->geno_827ATL_32outliers_12pops_sudGONB
popNames(geno_827ATL_32outliers_12pops_sudGONB) #  [1] "AQUI" "ASTU" "CANT" "CHAR" "GALI" "GONB" "LOIR" "PBRE" "PENI" "PORT" "SINE"
genpop_827ATL_32outliers_sudGONB <- genind2genpop(geno_827ATL_32outliers_12pops_sudGONB)
Dgen_outliers_sudGONB <- dist.genpop(genpop_827ATL_32outliers_sudGONB,method=1)
mantel_outliers_sudGONB <- mantel.randtest(Dgen_outliers_sudGONB,Dgeo_sudGONB,nrepet=10000)
# Simulated p-value:  0.1306869

# Gascogne outliers
repool(y$AQUI,y$ASTU,y$CANT,y$CHAR,y$GONB,y$LOIR,y$PBRE)->geno_827ATL_32outliers_7pops_Gascogne
genpop_827ATL_32outliers_Gascogne <- genind2genpop(geno_827ATL_32outliers_7pops_Gascogne)
Dgen_outliers_Gascogne <- dist.genpop(genpop_827ATL_32outliers_Gascogne,method=1)
read.csv("/Volumes/data/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/IBD/21localities/7Gascogne_geoDist.csv") -> Dgeo_Gascogne
dist(Dgeo_Gascogne) -> Dgeo_Gascogne
mantel_outliers_Gascogne <- mantel.randtest(Dgen_outliers_Gascogne,Dgeo_Gascogne,nrepet=10000)
# Simulated p-value: 0.520448

# Gascogne neutres
repool(x$AQUI,x$ASTU,x$CANT,x$CHAR,x$GONB,x$LOIR,x$PBRE)->geno_827ATL_980neutral_7pops_Gascogne
genpop_827ATL_980neutral_7pops_Gascogne <- genind2genpop(geno_827ATL_980neutral_7pops_Gascogne)
Dgen_neutres_Gascogne <- dist.genpop(genpop_827ATL_980neutral_7pops_Gascogne,method=1)
mantel_neutres_Gascogne <- mantel.randtest(Dgen_neutres_Gascogne,Dgeo_Gascogne,nrepet=10000)
# Simulated p-value: 0.1754825 

# Portugal outliers
repool(y$CORU,y$PENI,y$PORT,y$SINE,y$VIGO)->geno_827ATL_32outliers_5pops_Portugal
genpop_827ATL_32outliers_Portugal <- genind2genpop(geno_827ATL_32outliers_5pops_Portugal)
Dgen_outliers_Portugal <- dist.genpop(genpop_827ATL_32outliers_Portugal,method=1)
read.csv("/Volumes/data/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/IBD/21localities/5Portugal_geoDist.csv") -> Dgeo_Portugal
dist(Dgeo_Portugal) -> Dgeo_Portugal
mantel_outliers_Portugal <- mantel.randtest(Dgen_outliers_Portugal,Dgeo_Portugal,nrepet=10000)
# Simulated p-value: 0.1783822

# Portugal neutres
repool(x$CORU,x$PENI,x$PORT,x$SINE,x$VIGO)->geno_827ATL_980neutral_5pops_Portugal
genpop_827ATL_980neutral_5pops_Portugal <- genind2genpop(geno_827ATL_980neutral_5pops_Portugal)
Dgen_neutres_Portugal <- dist.genpop(genpop_827ATL_980neutral_5pops_Portugal,method=1)
mantel_neutres_Portugal <- mantel.randtest(Dgen_neutres_Portugal,Dgeo_Portugal,nrepet=10000)
# Simulated p-value: 0.9684032 
