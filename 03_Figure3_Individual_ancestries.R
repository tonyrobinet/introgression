## Individual admixture plot
############################

# The code below generates a plot of inferred individual Mediterranean ancestry. 
# The program "Admixture" (Alexander, Novembre & Lange, 2009), which infers individual ancestry proportions from K ancestral populations (here K was set to 2), was ran untill convergence using 100 iterations. 
# The standard error of individual Mediterranean ancestry was estimated using 1000 bootstrap resamplings, each consisting of 10 iterations.
# The ADMIXTURE software was run for K=2 on "Merged_827_ChipsATL_10_GenomesMED_1012_SNPs_NoRep.bed" for ATL-MED. The corresponding .map file is needed in the same directory.
$ ./admixture32 Merged_827_ChipsATL_10_GenomesMED_1012_SNPs_NoRep.bed 2

################
library(ggplot2)
library(dplyr)
mytheme <- theme_classic()
mytheme$axis.line.x <- mytheme$axis.line.y <- mytheme$axis.line


"~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/admixture/noms_837DLAB_ATLMED_ICESNAME_sreg_xy.csv" %>%
  read.csv() %>% as_tibble() -> noms
"~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/admixture/dist_SINE_long_lat_23reg_VIGO_CORO.csv" %>%
  read.csv() %>% as_tibble() -> dist_SINE_lat
"~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/admixture/modele_admixture_PAG/ATL-MED_admixture_100iterations/Merged_827_Chips_10_Genomes_1012_SNPs_NoRep.2.Q" %>%
  read.csv(sep=" ",header=T) %>% as_tibble() -> data
"~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/admixture/modele_admixture_PAG/ATL-MED_admixture_100iterations/Merged_827_Chips_10_Genomes_1012_SNPs_NoRep.2.Q_se" %>%
  read.csv(sep=" ",header=F) %>% as_tibble() %>% rename(se1=V1,se2=V2) -> data_se

merge(noms,dist_SINE_lat,all=T) %>% as_tibble() %>% arrange(rank) -> noms # 837 rows
data %>% mutate(id_specimen=noms$id_specimen) %>% mutate(CLST=noms$CLST) %>% mutate(dist_SINE=noms$dist_SINE) %>% mutate(group=noms$group) %>%
  mutate(se=data_se$se2) %>% mutate(MIN=Q2-se) %>% mutate(MAX=Q2+se) %>% mutate(MIN=ifelse(MIN<0,0,MIN)) %>% mutate(MAX=ifelse(MAX>1,1,MAX)) -> data
data %>% arrange(dist_SINE,CLST) %>% mutate(POS=c(1:837))-> data

data %>% group_by(CLST) %>% summarise_all(funs(mean)) %>% as_tibble() -> data_mean
merge(dist_SINE_lat,data_mean, by="CLST") %>% as_tibble() -> data_mean
data_mean %>% filter(dist_SINE.x>=0) %>% mutate(CLST=CLST, dist_SINE=dist_SINE.x, group=group.x, long=long, lat=lat, Q1=Q1, Q2=Q2, se=se, MIN=MIN, MAX=MAX, POS=POS) -> data_mean_ATL


# noms %>% arrange(CLST) %>% print(n = Inf)
# data %>% arrange(CLST) %>% print(n = Inf)
# data_mean %>% print(n = Inf)
# data_mean_ATL %>% print(n = Inf)
# dist_SINE_lat %>% print(n = Inf)


# Make the plot
h <- ggplot(data, aes(x=POS,y=Q2)) + geom_point(aes(colour=CLST))

MED_ancestry_individual <- h + geom_errorbar(aes(ymin=MIN,ymax=MAX,colour=CLST)) +
  labs(x="Individuals within localities ranked by latitude", y = "Mediterranean ancestry") +
  mytheme + theme(legend.justification=c(0.8,1), legend.position=c(0.8,1)) +
  theme(legend.title = element_blank(),legend.background=element_rect(fill="grey95"),
        legend.text=element_text(size=11), axis.text=element_text(size=12), axis.title=element_text(size=14)) + 
  guides(colour = guide_legend(nrow = 3)) +  scale_color_discrete(breaks=c("EMED","WMED","SINE","PENI","PORT","VIGO","CORU","ASTU","CANT","AQUI","CHAR","LOIR","PBRE","GONB","CRNW","ISWI","WALE","BSEI","SGEO","DOVE","ANGL","NEUK")) 

MED_ancestry_individual