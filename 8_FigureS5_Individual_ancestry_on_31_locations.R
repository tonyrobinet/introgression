library(ggplot2)
library(labeling)
library(tidyr)
library(dplyr)
library(grid)
library(gridExtra)
library(ggrepel)
library(cowplot)
mytheme <- theme_classic()
mytheme$axis.line.x <- mytheme$axis.line.y <- mytheme$axis.line

## The ADMIXTURE software was run for K=2 on "FID_IID_DLAB1012loci_827ind.bed" for ATL. The corresponding .map file is needed in the same directory.
$ ./admixture32 FID_IID_DLAB1012loci_827ind.bed 2

# sur 761 DLAB 31 CLST
"~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/admixture/sur_tous_marqueurs/761DLAB/res_admixture_FID_IID_1012loci_761ind_31CLST.csv" %>%
  read.csv() %>% as_tibble() -> res_admixture_FID_IID_1012loci_761ind_31CLST
"~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/fichiers_ped_gtx/31localities/dist_SINE_long_lat_31reg.csv" %>%
  read.csv() %>% as_tibble() -> dist_SINE_long_lat_31reg
"~/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/freq_alleliques/Freq_alleliq_761ATL_52outliersATLMED_31reg.csv" %>%
  read.csv() %>% as_tibble() -> MAF_761ATL_52ATLMED_outliers

merge(res_admixture_FID_IID_1012loci_761ind_31CLST,dist_SINE_long_lat_31reg, by = "CLST",  all=T) %>% as_tibble() -> data_761
MAF_761ATL_52ATLMED_outliers %>% group_by(CLST) %>% summarize(MAF_avg=mean(MAF), MAF_sd=sd(MAF)) -> MAF_761ATL_52ATLMED_outliers
merge(MAF_761ATL_52ATLMED_outliers,data_761,all=T) %>% as_tibble() -> data_761

data_761 %>% filter(BRANCH=="SOUTH") -> data_SOUTH
data_761 %>% filter(BRANCH=="IRISH") -> data_IRISH
data_761 %>% filter(BRANCH=="CHANNEL") -> data_CHANNEL
full_join(data_SOUTH,data_CHANNEL) -> data_SOUTH_CHANNEL
full_join(data_SOUTH,data_IRISH) -> data_SOUTH_IRISH

data_SOUTH_IRISH %>% group_by(CLST, BRANCH) %>% summarise(avg_Q2=mean(Q2), sd_Q2=sd(Q2), MAF_avg=mean(MAF_avg), MAF_sd=mean(MAF_sd), dist_SINE=mean(dist_SINE)) %>%
  filter(dist_SINE>=100) -> data_SOUTH_IRISH_avg

data_SOUTH_CHANNEL %>% group_by(CLST, BRANCH) %>% summarise(avg_Q2=mean(Q2), sd_Q2=sd(Q2), MAF_avg=mean(MAF_avg), MAF_sd=mean(MAF_sd), dist_SINE=mean(dist_SINE)) %>%
  filter(dist_SINE>=100) -> data_SOUTH_CHANNEL_avg

full_join(data_SOUTH_CHANNEL_avg, data_SOUTH_IRISH_avg) -> data_761_avg

r <- ggplot(data=data_761_avg,aes(x=dist_SINE, y=avg_Q2, color = BRANCH)) + geom_line(size=1, shape=21) + 
  geom_text_repel(aes(label=CLST)) #+ geom_errorbar(aes(ymin=avg_Q2-sd_Q2, ymax=avg_Q2+sd_Q2)) 
v <- r +  mytheme + labs(x="Distance to SINE by the plateau", y="mean MED ancestrality") +
  theme(legend.title = element_blank(),legend.background=element_rect(fill="grey95"),
        legend.position="right", legend.text = element_text(size=14),
        axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),
        axis.text.x  = element_text(size=14), axis.text.y  = element_text(size=14)) + 
  scale_x_continuous(limits=c(100,3050)) #+ scale_y_continuous(limits=c(0,0.1))

s <- ggplot(data=data_761_avg,aes(x=dist_SINE, y=MAF_avg, color = BRANCH)) + geom_line(size=1, shape=21) + 
  geom_text_repel(aes(label=CLST)) #+ geom_errorbar(aes(ymin=avg_Q2-sd_Q2, ymax=avg_Q2+sd_Q2)) 
t <- s +  mytheme + labs(x="Distance to SINE by the plateau", y="mean MAF outliers MED-ATL") +
  theme(legend.title = element_blank(),legend.background=element_rect(fill="grey95"),
        legend.position="right", legend.text = element_text(size=14),
        axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),
        axis.text.x  = element_text(size=14), axis.text.y  = element_text(size=14)) + 
  scale_x_continuous(limits=c(100,3050)) #+ scale_y_continuous(limits=c(0.025,0.05))

v <- plot_grid(v , labels = "a")
t <- plot_grid(t , labels = "b")

grid.arrange(v,t,nrow=2)

# dev.print(png, 
#           '~/sync/poissons/genetique_poissons/genotypage_bar/drafts/figures/MAF_MED_ancestr_dist_SINE_2_branches.png', 
#           width=800, height=800)

