###################################
# Patterns of Isolation By Distance
###################################

library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(wesanderson)
library(dplyr)
mytheme <- theme_classic()
mytheme$axis.line.x <- mytheme$axis.line.y <- mytheme$axis.line


setwd("/Volumes/data/sync/poissons/genetique_poissons/genotypage_bar/resultats/analyses_finales/IBD/21localities/")

dgeo_fst <- "IBD_neutres_outliers_1dataset.csv" %>% read.csv() %>% as_tibble
dgeo_fst_2datasets <- "IBD_neutres_outliers_2datasets.csv" %>% read.csv() %>% as_tibble
dgeo_fst_sud_GONB <- "IBD_sud_GONB.csv" %>% read.csv() %>% as_tibble
dgeo_fst_nord_GONB <- "IBD_nord_GONB.csv" %>% read.csv() %>% as_tibble
combi_nord_sud_GONB <- rbind(dgeo_fst_nord_GONB,dgeo_fst_sud_GONB) %>% as_tibble()

toutes_localites <- ggplot(dgeo_fst_2datasets, aes(x=dist, y=Fst_Rousset, col=dataset)) + geom_jitter(size=.7) +
  scale_color_manual(values=wes_palette(n=2, name="Darjeeling1")) + geom_smooth(method="lm", size=.7) + expand_limits(x=c(0,3050)) +
  mytheme + theme(axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), 
                  legend.title = element_blank(),legend.background=element_rect(fill="grey95"), 
                  legend.position="none", legend.text=element_text(size=11), axis.text.x=element_text(size=13), axis.text.y=element_text(size=13)) +
  labs(x="distance by the plateau", y="Fst / (1  - Fst)")

neutres_nord_sudGONB <- ggplot(combi_nord_sud_GONB, aes(x= dist, y=Fst_Rousset_980neutral, col=dataset), ylim=c(0,0.07)) +
  scale_color_manual(values=c("#330099", "#E69F00")) + geom_jitter(size=.7) + mytheme +
  geom_smooth(method="lm", size=.7) + theme_classic() + expand_limits(y=c(0,0.06), x=c(0,3000)) +
  theme(axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), 
        legend.title = element_blank(),legend.background=element_rect(fill="grey95"), 
        legend.position="none", legend.text=element_text(size=13), axis.text.x=element_text(size=13), axis.text.y=element_text(size=13)) +
  labs(y="Fst / (1-Fst)", x="distance by the plateau")

outliers_nord_sudGONB <- ggplot(combi_nord_sud_GONB, aes(x= dist, y=Fst_Rousset_32outliers, col=dataset)) + 
  scale_color_manual(values=c("#330099", "#E69F00")) + geom_jitter(size=.7) + mytheme +
  geom_smooth(method="lm", size=.7) + theme_classic() + expand_limits(y=c(0,0.06), x=c(0,3000)) +
  theme(axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), 
        legend.title = element_blank(),legend.background=element_rect(fill="grey95"), 
        legend.position="none", legend.text=element_text(size=13), axis.text.x=element_text(size=13), axis.text.y=element_text(size=13)) +
  labs(y="Fst / (1-Fst)", x="distance by the plateau")

a <- arrangeGrob(toutes_localites, top = textGrob("a", x = unit(0, "npc")
                                                  , y   = unit(1, "npc"), just=c("left","top"),
                                                  gp=gpar(col="black", fontsize=14, fontfamily="Times Roman")))
b <- arrangeGrob(neutres_nord_sudGONB, top = textGrob("b", x = unit(0, "npc")
                                                      , y   = unit(1, "npc"), just=c("left","top"),
                                                      gp=gpar(col="black", fontsize=14, fontfamily="Times Roman")))
c <- arrangeGrob(outliers_nord_sudGONB, top = textGrob("c", x = unit(0, "npc")
                                                       , y   = unit(1, "npc"), just=c("left","top"),
                                                       gp=gpar(col="black", fontsize=14, fontfamily="Times Roman")))
# Make the plot
grid.arrange(a,b,c,nrow=3,ncol=1)

# IBD models
model_outliers_all <- lm(dgeo_fst$Fst_Rousset_32outliers~dgeo_fst$dist)
model_neutres_all <- lm(dgeo_fst$Fst_Rousset_980neutral~dgeo_fst$dist)
model_neutres_nordGONB <- lm(dgeo_fst_nord_GONB$Fst_Rousset_980neutral~dgeo_fst_nord_GONB$dist)
model_neutres_sudGONB <- lm(dgeo_fst_sud_GONB$Fst_Rousset_980neutral~dgeo_fst_sud_GONB$dist)
model_outliers_nordGONB <- lm(dgeo_fst_nord_GONB$Fst_Rousset_32outliers~dgeo_fst_nord_GONB$dist)
model_outliers_sudGONB <- lm(dgeo_fst_sud_GONB$Fst_Rousset_32outliers~dgeo_fst_sud_GONB$dist)

summary(model_outliers_all) #Adjusted R-squared:  0.08902 / p-value: 6.474e-06
summary(model_neutres_all) #Adjusted R-squared:  0.07423  / p-value: 3.737e-05
summary(model_neutres_nordGONB) #Adjusted R-squared:  0.2359  / p-value: 0.001574
summary(model_outliers_nordGONB) # Adjusted R-squared:   0.245 / p-value: 0.001266
summary(model_neutres_sudGONB) # Adjusted R-squared:  -0.0109  / p-value: p-value: 0.5863
summary(model_outliers_sudGONB) # Adjusted R-squared:  0.08277 / p-value: 0.01096
