# Stevick et al 2020 Oyster Gut Microbiome Function in an Estuary
# Script to analyze gene-level functional expression from SAMSA
# Figures S5, S6, S7, S8
# updated 4/16/2020 for resubmission

# Metaranscriptomic functional data at gene level


library("ggplot2")
library("tidyr")
library("dplyr")
library("readxl")
library("gridExtra")
library("leaflet")
library(RColorBrewer)
library(viridis)
library(pheatmap)

data<-read_excel("Function/DESeq_subsys_results_level4_summary.xlsx", sheet=1)
hierarchy<-read_excel("Function/DESeq_subsys_results_level4_summary.xlsx", sheet=2)

data$padj_PVD<-as.numeric(data$padj_PVD)
data$padj_GB<-as.numeric(data$padj_GB)
data$padj_BIS<-as.numeric(data$padj_BIS)
data$padj_NAR<-as.numeric(data$padj_NAR)
data$padj_NIN<-as.numeric(data$padj_NIN)


datah<-merge(data, hierarchy)

# NITROGEN ------------------------------------

datahdenit<-filter(datah, Level2=="Nitrogen Metabolism")
datahdenitg<-gather(datahdenit, foldchange, foldvalue, "log2FoldChange_PVD",
                    "log2FoldChange_GB","log2FoldChange_BIS",
                    "log2FoldChange_NAR","log2FoldChange_NIN")
datahdenitgp<-gather(datahdenit, padj, padjvalue, "padj_PVD","padj_GB","padj_BIS",
                     "padj_NAR","padj_NIN")
datahdenitg$padj<-datahdenitgp$padj
datahdenitg$padjvalue<-datahdenitgp$padjvalue

datahdenitg$signif<-ifelse(datahdenitgp$padjvalue <=0.01, "**", 
                        ifelse(datahdenitgp$padjvalue <=0.05, "*",
                               ifelse("ns")))

datahdenitg$foldchange = factor(datahdenitg$foldchange, 
                             levels = c("log2FoldChange_PVD",
                    "log2FoldChange_GB","log2FoldChange_BIS",
                    "log2FoldChange_NAR","log2FoldChange_NIN"))

ggplot(datahdenitg, aes(Level4, foldvalue,fill=foldchange))+
  geom_col(position = "dodge")+coord_flip()+
  facet_grid(Level3~., scales="free",space="free")+
  theme(strip.text.y = element_text(angle=0))+
  scale_fill_manual(values=c("#253494","#0868ac","#43a2ca","#7bccc4","#bae4bc"))+
   geom_text(aes(label=signif), position = position_dodge(width = 1), size=6)

datahdenitg$Level3f<-factor(datahdenitg$Level3, 
                           levels=c("Nitric_oxide_synthase","Allantoin_Utilization",
                                    "Amidase_clustered_with_urea_and_nitrile_hydratase_functions",
                                    "Denitrification",
                                    "Dissimilatory_nitrite_reductase",
                                    "Nitrosative_stress",
                                    "Cyanate_hydrolysis","Ammonia_assimilation",
                                    "Nitrate_and_nitrite_ammonification"))



ggplot(datahdenitg, aes(Level4, foldchange, label=signif, 
                         fill=foldvalue,color=signif))+theme_minimal()+
  geom_tile(size=0.7)+facet_grid(Level3f~., scales="free", space="free")+
  scale_fill_viridis()+coord_flip()+
  labs(x=NULL, y=NULL)+
  geom_text(size=5, nudge_x=-0.3)+
  theme(legend.position="bottom",
        plot.background = element_rect(fill = "transparent", color = NA),rect = element_rect(fill = "transparent"))+
  scale_color_manual(values=c("darkred","red","black"))



### PHOSPHOROUS -----------------------


datahphos<-filter(datah, Level2=="Phosphorus Metabolism")
datahphosg<-gather(datahphos, foldchange, foldvalue, "log2FoldChange_PVD",
                    "log2FoldChange_GB","log2FoldChange_BIS",
                    "log2FoldChange_NAR","log2FoldChange_NIN")
datahphosgp<-gather(datahphos, padj, padjvalue, "padj_PVD","padj_GB","padj_BIS",
                     "padj_NAR","padj_NIN")
datahphosg$padj<-datahphosgp$padj
datahphosg$padjvalue<-datahphosgp$padjvalue

datahphosg$signif<-ifelse(datahphosgp$padjvalue <=0.01, "**", 
                           ifelse(datahphosgp$padjvalue <=0.05, "*",
                                  ifelse("ns")))

datahphosg$foldchange = factor(datahphosg$foldchange, 
                                levels = c("log2FoldChange_PVD",
                                           "log2FoldChange_GB","log2FoldChange_BIS",
                                           "log2FoldChange_NAR","log2FoldChange_NIN"))

ggplot(datahphosg, aes(Level4, foldvalue,fill=foldchange))+
  geom_col(position = "dodge")+coord_flip()+
  facet_grid(Level3~., scales="free",space="free")+
#  theme(strip.text.y = element_text(angle=0))+
  scale_fill_manual(values=c("#253494","#0868ac","#43a2ca","#7bccc4","#bae4bc"))+
  geom_text(aes(label=signif), position = position_dodge(width = 1), size=6)


ggplot(datahphosg,aes(Level4, foldchange, label=signif,
                      fill=foldvalue,color=signif))+theme_minimal()+
  geom_tile(size=0.7)+facet_grid(Level3~., scales="free", space="free")+
  scale_fill_viridis()+coord_flip()+
  labs(x=NULL, y=NULL)+
  geom_text(size=5, nudge_x=-0.3)+
  theme(legend.position="bottom",
        plot.background = element_rect(fill = "transparent", color = NA),rect = element_rect(fill = "transparent"))+
  scale_color_manual(values=c("darkred","red","black"))





### STRESS RESPONSE -----------------------


datahstress<-filter(datah,Level2=="Osmotic stress" | Level2=="Periplasmic Stress" | 
                    Level2=="Oxidative stress")
datahstressg<-gather(datahstress, foldchange, foldvalue, "log2FoldChange_PVD",
                   "log2FoldChange_GB","log2FoldChange_BIS",
                   "log2FoldChange_NAR","log2FoldChange_NIN")
datahstressgp<-gather(datahstress, padj, padjvalue, "padj_PVD","padj_GB","padj_BIS",
                    "padj_NAR","padj_NIN")
datahstressg$padj<-datahstressgp$padj
datahstressg$padjvalue<-datahstressgp$padjvalue

datahstressg$signif<-ifelse(datahstressgp$padjvalue <=0.01, "**", 
                          ifelse(datahstressgp$padjvalue <=0.05, "*",
                                 ifelse("ns")))

datahstressg$foldchange = factor(datahstressg$foldchange, 
                               levels = c("log2FoldChange_PVD",
                                          "log2FoldChange_GB","log2FoldChange_BIS",
                                          "log2FoldChange_NAR","log2FoldChange_NIN"))


data3gpfpe<-gather(datahstress, folderror, errorvalue, "lfcSE_PVD","lfcSE_GB",
                   "lfcSE_BIS","lfcSE_NAR","lfcSE_NIN")
datahstressg$folderror<-data3gpfpe$folderror
datahstressg$errorvalue<-data3gpfpe$errorvalue


ggplot(datahstressg, aes(Level4, foldvalue,fill=foldchange))+
  geom_col(position = "dodge")+coord_flip()+
  facet_grid(Level2~., scales="free", space="free")+
  geom_errorbar(aes(ymin=foldvalue-errorvalue, ymax=foldvalue+errorvalue), 
                size=0.8,position = position_dodge(width=0.9), width = 0, color="grey60")+
 #   theme(strip.text.y = element_text(angle=0))+
  scale_fill_manual(values=c("#253494","#0868ac","#43a2ca","#7bccc4","#bae4bc"))+
  geom_text(aes(label=signif), position = position_dodge(width = 1), size=6)



ggplot(datahstressg, aes(Level4, foldchange, label=signif,
                         fill=foldvalue,color=signif))+theme_minimal()+
  geom_tile(size=0.7)+facet_grid(Level2~., scales="free", space="free")+
  scale_fill_viridis()+coord_flip()+
  labs(x=NULL, y=NULL)+
  geom_text(size=5, nudge_x=-0.3)+
  theme(legend.position="bottom",
        plot.background = element_rect(fill = "transparent", color = NA),rect = element_rect(fill = "transparent"))+
  scale_color_manual(values=c("darkred","red","black"))





#write.csv(datahdenit,"nitrogen.csv")
#write.csv(datahphos,"phosphorus.csv")
#write.csv(datahstress,"stress.csv")
