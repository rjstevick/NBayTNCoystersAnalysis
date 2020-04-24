# Stevick et al 2020 Oyster Gut Microbiome Function in an Estuary
# High level functional groups from SAMSA, North compared to South
# Stress response, phosphorus, and nitrogen gene expression at each site compared to the mean of all the others
# Figures 5 and 6
# updated 4/16/2020 for resubmission

# Metaranscriptomic functional data at pathway levels
# Subset for stress ressponse, nitrogen and phosphorus metabolism


library("ggplot2")
library("dplyr")
library("tidyr")
library("readxl")



########################################################
# LEVEL 1 ANALYSIS - North compared to South
###########################################################


# level 1
data1<-read_xlsx("Function/Subsystems_DESeq_results_NvsS.xlsx", sheet="Subsystems_level-1_DESeq_result")
# level hierarchy
hierarchy<-read_xlsx("Function/Subsystems_DESeq_results_NvsS.xlsx", sheet="hierarchykey")


# add in some stars
data1$signif<-ifelse(data1$padj <=0.01, "**", 
                     ifelse(data1$padj <=0.05, "*",
                            ifelse("ns")))
data1sig<-filter(data1, padj<=0.05)
data1sig$type<-"level1"

# order from downreg to upreg
data1sigo<-data1sig[order(data1sig$log2FoldChange),]
order<-data1sigo$Level1


########################################################
# FIGURE 5A ---------------------------------------------

ggplot(data1sig, aes(Level1, log2FoldChange,fill=log2FoldChange))+
  geom_rect(data=NULL,aes(xmin=6.5,xmax=7.5,ymin=-Inf,ymax=Inf),
            fill="grey92", alpha=0.1)+ #grey92 is the ggplot default background!!
  scale_y_continuous(limits=c(-2.6,2.6))+
  scale_x_discrete(limits=order)+
  geom_col()+theme_minimal()+labs(x=NULL, y=NULL, fill="Fold Change")+
  theme(axis.text.y = element_text(size=14, color="grey40"), legend.position = "bottom",
        legend.text = element_text(size=10, color="grey20"), legend.key.width = unit(2, "cm"),
        legend.title = element_text(size=14, color="grey40"))+
  scale_fill_gradientn(limits=c(-2.5,2.5),colours=c("#0571b0","#92c5de","#f4a582","#ca0020"))+
  geom_errorbar(aes(ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), 
                size=0.8,width = 0, color="grey60") +
  geom_hline(yintercept=0, color="grey20")+coord_flip()





###########################################################
# STRESS RESPONSE LEVEL 2 - compared to mean
###########################################################

data<-read_excel("Function/DESeq_subsys_results_level2_summary.xlsx", sheet="selected")
datasig<-filter(data, Group=="Key")
datasigg<-gather(datasig, foldchange, foldvalue, "pvd_log2FoldChange","gb_log2FoldChange","bis_log2FoldChange",
                 "nar_log2FoldChange","nin_log2FoldChange")
datasiggp<-gather(datasig, padj, padjvalue, "pvd_padj","gb_padj","bis_padj","nar_padj","nin_padj")
datasigg$padj<-datasiggp$padj
datasigg$padjvalue<-datasiggp$padjvalue

datasigg$foldchange = factor(datasigg$foldchange, 
                             levels = c("pvd_log2FoldChange","gb_log2FoldChange","bis_log2FoldChange","nar_log2FoldChange","nin_log2FoldChange"))
datasigg$padjvalue<-as.numeric(datasigg$padjvalue)
datasigg$signif<-ifelse(datasigg$padjvalue <=0.01, "**", 
                        ifelse(datasigg$padjvalue <=0.05, "*",
                               ifelse("ns")))
datasiggpe<-gather(datasig, folderror, errorvalue, "pvd_lfcSE","gb_lfcSE","bis_lfcSE","nar_lfcSE","nin_lfcSE")
datasigg$folderror<-datasiggpe$datasiggpe
datasigg$errorvalue<-datasiggpe$errorvalue

datasigg$foldstar<-0
datasigg$foldstar[datasigg$foldvalue < 0] <- -0.1-datasigg$errorvalue
datasigg$foldstar[datasigg$foldvalue > 0] <- 0.05+datasigg$errorvalue


########################################################
## FIGURE 5B ----------------------------------
#2400x500
ggplot(datasigg, aes(Subsystem, foldvalue,fill=foldchange))+

  geom_col(position = "dodge")+
  geom_errorbar(aes(ymin=foldvalue-errorvalue, ymax=foldvalue+errorvalue), 
                size=0.8,position = position_dodge(width=0.9), width = 0, color="grey50")+
  labs(x=NULL, y=NULL)+scale_x_discrete(expand=c(0,0))+
  facet_grid(.~Subsystem, scales="free",space="free")+
  theme(axis.text.y = element_text(size=20, color="grey40"), 
        axis.text.x = element_blank(),strip.background = element_rect(fill="white", color="white"),
        strip.text = element_text(size=20, color="black"),
        legend.position = "none")+
  geom_hline(yintercept=0, color="navy")+
  scale_fill_manual(values=c("#253494","#0868ac","#43a2ca","#7bccc4","#bae4bc"))






#######################################################
# NUTRIENT METABOLISM LEVEL 3 - compared to mean
########################################################

data<-read_excel("Function/DESeq_subsys_results_level3_summary.xlsx")
datasigg<-gather(data, foldchange, foldvalue, "pvd_log2FoldChange","gb_log2FoldChange","bis_log2FoldChange",
                 "nar_log2FoldChange","nin_log2FoldChange")
datap<-gather(data, padj, padjvalue, "pvd_padj","gb_padj","bis_padj","nar_padj","nin_padj")
datasigg$padjvalue<-datap$padjvalue
datasigg$Level3<-datasigg$Subsystem

datasigg$padjvalue<-as.numeric(datasigg$padjvalue)
datasigg$signif<-ifelse(datasigg$padjvalue <=0.01, "**",
                              ifelse(datasigg$padjvalue <=0.05, "*",ifelse("ns")))

datasigg$foldchange = factor(datasigg$foldchange, 
                            levels = c("pvd_log2FoldChange","gb_log2FoldChange","bis_log2FoldChange","nar_log2FoldChange","nin_log2FoldChange"))

data3gpfpe<-gather(data, folderror, errorvalue, "pvd_lfcSE","gb_lfcSE","bis_lfcSE","nar_lfcSE","nin_lfcSE")
datasigg$folderror<-data3gpfpe$folderror
datasigg$errorvalue<-data3gpfpe$errorvalue



# level hierarchy
hierarchy<-read_excel("Function/Subsystems_DESeq_results_NvsS.xlsx", sheet="hierarchykey")
# add in hierarchy
hierarhcysimp3<-hierarchy
hierarhcysimp3$Level4<-NULL
hierarhcysimp3<-unique(hierarhcysimp3)
data3gpf<-full_join(datasigg, unique(hierarhcysimp3), by=c("Level3"))

data3gpfselect<-filter(data3gpf, 
                       Level1=="Phosphorus Metabolism" | Level1=="Nitrogen Metabolism")
data3gpfselect$Level3o<-factor(data3gpfselect$Level3, 
                           levels=c("Nitric_oxide_synthase",
                                    "Allantoin_Utilization",
                                    "Amidase_clustered_with_urea_and_nitrile_hydratase_functions",
                                    "Denitrification",
                                    "Dissimilatory_nitrite_reductase",
                                    "Putative_diaminopropionate_ammonia-lyase_cluster",
                                    "Nitrosative_stress",
                                    "Cyanate_hydrolysis",
                                    "Ammonia_assimilation",
                                    "Nitrate_and_nitrite_ammonification",
                                    "Alkylphosphonate_utilization",
                                    "Phosphate_metabolism",
                                    "Phosphoenolpyruvate_phosphomutase",
                                    "Phosphonate_metabolism"))


data3gpfselect$foldstar<-0
data3gpfselect$foldstar[data3gpfselect$foldvalue < 0] <- -0.1-data3gpfselect$errorvalue
data3gpfselect$foldstar[data3gpfselect$foldvalue > 0] <- 0.05+data3gpfselect$errorvalue





########################################################
## FIGURE 6 bottom ----------------------------------

sub<-
  ggplot(data3gpfselect, aes(Level3o, foldvalue))+theme_grey()+
  facet_grid(~Level1+foldchange, scales="free")+
  geom_col(aes(fill=Level3o),color="grey60", position = "dodge")+
  scale_fill_manual(values=c("#ffeda0","#fed976",
                             "#feb24c","#fd8d3c","#fc4e2a","firebrick1",
                             "#e31a1c","#bd0026","#800026",
                             "#c2a5cf","#9970ab","#762a83","#40004b"))+
  scale_color_manual(values=c("grey60","grey60",
                             "#feb24c","#fd8d3c","#fc4e2a","firebrick1",
                             "#e31a1c","#bd0026","#800026",
                             "#c2a5cf","#9970ab","#762a83","#40004b"))+
  geom_hline(yintercept=0, color="navy")+labs(x=NULL, y=NULL)+
  scale_y_continuous(limits=c(-2,2),labels = scales::number_format(accuracy=.1))+
  theme(axis.text.y = element_text(size=20, color="grey40"), 
        legend.position = "none",
        strip.background = element_blank(), strip.text = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_rect(color="grey80", fill=NA),
        axis.ticks.length = unit(0.2, "cm"), 
        axis.ticks.y = element_line(color="grey60"),
        axis.ticks.x = element_blank())+
  geom_errorbar(aes(ymin=foldvalue-errorvalue, ymax=foldvalue+errorvalue, color=Level3o), 
                size=0.8,position = position_dodge(width=0.9), width = 0)




#################################################
# Nitrogen and Phosphorous only level 2 ---------------------------------
#################################################

# level 2
data<-read_excel("Function/DESeq_subsys_results_level2_summary.xlsx", sheet=1)

datasigg<-gather(data, foldchange, foldvalue, "pvd_log2FoldChange","gb_log2FoldChange","bis_log2FoldChange",
                 "nar_log2FoldChange","nin_log2FoldChange")
datasiggp<-gather(data, padj, padjvalue, "pvd_padj","gb_padj","bis_padj","nar_padj","nin_padj")
datasigg$padjvalue<-datasiggp$padjvalue

datasigg$foldchange = factor(datasigg$foldchange, 
                             levels = c("pvd_log2FoldChange","gb_log2FoldChange","bis_log2FoldChange","nar_log2FoldChange","nin_log2FoldChange"))
datasigg$padjvalue<-as.numeric(datasigg$padjvalue)
datasigg$signif<-ifelse(datasigg$padjvalue <=0.01, "**", 
                        ifelse(datasigg$padjvalue <=0.05, "*",
                               ifelse("ns")))


data3gpfpe<-gather(data, folderror, errorvalue, "pvd_lfcSE","gb_lfcSE","bis_lfcSE","nar_lfcSE","nin_lfcSE")
datasigg$folderror<-data3gpfpe$folderror
datasigg$errorvalue<-data3gpfpe$errorvalue



datasiggselect<-filter(datasigg, 
                       Subsystem=="Phosphorus Metabolism" | Subsystem=="Nitrogen Metabolism")

########################################################
## FIGURE 6 top ----------------------------------
nitphos<- ggplot(datasiggselect, aes(Subsystem, foldvalue,fill=foldchange))+
  geom_col(color="grey60", position = "dodge")+
  scale_x_discrete(expand=c(0,0))+theme_minimal()+
  scale_y_continuous(labels = scales::number_format(accuracy=.1))+
  labs(x=NULL, y=NULL)+facet_grid(.~Subsystem+foldchange, scales="free",space="free")+
  theme(axis.text.y = element_text(size=20, color="grey40"), 
        axis.text.x = element_blank(),strip.background = element_rect(fill="white"),
        strip.text = element_blank(),axis.ticks.x = element_blank(),
        legend.position = "none", axis.ticks.length = unit(0.2, "cm"), axis.ticks.y = element_line(color="grey60"))+
  geom_hline(yintercept=0, color="navy")+
  scale_fill_manual(values=c("#253494","#0868ac","#43a2ca","#7bccc4","#bae4bc"))+
  geom_errorbar(aes(ymin=foldvalue-errorvalue, ymax=foldvalue+errorvalue), 
                size=0.8,position = position_dodge(width=0.9), width = 0, color="grey60")


########################################################
## FIGURE 6 all together ----------------------------------
cowplot::plot_grid(nitphos,sub,nrow=2, align="hv", axis="lr", rel_heights = c(2,5))
#1500x800
# add legend and labels in Inkscape


