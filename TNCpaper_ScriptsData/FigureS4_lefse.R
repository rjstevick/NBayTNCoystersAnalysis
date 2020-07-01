# Stevick et al 2020 Oyster Gut Microbiome Function in an Estuary
# 16S abundances using LefSe
# Figure S4
# updated 4/16/2020 for resubmission

# 16S Amplicon data at order level analyzed with LefSe
# Access from here: https://huttenhower.sph.harvard.edu/galaxy/


library(tidyr)
library(dplyr)
library(ggplot2)
library(forcats)

#import data - summarized outputs from Galaxy
datatype<-read_xlsx("Taxonomy/TNC_summaryLefSeResults.xlsx", sheet="Types_gutwater")
datasites<-read_xlsx("Taxonomy/TNC_summaryLefSeResults.xlsx", sheet="Sites")

count(datatype$Type)

# Plot for sample types
tplot<- datatype %>%
  mutate(Order=fct_reorder(Order, LDA_Score_log10)) %>%
  ggplot( aes(Order,LDA_Score_log10,fill=Type))+
  geom_bar(stat="identity")+coord_flip()+
  facet_grid(Type~., scales="free",space="free")+
  theme(legend.position = "none")+
  labs(y="LDA Score (log 10)",x=NULL)+
  ggtitle("Types")+
  scale_fill_manual(values=c("orange", "darkred"))

# Plot for sites
splot<- datasites %>%
  mutate(Order=fct_reorder(Order, LDA_Score_log10)) %>%
  ggplot( aes(Order,LDA_Score_log10,fill=Site))+
  geom_bar(stat="identity")+coord_flip()+
  facet_grid(Site~., scales="free",space="free")+
  theme(legend.position = "none")+
  labs(y="LDA Score (log 10)",x=NULL)+
  ggtitle("Sites")+
  scale_fill_manual(values=c("#253494","#0868ac","#43a2ca","#7bccc4","#bae4bc"))


# Plot together - FIGURE S4
cowplot::plot_grid(tplot,splot, ncol=2, align="hv", axis="tbr")
#1300x800




# Just top 30 taxa --------------------
# show only if they are significantly different

# get this variable from script for Figure 4
toptax

# Extact data for orders included in toptax
topdatatype <- datatype[datatype$Order %in% toptax$Order,]
topdatasites <- datasites[datasites$Order %in% toptax$Order,]


# Plot for sample types
tplot<- topdatatype %>%
  mutate(Order=fct_reorder(Order, LDA_Score_log10)) %>%
  ggplot( aes(Order,LDA_Score_log10,fill=Type))+
  geom_bar(stat="identity")+coord_flip()+
  facet_grid(Type~., scales="free",space="free")+
  theme(legend.position = "none")+
  labs(y="LDA Score (log 10)",x=NULL)+
  ggtitle("Types")+
  scale_fill_manual(values=c("orange", "darkred"))

# Plot for sites
splot<- topdatasites %>%
  mutate(Order=fct_reorder(Order, LDA_Score_log10)) %>%
  ggplot( aes(Order,LDA_Score_log10,fill=Site))+
  geom_bar(stat="identity")+coord_flip()+
  facet_grid(Site~., scales="free",space="free")+
  theme(legend.position = "none")+
  labs(y="LDA Score (log 10)",x=NULL)+
  ggtitle("Sites")+
  scale_fill_manual(values=c("#253494","#0868ac","#43a2ca","#7bccc4","#bae4bc"))
