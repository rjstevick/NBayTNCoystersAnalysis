# Stevick et al 2020 Oyster Gut Microbiome Function in an Estuary
# Map and environmental data
# Figures 1 & 2 and Table 1
# updated 20200903


# Figure 1 Map --------------------------------------

# load mapping packages
library(tidyverse)
library(sf)

# import metadata averaged per site and lat/longs for mapping
data <- readxl::read_xlsx("Metadata/EnvironmentalPCAmetadata.xlsx", sheet="SiteSummary")
# load RI coastline
outline <- st_read("Metadata/NSDE66796/CUSPLine.shp")

# Plot map and site points
ggplot() + theme_bw()+
  geom_sf(data = outline, fill= "grey20") +
  geom_point(data, mapping = aes(x = Long, y = Lat, fill = Station), size=6, shape=21)+
  scale_fill_manual(values=c("#253494","#0868ac","#43a2ca","#7bccc4","#bae4bc"))+ 
  coord_sf(xlim = c(-71.8,-71.0), ylim = c(41.2,41.9), expand = FALSE) + 
  theme(panel.background = element_rect(fill = "white"))+
  labs(x=NULL, y=NULL)+theme(axis.text = element_text(size=10),legend.position = "none")


# Figure 2 PCA -----------------------------------------

# load packages
library(gridExtra)
library(ggfortify)

# clean up data per site
pcadata<-as.data.frame(data)
rownames(pcadata)<-c(pcadata[,1]) # copy site names to rownames
pcadatahisto<-pcadata[,c(4:12)] # remove lat/long from PCA

# Calculate PCA 
pca_result <- prcomp(pcadatahisto, scale=TRUE)
pca_result$rotation <- -pca_result$rotation # rotate PC axes
pca_result$rotation # check on variable

pca_result$x <- - pca_result$x # rotate x-axis
autoplot(pca_result) # plot basic result

theme_set(theme_grey())
# generate basic plot
aplot<-autoplot(pca_result,
                size=8,  # size of points
                data=pcadata, # identify metadata to use for plotting
                fill="Station",shape="Station", # change point fill and shape by site
                loadings=TRUE, # include arrows
                loadings.label=TRUE, # include arrow labels
                loadings.label.label=c("Salinity (ppt)", "Temperature (\u00B0C)",
                                       "Ammonium (\U003BCM)","Phosphate (\U003BCM)", "Nitrite (\U003BCM)","Nitrate (\U003BCM)",
                                       "pH","Chlorophyll-a (\U003BCM)","Dissolved Oxygen \n(mg/L)\n"),
                loadings.label.fontface="bold",
                loadings.label.repel=TRUE, # don't overlap arrow labels
                loadings.label.size=3.5, 
                loadings.label.colour="black",
                loadings.colour=c("orange","orange",
                                  "cornflowerblue","cornflowerblue","cornflowerblue","cornflowerblue",
                                  "orange","orange","orange"))
# adjust labels, colors, point shapes and plot
aplot + 
  geom_text(aes(label=Station), nudge_x=0.03, nudge_y=-0.06, color="grey20") + 
  scale_fill_manual(values=c("#253494","#0868ac","#43a2ca","#7bccc4","#bae4bc"))+
  scale_shape_manual(values = c(21,22,23,24,25))+
  scale_y_continuous(limits=c(-.65, .65)) + scale_x_continuous(limits=c(-.7, .7))


# Table 1 correlations -----------------------------------------

data %>%
  pivot_longer("Salinity_ppt":"DO_mgL", names_to="parameter") %>% 
  group_by(parameter) %>% 
  summarise(cor(Lat, value), method=c("pearson", "spearman"))
# note: pearson and spearman give the same results
