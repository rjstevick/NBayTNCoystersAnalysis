# Stevick et al 2020 Oyster Gut Microbiome Function in an Estuary
# 16S diversity and rarefaction plots
# Figures 3 and S1
# updated 4/16/2020 for resubmission

# 16S Amplicon data at ASV level


library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(phyloseq)
library(vegan)
library(readxl)
library(ggpubr)
library(ggbiplot)
library(vegan)
library(ggpubr)
library(stringr)


#import data - ASV counts, **not normalized**
data<-read_xlsx("Taxonomy/16SallData_SILVAtaxa.xlsx", sheet="ASVcounts")
taxakey<-read_xlsx("Taxonomy/16SallData_SILVAtaxa.xlsx", sheet="Taxonomy")
metadata<-read_xlsx("Taxonomy/16SallData_SILVAtaxa.xlsx", sheet="Metadata")

#add in all taxa hierarchy to ASV counts
datafull<-full_join(data, taxakey)
#convert to long form
datafulllong<-gather(datafull, SampleID, count, "TNC01":"TNC64") 


####
### Reads per sample including chloroplast reads ----------------------------
####


# Add in percent-normalized data
# Determine number of reads per sample
SampleReadCounts<-
  datafulllong %>%
  group_by(SampleID) %>%
  dplyr::summarise(SampleTotalReads=(sum(count)))
# add read counts in to main table
metadata$SampleTotalReads<-SampleReadCounts$SampleTotalReads
datafulllong<-
  full_join(datafulllong, SampleReadCounts, by="SampleID")
# calculate percent abundances
datafulllong$percent<-
  datafulllong$count/datafulllong$SampleTotalReads
#add in the metadata
datafulllongmeta<-full_join(datafulllong, metadata, by = "SampleID",
                            copy=FALSE, suffix=c(".x",".y"))

#select only gut, water
metadatags<-filter(metadata, SampleType=="gut" | SampleType=="water")

# Number of reads per sample
ggplot(metadatags,aes(SampleName,SampleTotalReads,fill=Station))+geom_bar(stat="identity")+
  facet_grid(.~SampleType+Station, scales="free",space="free")+
  theme(legend.position = "none", axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.length = unit(0.6, 'lines'))+
  labs(y="Number of reads \nper sample",x=NULL)+
  scale_fill_manual(values=c("#253494","#0868ac","#43a2ca","#7bccc4","#bae4bc"), 
                    labels=c("1. Providence River", "2. Greenwich Bay", "3. Bissel Cove", "4. Narrow River", "5. Ninigret Pond"))



####
### 16S Rarefaction curve and coverage - Figure S1 ----------------------------------------------
####

datamat<-data[1:61] # remove controls
datamat$ASVID<-NULL
datamatt<-data.table::transpose(datamat)
raremax <- min(rowSums(datamatt))
S <- specnumber(datamatt) # observed number of ASVs
Srare <- rarefy(datamatt, raremax)


#define colors based on gut/site or water
colors<-c("#253494","#253494","#253494","#253494","#253494","#253494","#253494","#253494","#253494","#253494",
          "#0868ac","#0868ac","#0868ac","#0868ac","#0868ac","#0868ac","#0868ac","#0868ac","#0868ac","#0868ac",
          "#43a2ca","#43a2ca","#43a2ca","#43a2ca","#43a2ca","#43a2ca","#43a2ca","#43a2ca","#43a2ca","#43a2ca",
          "#7bccc4","#7bccc4","#7bccc4","#7bccc4","#7bccc4","#7bccc4","#7bccc4","#7bccc4","#7bccc4","#7bccc4",
          "#bae4bc","#bae4bc","#bae4bc","#bae4bc","#bae4bc","#bae4bc","#bae4bc","#bae4bc","#bae4bc","#bae4bc",
          "darkred","darkred","darkred","darkred","darkred","darkred","darkred","darkred","darkred","darkred")
#put plotting parameters into pars
pars <- expand.grid(col = colors, stringsAsFactors = FALSE)


#plot number of ASVs vs. rarefied
#550x400
plot(S, Srare, pch=21, cex=2, col="black", bg=colors,
     xlab = "Observed Number of ASVs", ylab = "Rarefied Number of ASVs")
abline(0, 1, lty="dashed")
legend(700,250, c("1.PVD gut","2.GB gut","3.BIS gut","4.NAR gut","5.NIN gut", "All water"),
       cex=0.8,pch=19,
       text.width = strwidth("60000000"),xjust = 1, yjust = 1,
       col=c("#253494","#0868ac","#43a2ca","#7bccc4","#bae4bc","darkred"))

#generate rarefaction curves
#NOTE: THIS STEP TAKES A LITTLE WHILE TO RUN.
outw <- with(pars, rarecurve(datamatt, step = 20, 
                             sample = raremax, col = colors, label = FALSE))
#determine limits of the graph
Smax <- sapply(outw, max)
#plot the rarefaction curves
#550x400
plot(c(1, 80000), c(1, max(Smax)), xlab = "Number of Reads",
     ylab = "Observed Number of ASVs", type = "n")
abline(v = raremax, lty="dashed")
for (i in seq_along(outw)) {
  N <- attr(outw[[i]], "Subsample")
  with(pars, lines(N, outw[[i]], col = col[i]))
}

# check slopes of curves at the end, make sure they're ~0
end<-raremax-1
slopes<-rareslope(datamatt, end)
plot(slopes, type="p", pch=21, bg = colors,
     cex=2, xlab = "16S Amplicon Samples", ylim=c(0,0.05))
# coverage proxy:
coverage<-100-100*rareslope(datamatt, end)
plot(coverage, type="p", pch=21, bg = colors,
     cex=2, xlab = "16S Amplicon Samples", 
     ylab="Coverage (%)", ylim=c(95,100))

mean(coverage)
sd(coverage)



### ASV level without chloro reads---------------------------------------------

####
### Remove chloroplast reads and recalculate percent abundances ------------
####

# remove chloroplast reads
datafulllong <- datafulllong %>%
  filter(str_detect(Taxon, "Chloroplast", negate=TRUE)) 

# calculate new number of reads per sample
SampleReadCounts<-
  datafulllong %>%
  group_by(SampleID) %>%
  dplyr::summarise(SampleTotalReads=(sum(count)))
# add read counts in to main table
metadata$SampleTotalReads<-SampleReadCounts$SampleTotalReads
datafulllong<-
  full_join(datafulllong, SampleReadCounts, by="SampleID")
# calculate percent abundances
datafulllong$percent<-
  datafulllong$count/datafulllong$SampleTotalReads
#add in the metadata
datafulllongmeta<-full_join(datafulllong, metadata, by = "SampleID",
                            copy=FALSE, suffix=c(".x",".y"))
#select only gut, water
dataASVmetagw<-filter(datafulllongmeta, SampleType=="gut" | SampleType=="water")




####
### Alpha diversity Plotting -----------------------------------------------------
####

#convert normalized data (with chloros removed) back to wide
dataASVwide<-select(dataASVmetagw, SampleID, ASVID, percent)
dataASVtable<-spread(dataASVwide, ASVID, percent) %>%
  tibble::column_to_rownames("SampleID") 

# calculate diversity
diversitytotal<-diversity(dataASVtable, index="simpson")
# add into metadata variable
metadatags$Simpsons<-diversitytotal

# make sub dataframes for stats if needed.
dataag_gut<-filter(metadatags, SampleType=="gut")
dataag_water<-filter(metadatags, SampleType=="water")


####
# Simpson's diversity plot - S1A
divplot<-ggplot(metadatags,aes(x=Station,y=Simpsons, fill=Station))+
  geom_boxplot()+
  geom_point(size=4, shape=23)+
  facet_grid(~SampleType, labeller=as_labeller(c(gut="Gut samples (n=10)", water="Water samples (n=2)"))) + 
  labs(x=NULL,y="Simpson's Index of Diversity",fill="Site")+
  scale_color_manual(values=c("#253494","#0868ac","#43a2ca","#7bccc4","#bae4bc")) +
  scale_fill_manual(values=c("#253494","#0868ac","#43a2ca","#7bccc4","#bae4bc")) + 
  #                    labels=c("1. Providence River", "2. Greenwich Bay", "3. Bissel Cove", "4. Narrow River", "5. Ninigret Pond"))+
  theme_bw()+scale_y_continuous(limits=c(0,1.0))+
  theme(legend.position = "none",
        strip.background = element_rect("grey90"))

# generate legend for plotting
divlegend<-ggplot(metadatags,aes(x=Station,y=Simpsons, fill=Station))+
  geom_point() + facet_grid(~SampleType) + theme_bw()+
  labs(x=NULL,y="Simpson's Index of Diversity",fill="Site")+
  scale_color_manual(values=c("#253494","#0868ac","#43a2ca","#7bccc4","#bae4bc")) +
  scale_fill_manual(values=c("#253494","#0868ac","#43a2ca","#7bccc4","#bae4bc")) + 
  theme(legend.position = c(.86,.33),
        legend.text = element_text(size=14, colour="gray20"), 
        strip.background = element_rect("grey90"),
        legend.background = element_rect(color="grey40"))




####
### Alpha diversity stats -----------------------------------------------------
####


compare_means(data=dataag_water, Simpsons ~ Station, method="kruskal")
kruskal.test(dataag_water$Simpsons, as.factor(dataag_water$Station))
# From the output of the Kruskal-Wallis test, we know that there is a significant difference between groups, but we don't know which pairs of groups are different.
# It's possible to use the function pairwise.wilcox.test() to calculate pairwise comparisons between group levels with corrections for multiple testing.
pairwise.wilcox.test(dataag_water$Simpsons, dataag_water$Station,
                     p.adjust.method = "BH")


compare_means(data=dataag_gut, Simpsons ~ Station, method="kruskal")
kruskal.test(dataag_gut$Simpsons, as.factor(dataag_gut$Station))
# From the output of the Kruskal-Wallis test, we know that there is a significant difference between groups, but we don't know which pairs of groups are different.
# It's possible to use the function pairwise.wilcox.test() to calculate pairwise comparisons between group levels with corrections for multiple testing.
pairwise.wilcox.test(dataag_gut$Simpsons, dataag_gut$Station,
                     p.adjust.method = "BH")

# all together now!
compare_means(data=metadatags, Simpsons~Station, group.by="SampleType",
              method="wilcox", p.adjust.method = "BH")




####
### Beta diversity: NMDS Bray-Curtis----------------------------------------
####

# ellipse function
theme_set(theme_bw())
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}



# Get spread of gut points based on Station
abund_table<-dataASVtable[c(1:50),]
meta_table<-dataag_gut

grouping_info<-data.frame(meta_table$SampleType, meta_table$Station, meta_table$OysterNum)
sol<-metaMDS(abund_table,distance = "bray", k = 2, trymax = 50)
stressplot(sol) # R2 = 0.964
NMDS=data.frame(x=sol$point[,1],y=sol$point[,2],Type=as.factor(grouping_info[,1]),Station=as.factor(grouping_info[,2]),OysterNum=as.factor(grouping_info[,3]))
plot.new()
ord<-ordiellipse(sol, as.factor(grouping_info[,2]), display = "sites", kind ="sd", conf = 0.95, label = T)
dev.off()
df_ell <- data.frame()
for(g in levels(NMDS$Station)){
  if(g!="" && (g %in% names(ord))){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$Station==g,],
                                                     veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale))),Station=g))}}
head(df_ell)
NMDS.mean=aggregate(NMDS[,1:2],list(group=NMDS$Station),mean)
#Now do the actual plotting
gutNMDS<-ggplot(data=NMDS,aes(x,y,colour=Station,fill=Station))+theme_bw() + 
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2, lty=Station), size=1) + 
  geom_point(size=4, alpha=0.9,aes(shape=Station))+scale_shape_manual(values = c(21,22,23,24,25))+
  annotate("text",x=NMDS.mean$x,y=NMDS.mean$y,label=NMDS.mean$group,size=5, color="gray40") + 
  #  annotate("text",x=NMDS$x+0.2,y=NMDS$y+0.05,label=paste(NMDS$Station,NMDS$OysterNum),size=4, color="gray40") + # to determine which point is which 
  scale_fill_manual(values=c("#253494","#0868ac","#43a2ca","#7bccc4","#bae4bc"))+
  scale_colour_manual(values=c("#253494","#0868ac","#43a2ca","#7bccc4","#bae4bc"))+
  scale_linetype_manual(values=c("solid","dotted","twodash","longdash", "solid"), labels=c("1. Providence River", "2. Greenwich  Bay", "3. Bissel Cove", "4. Narrow River", "5. Ninigret Pond"))+
  theme(legend.text = element_text(size=14, colour="gray20"), legend.position = "right",
        legend.title = element_blank(),legend.box="horizontal")

adonis2(abund_table~Station, data=NMDS, by=NULL,method="bray", k=2)



# Get spread of points based on Type
abund_table<-dataASVtable[c(1:60),]
meta_table<-metadatags

grouping_info<-data.frame(meta_table$SampleType, meta_table$Station, meta_table$OysterNum)
sol<-metaMDS(abund_table,distance = "bray", k = 2, trymax = 50)
stressplot(sol) # R2 = 0.964
NMDS=data.frame(x=sol$point[,1],y=sol$point[,2],Type=as.factor(grouping_info[,1]),
                Station=as.factor(grouping_info[,2]),OysterNum=as.factor(grouping_info[,3]))
plot.new()
ord<-ordiellipse(sol, as.factor(grouping_info[,1]), display = "sites", kind ="sd", conf = 0.95, label = T)
dev.off()
df_ell <- data.frame()
for(g in levels(NMDS$Type)){
  if(g!="" && (g %in% names(ord))){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$Type==g,],
                                                     veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale))),Type=g))}}
head(df_ell)
NMDS.mean=aggregate(NMDS[,1:2],list(group=NMDS$Type),mean)
head(NMDS.mean)
#Now do the actual plotting
typeNMDS<-ggplot(data=NMDS,aes(x,y,colour=Type,fill=Type))+theme_bw() + 
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2, lty=Type), size=1) + 
  geom_point(size=4, alpha=0.9,aes(shape=Station))+
  annotate("text",x=NMDS.mean$x,y=NMDS.mean$y,label=NMDS.mean$group,size=6, color="gray10") + 
  scale_fill_manual(values=c("orange","darkred"), labels=c("Gut","Water"))+
  scale_colour_manual(values=c("orange","darkred"), labels=c("Gut","Water"))+
  scale_linetype_manual(values=c("solid","twodash"), labels=c("Gut","Water"))+
  scale_shape_manual(values = c(21,22,23,24,25))+
  # , labels=c("1. Providence River", "2. Greenwich  Bay", "3. Bissel Cove", "4. Narrow River", "5. Ninigret Pond"))+
  theme(legend.text = element_text(size=14, colour="gray20"), legend.position = "right",
        legend.title = element_blank())

adonis2(abund_table~Type, data=NMDS, by=NULL,method="bray", k=2)




####
### Build Figure 3 with cowplots ------------------------------------------
####

gutleg<-cowplot::get_legend(gutNMDS+scale_linetype_discrete(guide = FALSE))
typeleg<-cowplot::get_legend(typeNMDS+scale_shape_discrete(guide = FALSE))
divleg<-cowplot::get_legend(divlegend+theme(legend.title = element_blank(), legend.position = "right", legend.background = element_rect(color="white")))
nmds<-cowplot::plot_grid(gutNMDS+theme(legend.position = "none"),
                         typeNMDS+theme(legend.position = "none"),
                         align="hv", axis="t",nrow=1)

plots<-cowplot::plot_grid(divplot+theme(legend.position = "none"),nmds, nrow=2, 
                          align="v", axis="l",rel_heights = c(40,60))
legends<-cowplot::plot_grid(divleg,typeleg,gutleg, ncol = 1,rel_heights = c(50,20,50))


cowplot::plot_grid(plots,legends, ncol=2, rel_widths = c(80,20), align="hv", axis="tbr")
# 900x550
# label y as "Bray-Curtis beta diversity (k=2)" in Inkscape
