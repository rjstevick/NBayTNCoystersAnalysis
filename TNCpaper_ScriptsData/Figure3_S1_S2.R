# Stevick et al 2020 Oyster Gut Microbiome Function in an Estuary
# 16S controls, rarefaction, and diversity plots
# Figures 3, S1, and S2
# RJS updated 6/30/2020 for resubmission

# 16S Amplicon data at ASV level

library(ggbiplot)
library(tidyverse)
library(RColorBrewer)
library(phyloseq)
library(vegan)
library(readxl)
library(ggpubr)
library(vegan)
library(stringr)
library(scales)

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
  full_join(datafulllong, SampleReadCounts, by="SampleID") %>%
  # calculate percent abundances
  mutate(percent=count/SampleTotalReads)
#add in the metadata
datafulllongmeta<-full_join(datafulllong, metadata, by = "SampleID",
                            copy=FALSE, suffix=c(".x",".y"))

#select only gut, water
metadatags<-filter(metadata, SampleType=="gut" | SampleType=="water")

# Number of QCd sequencing reads per sample
ggplot(metadata,aes(SampleName,SampleTotalReads,fill=Station))+
  geom_bar(stat="identity")+
  facet_grid(.~SampleType+Station, scales="free",space="free")+
  theme(legend.position = "none", axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.length = unit(0.6, 'lines'))+
  labs(y="Number of reads \nper sample",x=NULL)+
  scale_fill_manual(values=c("#253494","#0868ac","#43a2ca","#7bccc4","#bae4bc"),
                    labels=c("1. Providence River", "2. Greenwich Bay", "3. Bissel Cove", "4. Narrow River", "5. Ninigret Pond"))



####
### 16S Rarefaction curve and coverage - Figure S2 ----------------------------------------------
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
#put plotting colors into pars
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
#plot the rarefaction curves - Figure S2A
#550x400
par(mfrow=c(3,1))
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
# Figure S2B
plot(slopes, type="p", pch=21, bg = colors,
     cex=2, xlab = "16S Amplicon Samples", 
     ylab="Slopes at Raremax", ylim=c(0,0.05))
# coverage proxy:
coverage<-100-100*rareslope(datamatt, end)
# Figure S2C
plot(coverage, type="p", pch=21, bg = colors,
     cex=2, xlab = "16S Amplicon Samples",
     ylab="Estimated Coverage (%)", ylim=c(95,100))

mean(coverage)
sd(coverage)



### ASV level without chloro reads---------------------------------------------

####
### Remove chloroplast reads and recalculate percent abundances ------------
####

datafulllong<-gather(datafull, SampleID, count, "TNC01":"TNC64")

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
  full_join(datafulllong, SampleReadCounts, by="SampleID") %>%
  # calculate percent abundances
  mutate(percent=count/SampleTotalReads)

#add in the metadata
datafulllongmeta<-full_join(datafulllong, metadata, by = "SampleID",
                            copy=FALSE, suffix=c(".x",".y"))
#select only gut, water
dataASVmetagw<-filter(datafulllongmeta, SampleType=="gut" | SampleType=="water")


####
### Bar plots with Negative & Positive Controls - Figure S1 -----------------------------------------------------
####

# sum ASV percentages by phylum
phylumperc<- datafulllongmeta %>%
  group_by(Phylum,SampleID, SampleName, SampleType, Station, TypeStationGroup) %>%
  dplyr::summarise(physum=sum(percent)) %>%
  filter(SampleName != "NEG.CON" & SampleName != "NEG.CON.2") %>% ungroup() %>% 
  # show only the top 10 phyla, put the rest in "Other"
  mutate(PhylumOther=forcats::fct_lump(f=Phylum, w=physum, other_level="Others", n=10))

# define palette
palette<-c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
           "#FDBF6F", "#FF7F00", "#CAB2D6","gray48","gray30")

barp<-ggplot(phylumperc,
       (aes(x=SampleName, y=physum,
            fill=factor(PhylumOther, levels=c("Actinobacteria",  "Bacteroidetes","Chloroflexi","Cyanobacteria","Firmicutes","Fusobacteria","Planctomycetes","Proteobacteria", "Tenericutes","Verrucomicrobia","Unknown","Others")))))+
  facet_grid(.~TypeStationGroup, scales="free", space="free")+
  geom_col(position="fill", alpha=0.8)+theme_minimal()+
  theme(legend.text = element_text(size=12, colour="gray20"),
        legend.position = "right",
        axis.text.x = element_text(angle=60, hjust=1, vjust=1),
        axis.ticks = element_line(inherit.blank=FALSE, color="grey30"))+
  facet_grid(.~SampleType+Station, scales="free",space="free")+
  scale_fill_manual(values=c(palette))+
  labs(y="Percent abundance",x=NULL,fill="Phylum")+
  scale_y_continuous(labels = scales::percent, expand=c(0,0))


# mock control samples ASV percentages
controlASVperc<- datafulllongmeta %>%
  filter(SampleName == "MOCK.CON" | SampleName == "MOCK.CON.2") %>%
  group_by(ASVID, Taxon, SampleName) %>%
  dplyr::summarise(Taxonsum=sum(percent)) %>% ungroup() %>% 
  mutate(TaxonOther=forcats::fct_lump(f=Taxon, w=Taxonsum, other_level="Others", n=10))

palette2<-c("#8c48d6",
            "#5e9850",
            "#9bae34",
            "#c04fae",
            "#767ba7",
            "#b74f6f",
            "#cc5d32",
            "#906e49",
            "#509f8d",
            "#6259b2",'gray30')

# plot controls at ASV level
ctrlp<-ggplot(dplyr::arrange(controlASVperc,TaxonOther), (aes(x=SampleName, y=Taxonsum, fill=TaxonOther)))+
  geom_col(position="fill", color="white")+theme_minimal()+
  coord_flip()+
  theme(legend.text = element_text(size=9, colour="gray20"),
        legend.position = "bottom",
        axis.ticks = element_line(inherit.blank=FALSE, color="grey30"))+
  scale_fill_manual(values=c(palette2))+
  labs(y="Percent abundance",x=NULL,fill=NULL)+
  scale_y_continuous(labels = scales::label_percent(accuracy=1), breaks=breaks_pretty(n=10), expand=c(0,0))+
  guides(fill = guide_legend(ncol = 2))


# Figure S1
cowplot::plot_grid(barp, ctrlp, nrow=2, rel_heights = c(60,25), labels=c("A","B"))
ggsave("FigureS1.png", width = 15, height = 10, dpi=400)
ggsave("FigureS1.pdf", width = 15, height = 10)
ggsave("FigureS1.svg", width = 15, height = 10)





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
sol$stress
stressplot(sol)
# stress=0.19, R2=0.96

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
sol$stress
stressplot(sol)
# stress=0.19, R2=0.96

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





####
### Core Microbiome - Figure XXX ----------------------------------------------
####

coretaxa<-dataASVwide %>%
  left_join(metadatags) %>% 
  # select only the gut samples, not water
  filter(SampleType=="gut") %>% 
  # remove when the ASV is absent
  filter(percent!=0) %>%
  # group by ASV and the count the number of samples it occurs in
  group_by(ASVID) %>% count() %>% 
  # keep ASVs that occur in at least 80% of samples
  arrange(desc(n)) %>% filter(n>=40) %>% 
  # add back in the other metadata
  left_join(taxakey) %>% 
  left_join(dataASVwide) %>% 
  left_join(metadatags) %>% 
  filter(SampleType=="gut") %>% 
  # make a binary scale for when the ASV is present in each sample
  mutate(presence=case_when(percent==0 ~ "0", TRUE ~ "1")) %>%
  mutate(taxonlabel=paste(Phylum,Class,Order,Family,Genus,Species,sep="; "))

ggplot(coretaxa, 
       aes(x=SampleID, y=reorder(taxonlabel,n), fill=presence))+
  geom_tile(color="white")+
  facet_grid(.~Station,scales="free", space="free")+
  scale_fill_manual(values=c("grey80","aquamarine4"))+
  scale_y_discrete(labels=wrap_format(50))+
  theme_minimal()+
  theme(legend.position = "none", axis.text.x=element_blank())+
  labs(x=NULL, y=NULL)

