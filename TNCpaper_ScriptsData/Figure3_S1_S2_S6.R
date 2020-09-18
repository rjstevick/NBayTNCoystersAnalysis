# Stevick et al 2020 Oyster Gut Microbiome Function in an Estuary
# 16S controls, rarefaction, and diversity plots
# Figures 3, S1, S2, and S6. Tables S2 and S3
# RJS updated 09/18/2020 for reresubmission

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
library(cowplot)

#import data - ASV counts, **not normalized**
data<-read_xlsx("Taxonomy/16SallData_SILVAtaxa.xlsx", sheet="ASVcounts")
taxakey<-read_xlsx("Taxonomy/16SallData_SILVAtaxa.xlsx", sheet="Taxonomy")
metadata<-read_xlsx("Taxonomy/16SallData_SILVAtaxa.xlsx", sheet="Metadata")

#add in all taxa hierarchy to ASV counts
datafull<-full_join(data, taxakey)
#convert to long form
datafulllong<-gather(datafull, SampleID, count, "TNC01":"TNC64")


####
### + Reads per sample including chloroplast reads ----------------------------
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

# Number of QCd sequencing reads per sample
metadata %>% filter(SampleType=="gut" | SampleType=="water") %>% 
  ggplot(aes(SampleName,SampleTotalReads,fill=Station))+
  geom_bar(stat="identity")+
  facet_grid(.~SampleType+Station, scales="free",space="free")+
  theme(legend.position = "none", axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.length = unit(0.6, 'lines'))+
  labs(y="Number of reads \nper sample",x=NULL)+
  scale_fill_manual(values=c("#253494","#0868ac","#43a2ca","#7bccc4","#bae4bc"),
                    labels=c("1. Providence River", "2. Greenwich Bay", "3. Bissel Cove", "4. Narrow River", "5. Ninigret Pond"))



####
### Figure S2. 16S Rarefaction curve and coverage ----------------------------------------------
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



####
### + ASV level without chloroplast reads -----------------------------------------------
####

# Remove chloroplast reads and recalculate percent abundances
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

#convert normalized data (with chloros removed) back to wide
dataASVwide<-select(dataASVmetagw, SampleID, ASVID, percent)
dataASVtable<-spread(dataASVwide, ASVID, percent) %>%
  tibble::column_to_rownames("SampleID")

####
### Figure S1. Bar plots with Negative & Positive Controls -----------------------------------
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
            fill=factor(PhylumOther, levels=c("Actinobacteria",  "Bacteroidetes","Chloroflexi","Cyanobacteria","Firmicutes","Planctomycetes","Proteobacteria", "Tenericutes","Verrucomicrobia","Unknown","Others")))))+
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
#ggsave("FigureS1.png", width = 15, height = 10, dpi=400)
#ggsave("FigureS1.pdf", width = 15, height = 10)
#ggsave("FigureS1.svg", width = 15, height = 10)





####
### Figure 3A. Alpha diversity -----------------------------------------------------
####


#select only gut, water
metadatags<-filter(metadata, SampleType=="gut" | SampleType=="water")
# calculate diversity
diversitytotal<-diversity(dataASVtable, index="simpson")
# add into metadata variable
metadatags$Simpsons<-diversitytotal




#### Simpson's diversity plots
gutdiv<-metadatags %>% filter(SampleType=="gut") %>% 
  ggplot(aes(x=Station,y=Simpsons, fill=Station))+
  geom_jitter(width=0.15, size=3, shape=23, alpha=0.8)+
  geom_boxplot(alpha=0.8)+
  labs(title="Gut samples (n=10)", x=NULL,y="Simpson's Index of Diversity",fill="Site")+
  scale_color_manual(values=c("#253494","#0868ac","#43a2ca","#7bccc4","#bae4bc")) +
  scale_fill_manual(values=c("#253494","#0868ac","#43a2ca","#7bccc4","#bae4bc")) +
  theme_bw()+scale_y_continuous(limits=c(0,1.19), labels=c("0.00","0.25","0.50","0.75","1.00"," "))+
  theme(legend.position = "none", plot.title = element_text(hjust=0.5))+
  geom_signif(comparisons = list(c("3.BIS","4.NAR"),c("1.PVD", "4.NAR"), c("1.PVD","5.NIN")), 
              map_signif_level=c("*"=0.01, "*"=0.05), vjust=0.5, y_position = c(1.05,1.1,1.15))+
  stat_compare_means(label="p.format", label.y = 0.1, label.x="4.NAR", size=4)
  
waterdiv<-metadatags %>% filter(SampleType=="water") %>% 
  ggplot(aes(x=Station,y=Simpsons, fill=Station))+
  geom_jitter(width=0.15, size=3, shape=23)+
  geom_boxplot(alpha=0.8)+
  labs(title="Water samples (n=2)", x=NULL,y=NULL,fill="Site")+
  scale_color_manual(values=c("#253494","#0868ac","#43a2ca","#7bccc4","#bae4bc")) +
  scale_fill_manual(values=c("#253494","#0868ac","#43a2ca","#7bccc4","#bae4bc")) +
  theme_bw()+scale_y_continuous(limits=c(0,1.19), labels=c("0.00","0.25","0.50","0.75","1.00"," "))+
  theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y=element_blank(),
        plot.title = element_text(hjust=0.5))+
  stat_compare_means(label="p.format", label.y = 0.1, label.x="4.NAR", size=4)


plot_grid(gutdiv,waterdiv,rel_widths = c(60,50))

####
### + Alpha diversity stats -----------------------------------------------------
####

dataag_gut<-filter(metadatags, SampleType=="gut")
dataag_water<-filter(metadatags, SampleType=="water")

# within water samples
compare_means(data=dataag_water, Simpsons ~ Station, method="kruskal")
kruskal.test(dataag_water$Simpsons, as.factor(dataag_water$Station))
# From the output of the Kruskal-Wallis test, we know that there is no significant difference between groups. 
pairwise.wilcox.test(dataag_water$Simpsons, dataag_water$Station, p.adjust.method = "BH")

# within gut samples
compare_means(data=dataag_gut, Simpsons ~ Station, method="kruskal")
kruskal.test(dataag_gut$Simpsons, as.factor(dataag_gut$Station))
# From the output of the Kruskal-Wallis test, we know that there is a significant difference between groups, but we don't know which pairs of groups are different.
# It's possible to use the function pairwise.wilcox.test() to calculate pairwise comparisons between group levels with corrections for multiple testing.
pairwise.wilcox.test(dataag_gut$Simpsons, dataag_gut$Station, p.adjust.method = "BH")


# all together now!
compare_means(data=metadatags, Simpsons~Station, group.by="SampleType", method="wilcox", p.adjust.method = "BH")


####
### Figure 3B. Beta diversity: NMDS Bray-Curtis ----------------
####

# ellipse function
theme_set(theme_bw())
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100)
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

## NMDS with only gut samples, examine site effect
# subset data
abund_table_gut<-dataASVtable[c(1:50),]
# generate bray-curtis results
sol<-metaMDS(abund_table_gut,distance = "bray", k = 2, trymax = 50)
sol$stress
stressplot(sol) # stress=0.19, R2=0.96
# get spread of gut points based on Station and generate ellipses
NMDSgut=data.frame(x=sol$point[,1],y=sol$point[,2],Type=as.factor(dataag_gut$SampleType),Station=as.factor(dataag_gut$Station),OysterNum=as.factor(dataag_gut$OysterNum))
plot.new()
ordgut<-ordiellipse(sol, as.factor(dataag_gut$Station), display = "sites", kind ="sd", conf = 0.95, label = T)
dev.off()
# make a dataframe with all the information needed for plotting
df_ellgut <- data.frame()
for(g in levels(NMDSgut$Station)){
  if(g!="" && (g %in% names(ordgut))){
    df_ellgut <- rbind(df_ellgut, cbind(as.data.frame(with(NMDSgut[NMDSgut$Station==g,], veganCovEllipse(ordgut[[g]]$cov,ordgut[[g]]$center,ordgut[[g]]$scale))),Station=g))}}
head(df_ellgut)
# make a dataframe with the ellipse centers
NMDSgut.mean=aggregate(NMDSgut[,1:2],list(group=NMDSgut$Station),mean)
# calculate PERMANOVA significance with adonis2
adonis2gutsite<-adonis2(abund_table_gut~Station, data=dataag_gut, by=NULL, method="bray", k=2)
# plot the NMDS
gutNMDS<-ggplot(data=NMDSgut, aes(x,y))+
  # add the ellipses
  geom_path(data=df_ellgut, aes(x=NMDS1, y=NMDS2,colour=Station, lty=Station), size=1) +
  # add the points for each sample
  geom_point(size=4, alpha=0.9,aes(fill=Station,colour=Station, shape=Station))+
  # add the labels for each site at the center of the ellipses
  geom_text(data=NMDSgut.mean, aes(x=x,y=y,label=group),size=5, color="gray40") +
  # add the global p-value from the adonis2 function
  geom_text(data=adonis2gutsite, aes(x=1,y=-1.2,label=paste("Gut only:\np=",`Pr(>F)`[1],"**")))+
  # edit theme and colors, etc.
  scale_shape_manual(values = c(21,22,23,24,25))+
  scale_fill_manual(values=c("#253494","#0868ac","#43a2ca","#7bccc4","#bae4bc"))+
  scale_colour_manual(values=c("#253494","#0868ac","#43a2ca","#7bccc4","#bae4bc"))+
  scale_linetype_manual(values=c("solid","dotted","twodash","longdash", "solid"), 
                        labels=c("1. Providence River", "2. Greenwich  Bay", "3. Bissel Cove", "4. Narrow River", "5. Ninigret Pond"))+
  theme(legend.text = element_text(size=14, colour="gray20"), legend.position = "right", legend.title = element_blank(), legend.box="horizontal")


## NMDS with all samples, examine sample type effect
abund_table_all<-dataASVtable[c(1:60),]
sol<-metaMDS(abund_table_all,distance = "bray", k = 2, trymax = 50)
sol$stress
stressplot(sol) # stress=0.19, R2=0.96
NMDS=data.frame(x=sol$point[,1],y=sol$point[,2],SampleType=as.factor(metadatags$SampleType), Station=as.factor(metadatags$Station),OysterNum=as.factor(metadatags$OysterNum))
plot.new()
ord<-ordiellipse(sol, as.factor(metadatags$SampleType), display = "sites", kind ="sd", conf = 0.95, label = T)
dev.off()
df_ell <- data.frame()
for(g in levels(NMDS$SampleType)){
  if(g!="" && (g %in% names(ord))){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$SampleType==g,], veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale))),SampleType=g))}}
head(df_ell)
NMDS.mean=aggregate(NMDS[,1:2],list(group=NMDS$SampleType),mean)
head(NMDS.mean)
adonis2type<-adonis2(abund_table_all~SampleType, data=metadatags, by=NULL, method="bray", k=2)
typeNMDS<-ggplot(data=NMDS,aes(x,y))+
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,colour=SampleType, lty=SampleType), size=1) +
  geom_point(size=4, alpha=0.9,aes(shape=Station,colour=SampleType,fill=SampleType))+
  annotate("text",x=NMDS.mean$x,y=NMDS.mean$y,label=NMDS.mean$group,size=6, color="gray10") +
  scale_fill_manual(values=c("orange","darkred"), labels=c("Gut","Water"))+
  scale_colour_manual(values=c("orange","darkred"), labels=c("Gut","Water"))+
  scale_linetype_manual(values=c("solid","twodash"), labels=c("Gut","Water"))+
  scale_shape_manual(values = c(21,22,23,24,25))+
  theme(legend.text = element_text(size=14, colour="gray20"), legend.position = "right",legend.title = element_blank())+
  geom_text(data=adonis2type, aes(x=1,y=-2,label=paste("p=",`Pr(>F)`[1],"**")))




####
###  Figure 3C. Within site beta-diversity -----------------------------------------------------
####

stationids<-dataag_gut %>% select(Station, SampleID)

braycurtisdistances<-
  # calculate the dissimilarity matrix between each sample
  vegdist(abund_table_gut, method="bray", k=2) %>% as.matrix() %>% as.data.frame() %>% 
  # Add in the sites based on the rows (paired sample)
  rownames_to_column(var="PairedID") %>% mutate(PairedStation=dataag_gut$Station) %>% 
  # Make into longform based on the columns (OG SampleID)
  pivot_longer(cols="TNC01":"TNC50", names_to="SampleID") %>% 
  # Add in sites based on the columns (OG SampleID)
  left_join(stationids) %>% 
  # Remove rows where a sample is paired with a sample from another site
  filter(PairedStation==Station) %>% 
  # Remove the diagonal rows, where each sample was compared to itself
  filter(value!=0)

# stats
compare_means(data=braycurtisdistances, value~Station, method="wilcox", p.adjust.method = "BH")
compare_means(data=braycurtisdistances, value~Station, method="kruskal")
# filter out only significant p-values
compare_means(data=braycurtisdistances, value~Station, method="wilcox", p.adjust.method = "BH") %>% 
  filter(p.adj<=0.05)

# plot within site variation
braydistanceplot<-ggplot(braycurtisdistances, aes(x=Station, y=value, fill=Station))+
  geom_jitter(width=0.15, size=2, shape=23, alpha=0.7)+geom_boxplot(alpha=0.8)+
  labs(x=NULL,y="within site dissimilarity index",fill="Site")+
  scale_color_manual(values=c("#253494","#0868ac","#43a2ca","#7bccc4","#bae4bc")) +
  scale_fill_manual(values=c("#253494","#0868ac","#43a2ca","#7bccc4","#bae4bc")) +
  theme_bw()+scale_y_continuous(limits=c(0,1.19), labels=c("0.00","0.25","0.50","0.75","1.00"," "))+
  theme(legend.position = "none", strip.background = element_rect("grey90"))+
  # add statistics
  geom_signif(comparisons = list(c("4.NAR","5.NIN"),c("1.PVD", "2.GB")), 
              map_signif_level=TRUE, vjust=0.5, y_position = c(1.05,1.05))+
  geom_signif(annotations = c("***"), y_position = 1.18, xmin=c(2), vjust=0.5, xmax=c(4.5),tip_length = c(0.08,0.04))+
  geom_signif(annotations = c("", ""), y_position = 1.15, xmin=c(1,4), xmax=c(3,5))+
  # add global p-value
  stat_compare_means(label="p.format", label.y = 0.1, label.x="4.NAR", size=4)
  


####
### + Build Figure 3 with cowplots ------------------------------------------
####

legends<-plot_grid(get_legend(gutNMDS+scale_linetype_discrete(guide = FALSE)),
                   get_legend(typeNMDS+scale_shape_discrete(guide = FALSE)),
                   ncol = 2,rel_heights = c(20,50))

divplot<-plot_grid(gutdiv,waterdiv)
nmds<-plot_grid(gutNMDS+theme(legend.position = "none"),
                typeNMDS+theme(legend.position = "none"),
                align="hv", axis="t", nrow=1)
braydistance_legends<-plot_grid(braydistanceplot, legends)

plot_grid(divplot+theme(plot.margin = margin(0, 0, 0, 7)),
          nmds+draw_label("Bray-Curtis beta diversity (k=2)", size=11, x=0, y=0.5, angle=90), 
          braydistance_legends+draw_label("Bray-Curtis beta diversity", size=11, x=0, y=0.55, angle=90),
          nrow=3, align="v", axis="tl",rel_heights = c(50,60,50),
          labels = "AUTO")
#ggsave("Figure3.png", width = 8, height = 9, dpi=400)
#ggsave("Figure3.pdf", width = 8, height = 9)
#ggsave("Figure3.svg", width = 8, height = 9)



####
### Table S2. Beta diversity stats -----------------------------------------------------
####

beta<-NULL

# by station
beta$overallsites<-adonis2(abund_table_gut~Station, data=dataag_gut, by=NULL, method="bray", k=2)
# pvd vs gb
beta$pvdgb<-adonis2(abund_table_gut[c(1:20),]~Station, data=dataag_gut[c(1:20),], by=NULL, method="bray", k=2)
# pvd vs bis
beta$pvdbis<-adonis2(abund_table_gut[c(1:10,21:30),]~Station, data=dataag_gut[c(1:10,21:30),], by=NULL, method="bray", k=2)
# pvd vs nar
beta$pvdnar<-adonis2(abund_table_gut[c(1:10,31:40),]~Station, data=dataag_gut[c(1:10,31:40),], by=NULL, method="bray", k=2)
# pvd vs nin
beta$pvdnin<-adonis2(abund_table_gut[c(1:10,41:50),]~Station, data=dataag_gut[c(1:10,41:50),], by=NULL, method="bray", k=2)
# gb vs bis
beta$gbbis<-adonis2(abund_table_gut[c(11:30),]~Station, data=dataag_gut[c(11:30),], by=NULL, method="bray", k=2)
# gb vs nar
beta$gbnar<-adonis2(abund_table_gut[c(11:20,31:40),]~Station, data=dataag_gut[c(11:20,31:40),], by=NULL, method="bray", k=2)
# gb vs nin
beta$gbnin<-adonis2(abund_table_gut[c(11:20,41:50),]~Station, data=dataag_gut[c(11:20,41:50),], by=NULL, method="bray", k=2)
# bis vs nar
beta$bisnar<-adonis2(abund_table_gut[c(21:40),]~Station, data=dataag_gut[c(21:40),], by=NULL, method="bray", k=2)
# bis vs nin
beta$bisnin<-adonis2(abund_table_gut[c(21:30,41:50),]~Station, data=dataag_gut[c(21:30, 41:50),], by=NULL, method="bray", k=2)
# nar vs nin
beta$narnin<-adonis2(abund_table_gut[c(31:50),]~Station, data=dataag_gut[c(31:50),], by=NULL, method="bray", k=2)

# by sample type
beta$sampletypeonly <- adonis2(abund_table_all~SampleType, data=metadatags, by=NULL, method="bray", k=2)
beta$overallsampletype <- adonis2(abund_table_all~SampleType*Station, data=metadatags, by="terms", method="bray", k=2)

# export adonis2 results to use as table
betanames <- beta %>% enframe() %>% unnest(cols=c(value)) %>% pull(name)
beta %>% bind_rows() %>% rownames_to_column() %>% mutate(comparison=betanames) %>% write.csv("betadiversity.csv")



####
### Figure S6 and Table S3. Core Microbiome ----------------------------------------------
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
  left_join(taxakey) %>% # write.csv("coremicrobiome16s.csv")
  left_join(dataASVwide) %>% 
  left_join(metadatags) %>% 
  filter(SampleType=="gut") %>% 
  # make a binary scale for when the ASV is present in each sample
  mutate(presence=case_when(percent==0 ~ "0", TRUE ~ "1")) %>%
  mutate(taxonlabel=paste(Phylum,Class,Order,Family,Genus,Species,sep="; "))

core16splot<-ggplot(coretaxa, 
       aes(x=SampleID, y=reorder(taxonlabel,n), fill=presence))+
  geom_tile(color="white")+
  facet_grid(.~Station,scales="free", space="free")+
  scale_fill_manual(values=c("grey80","aquamarine4"))+
  scale_y_discrete(labels=wrap_format(50))+
  theme_minimal()+
  theme(legend.position = "none", axis.text.x=element_blank(),
        axis.text.y = element_text(lineheight=0.7))+
  labs(x=NULL, y=NULL)


