# Stevick et al 2020 Oyster Gut Microbiome Function in an Estuary
# Compare taxonomy generated from RefSeq/metatranscriptomes and 16S data
# Figures 4, S2, S3, S5, and S6. Table S3
# updated 9/17/2020 for reresubmission

# Metaranscriptomic taxonomy data
# 16S Amplicon data at order level

library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggpubr)
library(readxl)
library(vegan)
library(cowplot)

data<-read_xlsx("Taxonomy/SAMSA_metatranscriptomes_RefSeqtaxa_edit.xlsx", sheet=2)
taxakey<- read_xlsx("Taxonomy/SAMSA_metatranscriptomes_RefSeqtaxa_edit.xlsx", sheet=3)
data$Sample<-as.character(data$Sample)
datatax<-left_join(data, unique(taxakey), by="Taxa")

dataphylasum<-datatax %>% group_by(Phylum, Sample) %>%
  dplyr::summarise(sumPhylum=sum(Percent)) %>%
  group_by(Phylum) %>% dplyr::summarise(mean=mean(sumPhylum), stdev=sd(sumPhylum))

dataordsum<-datatax %>% group_by(Order, Sample) %>%
  dplyr::summarise(sumOrder=sum(Percent)) %>%
  group_by(Order) %>% dplyr::summarise(mean=mean(sumOrder), stdev=sd(sumOrder)) %>%
  arrange(desc(mean))

dataord <- datatax %>% group_by(Order, Sample) %>%
  dplyr::summarise(sumOrder=sum(Percent))

# clean up metarans data for later venn diagrams
condenseddataord<-datatax %>% select(Sample, Order, Percent, SampleType) %>%
  group_by(Sample, Order, SampleType) %>%
  dplyr::summarise(PercentSum=sum(Percent)) %>% ungroup()
dataordtab<-spread(condenseddataord, Order, PercentSum)
rownames(dataordtab) <- dataordtab$Sample

# make a ginormous palette
palette<-c(brewer.pal(12, "Set3"),brewer.pal(8, "Set2"),brewer.pal(12, "Set3"),
           brewer.pal(12, "Set3"),brewer.pal(8, "Set2"),brewer.pal(12, "Set3"),
           brewer.pal(12, "Set3"),brewer.pal(8, "Set2"),brewer.pal(12, "Set3"),
           brewer.pal(12, "Set3"),brewer.pal(8, "Set2"),brewer.pal(12, "Set3"),
           brewer.pal(12, "Set3"),brewer.pal(8, "Set2"),brewer.pal(12, "Set3"),
           brewer.pal(12, "Set3"),brewer.pal(8, "Set2"),brewer.pal(12, "Set3"),
           brewer.pal(12, "Set3"),brewer.pal(8, "Set2"),brewer.pal(12, "Set3"),
           brewer.pal(12, "Set3"),brewer.pal(8, "Set2"),brewer.pal(12, "Set3"),
           brewer.pal(12, "Set3"),brewer.pal(8, "Set2"),brewer.pal(12, "Set3"))

#metatranscriptomic taxa by species
ggplot(data,aes(x=as.factor(Sample),Percent,fill=Taxa))+
  geom_col(position="fill")+
  theme(legend.text = element_text(size=6, colour="gray20"),
        legend.position = "none",axis.ticks.length = unit(0.4, 'lines'))+
  facet_grid(.~Station, scales="free",space="free")+
  scale_fill_manual(values=palette)+
  labs(y="Percent \nabundance",x=NULL,fill=NULL)+
  scale_y_continuous(labels = scales::percent, expand=c(0,0))

#metatranscriptomic taxa by order
ggplot(datatax,aes(x=as.factor(Sample),Percent,fill=Order_Name))+
  geom_col(position="fill")+
  theme(legend.text = element_text(size=6, colour="gray20"),
        legend.position = "bottom",axis.ticks.length = unit(0.4, 'lines'))+
  facet_grid(.~Station, scales="free",space="free")+
  scale_fill_manual(values=palette)+
  labs(y="Percent \nabundance",x=NULL,fill=NULL)+
  scale_y_continuous(labels = scales::percent)

#metatranscriptomic taxa by phylum
ggplot(datatax,aes(x=as.factor(Sample),Percent,fill=Phylum))+
  geom_col(position="fill")+
  theme(legend.text = element_text(size=6, colour="gray20"),
        legend.position = "bottom",axis.ticks.length = unit(0.4, 'lines'))+
  facet_grid(.~Station, scales="free",space="free")+
  scale_fill_manual(values=palette)+
  labs(y="Percent \nabundance",x=NULL,fill=NULL)+
  scale_y_continuous(labels = scales::percent)



# Add in  16S data
data16t<-read_xlsx("Taxonomy/16SallData_SILVAtaxa.xlsx", sheet="Level4-normtrim")
data16tg<-tidyr::gather(data16t, Order, Percent, "Acidobacteria;Acidobacteriia;Solibacterales":"Unknown")
# add in the extra taxa information from the key
taxakey$Taxa<-NULL
data16tgftax<-inner_join(data16tg, unique(taxakey), by="Order")

# merge metatranscriptomes and 16S data
alldatatax<-full_join(data16tgftax, datatax, by=c("Method","Sample","Station","Order","Order_Name","Class","Phylum","Percent","SampleType"))


#plot it all together by Order
ggplot(alldatatax,aes(x=as.factor(Sample),Percent,fill=Order))+
  geom_col(position="fill")+
  theme(legend.text = element_text(size=6, colour="gray20"),
        legend.position = "none",axis.ticks.length = unit(0.4, 'lines'))+
  facet_grid(.~SampleType+Station+Method, scales="free",space="free")+
  scale_fill_manual(values=palette)+
  labs(y="Percent \nabundance",x=NULL,fill=NULL)+
  scale_y_continuous(labels = scales::percent)

# by class
ggplot(alldatatax,aes(x=as.factor(Sample),Percent,fill=Class))+
  geom_col(position="fill")+
  theme(legend.text = element_text(size=6, colour="gray20"),
        legend.position = "none",axis.ticks.length = unit(0.4, 'lines'),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  facet_grid(.~SampleType+Method+Station, scales="free",space="free")+
  scale_fill_manual(values=palette)+
  labs(y="Percent \nabundance",x=NULL,fill=NULL)+
  scale_y_continuous(labels = scales::percent, expand=c(0,0))

# by phylum
ggplot(alldatatax,aes(x=as.factor(Sample),Percent,fill=Phylum))+
  geom_col(position="fill")+
  theme(legend.text = element_text(size=6, colour="gray20"),
        legend.position = "bottom",axis.ticks.length = unit(0.4, 'lines'),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  facet_grid(.~SampleType+Method+Station, scales="free",space="free")+
  scale_fill_manual(values=palette)+
  labs(y="Percent \nabundance",x=NULL,fill=NULL)+
  scale_y_continuous(labels = scales::percent, expand=c(0,0))



####
### Figure 4B. Heatmap of most abundant taxa by site and sequencing type ----------------------------------------------
####

# squareroot transformation for easier visualization
library(plyr)
alldatatax_norm <- plyr::ddply(alldatatax, .(Sample), transform, rescale=sqrt(Percent))

# remove lowest abundance taxa
toptax<-alldatatax_norm %>% group_by(Order) %>%
  dplyr::summarise(sums=sum(Percent)) %>% # calculate sums
  top_n(30, wt=sums) %>% arrange(desc(sums)) %>% pull(Order) # extract top 30 Orders

topdatatax_norm <- filter(alldatatax_norm, Order %in% toptax)

heatmap<-ggplot(topdatatax_norm, aes(Sample, Order)) +
  geom_tile(aes(fill = rescale),colour = "white") + ggpubr::theme_transparent() +
  facet_grid(factor(Phylum, levels=c("Actinobacteria","Bacteroidetes","Cyanobacteria","Firmicutes","Fusobacteria",
                                     "Proteobacteria","Tenericutes","Verrucomicrobia","Unknown"))~
               factor(SampleType, levels = c("water","gut"))+Method+Station,
             space="free", scales="free", switch="x")+
  scale_fill_gradientn(na.value = "salmon",labels = scales::percent,limits=c(0,1),
                       colours=c("white","#fecc5c","#fd8d3c","#f03b20","#bd0026","darkred"))+
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  theme(legend.position = c(-0.24,-0.06),legend.direction = "horizontal",
        axis.ticks = element_blank(),
        axis.text.y = element_text(size=10, colour="grey40"),
        axis.text.x = element_blank(),
        legend.text = element_text(size=10, colour="grey40"),
        legend.key.size = unit(1.5, 'lines'),legend.title = element_text(size=11, colour="grey40"),
        panel.background = element_rect(fill="white", colour="white"),
        strip.background = element_rect(fill="white", colour="white"),
        strip.text.y = element_blank())+
  labs(y="Most abundant taxa (Order level)", x=NULL, fill="Relative Percent \nAbundance of Order")
# 1030x550



## VENN DIAGRAMS --------------------

library(VennDiagram)
library(gplots)
library(UpSetR)
library(reshape2)
library(ComplexHeatmap)
theme_set(theme_bw())

# subset data
gutdata16t<-data16t[data16t$SampleType=="gut", apply(data16t[data16t$SampleType=="gut",], MARGIN=2, function(x) any(x >0))]
gut16S=colnames(data16t[data16t$SampleType=="gut", apply(data16t[data16t$SampleType=="gut",], MARGIN=2, function(x) any(x >0))])
gut16S<-gut16S[1:296]
gut16S1=colnames(gutdata16t[gutdata16t$Station=="1.PVD", apply(gutdata16t[gutdata16t$Station=="1.PVD",], MARGIN=2, function(x) any(x >0))])
gut16S1<-gut16S1[1:179]
gut16S2=colnames(gutdata16t[gutdata16t$Station=="2.GB", apply(gutdata16t[gutdata16t$Station=="2.GB",], MARGIN=2, function(x) any(x >0))])
gut16S2<-gut16S2[1:167]
gut16S3=colnames(gutdata16t[gutdata16t$Station=="3.BIS", apply(gutdata16t[gutdata16t$Station=="3.BIS",], MARGIN=2, function(x) any(x >0))])
gut16S3<-gut16S3[1:198]
gut16S4=colnames(gutdata16t[gutdata16t$Station=="4.NAR", apply(gutdata16t[gutdata16t$Station=="4.NAR",], MARGIN=2, function(x) any(x >0))])
gut16S4<-gut16S4[1:198]
gut16S5=colnames(gutdata16t[gutdata16t$Station=="5.NIN", apply(gutdata16t[gutdata16t$Station=="5.NIN",], MARGIN=2, function(x) any(x >0))])
gut16S5<-gut16S5[1:172]
water16S=colnames(data16t[data16t$SampleType=="water", apply(data16t[data16t$SampleType=="water",], MARGIN=2, function(x) any(x >0))])
water16S<-water16S[1:143]
guttrans=colnames(dataordtab[dataordtab$SampleType=="gut", apply(dataordtab[dataordtab$SampleType=="gut",], MARGIN=2, function(x) any(x >0))])
guttrans<-guttrans[2:68]

# simple venn diagrams
venn(list(gut16S=gut16S, water16S=water16S))
sitestogether<-venn(list("1.PVD Gut 16S"=gut16S1,
          "2.GB Gut 16S"=gut16S2,
          "3.BIS Gut 16S"=gut16S3,
          "4.NAR Gut 16S"=gut16S4,
          "5.NIN Gut 16S"=gut16S5))
vennlist<- venn(list("gut 16S"=gut16S, "water 16S"=water16S, "gut metatranscriptome"=guttrans))
summaryvenn<-attr(vennlist, "intersections")
summaryvennsites<-attr(sitestogether, "intersections")

set1 <- guttrans
set2 <- gut16S
set3 <- water16S
read_sets = list(gutRNA = set1,gut16S = set2, water16S = set3)
m = make_comb_mat(read_sets)

upset(fromList(read_sets),
      number.angles = 0, sets=c("gutRNA", "gut16S", "water16S"), point.size = 5, line.size = 1.3,
      mainbar.y.label = NULL, sets.x.label = NULL,
      text.scale = c(1.5, 1.5, 1.25, 1.25, 1, 1.8), mb.ratio = c(0.6, 0.4),
      keep.order = TRUE, order.by = "degree")

UpSet(t(m), set_order = order(c("water16S","gutRNA","gut16S")),
      top_annotation = upset_top_annotation(t(m),bar_width = 0.9))



####
### Figure 4A. Upset plot of orders shared between sequencing types ----------------------------------------------
####


top<-UpSet(m, set_order = order(c("water16S","gutRNA","gut16S")),
      pt_size = unit(.35, "cm"),lwd=3,
      left_annotation = rowAnnotation(" " = anno_barplot(set_size(m), bar_width=0.7,
                            axis_param = list(direction = "reverse", side = "top",labels_rot = 0),
                            border = FALSE, annotation_name_side = "top",
                            gp = gpar(fill = "black"),
                            width = unit(4, "cm"))),
      right_annotation = NULL,row_names_side = "left",
      top_annotation = upset_top_annotation(m,bar_width = 0.9))
#1000x155

read_sets = list("1.PVD Gut 16S"=gut16S1,
                 "2.GB Gut 16S"=gut16S2,
                 "3.BIS Gut 16S"=gut16S3,
                 "4.NAR Gut 16S"=gut16S4,
                 "5.NIN Gut 16S"=gut16S5,
                 "All water 16S"=water16S)
mgutw = make_comb_mat(read_sets)


####
### Figure S3. Upset plot of orders shared in 16S samples  ----------------------------------------------
####

#900x400

UpSet(mgutw,comb_col = c("#253494","#0868ac","#43a2ca","#7bccc4","#bae4bc","grey40",
                         "black","black","black","grey40","black","black",
                         "black","grey40","black","black","grey40","black",
                         "grey40","grey40","black","black","black","black",
                         "grey40","black","grey40","black","black","grey40",
                         "black","grey40","black","grey40","grey40","black",
                         "black","grey40","black","grey40","black","grey40",
                         "grey40","black","grey40", "black", "grey40", "grey40",
                         "grey40", "grey40", "grey40", "grey40"),
      set_order = order(c("2.GB Gut 16S","3.BIS Gut 16S",
                          "4.NAR Gut 16S", "5.NIN Gut 16S","All water 16S","1.PVD Gut 16S")),
      pt_size = unit(.5, "cm"),lwd=3,
      right_annotation = rowAnnotation(" " = anno_barplot(set_size(mgutw), bar_width=0.7,
                                                         axis_param = list(side = "top",labels_rot = 0),
                                                         border = FALSE, annotation_name_side = "top",
                                                         gp = gpar(fill = c("#253494","#0868ac","#43a2ca","#7bccc4","#bae4bc","grey40")),
                                                         width = unit(4, "cm"))),
      row_names_side = "left",
      top_annotation = upset_top_annotation(mgutw,bar_width = 0.9, height = unit(6, "cm")))



### + NMDS of Metatranscriptome taxonomy -----------------

Site<-as.factor(c("1.PVD","1.PVD","1.PVD","1.PVD","1.PVD",
          "2.GB","2.GB","2.GB","2.GB","2.GB",
          "3.BIS","3.BIS","3.BIS","3.BIS","3.BIS",
          "4.NAR","4.NAR","4.NAR","4.NAR","4.NAR",
          "5.NIN","5.NIN","5.NIN","5.NIN","5.NIN"))

library(vegan)
theme_set(theme_bw())
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100)
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

####
### Figure S5A. Species level NMDS ----------------------------------------------
####

abund_tablesp<-data %>%
  select(Sample, Taxa, Percent) %>%
  spread(Taxa, Percent)
abund_tablesp[is.na(abund_tablesp)]<-0
rownames(abund_tablesp)<-abund_tablesp$Sample
abund_tablesp$Sample<-NULL

sol<-metaMDS(abund_tablesp,distance = "bray", k = 2, trymax = 50)
sol$stress
stressplot(sol)
# stress=0.0825, R2=0.993

NMDS=data.frame(x=sol$point[,1],y=sol$point[,2], Site)
plot.new()
ord<-ordiellipse(sol, as.factor(Site), display = "sites", kind ="sd", conf = 0.95, label = T)
dev.off()
df_ell <- data.frame()
for(g in levels(NMDS$Site)){
  if(g!="" && (g %in% names(ord))){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$Site==g,],
                    veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale))),Site=g))}}
head(df_ell)
NMDS.mean=aggregate(NMDS[,1:2],list(group=NMDS$Site),mean)

spadonis2<-adonis2(abund_tablesp~Site, data=NMDS, by=NULL,method="bray", k=2)

metatranssp<-ggplot(data=NMDS,aes(x,y)) +
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2, lty=Site,colour=Site), size=1) +
  geom_point(size=4, alpha=0.9,aes(shape=Site,colour=Site,fill=Site))+scale_shape_manual(values = c(21,22,23,24,25))+
  annotate("text",x=NMDS.mean$x,y=NMDS.mean$y,label=NMDS.mean$group,size=5, color="gray40") +
  scale_fill_manual(values=c("#253494","#0868ac","#43a2ca","#7bccc4","#bae4bc"))+
  scale_colour_manual(values=c("#253494","#0868ac","#43a2ca","#7bccc4","#bae4bc"))+
  scale_linetype_manual(values=c("solid","dotted","twodash","longdash", "solid"), labels=c("1. Providence River", "2. Greenwich  Bay", "3. Bissel Cove", "4. Narrow River", "5. Ninigret Pond"))+
  theme(legend.text = element_text(size=14, colour="gray20"),
        legend.position = "none",
        legend.title = element_blank(),legend.box="horizontal") +
  ggtitle("Species level annotation")+
  # add the global p-value from the adonis2 function
  geom_text(data=spadonis2, aes(x=1,y=-0.5,label=paste("p=",`Pr(>F)`[1],"*")))
  


####
### Figure S5B. Order level NMDS ----------------------------------------------
####

abund_tableord<-dataord %>%
  select(Sample, Order, sumOrder) %>%
  spread(Order, sumOrder)
abund_tableord[is.na(abund_tableord)]<-0
rownames(abund_tableord)<-abund_tableord$Sample
abund_tableord$Sample<-NULL

sol<-metaMDS(abund_tableord,distance = "bray", k = 2, trymax = 50)
sol$stress
stressplot(sol)
# stress=0.075, R2=0.994

NMDS=data.frame(x=sol$point[,1],y=sol$point[,2], Site)
plot.new()
ord<-ordiellipse(sol, as.factor(Site), display = "sites", kind ="sd", conf = 0.95, label = T)
dev.off()
df_ell <- data.frame()
for(g in levels(NMDS$Site)){
  if(g!="" && (g %in% names(ord))){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$Site==g,],
                                                     veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale))),Site=g))}}
head(df_ell)
NMDS.mean=aggregate(NMDS[,1:2],list(group=NMDS$Site),mean)

ordadonis2<-adonis2(abund_tableord~Site, data=NMDS, by=NULL,method="bray", k=2)

metatransord<-ggplot(data=NMDS,aes(x,y))+
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2, lty=Site, colour=Site), size=1) +
  geom_point(size=4, alpha=0.9,aes(shape=Site, colour=Site,fill=Site))+scale_shape_manual(values = c(21,22,23,24,25))+
  annotate("text",x=NMDS.mean$x,y=NMDS.mean$y,label=NMDS.mean$group,size=5, color="gray40") +
  scale_fill_manual(values=c("#253494","#0868ac","#43a2ca","#7bccc4","#bae4bc"))+
  scale_colour_manual(values=c("#253494","#0868ac","#43a2ca","#7bccc4","#bae4bc"))+
  scale_linetype_manual(values=c("solid","dotted","twodash","longdash", "solid"), labels=c("1. Providence River", "2. Greenwich  Bay", "3. Bissel Cove", "4. Narrow River", "5. Ninigret Pond"))+
  theme(legend.text = element_text(size=14, colour="gray20"),
        legend.position = "right",
        legend.title = element_blank(),legend.box="vertical")+
  ggtitle("Order level annotation")+
  # add the global p-value from the adonis2 function
  geom_text(data=ordadonis2, aes(x=0.6,y=-0.3,label=paste("p=",`Pr(>F)`[1],"**")))


####
### + Combine Figure S5 ----------------------------------------------
####

#1000x800
plot_grid(metatranssp,get_legend(metatransord), metatransord+theme(legend.position="none"),
          ncol=2, rel_widths = c(50,10))


####
### Figure S2. Metatranscriptome Rarefaction curve ----------------------------------------------
####

# Load in the raw annotation data
data_table<-read_xlsx("Taxonomy/SAMSA_metatranscriptomes_RefSeqtaxa_edit.xlsx", sheet="refseqspecies")

# Make into table
datamatt<-data_table %>%
  select(Taxa, contains("reads")) %>%
  tibble::column_to_rownames("Taxa") %>%
  data.table::transpose()

raremax <- min(rowSums(datamatt))
S <- specnumber(datamatt) # observed number of annotated species
Srare <- rarefy(datamatt, raremax)

#define colors based on site
colors<-c("#253494","#253494","#253494","#253494","#253494",
          "#0868ac","#0868ac","#0868ac","#0868ac","#0868ac",
          "#43a2ca","#43a2ca","#43a2ca","#43a2ca","#43a2ca",
          "#7bccc4","#7bccc4","#7bccc4","#7bccc4","#7bccc4",
          "#bae4bc","#bae4bc","#bae4bc","#bae4bc","#bae4bc")
#put plotting parameters into pars
pars <- expand.grid(col = colors, stringsAsFactors = FALSE)

#plot number of species annotated vs. rarefied
#550x400
plot(S, Srare, pch=21, cex=2, col="black", bg=colors,
     xlab = "Observed Number of Annotated Species", ylab = "Rarefied Number of Annotated Species")
abline(0, 1, lty="dashed")
legend(700,250, c("1.PVD gut","2.GB gut","3.BIS gut","4.NAR gut","5.NIN gut", "All water"),
       cex=0.8,pch=19,
       text.width = strwidth("60000000"),xjust = 1, yjust = 1,
       col=c("#253494","#0868ac","#43a2ca","#7bccc4","#bae4bc","darkred"))


#generate rarefaction curves
#NOTE: THIS STEP TAKES A WHILE TO RUN.
# Step size has been increased so it doesn't take a day or so
outw <- with(pars, rarecurve(datamatt, step = 800,
                             sample = raremax, col = colors, label = FALSE))
#determine limits of the graph
Smax <- sapply(outw, max)
#plot the rarefaction curves - Figure S2A
#550x400
par(mfrow=c(3,1))
plot(c(1, 1103294), c(1, max(Smax)), xlab = "Number of Reads",
     ylab = "Observed Number of Annotated Species", type = "n")
abline(v = raremax, lty="dashed")
for (i in seq_along(outw)) {
  N <- attr(outw[[i]], "Subsample")
  with(pars, lines(N, outw[[i]], col = col[i]))
}

# check slopes of curves at the end, make sure they're ~0
end<-raremax-1
slopes<-rareslope(datamatt, end)
# Figure S2B
plot(slopes, type="p", pch=23, bg = colors,
     cex=2, xlab = "Metatranscriptome Samples",
     ylab= "Slopes at Raremax", ylim=c(0,0.05))
# coverage proxy:
coverage<-100-100*rareslope(datamatt, end)
# Figure S2C
plot(coverage, type="p", pch=23, bg = colors,
     cex=2, xlab = "Metatranscriptome Samples",
     ylab="Estimated Coverage (%)", ylim=c(95,100))

mean(coverage)
sd(coverage)




####
### Figure S6 and Table S3. Core Microbiome ----------------------------------------------
####


# with order level
dataord %>%
  # remove when the order is absent
  filter(sumOrder!=0) %>%
  # group by order and the count the number of samples it occurs in
  group_by(Order) %>% count() %>%
  # keep orders that occur in at least 80% of samples
  arrange(desc(n)) %>% filter(n>=20) %>% # write.csv("coremicrobiome16s.csv")
  # add back in the other metadata
  left_join(dataord) %>%
  # make a binary scale for when the order is present in each sample
  mutate(presence=case_when(sumOrder==0 ~ "0", TRUE ~ "1")) %>%
  # now plot it
  ggplot(aes(x=Sample, y=reorder(Order,n), fill=presence))+
  geom_tile(color="white")+
  scale_fill_manual(values="aquamarine4", na.value="grey80")+
  theme_minimal()+
  theme(legend.position = "none")+
  labs(x=NULL, y=NULL)

# make a dataframe of all possible combos
all<-data %>%  filter(Percent!=0) %>% group_by(Taxa) %>% count() %>%
  arrange(desc(n)) %>% filter(n>=20) %>% left_join(data) %>% expand(Taxa, 26:50) %>%
  mutate(Sample=as.character(`26:50`)) %>%
  select(-`26:50`)

stations<-datatax %>% select(Sample, Station) %>%
  group_by(Sample, Station) %>% summarise() %>%
  full_join(all)

# with species level
coretaxa<-data %>%
  # remove when the species is absent
  filter(Percent!=0) %>%
  # group by Taxon and the count the number of samples it occurs in
  group_by(Taxa) %>% count() %>%
  # keep Taxa that occur in at least 80% of samples
  arrange(desc(n)) %>% filter(n>=20) %>%
  right_join(stations) %>%
  # add back in the other metadata
  left_join(taxakey) %>% # write.csv("coremicrobiomemetatrans.csv")
  left_join(data) %>%
  # make a binary scale for when the order is present in each sample
  mutate(presence=case_when(is.na(Percent) ~ "0", TRUE ~ "1")) %>%
  mutate(taxonlabel=paste(Phylum,Taxa,sep="; "))

coremetatransplot <-
  ggplot(coretaxa, aes(x=Sample, y=reorder(taxonlabel,n), fill=presence))+
  geom_tile(color="white")+
  facet_grid(.~Station,scales="free", space="free")+
  scale_fill_manual(values=c("grey80","aquamarine4"))+
  theme_minimal()+labs(x=NULL, y=NULL)+
  theme(legend.position = "none", axis.text.x=element_blank(),
        strip.text=element_blank())

plot_grid(core16splot, coremetatransplot,
          nrow=2, align="v", axis="l", rel_heights = c(15,30),labels="AUTO")+
  draw_label("DNA: 16S rRNA gene amplicons", size=11, x=1, y=0.82, angle=-90)+
  draw_label("RNA: Metatranscriptomes", size=11, x=1, y=0.3, angle=-90)+
  theme(plot.margin = margin(0, 10, 0, 0))

#ggsave("FigureS6.png", width = 8, height = 12, dpi=400)
#ggsave("FigureS6.pdf", width = 8, height = 12)
#ggsave("FigureS6.svg", width = 8, height = 12)
