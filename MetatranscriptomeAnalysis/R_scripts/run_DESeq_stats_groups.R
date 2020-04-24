# Created 6/22/16, updated 6/16/2017
# Run with --help flag for help.

suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--input"), type="character", default="./",
              help="Input directory", metavar="character"),
  make_option(c("-O", "--out"), type="character", default="DESeq_results.tab", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-R", "--raw_counts"), type="character", default=NULL,
              help="raw (total) read counts for this starting file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print("USAGE: $ run_DESeq_stats.R -I working_directory/ -O save.filename")

# check for necessary specs
if (is.null(opt$input)) {
  print ("WARNING: No working input directory specified with '-I' flag.")
  stop()
} else {  cat ("Working directory is ", opt$input, "\n")
  wd_location <- opt$input  
  setwd(wd_location)  }

if (is.null(opt$out)) {
  print ("WARNING: No save name for DESeq results specified; defaulting to 'DESeq_results.tab'.") 
  save_filename <- opt$out
} else { save_filename <- opt$out }

if (is.null(opt$raw_counts)) {
  print ("WARNING: no raw counts file specified, skipping this info for DESeq analysis.")
} else {
  counts_file <- opt$raw_counts
}

# import other necessary packages
suppressPackageStartupMessages({
  library(DESeq2)
})

# GET FILE NAMES
pvd_files <- list.files(
  pattern = "1PVD_*", full.names = T, recursive = FALSE)
pvd_names = ""
for (name in pvd_files) {
  pvd_names <- c(pvd_names, unlist(strsplit(name, split='_', fixed=TRUE))[2])}
pvd_names <- pvd_names[-1]
pvd_names_trimmed = ""
for (name in pvd_names) {
  pvd_names_trimmed <- c(pvd_names_trimmed, unlist(strsplit(name, split='.', fixed=TRUE))[1])}
pvd_names_trimmed <- pvd_names_trimmed[-1]

gb_files <- list.files(
  pattern = "2GB_*", full.names = T, recursive = FALSE)
gb_names = ""
for (name in gb_files) {
  gb_names <- c(gb_names, unlist(strsplit(name, split='_', fixed=TRUE))[2])}
gb_names <- gb_names[-1]
gb_names_trimmed = ""
for (name in gb_names) {
  gb_names_trimmed <- c(gb_names_trimmed, unlist(strsplit(name, split='.', fixed=TRUE))[1])}
gb_names_trimmed <- gb_names_trimmed[-1]

bis_files <- list.files(
  pattern = "3BIS_*", full.names = T, recursive = FALSE)
bis_names = ""
for (name in bis_files) {
  bis_names <- c(bis_names, unlist(strsplit(name, split='_', fixed=TRUE))[2])}
bis_names <- bis_names[-1]
bis_names_trimmed = ""
for (name in bis_names) {
  bis_names_trimmed <- c(bis_names_trimmed, unlist(strsplit(name, split='.', fixed=TRUE))[1])}
bis_names_trimmed <- bis_names_trimmed[-1]

nar_files <- list.files(
  pattern = "4NAR_*", full.names = T, recursive = FALSE)
nar_names = ""
for (name in nar_files) {
  nar_names <- c(nar_names, unlist(strsplit(name, split='_', fixed=TRUE))[2])}
nar_names <- nar_names[-1]
nar_names_trimmed = ""
for (name in nar_names) {
  nar_names_trimmed <- c(nar_names_trimmed, unlist(strsplit(name, split='.', fixed=TRUE))[1])}
nar_names_trimmed <- nar_names_trimmed[-1]

nin_files <- list.files(
  pattern = "5NIN_*", full.names = T, recursive = FALSE)
nin_names = ""
for (name in nin_files) {
  nin_names <- c(nin_names, unlist(strsplit(name, split='_', fixed=TRUE))[2])}
nin_names <- nin_names[-1]
nin_names_trimmed = ""
for (name in nin_names) {
  nin_names_trimmed <- c(nin_names_trimmed, unlist(strsplit(name, split='.', fixed=TRUE))[1])}
nin_names_trimmed <- nin_names_trimmed[-1]


# READ IN FILES
# loading the pvd table
y <- 0
for (x in pvd_files) {
  y <- y + 1
  if (y == 1) {
    pvd_table <- read.table(file = x, header = F, quote = "", sep = "\t")
    colnames(pvd_table) = c("DELETE", x, "V3")
    pvd_table <- pvd_table[,c(2,3)]      }
  if (y > 1) {
    temp_table <- read.table(file = x, header = F, quote = "", sep = "\t")
    colnames(temp_table) = c("DELETE", x, "V3")
    temp_table <- temp_table[,c(2,3)]
    pvd_table <- merge(pvd_table, temp_table, by = "V3", all = T)  }
}
pvd_table[is.na(pvd_table)] <- 0
rownames(pvd_table) = pvd_table$V3
pvd_table_trimmed <- pvd_table[,-1]

# loading the gb table
y <- 0
for (x in gb_files) {
  y <- y + 1
  if (y == 1) {
    gb_table <- read.table(file = x, header=F, quote = "", sep = "\t")
    colnames(gb_table) = c("DELETE", x, "V3")
    gb_table <- gb_table[,c(2,3)]  }
  if (y > 1) {
    temp_table <- read.table(file = x, header = F, quote = "", sep = "\t")
    colnames(temp_table) = c("DELETE", x, "V3")
    temp_table <- temp_table[,c(2,3)]
    gb_table <- merge(gb_table, temp_table, by = "V3", all = T)  }
}
gb_table[is.na(gb_table)] <- 0
rownames(gb_table) = gb_table$V3
gb_table_trimmed <- gb_table[,-1]


# loading the bis table
y <- 0
for (x in bis_files) {
  y <- y + 1
  if (y == 1) {
    bis_table <- read.table(file = x, header=F, quote = "", sep = "\t")
    colnames(bis_table) = c("DELETE", x, "V3")
    bis_table <- bis_table[,c(2,3)]  }
  if (y > 1) {
    temp_table <- read.table(file = x, header = F, quote = "", sep = "\t")
    colnames(temp_table) = c("DELETE", x, "V3")
    temp_table <- temp_table[,c(2,3)]
    bis_table <- merge(bis_table, temp_table, by = "V3", all = T)  }
}
bis_table[is.na(bis_table)] <- 0
rownames(bis_table) = bis_table$V3
bis_table_trimmed <- bis_table[,-1]

# loading the nar table
y <- 0
for (x in nar_files) {
  y <- y + 1
  if (y == 1) {
    nar_table <- read.table(file = x, header=F, quote = "", sep = "\t")
    colnames(nar_table) = c("DELETE", x, "V3")
    nar_table <- nar_table[,c(2,3)]  }
  if (y > 1) {
    temp_table <- read.table(file = x, header = F, quote = "", sep = "\t")
    colnames(temp_table) = c("DELETE", x, "V3")
    temp_table <- temp_table[,c(2,3)]
    nar_table <- merge(nar_table, temp_table, by = "V3", all = T)  }
}
nar_table[is.na(nar_table)] <- 0
rownames(nar_table) = nar_table$V3
nar_table_trimmed <- nar_table[,-1]

# loading the nin table
y <- 0
for (x in nin_files) {
  y <- y + 1
  if (y == 1) {
    nin_table <- read.table(file = x, header=F, quote = "", sep = "\t")
    colnames(nin_table) = c("DELETE", x, "V3")
    nin_table <- nin_table[,c(2,3)]  }
  if (y > 1) {
    temp_table <- read.table(file = x, header = F, quote = "", sep = "\t")
    colnames(temp_table) = c("DELETE", x, "V3")
    temp_table <- temp_table[,c(2,3)]
    nin_table <- merge(nin_table, temp_table, by = "V3", all = T)  }
}
nin_table[is.na(nin_table)] <- 0
rownames(nin_table) = nin_table$V3
nin_table_trimmed <- nin_table[,-1]



# getting the column names simplified
colnames(pvd_table_trimmed) = pvd_names_trimmed
colnames(gb_table_trimmed) = gb_names_trimmed
colnames(bis_table_trimmed) = bis_names_trimmed
colnames(nar_table_trimmed) = nar_names_trimmed
colnames(nin_table_trimmed) = nin_names_trimmed
head(pvd_table_trimmed)


# merging the tables
#complete_tab1<-Reduce(function(x,y) merge(x,y,all=TRUE) ,list(pvd_table_trimmed, gb_table_trimmed,bis_table_trimmed, nar_table_trimmed,nin_table_trimmed), accumulate=FALSE)
complete_tab1<-merge(pvd_table_trimmed, gb_table_trimmed, by=0, all=TRUE)
complete_tab1[is.na(complete_tab1)] <- 0

# reducing stuff down to avoid duplicates
# removing extra Row.names column
complete_tab1 <- aggregate(. ~  Row.names, data = complete_tab1, sum)
complete_tab1 <- complete_tab1[!(complete_tab1$Row.names == ""), ]

bis_table_trimmed$Row.names<-row.names(bis_table_trimmed)
complete_tab2<-merge(complete_tab1, bis_table_trimmed, by="Row.names", all=TRUE)
complete_tab2[is.na(complete_tab2)] <- 0
complete_tab2 <- aggregate(. ~  Row.names, data = complete_tab2, sum)
complete_tab2 <- complete_tab2[!(complete_tab2$Row.names == ""), ]

nar_table_trimmed$Row.names<-row.names(nar_table_trimmed)
complete_tab3<-merge(complete_tab2, nar_table_trimmed, by="Row.names", all=TRUE)
complete_tab3[is.na(complete_tab3)] <- 0
complete_tab3 <- aggregate(. ~  Row.names, data = complete_tab3, sum)
complete_tab3 <- complete_tab3[!(complete_tab3$Row.names == ""), ]

nin_table_trimmed$Row.names<-row.names(nin_table_trimmed)
complete_table<-merge(complete_tab3, nin_table_trimmed, by="Row.names", all=TRUE)
complete_table[is.na(complete_table)] <- 0
complete_table <- aggregate(. ~  Row.names, data = complete_table, sum)
complete_table <- complete_table[!(complete_table$Row.names == ""), ]
row.names(complete_table)<-complete_table$Row.names
complete_table <- complete_table[,-1]
head(complete_table, n=15L)
names(complete_table)
print("complete_table merge successful \n\n")

# OPTIONAL: importing the raw counts
if (is.null(opt$raw_counts) == FALSE) {
  raw_counts_table <- read.table(counts_file, header=FALSE, sep = "\t", quote = "")
  raw_counts_table <- data.frame(raw_counts_table, 
        do.call(rbind, strsplit(as.character(raw_counts_table$V1),'_')))
  raw_counts_table$X2 <- as.numeric(as.character(raw_counts_table$X2))
  raw_counts_table <- t(raw_counts_table[,c("X2", "V2")])
  row.names(raw_counts_table) <- c("SAMPLE","RAW TOTAL")
  colnames(raw_counts_table) <- raw_counts_table[1,]
  raw_counts_table <- as.data.frame(raw_counts_table)
  raw_counts_table <- raw_counts_table[-1,]

  # Need to subtract off the total number of annotations
  raw_counts_table["ANNOTATION COUNT",] <- colSums(complete_table)
  raw_counts_table["OTHER",] <- raw_counts_table[1,] - raw_counts_table[2,]

  names(raw_counts_table)<-names(complete_table)
  complete_table <- rbind(complete_table, raw_counts_table["OTHER",])
}

head(complete_table, n=15L)


# DESeq statistical calculations
completeCondition <- data.frame(condition=factor(c(
  rep(paste("pvd", 1:length(pvd_files), sep=".")), 
  rep(paste("gb", 1:length(gb_files), sep=".")),
  rep(paste("bis", 1:length(bis_files), sep=".")), 
  rep(paste("nar", 1:length(nar_files), sep=".")), 
  rep(paste("nin", 1:length(nin_files), sep=".")))))
completeCondition1 <- t(completeCondition)
colnames(complete_table) <- completeCondition1
completeCondition2 <- data.frame(condition=factor(c(
  rep("pvd", length(pvd_files)), 
  rep("gb", length(gb_files)),
  rep("bis", length(bis_files)),
  rep("nin", length(nar_files)),
  rep("nar", length(nin_files)))))


dds <- DESeqDataSetFromMatrix(complete_table, completeCondition2, ~condition)
dds <- DESeq(dds)
head(dds)

baseMeanPerLvl <- sapply( levels(dds$condition), 
                          function(lvl) rowMeans( counts(dds,normalized=TRUE)[,dds$condition == lvl] ) )
head(baseMeanPerLvl)

# This step creates the summary results output
#res <- results(dds, contrast=list("pvd",c("gb","bis","nar","nin")), 
#               listValues=c(1, -1/4))
res_pvd <- results(dds, contrast=list("conditionpvd",c("conditiongb","conditionbis","conditionnar","conditionnin")), 
                   listValues=c(1, -1/4))
res_gb <- results(dds, contrast=list("conditiongb",c("conditionpvd","conditionbis","conditionnar","conditionnin")), 
                   listValues=c(1, -1/4))
res_bis <- results(dds, contrast=list("conditionbis",c("conditiongb","conditionpvd","conditionnar","conditionnin")), 
                   listValues=c(1, -1/4))
res_nar <- results(dds, contrast=list("conditionnar",c("conditiongb","conditionbis","conditionpvd","conditionnin")), 
                   listValues=c(1, -1/4))
res_nin <- results(dds, contrast=list("conditionnin",c("conditiongb","conditionbis","conditionnar","conditionpvd")), 
                   listValues=c(1, -1/4))
head(res_pvd)


res1 <- results(dds, contrast=c("condition","pvd","gb"))
head(res1)
res2 <- results(dds, contrast=c("condition","pvd","bis"))
head(res2)
res3 <- results(dds, contrast=c("condition", "pvd","nar"))
head(res3)
res4 <- results(dds, contrast=c("condition", "pvd","nin"))
head(res4)
res5 <- results(dds, contrast=c("condition", "gb","bis"))
head(res5)
res6 <- results(dds, contrast=c("condition", "gb","nar"))
head(res6)
res7 <- results(dds, contrast=c("condition", "gb","nin"))
head(res7)
res8 <- results(dds, contrast=c("condition", "bis","nar"))
head(res8)
res9 <- results(dds, contrast=c("condition", "bis","nin"))
head(res9)
res10 <- results(dds, contrast=c("condition", "nar","nin"))
head(res10)


org_results_pvd <- data.frame(res_pvd)
org_results_pvd <- merge(org_results_pvd, baseMeanPerLvl, by="row.names")
org_results_pvd <- org_results_pvd[,c(1,2,8,9,3,4,5,6,7)]
colnames(org_results_pvd)[c(3,4)] <- c("pvd", "allothers")
sorted_org_results_pvd <- org_results_pvd[order(-org_results_pvd$baseMean),]
colnames(sorted_org_results_pvd)[1] <- "Organism Name"

org_results_gb <- data.frame(res_gb)
org_results_gb <- merge(org_results_gb, baseMeanPerLvl, by="row.names")
org_results_gb <- org_results_gb[,c(1,2,8,9,3,4,5,6,7)]
colnames(org_results_gb)[c(3,4)] <- c("gb", "allothers")
sorted_org_results_gb <- org_results_gb[order(-org_results_gb$baseMean),]
colnames(sorted_org_results_gb)[1] <- "Organism Name"

org_results_bis <- data.frame(res_bis)
org_results_bis <- merge(org_results_bis, baseMeanPerLvl, by="row.names")
org_results_bis <- org_results_bis[,c(1,2,8,9,3,4,5,6,7)]
colnames(org_results_bis)[c(3,4)] <- c("bis", "allothers")
sorted_org_results_bis <- org_results_bis[order(-org_results_bis$baseMean),]
colnames(sorted_org_results_bis)[1] <- "Organism Name"

org_results_nar <- data.frame(res_nar)
org_results_nar <- merge(org_results_nar, baseMeanPerLvl, by="row.names")
org_results_nar <- org_results_nar[,c(1,2,8,9,3,4,5,6,7)]
colnames(org_results_nar)[c(3,4)] <- c("nar", "allothers")
sorted_org_results_nar <- org_results_nar[order(-org_results_nar$baseMean),]
colnames(sorted_org_results_nar)[1] <- "Organism Name"

org_results_nin <- data.frame(res_nin)
org_results_nin <- merge(org_results_nin, baseMeanPerLvl, by="row.names")
org_results_nin <- org_results_nin[,c(1,2,8,9,3,4,5,6,7)]
colnames(org_results_nin)[c(3,4)] <- c("nin", "allothers")
sorted_org_results_nin <- org_results_nin[order(-org_results_nin$baseMean),]
colnames(sorted_org_results_nin)[1] <- "Organism Name"



org_results1 <- data.frame(res1)
org_results1 <- merge(org_results1, baseMeanPerLvl, by="row.names")
org_results1 <- org_results1[,c(1,2,8,9,3,4,5,6,7)]
colnames(org_results1)[c(3,4)] <- c("pvd", "gb")
sorted_org_results1 <- org_results1[order(-org_results1$baseMean),]
colnames(sorted_org_results1)[1] <- "Organism Name"

org_results2 <- data.frame(res2)
org_results2 <- merge(org_results2, baseMeanPerLvl, by="row.names")
org_results2 <- org_results2[,c(1,2,8,9,3,4,5,6,7)]
colnames(org_results2)[c(3,4)] <- c("pvd", "bis")
sorted_org_results2 <- org_results2[order(-org_results2$baseMean),]
colnames(sorted_org_results2)[1] <- "Organism Name"

org_results3 <- data.frame(res3)
org_results3 <- merge(org_results3, baseMeanPerLvl, by="row.names")
org_results3 <- org_results3[,c(1,2,8,9,3,4,5,6,7)]
colnames(org_results3)[c(3,4)] <- c("pvd", "nar")
sorted_org_results3 <- org_results3[order(-org_results3$baseMean),]
colnames(sorted_org_results3)[1] <- "Organism Name"

org_results4 <- data.frame(res4)
org_results4 <- merge(org_results4, baseMeanPerLvl, by="row.names")
org_results4 <- org_results4[,c(1,2,8,9,3,4,5,6,7)]
colnames(org_results4)[c(3,4)] <- c("pvd", "nin")
sorted_org_results4 <- org_results4[order(-org_results4$baseMean),]
colnames(sorted_org_results4)[1] <- "Organism Name"

org_results5 <- data.frame(res5)
org_results5 <- merge(org_results5, baseMeanPerLvl, by="row.names")
org_results5 <- org_results5[,c(1,2,8,9,3,4,5,6,7)]
colnames(org_results5)[c(3,4)] <- c("gb", "bis")
sorted_org_results5 <- org_results5[order(-org_results5$baseMean),]
colnames(sorted_org_results5)[1] <- "Organism Name"

org_results6 <- data.frame(res6)
org_results6 <- merge(org_results6, baseMeanPerLvl, by="row.names")
org_results6 <- org_results6[,c(1,2,8,9,3,4,5,6,7)]
colnames(org_results6)[c(3,4)] <- c("gb", "nar")
sorted_org_results6 <- org_results6[order(-org_results6$baseMean),]
colnames(sorted_org_results6)[1] <- "Organism Name"

org_results7 <- data.frame(res7)
org_results7 <- merge(org_results7, baseMeanPerLvl, by="row.names")
org_results7 <- org_results7[,c(1,2,8,9,3,4,5,6,7)]
colnames(org_results7)[c(3,4)] <- c("gb", "nin")
sorted_org_results7 <- org_results7[order(-org_results7$baseMean),]
colnames(sorted_org_results7)[1] <- "Organism Name"

org_results8 <- data.frame(res8)
org_results8 <- merge(org_results8, baseMeanPerLvl, by="row.names")
org_results8 <- org_results8[,c(1,2,8,9,3,4,5,6,7)]
colnames(org_results8)[c(3,4)] <- c("bis", "nar")
sorted_org_results8 <- org_results8[order(-org_results8$baseMean),]
colnames(sorted_org_results8)[1] <- "Organism Name"

org_results9 <- data.frame(res9)
org_results9 <- merge(org_results9, baseMeanPerLvl, by="row.names")
org_results9 <- org_results9[,c(1,2,8,9,3,4,5,6,7)]
colnames(org_results9)[c(3,4)] <- c("bis", "nin")
sorted_org_results9 <- org_results9[order(-org_results9$baseMean),]
colnames(sorted_org_results9)[1] <- "Organism Name"

org_results10 <- data.frame(res10)
org_results10 <- merge(org_results10, baseMeanPerLvl, by="row.names")
org_results10 <- org_results10[,c(1,2,8,9,3,4,5,6,7)]
colnames(org_results10)[c(3,4)] <- c("nar", "nin")
sorted_org_results10 <- org_results10[order(-org_results6$baseMean),]
colnames(sorted_org_results10)[1] <- "Organism Name"


# saving and finishing up
write.table(sorted_org_results_pvd, file = paste(save_filename,"_pvd"), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(sorted_org_results_gb, file = paste(save_filename,"_gb"), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(sorted_org_results_bis, file = paste(save_filename,"_bis"), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(sorted_org_results_nar, file = paste(save_filename,"_nar"), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(sorted_org_results_nin, file = paste(save_filename,"_nin"), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

write.table(sorted_org_results1, file = paste(save_filename,"_01"), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(sorted_org_results2, file = paste(save_filename,"_02"), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(sorted_org_results3, file = paste(save_filename,"_03"), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(sorted_org_results4, file = paste(save_filename,"_04"), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(sorted_org_results5, file = paste(save_filename,"_05"), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(sorted_org_results6, file = paste(save_filename,"_06"), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(sorted_org_results7, file = paste(save_filename,"_07"), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(sorted_org_results8, file = paste(save_filename,"_08"), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(sorted_org_results9, file = paste(save_filename,"_09"), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(sorted_org_results10, file = paste(save_filename,"_10"), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


cat ("\nSuccess!\nSaved results file as ", save_filename, "\n")



# heatmap ---------------------------------------------------------------
library(pheatmap)

rowVars<- function (x,na.rm = TRUE) {
  sqr = function(x) x * x
  n = rowSums(!is.na(x))
  n[n <= 1] = NA
  return(rowSums(sqr(x - rowMeans(x,na.rm = na.rm)), na.rm = na.rm)/(n - 1))
}

transformed_data <- rlog(dds, blind=FALSE)
head(transformed_data)
topVarGenes <- head(order(rowVars(assay(transformed_data)), decreasing=T), 50)
head(topVarGenes)

# making the PCA plot

# calculate euclidean distances from the variance-stabilized data
dists <- dist(t(assay(transformed_data)))
PCAplot <- plotPCA(transformed_data, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(PCAplot, "percentVar"))

# making the PCA plot -----------------------------
# calculate euclidean distances from the variance-stabilized data
matrix <- assay(transformed_data)[ topVarGenes, ]
matrix <- matrix - rowMeans(matrix)
write.table(matrix, file = "heatmapmatrix.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)


# Create heatmap of distances ---------------------
#colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap1 <- heatmap(matrix)
heatmap <- pheatmap(matrix)

# saving and finishing up
cat ("Saving heatmap as ", paste(save_filename,"_heat"), " now.\n")
pdf(file = save_filename, width=10, height=10)
heatmap
heatmap1
ggplot(PCAplot, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  #    geom_text(aes(label=name), hjust=1, vjust=-1) +
  ggtitle("PCA Plot of data") +
  theme(legend.position = "bottom") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()




# Diversity plots -------------------------------

# getting diversity statistics
flipped_complete_table <- data.frame(t(complete_table))

graphing_table <- data.frame(condition=factor(c(
  rep("pvd", length(pvd_files)), 
  rep("gb", length(gb_files)),
  rep("bis", length(bis_files)),
  rep("nin", length(nar_files)),
  rep("nar", length(nin_files)))))
graphing_table[,"order"] <- c(1:nrow(graphing_table))
graphing_table[,"Shannon"] <- diversity(flipped_complete_table, index = "shannon")
graphing_table[,"Simpson"] <- diversity(flipped_complete_table, index = "simpson")

shannon_plot <- ggplot(data = graphing_table, aes(x=order, y=Shannon,
                                                  color = condition, fill = condition)) +
  geom_bar(stat="identity", width = 0.8) +
  ggtitle("Shannon diversity of samples") +
  theme(legend.position = "bottom")

simpson_plot <- ggplot(data = graphing_table, aes(x=order, y=Simpson,
                                                  color = condition, fill = condition)) +
  geom_bar(stat="identity", width = 0.8) +
  ggtitle("Simpson diversity of csamples") +
  theme(legend.position = "bottom")

cat ("\nSuccess!\nSaving diversity graphs as ", paste(save_filename,"_diversity"), " now.\n")
pdf(file = save_filename, width=10, height=7)
grid.arrange(shannon_plot, simpson_plot, ncol=1)
dev.off()


