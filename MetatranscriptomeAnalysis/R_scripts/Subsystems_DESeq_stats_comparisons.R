# Subsystems_DESeq_stats.R
# Created 3/06/2017, by Sam Westreich
# Last updated 6/16/2017
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
              help="raw (total) read counts for this starting file", metavar="character"),
  make_option(c("-L", "--level"), type="integer", default=1,
              help="level of Subsystems hierarchy for DESeq stats [default=%default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print("USAGE: $ Subsystems_DESeq_stats.R -I input_directory/ -O save.filename -L level (1,2,3,4) [-R raw_counts_file]")

# check for necessary specs
if (is.null(opt$input)) {
  print ("WARNING: No working input directory specified with '-I' flag.")
  stop()
} else {  cat ("Working directory is ", opt$input, "\n")
  wd_location <- opt$input  
  setwd(wd_location)  }

cat ("Saving results as ", opt$out, "\n")
save_filename <- opt$out

if (is.null(opt$raw_counts)) {
  print ("WARNING: no raw counts file specified, skipping this info for DESeq analysis.")
} else {
  counts_file <- opt$raw_counts
}

cat ("Calculating DESeq results for hierarchy level ", opt$level, "\n")
levelname=paste("Level", opt$level, sep="")

# import other necessary packages
suppressPackageStartupMessages({
  library(DESeq2)
  library("data.table")
  library(plyr)
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



# loading the pvd files in
y <- 0
for (x in pvd_files) {
  y <- y + 1
  if (y == 1) {
    pvd_table <- read.table(file = x, header = F, quote = "", sep = "\t", fill = TRUE)
    colnames(pvd_table) = c("DELETE", x, "Level4", "Level3", "Level1", "Level2")
    if (opt$level == 1) {
      pvd_table = pvd_table[,c(2, 5)]
    } else if (opt$level == 2) {
      pvd_table = pvd_table[,c(2, 6)]
    } else if (opt$level == 3) {
      pvd_table = pvd_table[,c(2, 4)]
    } else {
      pvd_table = pvd_table[,c(2, 3)]
    }
    pvd_table <- ddply(pvd_table, colnames(pvd_table)[2], numcolwise(sum))
    rownames(pvd_table) <- pvd_table[,1] }
  if (y > 1) {
    temp_table <- read.table(file = x, header = F, quote = "", sep = "\t", fill = TRUE)
    colnames(temp_table) = c("DELETE", x, "Level4", "Level3", "Level1", "Level2")
    if (opt$level == 1) {
      temp_table = temp_table[,c(2, 5)]
    } else if (opt$level == 2) {
      temp_table = temp_table[,c(2, 6)]
    } else if (opt$level == 3) {
      temp_table = temp_table[,c(2, 4)]
    } else {
      temp_table = temp_table[,c(2, 3)]
    }
    temp_table <- ddply(temp_table, colnames(temp_table)[2], numcolwise(sum))
    rownames(temp_table) <- temp_table[,1]
    pvd_table <- merge(temp_table, pvd_table, by = colnames(temp_table)[1], all=T)
  }
}
pvd_table <- pvd_table[!is.na(names(pvd_table))]
pvd_table[is.na(pvd_table)] <- ""

# Need to convert NAs to 0s
pvd_data <- pvd_table[,c(2:(length(pvd_files)+1))]
pvd_data <- lapply(pvd_data, function(x) as.numeric(as.character(x)))
pvd_data <- as.data.frame(pvd_data)
pvd_data[is.na(pvd_data)] <- 0
pvd_table[,c(2:(length(pvd_files)+1))] <- pvd_data




# loading the gb files in
y <- 0
for (x in gb_files) {
  y <- y + 1
  if (y == 1) {
    gb_table <- read.table(file = x, header = F, quote = "", sep = "\t", fill = TRUE)
    colnames(gb_table) = c("DELETE", x, "Level4", "Level3", "Level1", "Level2")
    if (opt$level == 1) {
      gb_table = gb_table[,c(2, 5)]
    } else if (opt$level == 2) {
      gb_table = gb_table[,c(2, 6)]
    } else if (opt$level == 3) {
      gb_table = gb_table[,c(2, 4)]
    } else {
      gb_table = gb_table[,c(2, 3)]
    }
    gb_table <- ddply(gb_table, colnames(gb_table)[2], numcolwise(sum))
    rownames(gb_table) <- gb_table[,1] }
  if (y > 1) {
    temp_table <- read.table(file = x, header = F, quote = "", sep = "\t", fill = TRUE)
    colnames(temp_table) = c("DELETE", x, "Level4", "Level3", "Level1", "Level2")
    if (opt$level == 1) {
      temp_table = temp_table[,c(2, 5)]
    } else if (opt$level == 2) {
      temp_table = temp_table[,c(2, 6)]
    } else if (opt$level == 3) {
      temp_table = temp_table[,c(2, 4)]
    } else {
      temp_table = temp_table[,c(2, 3)]
    }
    temp_table <- ddply(temp_table, colnames(temp_table)[2], numcolwise(sum))
    rownames(temp_table) <- temp_table[,1]
    gb_table <- merge(temp_table, gb_table, by = colnames(temp_table)[1], all=T)
  }
}
gb_table <- gb_table[!is.na(names(gb_table))]
gb_table[is.na(gb_table)] <- ""

# converting NAs to 0s
gb_data <- gb_table[,c(2:(length(gb_files)+1))]
gb_data <- lapply(gb_data, function(x) as.numeric(as.character(x)))
gb_data <- as.data.frame(gb_data)
gb_data[is.na(gb_data)] <- 0
gb_table[,c(2:(length(gb_files)+1))] <- gb_data




# loading the bis files in
y <- 0
for (x in bis_files) {
  y <- y + 1
  if (y == 1) {
    bis_table <- read.table(file = x, header = F, quote = "", sep = "\t", fill = TRUE)
    colnames(bis_table) = c("DELETE", x, "Level4", "Level3", "Level1", "Level2")
    if (opt$level == 1) {
      bis_table = bis_table[,c(2, 5)]
    } else if (opt$level == 2) {
      bis_table = bis_table[,c(2, 6)]
    } else if (opt$level == 3) {
      bis_table = bis_table[,c(2, 4)]
    } else {
      bis_table = bis_table[,c(2, 3)]
    }
    bis_table <- ddply(bis_table, colnames(bis_table)[2], numcolwise(sum))
    rownames(bis_table) <- bis_table[,1] }
  if (y > 1) {
    temp_table <- read.table(file = x, header = F, quote = "", sep = "\t", fill = TRUE)
    colnames(temp_table) = c("DELETE", x, "Level4", "Level3", "Level1", "Level2")
    if (opt$level == 1) {
      temp_table = temp_table[,c(2, 5)]
    } else if (opt$level == 2) {
      temp_table = temp_table[,c(2, 6)]
    } else if (opt$level == 3) {
      temp_table = temp_table[,c(2, 4)]
    } else {
      temp_table = temp_table[,c(2, 3)]
    }
    temp_table <- ddply(temp_table, colnames(temp_table)[2], numcolwise(sum))
    rownames(temp_table) <- temp_table[,1]
    bis_table <- merge(temp_table, bis_table, by = colnames(temp_table)[1], all=T)
  }
}
bis_table <- bis_table[!is.na(names(bis_table))]
bis_table[is.na(bis_table)] <- ""

# converting NAs to 0s
bis_data <- bis_table[,c(2:(length(bis_files)+1))]
bis_data <- lapply(bis_data, function(x) as.numeric(as.character(x)))
bis_data <- as.data.frame(bis_data)
bis_data[is.na(bis_data)] <- 0
bis_table[,c(2:(length(bis_files)+1))] <- bis_data





# loading the nar files in
y <- 0
for (x in nar_files) {
  y <- y + 1
  if (y == 1) {
    nar_table <- read.table(file = x, header = F, quote = "", sep = "\t", fill = TRUE)
    colnames(nar_table) = c("DELETE", x, "Level4", "Level3", "Level1", "Level2")
    if (opt$level == 1) {
      nar_table = nar_table[,c(2, 5)]
    } else if (opt$level == 2) {
      nar_table = nar_table[,c(2, 6)]
    } else if (opt$level == 3) {
      nar_table = nar_table[,c(2, 4)]
    } else {
      nar_table = nar_table[,c(2, 3)]
    }
    nar_table <- ddply(nar_table, colnames(nar_table)[2], numcolwise(sum))
    rownames(nar_table) <- nar_table[,1] }
  if (y > 1) {
    temp_table <- read.table(file = x, header = F, quote = "", sep = "\t", fill = TRUE)
    colnames(temp_table) = c("DELETE", x, "Level4", "Level3", "Level1", "Level2")
    if (opt$level == 1) {
      temp_table = temp_table[,c(2, 5)]
    } else if (opt$level == 2) {
      temp_table = temp_table[,c(2, 6)]
    } else if (opt$level == 3) {
      temp_table = temp_table[,c(2, 4)]
    } else {
      temp_table = temp_table[,c(2, 3)]
    }
    temp_table <- ddply(temp_table, colnames(temp_table)[2], numcolwise(sum))
    rownames(temp_table) <- temp_table[,1]
    nar_table <- merge(temp_table, nar_table, by = colnames(temp_table)[1], all=T)
  }
}
nar_table <- nar_table[!is.na(names(nar_table))]
nar_table[is.na(nar_table)] <- ""

# converting NAs to 0s
nar_data <- nar_table[,c(2:(length(nar_files)+1))]
nar_data <- lapply(nar_data, function(x) as.numeric(as.character(x)))
nar_data <- as.data.frame(nar_data)
nar_data[is.na(nar_data)] <- 0
nar_table[,c(2:(length(nar_files)+1))] <- nar_data






# loading the nin files in
y <- 0
for (x in nin_files) {
  y <- y + 1
  if (y == 1) {
    nin_table <- read.table(file = x, header = F, quote = "", sep = "\t", fill = TRUE)
    colnames(nin_table) = c("DELETE", x, "Level4", "Level3", "Level1", "Level2")
    if (opt$level == 1) {
      nin_table = nin_table[,c(2, 5)]
    } else if (opt$level == 2) {
      nin_table = nin_table[,c(2, 6)]
    } else if (opt$level == 3) {
      nin_table = nin_table[,c(2, 4)]
    } else {
      nin_table = nin_table[,c(2, 3)]
    }
    nin_table <- ddply(nin_table, colnames(nin_table)[2], numcolwise(sum))
    rownames(nin_table) <- nin_table[,1] }
  if (y > 1) {
    temp_table <- read.table(file = x, header = F, quote = "", sep = "\t", fill = TRUE)
    colnames(temp_table) = c("DELETE", x, "Level4", "Level3", "Level1", "Level2")
    if (opt$level == 1) {
      temp_table = temp_table[,c(2, 5)]
    } else if (opt$level == 2) {
      temp_table = temp_table[,c(2, 6)]
    } else if (opt$level == 3) {
      temp_table = temp_table[,c(2, 4)]
    } else {
      temp_table = temp_table[,c(2, 3)]
    }
    temp_table <- ddply(temp_table, colnames(temp_table)[2], numcolwise(sum))
    rownames(temp_table) <- temp_table[,1]
    nin_table <- merge(temp_table, nin_table, by = colnames(temp_table)[1], all=T)
  }
}
nin_table <- nin_table[!is.na(names(nin_table))]
nin_table[is.na(nin_table)] <- ""

# converting NAs to 0s
nin_data <- nin_table[,c(2:(length(nin_files)+1))]
nin_data <- lapply(nin_data, function(x) as.numeric(as.character(x)))
nin_data <- as.data.frame(nin_data)
nin_data[is.na(nin_data)] <- 0
nin_table[,c(2:(length(nin_files)+1))] <- nin_data




# reducing stuff down to avoid duplicates
colnames(pvd_table) <- c(colnames(pvd_table)[1], pvd_names_trimmed)
colnames(gb_table) <- c(colnames(gb_table)[1], gb_names_trimmed)
colnames(bis_table) <- c(colnames(bis_table)[1], bis_names_trimmed)
colnames(nar_table) <- c(colnames(nar_table)[1], nar_names_trimmed)
colnames(nin_table) <- c(colnames(nin_table)[1], nin_names_trimmed)


l1_table<-Reduce(function(x,y) merge(x,y, by=colnames(pvd_table)[1], all.x = T) ,list(pvd_table,gb_table,bis_table,nar_table,nin_table), accumulate=FALSE)
l1_table[,levelname][is.na(l1_table[,levelname])] <- ""
l1_table[,levelname] <- sub("^$", "NO HIERARCHY", l1_table[,levelname])
l1_table <- ddply(l1_table, colnames(l1_table)[1], numcolwise(sum))
rownames(l1_table) <- l1_table[,levelname]
l1_names <- l1_table[,levelname]
l1_table[,levelname] <- NULL
l1_table[is.na(l1_table)] <- 0
print("complete_table merge successful \n\n")
head(l1_table)






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
  raw_counts_table["ANNOTATION COUNT",] <- colSums(l1_table)
  raw_counts_table["OTHER",] <- raw_counts_table[1,] - raw_counts_table[2,]
  
  names(raw_counts_table)<-names(l1_table)
  l1_table <- rbind(l1_table, raw_counts_table["OTHER",])
}

head(l1_table, n=15L)



# now the DESeq stuff ---------------------------------------------------------------
cat ("Now running DESeq.\n")
completeCondition <- data.frame(condition=factor(c(
  rep("pvd", length(pvd_files)), 
  rep("gb", length(gb_files)),
  rep("bis", length(bis_files)),
  rep("nin", length(nar_files)),
  rep("nar", length(nin_files)))))
dds <- DESeqDataSetFromMatrix(l1_table, completeCondition, ~ condition)
dds <- DESeq(dds)
baseMeanPerLvl <- sapply( levels(dds$condition), function(lvl) rowMeans( 
  counts(dds,normalized=TRUE)[,dds$condition == lvl] ) )


res_pvd <- results(dds, contrast=list("conditionpvd",c("conditionnin")), 
                   listValues=c(1, -1))
res_gb <- results(dds, contrast=list("conditiongb",c("conditionnin")), 
                  listValues=c(1, -1))
res_bis <- results(dds, contrast=list("conditionbis",c("conditionnin")), 
                   listValues=c(1, -1))
res_nar <- results(dds, contrast=list("conditionnar",c("conditionnin")), 
                   listValues=c(1, -1))
res_nin <- results(dds, contrast=list(c("conditiongb","conditionbis","conditionpvd"),c("conditionnar","conditionnin")), 
                   listValues=c(1/2, -1/2))
head(res_nin)

l1_results <- data.frame(res_pvd)
head(l1_results)
rownames(l1_results) <- rownames(l1_table)
rownames(baseMeanPerLvl) <- rownames(l1_table)
l1_results <- merge(l1_results, baseMeanPerLvl, by="row.names")
l1_results <- l1_results[,c(1,2,8,9,3,4,5,6,7)]
#colnames(l1_results)[c(3,4)] <- c("controlMean", "experimentalMean")
l1_results_pvd <- l1_results[order(-l1_results$baseMean),]

l1_results <- data.frame(res_gb)
head(l1_results)
rownames(l1_results) <- rownames(l1_table)
rownames(baseMeanPerLvl) <- rownames(l1_table)
l1_results <- merge(l1_results, baseMeanPerLvl, by="row.names")
l1_results <- l1_results[,c(1,2,8,9,3,4,5,6,7)]
#colnames(l1_results)[c(3,4)] <- c("controlMean", "experimentalMean")
l1_results_gb <- l1_results[order(-l1_results$baseMean),]

l1_results <- data.frame(res_bis)
head(l1_results)
rownames(l1_results) <- rownames(l1_table)
rownames(baseMeanPerLvl) <- rownames(l1_table)
l1_results <- merge(l1_results, baseMeanPerLvl, by="row.names")
l1_results <- l1_results[,c(1,2,8,9,3,4,5,6,7)]
#colnames(l1_results)[c(3,4)] <- c("controlMean", "experimentalMean")
l1_results_bis <- l1_results[order(-l1_results$baseMean),]

l1_results <- data.frame(res_nar)
head(l1_results)
rownames(l1_results) <- rownames(l1_table)
rownames(baseMeanPerLvl) <- rownames(l1_table)
l1_results <- merge(l1_results, baseMeanPerLvl, by="row.names")
l1_results <- l1_results[,c(1,2,8,9,3,4,5,6,7)]
#colnames(l1_results)[c(3,4)] <- c("controlMean", "experimentalMean")
l1_results_nar <- l1_results[order(-l1_results$baseMean),]

l1_results <- data.frame(res_nin)
head(l1_results)
rownames(l1_results) <- rownames(l1_table)
rownames(baseMeanPerLvl) <- rownames(l1_table)
l1_results <- merge(l1_results, baseMeanPerLvl, by="row.names")
l1_results <- l1_results[,c(1,2,8,9,3,4,5,6,7)]
#colnames(l1_results)[c(3,4)] <- c("controlMean", "experimentalMean")
l1_results_nin <- l1_results[order(-l1_results$baseMean),]

# saving and finishing up
cat ("\nSuccess!\nSaving results file as ", save_filename, "\n")
write.table(l1_results_pvd, file = paste(save_filename,"_pvd"),  append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(l1_results_gb, file = paste(save_filename,"_gb"),  append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(l1_results_bis, file = paste(save_filename,"_bis"),  append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(l1_results_nar, file = paste(save_filename,"_nar"),  append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(l1_results_nin, file = paste(save_filename,"_south"),  append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
# 
# 
# # heatmap ---------------------------------------------------------------
# library(pheatmap)
# 
# rowVars<- function (x,na.rm = TRUE) {
#   sqr = function(x) x * x
#   n = rowSums(!is.na(x))
#   n[n <= 1] = NA
#   return(rowSums(sqr(x - rowMeans(x,na.rm = na.rm)), na.rm = na.rm)/(n - 1))
# }
# 
# transformed_data <- rlog(dds, blind=FALSE)
# head(transformed_data)
# topVarGenes <- head(order(rowVars(assay(transformed_data)), decreasing=T), 50)
# head(topVarGenes)
# 
# # making the PCA plot
# 
# # calculate euclidean distances from the variance-stabilized data
# dists <- dist(t(assay(transformed_data)))
# PCAplot <- plotPCA(transformed_data, intgroup = "condition", returnData = TRUE)
# percentVar <- round(100 * attr(PCAplot, "percentVar"))
# 
# # making the PCA plot -----------------------------
# # calculate euclidean distances from the variance-stabilized data
# matrix <- assay(transformed_data)[ topVarGenes, ]
# matrix <- matrix - rowMeans(matrix)
# write.table(matrix, file = "heatmapmatrix.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
# 
# 
# # Create heatmap of distances ---------------------
# #colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# heatmap1 <- heatmap(matrix)
# heatmap <- pheatmap(matrix)
# 
# # saving and finishing up
# cat ("Saving heatmap as ", paste(save_filename,"_heat"), " now.\n")
# pdf(file = save_filename, width=10, height=10)
# heatmap
# heatmap1
# plot(PCAplot$PC1, PCAplot$PC2)
# dev.off()
# 
# cat ("Saving heatmap as png", paste(save_filename,"_heat"), " now.\n")
# png(file = save_filename)
# heatmap
# heatmap1
# plot(PCAplot$PC1, PCAplot$PC2)
# dev.off()



