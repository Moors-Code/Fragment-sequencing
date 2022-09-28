########## Part 1.1: Pre-processing - MULTI-seq demultiplexing ##########
#This part uses the workflow from https://github.com/chris-mcginnis-ucsf/MULTI-seq to demultiplex FASTQ files of MULTIseq libraries 
#to allocate single cells to their MULTI-seq BC and thereby sphere of origin 

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")
source("./functions_and_packages/2.Functions_preprocessing.R")

########## Liver samples ##########
###Sample 1
##UMI count table generation 
bar.table.start <- MULTIseq_UMI_count_matrix_generation_BD(
  "./data/MULTI-seq_3Prime_3Plates_barcodes.xlsx",
  "./data/MetsM1kept_barcodes.txt",
  "/home/khandler/NAS/Kristina/Data_paper/LiverM1/Mets1MS_merged_R1.fastq.gz",
  "/home/khandler/NAS/Kristina/Data_paper/LiverM1/Mets1MS_merged_R2.fastq.gz")

##Extract only columns with barcode information
bar.table.full <- bar.table.start[, 1:288] 

##Visualize UMI counts in UMAP space for each individual MULTI-seq barcode and extract only positive ones (indicated by clusters in color red)
bar.tsne <- barTSNE(bar.table.start[,1:288]) 
for (i in 3:ncol(bar.tsne)) {
  g <- ggplot(bar.tsne, aes(x = TSNE1, y = TSNE2, color = bar.tsne[,i])) +
    geom_point() +
    scale_color_gradient(low = "black", high = "red") +
    ggtitle(colnames(bar.tsne)[i]) +
    theme(legend.position = "none") 
  print(g)
}

good.bars <- paste("Bar",c(2:7,9:14,17,19:23,26,28,30:36,39,41:43,45:48,62,64,65,66,70,72,73,
                           74,76:84,86:96,99:102,104,106,111,112,116,117,119,120,122,124,125,
                           127:129,134,136,137,139,142,144:149,151,152,154,156:161,163:168,
                           171:174,176,177,180:184,186:191,193,196,197,201:203,206,208:210,213,
                           214,216,222,226,227,234,237:240,242,243,245:252,258,261,265,268,272,
                           275,276,285,288),sep="") 
bar.table <- bar.table.start[, good.bars]  

##MULTI-seq BC Classification 
#1st round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- names(round.calls)[which(round.calls == "Negative")]
#Remove negative cells 
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#2nd round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#3rd round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#4th round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#5th round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#6th round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls) #no negative cells left 
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])

##Combine final calls and merge with Cell IDs in bar.table.full 
final.calls <- c(round.calls, rep("Negative", length(neg.cells)))
names(final.calls) <- c(names(round.calls), neg.cells)
#remove duplicates 
final.calls <- final.calls[!duplicated(names(final.calls))]

#make final table with matched Cell ID and MULTI-seq BC information 
final_round_classification <- as.data.frame(final.calls, row.names = names(final.calls))
classified_table <- merge(bar.table.full, final_round_classification, by = "row.names")
rownames(classified_table) <- rownames(bar.table.full)
classified_table$Barcode <- classified_table$final.calls

##save table 
saveRDS(classified_table, file = "./data_files_generated/classified_table_Liver_sample1.Rda")


###Sample 2
##UMI count table generation 
bar.table.start <- MULTIseq_UMI_count_matrix_generation_BD(
  "./data/MULTI-seq_3Prime_3Plates_barcodes.xlsx",
  "./data/MetsM3kept_barcodes.txt",
  "/home/khandler/NAS/Kristina/Data_paper/LiverM3/Mets3MS_merged_R1.fastq.gz",
  "/home/khandler/NAS/Kristina/Data_paper/LiverM3/Mets3MS_merged_R2.fastq.gz")

##Extract only columns with barcode information 
bar.table.full <- bar.table.start[, 1:288] 

##Visualize UMI counts in UMAP space for each individual MULTI-seq barcode and extract only positive ones (indicated by clusters in color red)
bar.tsne <- barTSNE(bar.table.start[,1:288]) 
for (i in 3:ncol(bar.tsne)) {
  g <- ggplot(bar.tsne, aes(x = TSNE1, y = TSNE2, color = bar.tsne[,i])) +
    geom_point() +
    scale_color_gradient(low = "black", high = "red") +
    ggtitle(colnames(bar.tsne)[i]) +
    theme(legend.position = "none") 
  print(g)
}

good.bars <- paste("Bar",c(2,3,5:14,18,20:23,25,26,28,30:33,35,36,38,41:42,44:48,50:54,56:60,62,63,65,67:76,
                           78:80,83,84,86,87,89,92:96,98,101:103,105:110,113,118:120,123:124,123,
                           124,126,127,130:133,135,137:142,144:147,151:153,155,157,161,163:165,168,
                           170:172,174:176,179,181,182,184,186:188,190,191,193,195,197:208,210:212,
                           214:229,231:239,241:244,246,248:255,258:268,270,273,275:277,279,280,282:285,287),sep="") 
bar.table <- bar.table.start[, good.bars]  

##MULTIseq BC Classification 
#1st round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- names(round.calls)[which(round.calls == "Negative")]
#Remove negative cells 
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#2nd round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#3rd round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls) #no negative cells left 
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])

##Combine final calls and merge with Cell IDs in bar.table.full 
final.calls <- c(round.calls, rep("Negative", length(neg.cells)))
names(final.calls) <- c(names(round.calls), neg.cells)
#remove duplicates 
final.calls <- final.calls[!duplicated(names(final.calls))]

#make final table with matched Cell ID and MULTI-seq BC information 
final_round_classification <- as.data.frame(final.calls, row.names = names(final.calls))
classified_table <- merge(bar.table.full, final_round_classification, by = "row.names")
rownames(classified_table) <- rownames(bar.table.full)
classified_table$Barcode <- classified_table$final.calls

##save table 
saveRDS(classified_table, file = "./data_files_generated/classified_table_Liver_sample2.Rda")

###Sample 3
##UMI count table generation 
bar.table.start <- MULTIseq_UMI_count_matrix_generation_BD(
  "./data/MULTI-seq_3Prime_3Plates_barcodes.xlsx",
  "./data/MetsM4kept_barcodes.txt",
  "/home/khandler/NAS/Kristina/Data_paper/LiverM4/Mets4MS_merged_R1.fastq.gz",
  "/home/khandler/NAS/Kristina/Data_paper/LiverM4/Mets4MS_merged_R2.fastq.gz")

##Extract only columns with barcode information
bar.table.full <- bar.table.start[, 1:288] 

##Visualize UMI counts in UMAP space for each individual MULTI-seq barcode and extract only positive ones (indicated by clusters in color red)
bar.tsne <- barTSNE(bar.table.start[,1:288]) 
for (i in 3:ncol(bar.tsne)) {
  g <- ggplot(bar.tsne, aes(x = TSNE1, y = TSNE2, color = bar.tsne[,i])) +
    geom_point() +
    scale_color_gradient(low = "black", high = "red") +
    ggtitle(colnames(bar.tsne)[i]) +
    theme(legend.position = "none") 
  print(g)
}

good.bars <- paste("Bar",c(2:4,7:9,11,12,14:16,18:22,25:27,30,31,35,36,40,42,46:60,62:77,79,80,82:84,
                           86,87,89,91,93:95,97:100,102:106,108,111,113:116,118,119,121:124,126,127,
                           129,131,135,136,138,140:142,144,147:149,151:156,160,161,164,166,168,170,
                           171,174:176,181:185,190:203,205,207,209:219,221,222,225:229,233,234,236:242,
                           245:247,250,253,256,259,260,262:265,267,268,271,272,274,276,279,280,283:288),sep="") 
bar.table <- bar.table.start[, good.bars]  

##MULTIseq BC Classification 
#1st round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- names(round.calls)[which(round.calls == "Negative")]
#Remove negative cells 
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#2nd round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#3rd round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#4th round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#5th round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls) #no negative cells left 
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])

##Combine final calls and merge with Cell IDs in bar.table.full 
final.calls <- c(round.calls, rep("Negative", length(neg.cells)))
names(final.calls) <- c(names(round.calls), neg.cells)
#remove duplicates 
final.calls <- final.calls[!duplicated(names(final.calls))]

#make final table with matched Cell ID and MULTI-seq BC information
final_round_classification <- as.data.frame(final.calls, row.names = names(final.calls))

classified_table <- merge(bar.table.full, final_round_classification, by = "row.names")
rownames(classified_table) <- rownames(bar.table.full)
classified_table$Barcode <- classified_table$final.calls

##save table 
saveRDS(classified_table, file = "./data_files_generated/classified_table_Liver_sample3.Rda")

###Sample 4
##UMI count table generation 
bar.table.start <- MULTIseq_UMI_count_matrix_generation_BD(
  "./data/MULTI-seq_3Prime_3Plates_barcodes.xlsx",
  "./data/Mets4M1kept_barcodes.txt",
  "/home/khandler/NAS/Kristina/Data_paper/Liver4M1/Mets4_M1_MS_S2_R1_001.fastq.gz",
  "/home/khandler/NAS/Kristina/Data_paper/Liver4M1/Mets4_M1_MS_S2_R2_001.fastq.gz")

##Extract only columns with barcode information 
bar.table.full <- bar.table.start[, 1:288] 

##Visualize UMI counts in UMAP space for each individual MULTI-seq barcode and extract only positive ones (indicated by clusters in color red)
bar.tsne <- barTSNE(bar.table.start[,1:288]) 
for (i in 3:ncol(bar.tsne)) {
  g <- ggplot(bar.tsne, aes(x = TSNE1, y = TSNE2, color = bar.tsne[,i])) +
    geom_point() +
    scale_color_gradient(low = "black", high = "red") +
    ggtitle(colnames(bar.tsne)[i]) +
    theme(legend.position = "none") 
  print(g)
}

good.bars <- paste("Bar",c(2,6,12:16,18:21,23,26,30,34:36,38,41,44,46,49,53:56,59:62,66,
                           68:72,74,76,77,79,81,83,88,89,92,94,98:104,106,107,108,110:115,117,
                           120,121,124,125,126,128:133,136:143,145,146,148:151,153,154:161,
                           163,165,166,168,169,171:189,195:198,202:210,212:222,225:227,231,
                           239,242,244,245,248:250,254,259,260,262,268:271,273,274,282,284,285),sep="") 
bar.table <- bar.table.start[, good.bars]  

##MULTIseq BC Classification 
#1st round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- names(round.calls)[which(round.calls == "Negative")]
#Remove negative cells 
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#2nd round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#3rd round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#4th round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls) #no negative cells left 
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])

##Combine final calls and merge with Cell IDs in bar.table.full 
final.calls <- c(round.calls, rep("Negative", length(neg.cells)))
names(final.calls) <- c(names(round.calls), neg.cells)
#remove duplicates 
final.calls <- final.calls[!duplicated(names(final.calls))]

#make final table with matched Cell ID and MULTI-seq BC information
final_round_classification <- as.data.frame(final.calls, row.names = names(final.calls))
classified_table <- merge(bar.table.full, final_round_classification, by = "row.names")
rownames(classified_table) <- rownames(bar.table.full)
classified_table$Barcode <- classified_table$final.calls

##save table  
saveRDS(classified_table, file = "./data_files_generated/classified_table_Liver_sample4.Rda")

###Sample 5
##UMI count table generation 
bar.table.start <- MULTIseq_UMI_count_matrix_generation_BD(
  "./data/MULTI-seq_3Prime_3Plates_barcodes.xlsx",
  "./data/Mets4M2kept_barcodes.txt",
  "/home/khandler/NAS/Kristina/Data_paper/Liver4M2/Mets4_M2_MS_S4_R1_001.fastq.gz",
  "/home/khandler/NAS/Kristina/Data_paper/Liver4M2/Mets4_M2_MS_S4_R2_001.fastq.gz")

##Extract only columns with barcode information 
bar.table.full <- bar.table.start[, 1:288] 

##Visualize UMI counts in UMAP space for each individual MULTI-seq barcode and extract only positive ones (indicated by clusters in color red)
bar.tsne <- barTSNE(bar.table.start[,1:288]) 
for (i in 3:ncol(bar.tsne)) {
  g <- ggplot(bar.tsne, aes(x = TSNE1, y = TSNE2, color = bar.tsne[,i])) +
    geom_point() +
    scale_color_gradient(low = "black", high = "red") +
    ggtitle(colnames(bar.tsne)[i]) +
    theme(legend.position = "none") 
  print(g)
}

good.bars <- paste("Bar",c(8,10,13:15,20,22,23,26:28,34,39,40,45,48,50,52,56,60,61,63,
                           68,69,75,83,91,96,101,104,105,112,115,117,122,124,126:128,133,
                           137,138,139,141,143,145,148,152,157:159,166,168,171,173,175,
                           179,182,186,188,190:209,211:215,218:222,223,225,227,228,230,
                           231,233:258,260:265,267:272,274:288),sep="") 
bar.table <- bar.table.start[, good.bars]  


##MULTIseq BC Classification 
#1st round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- names(round.calls)[which(round.calls == "Negative")]
#Remove negative cells 
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#2nd round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#3rd round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#4th round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#5th round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#6th round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls) #no negative cells left 
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])

##Combine final calls and merge with Cell IDs in bar.table.full 
final.calls <- c(round.calls, rep("Negative", length(neg.cells)))
names(final.calls) <- c(names(round.calls), neg.cells)
#remove duplicates 
final.calls <- final.calls[!duplicated(names(final.calls))]

#make final table with matched Cell ID and MULTI-seq BC information
final_round_classification <- as.data.frame(final.calls, row.names = names(final.calls))
classified_table <- merge(bar.table.full, final_round_classification, by = "row.names")
rownames(classified_table) <- rownames(bar.table.full)
classified_table$Barcode <- classified_table$final.calls

##save table 
saveRDS(classified_table, file = "./data_files_generated/classified_table_Liver_sample5.Rda")

###Sample 6
##UMI count table generation 
bar.table.start <- MULTIseq_UMI_count_matrix_generation_BD(
  "./data/MULTI-seq_3Prime_3Plates_barcodes.xlsx",
  "./data/MetsM5kept_barcodes.txt",
  "/home/khandler/NAS/Kristina/Data_paper/LiverM5/Mets5_MS_S6_R1_001.fastq.gz",
  "/home/khandler/NAS/Kristina/Data_paper/LiverM5/Mets5_MS_S6_R2_001.fastq.gz")

##Extract only columns with barcode information 
bar.table.full <- bar.table.start[, 1:288] 

##Visualize UMI counts in UMAP space for each individual MULTI-seq barcode and extract only positive ones (indicated by clusters in color red)
bar.tsne <- barTSNE(bar.table.start[,1:288]) 
for (i in 3:ncol(bar.tsne)) {
  g <- ggplot(bar.tsne, aes(x = TSNE1, y = TSNE2, color = bar.tsne[,i])) +
    geom_point() +
    scale_color_gradient(low = "black", high = "red") +
    ggtitle(colnames(bar.tsne)[i]) +
    theme(legend.position = "none") 
  print(g)
}

good.bars <- paste("Bar",c(2,3,6,8:12,14:17,19,23,26,27,30:35,37,38,40,42,43,47,50:52,54,55,
                           60:65,70:79,81,82,84,87:90,92,93,97:99,101:103,107:116,120,122:124,
                           127,128,133,139,140,143:146,151,154,155,162,170,172,177,181,182,186,
                           188,190,192:194,197:199,204,206,207,211,212,214,215,217,220,227,230,
                           232,233,235,236,240,242,243,248,249,253,254,256,258,263,265,270,273,276,283,284,286,287),sep="") 
bar.table <- bar.table.start[, good.bars]  

##MULTIseq BC Classification 
#1st round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- names(round.calls)[which(round.calls == "Negative")]
#Remove negative cells 
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#2nd round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#3rd round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#4th round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls) #no negative cells left 
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])

##Combine final calls and merge with Cell IDs in bar.table.full 
final.calls <- c(round.calls, rep("Negative", length(neg.cells)))
names(final.calls) <- c(names(round.calls), neg.cells)
#remove duplicates 
final.calls <- final.calls[!duplicated(names(final.calls))]

#make final table with matched Cell ID and MULTI-seq BC information
final_round_classification <- as.data.frame(final.calls, row.names = names(final.calls))
classified_table <- merge(bar.table.full, final_round_classification, by = "row.names")
rownames(classified_table) <- rownames(bar.table.full)
classified_table$Barcode <- classified_table$final.calls

##save table 
saveRDS(classified_table, file = "./data_files_generated/classified_table_Liver_sample6.Rda")

###Sample 7
##UMI count table generation 
bar.table.start <- MULTIseq_UMI_count_matrix_generation_BD(
  "./data/MULTI-seq_3Prime_3Plates_barcodes.xlsx",
  "./data/Mets6M1kept_barcodes.txt",
  "/home/khandler/NAS/Kristina/Data_paper/Liver6M1/Mets6_M1_MS_S2_L004_R1_001.fastq.gz",
  "/home/khandler/NAS/Kristina/Data_paper/Liver6M1/Mets6_M1_MS_S2_L004_R2_001.fastq.gz")

##Extract only columns with barcode information 
bar.table.full <- bar.table.start[, 1:288] 

##Visualize UMI counts in UMAP space for each individual MULTI-seq barcode and extract only positive ones (indicated by clusters in color red)
bar.tsne <- barTSNE(bar.table.start[,1:288]) 
for (i in 3:ncol(bar.tsne)) {
  g <- ggplot(bar.tsne, aes(x = TSNE1, y = TSNE2, color = bar.tsne[,i])) +
    geom_point() +
    scale_color_gradient(low = "black", high = "red") +
    ggtitle(colnames(bar.tsne)[i]) +
    theme(legend.position = "none") 
  print(g)
}

good.bars <- paste("Bar",c(6,8,9,11,15:17,20,21,24,27,32,43,77,78,90,98:102,107,108,110:113,
                           115,120,121,125:128,130,132:135,139,140,142,145,148,149,151,154,
                           155,157,160,161,163,164,167:171,173,175,178,179,181:186,188,189,
                           191,193,194,197,199,201:204,207:215,217:219,224,230,232:235,237:239,
                           242,248,251,252,257,273,282),sep="") 
bar.table <- bar.table.start[, good.bars]  

##MULTIseq BC Classification 
#1st round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- names(round.calls)[which(round.calls == "Negative")]
#Remove negative cells 
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#2nd round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#3rd round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#4th round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls) #no negative cells left 
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])

##Combine final calls and merge with Cell IDs in bar.table.full 
final.calls <- c(round.calls, rep("Negative", length(neg.cells)))
names(final.calls) <- c(names(round.calls), neg.cells)
#remove duplicates 
final.calls <- final.calls[!duplicated(names(final.calls))]

#make final table with matched Cell ID and MULTI-seq BC information
final_round_classification <- as.data.frame(final.calls, row.names = names(final.calls))
classified_table <- merge(bar.table.full, final_round_classification, by = "row.names")
rownames(classified_table) <- rownames(bar.table.full)
classified_table$Barcode <- classified_table$final.calls

##save table 
saveRDS(classified_table, file = "./data_files_generated/classified_table_Liver_sample7.Rda")

###Sample 8
##UMI count table generation 
bar.table.start <- MULTIseq_UMI_count_matrix_generation_BD(
  "./data/MULTI-seq_3Prime_3Plates_barcodes.xlsx",
  "./data/Mets6M2kept_barcodes.txt",
  "/home/khandler/NAS/Kristina/Data_paper/Liver6M2/Mets6_M2_MS_S4_L004_R1_001.fastq.gz",
  "/home/khandler/NAS/Kristina/Data_paper/Liver6M2/Mets6_M2_MS_S4_L004_R2_001.fastq.gz")

##Extract only columns with barcode information 
bar.table.full <- bar.table.start[, 1:288] 

##Visualize UMI counts in UMAP space for each individual MULTI-seq barcode and extract only positive ones (indicated by clusters in color red)
bar.tsne <- barTSNE(bar.table.start[,1:288]) 
for (i in 3:ncol(bar.tsne)) {
  g <- ggplot(bar.tsne, aes(x = TSNE1, y = TSNE2, color = bar.tsne[,i])) +
    geom_point() +
    scale_color_gradient(low = "black", high = "red") +
    ggtitle(colnames(bar.tsne)[i]) +
    theme(legend.position = "none") 
  print(g)
}

good.bars <- paste("Bar",c(2,10,11,13,15,16,19,23:25,32,39,46,48,49,51:53,54,58,63,65,66,68,81,
                           84,86:88,92,97,100:102,105:106,111,116,117,119,124,127:129,131:132,142,143,145,148,150,
                           158,159,165,166,168:171,173,178,179,181,184,186,188,189,190,193:196,199,
                           201:203,205,206,210:212,215,216,218,220:224,227,228,231,235,240,246,251:252,255,
                           256,259,262,263,266,268,271,276,277,280,282:284,286,288),sep="") 

bar.table <- bar.table.start[, good.bars]  

##MULTIseq BC Classification 
#1st round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- names(round.calls)[which(round.calls == "Negative")]
#Remove negative cells 
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#2nd round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#3rd round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#4th round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])

##Combine final calls and merge with Cell IDs in bar.table.full 
final.calls <- c(round.calls, rep("Negative", length(neg.cells)))
names(final.calls) <- c(names(round.calls), neg.cells)
#remove duplicates 
final.calls <- final.calls[!duplicated(names(final.calls))]

#make final table with matched Cell ID and MULTI-seq BC information
final_round_classification <- as.data.frame(final.calls, row.names = names(final.calls))

classified_table <- merge(bar.table.full, final_round_classification, by = "row.names")
rownames(classified_table) <- rownames(bar.table.full)
classified_table$Barcode <- classified_table$final.calls

##save table 
saveRDS(classified_table, file = "./data_files_generated/classified_table_Liver_sample8.Rda")

###Sample 9
##UMI count table generation 
bar.table.start <- MULTIseq_UMI_count_matrix_generation_BD(
  "./data/MULTI-seq_3Prime_3Plates_barcodes.xlsx",
  "./data/Mets6M3kept_barcodes.txt",
  "/home/khandler/NAS/Kristina/Data_paper/Liver6M3/Mets6_M3_MS_S6_L004_R1_001.fastq.gz",
  "/home/khandler/NAS/Kristina/Data_paper/Liver6M3/Mets6_M3_MS_S6_L004_R2_001.fastq.gz")

##Extract only columns with barcode information 
bar.table.full <- bar.table.start[, 1:288] 

##Visualize UMI counts in UMAP space for each individual MULTI-seq barcode and extract only positive ones (indicated by clusters in color red)
bar.tsne <- barTSNE(bar.table.start[,1:288]) 
for (i in 3:ncol(bar.tsne)) {
  g <- ggplot(bar.tsne, aes(x = TSNE1, y = TSNE2, color = bar.tsne[,i])) +
    geom_point() +
    scale_color_gradient(low = "black", high = "red") +
    ggtitle(colnames(bar.tsne)[i]) +
    theme(legend.position = "none") 
  print(g)
}

good.bars <- paste("Bar",c(2:4,8,10:24,26:28,30,32:52,54,56,57,59:63,65:75,77:82,84,87:97,99,101:112,
                           114:121,123:127,131,132,134:141,143:145,147:151,153:159,161,163:170,
                           172:174,176:180,182,184:191,193:207,209:220,222:234,236:245,247,248,
                           252,253,255:262,264,265,268:270,272:284,286,287),sep="") 
bar.table <- bar.table.start[, good.bars]  

##MULTIseq BC Classification 
#1st round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- names(round.calls)[which(round.calls == "Negative")]
#Remove negative cells 
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#2nd round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#3rd round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#4th round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#5th round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])

##Combine final calls and merge with Cell IDs in bar.table.full 
final.calls <- c(round.calls, rep("Negative", length(neg.cells)))
names(final.calls) <- c(names(round.calls), neg.cells)
#remove duplicates 
final.calls <- final.calls[!duplicated(names(final.calls))]

#make final table with matched Cell ID and MULTI-seq BC information
final_round_classification <- as.data.frame(final.calls, row.names = names(final.calls))

classified_table <- merge(bar.table.full, final_round_classification, by = "row.names")
rownames(classified_table) <- rownames(bar.table.full)
classified_table$Barcode <- classified_table$final.calls

##save table 
saveRDS(classified_table, file = "./data_files_generated/classified_table_Liver_sample9.Rda")

###Sample 10
##UMI count table generation 
bar.table.start <- MULTIseq_UMI_count_matrix_generation_BD(
  "./data/MULTI-seq_3Prime_3Plates_barcodes.xlsx",
  "./data/WTkept_barcodes.txt",
  "/home/khandler/NAS/Kristina/Data_paper/LiverWTCre/WT_MS_S8_L004_R1_001.fastq.gz",
  "/home/khandler/NAS/Kristina/Data_paper/LiverWTCre/WT_MS_S8_L004_R2_001.fastq.gz")

##Extract only columns with barcode information
bar.table.full <- bar.table.start[, 1:288] 

##Visualize UMI counts in UMAP space for each individual MULTI-seq barcode and extract only positive ones (indicated by clusters in color red)
bar.tsne <- barTSNE(bar.table.start[,1:288]) 
for (i in 3:ncol(bar.tsne)) {
  g <- ggplot(bar.tsne, aes(x = TSNE1, y = TSNE2, color = bar.tsne[,i])) +
    geom_point() +
    scale_color_gradient(low = "black", high = "red") +
    ggtitle(colnames(bar.tsne)[i]) +
    theme(legend.position = "none") 
  print(g)
}

good.bars <- paste("Bar",c(2,3,14,19:22,24:28,30:32,34,36,40,42,43,45,47,48,50,51,
                           53,62:64,66:69,73,76:78,88,92,93,95,97,101,102,105:111,
                           113:118,120:127,132:138,141:146,154,155,161,165,166,168,170,
                           174,175,179,181:183,188,195:220,222:225,227:244,246:251,253,254,
                           257:259,262,263,265:269,271,274:288),sep="") 
bar.table <- bar.table.start[, good.bars]  

##MULTIseq BC Classification 
#1st round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- names(round.calls)[which(round.calls == "Negative")]
#Remove negative cells 
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#2nd round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#3rd round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#4th round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#5th round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])

##Combine final calls and merge with Cell IDs in bar.table.full 
final.calls <- c(round.calls, rep("Negative", length(neg.cells)))
names(final.calls) <- c(names(round.calls), neg.cells)
#remove duplicates 
final.calls <- final.calls[!duplicated(names(final.calls))]

#make final table with matched Cell ID and MULTI-seq BC information
final_round_classification <- as.data.frame(final.calls, row.names = names(final.calls))

classified_table <- merge(bar.table.full, final_round_classification, by = "row.names")
rownames(classified_table) <- rownames(bar.table.full)
classified_table$Barcode <- classified_table$final.calls

##save table
saveRDS(classified_table, file = "./data_files_generated/classified_table_Liver_sample10.Rda")

########## CRC organoids mixing species sample ##########
##UMI count table generation 
bar.table.start <- MULTIseq_UMI_count_matrix_generation_BD(
  "./data/MULTI-seq_3Prime_3Plates_barcodes.xlsx",
  "./data/OrgMixingkept_barcodes.txt",
  "/home/khandler/NAS/Kristina/Data_paper/Org_mix/Org_mix_MS_S2_R1_001.fastq.gz",
  "/home/khandler/NAS/Kristina/Data_paper/Org_mix/Org_mix_MS_S2_R2_001.fastq.gz")

##Extract only columns with barcode information 
bar.table.full <- bar.table.start[, 1:288] 

##Visualize UMI counts in UMAP space for each individual MULTI-seq barcode and extract only positive ones (indicated by clusters in color red)
bar.tsne <- barTSNE(bar.table.start[,1:288]) 
for (i in 3:ncol(bar.tsne)) {
  g <- ggplot(bar.tsne, aes(x = TSNE1, y = TSNE2, color = bar.tsne[,i])) +
    geom_point() +
    scale_color_gradient(low = "black", high = "red") +
    ggtitle(colnames(bar.tsne)[i]) +
    theme(legend.position = "none") 
  print(g)
}

good.bars <- paste("Bar",c(2,7,9,14,15,17:20,23,24,26,28:37,42,43,47,49,51:53,56:63,65,68,70,
                           74:77,80,82,83,87,90,95:97,99,100,102,105,106,111,114,117:120,124,
                           127,130:132,134,136,138,140:143,148,150:153,155,156,159,160,165:167,
                           170,174,177,179,180,182,184:188,191,194,196,198,201,205:209,211,212,
                           214,217,221,222,224,225,232,236,241,243:245,247,251,254,255,260:264,
                           267:270,273,274,276,281,282,283,286),sep="") 

bar.table <- bar.table.start[, good.bars]  

##MULTIseq BC Classification 
#1st round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- names(round.calls)[which(round.calls == "Negative")]
#Remove negative cells 
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#2nd round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#3rd round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#4th round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])

##Combine final calls and merge with Cell IDs in bar.table.full 
final.calls <- c(round.calls, rep("Negative", length(neg.cells)))
names(final.calls) <- c(names(round.calls), neg.cells)
#remove duplicates 
final.calls <- final.calls[!duplicated(names(final.calls))]

#make final table with matched Cell ID and MULTI-seq BC information
final_round_classification <- as.data.frame(final.calls, row.names = names(final.calls))
classified_table <- merge(bar.table.full, final_round_classification, by = "row.names")
rownames(classified_table) <- rownames(bar.table.full)
classified_table$Barcode <- classified_table$final.calls

##save table 
saveRDS(classified_table, file = "./data_files_generated/classified_table_Organoid_mixing.Rda")

########## Spleen samples ##########
###Sample 1
##UMI count table generation 
bar.table.startS1 <- MULTIseq_bar_table_generation10X(
  "./data/MULTI-seq_McGinnis_96bc.xlsx",
  "./data/Spleen_sample1_filtered_feature_bc_matrix",
  "/home/khandler/NAS/Kristina/Data_paper/Spleen5/Spleen25_MS_S1_L001_R1_001.fastq.gz",
  "/home/khandler/NAS/Kristina/Data_paper/Spleen5/Spleen25_MS_S1_L001_R2_001.fastq.gz")

##Extract only columns with barcode information 
bar.table.full <- bar.table.startS1[, 1:96] 

##Visualize UMI counts in UMAP space for each individual MULTI-seq barcode and extract only positive ones (indicated by clusters in color red)
bar.tsne <- barTSNE(bar.table.start[,1:96]) 
for (i in 3:ncol(bar.tsne)) {
  g <- ggplot(bar.tsne, aes(x = TSNE1, y = TSNE2, color = bar.tsne[,i])) +
    geom_point() +
    scale_color_gradient(low = "black", high = "red") +
    ggtitle(colnames(bar.tsne)[i]) +
    theme(legend.position = "none") 
  print(g)
}

good.bars <- paste("Bar",c(3,4,6:9,13,14,16,19:21,23,24,29,34:36,38,41,51,54:57,59,60,
                           66,68,70,72,78,79,80,82,83,85,86,87,89,93:96),sep="") 
bar.table <- bar.table.startS1[, good.bars]  

##MULTIseq BC Classification 
#1st round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- names(round.calls)[which(round.calls == "Negative")]
#Remove negative cells 
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#2nd round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#3rd round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])

##Combine final calls and merge with Cell IDs in bar.table.full 
final.calls <- c(round.calls, rep("Negative", length(neg.cells)))
names(final.calls) <- c(names(round.calls), neg.cells)
#remove duplicates 
final.calls <- final.calls[!duplicated(names(final.calls))]

#make final table with matched Cell ID and MULTI-seq BC information
final_round_classification <- as.data.frame(final.calls, row.names = names(final.calls))
classified_table <- merge(bar.table.full, final_round_classification, by = "row.names")
rownames(classified_table) <- rownames(bar.table.full)
classified_table$Barcode <- classified_table$final.calls

##save table
saveRDS(classified_table, file = "./data_files_generated/classified_table_Spleen_sample1.Rda")

###Sample 2
##UMI count table generation 
bar.table.start <- MULTIseq_bar_table_generation10X(
  "./data/MULTI-seq_McGinnis_96bc.xlsx",
  "./data/Spleen_sample2_filtered_feature_bc_matrix",
  "/home/khandler/NAS/Kristina/Data_paper/Spleen6/Spleen50_L1_MS_S2_L001_R1_001.fastq.gz",
  "/home/khandler/NAS/Kristina/Data_paper/Spleen6/Spleen50_L1_MS_S2_L001_R2_001.fastq.gz")

##Extract only columns with barcode information 
bar.table.full <- bar.table.start[, 1:96] 

##Visualize UMI counts in UMAP space for each individual MULTI-seq barcode and extract only positive ones (indicated by clusters in color red)
bar.tsne <- barTSNE(bar.table.start[,1:96]) 
for (i in 3:ncol(bar.tsne)) {
  g <- ggplot(bar.tsne, aes(x = TSNE1, y = TSNE2, color = bar.tsne[,i])) +
    geom_point() +
    scale_color_gradient(low = "black", high = "red") +
    ggtitle(colnames(bar.tsne)[i]) +
    theme(legend.position = "none") 
  print(g)
}

good.bars <- paste("Bar",c(3,4:8,11,14:17,21:24,26:30,34:36,38,41:43,45,46,48,50,53,60,70,72,86,87,95,96),sep="") 
bar.table <- bar.table.start[, good.bars]  

##MULTIseq BC Classification 
#1st round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- names(round.calls)[which(round.calls == "Negative")]
#Remove negative cells 
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#2nd round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#3rd round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#4th round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])

##Combine final calls and merge with Cell IDs in bar.table.full 
final.calls <- c(round.calls, rep("Negative", length(neg.cells)))
names(final.calls) <- c(names(round.calls), neg.cells)
#remove duplicates 
final.calls <- final.calls[!duplicated(names(final.calls))]

#make final table with matched Cell ID and MULTI-seq BC information
final_round_classification <- as.data.frame(final.calls, row.names = names(final.calls))
classified_table <- merge(bar.table.full, final_round_classification, by = "row.names")
rownames(classified_table) <- rownames(bar.table.full)
classified_table$Barcode <- classified_table$final.calls

##save table 
saveRDS(classified_table, file = "./data_files_generated/classified_table_Spleen_sample2.Rda")

########## Crohn's biopsy samples ##########
###P387
##UMI count table generation 
bar.table.start <- MULTIseq_UMI_count_matrix_generation_BD(
  "./data/MULTI-seq_3Prime_3Plates_barcodes.xlsx",
  "./data/CrohnP387kept_barcodes.txt",
  "/home/khandler/NAS/Kristina/Data_paper/CrohnP387/CrohnP387MS_S2_R1_001.fastq.gz",
  "/home/khandler/NAS/Kristina/Data_paper/CrohnP387/CrohnP387MS_S2_R2_001.fastq.gz")

##Extract only columns with barcode information 
bar.table.full <- bar.table.start[, 1:288] 

##Visualize UMI counts in UMAP space for each individual MULTI-seq barcode and extract only positive ones (indicated by clusters in color red)
bar.tsne <- barTSNE(bar.table.start[,1:288]) 
for (i in 3:ncol(bar.tsne)) {
  g <- ggplot(bar.tsne, aes(x = TSNE1, y = TSNE2, color = bar.tsne[,i])) +
    geom_point() +
    scale_color_gradient(low = "black", high = "red") +
    ggtitle(colnames(bar.tsne)[i]) +
    theme(legend.position = "none") 
  print(g)
}

good.bars <- paste("Bar",c(12,14,28,35,43,55,63,66,69,74,78,81,90,96,
                           100,103,109,111,114,116,129,131,138,149,160,178,
                           211,222,224,227,229,236,242,258,266,267,274,277,279,
                           280,283,284,285,286,287),sep="") 
bar.table <- bar.table.start[, good.bars]  

##MULTIseq BC Classification 
#1st round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- names(round.calls)[which(round.calls == "Negative")]
#Remove negative cells 
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#2nd round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#3rd round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])

##Combine final calls and merge with Cell IDs in bar.table.full 
final.calls <- c(round.calls, rep("Negative", length(neg.cells)))
names(final.calls) <- c(names(round.calls), neg.cells)
#remove duplicates 
final.calls <- final.calls[!duplicated(names(final.calls))]

#make final table with matched Cell ID and MULTI-seq BC information
final_round_classification <- as.data.frame(final.calls, row.names = names(final.calls))
classified_table <- merge(bar.table.full, final_round_classification, by = "row.names")
rownames(classified_table) <- rownames(bar.table.full)
classified_table$Barcode <- classified_table$final.calls

##save table 
saveRDS(classified_table, file = "./data_files_generated/classified_table_Crohn_P387.Rda")

###P393
##UMI count table generation 
bar.table.start <- MULTIseq_UMI_count_matrix_generation_BD(
  "./data/MULTI-seq_3Prime_3Plates_barcodes.xlsx",
  "./data/CrohnP393kept_barcodes.txt",
  "/home/khandler/NAS/Kristina/Data_paper/CrohnP393/Crohn_P393_MS_S2_R1_001.fastq.gz",
  "/home/khandler/NAS/Kristina/Data_paper/CrohnP393/Crohn_P393_MS_S2_R2_001.fastq.gz")

##Extract only columns with barcode information 
bar.table.full <- bar.table.start[, 1:288] 

##Visualize UMI counts in UMAP space for each individual MULTI-seq barcode and extract only positive ones (indicated by clusters in color red)
bar.tsne <- barTSNE(bar.table.start[,1:288]) 
for (i in 3:ncol(bar.tsne)) {
  g <- ggplot(bar.tsne, aes(x = TSNE1, y = TSNE2, color = bar.tsne[,i])) +
    geom_point() +
    scale_color_gradient(low = "black", high = "red") +
    ggtitle(colnames(bar.tsne)[i]) +
    theme(legend.position = "none") 
  print(g)
}

good.bars <- paste("Bar",c(5,6,8,9,18,24,25,28,38,40,58,64,74,87,93,98,99:103,106,107,
                           109,110,115,116,118:120,122,124,125,128,131,132,134,136,
                           137,144:146,152,157,160,165,167,168,171,172,174,175,176,
                           181,182,185,189,199:201,212,218,223,239,241,253,271),sep="") 
bar.table <- bar.table.start[, good.bars]  

##MULTIseq BC Classification 
#1st round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- names(round.calls)[which(round.calls == "Negative")]
#Remove negative cells 
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#2nd round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#3rd round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#4th round 
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
threshold.results <- findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "top") +
  geom_vline(xintercept=threshold.results$extrema, lty=2) + scale_color_manual(values=c("red","black","blue")) + 
  xlab("Inter-maxima Quantile")
round.calls <- classifyCells(bar.table, q=findQ(threshold.results$res, threshold.results$extrema))
table(round.calls)
neg.cells <- c(neg.cells, names(round.calls)[which(round.calls == "Negative")])

##Combine final calls and merge with Cell IDs in bar.table.full 
final.calls <- c(round.calls, rep("Negative", length(neg.cells)))
names(final.calls) <- c(names(round.calls), neg.cells)
#remove duplicates 
final.calls <- final.calls[!duplicated(names(final.calls))]

#make final table with matched Cell ID and MULTI-seq BC information
final_round_classification <- as.data.frame(final.calls, row.names = names(final.calls))
classified_table <- merge(bar.table.full, final_round_classification, by = "row.names")
rownames(classified_table) <- rownames(bar.table.full)
classified_table$Barcode <- classified_table$final.calls

##save table 
saveRDS(classified_table, file = "./data_files_generated/classified_table_Crohn_P393.Rda")

