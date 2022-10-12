#This script just puts everything together
library(SingleCellExperiment)
fls <- list.files("data/resolve_out",full.names=TRUE)
cnts <- fls[grepl("counts",fls)]
cnts_nms <- gsub("data/resolve_out/|_counts.csv","",cnts)
pDs <- fls[grepl("centroids",fls)]
pDs_nms <- gsub("data/resolve_out/|_centroids.csv","",pDs)

m.list <- lapply(cnts,function(FILE) as.matrix(read.csv(FILE,header=TRUE,row.names=1)))
gns <- rownames(m.list[[1]])

#Add slide names to cells to avoid non-unique barcodes
m.list <- lapply(1:length(m.list), function(SLIDE) {
		     tmp <- m.list[[SLIDE]][gns,] ## !! Enforce the same row order, cbind is rowname unaware
		     nm <- cnts_nms[SLIDE]
		     colnames(tmp) <- paste0(nm,"_",colnames(tmp))
		     return(tmp)
})

pD.list <- lapply(pDs,function(FILE) read.csv(FILE,header=TRUE,row.names=1))

pD.list <- lapply(1:length(pD.list), function(SLIDE) {
		     tmp <- pD.list[[SLIDE]]
		     nm <- pDs_nms[SLIDE]
		     tmp$Cell <- paste0(nm,"_",tmp$Cell)
		     tmp$Slide <- nm
		     return(tmp)
})

pD <- Reduce(rbind,pD.list)
m <- Reduce(cbind,m.list)
rownames(pD) <- pD$Cell
sce <- SingleCellExperiment(assays=SimpleList(counts=m),colData=DataFrame(pD))
saveRDS(sce,"data/SCE.RDS")
