########## Part 4.3: Veins cell type co-localisation  ##########
#This part does cell type co-localisation analysis comparing CV and PV areas 
#Code by Karsten Bach 

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Sphere-sequencing/Highly_multiplexed_FISH_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")
source("./functions_and_packages/5.Functions_cell_co_localisation.R")

###Load data 
merged <- readRDS(file = "./data_files_generated/Resolve_seurat_anno.rds")

#only consider CV and PV areas 
Idents(merged) <- "vein"
veins <- subset(merged, idents = c("CV","PV"))
srt <- veins

########## Co-localisation analysis #########
#Extract meta data as df
pD <- data.frame("X"=srt$X,
                 "Y"=srt$Y,
                 "sampleID"=srt$sampleID,
                 "Area"=srt$vein, 
                 "Slide"=srt$Slide,
                 "annotation"=srt$annotation,
                 "Mets"=srt$Mets,
                 "Barcode"=srt$Cell
)
#We don't want to distinguish hepatocytes
pD$annotation <- gsub("_CV|_PV","",pD$annotation) 

#Define ROI
pD$ROI <- paste0(pD$Slide,"_",pD$Area)

###Parameter settings
#Use a knn function to look for the K nearest neighbors in a pre-defined radius (distance), 
#the k should thus be sufficiently large to include most neighbors, 
#that's why I set it to 41 (+1 as the cell itself is always the closest neighbor)
k <- 41 
slide.scale <- 0.138 # Scale
contactdist <- 10 # Distance to look for neighbors in um
contactdist.pix <- floor(contactdist/slide.scale) # pixel distance used in the functions
n.it <- 1000 # Number of iterations (ideally this should be going up to 1k but then it takes quite a while)
BPPARAM <- MulticoreParam(workers = 4, tasks = 20)# Parallelizing this, if your computer suffers reduce workers

###Here you set the Areas
#Defining the two Areas (result is area1-area2)
area1 <- "PV"
area2 <- "CV"

cons.list <- lapply(unique(pD$ROI), function(ROI) {
  cons <- getContacts(input=pD[pD$ROI==ROI,],rad=contactdist.pix,knn=k,expected=unique(pD$annotation))
  cons$ROI <- ROI
  cons$Slide <- unique(pD[pD$ROI==ROI,"Slide"])
  cons$Area <- unique(pD[pD$ROI==ROI,"Area"])
  return(cons)
})

#This object contains number of edges between each cell type for all slides and areas
cons <- do.call(rbind,cons.list)

###Now this is reduced to only the PV area and contains a column with the difference with regards to the CV area
cons.diff <- getDiff(input=cons,diffFrom=area1,diffTo=area2)

###Now we test for significance across all slides by running the permutation test
res.list <- lapply(unique(pD$Slide), function(SLIDE) {
  print(SLIDE)
  res <- runTest(input=pD[pD$Slide==SLIDE,], con=cons.diff[cons.diff$Slide==SLIDE,],
                 k=k,rad=contactdist.pix,n_it=n.it,BPPARAM=BPPARAM,expected=unique(pD$annotation),
                 diffFrom=area1,diffTo=area2)
  res$Slide <- SLIDE
  return(res)})

res <- do.call(rbind,res.list)

###We can visualize this now by summing across the slides the ones that were significant at alpha = 0.01 & more than 30% difference
res$FC <- (res$Edges + res$EdgeDiff_To_CV) / res$Edges
res$FC[is.nan(res$FC)] <- 0
res.sum <- group_by(res, From, To) %>%
  summarize(Score = sum((Pval <= 0.01 & abs(FC) > 1.3) * sign(EdgeDiff_To_CV))/n())

###Triangle plotting
pmat <- xtabs(Score ~ From + To,data=res.sum)
cellCellContactMap(pmat)
cellCellContactMap(pmat,order=levels(factor(res.sum$From)))
ggsave("./figures/4.3/Contacts_PV-CV.pdf",width = 12, height = 10)
ggsave("./figures/4.3/Contacts_PV-CV.svg",width = 12, height = 10)
