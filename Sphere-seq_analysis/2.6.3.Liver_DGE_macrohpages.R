########## Part 2.6.3: DGE analysis of macrophage subsets  ##########
#This part does DGE analysis between C1q+ and Ly6c+ macrophages 

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")

###load R object
metastasis <- readRDS("./data_files_generated/LiverMerged_afterBC_anno_BS_5cells_zC_lobules_mets_distance.Rda")

########## DGE analysis of different macrophages subtypes ##########
Idents(metastasis) <- "annotation"
mac <- subset(metastasis, idents = c("Mac_C1q","Mac_Ly6c"))

###Prepare data
mac$sample_mac_type <- paste0(mac$orig.ident,"_",mac$annotation)
seurt <- mac

#convert Seurat object to SCE object with required meta data information 
m <- GetAssayData(mac, assay = "RNA", slot = "counts")
pD <- data.frame("barcode"=colnames(m),
                 "sample"=seurt$sample_mac_type,
                 "annotation"=seurt$annotation,
                 "Batch"=seurt$orig.ident)

sce <- SingleCellExperiment(assays=list(counts=m),colData=DataFrame(pD))

#Sum across the samples
sumd <- aggregateAcrossCells(sce,ids=sce$sample)
#Only consider samples with at least 5 cells
sumd <- sumd[,sumd$ncells >= 5]

###Set up edgeR object
y <- DGEList(counts=counts(sumd),samples=data.frame(colData(sumd)))

#Remove lowly expressed genes
keep <- filterByExpr(y, group=y$samples$annotation)
y <- y[keep,]

#Normalization
y <- calcNormFactors(y)

###Defining the model matrix
y$samples$annotation <- factor(y$samples$annotation)

mdl <- model.matrix(~Batch + annotation,y$samples)

###Follow the standard edgeR workflow
y <- estimateDisp(y, mdl)
fit <- glmQLFit(y, mdl, robust=TRUE)

res <- glmQLFTest(fit)

###save results for plotting
top <- topTags(res,n=nrow(y))$table
top$Gene <- rownames(top)
write.csv(top,"./figures/2.6.3/MacC1q_vs_MacLy6c_DGE.csv")

###Plot result of DGE analysis in vulcano plot 
#read in DGE analysis results 
top <- read.csv("./figures/2.6.3/MacC1q_vs_MacLy6c_DGE.csv")

top <- top %>% 
  mutate(
    Expression = case_when(logFC >= 0.5 & FDR <= 0.05 ~ "High in Mac_Ly6c_high",
                           logFC <= -0.5 & FDR <= 0.05 ~ "High in Mac_C1q",
                           TRUE ~ "Non sig.")
  )

#highlight specific genes 
top$genelabels <- ifelse(top$Gene == "Spp1"
                         |top$Gene == "C1qbp"
                         |top$Gene == "Dab2"
                         |top$Gene == "Mrc1"
                         |top$Gene == "Folr2"
                         |top$Gene == "Apoe"
                         |top$Gene == "Vcan"
                         |top$Gene == "Ly6c2"
                         |top$Gene == "Chil3"
                         |top$Gene == "Ccr2"
                         |top$Gene == "Thbs1"
                         |top$Gene == "Lyz2"
                         |top$Gene == "Il6ra"
                         |top$Gene == "Tgfbi"
                         |top$Gene == "Cd36"
                         |top$Gene == "Lpl"
                         |top$Gene == "Vegfa"
                         |top$Gene == "Tnfrsf1a"
                         |top$Gene == "Tnfrsf1b"
                         |top$Gene == "Trem2",
                         TRUE,FALSE)

p <- ggplot(top, aes(x=logFC, y=-log10(FDR))) +
  geom_point(aes(color = Expression),size=3) +
  geom_text_repel(aes(label=ifelse(top$genelabels, top$Gene,"")),size=8,max.overlaps = Inf) +
  xlab("logFC") + 
  ylab("-log10(FDR)") + ggtitle("MacC1q vs. MacLy6c (FDR ≤ 0.05, logFC >0.5") + 
  scale_color_manual(values = c("dodgerblue3", "firebrick3", "gray50"),guide = "none") + theme_classic() + 
  theme(axis.title= element_text(size = 25)) + theme(axis.text = element_text(size = 30))  + 
  theme(plot.title = element_text(size = 25, face = "bold"))  + 
  guides(colour = guide_legend(override.aes = list(size=1.5))) + 
  theme(legend.title = element_text(size = 15), legend.text = element_text(size = 15)) +  
  geom_vline(xintercept = c(-0.5,0.5), linetype = "dashed", color = "gray50") 
p + ggsave("./figures/2.6.3/MacC1q_vs_MacLy6c_vulcano.pdf",width = 12, height = 10)
p + ggsave("./figures/2.6.3/MacC1q_vs_MacLy6c_vulcano.svg",width = 12, height = 10)  








