########## Part 6: Fragment size and cell counts bias ##########
#This part investigates potential biases introduced by different cutoffs of fragment-size and cell counts per fragment 

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Fragment-sequencing/Fragment-seq_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")
source("./functions_and_packages/1.Packages.R")
source("./functions_and_packages/5.Functions_lobule_layers.R")
source("./functions_and_packages/4.Functions_cells_per_fragment_cutoff.R")
source("./functions_and_packages/3.Functions_fragment_size_GFP_integration.R")
source("./functions_and_packages/7.Functions_ligand_receptor_analysis.R")

###load R object
liverSpS5C <- readRDS("./data_files_generated/LiverMerged_afterBC_anno_BS_5cells_zC_lobules.Rda")
metastasis <- readRDS("./data_files_generated/LiverMerged_afterBC_anno_BS_5cells_zC_lobules_mets_distance.Rda")

#only use injected samples for comparison analysis 
Idents(liverSpS5C) <- "orig.ident"
Injected <- subset(liverSpS5C, idents = c("M1","M3","M4","4M1","4M2","M5","6M1","6M2","6M3"))

########## Investigation of bias introduced by fragment size #########
##### Size distributions between conditions 
###distal cs. proximal 
Idents(metastasis) <- "Mets_distance"
subS1 <- subset(metastasis, idents = c("distal"))
subS2 <- subset(metastasis, idents = c("proximal"))

##Make data frames of fragment sizes per fragment and sample
dfS1_size <- BS_df_for_boxplot_per_sample(subS1,"fragment","fragment_size","distal")
dfS2_size <- BS_df_for_boxplot_per_sample(subS2,"fragment","fragment_size","proximal")

#merge data frames
size_merged_df <- rbind(dfS1_size,dfS2_size)

##Plot in violinplot 
p <- ggplot(size_merged_df,aes(x = sample,y = fragment_size, fill = sample)) +theme_classic() +
  geom_violin() +
  geom_jitter(position = position_jitter(seed = 1, width =0.4),size = 1) + 
  theme(axis.text = element_text(size = 30))  +
  ggtitle("Fragment size per fragment - Mets cond") + xlab("Sample") + 
  ylab("Size (µm)") + theme(axis.title= element_text(size = 25)) + 
  theme(plot.title = element_text(size = 25, face = "bold")) + 
  ggsignif::geom_signif(comparisons = list(c("distal", "proximal")), textsize=7,test = "wilcox.test",map_signif_level = c("***"=0.001,"**"=0.01,"*"=0.05)) + 
  theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30)) 
p + ggsave("./figures/6/Size_violinplot_Mets.pdf",width = 15, height = 10)
p + ggsave("./figures/6/Size_violinplot_Mets.svg",width = 15, height = 10)

###CV vs. PV 
Idents(Injected) <- "vein"
subS1 <- subset(Injected, idents = c("CV"))
subS2 <- subset(Injected, idents = c("PV"))

##Make data frames of fragment sizes per sample
dfS1_size <- BS_df_for_boxplot_per_sample(subS1,"fragment","fragment_size","CV")
dfS2_size <- BS_df_for_boxplot_per_sample(subS2,"fragment","fragment_size","PV")

#merge data frames
size_merged_df <- rbind(dfS1_size,dfS2_size)

##Plot in violinplot 
p <- ggplot(size_merged_df,aes(x = sample,y = fragment_size, fill = sample)) +theme_classic() +
  geom_violin() +
  geom_jitter(position = position_jitter(seed = 1, width =0.4),size = 1) + 
  theme(axis.text = element_text(size = 30))  +
  ggtitle("Fragment size per fragment - vein") + xlab("Sample") + 
  ylab("Size (µm)") + theme(axis.title= element_text(size = 25)) + 
  theme(plot.title = element_text(size = 25, face = "bold")) + 
  ggsignif::geom_signif(comparisons = list(c("CV", "PV")),textsize=7,test = "wilcox.test",map_signif_level = c("***"=0.001,"**"=0.01,"*"=0.05)) +
  theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30)) 
p + ggsave("./figures/6/Size_violinnplot_vein.pdf",width = 15, height = 10)
p + ggsave("./figures/6/Size_violinplot_vein.svg",width = 15, height = 10)

#####Comparison of large and small fragments 
#add small and large identity based on fragment-size, split fragment sizes in half 
a <- Injected@meta.data
a <- a[,c(4,14)]
#remove duplicated
a <- a[!duplicated(a$fragment), ]

small <- 211:325
fragments_small <- a[a$fragment_size %in% small, ]

large <- 326:457
fragments_large <- a[a$fragment_size %in% large, ]

#add identity in meta data 
Injected$Size_cond <- NA
Injected@meta.data <- Injected@meta.data %>%
  mutate(Size_cond = case_when(
    fragment %in% fragments_small$fragment ~ "small",
    fragment %in% fragments_large$fragment ~ "large",
    TRUE ~ NA_character_))

##colours 
#B cells (orange): B_mem: #EF975B; B_plasma: #F4741E    
#Granulocytes (yellow): Basophils:#56595B; Neutrophils: #AEAEAF  
#Endothelial (grey): LVECs: #F4CB1C; LSECs:#F4E740
#Stromal cells (brown): Fibroblasts: #8E5229; Stellate cells: #C1A08A
#Hepatocytes (redish): #F46042
#DCs (turquoise): cDC1: #91C6C4; cDC2: #315B5A; pDC: #7AEDD9  
#KC (purple): Kupffer: #E20FE8; Kupffer_Endo: #C491C6  
#Monocytes (green): M_C1q: #19E80F; M_patrolling: #B5EDB2; Mac_Ly6c: #566B44  
#T (blue): CD4: #2323E5 ; CD8: #9A9AE5; NKT: #2C91E5; NK:#95A0C6   
#Metastatic cells (red): #F90606
#Cholangiocytes (pink): #EDB2D4

p <- DimPlot(Injected, split.by = "Size_cond", label = FALSE,label.size = 5, group.by = "annotation", pt.size = 0.5, 
             cols = c("#EF975B","#F4741E","#56595B","#91C6C4",
                      "#315B5A","#EDB2D4",
                      "#8E5229","#F46042","#E20FE8",
                      "#C491C6","#F4E740", "#F4CB1C", "#19E80F",
                      "#566B44","#F90606","#B5EDB2","#AEAEAF",
                      "#95A0C6","#2C91E5",  "#7AEDD9","#C1A08A", "#2323E5","#9A9AE5")) + 
  theme(legend.title = element_text(size = 22), legend.text = element_text(size = 22)) + 
  theme(title = element_text(size = 25))+ theme(axis.text = element_text(size = 30)) +
  ggtitle("Annotation") + NoLegend()
p + ggsave("./figures/6/Annotated_umap_split_size_cond.pdf", width = 15, height = 10)
p + ggsave("./figures/6/Annotated_umap_split_size_cond.svg", width = 15, height = 10)

#####DEG between LECs of fragments of different sizes
##prepare data
#subset large and small and then merge again to remove NA 
Idents(Injected) <- "Size_cond"
smallS <- subset(Injected, idents = "small")
largeS <- subset(Injected, idents = "large")

#merge
seurt <- merge(smallS, largeS)
#subset LECs only 
Idents(seurt) <- "annotation.broad"
seurt <- subset(seurt, idents = "LECs")

#convert Seurat object to SCE object with required meta data information 
m <- GetAssayData(seurt, assay = "RNA", slot = "counts")
pD <- data.frame("barcode"=colnames(m),
                 "Fragment"=seurt$fragment,
                 "Layer"=seurt$lobule_layer,
                 "Size_cond"=seurt$Size_cond,
                 "Batch"=seurt$orig.ident)

sce <- SingleCellExperiment(assays=list(counts=m),colData=DataFrame(pD))

#Sum across the fragments
sumd <- aggregateAcrossCells(sce,ids=sce$Fragment)
#Only consider fragments with at least 5 LECs  
sumd <- sumd[,sumd$ncells >= 5]

###Set up edgeR object 
y <- DGEList(counts=counts(sumd),samples=data.frame(colData(sumd)))

#Remove lowly expressed genes
keep <- filterByExpr(y, group=y$samples$Size_cond)
y <- y[keep,]

#Normalization
y <- calcNormFactors(y)

##Defining the model matrix
y$samples$Size_cond <- factor(y$samples$Size_cond)
y$samples$Layer <- factor(y$samples$Layer)
y$samples$Layer <- ordered(y$samples$Layer)

#account for variability in samples with '~ Batch + ...' and for lobule layers with '~ Layer'
mdl <- model.matrix(~Batch + Size_cond + Layer,y$samples)

##Follow the standard edgeR workflow
y <- estimateDisp(y, mdl)
fit <- glmQLFit(y, mdl, robust=TRUE)

res <- glmQLFTest(fit, coef = "Size_condsmall")

##plot DEGs
top <- topTags(res,n=nrow(y))$table
top$Gene <- rownames(top)

top <- top %>% 
  mutate(
    Expression = case_when(logFC >= 0.5 & FDR <= 0.05 ~ "High in small",
                           logFC <= -0.5 & FDR <= 0.05 ~ "High in large",
                           TRUE ~ "Non sig.")
  )

#label the top genes 
top$genelabels <- ifelse(top$Gene == "Zbtb20"
                         |top$Gene == "Slfn5"
                         |top$Gene == "mt-Nd2"
                         |top$Gene == "Nptn"
                         |top$Gene == "mt-Nd6"
                         |top$Gene == "Degs1",
                         TRUE,FALSE)

p <- ggplot(top, aes(x=logFC, y=-log10(FDR))) +
  geom_point(aes(color = Expression),size=5) +
  geom_text_repel(aes(label=ifelse(top$genelabels, top$Gene,"")),size=8) +
  xlab("logFC") + 
  ylab("-log10(FDR)") + ggtitle("LECs small vs. large fragment size") + 
  scale_color_manual(values = c( "gray50"),guide = "none") + theme_classic() + 
  theme(axis.title= element_text(size = 25)) + theme(axis.text = element_text(size = 30))  + 
  theme(plot.title = element_text(size = 25, face = "bold"))  + 
  guides(colour = guide_legend(override.aes = list(size=1.5))) + 
  theme(legend.title = element_text(size = 15), legend.text = element_text(size = 15)) +  
  geom_vline(xintercept = c(-0.5,0.5), linetype = "dashed", color = "gray50") + 
  coord_cartesian(ylim = c(0,2)) + 
  #plot dashed line at 0.05 p-value 
  geom_hline(yintercept = c(1.30103), linetype = "dashed", color = "gray50") 
p + ggsave("./figures/6/Vulcano_small_large_LECs.pdf",width = 12, height = 10)
p + ggsave("./figures/6/Vulcano_small_large_LECs.svg",width = 12, height = 10)  

#####DEG between KCs of fragments of different sizes 
##prepare data 
seurt <- merge(smallS, largeS)
Idents(seurt) <- "annotation.broad"
seurt <- subset(seurt, idents = "Kupffer")

#convert Seurat object to SCE object with required meta data information 
m <- GetAssayData(seurt, assay = "RNA", slot = "counts")
pD <- data.frame("barcode"=colnames(m),
                 "Fragment"=seurt$fragment,
                 "Layer"=seurt$lobule_layer,
                 "Size_cond"=seurt$Size_cond,
                 "Batch"=seurt$orig.ident)

sce <- SingleCellExperiment(assays=list(counts=m),colData=DataFrame(pD))

#Sum across the fragments
sumd <- aggregateAcrossCells(sce,ids=sce$Fragment)
#Only consider fragments with at least 5 cells 
sumd <- sumd[,sumd$ncells >= 5]

##Set up edgeR object 
y <- DGEList(counts=counts(sumd),samples=data.frame(colData(sumd)))

#Remove lowly expressed genes
keep <- filterByExpr(y, group=y$samples$Size_cond)
y <- y[keep,]

#Normalization
y <- calcNormFactors(y)

##Defining the model matrix
y$samples$Size_cond <- factor(y$samples$Size_cond)
y$samples$Layer <- factor(y$samples$Layer)
y$samples$Layer <- ordered(y$samples$Layer)

#account for variability in samples with '~ Batch + ...' and for lobule layers with '~ Layer'
mdl <- model.matrix(~Batch + Size_cond + Layer,y$samples)

##Follow the standard edgeR workflow
y <- estimateDisp(y, mdl)
fit <- glmQLFit(y, mdl, robust=TRUE)

res <- glmQLFTest(fit, coef = "Size_condsmall")

##plot DEGs
top <- topTags(res,n=nrow(y))$table
top$Gene <- rownames(top)

top <- top %>% 
  mutate(
    Expression = case_when(logFC >= 0.5 & FDR <= 0.05 ~ "High in small",
                           logFC <= -0.5 & FDR <= 0.05 ~ "High in large",
                           TRUE ~ "Non sig.")
  )

#label the top genes 
top$genelabels <- ifelse(top$Gene == "Rab10"
                         |top$Gene == "Tmsb4x"
                         |top$Gene == "Tram1",
                         TRUE,FALSE)

p <- ggplot(top, aes(x=logFC, y=-log10(FDR))) +
  geom_point(aes(color = Expression),size=5) +
  geom_text_repel(aes(label=ifelse(top$genelabels, top$Gene,"")),size=8) +
  xlab("logFC") + 
  ylab("-log10(FDR)") + ggtitle("KCs small vs. large fragment sizes") + 
  scale_color_manual(values = c( "gray50"),guide = "none") + theme_classic() + 
  theme(axis.title= element_text(size = 25)) + theme(axis.text = element_text(size = 30))  + 
  theme(plot.title = element_text(size = 25, face = "bold"))  + 
  guides(colour = guide_legend(override.aes = list(size=1.5))) + 
  theme(legend.title = element_text(size = 15), legend.text = element_text(size = 15)) +  
  geom_vline(xintercept = c(-0.5,0.5), linetype = "dashed", color = "gray50") + 
  coord_cartesian(ylim = c(0,2)) + 
  #plot dashed line at 0.05 p-value 
  geom_hline(yintercept = c(1.30103), linetype = "dashed", color = "gray50") 
p + ggsave("./figures/6/Vulcano_small_large_KCs.pdf",width = 12, height = 10)
p + ggsave("./figures/6/Vulcano_small_large_KCs.svg",width = 12, height = 10)  

#####plot zonated DEGs of LECs between different fragment sizes in a combined vulcano plot 
#run DEG of small and large
DE_zonated_genes(smallS,"annotation.broad","LECs","./figures/6/","small")
DE_zonated_genes(largeS,"annotation.broad","LECs","./figures/6/","large")

#read in results 
topSmall <- read.csv("./figures/6/small_LECs_LM_zonated_genes.csv")
topLarge <- read.csv("./figures/6/large_LECs_LM_zonated_genes.csv")

#apply significant cutoff of FDR and logFC to highlight significant ones in vulcano plot  
topSmall <- topSmall %>% 
  mutate(
    Expression = case_when(logFC >= 0.5 & FDR <= 0.05 ~ "High in PV small",
                           logFC <= -0.5 & FDR <= 0.05 ~ "High in CV small",
                           TRUE ~ "Non sig.")
  )

topLarge <- topLarge %>% 
  mutate(
    Expression = case_when(logFC >= 0.5 & FDR <= 0.05 ~ "High in PV large",
                           logFC <= -0.5 & FDR <= 0.05 ~ "High in CV large",
                           TRUE ~ "Non sig.")
  )

#extract the genes we describe in more detail 
topSmall2 <- topSmall[topSmall$X %in% c("Galnt15","Mecom","Cd34","Ifi207",
                                        "Cd36","Lhx6","Plpp1","Jup","Slco2a1","Car8",
                                        "Acer2"),]
topLarge2 <- topLarge[topLarge$X %in% c("Galnt15","Mecom","Cd34","Ifi207",
                                        "Cd36","Lhx6","Plpp1","Jup","Slco2a1","Car8",
                                        "Acer2"),]

top <- rbind(topSmall2, topLarge2)

p <- ggplot(top, aes(x=logFC, y=-log10(FDR))) +
  geom_point(aes(color = Expression),size=5) +
  geom_text(data=top[top$FDR<1 & abs(top$logFC) > 0,], aes(label=Gene),size=3) +
  xlab("logFC") + 
  ylab("-log10(FDR)") + ggtitle("LECs small vs. large zonated genes") + 
  scale_color_manual(values = c("darkgoldenrod2","blue","red","purple"),guide = "none") + theme_classic() + 
  theme(axis.title= element_text(size = 25)) + theme(axis.text = element_text(size = 30))  + 
  theme(plot.title = element_text(size = 25, face = "bold"))  + 
  guides(colour = guide_legend(override.aes = list(size=1.5))) + 
  theme(legend.title = element_text(size = 15), legend.text = element_text(size = 15)) +  
  geom_vline(xintercept = c(-0.5,0.5), linetype = "dashed", color = "gray50") 
p + ggsave("./figures/6/Vulcano_sig_zonated_genes_size_comp.pdf",width = 12, height = 10)
p + ggsave("./figures/6/Vulcano_sig_Zonated_genes_size_comp.svg",width = 12, height = 10)  

#####L-R interaction analysis 
#split dataset 
Idents(Injected) <- "Size_cond"
small_LR <- subset(Injected, ident = "small")
large_LR <- subset(Injected, ident = "large")

###split veins for each sample and produce CellPhoneDB input files 
#load human and mouse ensemble symbols
human <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
human <- useMart(host='apr2019.archive.ensembl.org', dataset="hsapiens_gene_ensembl", biomart='ENSEMBL_MART_ENSEMBL')

mouse <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
mouse <- useMart(host='apr2019.archive.ensembl.org', dataset="mmusculus_gene_ensembl", biomart='ENSEMBL_MART_ENSEMBL')

###small
#subset veins 
Idents(small_LR) <- "vein"
cv <- subset(small_LR, idents = "CV")
pv <- subset(small_LR, idents = "PV")

#input file generation 
Input_files_CellPhoneDB_generation(cv,'annotation.broad',"cv","./figures/6/CellPhoneDB/small/") 
Input_files_CellPhoneDB_generation(pv,'annotation.broad',"pv","./figures/6/CellPhoneDB/small/") 

###large
#subset veins 
Idents(large_LR) <- "vein"
cv <- subset(large_LR, idents = "CV")
pv <- subset(large_LR, idents = "PV")

#input file generation 
Input_files_CellPhoneDB_generation(cv,'annotation.broad',"cv","./figures/6/CellPhoneDB/large/") 
Input_files_CellPhoneDB_generation(pv,'annotation.broad',"pv","./figures/6/CellPhoneDB/large/") 

###run CellPhoneDB in 6.1.CellPHoneDB_veins_bias_fragment_size_counts.sh 

###PLotting Kupffer and T in different fragment sizes 
##small
file1_mean <- "./figures/6/CellPhoneDB/small/cv/significant_means.txt"
file1_pval <- "./figures/6/CellPhoneDB/small/cv/pvalues.txt"

file2_mean <- "./figures/6/CellPhoneDB/small/pv/significant_means.txt"
file2_pval <- "./figures/6/CellPhoneDB/small/pv/pvalues.txt"

L_R_CellPhoneDB_comp_2samples("Kupffer.T",file1_mean,file1_pval,file2_mean,file2_pval,"./figures/6/CellPhoneDB/","small")

##large 
file1_mean <- "./figures/6/CellPhoneDB/large/cv/significant_means.txt"
file1_pval <- "./figures/6/CellPhoneDB/large/cv/pvalues.txt"

file2_mean <- "./figures/6/CellPhoneDB/large/pv/significant_means.txt"
file2_pval <- "./figures/6/CellPhoneDB/large/pv/pvalues.txt"

L_R_CellPhoneDB_comp_2samples("Kupffer.T",file1_mean,file1_pval,file2_mean,file2_pval,"./figures/6/CellPhoneDB/","large")

###read in combined L-R, inculde original analysis with all fragments 
df_vein_small <- read.csv("./figures/6/CellPhoneDB/interactions_comp_two_cond_Kupffer.T_small.csv")
df_vein_large <- read.csv("./figures/6/CellPhoneDB/interactions_comp_two_cond_Kupffer.T_large.csv")
df_vein_all <- read.csv("./figures/2.5.3/interactions_comp_two_cond_Kupffer.T_vein.csv")

###extract top 10 based on all analysis 
#calculate difference between distal and proximal to order interacting pairs in decreasing order
#make values positive 
df_vein_all$diff <- df_vein_all$sample1_mean - df_vein_all$sample2_mean
df_vein_all$diff <- abs(df_vein_all$diff)

#order dataframe decreasing based on the difference of interaction scores 
df_vein_all <- df_vein_all[order(df_vein_all$diff,decreasing = TRUE),]

##split the plot so you can put cv and pv scores in one column for interaction_score and one for pvals
df_vein_all_cv <- df_vein_all[,c("interacting_pair", "sample1_mean","sample1_pval","diff")]
df_vein_all_cv$vein <- "all"
colnames(df_vein_all_cv) <- c("interacting_pair","interaction_score","pvalue","diff", "condition")
#make minus to put on one axis, has to be changed to + scores later in Illustrator 
df_vein_all_cv$interaction_score <- -(df_vein_all_cv$interaction_score)
#take only the top 15 interactions for plotting
df_vein_all_cv <- df_vein_all_cv[1:10,]

df_vein_all_pv <- df_vein_all[,c("interacting_pair", "sample2_mean","sample2_pval","diff")]
df_vein_all_pv$vein <- "all"
colnames(df_vein_all_pv) <- c("interacting_pair","interaction_score","pvalue", "diff", "condition")
#take only the top 15 interactions for plotting
df_vein_all_pv <- df_vein_all_pv[1:10,]

df_vein_all1 <- rbind(df_vein_all_cv,df_vein_all_pv)

df_vein_all1 <- df_vein_all1[df_vein_all1$pval < 0.05,] 

##these are the interactions for plotting 
interactions <- df_vein_all1$interacting_pair

#extract from small 
df_vein_small_cv <- df_vein_small[,c("interacting_pair", "sample1_mean","sample1_pval")]
df_vein_small_cv$vein <- "small"
colnames(df_vein_small_cv) <- c("interacting_pair","interaction_score","pvalue", "condition")
df_vein_small_cv$interaction_score <- -(df_vein_small_cv$interaction_score)

df_vein_small_pv <- df_vein_small[,c("interacting_pair", "sample2_mean","sample2_pval")]
df_vein_small_pv$vein <- "small"
colnames(df_vein_small_pv) <- c("interacting_pair","interaction_score","pvalue", "condition")

df_vein_small1 <- rbind(df_vein_small_cv, df_vein_small_pv)
df_vein_small2 <- df_vein_small1[df_vein_small1$interacting_pair %in% interactions,]

#extract from large 
df_vein_large_cv <- df_vein_large[,c("interacting_pair", "sample1_mean","sample1_pval")]
df_vein_large_cv$vein <- "large"
colnames(df_vein_large_cv) <- c("interacting_pair","interaction_score","pvalue", "condition")
df_vein_large_cv$interaction_score <- -(df_vein_large_cv$interaction_score)

df_vein_large_pv <- df_vein_large[,c("interacting_pair", "sample2_mean","sample2_pval")]
df_vein_large_pv$vein <- "large"
colnames(df_vein_large_pv) <- c("interacting_pair","interaction_score","pvalue", "condition")

df_vein_large1 <- rbind(df_vein_large_cv, df_vein_large_pv)
df_vein_large2 <- df_vein_large1[df_vein_large1$interacting_pair %in% interactions,]

#merge
df_vein_all1$diff <- NULL
df1 <- rbind(df_vein_all1, df_vein_small2)
df1 <- rbind(df1, df_vein_large2)

#make barplot 
p <- ggplot(df1, aes(x=interaction_score, y=reorder(interacting_pair,+interaction_score),fill = pvalue, groups = condition)) + theme_classic() +
  geom_bar(stat = "identity",position = "dodge",aes(colour = condition), linewidth = 0.5) + 
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 6, name ="RdBu")))(256), values = c(0,0.001,0.003,0.005,0.01,0.05), 
                       breaks = c(0,0.001,0.003,0.005,0.01,0.05)) +
  theme(axis.text = element_text(size = 15))  + theme(axis.text.x = element_text(angle = 90)) +
  scale_color_manual(values = c("purple","red","yellow", "purple","red","yellow")) +
  theme(axis.title= element_text(size = 25)) + 
  ggtitle("L-R interactions - CellPhoneDB (Kupffer|T) small vs. large") + xlab("Interaction score") + 
  ylab("L-R interaction")  +
  theme(plot.title = element_text(size = 15, face = "bold")) + theme(legend.text = element_text(size = 22),
                                                                     legend.title= element_text(size = 22)) + 
  guides(fill=guide_legend(title="P-value")) 
p + ggsave("./figures/6/Interaction_KC_T_vein_all_small_large.pdf", width = 12, height = 5)
p + ggsave("./figures/6/Interaction_KC_T_vein_all_small_large.svg", width = 12, height = 5)

########## Investigation of bias introduced by cell counts per fragment cutoffs #########
#####different cell number cutoffs, compare > 5 cells/fragment with > 20 cells/fragment
Fragment_cell_cutoff("/mnt/khandler/R_projects/Fragment-sequencing/Fragment-seq_analysis/data_files_generated/LiverMerged_afterBC_anno_BS_5cells_zC_lobules.Rda",
                   20,"/mnt/khandler/R_projects/Fragment-sequencing/Fragment-seq_analysis/data_files_generated/LiverMerged_afterBC_anno_BS_20cells_zC_lobules.Rda")

liverSpS20C <- readRDS("/mnt/khandler/R_projects/Fragment-sequencing/Fragment-seq_analysis/data_files_generated/LiverMerged_afterBC_anno_BS_20cells_zC_lobules.Rda")
liverSpS5C <- readRDS("/mnt/khandler/R_projects/Fragment-sequencing/Fragment-seq_analysis/data_files_generated/LiverMerged_afterBC_anno_BS_5cells_zC_lobules.Rda")

#only take injected samples into consideration 
Idents(liverSpS5C) <- "orig.ident"
liverSpS5C <- subset(liverSpS5C, idents = c("M1","M3","M4","4M1","4M2","M5","6M1","6M2","6M3"))

Idents(liverSpS20C) <- "orig.ident"
liverSpS20C <- subset(liverSpS20C, idents = c("M1","M3","M4","4M1","4M2","M5","6M1","6M2","6M3"))

#####plot umap with > 5 cells/fragment and with > 20 cells/fragment
###colors 
#B cells (orange): B_mem: #EF975B; B_plasma: #F4741E    
#Granulocytes (yellow): Basophils:#56595B; Neutrophils: #AEAEAF  
#Endothelial (grey): LVECs: #F4CB1C; LSECs:#F4E740
#Stromal cells (brown): Fibroblasts: #8E5229; Stellate cells: #C1A08A
#Hepatocytes (redish): #F46042
#DCs (turquoise): cDC1: #91C6C4; cDC2: #315B5A; pDC: #7AEDD9  
#KC (purple): Kupffer: #E20FE8; Kupffer_Endo: #C491C6  
#Monocytes (green): M_C1q: #19E80F; M_patrolling: #B5EDB2; Mac_Ly6c: #566B44  
#T (blue): CD4: #2323E5 ; CD8: #9A9AE5; NKT: #2C91E5; NK:#95A0C6   
#Metastatic cells (red): #F90606
#Cholangiocytes (pink): #EDB2D4

p <- DimPlot(liverSpS20C, label = FALSE,label.size = 5, group.by = "annotation", pt.size = 0.5, 
             cols = c("#EF975B","#F4741E","#56595B","#91C6C4",
                      "#315B5A","#EDB2D4",
                      "#8E5229","#F46042","#E20FE8",
                      "#C491C6","#F4E740", "#F4CB1C", "#19E80F",
                      "#566B44","#F90606","#B5EDB2","#AEAEAF",
                      "#95A0C6","#2C91E5",  "#7AEDD9","#C1A08A", "#2323E5","#9A9AE5")) + 
  theme(legend.title = element_text(size = 22), legend.text = element_text(size = 22)) + 
  theme(title = element_text(size = 25))+ theme(axis.text = element_text(size = 30)) +
  ggtitle("Annotation") + NoLegend()
p + ggsave("./figures/6/Annotated_umap_C20.pdf", width = 15, height = 10)
p + ggsave("./figures/6/Annotated_umap_C20.svg", width = 15, height = 10)

p <- DimPlot(liverSpS5C, label = FALSE,label.size = 5, group.by = "annotation", pt.size = 0.5, 
             cols = c("#EF975B","#F4741E","#56595B","#91C6C4",
                      "#315B5A","#EDB2D4",
                      "#8E5229","#F46042","#E20FE8",
                      "#C491C6","#F4E740", "#F4CB1C", "#19E80F",
                      "#566B44","#F90606","#B5EDB2","#AEAEAF",
                      "#95A0C6","#2C91E5",  "#7AEDD9","#C1A08A", "#2323E5","#9A9AE5")) + 
  theme(legend.title = element_text(size = 22), legend.text = element_text(size = 22)) + 
  theme(title = element_text(size = 25))+ theme(axis.text = element_text(size = 30)) +
  ggtitle("Annotation") + NoLegend()
p + ggsave("./figures/6/Annotated_umap_C5.pdf", width = 15, height = 10)
p + ggsave("./figures/6/Annotated_umap_C5.svg", width = 15, height = 10)

#####assess cell counts per fragment between conditions 
### CV vs. PV 
Idents(liverSpS5C) <- "vein"
subS1 <- subset(liverSpS5C, idents = c("CV"))
subS2 <- subset(liverSpS5C, idents = c("PV"))

df1 <- as.data.frame(table(subS1$fragment))
df1$fragment <- df1$Var1
df1$sample <- "CV"
df1$Var1 <- NULL
colnames(df1) <- c("fragment_size","fragment", "sample")

df2 <- as.data.frame(table(subS2$fragment))
df2$fragment <- df2$Var1
df2$sample <- "PV"
df2$Var1 <- NULL
colnames(df2) <- c("fragment_size","fragment", "sample")

#merge data frames
size_merged_df <- rbind(df1,df2)

#Plot in violinplot 
p <- ggplot(size_merged_df,aes(x = sample,y = fragment_size, fill = sample)) +theme_classic() +
  geom_violin() +
  geom_jitter(position = position_jitter(seed = 1, width =0.4),size = 1) + 
  theme(axis.text = element_text(size = 30))  +
  ggtitle("Cell counts per fragment - Vein") + xlab("Sample") + 
  ylab("Cell number per fragment") + theme(axis.title= element_text(size = 25)) + 
  theme(plot.title = element_text(size = 25, face = "bold")) +
  ggsignif::geom_signif(comparisons = list(c("CV", "PV")),textsize=7,test = "wilcox.test",map_signif_level = c("***"=0.001,"**"=0.01,"*"=0.05)) +
  theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30)) 
p + ggsave("./figures/6/Cell_numbers_per_fragment_vein.pdf",width = 15, height = 10)
p + ggsave("./figures/6/Cell_numbers_per_fragment_vein.svg",width = 15, height = 10)

###distal vs. proximal
Idents(metastasis) <- "Mets_distance"
subS1 <- subset(metastasis, idents = c("distal"))
subS2 <- subset(metastasis, idents = c("proximal"))

df1 <- as.data.frame(table(subS1$fragment))
df1$fragment <- df1$Var1
df1$sample <- "distal"
df1$Var1 <- NULL
colnames(df1) <- c("fragment_size","fragment", "sample")

df2 <- as.data.frame(table(subS2$fragment))
df2$fragment <- df2$Var1
df2$sample <- "proximal"
df2$Var1 <- NULL
colnames(df2) <- c("fragment_size","fragment", "sample")

#merge data frames
size_merged_df <- rbind(df1,df2)

#Plot in violinplot 
p <- ggplot(size_merged_df,aes(x = sample,y = fragment_size, fill = sample)) +theme_classic() +
  geom_violin() +
  geom_jitter(position = position_jitter(seed = 1, width =0.4),size = 1) + 
  theme(axis.text = element_text(size = 30))  +
  ggtitle("Cell number per fragment Mets") + xlab("Sample") + 
  ylab("Cell number per fragment") + theme(axis.title= element_text(size = 25)) + 
  theme(plot.title = element_text(size = 25, face = "bold")) +
  ggsignif::geom_signif(comparisons = list(c("distal", "proximal")),textsize=7,test = "wilcox.test",map_signif_level = c("***"=0.001,"**"=0.01,"*"=0.05)) +
  theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30)) 
p + ggsave("./figures/6/Cell_numbers_per_fragment_Mets.pdf",width = 15, height = 10)
p + ggsave("./figures/6/Cell_numbers_per_fragment_Mets.svg",width = 15, height = 10)

#####DEG between LECs of fragments from different cell count conditions 
###prepare data 
#add metadata column for C5 and C20
liverSpS5C@meta.data$cell_numb_cond <- "C5"
liverSpS20C@meta.data$cell_numb_cond <- "C20"

liverSpS5C$fragment <- paste(liverSpS5C$fragment, "_C5", sep="")
liverSpS20C$fragment <- paste(liverSpS20C$fragment, "_C20", sep="")

seurt <- merge(liverSpS5C, liverSpS20C)
Idents(seurt) <- "annotation.broad"
seurt <- subset(seurt, idents = "LECs")

#convert Seurat object to SCE object with required meta data information 
m <- GetAssayData(seurt, assay = "RNA", slot = "counts")
pD <- data.frame("barcode"=colnames(m),
                 "Fragment"=seurt$fragment,
                 "Layer"=seurt$lobule_layer,
                 "Size_cond"=seurt$cell_numb_cond,
                 "Batch"=seurt$orig.ident)

sce <- SingleCellExperiment(assays=list(counts=m),colData=DataFrame(pD))

#Sum across the fragments
sumd <- aggregateAcrossCells(sce,ids=sce$Fragment)
#Only consider fragments with at least 5 cells 
sumd <- sumd[,sumd$ncells >= 5]

###Set up edgeR object 
y <- DGEList(counts=counts(sumd),samples=data.frame(colData(sumd)))

#Remove lowly expressed genes
keep <- filterByExpr(y, group=y$samples$Size_cond)
y <- y[keep,]

#Normalization
y <- calcNormFactors(y)

###Defining the model matrix
y$samples$Size_cond <- factor(y$samples$Size_cond)
y$samples$Layer <- factor(y$samples$Layer)
y$samples$Layer <- ordered(y$samples$Layer)

#account for variability in samples with '~ Batch + ...' and for lobule layers with '~ Layer'
mdl <- model.matrix(~Batch + Size_cond + Layer,y$samples)

###Follow the standard edgeR workflow
y <- estimateDisp(y, mdl)
fit <- glmQLFit(y, mdl, robust=TRUE)

res <- glmQLFTest(fit, coef = "Size_condC5")

###plot DEGs
top <- topTags(res,n=nrow(y))$table
top$Gene <- rownames(top)

top <- top %>% 
  mutate(
    Expression = case_when(logFC >= 0.5 & FDR <= 0.05 ~ "High in C5",
                           logFC <= -0.5 & FDR <= 0.05 ~ "High in C20",
                           TRUE ~ "Non sig.")
  )

p <- ggplot(top, aes(x=logFC, y=-log10(FDR))) +
  geom_point(aes(color = Expression),size=5)  +
  xlab("logFC") + 
  ylab("-log10(FDR)") + ggtitle("C5 vs. C20 LECs") + 
  scale_color_manual(values = c( "gray50"),guide = "none") + theme_classic() + 
  theme(axis.title= element_text(size = 25)) + theme(axis.text = element_text(size = 30))  + 
  theme(plot.title = element_text(size = 25, face = "bold"))  + 
  guides(colour = guide_legend(override.aes = list(size=1.5))) + 
  theme(legend.title = element_text(size = 15), legend.text = element_text(size = 15)) +  
  geom_vline(xintercept = c(-0.5,0.5), linetype = "dashed", color = "gray50") + 
  coord_cartesian(ylim = c(0,2)) + 
  geom_hline(yintercept = c(1.30103), linetype = "dashed", color = "gray50") 
p + ggsave("./figures/6/Vulcano_C5_C20_LECs.pdf",width = 12, height = 10)
p + ggsave("./figures/6/Vulcano_C5_C20_LECs.svg",width = 12, height = 10)  

#####DEG between KCs of fragments from different cell counts conditions 
###prepare data 
liverSpS5C@meta.data$cell_numb_cond <- "C5"
liverSpS20C@meta.data$cell_numb_cond <- "C20"

liverSpS5C$fragment <- paste(liverSpS5C$fragment, "_C5", sep="")
liverSpS20C$fragment <- paste(liverSpS20C$fragment, "_C20", sep="")

seurt <- merge(liverSpS5C, liverSpS20C)
Idents(seurt) <- "annotation.broad"
seurt <- subset(seurt, idents = "Kupffer")

#convert Seurat object to SCE object with required meta data information 
m <- GetAssayData(seurt, assay = "RNA", slot = "counts")
pD <- data.frame("barcode"=colnames(m),
                 "Fragment"=seurt$fragment,
                 "Layer"=seurt$lobule_layer,
                 "Size_cond"=seurt$cell_numb_cond,
                 "Batch"=seurt$orig.ident)

sce <- SingleCellExperiment(assays=list(counts=m),colData=DataFrame(pD))

#Sum across the fragments
sumd <- aggregateAcrossCells(sce,ids=sce$Fragment)
#Only consider fragments with at least 5 cells 
sumd <- sumd[,sumd$ncells >= 5]

###Set up edgeR object 
y <- DGEList(counts=counts(sumd),samples=data.frame(colData(sumd)))

#Remove lowly expressed genes
keep <- filterByExpr(y, group=y$samples$Size_cond)
y <- y[keep,]

#Normalization
y <- calcNormFactors(y)

###Defining the model matrix
y$samples$Size_cond <- factor(y$samples$Size_cond)
y$samples$Layer <- factor(y$samples$Layer)
y$samples$Layer <- ordered(y$samples$Layer)

#account for variability in samples with '~ Batch + ...' and for lobule layers  with '~ Layer'
mdl <- model.matrix(~Batch + Size_cond + Layer,y$samples)

###Follow the standard edgeR workflow
y <- estimateDisp(y, mdl)
fit <- glmQLFit(y, mdl, robust=TRUE)

res <- glmQLFTest(fit, coef = "Size_condC5")

###plot DEGs
top <- topTags(res,n=nrow(y))$table
top$Gene <- rownames(top)

top <- top %>% 
  mutate(
    Expression = case_when(logFC >= 0.5 & FDR <= 0.05 ~ "High in C5",
                           logFC <= -0.5 & FDR <= 0.05 ~ "High in C20",
                           TRUE ~ "Non sig.")
  )

p <- ggplot(top, aes(x=logFC, y=-log10(FDR))) +
  geom_point(aes(color = Expression),size=5)  +
  xlab("logFC") + 
  ylab("-log10(FDR)") + ggtitle("C5 vs. C20 KCs") + 
  scale_color_manual(values = c( "gray50"),guide = "none") + theme_classic() + 
  theme(axis.title= element_text(size = 25)) + theme(axis.text = element_text(size = 30))  + 
  theme(plot.title = element_text(size = 25, face = "bold"))  + 
  guides(colour = guide_legend(override.aes = list(size=1.5))) + 
  theme(legend.title = element_text(size = 15), legend.text = element_text(size = 15)) +  
  geom_vline(xintercept = c(-0.5,0.5), linetype = "dashed", color = "gray50") + 
  coord_cartesian(ylim = c(0,2)) + 
  #dashed line at p value 0.05
  geom_hline(yintercept = c(1.30103), linetype = "dashed", color = "gray50") 
p + ggsave("./figures/6/Vulcano_C5_C20_KCs.pdf",width = 12, height = 10)
p + ggsave("./figures/6/Vulcano_C5_C20_KCs.svg",width = 12, height = 10)  

#####plot zonated DEGs of LECs between different cutoffs (>5 cells/fragment and >20 cells/fragment) in one vulcano
#run DEG of C5 and C20
DE_zonated_genes(liverSpS5C,"annotation.broad","LECs","./figures/6/","C5")
DE_zonated_genes(liverSpS20C,"annotation.broad","LECs","./figures/6/","C20")

#read in results 
topC5 <- read.csv("./figures/6/C5_LECs_LM_zonated_genes.csv")
topC20 <- read.csv("./figures/6/C20_LECs_LM_zonated_genes.csv")

#apply significant cutoff of FDR and logFC to highlight significant ones in vulcano plot  
topC5 <- topC5 %>% 
  mutate(
    Expression = case_when(logFC >= 0.5 & FDR <= 0.05 ~ "High in PV C5",
                           logFC <= -0.5 & FDR <= 0.05 ~ "High in CV C5",
                           TRUE ~ "Non sig.")
  )

topC20 <- topC20 %>% 
  mutate(
    Expression = case_when(logFC >= 0.5 & FDR <= 0.05 ~ "High in PV C20",
                           logFC <= -0.5 & FDR <= 0.05 ~ "High in CV C20",
                           TRUE ~ "Non sig.")
  )

#extract the genes we describe in more detail 
topC52 <- topC5[topC5$X %in% c("Galnt15","Mecom","Cd34","Ifi207",
                                        "Cd36","Lhx6","Plpp1","Jup","Slco2a1","Car8",
                                        "Acer2"),]
topC202 <- topC20[topC20$X %in% c("Galnt15","Mecom","Cd34","Ifi207",
                                        "Cd36","Lhx6","Plpp1","Jup","Slco2a1","Car8",
                                        "Acer2"),]

top <- rbind(topC52, topC202)

p <- ggplot(top, aes(x=logFC, y=-log10(FDR))) +
  geom_point(aes(color = Expression),size=5) +
  geom_text(data=top[top$FDR<1 & abs(top$logFC) > 0,], aes(label=Gene),size=3) +
  xlab("logFC") + 
  ylab("-log10(FDR)") + ggtitle("LECs DEGs zonated C5 vs. C20") + 
  scale_color_manual(values = c("darkgoldenrod2","blue","red","purple"),guide = "none") + theme_classic() + 
  theme(axis.title= element_text(size = 25)) + theme(axis.text = element_text(size = 30))  + 
  theme(plot.title = element_text(size = 25, face = "bold"))  + 
  guides(colour = guide_legend(override.aes = list(size=1.5))) + 
  theme(legend.title = element_text(size = 15), legend.text = element_text(size = 15)) +  
  geom_vline(xintercept = c(-0.5,0.5), linetype = "dashed", color = "gray50") 
p + ggsave("./figures/6/Vulcano_sig_genes_C5_C20_LEC_comp.pdf",width = 12, height = 10)
p + ggsave("./figures/6/Vulcano_sig_genes_C5_C20_LEC_comp.svg",width = 12, height = 10)  

#####L-R interaction analysis
liverSpS20C <- readRDS("/mnt/khandler/R_projects/Fragment-sequencing/Fragment-seq_analysis/data_files_generated/LiverMerged_afterBC_anno_BS_20cells_zC_lobules.Rda")

#only take injected samples into consideration 
Idents(liverSpS20C) <- "orig.ident"
liverSpS20C <- subset(liverSpS20C, idents = c("M1","M3","M4","4M1","4M2","M5","6M1","6M2","6M3"))

##C20
#subset veins 
Idents(liverSpS20C) <- "vein"
cv <- subset(liverSpS20C, idents = "CV")
pv <- subset(liverSpS20C, idents = "PV")

#input file generation 
Input_files_CellPhoneDB_generation(cv,'annotation.broad',"cv","./figures/6/CellPhoneDB/C20/") 
Input_files_CellPhoneDB_generation(pv,'annotation.broad',"pv","./figures/6/CellPhoneDB/C20/") 

###run CellPhoneDB in 6.1.CellPHoneDB_veins_bias_fragment_size_counts

###PLotting of Kupffer and T in different cell count cutoffs 
###C20
file1_mean <- "./figures/6/CellPhoneDB/C20/cv/significant_means.txt"
file1_pval <- "./figures/6/CellPhoneDB/C20/cv/pvalues.txt"

file2_mean <- "./figures/6/CellPhoneDB/C20/pv/significant_means.txt"
file2_pval <- "./figures/6/CellPhoneDB/C20/pv/pvalues.txt"

L_R_CellPhoneDB_comp_2samples("Kupffer.T",file1_mean,file1_pval,file2_mean,file2_pval,"./figures/6/CellPhoneDB/","C20")

###read in combined L-R 
df_vein_C20 <- read.csv("./figures/6/CellPhoneDB/interactions_comp_two_cond_Kupffer.T_C20.csv")
df_vein_C5 <- read.csv("./figures/2.5.3/interactions_comp_two_cond_Kupffer.T_vein.csv")

###extract top 10 based on all analysis 
#calculate difference between distal and proximal to order interacting pairs in decreasing order
#make values positive 
df_vein_C5$diff <- df_vein_C5$sample1_mean - df_vein_C5$sample2_mean
df_vein_C5$diff <- abs(df_vein_C5$diff)

#order dataframe decreasing based on the difference of interaction scores 
df_vein_C5 <- df_vein_C5[order(df_vein_C5$diff,decreasing = TRUE),]

##split the plot so you can put cv and pv scores in one column for interaction_score and one for pvals
df_vein_C5_cv <- df_vein_C5[,c("interacting_pair", "sample1_mean","sample1_pval","diff")]
df_vein_C5_cv$vein <- "C5"
colnames(df_vein_C5_cv) <- c("interacting_pair","interaction_score","pvalue","diff", "condition")
#make minus to put on one axis, has to be changed to + scores later in Illustrator 
df_vein_C5_cv$interaction_score <- -(df_vein_C5_cv$interaction_score)
#take only the top 15 interactions for plotting
df_vein_C5_cv <- df_vein_C5_cv[1:10,]

df_vein_C5_pv <- df_vein_C5[,c("interacting_pair", "sample2_mean","sample2_pval","diff")]
df_vein_C5_pv$vein <- "C5"
colnames(df_vein_C5_pv) <- c("interacting_pair","interaction_score","pvalue", "diff", "condition")
#take only the top 15 interactions for plotting
df_vein_C5_pv <- df_vein_C5_pv[1:10,]

df_vein_C51 <- rbind(df_vein_C5_cv,df_vein_C5_pv)

df_vein_C51 <- df_vein_C51[df_vein_C51$pval < 0.05,] #only 13 with p-values smaller than 0.05

##these are the interactions for plotting 
interactions <- df_vein_C51$interacting_pair

#extract from C20
df_vein_C20_cv <- df_vein_C20[,c("interacting_pair", "sample1_mean","sample1_pval")]
df_vein_C20_cv$vein <- "C20"
colnames(df_vein_C20_cv) <- c("interacting_pair","interaction_score","pvalue", "condition")
df_vein_C20_cv$interaction_score <- -(df_vein_C20_cv$interaction_score)

df_vein_C20_pv <- df_vein_C20[,c("interacting_pair", "sample2_mean","sample2_pval")]
df_vein_C20_pv$vein <- "C20"
colnames(df_vein_C20_pv) <- c("interacting_pair","interaction_score","pvalue", "condition")

df_vein_C201 <- rbind(df_vein_C20_cv, df_vein_C20_pv)
df_vein_C202 <- df_vein_C201[df_vein_C201$interacting_pair %in% interactions,]

#merge
df_vein_C51$diff <- NULL
df1 <- rbind(df_vein_C51, df_vein_C202)

#make barplot 
p <- ggplot(df1, aes(x=interaction_score, y=reorder(interacting_pair,+interaction_score),fill = pvalue, groups = condition)) + theme_classic() +
  geom_bar(stat = "identity",position = "dodge",aes(colour = condition), linewidth = 0.5) + 
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 6, name ="RdBu")))(256), values = c(0,0.001,0.003,0.005,0.01,0.05), 
                       breaks = c(0,0.001,0.003,0.005,0.01,0.05)) +
  theme(axis.text = element_text(size = 15))  + theme(axis.text.x = element_text(angle = 90)) +
  scale_color_manual(values = c("black","red", "black","red")) +
  theme(axis.title= element_text(size = 25)) + 
  ggtitle("L-R interactions - CellPhoneDB (Kupffer|T)") + xlab("Interaction score") + 
  ylab("L-R interaction")  +
  theme(plot.title = element_text(size = 15, face = "bold")) + theme(legend.text = element_text(size = 22),
                                                                     legend.title= element_text(size = 22)) + 
  guides(fill=guide_legend(title="P-value")) 
p + ggsave("./figures/6/Interaction_KC_T_vein_C5_C20.pdf", width = 12, height = 5)
p + ggsave("./figures/6/Interaction_KC_T_vein_C5_C20.svg", width = 12, height = 5)






