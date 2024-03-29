########## Part 5.2 DGE analysis between Monocytes of different Mets distances ##########
#This part validates findings of DGE analysis between macrophage subtypes in Monocytes of highly multiplexed FISH 

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Fragment-sequencing/Highly_multiplexed_FISH_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")
source("./functions_and_packages/3.Functions_DGE.R")

###Load data 
merged <- readRDS(file = "./data_files_generated/Resolve_seurat_anno.rds")

#subset only samples with visible metastasis
Idents(merged) <- "samples"
mets_samples <- subset(merged, idents = "mets")

########## DGE analysis ##########
DGE_between_Mets_distance_hmFISH(mets_samples, "annotation","Monocytes","Mets_distance","./figures/5.2/") 

top <- read.csv("./figures/5.2/Mets_distance_Monocytes_DGE_hmFISH.csv")

top <- top %>% 
  mutate(
    Expression = case_when(logFC >=0.4 & FDR <= 0.05 ~ "High in proximal",
                           logFC <= -0.4 & FDR <= 0.05 ~ "High in distal",
                           TRUE ~ "Non sig.")
  )

#highlight specific genes 
top$genelabels <- ifelse(top$Gene == "C1qc"
                         |top$Gene == "C1qb"
                         |top$Gene == "Dab2"
                         |top$Gene == "Cd5l"
                         |top$Gene == "Il6ra"
                         |top$Gene == "Tgfbi"
                         |top$Gene == "Tnfrsf1b",
                         TRUE,FALSE)

p <- ggplot(top, aes(x=logFC, y=-log10(FDR))) +
  geom_point(aes(color = Expression),size=3) +
  geom_text_repel(aes(label=ifelse(top$genelabels, top$Gene,"")),size=8) +
  xlab("logFC") + 
  ylab("-log10(FDR)") + ggtitle("linear model Monocytes (FDR 0.05, logFC 0.5)- hmFISH") + 
  scale_color_manual(values = c( "firebrick3","dodgerblue3", "gray50"),guide = "none") + theme_classic() + 
  theme(axis.title= element_text(size = 25)) + theme(axis.text = element_text(size = 30))  + 
  theme(plot.title = element_text(size = 25, face = "bold"))  + guides(colour = guide_legend(override.aes = list(size=1.5))) + 
  theme(legend.title = element_text(size = 15), legend.text = element_text(size = 15)) +  
  geom_vline(xintercept = c(-0.4,.4), linetype = "dashed", color = "gray50") 
p + ggsave("./figures/5.2/Mono_LM_vulcano_proximal_distal_Resolve.pdf",width = 12, height = 10)
p + ggsave("./figures/5.2/Mono_LM_vulcano_proximal_distal_Resolve.svg",width = 12, height = 10)  


