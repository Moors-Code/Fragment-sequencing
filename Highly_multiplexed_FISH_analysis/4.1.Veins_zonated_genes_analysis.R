########## Part 4.1: Veins zonated genes analysis ##########
#This part validates findings of new zonated genes from sphere-seq analysis 

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Sphere-sequencing/Highly_multiplexed_FISH_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")
source("./functions_and_packages/3.Functions_DGE.R")

###Load data 
merged <- readRDS(file = "./data_files_generated/Resolve_seurat_anno.rds")

#only consider CV and PV areas 
Idents(merged) <- "vein"
veins <- subset(merged, idents = c("CV","PV"))

######### LECs ##########
###DGE analysis and vulcano plotting 
DGE_between_veins_hmFISH(veins, "annotation","LECs","vein", "./figures/4.1/") 
top <- read.csv("./figures/4.1/vein_LECs_DGE_hmFISH.csv")

top <- top %>% 
  mutate(
    Expression = case_when(logFC >= 0.5 & FDR <= 0.05 ~ "High in PV",
                           logFC <= -0.5 & FDR <= 0.05 ~ "High in CV",
                           TRUE ~ "Non sig.")
  )

#highlight specific genes 
top$genelabels <- ifelse(top$Gene == "Galnt15"
                         |top$Gene == "Mecom"
                         |top$Gene == "Cd34"
                         |top$Gene == "Ifi207"
                         |top$Gene == "Cd36"
                         |top$Gene == "Lhx6"
                         |top$Gene == "Plpp1"
                         |top$Gene == "Jup"
                         |top$Gene == "Slco2a1"
                         |top$Gene == "Car8"
                         |top$Gene == "Acer2",
                         TRUE,FALSE)


p <- ggplot(top, aes(x=logFC, y=-log10(FDR))) +
  geom_point(aes(color = Expression),size=3) +
  geom_text_repel(aes(label=ifelse(top$genelabels, top$Gene,"")),size=8) +
  xlab("logFC") + 
  ylab("-log10(FDR)") + ggtitle("LECs - Resolve (FDR ≤ 0.05, logFC >0.5") + 
  scale_color_manual(values = c("red", "darkgoldenrod2", "gray50"),guide = "none") + theme_classic() + 
  theme(axis.title= element_text(size = 25)) + theme(axis.text = element_text(size = 30))  + 
  theme(plot.title = element_text(size = 25, face = "bold"))  + 
  guides(colour = guide_legend(override.aes = list(size=1.5))) + 
  theme(legend.title = element_text(size = 15), legend.text = element_text(size = 15))  + 
  geom_vline(xintercept = c(-0.5,0.5), linetype = "dashed", color = "gray50")
p + ggsave("./figures/4.1/LECs_CV_PV_Resolve_vulcano.pdf",width = 12, height = 10)
p + ggsave("./figures/4.1/LECs_CV_PV_Resolve_vulcano.svg",width = 12, height = 10)  

###Boxplot of genes of interest 
#CV zonated
boxplot_LECs_per_spatial_area(veins,"Plpp1","red","./figures/4.1/",".pdf")
boxplot_LECs_per_spatial_area(veins,"Plpp1","red","./figures/4.1/",".svg")

#PV zonated 
boxplot_LECs_per_spatial_area(veins,"Galnt15","darkgoldenrod2","./figures/4.1/",".pdf")
boxplot_LECs_per_spatial_area(veins,"Galnt15","darkgoldenrod2","./figures/4.1/",".svg")


########## Kupffer cells ##########
###DGE analysis and vulcano plotting 
DGE_between_veins_hmFISH(veins, "annotation","Kupffer","vein", "./figures/4.1/") 
top <- read.csv("./figures/4.1/vein_Kupffer_DGE_hmFISH.csv")

###plot Vcam1 in boxplot 
boxplot_KCs_per_spatial_area(veins,"Vcam1","purple","./figures/4.1/",".pdf")
boxplot_KCs_per_spatial_area(veins,"Vcam1","purple","./figures/4.1/",".svg")

