########## Part 2.5.2: Liver lobule layer zonated gene expression analysis##########
#This part investigates zonation specific gene expression in LECs and KCs

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")
source("./functions_and_packages/5.Functions_lobule_layers.R")

###load R object
liverSpS5C <- readRDS("./data_files_generated/LiverMerged_afterBC_anno_BS_5cells_zC_lobules.Rda")

#WT sample 
Idents(liverSpS5C) <- "orig.ident"
WT <- subset(liverSpS5C, idents = "WTCre")

#Injected samples  
Idents(liverSpS5C) <- "orig.ident"
Injected <- subset(liverSpS5C, idents = c("M1","M3","M4","4M1","4M2","M5","6M1","6M2","6M3"))

########## DE analysis along lobule layers ##########
###LECs 
DE_zonated_genes(Injected,"annotation.broad","LECs","./data_files_generated/","Injected")

##plot results of DE analysis in vulcano plot, only consider genes that were not in the landmark gene panel for ZC calculation 
#read in results 
top <- read.csv("./data_files_generated/Injected_LECs_LM_zonated_genes.csv")

#genes used in ZC calculation 
PV_genes <- c("Sh3rf1","Ppp1r9b","Ctdspl","Serpina3f", "Impdh1","Meis1","Lama4", "Armcx1",
              "Txndc16","Klhl3","Lpcat1","Tnik","Gdf2","Trim30a", "Rnf157", "Lyve1","Pcolce",
              "Osmr","Enpp6","Pla2r1","Pygb", "Pld1","Ntm","Il33","Igfbp3","Hoxb5","Art3","Tgfb1",
              "B4galt4","Klf4","Sod3","Fgfr1", "Ly6a","Hlx","Slc43a2", "Rapgef3","Sept8","Itgb3",
              "Efnb2","Twist1","Sox17","Itga9", "Rasd1", "Slc41a1","Spats2l","Gata2", "Col1a2",
              "Chst2","Nid2","Tmem44","Adam23","Prkcq","Tomm40l", "Dll4", "Il1a","Ntn4", "Msr1")

CV_genes <- c("Wnt9b", "Rspo3", "Cdh13","Thbd","Lmcd1","Ar", "Ier3","Cc2d2a", "Dkk3",
              "Fabp4","Dennd2d","Cln6","Lfng","Cebpd","Tmed8","Rgp1", "Pfkfb3","Wnt2", "Robo1",
              "Fam84b","Rab3b", "Ptgs1","P2ry1","Olfm1", "Amigo2","Trappc2", "Kit",
              "Mpp1","Ap4m1","Nkd1","Sgsh","Zxdb","Mbd5","Atr","Slc48a1", "Gen1","Pigh","Unkl","Itgb3bp",
              "Zmat1","Med9","Anxa3","Kcnb1","Gja4","Eid1","Arl16","Orc2","Hs1bp3","Gorasp1",
              "Skp1a", "Med21", "Ptcd1","Rpa1","Dbp","Tsr2","Diablo")

genes_known <- c(PV_genes,CV_genes)

#remove known zonation specific genes 
top <- top[!(top$Gene %in% genes_known),]

#apply significant cutoff of FDR and logFC to highlight significant ones in vulcano plot  
top <- top %>% 
  mutate(
    Expression = case_when(logFC >= 0.5 & FDR <= 0.05 ~ "High in PV",
                           logFC <= -0.5 & FDR <= 0.05 ~ "High in CV",
                           TRUE ~ "Non sig.")
  )

#highlight specific genes that are also present in panel of highly multiplexed FISH 
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
  ylab("-log10(FDR)") + ggtitle("LECs - Sphere-seq (FDR ≤ 0.05, logFC >0.5") + 
  scale_color_manual(values = c("red", "darkgoldenrod2", "gray50"),guide = "none") + theme_classic() + 
  theme(axis.title= element_text(size = 25)) + theme(axis.text = element_text(size = 30))  + 
  theme(plot.title = element_text(size = 25, face = "bold"))  + 
  guides(colour = guide_legend(override.aes = list(size=1.5))) + 
  theme(legend.title = element_text(size = 15), legend.text = element_text(size = 15)) + 
  geom_vline(xintercept = c(-0.5,0.5), linetype = "dashed", color = "gray50")
p + ggsave("./figures/2.5.2/LECs_CV_PV_SpS_vulcano.pdf",width = 12, height = 10)
p + ggsave("./figures/2.5.2/LECs_CV_PV_SpS_vulcano.svg",width = 12, height = 10)  

##plot some genes in boxplot 
#CV zonated 
boxplot_zonation_genes_LECs_injected(Injected,"Plpp1","red","./figures/2.5.2/","Injected",".pdf")
boxplot_zonation_genes_LECs_injected(Injected,"Plpp1","red","./figures/2.5.2/","Injected",".svg")
boxplot_zonation_genes_LECs_injected(Injected,"Lhx6","red","./figures/2.5.2/","Injected",".pdf")
boxplot_zonation_genes_LECs_injected(Injected,"Lhx6","red","./figures/2.5.2/","Injected",".svg")

#PV zonated 
boxplot_zonation_genes_LECs_injected(Injected,"Cd36","darkgoldenrod2","./figures/2.5.2/","Injected",".pdf")
boxplot_zonation_genes_LECs_injected(Injected,"Cd36","darkgoldenrod2","./figures/2.5.2/","Injected",".svg")
boxplot_zonation_genes_LECs_injected(Injected,"Galnt15","darkgoldenrod2","./figures/2.5.2/","Injected",".pdf")
boxplot_zonation_genes_LECs_injected(Injected,"Galnt15","darkgoldenrod2","./figures/2.5.2/","Injected",".svg")

###Kupffer cells 
DE_zonated_genes(Injected,"annotation.broad","Kupffer","./data/","Injected")

##plot results of DE analysis in vulcano plot
#read in results 
top <- read.csv("./data/Injected_Kupffer_LM_zonated_genes.csv")

top <- top %>% 
  mutate(
    Expression = case_when(logFC >= 0.5 & FDR <= 0.05 ~ "High in PV",
                           logFC <= -0.5 & FDR <= 0.05 ~ "High in CV",
                           TRUE ~ "Non sig.")
  )

#highlight interesting genes 
top$genelabels <- ifelse(top$Gene == "Vcam1"
                         |top$Gene == "Ccl5"
                         |top$Gene == "Asb2"
                         |top$Gene == "Ccl24"
                         |top$Gene == "Marco"
                         |top$Gene == "Ccl2"
                         |top$Gene == "Cxcl13"
                         |top$Gene == "Il6"
                         |top$Gene == "Cxcl1"
                         |top$Gene == "Cd207",
                         TRUE,FALSE)

p <- ggplot(top, aes(x=logFC, y=-log10(FDR))) +
  geom_point(aes(color = Expression),size=3) +
  geom_text_repel(aes(label=ifelse(top$genelabels, top$Gene,"")),size=8, max.overlaps = 100) +
  xlab("logFC") + 
  ylab("-log10(FDR)") + ggtitle("Kupffer - Sphere-seq (FDR ≤ 0.05, logFC >0.5") + 
  scale_color_manual(values = c("red", "darkgoldenrod2", "gray50"),guide = "none") + theme_classic() + 
  theme(axis.title= element_text(size = 25)) + theme(axis.text = element_text(size = 30))  + 
  theme(plot.title = element_text(size = 25, face = "bold"))  + 
  guides(colour = guide_legend(override.aes = list(size=1.5))) + 
  theme(legend.title = element_text(size = 15), legend.text = element_text(size = 15)) + 
  geom_vline(xintercept = c(-0.5,0.5), linetype = "dashed", color = "gray50")
p + ggsave("./figures/2.5.2/Kupffer_CV_PV_SpS_vulcano.pdf",width = 12, height = 10)
p + ggsave("./figures/2.5.2/Kupffer_CV_PV_SpS_vulcano.svg",width = 12, height = 10)  

##plot Vcam1 in boxplot 
boxplot_zonation_genes_KCs_injected(Injected,"Vcam1","purple","./figures/2.5.2/","Injected",".pdf")
boxplot_zonation_genes_KCs_injected(Injected,"Vcam1","purple","./figures/2.5.2/","Injected",".svg")


########## WT sample analysis ##########
#only one sample, therefore I applied a wilcox rank sum test to compare CV and PV gene expression 

###LECs 
boxplot_zonation_genes_LECs_WT(WT,"Plpp1","red","./figures/2.5.2/","WT",".pdf")
boxplot_zonation_genes_LECs_WT(WT,"Plpp1","red","./figures/2.5.2/","WT",".svg")
boxplot_zonation_genes_LECs_WT(WT,"Lhx6","red","./figures/2.5.2/","WT",".pdf")
boxplot_zonation_genes_LECs_WT(WT,"Lhx6","red","./figures/2.5.2/","WT",".svg")

boxplot_zonation_genes_LECs_WT(WT,"Cd36","darkgoldenrod2","./figures/2.5.2/","WT",".pdf")
boxplot_zonation_genes_LECs_WT(WT,"Cd36","darkgoldenrod2","./figures/2.5.2/","WT",".svg")
boxplot_zonation_genes_LECs_WT(WT,"Galnt15","darkgoldenrod2","./figures/2.5.2/","WT",".pdf")
boxplot_zonation_genes_LECs_WT(WT,"Galnt15","darkgoldenrod2","./figures/2.5.2/","WT",".svg")

