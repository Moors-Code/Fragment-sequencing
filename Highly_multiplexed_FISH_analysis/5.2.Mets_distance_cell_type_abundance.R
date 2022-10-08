########## Part 5.2 Mets distance cell type abundance analysis ##########
#This part validates findings of cell type abundance analysis of sphere-seq 

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Sphere-sequencing/Highly_multiplexed_FISH_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")
source("./functions_and_packages/4.Functions_cell_type_abundance.R")

###Load data 
merged <- readRDS(file = "./data_files_generated/Resolve_seurat_anno.rds")

#subset only samples with visible metastasis 
Idents(merged) <- "samples"
mets_samples <- subset(merged, idents = "mets")

########## DA abundance analysis ##########
DA_analysis_cell_type_abundance_Mets_distance_hmFISH(mets_samples,"././figures/5.2/")

top <- read.csv("./figures/5.2//Cell_type_abundance_Mets_distance_hmFISH.csv")

top <- top %>% 
  mutate(
    Expression = case_when(logFC >= 0.5 & PValue <= 0.05 ~ "High in proximal",
                           logFC <= -0.5 & PValue <= 0.05 ~ "High in distal",
                           TRUE ~ "Non sig.")
  )

p <- ggplot(top, aes(x=logFC, y=-log10(PValue))) +
  geom_point(aes(color = Expression),size=5) +
  geom_text(data=top[top$PValue<1 & abs(top$logFC) > 0,], aes(label=Gene),size=8) +
  xlab("logFC") + 
  ylab("-log10(PValue)") + ggtitle("Cell type prop - hmFISH (PValue ≤ 0.05, logFC >0.5") + 
  scale_color_manual(values = c("firebrick3","dodgerblue3", "gray50"),guide = "none") + theme_classic() + 
  theme(axis.title= element_text(size = 25)) + theme(axis.text = element_text(size = 30))  + 
  theme(plot.title = element_text(size = 25, face = "bold"))  + 
  guides(colour = guide_legend(override.aes = list(size=1.5))) + 
  theme(legend.title = element_text(size = 15), legend.text = element_text(size = 15)) +  
  geom_vline(xintercept = c(-0.5,0.5), linetype = "dashed", color = "gray50") 
p + ggsave("./figures/5.2/Cell_type_prop_metsDistance_hmFISH_vulcano.pdf",width = 12, height = 10)
p + ggsave("./figures/5.2/Cell_type_prop_metsDistance_hmFISH_vulcano.svg",width = 12, height = 10)  

###Boxplots of some cell types 
spatial_feature_number_ids <- as.data.frame(table(mets_samples$spatial_feature_number))$Var1
for (i in spatial_feature_number_ids) {
  Idents(mets_samples) <- "spatial_feature_number"
  j <- subset(mets_samples, idents = i)
  create_table_cell_type_prop_resolve(j, "Mets_distance","annotation","./figures/5.2/", i)
}

#read in S1 areas
df1 <- read.csv("./figures/5.2/s1_cv_nm_1_proportions_Mets_distance_annotation.csv")
df2 <- read.csv("./figures/5.2/s1_cv_nm_2_proportions_Mets_distance_annotation.csv")
df3 <- read.csv("./figures/5.2/s1_cv_nm_3_proportions_Mets_distance_annotation.csv")
df4 <- read.csv("./figures/5.2/s1_cv_nm_4_proportions_Mets_distance_annotation.csv")
df5 <- read.csv("./figures/5.2/s1_cv_nm_5_proportions_Mets_distance_annotation.csv")
df6 <- read.csv("./figures/5.2/s1_cv_nm_6_proportions_Mets_distance_annotation.csv")
df7 <- read.csv("./figures/5.2/s1_cv_nm_7_proportions_Mets_distance_annotation.csv")
df8 <- read.csv("./figures/5.2/s1_cv_nm_8_proportions_Mets_distance_annotation.csv")
df9 <- read.csv("./figures/5.2/s1_cv_nm_9_proportions_Mets_distance_annotation.csv")
df10 <- read.csv("./figures/5.2/s1_cv_nm_10_proportions_Mets_distance_annotation.csv")

df11 <- read.csv("./figures/5.2/s1_cv_sm_1_proportions_Mets_distance_annotation.csv")
df12 <- read.csv("./figures/5.2/s1_cv_sm_2_proportions_Mets_distance_annotation.csv")
df13 <- read.csv("./figures/5.2/s1_cv_sm_3_proportions_Mets_distance_annotation.csv")
df14 <- read.csv("./figures/5.2/s1_cv_sm_4_proportions_Mets_distance_annotation.csv")
df15 <- read.csv("./figures/5.2/s1_cv_sm_5_proportions_Mets_distance_annotation.csv")
df16 <- read.csv("./figures/5.2/s1_cv_sm_6_proportions_Mets_distance_annotation.csv")
df17 <- read.csv("./figures/5.2/s1_cv_sm_7_proportions_Mets_distance_annotation.csv")
df18 <- read.csv("./figures/5.2/s1_cv_sm_8_proportions_Mets_distance_annotation.csv")
df19 <- read.csv("./figures/5.2/s1_cv_sm_9_proportions_Mets_distance_annotation.csv")
df20 <- read.csv("./figures/5.2/s1_cv_sm_10_proportions_Mets_distance_annotation.csv")
df21 <- read.csv("./figures/5.2/s1_cv_sm_11_proportions_Mets_distance_annotation.csv")
df22 <- read.csv("./figures/5.2/s1_cv_sm_12_proportions_Mets_distance_annotation.csv")
df23 <- read.csv("./figures/5.2/s1_cv_sm_13_proportions_Mets_distance_annotation.csv")

df24 <- read.csv("./figures/5.2/s1_m1_proportions_Mets_distance_annotation.csv")
df25 <- read.csv("./figures/5.2/s1_m2_proportions_Mets_distance_annotation.csv")
df26 <- read.csv("./figures/5.2/s1_m3_proportions_Mets_distance_annotation.csv")
df27 <- read.csv("./figures/5.2/s1_m4_proportions_Mets_distance_annotation.csv")
df28 <- read.csv("./figures/5.2/s1_m5_proportions_Mets_distance_annotation.csv")

df29 <- read.csv("./figures/5.2/s1_pv_nm_1_proportions_Mets_distance_annotation.csv")
df31 <- read.csv("./figures/5.2/s1_pv_nm_3_proportions_Mets_distance_annotation.csv")
df32 <- read.csv("./figures/5.2/s1_pv_nm_4_proportions_Mets_distance_annotation.csv")
df33 <- read.csv("./figures/5.2/s1_pv_nm_5_proportions_Mets_distance_annotation.csv")
df34 <- read.csv("./figures/5.2/s1_pv_nm_6_proportions_Mets_distance_annotation.csv")
df35 <- read.csv("./figures/5.2/s1_pv_nm_7_proportions_Mets_distance_annotation.csv")
df36 <- read.csv("./figures/5.2/s1_pv_nm_8_proportions_Mets_distance_annotation.csv")
df37 <- read.csv("./figures/5.2/s1_pv_nm_9_proportions_Mets_distance_annotation.csv")
df38 <- read.csv("./figures/5.2/s1_pv_nm_10_proportions_Mets_distance_annotation.csv")
df39 <- read.csv("./figures/5.2/s1_pv_nm_11_proportions_Mets_distance_annotation.csv")

df40 <- read.csv("./figures/5.2/s1_pv_sm_1_proportions_Mets_distance_annotation.csv")
df41 <- read.csv("./figures/5.2/s1_pv_sm_2_proportions_Mets_distance_annotation.csv")
df42 <- read.csv("./figures/5.2/s1_pv_sm_3_proportions_Mets_distance_annotation.csv")
df43 <- read.csv("./figures/5.2/s1_pv_sm_4_proportions_Mets_distance_annotation.csv")
df44 <- read.csv("./figures/5.2/s1_pv_sm_5_proportions_Mets_distance_annotation.csv")
df45 <- read.csv("./figures/5.2/s1_pv_sm_6_proportions_Mets_distance_annotation.csv")

#S4 samples 
df46 <- read.csv("./figures/5.2/s4_cv_nm_1_proportions_Mets_distance_annotation.csv")
df47 <- read.csv("./figures/5.2/s4_cv_nm_2_proportions_Mets_distance_annotation.csv")
df48 <- read.csv("./figures/5.2/s4_cv_nm_3_proportions_Mets_distance_annotation.csv")
df49 <- read.csv("./figures/5.2/s4_cv_nm_4_proportions_Mets_distance_annotation.csv")
df50 <- read.csv("./figures/5.2/s4_cv_nm_5_proportions_Mets_distance_annotation.csv")
df51 <- read.csv("./figures/5.2/s4_cv_nm_6_proportions_Mets_distance_annotation.csv")

df52 <- read.csv("./figures/5.2/s4_cv_sm_1_proportions_Mets_distance_annotation.csv")
df53 <- read.csv("./figures/5.2/s4_cv_sm_2_proportions_Mets_distance_annotation.csv")
df54 <- read.csv("./figures/5.2/s4_cv_sm_3_proportions_Mets_distance_annotation.csv")
df55 <- read.csv("./figures/5.2/s4_cv_sm_4_proportions_Mets_distance_annotation.csv")
df56 <- read.csv("./figures/5.2/s4_cv_sm_5_proportions_Mets_distance_annotation.csv")
df57 <- read.csv("./figures/5.2/s4_cv_sm_6_proportions_Mets_distance_annotation.csv")
df58 <- read.csv("./figures/5.2/s4_cv_sm_7_proportions_Mets_distance_annotation.csv")

df59 <- read.csv("./figures/5.2/s4_m1_proportions_Mets_distance_annotation.csv")
df60 <- read.csv("./figures/5.2/s4_m2_proportions_Mets_distance_annotation.csv")
df61 <- read.csv("./figures/5.2/s4_m3_proportions_Mets_distance_annotation.csv")

df62 <- read.csv("./figures/5.2/s4_pv_nm_1_proportions_Mets_distance_annotation.csv")
df64 <- read.csv("./figures/5.2/s4_pv_nm_3_proportions_Mets_distance_annotation.csv")
df65 <- read.csv("./figures/5.2/s4_pv_nm_4_proportions_Mets_distance_annotation.csv")
df66 <- read.csv("./figures/5.2/s4_pv_nm_5_proportions_Mets_distance_annotation.csv")
df67 <- read.csv("./figures/5.2/s4_pv_nm_6_proportions_Mets_distance_annotation.csv")

df68 <- read.csv("./figures/5.2/s4_pv_sm_1_proportions_Mets_distance_annotation.csv")
df69 <- read.csv("./figures/5.2/s4_pv_sm_2_proportions_Mets_distance_annotation.csv")
df70 <- read.csv("./figures/5.2/s4_pv_sm_3_proportions_Mets_distance_annotation.csv")

#S5 samples 
df71 <- read.csv("./figures/5.2/s5_cv_nm_1_proportions_Mets_distance_annotation.csv")
df72 <- read.csv("./figures/5.2/s5_cv_nm_2_proportions_Mets_distance_annotation.csv")

df73 <- read.csv("./figures/5.2/s5_cv_sm_1_proportions_Mets_distance_annotation.csv")
df74 <- read.csv("./figures/5.2/s5_cv_sm_2_proportions_Mets_distance_annotation.csv")
df75 <- read.csv("./figures/5.2/s5_cv_sm_3_proportions_Mets_distance_annotation.csv")
df76 <- read.csv("./figures/5.2/s5_cv_sm_4_proportions_Mets_distance_annotation.csv")
df77 <- read.csv("./figures/5.2/s5_cv_sm_5_proportions_Mets_distance_annotation.csv")

df78 <- read.csv("./figures/5.2/s5_m1_proportions_Mets_distance_annotation.csv")
df79 <- read.csv("./figures/5.2/s5_m2_proportions_Mets_distance_annotation.csv")
df80 <- read.csv("./figures/5.2/s5_m3_proportions_Mets_distance_annotation.csv")

df81 <- read.csv("./figures/5.2/s5_pv_nm_1_proportions_Mets_distance_annotation.csv")

df82 <- read.csv("./figures/5.2/s5_pv_sm_1_proportions_Mets_distance_annotation.csv")
df83 <- read.csv("./figures/5.2/s5_pv_sm_2_proportions_Mets_distance_annotation.csv")
df84 <- read.csv("./figures/5.2/s5_pv_sm_3_proportions_Mets_distance_annotation.csv")

#combine all dataframes 
df_all <- rbind.fill(df1, df2)
df_all <- rbind.fill(df_all,df3)
df_all <- rbind.fill(df_all,df4)
df_all <- rbind.fill(df_all,df5)
df_all <- rbind.fill(df_all,df6)
df_all <- rbind.fill(df_all,df7)
df_all <- rbind.fill(df_all,df8)
df_all <- rbind.fill(df_all,df9)
df_all <- rbind.fill(df_all,df10)
df_all <- rbind.fill(df_all,df11)
df_all <- rbind.fill(df_all,df12)
df_all <- rbind.fill(df_all,df13)
df_all <- rbind.fill(df_all,df14)
df_all <- rbind.fill(df_all,df15)
df_all <- rbind.fill(df_all,df16)
df_all <- rbind.fill(df_all,df17)
df_all <- rbind.fill(df_all,df18)
df_all <- rbind.fill(df_all,df19)
df_all <- rbind.fill(df_all,df20)
df_all <- rbind.fill(df_all,df21)
df_all <- rbind.fill(df_all,df22)
df_all <- rbind.fill(df_all,df23)
df_all <- rbind.fill(df_all,df24)
df_all <- rbind.fill(df_all,df25)
df_all <- rbind.fill(df_all,df26)
df_all <- rbind.fill(df_all,df27)
df_all <- rbind.fill(df_all,df28)
df_all <- rbind.fill(df_all,df29)
df_all <- rbind.fill(df_all,df31)
df_all <- rbind.fill(df_all,df32)
df_all <- rbind.fill(df_all,df33)
df_all <- rbind.fill(df_all,df34)
df_all <- rbind.fill(df_all,df35)
df_all <- rbind.fill(df_all,df36)
df_all <- rbind.fill(df_all,df37)
df_all <- rbind.fill(df_all,df38)
df_all <- rbind.fill(df_all,df39)
df_all <- rbind.fill(df_all,df40)
df_all <- rbind.fill(df_all,df41)
df_all <- rbind.fill(df_all,df42)
df_all <- rbind.fill(df_all,df43)
df_all <- rbind.fill(df_all,df44)
df_all <- rbind.fill(df_all,df45)
df_all <- rbind.fill(df_all,df46)
df_all <- rbind.fill(df_all,df47)
df_all <- rbind.fill(df_all,df48)
df_all <- rbind.fill(df_all,df49)
df_all <- rbind.fill(df_all,df50)
df_all <- rbind.fill(df_all,df51)
df_all <- rbind.fill(df_all,df52)
df_all <- rbind.fill(df_all,df53)
df_all <- rbind.fill(df_all,df54)
df_all <- rbind.fill(df_all,df55)
df_all <- rbind.fill(df_all,df56)
df_all <- rbind.fill(df_all,df57)
df_all <- rbind.fill(df_all,df58)
df_all <- rbind.fill(df_all,df59)
df_all <- rbind.fill(df_all,df60)
df_all <- rbind.fill(df_all,df61)
df_all <- rbind.fill(df_all,df62)
df_all <- rbind.fill(df_all,df64)
df_all <- rbind.fill(df_all,df65)
df_all <- rbind.fill(df_all,df66)
df_all <- rbind.fill(df_all,df67)
df_all <- rbind.fill(df_all,df68)
df_all <- rbind.fill(df_all,df69)
df_all <- rbind.fill(df_all,df70)
df_all <- rbind.fill(df_all,df71)
df_all <- rbind.fill(df_all,df72)
df_all <- rbind.fill(df_all,df73)
df_all <- rbind.fill(df_all,df74)
df_all <- rbind.fill(df_all,df75)
df_all <- rbind.fill(df_all,df76)
df_all <- rbind.fill(df_all,df77)
df_all <- rbind.fill(df_all,df78)
df_all <- rbind.fill(df_all,df79)
df_all <- rbind.fill(df_all,df80)
df_all <- rbind.fill(df_all,df81)
df_all <- rbind.fill(df_all,df82)
df_all <- rbind.fill(df_all,df83)
df_all <- rbind.fill(df_all,df84)


#Kupffer
boxplot_cell_prop(df_all,"Kupffer","#E20FE8","./figures/5.2/",".pdf")
boxplot_cell_prop(df_all,"Kupffer","#E20FE8","./figures/5.2/",".svg")

#LECs
boxplot_cell_prop(df_all,"LECs","#F4E740","./figures/5.2/",".pdf")
boxplot_cell_prop(df_all,"LECs","#F4E740","./figures/5.2/",".svg")

#Monoytes
boxplot_cell_prop(df_all,"Monocytes","#19E80F","./figures/5.2/",".pdf")
boxplot_cell_prop(df_all,"Monocytes","#19E80F","./figures/5.2/",".svg")

#Metastatic cells
boxplot_cell_prop(df_all,"Metastasis","#F90606","./figures/5.2/",".pdf")
boxplot_cell_prop(df_all,"Metastasis","#F90606","./figures/5.2/",".svg")

