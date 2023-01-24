#
library(reticulate)
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(tidyr)
library(nichenetr)
library(tidyverse)
library(stringr)
library(readxl)
library(plyr)
library(ComplexHeatmap)
library(circlize)


location = "mac"
cancer_name = c("Lung","Liver","Colon","Ovary","Breast")

if (location == "mac"){  wd="/Users/cynthia_ye/Documents/OneDrive - Linköpings universitet"}
if (location == "omika"){ wd="/home/yelzh67/Projects"}

outdir= paste0(wd,"/cancer/datafile/outputs_5cancers/new_outputs/Nichenet")
if (dir.exists(outdir)==F){dir.create(outdir)
  print("outdir created")}
input.dir=paste0(wd,"/cancer/datafile/outputs_5cancers/new_outputs/Seurat_subset_celltyping/2022Feb12")
dir.create(input.dir)
############################################################################################
#Step 0:  NicheNet’s ligand-target prior model
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds")) #ligand_target_matrix
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns

lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
head(lr_network)

weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))
head(weighted_networks$lr_sig)

DEG_method = "MAST"

for (cancer in cancer_name){
  # #1.import DEGs####
  # DEGs = read.table(file = paste0(wd,"/cancer/datafile/outputs_5cancers/new_outputs/DEGs/DEGs_of_",cancer,"/DEGs_merged_",cancer,"_",DEG_method,".txt"), header= T, sep="\t", stringsAsFactor = T)  
  # DEGs = DEGs[, apply(DEGs, 2, function(y) any(!is.na(y)))]
  # head(DEGs)
  # 
  # #2.background genes####
  # BGs = read.table(file = paste0(wd,"/cancer/datafile/outputs_5cancers/new_outputs/DEGs/Background_genes/background_genes_",cancer,".csv"), header= T, sep=",", stringsAsFactor = T)
  # head(BGs)
  # BGs = BGs[,2:dim(BGs)[2]]
  # 
  #3.setup for nichenet####
  outputdir <- paste0(outdir,"/",cancer)
  if (dir.exists(outputdir)==F){dir.create(outputdir)
    print("outputdir created")}
  # 
  # all_expressed_genes <- BGs
  # colnames(all_expressed_genes)
  # dim(all_expressed_genes)
  # 
  # DEGs_human <- DEGs
  # colnames(DEGs_human)
  # DEGs_human[1:5,1:7]
  # 
  # ##3.1 looping for ligands activity for single dataset ####
  # source(paste0(wd,"/cancer/scripts/scripts_new_remerge/3.0 Looping for ligands activity.R"))
  # 
  # write.table(all_ligand_activity_Yelin, 
  #             file=paste0(outputdir,"/all_ligand_activity.csv"), 
  #             sep=",", quote=F, col.names=T,row.names=F)
  # 
  # ##3.2 create matrix for all_ligand_activity_Yelin
  # all_ligand_activity_Yelin <- read.table(file = paste0(outputdir,"/all_ligand_activity.csv"), header= T, sep=",",stringsAsFactor = F)
  # all_ligand_activity_Yelin <- as.data.frame(all_ligand_activity_Yelin)
  # is.atomic(all_ligand_activity_Yelin)
  # dim(all_ligand_activity_Yelin)
  # length(unique(all_ligand_activity_Yelin$test_ligand))
  # table(all_ligand_activity_Yelin$Sender,all_ligand_activity_Yelin$Target)
  # dim(table(all_ligand_activity_Yelin$Sender,all_ligand_activity_Yelin$Target))
  # 
  # #4. make cluster interaction matrix - for cytoscape####
  # cluster_intera_ye <- data.frame(table(all_ligand_activity_Yelin$Sender,all_ligand_activity_Yelin$Target))
  # colnames(cluster_intera_ye) = c("Sender", "Target", "No.Interactions")
  # 
  # write.table(cluster_intera_ye, 
  #             file=paste0(outputdir,"/cluster_interaction.txt"), 
  #             sep="\t", quote=F, col.names=T,row.names=F)
  # dim(cluster_intera_ye)
  # 
  # #5. calculate the most frequent ligands from all clusters####
  # ligand_freq_ranked  <- data.frame(table(all_ligand_activity_Yelin$test_ligand))
  # ligand_freq_ranked <- ligand_freq_ranked %>% arrange(-Freq)
  # write.table(ligand_freq_ranked, 
  #             file=paste0(outputdir,"/ligand_freq_ranked.txt"),
  #             sep="\t", quote=F, col.names=T,row.names=F)
  # 
  # #6. positive PCC only####
  # #6.1  get PCC positive ligands
  # PCC_positive_ligand  = all_ligand_activity_Yelin %>% filter(pearson > 0)
  # length(unique(PCC_positive_ligand$test_ligand))
  # dim(PCC_positive_ligand)
  # write.table(PCC_positive_ligand, 
  #             file=paste0(outputdir,"/PCC_positive_ligand.txt"), 
  #             sep="\t", quote=F, col.names=T,row.names=F)
  PCC_positive_ligand = read.table(file =paste0(outputdir,"/PCC_positive_ligand.txt"), header= T, sep=",",stringsAsFactor = F)
  # 
  # #6.2 calculate the most frequent ligands from all clusters####
  # PCC_positive_ligand_freq_ranked <- data.frame(table(PCC_positive_ligand$test_ligand))
  # PCC_positive_ligand_freq_ranked <- PCC_positive_ligand_freq_ranked %>% arrange(-Freq)
  # 
  # write.table(PCC_positive_ligand_freq_ranked, 
  #             file=paste0(outputdir,"/PCC_positive_ligand_freq_ranked.txt"), 
  #             sep="\t", quote=F, col.names=T,row.names=F)
  # 
  # #6.3 make cluster interaction matrix - for cytoscape
  # PCC_positive_cluster_intera_ye <- data.frame(table(PCC_positive_ligand$Sender,PCC_positive_ligand$Target))
  # colnames(PCC_positive_cluster_intera_ye) = c("Sender", "Target", "No.Interactions")
  # 
  # write.table(PCC_positive_cluster_intera_ye, 
  #             file=paste0(outputdir,"/cluster_interaction_PCC_positive.txt"), 
  #             sep="\t", quote=F, col.names=T,row.names=F)
  # dim(PCC_positive_cluster_intera_ye)
  # 
  # positive_ligands_counts = matrix(NA,1,4)
  # rownames(positive_ligands_counts) = cancer
  # colnames(positive_ligands_counts) = c("N_positive_ligands","N_all_ligands","N_positive_interactions","N_all_interactions")
  # positive_ligands_counts[1,1] = length(unique(PCC_positive_ligand$test_ligand))
  # positive_ligands_counts[1,2] = length(unique(all_ligand_activity_Yelin$test_ligand))
  # positive_ligands_counts[1,3] = dim(PCC_positive_ligand)[1]
  # positive_ligands_counts[1,4] = dim(all_ligand_activity_Yelin)[1]
  # write.table(positive_ligands_counts, 
  #             file=paste0(outputdir,"/positive_ligands_counts.txt"), 
  #             sep="\t", quote=F, col.names=T,row.names=T) 
  
  #7. Centrality####
  source(paste0(wd,"/cancer/scripts/Cancer_Project_shared_unique/3.0 NicheNet_centrality_single_dataset.R"))

  print(paste0(cancer, "NicheNet done!"))
}

##7.1 Centrality heatmap ####
centrality_merge <- matrix(NA, ncol=11) %>% as.data.frame()
colnames(centrality_merge) = c("node_degree_all","node_degree_in","node_degree_out","closeness","eigenvector centralities","Kleinberg's hub centrality scores","Laplacian Centrality","Leverage Centrality","Local Bridging Centrality","cancer","celltype")
for (cancer in cancer_name){
  a = read.csv(paste0(wd,"/cancer/datafile/outputs_5cancers/new_outputs/Nichenet/",cancer,"/Cell-cell_centrality_summary_PCC_positive.txt"),sep = "\t",row.names = "X") %>% as.data.frame()
  colnames(a) = lapply(colnames(a), function(x) {strsplit(x,"_") %>% unlist %>% .[2]})  %>% unlist  
  a = t(a)  %>% as.data.frame()
  a$cancer = cancer
  a$celltype = rownames(a)
  centrality_merge = rbind(centrality_merge,a)
}
centrality_merge = centrality_merge[2:dim(centrality_merge)[1],]

library(reshape2)
df = centrality_merge
head(df)
df.wide = dcast(df,celltype~cancer, value.var ="eigenvector centralities")
rowname = df.wide$celltype
df.wide = df.wide[,2:dim(df.wide)[2]]
df.wide[1:5,1:4]
df.wide[is.na(df.wide)==T] <- 0
df.wide = apply(df.wide,2,as.numeric)
rownames(df.wide) = rowname

dpi = 300
png(file=paste0(outdir,"/heatmap_eigenvector_centrality.png"), width = dpi*10, height = dpi*5, units = "px",res = dpi)
p = Heatmap(t(df.wide),
          # split
          # row_split = 3,
          # column_split = 2,  #split by dendrograms
          # row_km = 2, 
          column_km = 2 , #split rows and colums
          # split = data.frame(cyl = mtcars$cyl, am = mtcars$am), #根据数据进行分割
          # column_split = anno_col$Group,
          #设置title,legend,font
          name = "Eigenvector centrality", #title of legend
          column_title = "Cell type",
          column_title_side = "bottom",  #设置列标题的位置，可选"top"或"bottom"
          column_title_gp = gpar(fontsize = 15, fontface = "bold"),  #更改列文本的字体
          column_names_gp = gpar(fontsize = 15),
          row_title = "Cancers",
          row_title_side ="right",  #设置行标题的位置，可选"left"或者"right"
          row_title_gp = gpar(fontsize = 15) ,
          row_names_gp = gpar(fontsize = 15), # Text size for row names
          # heatmap_legend_param = list(title="legend"),
          # top_annotation = ha, #put annotation
          cluster_columns = T,
          cluster_rows = T,
          # cluster_columns = cluster_within_group(df,anno_col$Celltype),  #组内聚类
          #设置颜色
          col = colorRamp2(c(0,1), c("white", "red")),
          # col = colorpanel(50,"blue","white","red"),
          # col = colorRamp2(c(0, 10, 50), c("blue", "white", "dark red"))  #设置绘图的颜色，可以根据数据实际情况调整
          # na_col = "grey",
          # clustering_distance_rows = "pearson"
          border = T
) 
print(p)
dev.off()

##7.1 Centrality ranking heatmap #### 
library(reshape2)
df = read.csv(paste0(wd,"/cancer/datafile/outputs_5cancers/new_outputs/GWAS/GWAS_centrality_corr/centrality_rank.csv"),row.names = 1)

head(df)
df$sum = apply(df,1,sum)
df = df %>% arrange(sum)

p = Heatmap(t(df[,1:5]),
            # split
            # row_split = 3,
            # column_split = 2,  #split by dendrograms
            # row_km = 2, 
            column_km = 2 , #split rows and colums
            # split = data.frame(cyl = mtcars$cyl, am = mtcars$am), #根据数据进行分割
            # column_split = anno_col$Group,
            #设置title,legend,font
            rect_gp = gpar(col = "white"),
            name = "Rank by centrality", #title of legend
            column_title = "Cell type",
            column_title_side = "bottom",  #设置列标题的位置，可选"top"或"bottom"
            column_title_gp = gpar(fontsize = 15, fontface = "bold"),  #更改列文本的字体
            column_names_gp = gpar(fontsize = 15),
            row_title = "Cancers",
            row_title_side ="right",  #设置行标题的位置，可选"left"或者"right"
            row_title_gp = gpar(fontsize = 15) ,
            row_names_gp = gpar(fontsize = 15), # Text size for row names
            # heatmap_legend_param = list(title="legend"),
            # top_annotation = ha, #put annotation
            cluster_columns = T,
            cluster_rows = T,
            # cluster_columns = cluster_within_group(df,anno_col$Celltype),  #组内聚类
            #设置颜色
            col = colorRamp2(c(13,1), c("white", "red")),
            # col = colorpanel(50,"blue","white","red"),
            # col = colorRamp2(c(0, 10, 50), c("blue", "white", "dark red"))  #设置绘图的颜色，可以根据数据实际情况调整
            # na_col = "grey",
            # clustering_distance_rows = "pearson"
            border = T
) 

pdf(file=paste0(outdir,"/heatmap_eigenvector_centrality_ranking.pdf"), width = 10, height = 5)
print(p)
dev.off()





