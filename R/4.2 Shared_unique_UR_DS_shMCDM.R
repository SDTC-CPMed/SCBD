
# devtools::install_github('satijalab/seurat-data')
#library(SeuratData)
library(Matrix)
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

location = "mac"
cancer_name = c("Lung","Liver","Colon","Ovary","Breast")

if (location == "mac"){  wd="/Users/cynthia_ye/Documents/OneDrive - Linköpings universitet"}
if (location == "omika"){ wd="/home/yelzh67/Projects"}

outdir= paste0(wd,"/cancer/datafile/outputs_5cancers/new_outputs/MCDM")
if (dir.exists(outdir)==F){dir.create(outdir)
  print("outdir created")}

outdir_shDS_of_shUR = paste0(outdir,"/Ligand_target_shared_DS_of_shared_UR/")
if (dir.exists(outdir_shDS_of_shUR)==F){dir.create(outdir_shDS_of_shUR)} else {print("outdir_shDS_of_shUR exsists")}

outdir_allDS_of_allUR = paste0(outdir,"/Ligand_target_allDS_of_allUR/")
if (dir.exists(outdir_allDS_of_allUR)==F){dir.create(outdir_allDS_of_allUR)} else {print("outdir_allDS_of_allUR exsists")}

#1.Find shared Ligands###########################
ligand_list = list()
for (cancer in cancer_name){
  ligands <- read.csv(paste0(wd,"/cancer/datafile/outputs_5cancers/new_outputs/Nichenet/",cancer,"/PCC_positive_ligand.txt"), sep="\t",header= T)
  ligand_list = c(ligand_list,list(ligands$test_ligand %>% unique ) )
}
names(ligand_list) = cancer_name
ligand_list = c(ligand_list['Breast'],ligand_list['Colon'],ligand_list['Liver'],ligand_list['Lung'],ligand_list['Ovary'])

shared_ligands <- intersect(unlist(ligand_list[1]), unlist(ligand_list[2])) %>% 
  intersect(.,unlist(ligand_list[3])) %>% 
  intersect(.,unlist(ligand_list[4])) %>% 
  intersect(.,unlist(ligand_list[5]))
length(shared_ligands)                     
write.table(shared_ligands, file=paste0(outdir,"/shared_UR.txt"), sep="\t", col.names = F, row.names = F)

shared_ov_breast_ligands <- intersect(unlist(ligand_list[4]), unlist(ligand_list[5]))
shared_lung_liver_colon_ligands <- intersect(unlist(ligand_list[1]), unlist(ligand_list[2])) %>% 
  intersect(.,unlist(ligand_list[3])) 

`%!in%` <- Negate(`%in%`)
shared_ov_breast_ligands[shared_ov_breast_ligands %!in% shared_lung_liver_colon_ligands]
shared_lung_liver_colon_ligands[shared_lung_liver_colon_ligands %!in% shared_ov_breast_ligands]

old_ligands <- read.csv("/Users/cynthia_ye/Documents/OneDrive/5. Linkoping University/2.5 Cancer/outputs/2021Aug20/shared_positive_ligand.txt",sep = "\t",header=F)
intersect(shared_ligands,as.vector(old_ligands$V1)) %>% length #compare with the results from the previous celltyping 

##1.1 venn plot of all URs####
library(VennDiagram)
library(gplots)
library(RColorBrewer)

#color
Chevalier1<-c("#355243","#fbca50","#c9d5d4","#baa28a")
FantasticFox1<-c("#d37a20","#dbcb09","#3a9cbc","#dd7208","#a30019")
Moonrise3<-c("#75cbdc","#f0a4af","#8a863a","#c2b479","#f8d068")
Darjeeling2<-c("#e6c09e","#0d5888","#cb8b3e","#9cd6d6","#000000")
Darjeeling1<-c("#fb0007","#139177","#ed9e08","#f56f08","#4caecc")
Royal2<-c("#e4c9b2","#f1c2a5","#f49d98","#fcd68f","#629076")
IsleofDogs2<-c("#e4c9b2","#998273","#a6723d","#2b2523","#151213")
colorset1 <- c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3")
colorset2 <- brewer.pal(5, "RdYlBu")
colorset3 <- paletteer_d("ggthemes::Classic_Green_Orange_6",5)
colorset4 <- c("#e97d72","#b49e33","#53b54c","#49acf4","#ed739b")
colorset5 <- c("#3e0d51","#42246b","#3f5086","#62b67c","#fae655")
scale_color_viridis_c()
cancer_color <- c("#ABDDA4","#FFFFBF","#FDAE61","#D7191C","#2B83BA")

venn_list <- ligand_list
venn.plot <- venn.diagram(x=venn_list, 
             col="black",fill=cancer_color , # #paletteer_d("ggthemes::Classic_Green_Orange_6",5)
             # filename=paste0(outdir,"/Venn_plot/Venn_all_URs.tiff"),
             filename = NULL,
             cex=1,cat.cex = 1,cat.fontface = "bold",alpha=0.8,
             # fontfamily=3,
             margin = 0.05)
pdf(paste0(outdir,"/Venn_plot/Venn_all_URs.pdf"))
grid.draw(venn.plot)
dev.off()

#export venn result
inter<-get.venn.partitions(x=venn_list)
venn_results<- data.frame()
for (i in 1:nrow(inter)) {
  venn_results[i,'set'] <- paste(inter[[i,'..set..']], collapse = ', ')
  venn_results[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
  venn_results[i,'counts'] <- paste(inter[[i,'..count..']], collapse = ', ')
}
write.table(venn_results,file =paste0(outdir,"/Venn_plot/venn_result_all_URs.txt"), sep ="\t" , append = F, row.names = T)


#2.Find DSs of shared Ligands ###########################
ligand_target_merged_allDS <- read.csv(paste0(outdir_allDS_of_allUR,"/ligand_target_merged_withFC.csv"), sep=",",header=T)
unique(ligand_target_merged_allDS$target) %>% length()
unique(ligand_target_merged_allDS$ligand) %>% length()
shared_ligands = read.csv(paste0(outdir,"/shared_UR.txt"),sep = "\t",header = F)
shared_ligands = shared_ligands$V1  %>% as.vector
  
length(shared_ligands)
ligand_target_merged_DS_of_shUR = ligand_target_merged_allDS[ligand_target_merged_allDS$ligand %in% shared_ligands,]
dim(ligand_target_merged_DS_of_shUR)
unique(ligand_target_merged_DS_of_shUR$target) %>% length()
write.table(ligand_target_merged_DS_of_shUR,file = paste0(outdir_shDS_of_shUR,"/ligand_target_merged_withFC.csv"),sep = ",") 

##2.3 shared DS of shared UR####
shared_DS <- intersect(ligand_target_merged_DS_of_shUR[ligand_target_merged_DS_of_shUR$cancer==cancer_name[1],"target"] %>% unique,
                       ligand_target_merged_DS_of_shUR[ligand_target_merged_DS_of_shUR$cancer==cancer_name[2],"target"] %>% unique) %>% 
  intersect(.,ligand_target_merged_DS_of_shUR[ligand_target_merged_DS_of_shUR$cancer==cancer_name[3],"target"] %>% unique) %>% 
  intersect(.,ligand_target_merged_DS_of_shUR[ligand_target_merged_DS_of_shUR$cancer==cancer_name[4],"target"] %>% unique) %>% 
  intersect(.,ligand_target_merged_DS_of_shUR[ligand_target_merged_DS_of_shUR$cancer==cancer_name[5],"target"] %>% unique)
length(shared_DS)                     
write.table(shared_DS, file=paste0(outdir,"/shared_DS_of_shared_UR.txt"), sep="\t", col.names = F, row.names = F)

##2.4 Venn plot for DSs of all shUR####
library(VennDiagram)
library(gplots)
library(RColorBrewer)

venn_list = list(Breast=ligand_target_merged_DS_of_shUR[ligand_target_merged_DS_of_shUR$cancer=="Breast","target"]  %>% unique,
                 Colon=ligand_target_merged_DS_of_shUR[ligand_target_merged_DS_of_shUR$cancer=="Colon","target"]  %>% unique,
                 Liver=ligand_target_merged_DS_of_shUR[ligand_target_merged_DS_of_shUR$cancer=="Liver","target"]  %>% unique,
                 Lung=ligand_target_merged_DS_of_shUR[ligand_target_merged_DS_of_shUR$cancer=="Lung","target"]  %>% unique,
                 Ovary=ligand_target_merged_DS_of_shUR[ligand_target_merged_DS_of_shUR$cancer=="Ovary","target"]  %>% unique)
venn.plot <- venn.diagram(x=venn_list, 
             col="black",fill=cancer_color, # #paletteer_d("ggthemes::Classic_Green_Orange_6",5)
             filename=NULL,
             cex=1,cat.cex = 1,cat.fontface = "bold",alpha=0.9,
             margin = 0.05)
pdf(paste0(paste0(wd,"/cancer/datafile/outputs_5cancers/new_outputs/MCDM/Venn_plot/Venn_DS_of_shared_UR.pdf")))
grid.draw(venn.plot)
dev.off()

#export venn result
inter<-get.venn.partitions(x=venn_list)
venn_results<- data.frame()
for (i in 1:nrow(inter)) {
  venn_results[i,'set'] <- paste(inter[[i,'..set..']], collapse = ', ')
  venn_results[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
  venn_results[i,'counts'] <- paste(inter[[i,'..count..']], collapse = ', ')
}
write.table(venn_results,file =paste0(wd,"/cancer/datafile/outputs_5cancers/new_outputs/MCDM/Venn_plot/Venn_result_DS_of_shared_UR.txt"), sep ="\t" , append = F, row.names = T)

##################################################################################
#3.scUR - Shared URs that has the same expression change in > 4 cancers ####
shared_ligands <- read.csv(paste0(outdir,"/shared_UR.txt"), sep="\t",header=F) %>% as.vector
all_ligands <- read.csv(paste0(outdir,"/all_positive_ligand.txt"), sep="\t",header=F) %>% as.vector

##3.1. import ligands_FC and merge####
allUR_FC <- list.files(paste0(outdir,"/UR_DS_FC/"),pattern = "*_allUR_FC.txt")
name <- sapply(strsplit(allUR_FC,"_"),'[[',1)

cancer = c()
shUR_FC_list <- list()
for (i in 1:length(allUR_FC)){
  a = read.table(file = paste0(outdir,"/UR_DS_FC/",allUR_FC[i]),header = T,row.names = 1)
  a = a[rownames(a) %in% shared_ligands$V1,] 
  cancer = c(cancer,rep(name[i],(dim(a)[2])))
  # colnames(a) = paste0(colnames(a),'_',name[i])
  a = a %>% list
  names(a) = name[i]
  shUR_FC_list = c(shUR_FC_list,a)
}

shUR_FC_merged <- cbind(shUR_FC_list[[1]],shUR_FC_list[[2]],shUR_FC_list[[3]],shUR_FC_list[[4]],shUR_FC_list[[5]])
shUR_FC_merged[is.na(shUR_FC_merged)] <- 0
# write.table(shUR_FC_merged, file = paste0(outdir,"/UR_DS_FC/shared_UR_FC_cell_cancer_merged.txt"), sep ="\t" , append = F, row.names = T,header=T)

write.table(shUR_FC_merged, file = paste0(outdir,"/UR_DS_FC/shared_UR_FC_merged.txt"), sep ="\t" , append = F, row.names = T)
dim(shUR_FC_merged)

#### 3.1.1 calculate mean FC in each cell type by cancer ####
a = t(shUR_FC_merged) 
cell = rownames(a)
a = a %>% as.data.frame()
a$cell = cell
a$cancer = cancer
a[1:20,-1:-5]

library(reshape2)
library(tidyr)
a.long <- melt(a, id.vars=c("cell","cancer"),
               variable.name="shUR",
               value.name="log2FC")
a.long$gene_cell = paste0(a.long$shUR,"_",a.long$cell)
a.long <- a.long %>% distinct()
a.wide = spread(a.long,key=cancer, value=log2FC) %>% as.data.frame()
a.wide[1:5,]
a.wide$mean_log2FC <- apply(a.wide[,4:8],1,function(x){mean(x,na.rm = TRUE)})
write.table(a.wide, file = paste0(outdir,"/UR_DS_FC/shared_UR_FC_merged_bycancer_mean.txt"), sep ="\t" , append = F, row.names = F)


##3.2 heatmap of shUR####
library(ComplexHeatmap)
library(circlize)
df <- shUR_FC_merged %>% as.matrix
Celltype = colnames(df)

unique(Celltype)
anno_col <- data.frame(Group=cancer)
anno_col$Group <- factor(anno_col$Group,levels=unique(anno_col$Group))

ha <- HeatmapAnnotation(Cancer_type = anno_col$Group,
                        #= c(rep("Lung",11),rep("Liver",11),rep("Colon_bel",12),rep("Breast_ER",10),rep("Ovary",13)),
                        # Celltype = Celltype
                        col = list(Cancer_type = c("Lung" = "green", "Liver" = "yellow", "Colon" = "blue", "Breast" = "darkred", "Ovary" = "red")))

dpi = 150
png(file=paste0(outdir,"/UR_DS_FC/heatmap_shared_UR_FC.png"), width = dpi*12, height = dpi*18, units = "px",res = dpi)
p = Heatmap(df,
            # split
            # row_split = 2, column_split = 3,  #split by dendrograms    
            # row_km = 2, column_km = 3 , #split rows and colums
            # split = data.frame(cyl = mtcars$cyl, am = mtcars$am), #根据数据进行分割
            column_split = anno_col$Group,
            #设置title,legend,font
            name = "Cancers", #title of legend
            column_title = "Cell type", 
            column_title_side = "bottom",  #设置列标题的位置，可选"top"或"bottom"
            column_title_gp = gpar(fontsize = 15, fontface = "bold"),  #更改列文本的字体
            column_names_gp = gpar(fontsize = 8),
            row_title = "Ligands",
            row_title_side ="right",  #设置行标题的位置，可选"left"或者"right"
            row_title_gp = gpar(fontsize = 15) ,
            row_names_gp = gpar(fontsize = 8), # Text size for row names
            # heatmap_legend_param = list(title="legend"),
            top_annotation = ha, #put annotation
            cluster_columns = F,
            # cluster_columns = cluster_within_group(df,anno_col$Group) #组内聚类
            #设置颜色
            col = colorRamp2(c(-4,0,4), c("blue", "white","red"))
            # col = colorpanel(50,"blue","white","red"),
            # col = colorRamp2(c(0, 10, 50), c("blue", "white", "red"))  #设置绘图的颜色，可以根据数据实际情况调整
            # na_col = "grey",
            # clustering_distance_rows = "pearson"
)
print(p)
dev.off()

##3.3 select scUR####
#1. change df to 0, -1, 1
df_log2fc_binary <- df 
# df_log2fc_binary[df_log2fc_binary<=1 & df_log2fc_binary>=-1] <- 0
df_log2fc_binary[df_log2fc_binary > 0] <- 1
df_log2fc_binary[df_log2fc_binary < 0] <- -1
df_log2fc_binary <- as.matrix(df_log2fc_binary,rownames.force = NA)
class(df_log2fc_binary)
unique(Celltype)

#2.identify shared ligands of each cell type
URs_same_FC_in_all_cancers = data.frame(matrix(NA, ncol=2))
names(URs_same_FC_in_all_cancers) = c("URs_same_fc","cell_type")
for (i in 1:length(unique(Celltype))){
  df_temp = df_log2fc_binary[,Celltype == unique(Celltype)[i]] %>% as.matrix()
  group_temp = anno_col$Group[Celltype == unique(Celltype)[i]]
  sum = apply(df_temp, 1, sum)
  URs_same_fc = rownames(df_temp[sum>=4|sum<=-4,])
  cell_type = rep(unique(Celltype)[i],length(URs_same_fc))
  a = data.frame(URs_same_fc,cell_type)
  print(a)
  URs_same_FC_in_all_cancers=rbind(URs_same_FC_in_all_cancers,a)
}

URs_same_FC_in_all_cancers <- URs_same_FC_in_all_cancers[2:dim(URs_same_FC_in_all_cancers)[1],]
scUR = URs_same_FC_in_all_cancers
a = table(URs_same_FC_in_all_cancers$URs_same_fc,URs_same_FC_in_all_cancers$cell_type) 
write.csv(scUR,paste0(outdir,"/scUR_celltype.csv"),row.names = F, col.names = T,quote = F)
write.csv(a,paste0(outdir,"/UR_DS_FC/No.URs_same_FC_in_all_cancers.csv"),row.names = T, col.names = T,quote = F)
write.csv(URs_same_FC_in_all_cancers,paste0(outdir,"/UR_DS_FC/URs_same_FC_in_all_cancers.csv"),row.names = T, col.names = T,quote = F)

##3.4 heatmap-of scUR####
Ligands_of_interests <- scUR$URs_same_fc %>% unique
df_log2fc_binary_candid <- df_log2fc_binary[rownames(df_log2fc_binary) %in% Ligands_of_interests,]

anno_col <- data.frame(Group=cancer)
anno_col$Group <- factor(anno_col$Group,levels=unique(anno_col$Group))

ha <- HeatmapAnnotation(Cancer_type = anno_col$Group,
                        #= c(rep("Lung",11),rep("Liver",11),rep("Colon_bel",12),rep("Breast_ER",10),rep("Ovary",13)),
                        # Celltype = Celltype
                        col = list(Cancer_type = c("Lung" = "green", "Liver" = "yellow", "Colon" = "blue", "Breast" = "darkred", "Ovary" = "red")))

dpi = 150
png(file=paste0(outdir,"/UR_DS_FC/heatmap_scUR_FC_binary_bycancer.png"), width = dpi*12, height = dpi*18, units = "px",res = dpi)
p=Heatmap(df_log2fc_binary_candid,
          # row_split = 2, column_split = 3,  #split by dendrograms    
          # row_km = 2, column_km = 3 , #split rows and colums
          # split = data.frame(cyl = mtcars$cyl, am = mtcars$am), #根据数据进行分割
          column_split = anno_col$Group, #anno_col$Group
          # row_split = 3,
          # row_dend_reorder = TRUE,
          #设置title,legend,font
          name = "Cancers", #title of legend
          column_title = "Cell type", 
          column_title_side = "bottom",  #设置列标题的位置，可选"top"或"bottom"
          column_title_gp = gpar(fontsize = 15, fontface = "bold"),  #更改列文本的字体
          column_names_gp = gpar(fontsize = 8),
          row_title = "Ligands",
          row_title_side ="right",  #设置行标题的位置，可选"left"或者"right"
          row_title_gp = gpar(fontsize = 15) ,
          row_names_gp = gpar(fontsize = 8), # Text size for row names
          # heatmap_legend_param = list(title="legend"),
          top_annotation = ha, #put annotation
          cluster_columns = T,
          cluster_rows = T,
          # cluster_rows = cluster_within_group(df,anno_col$Group), #组内聚类
          clustering_distance_rows = "euclidean", #"pearson","euclidean"
          clustering_method_rows  = "complete",
          #设置颜色
          border = T,
          col = colorRamp2(c(-1, 0,1), c("blue","white", "red"))
          # col = colorpanel(50,"blue","white","red"),
          # col = colorRamp2(c(0, 10, 50), c("blue", "white", "red"))  #设置绘图的颜色，可以根据数据实际情况调整
          # na_col = "grey"
)
print(p)
dev.off()

png(file=paste0(outdir,"/UR_DS_FC/heatmap_scUR_FC_binary_byCelltype.png"), width = dpi*12, height = dpi*18, units = "px",res = dpi)
p=Heatmap(df_log2fc_binary_candid,
          # row_split = 2, column_split = 3,  #split by dendrograms    
          # row_km = 2, column_km = 3 , #split rows and colums
          # split = data.frame(cyl = mtcars$cyl, am = mtcars$am), #根据数据进行分割
          column_split = anno_col$Celltype, #anno_col$Group
          # row_split = 11, #anno_col$Group
          # row_split = 2,
          # row_title = "cluster_%s",
          #设置title,legend,font
          name = "Cancers", #title of legend
          column_title = "Cell type", 
          column_title_side = "bottom",  #设置列标题的位置，可选"top"或"bottom"
          column_title_gp = gpar(fontsize = 15, fontface = "bold"),  #更改列文本的字体
          column_names_gp = gpar(fontsize = 8),
          row_title = "Ligands",
          row_title_side ="right",  #设置行标题的位置，可选"left"或者"right"
          row_title_gp = gpar(fontsize = 15) ,
          row_names_gp = gpar(fontsize = 8), # Text size for row names
          # heatmap_legend_param = list(title="legend"),
          top_annotation = ha, #put annotation
          cluster_columns = T,
          cluster_rows = T,
          border = T,
          # cluster_columns = cluster_within_group(df,anno_col$Group) #组内聚类
          # clustering_distance_rows = "euclidean", #"pearson","euclidean"
          # clustering_method_rows = "complete",
          #设置颜色
          col = colorRamp2(c(-1, 0,1), c("blue","white", "red"))
          # col = colorpanel(50,"blue","white","red"),
          # col = colorRamp2(c(0, 10, 50), c("blue", "white", "red"))  #设置绘图的颜色，可以根据数据实际情况调整
          # na_col = "grey",
)
print(p)
dev.off()

clustersCelltype = cutree(p$tree_col,k=15) %>% as.data.frame()   #cut tree_col
table(clustersCelltype$.)
head(clustersCelltype)

#4.scDS - downstreamers of scUR that have the same expression pattern in >4 cancers ###############
ligand_target_merged_allDS <- read.csv(paste0(outdir_allDS_of_allUR,"/ligand_target_merged_withFC.csv"), sep=",",header=T)
unique(ligand_target_merged_allDS$target) %>% length()
scUR <- read.csv(paste0(outdir,"/scUR_celltype.csv"),header = T)
scUR <- scUR$URs_same_fc %>% unique
length(scUR)
ligand_target_merged_DS_of_scUR = ligand_target_merged_allDS[ligand_target_merged_allDS$ligand %in% scUR,]
dim(ligand_target_merged_DS_of_scUR)
unique(ligand_target_merged_DS_of_scUR$target) %>% length()

##4.1 Venn plot of DS_of_scUR ####
library(VennDiagram)
library(gplots)
library(RColorBrewer)

colorset4 <- c("#e97d72","#b49e33","#53b54c","#49acf4","#ed739b")

venn_list = list(Ovary=ligand_target_merged_DS_of_scUR[ligand_target_merged_DS_of_scUR$cancer=="Ovary","target"]  %>% unique,
                 Breast=ligand_target_merged_DS_of_scUR[ligand_target_merged_DS_of_scUR$cancer=="Breast","target"]  %>% unique,
                 Colon=ligand_target_merged_DS_of_scUR[ligand_target_merged_DS_of_scUR$cancer=="Colon","target"]  %>% unique,
                 Liver=ligand_target_merged_DS_of_scUR[ligand_target_merged_DS_of_scUR$cancer=="Liver","target"]  %>% unique,
                 Lung=ligand_target_merged_DS_of_scUR[ligand_target_merged_DS_of_scUR$cancer=="Lung","target"]  %>% unique)
venn.diagram(x=venn_list, 
             col="white",fill=colorset4, # #paletteer_d("ggthemes::Classic_Green_Orange_6",5)
             filename=paste0(wd,"/cancer/datafile/outputs_5cancers/new_outputs/MCDM/Venn_plot/Venn_DS_of_scUR.png"),
             cex=1,cat.cex = 1,cat.fontface = "bold",alpha=0.9,
             margin = 0.05)
#export venn result
inter<-get.venn.partitions(x=venn_list)
venn_results<- data.frame()
for (i in 1:nrow(inter)) {
  venn_results[i,'set'] <- paste(inter[[i,'..set..']], collapse = ', ')
  venn_results[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
  venn_results[i,'counts'] <- paste(inter[[i,'..count..']], collapse = ', ')
}
write.table(venn_results,file =paste0(wd,"/cancer/datafile/outputs_5cancers/new_outputs/MCDM/Venn_plot/Venn_result_DS_of_scUR.txt"), sep ="\t" , append = F, row.names = T)

##4.2 select scDS based on the FC, same change in >4 cancers ####
#make it to wide format
library("tidyr")
df <- ligand_target_merged_DS_of_scUR
df <-df[,c(1,2,3,7)]
df <- df[!duplicated(df),] #same target of different ligands
head(df)
df <- spread(df,                                  
             key = cancer,
             value = avg_log2FC)
head(df)
df[is.na(df)] <- 0

#choose scDS based on FC
targets_shared_4cancers <- data.frame(matrix(NA,ncol=2))
colnames(targets_shared_4cancers) = c("shared_target","receiver_cell")

for (i in 1:length(unique(df$receiver))){
  df_temp = df[df$receiver == unique(df$receiver)[i],]
  print(paste0("Targets(DSs) in ",unique(df$receiver)[i]))
  
  df_temp = t(df_temp)
  df_temp[1:4,1:2]
  colnames(df_temp) = df_temp[1,]
  df_temp = df_temp[3:dim(df_temp)[1],]
  rowname = rownames(df_temp)
  df_temp = apply(df_temp,2,as.numeric)
  rownames(df_temp) = rowname
  
  df_temp[1:5,1:3]
  df_temp %>% max()
  
  df_temp_binary <- df_temp 
  df_temp_binary[is.na(df_temp_binary)] <- 0
  # df_temp_binary[df_temp_binary<=1 & df_temp_binary>=-1] <- 0
  df_temp_binary[df_temp_binary > 0] <- 1
  df_temp_binary[df_temp_binary < 0] <- -1
  df_temp_binary <- as.matrix(df_temp_binary,rownames.force = NA)
  class(df_temp_binary)
  df_temp_binary[1:5,1:5]
  
  ## choose genes have the same expression change in >4 cancers
  sum = apply(df_temp_binary,2,sum)
  if (length(sum[sum>=4|sum<=-4]) ==0){ next }
  if (length(sum[sum>=4|sum<=-4]) ==1){
    df_temp_select <- matrix(df_temp_binary[,sum>=4|sum<=-4],)
    rownames(df_temp_select) <- rownames(df_temp_binary)
    colnames(df_temp_select) <- names(sum[sum>=4|sum<=-4])
  } else {
    df_temp_select <- df_temp_binary[,sum>=4|sum<=-4]
  }
  
  a <- data.frame(colnames(df_temp_select),unique(df$receiver)[i])
  colnames(a) <- c("shared_target","receiver_cell")
  targets_shared_4cancers <- rbind(targets_shared_4cancers,a)
  
  # dpi = 300
  # png(file=paste0(outdir,"/UR_DS_FC/heatmap_ligand_target_",unique(df$receiver)[i],"_logFC.png"), width = dpi*10, height = dpi*5, units = "px",res = dpi)
  # p = Heatmap(df_temp_select,
  #             # clustering_distance_columns = "euclidean", #specify distance as a pre-defined option. The valid values are the supported methods in dist() function and in "pearson", "spearman" and "kendall"
  #             clustering_method_columns = "complete", # one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).)
  #             # split
  #             # column_split = 2,  #split by dendrograms
  #             # row_km = 2, column_km = 3 , #split rows and colums
  #             # split = data.frame(cyl = mtcars$cyl, am = mtcars$am), #根据数据进行分割
  #             # column_split = anno_col$Group,
  #             #设置title,legend,font
  #             name = "log2FC", #title of legend
  #             column_title = paste0("Down stream targets in ",unique(df$receiver)[i]),
  #             column_title_side = "bottom",  #设置列标题的位置，可选"top"或"bottom"
  #             column_title_gp = gpar(fontsize = 15, fontface = "bold"),  #更改列文本的字体
  #             column_names_gp = gpar(fontsize = 6),
  #             row_title = "Cancers",
  #             row_title_side ="right",  #设置行标题的位置，可选"left"或者"right"
  #             row_title_gp = gpar(fontsize = 15) ,
  #             row_names_gp = gpar(fontsize = 10), # Text size for row names
  #             # heatmap_legend_param = list(title="legend"),
  #             # top_annotation = ha, #put annotation
  #             cluster_columns = T,
  #             cluster_rows = T,
  #             # cluster_columns = cluster_within_group(df,anno_col$Target_cell_type),  #组内聚类
  #             #设置颜色
  #             col = colorRamp2(c(-1,0,1), c("blue","white","red"))
  #             # col = colorpanel(50,"blue","white","red"),
  #             # col = colorRamp2(c(0, 10, 50), c("blue", "white", "dark red"))  #设置绘图的颜色，可以根据数据实际情况调整
  #             # na_col = "grey",
  #             # clustering_distance_rows = "pearson"
  # )
  # print(p)
  # dev.off()
}

targets_shared_4cancers <- targets_shared_4cancers[2:dim(targets_shared_4cancers)[1],]
write.csv(targets_shared_4cancers,file = paste0(outdir,"/scDS_celltype.csv"),row.names = T)
scDS = unique(targets_shared_4cancers$shared_target) 
write.csv(scDS,paste0(outdir,"/scDS.csv"),row.names = F, col.names = F,quote = F)

#5. Identify MCDM####
scUR <- read.csv(paste0(outdir,"/scUR_celltype.csv"),header = T)
scUR$ligands_sender <- paste0(scUR$URs_same_fc,"_",scUR$cell_type)
scUR$URs_same_fc %>% unique %>% length
# targets_shared_4cancers = read.csv(paste0(outdir,"/scDS_celltype.csv"),header=T)

targets_shared_4cancers$target_receiver <- paste0(targets_shared_4cancers$shared_target,"_",targets_shared_4cancers$receiver_cell) 
targets_shared_4cancers %>% head
ligand_target_merged_DS_of_scUR$target_receiver  <- paste0(ligand_target_merged_DS_of_scUR$target,"_",ligand_target_merged_DS_of_scUR$receiver) 
ligand_target_merged_DS_of_scUR %>% head
ligand_target_merged_DS_of_scUR %>% dim

a <- ligand_target_merged_DS_of_scUR[ligand_target_merged_DS_of_scUR$target_receiver %in% unique(targets_shared_4cancers$target_receiver),] # extract DS_celltype info of all DS of scURs
b <- table(scUR$URs_same_fc,scUR$cell_type) %>% as.data.frame.array()
b$ligand = rownames(b)
c.wide <- merge(a,b,by="ligand",all.y=T) #add information of scUR_celltype for each DS
dim(a)
head(c.wide)
library(reshape2)
c.long <- melt(c.wide, id.vars=c(colnames(c.wide)[1:8]), #change from wide to long
               variable.name="sender",
               value.name="is_sender"
)
head(c.long) 
dim(c.long)
d <- c.long[c.long$is_sender==1,] #only keep the interactions from shUR_celltype
d <- d[,c(1,4,5,9)]
d<- d[!duplicated(d),]  
shared_MCDM <- d
head(shared_MCDM)
dim(shared_MCDM)
shared_MCDM$sender %>% table
write.table(shared_MCDM, file = paste0(outdir,"/shared_MCDM.txt"), sep ="\t" ,  row.names = F, quote=F)

e <- d %>% group_by_at(vars(ligand,receiver,sender)) %>%  
  summarize_all(paste, collapse=",")
e$no.targets <- apply(e[,4],1,function(x){
  strsplit(as.character(x),",")[[1]] %>% length})
MCDM_file <- as.data.frame(e)
head(MCDM_file)
dim(MCDM_file)
MCDM_file$ligand %>% unique %>% length

write.table(MCDM_file, file = paste0(outdir,"/shared_MCDM_for_cytoscape.txt"), sep ="\t" ,  row.names = F, quote=F)

##make file for pathway analysis - shared interactions####
list_for_IPA <- c()
names <- c()
for (i in 1:length(e$sender %>% unique)){
  for (j in 1:length(e$receiver %>% unique)){
    sub = e[(e$sender==(e$sender %>% unique)[i]) & (e$receiver==(e$receiver %>% unique)[j]),]
    if(dim(sub)[1]==0){next}
    subURgenes = sub[,1]  %>% unique %>% as.vector
    subDSgenes = apply(sub[,4],1,function(x){strsplit(as.character(x),",")[[1]]}) %>% unique 
    all = c(subURgenes,subDSgenes) %>% unique %>% unlist 
    list_for_IPA = c(list_for_IPA,list(all) )
    names <- c(names,paste0((e$sender %>% unique)[i]," to ",(e$receiver %>% unique)[j])) 
  }
}
names(list_for_IPA) <- names

list_for_IPA <- do.call(cbind,         #merge vectors with different length
                        lapply(lapply(list_for_IPA, unlist), `length<-`,
                               max(lengths(list_for_IPA))))
list_for_IPA[is.na(list_for_IPA)] = "NA"

write.table(list_for_IPA, file = paste0(outdir,"/shared_MCDM_for_pathway.txt"), sep ="\t" ,  row.names = F, quote=F)

#cytoscape####
library(RCy3)
cytoscapePing()
#node features
setNodeShapeDefault ('ELLIPSE') #ELLIPSE,OCTAGON #change node shape
setNodeSizeDefault(15)
setNodeFontSizeDefault (15)
lockNodeDimensions(TRUE)

nodedata <- getTableColumns("node")
edgedata <- getTableColumns("edge")
column <- 'name'
# values <- c("Neurons","B","Plasma","CD4T","CD4T.gammadeltaT","CD8T","CD8T.GammadeltaT","CD8T.NK","NK","Mast","DC","Monocyte","Macrophage","Endothelial","Pericyte","Fibroblast","Epithelial")
# colors <- c("#4F7CBA","#1BA3C6", "#2CB5C0","#21B087", "#33A65C", "#57A337", "#A2B627", "#D5BB21", "#F8B620", "#F89217", "#F06719", "#E03426", "#F64971", "#FC719E", "#EB73B3", "#CE69BE", "#A26DC2" ) #brewer.pal(11, "RdYlBu")
# values <- nodedata$name
# colors <- brewer.pal(n = 7, name = "Set3")
values <- c("Neurons","B","Plasma","CD4T","CD4T.gammadeltaT","CD8T","CD8T.GammadeltaT","CD8T.NK","NK","Mast","DC",
            "NKT","Monocyte","Macrophage","Endothelial","Pericyte","Fibroblast","Epithelial")
colors <- c("#d0e562","#ff579f", "#ff62c3","#f663e3", "#33A65C", "#db73fc", "#A2B627", "#D5BB21", "#728ff0", "#39b24b", "#4ba8db", 
            "#ae88ff","#3ac4d0", "#38c094", "#f4a666", "#93c05f", "#eed163", "#f9766d" ) #brewer.pal(11, "RdYlBu")

setNodeColorMapping (column, values, colors , mapping.type = "d") 

setNodeSizeMapping(table.column = "no.ligands",
                   table.column.values = c(0, 15),
                   sizes = c(20,50),
                   mapping.type = "c")

# setNodeCustomRadialGradient(
#   colors = c("#DDDDDD", "#DDDDDD"),
#   anchors = c(0, 1),
#   xCenter = 0.5,
#   yCenter = 0.5,
#   slot = 1,
#   style.name = NULL,
#   base.url = .defaultBaseUrl
# )


#edges
colnames(edgedata)
setEdgeLineWidthMapping(
  table.column='no.targets',
  table.column.values =c(1, mean(as.numeric(edgedata$no.targets)), max(edgedata$no.targets)),
  # widths = c(0.2,0.9,1.8),
  widths = c(0.01,0.5,1),
  mapping.type = "c",
  default.width = NULL,
  style.name = NULL,
  network = NULL)

setEdgeColorMapping(
  table.column='edge_color',
  table.column.values = values,
  colors = colors,
  mapping.type = "d",
  default.color = NULL,
  style.name = NULL,
  network = NULL)

setEdgeLabelMapping(
  table.column="ligand",
  style.name = NULL
)

setEdgeLabelMapping('ligand')
matchArrowColorToEdge(TRUE)

exportImage("/Users/cynthia_ye/Documents/OneDrive - Linköpings universitet/cancer/datafile/outputs_5cancers/new_outputs/MCDM/sharedMCDM2",'pdf', zoom=300) #.png scaled by 200%





#################
DS_plau = ligand_target_merged_allDS[ligand_target_merged_allDS$ligand == "PLAU",]
DS_plau = DS_plau[DS_plau$receiver == "Macrophage",]
DS_plau = unique(DS_plau$target)
write.table(DS_plau,file="/Users/cynthia_ye/Documents/OneDrive - Linköpings universitet/cancer/datafile/outputs_5cancers/new_outputs/Pathway_analysis/clusterProfiler/DS_PLAU_Macrop/macrophage_DS_plau.txt",quote=F)

DS_PLAUR = read.table(paste0(wd,"/cancer/datafile/outputs_5cancers/new_outputs/Pathway_analysis/clusterProfiler/DS_PLAU_Macrop/downstreams of PLAUR.txt"),sep="\t",header=T)
DS_PLAUR = DS_PLAUR$Symbol

DS_PLAU_PLAUR = intersect(DS_plau,DS_PLAUR)


paste(DS_PLAU_PLAUR,collapse=",")











