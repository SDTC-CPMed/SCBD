#When it says: Warning: namespace ‘Seurat’ is not available and has been replaced by .GlobalEnv when processing object ‘p1’
myPaths <- .libPaths() 
myPaths
# myPaths <- c(myPaths, "/home/yelzh67/R/x86_64-pc-linux-gnu-library/4.0") # add new path
myPaths <- c(myPaths[2], myPaths[1])  # switch them
.libPaths(myPaths)  

rm(list = ls())
library(tidyverse)
library(nichenetr)
library(Seurat) # Seurat V4
library(MAST)
library(patchwork)
library(ggplot2)
library(SingleR)
library(SummarizedExperiment)
# sessionInfo()

#mac or omika
location = "mac"

cancer_name = c("Lung","Liver","Colon","Ovary","Breast")

if (location == "mac"){
  cancer_seuratobj_list = c("/Users/cynthia_ye/Documents/OneDrive/5. Linkoping University/2.5 Cancer/1. GSE123902 LA/datafile/DEGs/Aug16_seurat_int/all_integrated_annotated.RData",
                            "/Users/cynthia_ye/Documents/OneDrive/5. Linkoping University/2.5 Cancer/2. GSE138709 ICC/datafile/Seurat_int/Aug20/OneDrive_2_23-08-2021/all_integrated_annotated.RData",
                            "/Users/cynthia_ye/Documents/OneDrive/5. Linkoping University/2.5 Cancer/5. GSE144735 CRC_belg/datafile/yl_seurat_DEGs/Aug17_seurat_inte_DEGs/OneDrive_1_24-08-2021/all.integrated_annotated.RData",
                            "/Users/cynthia_ye/Documents/OneDrive/5. Linkoping University/2.5 Cancer/6. PMID32561858_Ovary_C/datafile/DEGs/DEGs_MAST_Aug22/OneDrive_5_24-08-2021/all.integrated_annotated.RData",
                            "/Users/cynthia_ye/Documents/OneDrive/5. Linkoping University/2.5 Cancer/4. GSE161529 BC /datafile/DEGs_seurat/Aug17/all.integrated_annotated.RData"
  )
  # outdir = "/Users/cynthia_ye/Documents/OneDrive - Linköpings universitet/cancer/datafile/outputs_5cancers/new_outputs/Seurat_subset_celltyping/2022Feb12"
  outdir = "/Users/cynthia_ye/Documents/OneDrive/cancer/datafile/outputs_5cancers/new_outputs/Seurat_subset_celltyping/2022Feb12"
  dir.create(outdir)
}

if (location == "omika"){
  cancer_seuratobj_list = c("/home/yelzh67/Projects/cancer/datafile/GSE123902_LA/outputs/Seurat_int/Aug16/all_integrated_annotated.RData",
                            "/home/yelzh67/Projects/cancer/datafile/GSE138709_ICC/outputs/Seurat_int/Aug19/all_integrated_annotated.RData",
                            "/home/yelzh67/Projects/cancer/datafile/GSE144735_CRC_Belgian/outputs/seurat_intergra/Jul17_res1/all.integrated_annotated.RData",
                            "/home/yelzh67/Projects/cancer/datafile/PMID32561858_Ovary_C/outputs/Seurat_int/Aug21/all.integrated_annotated.RData",
                            "/home/yelzh67/Projects/cancer/datafile/GSE161529_BC/outputs/seurat_intergra/Jul20_ER_N_qcncount500_hv2000/all.integrated_annotated.RData"
  )
  outdir = "/home/yelzh67/Projects/cancer/datafile/outputs_5cancers/new_outputs/Seurat_subset_celltyping/2022Feb12"
  dir.create(outdir)
}

outdir_PCA <- paste0(outdir,"/PCA")
if (dir.exists(outdir_PCA)==F){dir.create(outdir_PCA)}

outdir_singleR=paste0(outdir,"/singleR")
if (dir.exists(outdir_singleR) != 1) {  dir.create(outdir_singleR)}

outdir_markers=paste0(outdir,"/markers")
if (dir.exists(outdir_markers) != 1) { dir.create(outdir_markers)}

outdir_subset = paste0(outdir,"_subset")
if (dir.exists(outdir_subset)!=T){dir.create(outdir_subset)}

outdir_subset_PCA = paste0(outdir,"_subset/subset_PCA")
if (dir.exists(outdir_subset_PCA)!=T){dir.create(outdir_subset_PCA)}

outdir_subset_singleR = paste0(outdir,"_subset/subset_singleR")
if (dir.exists(outdir_subset_singleR)!=T){dir.create(outdir_subset_singleR)}

outdir_subset_markers = paste0(outdir,"_subset/subset_markers")
if (dir.exists(outdir_subset_markers)!=T){dir.create(outdir_subset_markers)}

######################################################################################## 
#1. subset celltype from each cancers ####
Stromal_sublist = list()
Lymphocytes_sublist = list()
Myeloids_sublist = list()
Epithelial_sublist = list()

outdir = "/home/yelzh67/Projects/cancer/datafile/outputs_5cancers/new_outputs/Seurat_subset_celltyping/2022Feb12"
file_list= list.files(path=outdir, pattern="all_integrated_singleR_main_")
cancer_name= sapply(strsplit(file_list,".RData"),'[[',1)
cancer_name= sapply(strsplit(cancer_name,"_"),'[[',5)

for (i in 1:length(file_list)) {
  load(paste0(outdir,"/",file_list[i]))
  print(i)
  DefaultAssay(all.integrated) <- "RNA"
  seuratObj = all.integrated
  seuratObj@meta.data %>% head()
  
  seuratObj@meta.data$Celltype_main %>% table()
  seuratObj@meta.data$cancer_name = cancer_name[i]
  
  Stromal_sub <- subset(x = seuratObj, 
                        subset = Celltype_main == "Stromal")
  Epi_sub <- subset(x = seuratObj, 
                    subset = Celltype_main == "Epithelial")
  Myeloids_sub <- subset(x = seuratObj, 
                         subset = Celltype_main == "Myeloids")  
  Lymphocytes_sub <- subset(x = seuratObj, 
                            subset = Celltype_main == "Lymphocytes")  
  
  Stromal_sublist = c(Stromal_sublist,Stromal_sub)
  Lymphocytes_sublist = c(Lymphocytes_sublist,Lymphocytes_sub)
  Myeloids_sublist = c(Myeloids_sublist,Myeloids_sub)
  Epithelial_sublist = c(Epithelial_sublist,Epi_sub)
  
  print(paste0(cancer_name[i], " done!"))
  rm(all.integrated,seuratObj)
}

save(Stromal_sublist,file = paste0(outdir_subset,"/Stromal_sublist.RData"))
save(Lymphocytes_sublist,file = paste0(outdir_subset,"/Lymphocytes_sublist.RData"))
save(Myeloids_sublist,file = paste0(outdir_subset,"/Myeloids_sublist.RData"))
save(Epithelial_sublist,file = paste0(outdir_subset,"/Epithelial_sublist.RData"))
Epithelial_sublist[[1]]$Celltype_main %>% table
Stromal_sublist[[1]]@meta.data %>% head()

######################################################################################## 
#2. merge the same celltype across cancers - seurat integration####
# library(SeuratDisk)
# library(SeuratWrappers)
##2.1 select integration features####
Stromal_sublist <- lapply(X = Stromal_sublist, FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = Stromal_sublist)
save(features,file = paste0(outdir_subset,"/features_Stromal.RData"))

Epithelial_sublist <- lapply(X = Epithelial_sublist, FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = Epithelial_sublist)
save(features,file = paste0(outdir_subset,"/features_Epithelial.RData"))

Myeloids_sublist <- lapply(X = Myeloids_sublist, FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = Myeloids_sublist)
save(features,file = paste0(outdir_subset,"/features_Myeloids.RData"))

Lymphocytes_sublist <- lapply(X = Lymphocytes_sublist, FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = Lymphocytes_sublist)
save(features,file = paste0(outdir_subset,"/features_Lymphocytes.RData"))

##2.2 integration  ####
load(paste0(outdir_subset,"/features_Lymphocytes.RData"))
anchors <- FindIntegrationAnchors(object.list = Lymphocytes_sublist, anchor.features = features)
all.integrated <- IntegrateData(anchorset = anchors)
save(all.integrated,file = paste0(outdir_subset,"/all_integrated_Lymphocytes.RData"))
rm(features)

load(paste0(outdir_subset,"/features_Myeloids.RData"))
anchors <- FindIntegrationAnchors(object.list = Myeloids_sublist, anchor.features = features)
all.integrated <- IntegrateData(anchorset = anchors)
save(all.integrated,file = paste0(outdir_subset,"/all_integrated_Myeloids.RData"))
rm(features)

load(paste0(outdir_subset,"/features_Stromal.RData"))
anchors <- FindIntegrationAnchors(object.list = Stromal_sublist, anchor.features = features)
all.integrated <- IntegrateData(anchorset = anchors)
save(all.integrated,file = paste0(outdir_subset,"/all_integrated_Stromal.RData"))
rm(features)

load(paste0(outdir_subset,"/features_Epithelial.RData"))
anchors <- FindIntegrationAnchors(object.list = Epithelial_sublist, anchor.features = features)
all.integrated <- IntegrateData(anchorset = anchors)
save(all.integrated,file = paste0(outdir_subset,"/all_integrated_Epithelial.RData"))
rm(features)

##2.3 Clustering after integration####
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
for (i in c("Myeloids","Lymphocytes","Stromal","Epithelial")){
  # i="Myeloids"
  load(paste0(outdir_subset,"/all_integrated_",i,".RData"))
  all.integrated@meta.data <- all.integrated@meta.data[colnames(all.integrated@assays$RNA@counts),]
  DefaultAssay(all.integrated) <- "integrated"
  seed = 42
  nPCs = 30
  res= 0.8
  reduction= "umap"
  
  all.integrated <- ScaleData(all.integrated, verbose = T)
  all.integrated <- RunPCA(all.integrated,seed.use = seed, verbose = T)
  p <- ElbowPlot(all.integrated,ndims = 30)
  ggsave(plot=p,  path = outdir_subset_PCA, dpi = 300, width = 25, height = 20, units = "cm" ,filename=paste0(i,"_PCA-ElbowPlot.pdf"))
  all.integrated <- RunUMAP(all.integrated, seed.use = seed, reduction = "pca", dims = seq(nPCs))
  all.integrated <- RunTSNE(all.integrated, seed.use = seed, reduction = "pca", dims=seq(nPCs))
  all.integrated <- FindNeighbors(all.integrated, reduction = "pca", dims = seq(nPCs))
  all.integrated <- FindClusters(all.integrated, random.seed = seed, resolution = res)
  
  # save(all.integrated, file = paste0(outdir_subset,"/all_integrated_",i,".RData"))
  
  ##2.4 singleR  ####
  sce <- all.integrated
  sce@meta.data = sce@meta.data[colnames(sce@assays$RNA@counts),]
  # rm(all.integrated)
  sce_for_singleR <- GetAssayData(sce, slot="data") #或者从seuratobject读取
  clusters=sce@meta.data$seurat_clusters
  head(sce@meta.data)
  table(sce$celltype_by_marker,sce$Celltype_main)
  table(sce$celltype_by_marker,sce$cancer_name)
  
  # hpca.se=HumanPrimaryCellAtlasData() #SingleR自带的reference，多为bulk数据或者microarray数据
  # save(hpca.se,file=paste0(outdir_singleR,"/hpca.se.RData"))
  # load(paste0(outdir_singleR,"/hpca.se.RData"))
  common_hpca <- intersect(rownames(sce), rownames(hpca.se))
  
  #小类注释
  pred.fine.hpca <- SingleR(test = sce_for_singleR, ref = hpca.se, labels = hpca.se$label.fine,
                            clusters = clusters, 
                            assay.type.test = "logcounts", assay.type.ref = "logcounts")
  
  cellType=data.frame(ClusterID=levels(as.factor(sce@meta.data$seurat_clusters)),
                      fine.hpca=pred.fine.hpca$labels)
  
  write.table(cellType, file =paste0(outdir_subset_singleR, "/cellType.SingleR.fine","_",i,"_",reduction,"_seed",seed,"_res",res,".txt"), sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE) #保存下来，方便以后调用
  
  #4 sce@meta.data中，然后画tsne/umap展示一下
  sce@meta.data$CB=rownames(sce@meta.data)
  sce@meta.data=merge(sce@meta.data,cellType,by.x="seurat_clusters",by.y="ClusterID")
  rownames(sce@meta.data)=sce@meta.data$CB
  
  all.integrated@meta.data <- sce@meta.data
  DefaultAssay(all.integrated) <- "integrated"
  p5 <- DimPlot(all.integrated, reduction = reduction, group.by = "main.hpca.y", pt.size=0.5)+theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),axis.text = element_blank())
  p6 <- DimPlot(all.integrated, reduction = reduction, group.by = "ident",   pt.size=0.5, label = TRUE,repel = TRUE)+theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),axis.text = element_blank())
  # fig_tsne <- plot_grid(p6, p5, labels = c('ident','HPCA_Main'),rel_widths = c(2,3))
  p7 <- DimPlot(all.integrated, reduction = reduction, group.by = "cancer_name",   pt.size=0.5, label = TRUE,repel = TRUE)+theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),axis.text = element_blank())
  fig <- p5 + p6 + p7
  ggsave(filename = paste0(outdir_subset_singleR, "/singleR_fine","_",i,"_",reduction,"_seed",seed,"_PC",nPCs,"_res",res,".pdf"), plot = fig, device = 'pdf', width = 48, height = 12, units = 'cm')
  
  p <- DimPlot(all.integrated, reduction = reduction, group.by = "seurat_clusters", label = TRUE)
  ggsave(plot=p,  path = outdir_subset_PCA, dpi = 300, width = 45, height = 40, units = "cm" ,filename=paste0("DimPlot_by_celltype","_",i,"_",reduction,"_seed",seed,"_PC",nPCs,"_res",res,".png"))
  
  DefaultAssay(all.integrated) <- "RNA"
  p <- FeaturePlot(all.integrated, features = c('PLAU'),split.by = "group", max.cutoff = 3, label = TRUE, cols = c("grey", "red"))
  ggsave(plot=p,  path = outdir_subset_PCA, dpi = 300, width = 45, height = 40, units = "cm" ,filename=paste0("VlnPlot","_",i,"_",reduction,"_seed",seed,"_PC",nPCs,"_res",res,".png"))
  
  p <- DotPlot(all.integrated, features = unique(genes_to_check),
               assay='RNA' , ##can use this to check orignial or normalized expression
               group.by = 'seurat_clusters' ##this can change to other group method like 'seurat_clusters'
  )  + scale_color_viridis_c() + coord_flip() 
  ggsave(plot=p,  path = outdir_subset_PCA, dpi = 300, width = 45, height = 40, units = "cm" ,filename=paste0("check_markers_",i,"_seed",seed,"_PC",nPCs,"_res",res,".pdf"))
  
  DefaultAssay(all.integrated) <- "RNA"
  p <- FeaturePlot(all.integrated, features = unique(genes_to_check), max.cutoff = 3, label = TRUE,
                   cols = c("grey", "red"))
  ggsave(plot=p,  path = outdir_subset_PCA, device = "png",dpi = 150, width = 18, height = 45, units = "in" ,filename=paste0("Feature_Plot","_",i,"_seed",seed,"_PC",nPCs,"_res",res,".png"))
  
  all.integrated@meta.data = all.integrated@meta.data[colnames(all.integrated@assays$RNA@counts),]
  
  table(all.integrated$celltype_by_marker,all.integrated$seurat_clusters)
  table(all.integrated$seurat_clusters,all.integrated$integrated_snn_res.0.8)
  # table(sce$celltype_by_marker,sce$integrated_snn_res.0.5)
  
  save(all.integrated,file= paste0(outdir_subset,"/all_integrated_",i,".RData"))
}

## 2.5 find marker genes ####
for (i in c("Myeloids","Lymphocytes","Stromal","Epithelial")){ #findallmarkers
  load(paste0(outdir_subset,"/all_integrated_",i,".RData"))
  print(i)
  DefaultAssay(all.integrated) <- "RNA"
  all.integrated@assays$RNA@counts[1,1:5]
  
  all.integrated@meta.data = all.integrated@meta.data[colnames(all.integrated@assays$RNA@counts),]
  Idents(all.integrated) <- all.integrated$seurat_clusters
  head(Idents(all.integrated))
  
  #这一步计算的时候可以把min.pct以及logfc.threshold调的比较低，然后再基于结果手动筛选
  all.markers <- FindAllMarkers(all.integrated, only.pos = TRUE, 
                                min.pct = 0.1, logfc.threshold = 0.25)
  #把每个细胞群都计算一遍
  write.table(all.markers,
              file=paste0(outdir_subset_markers,"/",i,"_total_marker_genes_",nPCs,"PC.txt"),
              sep="\t",quote = F,row.names = F)
  
  # 遍历每一个cluster然后展示其中前10个基因
  marker.sig <- all.markers %>% 
    mutate(Ratio = round(pct.1/pct.2,3)) %>%
    filter(p_val_adj <= 0.05)  # 本条件为过滤统计学不显著的基因
  
  for(cluster_id in unique(marker.sig$cluster)){
    # cluster.markers <- FindMarkers(experiment.aggregate, ident.1 = cluster, min.pct = 0.3)
    # cluster.markers <- as.data.frame(cluster.markers) %>%
    #   mutate(Gene = rownames(cluster.markers))
    cl4.genes <- marker.sig %>%
      filter(cluster == cluster_id) %>%
      arrange(desc(avg_log2FC))
    cl4.genes <- cl4.genes[1:min(nrow(cl4.genes),10),"gene"]
    
    # #VlnPlot
    # pvn <- VlnPlot(all.integrated, features = cl4.genes,ncol = 5)
    # pdf(paste0(outdir_subset_marker,"/",i,'_',cluster_id,"_vlnplot_",nPCs,"PC.pdf"),width = 30,height = 15)
    # print(pvn)
    # dev.off()
    # 
    # #feature plot 
    # pvn <- FeaturePlot(all.integrated,features=cl4.genes,ncol = 5,label = TRUE,)
    # pdf(paste0(outdir_subset_marker,"/",i,'_',cluster_id,"_feature_tsne_",nPCs,"PC.pdf"),width = 30,height = 15)
    # print(pvn)
    # dev.off()
    # 
    #RidgePlot
    # pvn<-RidgePlot(all.integrated, features = cl4.genes, ncol = 5)
    # pdf(paste0(outdir_subset_marker,"/",i,'_',cluster_id,"_ridge_tsne_",nPCs,"PC.pdf"),width = 30, height = 15)
    # print(pvn)
    # dev.off()
  }
  # rm(cl4.genes,cluster_id,pvn)
  
  #热图展示Top marker基因
  #筛选top的marker基因，可以通过参数改为其他数值
  top <- marker.sig %>% group_by(cluster) %>% 
    top_n(n = 10, wt = avg_log2FC)
  
  #top-marker基因dotplot
  pdf(paste0(outdir_subset_markers,"/MarkerGene-DotPlot_all_cluster_",i,"_","res",res,"_PC",nPCs,"PC.pdf"),width = 50,height = 10)
  DotPlot(all.integrated, features = unique(top$gene))+scale_color_viridis_c()+RotatedAxis()
  dev.off()
  
  #top-marker基因热图
  DefaultAssay(all.integrated)  <- 'integrated'
  pdf(paste0(outdir_subset_markers,"/MarkerGene-Heatmap_all_cluster_",i,"_","res",res,"_PC",nPCs,"PC.pdf"),width= 10, height= 15 )
  DoHeatmap(all.integrated, features = top$gene,size = 2) +
    theme(legend.position = "none", 
          axis.text.y = element_text(size = 5))
  dev.off()
  
  p <- DimPlot(all.integrated, reduction = reduction, group.by = "group",   pt.size=0.5, label = TRUE,repel = TRUE)+theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),axis.text = element_blank())
  ggsave(plot=p,  path = outdir_subset_PCA, dpi = 300, width = 36, height = 24, units = "cm" ,filename=paste0("/DimPlot_by_group_",i,"_",reduction,"_seed",seed,"_res",res,".pdf"))
}

##2.6 re-identify celltype####
for (i in c("Myeloids","Lymphocytes","Stromal","Epithelial")){
  load(paste0(outdir_subset,"/all_integrated_",i,".RData"))
  print(i)
  DefaultAssay(all.integrated) <- "integrated"
  all.integrated@assays$RNA@counts[1,1:5]
  
  all.integrated@meta.data = all.integrated@meta.data[colnames(all.integrated@assays$RNA@counts),]
  Idents(all.integrated) <- all.integrated$seurat_clusters
  head(Idents(all.integrated))
  
  if (i=="Myeloids"){ #done
    all.integrated <- RenameIdents(all.integrated, 
                                   `0` = "Macrophage",# 
                                   `1` = "Macrophage", #
                                   `2` = "Macrophage",
                                   `3` = "Macrophage", 
                                   `4` = "DC", 
                                   `5` = "Macrophage", 
                                   `6` = "Macrophage", 
                                   `7` = "Macrophage", 
                                   `8` = "Macrophage",
                                   `9` = "Macrophage",
                                   `10` = "Monocyte",
                                   `11` = "Endothelial", 
                                   `12` = "Mast", 
                                   `13` = "NKT",#CD8T?
                                   `14` = "Monocyte", 
                                   `15` = "B", 
                                   `16` = "Macrophage", #?epi
                                   `17` = "Macrophage",
                                   `18` = "Macrophage", 
                                   `19` = "DC",
                                   `20` = "DC", 
                                   `21` = "DC", 
                                   `22` = "Macrophage",
                                   `23` = "Fibroblast",  #?
                                   `24` = "DC"
    )}
  
  if (i=="Lymphocytes"){#done
    all.integrated <- RenameIdents(all.integrated, 
                                   `0` = "CD8T", 
                                   `1` = "CD4T", 
                                   `2` = "CD4T",
                                   `3` = "NKT",  
                                   `4` = "CD4T", 
                                   `5` = "B", 
                                   `6` = "CD8T", 
                                   `7` = "Plasma",
                                   `8` = "CD4T", 
                                   `9` = "NK",
                                   `10` = "CD8T", 
                                   `11` = "NK", 
                                   `12` = "deadcell-remove", 
                                   `13` = "NK",
                                   `14` = "CD4T",
                                   `15` = "CD8T", 
                                   `16` = "CD4T", 
                                   `17` = "CD4T",  
                                   `18` = "CD8T", 
                                   `19` = "CD4T", 
                                   `20` = "Epithelial", 
                                   `21` = "CD4T", 
                                   `22` = "NK", 
                                   `23` = "Endothelial",
                                   `24` = "Plasma", 
                                   `25` = "NK"
    )
    if (i=="Stromal"){#done
      all.integrated <- RenameIdents(all.integrated, 
                                     `0` = "Fibroblast", 
                                     `1` = "Endothelial", 
                                     `2` = "Fibroblast",
                                     `3` = "Fibroblast",  
                                     `4` = "Fibroblast", 
                                     `5` = "Epithelial", 
                                     `6` = "Fibroblast", 
                                     `7` = "Pericyte",
                                     `8` = "Pericyte", 
                                     `9` = "Epithelial",
                                     `10` = "Endothelial", 
                                     `11` = "Fibroblast", 
                                     `12` = "Fibroblast", 
                                     `13` = "Pericyte",
                                     `14` = "Fibroblast",
                                     `15` = "Fibroblast", 
                                     `16` = "Endothelial", 
                                     `17` = "Fibroblast",  
                                     `18` = "Endothelial", 
                                     `19` = "Fibroblast", 
                                     `20` = "Endothelial", 
                                     `21` = "Fibroblast", 
                                     `22` = "Fibroblast", 
                                     `23` = "Endothelial",
                                     `24` = "Fibroblast", 
                                     `25` = "Epithelial",
                                     `26` = "Epithelial"
      )}
  }
  
  if (i=="Epithelial"){all.integrated$Celltype_final = "Epithelial"}
  else {all.integrated$Celltype_final <-Idents(all.integrated)}
  # table(all.integrated$Cell_type_paper,all.integrated$Celltype_main)
  table(all.integrated$celltype_by_marker,all.integrated$Celltype_final)
  table(all.integrated$Celltype_final,all.integrated$seurat_clusters)
  save(all.integrated,file = paste0(outdir_subset,"/all_integrated_",i,".RData"))
}

#Visualization 
for (i in c("Myeloids","Lymphocytes","Stromal","Epithelial")){
  load(paste0(outdir_subset,"/all_integrated_",i,".RData"))
  
  DefaultAssay(all.integrated) <- "integrated"
  
  p1 <- DimPlot(all.integrated, reduction = reduction, group.by = "group")
  p2 <- DimPlot(all.integrated, reduction = reduction, group.by = "Celltype_final",label = TRUE)
  # p3 <- DimPlot(all.integrated, reduction = reduction, split.by = "group") #To visualize the two conditions side-by-side, we can use the split.by argument to show each condition colored by cluster.
  p3 <- DimPlot(all.integrated, reduction = reduction, group.by = "cancer_name")
  p4 <- DimPlot(all.integrated, reduction = reduction, split.by = "group" )
  p =  (p1 + p2 + p3)
  ggsave(plot=p,  path = outdir_subset_PCA, dpi = 100, width = 16, height = 6, units = "in" ,filename=paste0("/DimPlot_celltypefinal_",i,"_",reduction,"_seed",seed,"_res",res,".pdf"))
  
  p <- FeaturePlot(all.integrated, features = c("MDK","SPP1","PLAU"), max.cutoff = 3, label = TRUE, cols = c("grey", "red"))
  ggsave(plot=p,  path = outdir_subset_PCA, device = "pdf",dpi = 150, width = 12, height = 6 ,filename=paste0("Feature_Plot","_",i,"_",reduction,"_seed",seed,"_res",res,".pdf"))
  
  table(all.integrated$seurat_clusters)
  table(all.integrated$seurat_clusters,all.integrated$group)
  
  DefaultAssay(all.integrated) <- "RNA"
  # p <- FeaturePlot(all.integrated, features = unique(genes_to_check), max.cutoff = 3, label = TRUE,
  #                  cols = c("grey", "red"))
  # ggsave(plot=p,  path = outdir_subset_PCA, device = "png",dpi = 150, width = 18, height = 45, units = "in" ,filename=paste0("Feature_Plot","_",i,"_seed",seed,"_res",res,".png"))
  
  p <- DotPlot(all.integrated, features = unique(genes_to_check),
               assay='RNA' , ##can use this to check orignial or normalized expression
               group.by = 'Celltype_final' ##this can change to other group method like 'seurat_clusters'
  )  + scale_color_viridis_c()  + coord_flip() 
  ggsave(plot=p,  path = outdir_subset_PCA, dpi = 300, width = 24, height = 45, units = "cm" ,filename=paste0("check_markers_celltypefinal_",i,"_seed",seed,"_res",res,".pdf"))
  
  print(paste0(i, " done!"))
}

######################################################################################## 
#3. Assign cell type to each cancer data set ####
#obtain celltype info from each main type file
for (cancer in cancer_name){
  anno_final_all <- as.data.frame(matrix(ncol=4)) 
  colnames(anno_final_all) = c("BC","cancer_name","Celltype_final","celltype_by_marker")  
  for (maintype in c("Myeloids","Lymphocytes","Stromal","Epithelial")){
    load(paste0(outdir_subset,"/all_integrated_",maintype,".RData"))
    print(maintype)
    anno_final <- all.integrated@meta.data
    anno_final %>% head
    anno_final <- anno_final[anno_final$cancer_name==cancer,]
    anno_final <- anno_final[,c("BC","cancer_name","Celltype_final","celltype_by_marker")]
    anno_final_all <- rbind(anno_final_all,anno_final)
  }
  anno_final_all <- anno_final_all[2:dim(anno_final_all)[1],]
  write.table(anno_final_all, file =paste0(outdir, "/celltyped/Celltype_final","_",cancer,".txt"), sep = '\t', row.names = T, col.names = TRUE, quote = FALSE) 
  print(paste0(cancer," done!"))
  table(anno_final_all$celltype_by_marker,anno_final_all$Celltype_final)
}

file_list= list.files(path=outdir, pattern="all_integrated_singleR_main_*")
cancer_name= sapply(strsplit(file_list,".RData"),'[[',1)
cancer_name= sapply(strsplit(cancer_name,"_"),'[[',5)

for (i in 1:length(file_list)){
  load(paste0(outdir,"/",file_list[i]))
  
  #add final celltype
  anno_final_all <- read.csv(paste0(outdir, "/celltyped/Celltype_final","_",cancer_name[i],".txt"),sep="\t",header = T)
  anno_final_all <- anno_final_all[,c("BC","cancer_name","Celltype_final")]
  anno_final_all$BC = rownames(anno_final_all)
  head(anno_final_all)
  colnames(all.integrated@meta.data)
  all.integrated$BC = rownames(all.integrated@meta.data)
  all.integrated@meta.data <- merge(all.integrated@meta.data,anno_final_all, by = "BC",all.x=TRUE)
  rownames(all.integrated@meta.data) <- all.integrated$BC
  head(all.integrated@meta.data)
  
  #reorder metadata
  all.integrated@meta.data = all.integrated@meta.data[colnames(all.integrated@assays$RNA@counts),]
  
  #add Neurons from Colon cancer
  if(cancer_name[i] =="Colon"){
    all.integrated$Celltype_final[which(all.integrated$Celltype_main == "Neurons")] = "Neurons"
  }
  
  if(cancer_name[i] =="Lung"){
    all.integrated$Celltype_final[which(all.integrated$seurat_clusters == 11)] = "Fibroblast"
  }
  table(all.integrated$Celltype_final)
  
  table(all.integrated$Celltype_final,all.integrated$group)
  table(all.integrated$celltype_by_marker,all.integrated$Celltype_final)
  
  #remove dead cells
  Idents(all.integrated) <- "Celltype_final"
  all.integrated <- subset(all.integrated, subset = Celltype_final != "deadcell-remove") 
  save(all.integrated,file= paste0(outdir,"/all_integrated_Celltype_final_",cancer_name[i],".RData"))
}

#4.visualization ####
genes_to_check = c('PTPRC','CD45', ## immune cell marker ('PTPRC'='CD45')
                   'PTPRC','CD3D', 'CD3E','CD4', 'CD8A','FOXP3','KLRD1', ## Tcells
                   'CD8A','GZMK','CD8B','GZMB','TRAC',  ## CD8+T cells
                   # 'PDCD1','LAG3','HAVCR2','TOX',# CD8+T exhausted (only high with these markers) and Proliferating
                   # 'MKI67','TOP2A','STMN2', #CD8A+ T proliferating -except above, also highly express these genes
                   # 'JUNB', 'FOS', 'ATF3', 'HSPA1A','DNAJB1', 'DNAJB1',#CD8A+T exhausted IEG
                   'CD4','IL7R', ##CD4+T
                   'FOXP3','CD4','IL2RA','IKZF2','TNFRSF4',  ## Treg
                   'CD40LG','IL7R', ## T memory
                   'CD40LG','IL7R','STAT4','CD3G', ## T helper
                   'CD16','CD16A','CD16B','CCR5','CCR6','CD27','TNFRSF7','TNFSF7', #T gamma delta
                   'GNLY','KLRD1','KLRC1','KLRF1','GZMB','NKG7','NCAM1','HAVR2C','NCR1','FGFBP2','FCGR3A','PRF1',  ##NK(CD56=NCAM1) (but no cd3d,cd4,cd8a)
                   'ZNF683','CD8A','CD8B','CD3','NCAM1','CD56', ##NKT  
                   'CD19', 'CD79A','CD79B','MS4A1' , # B cells
                   'IGHG1', 'MZB1', 'SDC1','IGHG4','CD38',  # plasma 
                   'CD68', 'LYZ', 'AIF1', #myeloid
                   'VCAN','FCN1','S100A12', ## monocyte
                   'CD68', 'CD163', 'CD14',  'CD86', 'CCL22','S100A4','CD207','CCL17','XCR1', ## DC(belong to monocyte)
                   'CD40','CD80', 'HLA-DOB','DOB', # DC-activated
                   'CD1C','CLEC9A', ## mDCs
                   'LILRA4','IL3RA','TCF4','TCL1A','CLEC4C', ## pDCs
                   'CD68','CD163', 'LAMP2','LGALS3', 'MRC1','MSR1' ,'S100A8','CD14','CD11B','APOE','C1QA','C1QB','ITGMX','CD11C','ITGAM','CD11B', ## Macrophage (belong to monocyte)
                   'FCGR3B', ## Neutrophil
                   'CD33','KIT','VIM','MS4A2','TPSAB1','CPA3','TPBS2','ENPP3','SLC18A2',  ##Mast cells
                   'EPCAM', 'KRT19','KRT7','KRT8','KRT18','PROM1',  ## epi or tumor
                   'FGF7', 'MME','COL1A1','ACTA2','PECAM1','VWF' ,'PROX1','PDGFRA', ## Fibroblasts,Endothelial
                   'MME','CD10','ACTA2', 'COL1A1', 'FN1','BGN','DCN','FAP',     ##stromal_fibroblasts('MME'='CD10')
                   'PECAM1',"CD31", 'ENG','VWF', 'CD36',  ##stromal_endo('PECAM1'="CD31")
                   'MCAM','RGS5','NDUFA4','KCNE4','CD31','NG2','PDGFRB','CD146' #Pericytes
                   # 'ASGR1', # hepatocyte
                   # 'FXYD2',  # cholangiocyte
                   # 'CLDN18', 'SFTPA1', 'SFTPA2', 'SFTPC', #alveolar cell 
                   # 'S100B', 'PLP1' #enteric glia- enteric nervous system
)

#shortened known marker gene list
genes_to_check = c('CD19', 'CD79A','CD79B','MS4A1',
                   'IGHG1', 'MZB1', 'IGHG4',
                   'PTPRC','CD4','IL7R', 'FOXP3','CD3D', 'CD3E', #CD4T
                   'PTPRC','CD8A','CD8B','GZMK','CD3D', 'CD3E', #CD8T
                   'ZNF683','CD8A','CD8B','CD3','CD56', #NKT
                   'GNLY','KLRD1','KLRC1','KLRF1','GZMB','NKG7', #NK
                   'FCN1','S100A12','CD14', #mono
                   'CD68','CD163', 'LAMP2','APOE','C1QA','C1QB','ITGMX','CD11C','CD11B', #macro ,'LGALS3', 'MRC1','MSR1' ,'CD14','CD11B'
                   'CD68', 'CD163', 'CD14',  'CD86', 'CCL22','CD207','CCL17', 'LILRA4','IL3RA',#DC 'TCL1A','CLEC4C','CD40','HLA-DOB','DOB'
                   'KIT','MS4A2','TPSAB1', 'CPA3', #Mast
                   'ENG','VWF', 'CD36', #Endo
                   'CD10', 'COL1A1', 'FN1','ACTA2',#Fib  'BGN', 'DCN', 'MME',
                   'MCAM','RGS5',
                   'EPCAM', 'KRT19','KRT7','KRT8','KRT18',#epi 'PROM1'
                   'S100B', 'PLP1' #enteric glia- enteric nervous system
                   ) %>% unique

for (i in 1:length(cancer_name)){
  # i=1
  res=0.2
  reduction="umap"
  load(paste0(outdir,"/all_integrated_Celltype_final_",cancer_name[i],".RData"))
  #plot
  # cols = c("B" = "#4F7CBA","Monocyte" = "#1BA3C6","Plasma"="#2CB5C0","Macrophage" = "#21B087","CD4T" = "#33A65C",
  #          "Fibroblast" = "#57A337",  "CD8T" = "#A2B627","Epithelial" = "#D5BB21",  "NKT" = "#F8B620","Endothelial" = "#F89217",
  #          "NK" = "#F06719","DC" = "#E03426","Mast" = "#F64971","Pericyte" = "#FC719E","Neurons" = "#EB73B3") # color used in heatmap
  cols = c("B" = "#ff579f","Monocyte" = "#3ac4d0","Plasma"="#ff62c3","Macrophage" = "#38c094","CD4T" = "#f663e3",
           "Fibroblast" = "#eed163",  "CD8T" = "#db73fc","Epithelial" = "#f9766d",  "NKT" = "#ae88ff","Endothelial" = "#f4a666",
           "NK" = "#728ff0","DC" = "#4ba8db","Mast" = "#39b24b","Pericyte" = "#93c05f","Neurons" = "#d0e562") 

  # DefaultAssay(all.integrated) <- "integrated"
  # p1 <- DimPlot(all.integrated, reduction = reduction, group.by = "celltype_by_marker",label = F,raster=FALSE )
  # p2 <- DimPlot(all.integrated, reduction = reduction, group.by = "Celltype_final",label = F,raster=FALSE )
  # p3 <- DimPlot(all.integrated, reduction = reduction, group.by = "seurat_clusters",label = TRUE,raster=FALSE)
  # p4 <- DimPlot(all.integrated, reduction = reduction, group.by = "group" ,label = TRUE,raster=FALSE)
  # p =  (p1 + p2)/(p3 + p4)
  # ggsave(plot=p,  path = paste0(outdir, "/celltyped"), dpi = 100, width = 12, height = 6, units = "in" ,filename=paste0("DimPlot_celltypefinal_",cancer_name[i],"_",reduction,"_res",res,".pdf"))
  # 
  # p2 <- DimPlot(all.integrated, reduction = reduction, group.by = "Celltype_final",label = F,raster=FALSE,cols=cols)
  # ggsave(plot=p2,  path = paste0(outdir, "/celltyped"), dpi = 100, width = 8, height = 4 , units = "in" ,filename=paste0("DimPlot_celltypefinal_",cancer_name[i],"_",reduction,"_res",res,"_2.pdf"))

  DefaultAssay(all.integrated) <- "RNA"
  Idents(all.integrated) <- factor(Idents(all.integrated), levels= c('B','Plasma','CD4T','CD8T','NKT','NK','Monocyte','Macrophage','DC','Mast','Endothelial','Fibroblast','Pericyte','Epithelial','Neurons'))
  p= DotPlot(all.integrated, features = unique(genes_to_check),
               assay='RNA' , ##can use this to check original or normalized expression
               # group.by = 'Celltype_final' ##this can change to other group method like 'seurat_clusters'
  )  +    #scale_color_viridis_c() +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0) +
    theme(axis.text.x=element_text(angle=90,hjust=1))
  ggsave(plot=p,  path = paste0(outdir, "/celltyped"), width = 25, height = 8 ,filename=paste0("check_markers_celltypefinal_",cancer_name[i],"_res",res,".pdf"))

  # p = DotPlot(all.integrated, features = c('COL1A1','COL4A1','COL18A1','CLEC11A',"FN1","PLAU",'SPP1','MDK'),
  #              assay='RNA' , ##can use this to check original or normalized expression
  #              group.by = 'Celltype_final' ##this can change to other group method like 'seurat_clusters'
  # )  + theme(axis.text.x=element_text(angle=90,hjust=1)) +
  #   scale_colour_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0) +
  #   guides(color = guide_colorbar(title = 'Scaled Average Expression')) +
  #   ggtitle(cancer_name[i]) # +
  # # geom_point(mapping = aes_string(size = 'pct.exp', color = 'avg.exp')) #add if want to use average expression instead scale
  # ggsave(plot=p,  path = paste0(outdir, "/celltyped"), width = 10, height = 10 ,filename=paste0("Dotplot_check_topUR_scaled_",cancer_name[i],"_res",res,".pdf"))
  # 
  # all.integrated$new_id = paste0(all.integrated$Celltype_final,'_',all.integrated$group)
  # p = DotPlot(all.integrated, features = c('COL1A1','COL4A1','COL18A1','CLEC11A',"FN1","PLAU",'SPP1','MDK'),
  #             assay='RNA' , ##can use this to check original or normalized expression
  #             # split.by='group',
  #             group.by = 'new_id' ##this can change to other group method like 'seurat_clusters'
  # )  + theme(axis.text.x=element_text(angle=90,hjust=1)) +
  #   guides(color = guide_colorbar(title = 'Scaled Average Expression')) + 
  #   scale_colour_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0) + 
  #   ggtitle(cancer_name[i]) 
  # ggsave(plot=p,  path = paste0(outdir, "/celltyped"), width = 8, height = 10 ,filename=paste0("Dotplot_check_topUR_scaled_bygroup_",cancer_name[i],"_res",res,".pdf"))
 }

#plot the expression of DS of fib,wound healing, and hotair
for (i in 1:length(cancer_name)){
  res=0.2
  reduction="umap"
  load(paste0(outdir,"/all_integrated_Celltype_final_",cancer_name[i],".RData"))
  
  p = DotPlot(all.integrated, features = c('NFKB1', 'VIM', 'TGFB1', 'COL1A1', 'COL1A2', 'COL3A1', 'MMP1', 'MMP10', 'MMP9'),
              assay='RNA' , ##can use this to check original or normalized expression
              # cols = c("white", "red"),
              group.by = 'Celltype_final' ##this can change to other group method like 'seurat_clusters'
  )  + theme(axis.text.x=element_text(angle=90,hjust=1)) +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0) +
    # scale_color_viridis_c() +
    guides(color = guide_colorbar(title = 'Scaled Average Expression')) +
    ggtitle(cancer_name[i]) # +
  # geom_point(mapping = aes_string(size = 'pct.exp', color = 'avg.exp')) #add if want to use average expression instead scale
  ggsave(plot=p,  path = paste0(outdir, "/celltyped"), width = 10, height = 10 ,filename=paste0("Dotplot_check_sharedDSof3pathways_scaled",cancer_name[i],"_res",res,".pdf"))
  

  all.integrated$new_id = paste0(all.integrated$Celltype_final,'_',all.integrated$group)
  p = DotPlot(all.integrated, features = c('NFKB1', 'VIM', 'TGFB1', 'COL1A1', 'COL1A2', 'COL3A1', 'MMP1', 'MMP10', 'MMP9','PLAUR'),
              assay='RNA' , ##can use this to check original or normalized expression
              # cols = c("white", "red"),
              # split.by='group',
              group.by = 'new_id' ##this can change to other group method like 'seurat_clusters'
  )  + theme(axis.text.x=element_text(angle=90,hjust=1)) +
    guides(color = guide_colorbar(title = 'Scaled Average Expression')) +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0) +
    # scale_color_viridis_c() +
    ggtitle(cancer_name[i])
  ggsave(plot=p,  path = paste0(outdir, "/celltyped"), width = 8, height = 10 ,filename=paste0("Dotplot_check_sharedDSof3pathways_scaled_bygroup_",cancer_name[i],"_res",res,".pdf"))
}


#plot the expression of DS in macrophage of plau 
for (i in 1:length(cancer_name)){
  res=0.2
  reduction="umap"
  load(paste0(outdir,"/all_integrated_Celltype_final_",cancer_name[i],".RData"))
  
  p = DotPlot(all.integrated, features = c('CCND1','FN1','JUNB','MMP9','VEGFA','ZFP36','C5AR1','ACTA2','MYC','BCL2L1','FOSL1','JUN','SRC','STAT3'),
              assay='RNA' , ##can use this to check original or normalized expression
              # cols = c("white", "red"),
              group.by = 'Celltype_final' ##this can change to other group method like 'seurat_clusters'
  )  + theme(axis.text.x=element_text(angle=90,hjust=1)) +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0) +
    # scale_color_viridis_c() +
    guides(color = guide_colorbar(title = 'Scaled Average Expression')) +
    ggtitle(cancer_name[i]) # +
  # geom_point(mapping = aes_string(size = 'pct.exp', color = 'avg.exp')) #add if want to use average expression instead scale
  ggsave(plot=p,  path = paste0(outdir, "/celltyped"), width = 10, height = 10 ,filename=paste0("Dotplot_check_DS_PLAU_PLAUR_scaled",cancer_name[i],"_res",res,".pdf"))
  
  
  all.integrated$new_id = paste0(all.integrated$Celltype_final,'_',all.integrated$group)
  p = DotPlot(all.integrated, features = c('CCND1','FN1','JUNB','MMP9','VEGFA','ZFP36','C5AR1','ACTA2','MYC','BCL2L1','FOSL1','JUN','SRC','STAT3'),
              assay='RNA' , ##can use this to check original or normalized expression
              # cols = c("white", "red"),
              # split.by='group',
              group.by = 'new_id' ##this can change to other group method like 'seurat_clusters'
  )  + theme(axis.text.x=element_text(angle=90,hjust=1)) +
    guides(color = guide_colorbar(title = 'Scaled Average Expression')) +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0) +
    # scale_color_viridis_c() +
    ggtitle(cancer_name[i])
  ggsave(plot=p,  path = paste0(outdir, "/celltyped"), width = 8, height = 10 ,filename=paste0("Dotplot_check_DS_PLAU_PLAUR_scaled_scaled_bygroup_",cancer_name[i],"_res",res,".pdf"))
}

##plot scMCTM####
library(seriation)
allscUR = read.csv(paste0("/home/yelzh67/Projects/cancer/datafile/outputs_5cancers/new_outputs/MCDM/scUR_celltype.csv"),sep="," ) %>% .$URs_same_fc  %>% unique
allscDS = read.csv(paste0("/home/yelzh67/Projects/cancer/datafile/outputs_5cancers/new_outputs/MCDM/scDS.csv"),sep="\t") %>% .$x  %>% unique
genes_to_check_scMCTM = c(allscUR,allscDS) %>% unique()

for (i in 1:length(cancer_name)){
  res=0.2
  reduction="umap"
  load(paste0(outdir,"/all_integrated_Celltype_final_",cancer_name[i],".RData"))
  
  p = DotPlot(all.integrated, features = genes_to_check_scMCTM,
              assay='RNA' , ##can use this to check original or normalized expression
              # cols = c("white", "red"),
              group.by = 'Celltype_final' ##this can change to other group method like 'seurat_clusters'
  )  + theme(axis.text.x=element_text(angle=90,hjust=1)) +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0) +
    # scale_color_viridis_c() +
    guides(color = guide_colorbar(title = 'Scaled Average Expression')) +
    ggtitle(cancer_name[i]) # +
  # geom_point(mapping = aes_string(size = 'pct.exp', color = 'avg.exp')) #add if want to use average expression instead scale
  ggsave(plot=p,  path = paste0(outdir, "/celltyped"), width = 30, height = 10 ,filename=paste0("Dotplot_check_scMCTM_scaled",cancer_name[i],"_res",res,".pdf"))
  data_for_plot = p$data # get data of average expression and percentage expression
  write.csv(data_for_plot, file = paste0(outdir,"/celltyped/averageExpression_percent_of_shMCTM_",cancer_name[i],".csv"),sep="\t",quote=F)
  # 
  # #extract the two matrices: avg.exp.scaled and pct.exp
  # df_avg = dcast(data_for_plot,features.plot~id, value.var ="avg.exp.scaled")
  # rownames(df_avg) = df_avg$features.plot
  # df_avg = df_avg[,2:dim(df_avg)[2]]
  # #for pct.exp
  # df_pct = dcast(data_for_plot,features.plot~id, value.var ="pct.exp")
  # rownames(df_pct) = df_pct$features.plot
  # df_pct = df_pct[,2:dim(df_pct)[2]]
  # 
  # #get order by euclidean distance for columns/clusters
  # coldist1 <- dist(t(scale(df_avg)))
  # coldist2 <- dist(t(scale(df_pct)))
  # colorder_ser <- seriate((coldist1+coldist2)/2, "OLO")
  # colidx <- get_order(colorder_ser)
  # ID_new_order <- colnames(df_avg)[colidx]
  # #get order by euclidean distance for rows/genes
  # rowdist1 <- dist(scale(df_avg))
  # rowdist2 <- dist(scale(df_pct))
  # roworder_ser <- seriate((rowdist1+rowdist2)/2, "OLO")
  # rowidx <- get_order(roworder_ser)
  # genes_to_check_new_order <- rownames(df_avg)[rowidx]
  # 
  # # Idents(all.integrated) <- factor(Idents(all.integrated), levels= ID_new_order)
  # #plot
  # all.integrated$new_id = paste0(all.integrated$Celltype_final,'_',all.integrated$group)
  # p= DotPlot(object = all.integrated, features = genes_to_check_new_order, assay="RNA", group.by = 'new_id')  +
  #   guides(color = guide_colorbar(title = 'Scaled Average Expression')) +RotatedAxis() +
  #   scale_colour_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0)
  # ggsave(plot=p,  path = paste0(outdir, "/celltyped"), device = "pdf", width = 30, height = 7 ,filename=paste0("DotPlot_shMCTM_bygroup.ordered_",cancer_name[i],"_res",res,".pdf"))
  # 
  # p = DotPlot(all.integrated, features = genes_to_check_new_order,assay='RNA' , ##can use this to check original or normalized expression
  #             group.by = 'Celltype_final',cluster.idents=T
  # )  + theme(axis.text.x=element_text(angle=90,hjust=1)) +
  #   scale_colour_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0) +
  #   # scale_color_viridis_c() +
  #   guides(color = guide_colorbar(title = 'Scaled Average Expression')) +
  #   ggtitle(cancer_name[i]) # +
  # # geom_point(mapping = aes_string(size = 'pct.exp', color = 'avg.exp')) #add if want to use average expression instead scale
  # ggsave(plot=p,  path = paste0(outdir, "/celltyped"), width = 30, height = 10 ,filename=paste0("Dotplot_scMCTM.orderd_",cancer_name[i],"_res",res,".pdf"))
  }
