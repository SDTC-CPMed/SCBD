
#analysis based on each subset/Fibroblast
rm(list = ls())
library(tidyverse)
library(nichenetr)
library(Seurat) # Seurat V4
library(MAST)
library(patchwork)
library(ggplot2)
library(SingleR)
library(SummarizedExperiment)
library(reshape2)

#mac or omika
location = "mac"

cancer_name = c("Lung","Liver","Colon","Ovary","Breast")
if (location == "mac"){  wd="/Users/cynthia_ye/Documents/OneDrive - Linköpings universitet"}
if (location == "omika"){ wd="/home/yelzh67/Projects"}

if (location == "mac"){
  cancer_seuratobj_list = c("/Users/cynthia_ye/Documents/OneDrive/5. Linkoping University/2.5 Cancer/1. GSE123902 LA/datafile/DEGs/Aug16_seurat_int/all_integrated_annotated.RData",
                            "/Users/cynthia_ye/Documents/OneDrive/5. Linkoping University/2.5 Cancer/2. GSE138709 ICC/datafile/Seurat_int/Aug20/OneDrive_2_23-08-2021/all_integrated_annotated.RData",
                            "/Users/cynthia_ye/Documents/OneDrive/5. Linkoping University/2.5 Cancer/5. GSE144735 CRC_belg/datafile/yl_seurat_DEGs/Aug17_seurat_inte_DEGs/OneDrive_1_24-08-2021/all.integrated_annotated.RData",
                            "/Users/cynthia_ye/Documents/OneDrive/5. Linkoping University/2.5 Cancer/6. PMID32561858_Ovary_C/datafile/DEGs/DEGs_MAST_Aug22/OneDrive_5_24-08-2021/all.integrated_annotated.RData",
                            "/Users/cynthia_ye/Documents/OneDrive/5. Linkoping University/2.5 Cancer/4. GSE161529 BC /datafile/DEGs_seurat/Aug17/all.integrated_annotated.RData"
  )
  outdir = "/Users/cynthia_ye/Documents/OneDrive - Linköpings universitet/cancer/datafile/outputs_5cancers/new_outputs/Seurat_subset_celltyping/2022Feb12"
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

outdir_subset_final = paste0(outdir,"_subset/subset_final")
if (dir.exists(outdir_subset_final)!=T){dir.create(outdir_subset_final)}
######################################################################################## 
#1 pilot - look all cell types####
seed = 42
nPCs = 30
res= 0.8
reduction= "umap"
genes_to_check = c('PTPRC','CD45', ## immune cell marker ('PTPRC'='CD45')
                   'PTPRC','CD3D', 'CD3E','CD4', 'CD8A','FOXP3','KLRD1', ## Tcells
                   'CD8A','GZMK','CD8B',  ## CD8+T cells
                   # 'PDCD1','LAG3','HAVCR2','TOX',# CD8+T exhausted (only high with these markers) and Proliferating
                   # 'MKI67','TOP2A','STMN2', #CD8A+ T proliferating -except above, also highly express these genes
                   # 'JUNB', 'FOS', 'ATF3', 'HSPA1A','DNAJB1', 'DNAJB1',#CD8A+T exhausted IEG
                   'GZMB', # shared by CTL and NK 
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
                   'CD68', 'CD163', 'CD14',  'CD86', 'CCL22','S100A4','CD207','CCL17', ## DC(belong to monocyte)
                   'CD40','CD80', 'HLA-DOB','DOB', # DC-activated
                   'CD1C','CLEC9A', ## mDCs
                   'LILRA4','IL3RA','TCF4','TCL1A','CLEC4C', ## pDCs
                   'CD68','CD163', 'LAMP2','LGALS3', 'MRC1','MSR1' ,'S100A8','CD14','CD11B','APOE','C1QA','C1QB','ITGMX','CD11C','ITGAM','CD11B', ## Macrophage (belong to monocyte)
                   'FCGR3B', ## Neutrophil
                   'CD33','KIT','VIM','MS4A2','TPSAB1', 'CPA3',  ##Mast cells
                   'EPCAM', 'KRT19','KRT7','KRT8','KRT18', 'PROM1',  ## epi or tumor
                   'FGF7', 'MME','COL1A1','ACTA2','PECAM1', 'VWF' ,'PROX1','PDGFRA', ## Fibroblasts,Endothelial
                   'MME','CD10','ACTA2', 'COL1A1', 'FN1', 'BGN', 'DCN',     ##stromal_fibroblasts('MME'='CD10')
                   'PECAM1',"CD31", 'ENG','VWF', 'CD36',  ##stromal_endo('PECAM1'="CD31")
                   'MCAM','RGS5' #Pericytes
                   # 'ASGR1', # hepatocyte
                   # 'FXYD2',  # cholangiocyte
                   # 'CLDN18', 'SFTPA1', 'SFTPA2', 'SFTPC', #alveolar cell 
                   # 'S100B', 'PLP1' #enteric glia- enteric nervous system
)

##. Visualization for major cell type####
for (i in c("Myeloids","Lymphocytes","Stromal","Epithelial")){
  load(paste0(outdir_subset,"/all_integrated_",i,".RData"))
  
  DefaultAssay(all.integrated) <- "integrated"
  # 
  # p1 <- DimPlot(all.integrated, reduction = reduction, group.by = "group")
  # p2 <- DimPlot(all.integrated, reduction = reduction, group.by = "Celltype_final",label = TRUE)
  # # p3 <- DimPlot(all.integrated, reduction = reduction, split.by = "group") #To visualize the two conditions side-by-side, we can use the split.by argument to show each condition colored by cluster.
  # p3 <- DimPlot(all.integrated, reduction = reduction, group.by = "cancer_name")
  # p4 <- DimPlot(all.integrated, reduction = reduction, split.by = "group" )
  # p =  (p1 + p2 + p3)
  # ggsave(plot=p,  path = outdir_subset_PCA, dpi = 100, width = 16, height = 6, units = "in" ,filename=paste0("/DimPlot_celltypefinal_",i,"_",reduction,"_seed",seed,"_res",res,".pdf"))
  # 
  FeaturePlot(all.integrated, features = c("MDK","SPP1","PLAU",'HOXAS'),split.by = "group", max.cutoff = 3, label = TRUE, cols = c("grey", "red"))
  p <- FeaturePlot(all.integrated, features = c("MDK","SPP1","PLAU",'IL6','LIF','CCL5','CSF3','ACTA2','FAP','CTGF'),split.by = "group", max.cutoff = 3, label = TRUE, cols = c("grey", "red"))
  ggsave(plot=p,  path = outdir_subset_PCA, device = "pdf", width = 12, height = 45 ,filename=paste0("Feature_Plot","_",i,"_bygroup_",reduction,"_seed",seed,"_res",res,".pdf"))
  
  p <- FeaturePlot(all.integrated, features = c("MDK","SPP1","PLAU",'IL6','LIF','CCL5','CSF3','ACTA2','FAP','CTGF'),split.by = "cancer_name", max.cutoff = 3, label = TRUE)
  ggsave(plot=p,  path = outdir_subset_PCA, device = "pdf",width = 25, height = 45 ,filename=paste0("Feature_Plot","_",i,"_bycancer_",reduction,"_seed",seed,"_res",res,".pdf"))
  
  
  # table(all.integrated$seurat_clusters)
  # table(all.integrated$seurat_clusters,all.integrated$group)
  # 
  # DefaultAssay(all.integrated) <- "RNA"
  # # p <- FeaturePlot(all.integrated, features = unique(genes_to_check), max.cutoff = 3, label = TRUE,
  # #                  cols = c("grey", "red"))
  # # ggsave(plot=p,  path = outdir_subset_PCA, device = "png",dpi = 150, width = 18, height = 45, units = "in" ,filename=paste0("Feature_Plot","_",i,"_seed",seed,"_res",res,".png"))
  # 
  # p <- DotPlot(all.integrated, features = unique(genes_to_check),
  #              assay='RNA' , ##can use this to check orignial or normalized expression
  #              group.by = 'Celltype_final' ##this can change to other group method like 'seurat_clusters'
  # )  + scale_color_viridis_c()  + coord_flip() 
  # ggsave(plot=p,  path = outdir_subset_PCA, dpi = 300, width = 30, height = 45, units = "cm" ,filename=paste0("check_markers_celltypefinal_",i,"_seed",seed,"_res",res,".pdf"))
  # 
  print(paste0(i, " done!"))
}

######################################################################################## 
#2. re-sub merging for Fibroblast####
##2.1. subset celltype from each cancers ####
Fibroblast_sublist = list()

outdir = "/home/yelzh67/Projects/cancer/datafile/outputs_5cancers/new_outputs/Seurat_subset_celltyping/2022Feb12"
file_list= list.files(path=outdir, pattern="all_integrated_Celltype_final_")
cancer_name= sapply(strsplit(file_list,".RData"),'[[',1)
cancer_name= sapply(strsplit(cancer_name,"_"),'[[',5)

for (i in 1:length(file_list)) {
  load(paste0(outdir,"/",file_list[i]))
  print(i)
  DefaultAssay(all.integrated) <- "RNA"
  seuratObj = all.integrated
  seuratObj@meta.data %>% head()
  
  seuratObj@meta.data$Celltype_final %>% table()
  seuratObj@meta.data$cancer_name = cancer_name[i]
  
  Fibroblast_sub <- subset(x = seuratObj, 
                        subset = Celltype_final == "Fibroblast")
  
  Fibroblast_sublist = c(Fibroblast_sublist,Fibroblast_sub)
  
  print(paste0(cancer_name[i], " done!"))
  rm(all.integrated,seuratObj)
}

save(Fibroblast_sublist,file = paste0(outdir_subset_final,"/Fibroblast_sublist.RData"))

Fibroblast_sublist[[1]]$Celltype_final %>% table
Fibroblast_sublist[[1]]@meta.data %>% head()


##2.2. merge the same celltype across cancers - seurat integration####
library(SeuratDisk)
# remotes::install_github("lyc-1995/MySeuratWrappers")
library(MySeuratWrappers)
library(SeuratWrappers)
###2.2.1 select integration features####
Fibroblast_sublist <- lapply(X = Fibroblast_sublist, FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = Fibroblast_sublist)
save(features,file = paste0(outdir_subset_final,"/features_Fibroblast.RData"))

###2.2.2 integration  #### 
load(paste0(outdir_subset_final,"/features_Fibroblast.RData"))
anchors <- FindIntegrationAnchors(object.list = Fibroblast_sublist, anchor.features = features)
all.integrated <- IntegrateData(anchorset = anchors)
save(all.integrated,file = paste0(outdir_subset_final,"/all_integrated_Fibroblast_beforePCA.RData"))
rm(features)

###2.2.3 Clustering after integration #### 
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
for (i in c("Fibroblast")){
  i="Fibroblast"
  load(paste0(outdir_subset_final,"/all_integrated_",i,"_beforePCA.RData"))
  all.integrated <- subset(x = all.integrated, 
                           subset = (cancer_name != "nopericyteLung" & cancer_name !="withpericyteLung"))
  all.integrated@meta.data <- all.integrated@meta.data[colnames(all.integrated@assays$RNA@counts),]
  DefaultAssay(all.integrated) <- "integrated"
  seed = 42
  nPCs = 10
  res= 0.1
  reduction= "umap"
  
  all.integrated <- ScaleData(all.integrated, verbose = T)
  all.integrated <- RunPCA(all.integrated,seed.use = seed, verbose = T)
  p <- ElbowPlot(all.integrated,ndims = 30)
  ggsave(plot=p,  path = outdir_subset_final, dpi = 300, width = 25, height = 20, units = "cm" ,filename=paste0(i,"_PCA-ElbowPlot.pdf"))
  all.integrated <- RunUMAP(all.integrated, seed.use = seed, reduction = "pca", dims = seq(nPCs))
  all.integrated <- RunTSNE(all.integrated, seed.use = seed, reduction = "pca", dims=seq(nPCs))
  all.integrated <- FindNeighbors(all.integrated, reduction = "pca", dims = seq(nPCs))
  all.integrated <- FindClusters(all.integrated, random.seed = seed, resolution = res)
  table(all.integrated$seurat_clusters)
  
  # DefaultAssay(all.integrated) <- "RNA"
  # all.integrated <- NormalizeData(all.integrated)
  # all.integrated <- ScaleData(all.integrated, verbose = T)
  all.integrated@assays$RNA@counts[1:70,1:2]
  all.integrated@assays$RNA@data[1:70,1:2]
  all.integrated@assays$RNA@scale.data[1:70,1:2]
  
  all.integrated@assays$integrated@counts[1:70,1:2]
  all.integrated@assays$integrated@data[1:70,1:2]
  all.integrated@assays$integrated@scale.data[1:70,1:2]
  
  save(all.integrated, file = paste0(outdir_subset_final,"/all_integrated_",i,".RData"))
}

#visualization
for (i in c("Fibroblast")){
  i="Fibroblast"
  load(paste0(outdir_subset_final,"/all_integrated_",i,".RData"))
  all.integrated@meta.data <- all.integrated@meta.data[colnames(all.integrated@assays$RNA@counts),]
  DefaultAssay(all.integrated) <- "integrated"
  seed = 42
  nPCs = 10
  res= 0.1
  reduction= "umap"
  
  #DimPlot
  p1 <- DimPlot(all.integrated, reduction = reduction, group.by = "group")
  p2 <- DimPlot(all.integrated, reduction = reduction, label = TRUE)
  p3 <- DimPlot(all.integrated, reduction = reduction, group.by = "cancer_name")
  p =  (p1 + p2 + p3)
  ggsave(plot=p,  path = outdir_subset_final, width = 9, height = 3,filename=paste0("/DimPlot_celltypefinal_",i,"_",reduction,"_seed",seed,"_res",res,".pdf"))
  
  p <- DimPlot(all.integrated, reduction = reduction, split.by = "group",label=T) #To visualize the two conditions side-by-side, we can use the split.by argument to show each condition colored by cluster.
  ggsave(plot=p,  path = outdir_subset_final, width = 6, height = 3,filename=paste0("/DimPlot_celltypefinal_bygroup_",i,"_",reduction,"_seed",seed,"_res",res,".pdf"))
  
  p <- DimPlot(all.integrated, reduction = reduction, split.by = "cancer_name",label=T) #To visualize the two conditions side-by-side, we can use the split.by argument to show each condition colored by cluster.
  ggsave(plot=p,  path = outdir_subset_final, width = 9, height = 2,filename=paste0("/DimPlot_celltypefinal_bycancer_",i,"_",reduction,"_seed",seed,"_res",res,".pdf"))
  
  #Dotplot
  DefaultAssay(all.integrated) <- "RNA"
  #C5AR2 = GPR77, CD10=MME,CD146=MCAM/MUC18, CAV=CAV1,NID2=Nidogen2
  markercluster0 = c("COL1A1" ,"FN1" ,"MMP11","CTHRC1" ,  "COL1A2","COL3A1,SPARC","COL5A2","POSTN")
  top8UR = c("FN1","COL1A1","PLAU" ,"CLEC11A","MDK")
  DS_pw = c('NFKB1', 'VIM', 'TGFB1', 'COL1A1', 'COL1A2', 'COL3A1', 'MMP1', 'MMP10', 'MMP9')
  UR_pw = c('CYR61','RARRES2' )
  CAF_markers = c('ACTA2','IL6','IL10','LIF','FBLN1','PDGFRA','POSTN','KRT19','CD74')
  genes_to_check = c(CAF_markers,markercluster0,top8UR) %>% unique
  
  # genes_to_check = c("PLAU",'PLAUR','IL6','IL11','LIF','ACTA2','COL11A1','FN1','COL1A1','SPP1', 'COL4A1', 'COL18A1',  'CLEC11A', 'MDK','MMP3','IL17A','FOS','JUN','MMP11','CTHRC1','FAP','PDGFRB','CD74','POSTN','COL5A1','KRT19','KRT17','EPCAM','NID2','FBLN1','PDGFRA','MMP2','DCN','COL1A2','TAGLN')
  # genes_to_check = c('COL1A1','CLEC11A','FN1',"PLAU",'MDK','MMP11','CTHRC1','COL1A2','ACTA2','IL6','LIF','FBLN1','PDGFRA','CD74')
  
  all.integrated$new_id = paste0(all.integrated$seurat_clusters,'_',all.integrated$group)
  p=DotPlot(object = all.integrated, features = CAF_markers, 
            assay="RNA",group.by='new_id') + guides(color = guide_colorbar(title = 'Scaled Average Expression'))+
    scale_colour_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0)+
    theme(axis.text.x=element_text(angle=90,hjust=1))
  ggsave(plot=p,  path = outdir_subset_final, device = "pdf", width = 7, height = 7 ,filename=paste0("DotPlot_CAF_markers","_",i,"_bygroup_",reduction,"_seed",seed,"_res",res,".pdf"))
  
  p=DotPlot(object = all.integrated, features = CAF_markers, 
            assay="RNA",cols = c("white", "red"))+RotatedAxis()  + 
    scale_colour_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0) +
    #scale_color_viridis_c(option = "cividis")+
    # geom_point(mapping = aes_string(size = 'pct.exp', color = 'avg.exp'))+ #add if want to use average expression instead scale
    guides(color = guide_colorbar(title = 'Scaled Average Expression'))+ 
    theme(axis.text.x=element_text(angle=90,hjust=1)) 
  ggsave(plot=p,  path = outdir_subset_final, device = "pdf", width = 7, height = 7 ,filename=paste0("DotPlot_CAF_markers","_",i,"_",reduction,"_seed",seed,"_res",res,".pdf"))
  # 
  # #VlnPlot
  # p=VlnPlot(all.integrated, assay="RNA",slot="data",
  #          group.by= "seurat_clusters",
  #         # split.by='group', #cancer_name
  #         # split.plot = T,
  #         features = c("PLAU",'PLAUR'),ncol=1,pt.size = 0.1) + theme(legend.position = "right")  #+stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.5)
  # ggsave(plot=p,  path = outdir_subset_final, device = "pdf", width = 12, height = 9 ,filename=paste0("VlnPlot","_",i,"_",reduction,"_seed",seed,"_res",res,".pdf"))
  # 
  # 
  # table(all.integrated$cancer_name,all.integrated$seurat_clusters)
  # 
  #FeaturePlot
  DefaultAssay(all.integrated) <- "RNA"
  p = FeaturePlot(all.integrated, features = c('IL6','ACTA2','FN1', 'COL1A1', 'MMP11', 'CTHRC1', 'COL1A2','POSTN','CD74', 'HLA-DRA','KRT19'),
                  # split.by = "group", 
                  ncol = 6,
                  max.cutoff = 3, label = TRUE, cols = c("grey", "red"))
  ggsave(plot=p,  path = outdir_subset_final, device = "pdf", width = 20, height = 8 ,filename=paste0("Feature_Plot","_",i,"_",reduction,"_seed",seed,"_res",res,".pdf"))

  p = FeaturePlot(all.integrated, features = c('ACTA2','FN1','IL6', 'HSPA1A', 'KRT7','CD74'),
                  # split.by = "group", 
                  ncol = 3,
                  max.cutoff = 3, label = TRUE, cols = c("grey", "red"))
  ggsave(plot=p,  path = outdir_subset_final, device = "pdf", width = 16, height = 8 ,filename=paste0("Feature_Plot","_",i,"_",reduction,"_seed",seed,"_res",res,".pdf"))
  
  p = FeaturePlot(all.integrated, features = c('IL6','ACTA2','FN1', 'COL1A1', 'MMP11', 'CTHRC1', 'COL1A2','POSTN','CD74', 'HLA-DRA','KRT19'),
                  split.by = "group",
                  ncol = 6,
                  max.cutoff = 3, label = TRUE, cols = c("grey", "red"))
  ggsave(plot=p,  path = outdir_subset_final, device = "pdf", width = 8, height = 22 ,filename=paste0("Feature_Plot","_",i,"_bygroup_",reduction,"_seed",seed,"_res",res,".pdf"))
  
  # p = FeaturePlot(all.integrated, features = c("PLAU",'IL6','ACTA2'),split.by = "cancer_name", max.cutoff = 3, label = TRUE,cols = c("grey", "red"))
  # ggsave(plot=p,  path = outdir_subset_final, device = "pdf",width = 15, height = 9 ,filename=paste0("Feature_Plot","_",i,"_bycancer_",reduction,"_seed",seed,"_res",res,".pdf"))
  # 
  # iCAF: 'IL6','LIF','CCL5','CSF3', #https://www.nature.com/articles/s41420-021-00410-6 #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5339682/
  # myCAF: 'ACTA2','FAP' #https://www.nature.com/articles/s41420-021-00410-6
  # tumor-promoting CAFs (pCAFs): 'FAP'
  # tumor-inhibiting CAFs (rCAFs):
  # neutral CAFs (nCAFs):
  
}

###2.2.4 cell count####
all.integrated$cancer_group = paste0(all.integrated$cancer_name,'_',all.integrated$group)
count = table(all.integrated$seurat_clusters,all.integrated$cancer_group) %>% as.matrix() %>% as.data.frame.matrix()
write.csv(count,file = paste0(outdir_subset_final,"/Cluster_cellcount.csv"),sep=",")
percent = apply(count,2,function(x){x/sum(x)})
percent.long = melt(percent, value.name = "Percentage")
colnames(percent.long) = c("Cluster","Cancer_group","Percentage")
percent.long$Cluster = factor(percent.long$Cluster )
percent.long$Cancer = lapply(as.character(percent.long$Cancer_group), function(x){strsplit(x,"_") %>% unlist %>% .[1]}) %>% unlist
percent.long$Group = lapply(as.character(percent.long$Cancer_group), function(x){strsplit(x,"_") %>% unlist %>% .[2]}) %>% unlist

head(percent.long)

p = ggplot(percent.long, aes(x = Group, y = Percentage,fill= Cluster)) + 
  geom_bar(stat = "identity", color = "black") + scale_color_viridis_c() +
  facet_grid(~Cancer,scales='fixed')
pdf(file=paste0(outdir_subset_final,"/Cluster_percentage.pdf"), width = 5, height = 5)
print(p)
dev.off()

##2.3 find markers for clusters ####
for (i in c("Fibroblast")){ #findallmarkers
  load(paste0(outdir_subset_final,"/all_integrated_",i,".RData"))
  print(i)
  DefaultAssay(all.integrated) <- "RNA"
  
  all.integrated@meta.data = all.integrated@meta.data[colnames(all.integrated@assays$RNA@counts),]
  Idents(all.integrated) <- all.integrated$seurat_clusters
  head(Idents(all.integrated))
  
  all.markers <- FindAllMarkers(all.integrated, only.pos = TRUE, 
                                min.pct = 0.1, logfc.threshold = 0.25)
  #把每个细胞群都计算一遍
  write.table(all.markers,
              file=paste0(outdir_subset_final,"/ClusterMarkers/",i,"_total_marker_genes_",nPCs,"PC.txt"),
              sep="\t",quote = F,row.names = F)
  
  all.markers = read.csv2(paste0(outdir_subset_final,"/ClusterMarkers/",i,"_total_marker_genes_",nPCs,"PC.txt"),
                          sep="\t")
  
  all.markers.ranked = all.markers[order(all.markers$cluster,-as.numeric(all.markers$avg_log2FC)),]
  
  # 遍历每一个cluster然后展示其中前10个基因
  marker.sig <- all.markers %>% 
    mutate(Ratio = round(as.numeric(pct.1)/as.numeric(pct.2),3)) %>%
    filter(as.numeric(all.markers$p_val_adj) <= 0.05)  # 本条件为过滤统计学不显著的基因
  
  for(cluster_id in unique(marker.sig$cluster)){
    # cluster.markers <- FindMarkers(experiment.aggregate, ident.1 = cluster, min.pct = 0.3)
    # cluster.markers <- as.data.frame(cluster.markers) %>%
    #   mutate(Gene = rownames(cluster.markers))
    cl4.genes <- marker.sig %>%
      filter(cluster == cluster_id) %>%
      arrange(desc(avg_log2FC))
    cl4.genes <- cl4.genes[1:min(nrow(cl4.genes),10),"gene"]
  }
  # rm(cl4.genes,cluster_id,pvn)
  
  #热图展示Top marker基因
  #筛选top的marker基因，可以通过参数改为其他数值
  top <- marker.sig %>% group_by(cluster) %>% 
    top_n(n = 30, wt = avg_log2FC)
  
  #top-marker基因dotplot
  pdf(paste0(outdir_subset_final,"/ClusterMarkers/MarkerGene-DotPlot_all_cluster_",i,"_","res",res,"_PC",nPCs,"PC.pdf"),width = 20,height = 5)
  DotPlot(all.integrated, features = unique(top$gene))+
    # scale_color_viridis_c()+
    RotatedAxis() +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0)
  dev.off()
  
  #top-marker基因热图
  DefaultAssay(all.integrated)  <- 'integrated'
  pdf(paste0(outdir_subset_final,"/ClusterMarkers/MarkerGene-Heatmap_all_cluster_",i,"_","res",res,"_PC",nPCs,"PC.pdf"),width= 10, height= 10 )
  DoHeatmap(all.integrated, features = top$gene,size = 2) +
    theme(legend.position = "none", 
          axis.text.y = element_text(size = 7))
  dev.off()
}

###2.3.0 dotPlot all scMCTM genes ####
load(paste0(outdir_subset_final,"/all_integrated_Fibroblast.RData"))
all.integrated@meta.data <- all.integrated@meta.data[colnames(all.integrated@assays$RNA@counts),]
# DefaultAssay(all.integrated) <- "integrated"
DefaultAssay(all.integrated) <- "RNA"
table(all.integrated$seurat_clusters)
  
  # # Genes from scMCTM
  # allscUR = read.csv("/Users/cynthia_ye/Documents/OneDrive - Linköpings universitet/cancer/datafile/outputs_5cancers/new_outputs/MCDM/scUR_celltype.csv",sep="," ) %>% .$URs_same_fc  %>% unique
  # allscDS = read.csv("/Users/cynthia_ye/Documents/OneDrive - Linköpings universitet/cancer/datafile/outputs_5cancers/new_outputs/MCDM/scDS.csv",sep="\t") %>% .$x  %>% unique
  # genes_to_check_scMCTM = c(allscUR,allscDS) %>% unique()
  # 
  # # select scUR scDS which have 1.5 scaled expression in Fib
  # scMCTM_exp.list = list.files(path= paste0(wd,"/cancer/datafile/outputs_5cancers/new_outputs/Seurat_subset_celltyping/2022Feb12/celltyped"),pattern="averageExpression_percent_of_shMCTM_") 
  # scMCTM_Fib.high = c()
  # for (x in 1:5){
  # a = read.csv2(paste0(wd,"/cancer/datafile/outputs_5cancers/new_outputs/Seurat_subset_celltyping/2022Feb12/celltyped/",scMCTM_exp.list[x]),sep=",")
  # a = a[(a$pct.exp >10 & a$avg.exp.scaled >1.5 & a$id =="Fibroblast"),"features.plot"]
  # scMCTM_Fib.high = c(scMCTM_Fib.high,a)
  # }
  # scMCTM_Fib.high = scMCTM_Fib.high %>% unique
  # 
  # allscUR.fib=intersect(scMCTM_Fib.high,allscUR)
  # allscDS.fib=intersect(scMCTM_Fib.high,allscDS)
  # write.csv(allscUR.fib,file="/Users/cynthia_ye/Downloads/shUR_Fib_high.csv",sep="\t",quote=F)
  # write.csv(allscDS.fib,file="/Users/cynthia_ye/Downloads/shDS_Fib_high.csv",sep="\t",quote=F)
  
# shMCTM genes that highly expressed in fibroblast
allscUR.fib = read.csv(paste0(wd, "/cancer/datafile/outputs_5cancers/new_outputs/MCDM/CT_high_shUR_shDS/Findallmarkers_celltype/marker.high.shUR.csv"))
allscUR.fib = allscUR.fib[allscUR.fib$type =="shUR" & allscUR.fib$celltype == "Fibroblast","gene"]
allscDS.fib = read.csv(paste0(wd,"/cancer/datafile/outputs_5cancers/new_outputs/MCDM/CT_high_shUR_shDS/Findallmarkers_celltype/marker.high.shDS.csv"))
allscDS.fib = allscDS.fib[allscDS.fib$type =="shDS" & allscDS.fib$celltype == "Fibroblast","gene"]
  
#order by cluster marker
all.markers = read.csv2(paste0(outdir_subset_final,"/ClusterMarkers/Fibroblast_total_marker_genes_",nPCs,"PC.txt"),
                          sep="\t")
all.markers.ranked = all.markers[order(all.markers$cluster,-as.numeric(all.markers$avg_log2FC)),]
  
allscUR.ranked = all.markers.ranked$gene[all.markers.ranked$gene %in% allscUR.fib] %>% unique
allscDS.ranked = all.markers.ranked$gene[all.markers.ranked$gene %in% allscDS.fib] %>% unique
# genes_to_check_scMCTM_FibHigh.ranked = all.markers.ranked$gene[all.markers.ranked$gene %in%  scMCTM_Fib.high] %>% unique
  
#plot - order based on cluster marker order
p= DotPlot(object = all.integrated, features = allscUR.ranked, assay="RNA")  +
    guides(color = guide_colorbar(title = 'Scaled Average Expression')) +RotatedAxis() +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0)
  data_for_plot = p$data # get data of average expression and percentage expression
ggsave(plot=p,  path = outdir_subset_final, device = "pdf", width = 8, height = 7 ,filename=paste0("DotPlot","_fibroblast_allscUR.ordered","_seed",seed,"_res",res,".pdf"))
  
all.integrated$new_id = paste0(all.integrated$seurat_clusters,'_',all.integrated$group)
p= DotPlot(object = all.integrated, features = allscUR.ranked, assay="RNA",group.by = "new_id")  +
    guides(color = guide_colorbar(title = 'Scaled Average Expression')) +RotatedAxis() +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0)
  data_for_plot = p$data # get data of average expression and percentage expression
ggsave(plot=p,  path = outdir_subset_final, device = "pdf", width = 8, height = 7 ,filename=paste0("DotPlot","_fibroblast_allscUR_bygroup.ordered","_seed",seed,"_res",res,".pdf"))
  
p= DotPlot(object = all.integrated, features = allscDS.ranked, assay="RNA")  +
    guides(color = guide_colorbar(title = 'Scaled Average Expression')) +RotatedAxis() +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0)
data_for_plot = p$data # get data of average expression and percentage expression
ggsave(plot=p,  path = outdir_subset_final, device = "pdf", width = 16, height = 7 ,filename=paste0("DotPlot","_fibroblast_allscDS.ordered","_seed",seed,"_res",res,".pdf"))
  
p= DotPlot(object = all.integrated, features = allscDS.ranked, assay="RNA",group.by = "new_id")  +
    guides(color = guide_colorbar(title = 'Scaled Average Expression')) +RotatedAxis() +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0)
data_for_plot = p$data # get data of average expression and percentage expression
ggsave(plot=p,  path = outdir_subset_final, device = "pdf", width = 16, height = 7 ,filename=paste0("DotPlot","_fibroblast_allscDS_bygroup.ordered","_seed",seed,"_res",res,".pdf"))
  
  # p= DotPlot(object = all.integrated, features = genes_to_check_scMCTM_FibHigh.ranked, assay="RNA")  +
  #   guides(color = guide_colorbar(title = 'Scaled Average Expression')) +RotatedAxis() +
  #   scale_colour_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0)
  # data_for_plot = p$data # get data of average expression and percentage expression
  # ggsave(plot=p,  path = outdir_subset_final, device = "pdf", width = 30, height = 7 ,filename=paste0("DotPlot","_fibroblast_allscMCTM.ordered","_seed",seed,"_res",res,".pdf"))
  # 
  # p= DotPlot(object = all.integrated, features = genes_to_check_scMCTM_FibHigh.ranked, assay="RNA",group.by = "new_id")  +
  #   guides(color = guide_colorbar(title = 'Scaled Average Expression')) +RotatedAxis() +
  #   scale_colour_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0)
  # data_for_plot = p$data # get data of average expression and percentage expression
  # ggsave(plot=p,  path = outdir_subset_final, device = "pdf", width = 30, height = 7 ,filename=paste0("DotPlot","_fibroblast_allscMCTM_bygroup.ordered","_seed",seed,"_res",res,".pdf"))
  # 
  # 
  # # order - based on distance
  # library(seriation)
  # p= DotPlot(object = all.integrated, features = c(allscUR.fib,allscDS.fib), assay="RNA",group.by = "new_id")  +
  #   guides(color = guide_colorbar(title = 'Scaled Average Expression')) +RotatedAxis() +
  #   scale_colour_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0)
  # data_for_plot = p$data # get data of average expression and percentage expression
  # 
  # for (x in c("allscUR.fib","allscDS.fib")){
  #   data_for_plot.x = data_for_plot[data_for_plot$features.plot %in% eval(parse(text = x)),]
  #   
  #   #extract the two matrices: avg.exp.scaled and pct.exp
  #   df_avg = dcast(data_for_plot.x,features.plot~id, value.var ="avg.exp.scaled")
  #   rownames(df_avg) = df_avg$features.plot
  #   df_avg = df_avg[,2:dim(df_avg)[2]]
  #   #for pct.exp
  #   df_pct = dcast(data_for_plot.x,features.plot~id, value.var ="pct.exp")
  #   rownames(df_pct) = df_pct$features.plot
  #   df_pct = df_pct[,2:dim(df_pct)[2]]
  #   
  #   #get order by euclidean distance for columns/clusters
  #   coldist1 <- dist(t(scale(df_avg)))
  #   coldist2 <- dist(t(scale(df_pct)))
  #   colorder_ser <- seriate((coldist1+coldist2)/2, "OLO_complete")
  #   colidx <- get_order(colorder_ser)
  #   ID_new_order <- colnames(df_avg)[colidx]
  #   #get order by euclidean distance for rows/genes
  #   rowdist1 <- dist(scale(df_avg))
  #   rowdist2 <- dist(scale(df_pct))
  #   roworder_ser <- seriate((rowdist1+rowdist2)/2, "OLO_complete")
  #   rowidx <- get_order(roworder_ser)
  #   genes_to_check_new_order <- rownames(df_avg)[rowidx]
  #   
  #   #plot
  #   # Idents(all.integrated) <- factor(Idents(all.integrated), levels= ID_new_order)
  #   p= DotPlot(object = all.integrated, features = genes_to_check_new_order, assay="RNA")  +
  #     guides(color = guide_colorbar(title = 'Scaled Average Expression')) +RotatedAxis() +
  #     scale_colour_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0)
  #   ggsave(plot=p,  path = outdir_subset_final, device = "pdf", width = 12, height = 7 ,filename=paste0("DotPlot","_fibroblast_",x,".dis.ordered","_seed",seed,"_res",res,".pdf"))
  # }

### 2.3.1. Check whether top 10 cluster0 marker genes were DEGs against control ####
#order by cluster marker
all.markers = read.csv2(paste0(outdir_subset_final,"/ClusterMarkers/Fibroblast_total_marker_genes_",nPCs,"PC.txt"),
                        sep="\t")
all.markers.ranked = all.markers[order(all.markers$cluster,-as.numeric(all.markers$avg_log2FC)),]
all.markers.ranked.C0 = all.markers.ranked[all.markers.ranked$cluster ==0,"gene"][1:10]

#DEGs tumor vs. control - all Fibroblast
DEGs_of_all_celltype_cancer = read.csv2(paste0(wd,"/cancer/datafile/outputs_5cancers/new_outputs/DEGs/merge_DEGs_of_all_celltype_cancer.csv"),sep=",",row.names = 1)
DEGs_of_all_celltype_cancer.fib = DEGs_of_all_celltype_cancer[DEGs_of_all_celltype_cancer$cell == "Fibroblast" & DEGs_of_all_celltype_cancer$X %in% all.markers.ranked.C0,]
table(DEGs_of_all_celltype_cancer.fib$X,DEGs_of_all_celltype_cancer.fib$cancer)

write.table(DEGs_of_all_celltype_cancer.fib,
            file=paste0(outdir_subset_final,"/ClusterMarkers/top10C0_DEG_TvsN_Fib.txt"),
            sep="\t",quote = F,row.names = F)

#DEGs tumor vs. control - in C0
load(paste0(outdir_subset_final,"/all_integrated_Fibroblast.RData"))
DefaultAssay(all.integrated) <- "RNA"

all.integrated@meta.data = all.integrated@meta.data[colnames(all.integrated@assays$RNA@counts),]
Idents(all.integrated) <- all.integrated$seurat_clusters
head(Idents(all.integrated))

DEG_method <- "MAST"  #wilcox or MAST和DESeq2
DEG = data.frame()
for (i in 0:(length(unique(all.integrated$seurat_clusters))-1)) {
    print(i)
  # for (j in 1:5){
  #   print(unique(temp$cancer_name)[j])
    Idents(all.integrated) <- "seurat_clusters"
    temp <- subset(all.integrated, idents = i)
    # temp <- subset(temp, cancer_name == unique(temp$cancer_name)[j])
    Idents(temp) <- "group"
    
    if(dim(table(temp$Celltype_final,temp$group))[2] ==1){next} #skip if cell type only exists in one group
    
    iDEGs <- FindMarkers(temp, ident.1 = "T", ident.2 = "N", verbose = FALSE,
                         min.pct = 0.1, #only select genes that expressed in more than 10% of cells 
                         test.use = DEG_method #"wilcox"(default), can change to wilcox or MAST和DESeq2
    )
    head(iDEGs, n = 15)
    iDEGs <- iDEGs %>% filter(p_val_adj <= 0.05)
    iDEGs$cluster = i
    DEG = rbind(DEG,iDEGs)
    # }
}
write.csv(DEG,paste0(outdir_subset_final,"/ClusterMarkers/DEG_TvsN_allclusters","_",DEG_method,".csv"),row.names = T, col.names = T, sep = ',',quote = F)

DEGs = read.csv2(paste0(outdir_subset_final,"/ClusterMarkers/DEG_TvsN_allclusters","_",DEG_method,".csv"),sep=",")
DEGs_of_all_celltype_cancer.fib = DEGs[DEGs$X %in% all.markers.ranked.C0,]

write.table(DEGs_of_all_celltype_cancer.fib,
            file=paste0(outdir_subset_final,"/ClusterMarkers/top10C0_DEG_TvsN_C0.txt"),
            sep="\t",quote = F,row.names = F)

### 2.3.2. GO,KEGG for cluster0 marker ####
  library(org.Hs.eg.db)
  library(org.Hm.eg.db)  
  library(clusterProfiler)  
  library(dplyr)  
  library(ggplot2)  
  library(DOSE)
  library(topGO)
  library(pathview)
  #Write function to create all GO, KEGG output####
  GO_Kegg_Enrichment_Vis <- function(x,datasource,nshow=10,fromtype="SYMBOL",totype=c("ENSEMBL","ENTREZID"),orgdb="org.Hs.eg.db", out.dir){
    KEGG_results <- matrix(ncol=9)
    colnames(KEGG_results) <- c('ID','Description','GeneRatio','BgRatio','pvalue','p.adjust','qvalue','geneID','Count')
    ##translate####
    test = bitr(x, #数据集
                fromType=fromtype,  
                toType=totype,   
                OrgDb=orgdb) 
    head(test)
    
    ##KEGG enrichment####
    kk <- enrichKEGG(gene = test$ENTREZID,
                     organism = 'hsa', #human ='hsa', mouse= 'mmu' #https://www.genome.jp/kegg/catalog/org_list.html 
                     pvalueCutoff = 1)
    head(kk,2)
  
    #show results
    write.csv(summary(kk),paste0(out.dir,"/KEGG-enrich_",datasource,".csv"),row.names =FALSE)
    write.csv(kk,file = paste0(out.dir,"/KEGG-enrich_all_",datasource,".csv"),row.names =T)
    KEGG_results <- kk
    barplot
    p = dotplot(kk,orderBy = "x",title="Enrichment KEGG_dot",showCategory=nshow)
    ggsave(plot=p,path = paste0(out.dir), filename = paste0("KEGG-enrich_dotplot_",datasource,".pdf"),dpi = 300, width = 15, height = 5, units = "cm")
    p = barplot(kk,orderBy = "x",title="Enrichment KEGG_dot",showCategory=nshow)
    ggsave(plot=p,path = paste0(out.dir), filename = paste0("KEGG-enrich_barplot_",datasource,".pdf"),dpi = 300, width = 15, height = 5, units = "cm")

    ##GO enrichment ####
    ggo <- groupGO(gene = test$ENTREZID, OrgDb = org.Hs.eg.db, ont = "CC",level = 3,readable = TRUE)
    ego_ALL <- enrichGO(gene = test$ENTREZID,  #基因列表(转换的ID)
                        keyType = "ENTREZID",  #指定的基因ID类型，默认为ENTREZID
                        OrgDb=org.Hs.eg.db,  #物种对应的org包
                        # universe = names(geneList), #背景基因集  If missing, the all genes listed in the database (eg TERM2GENE table) will be used as background.
                        ont = "ALL",   #CC细胞组件，MF分子功能，BP生物学过程
                        pvalueCutoff = 0.01,  #p值阈值
                        pAdjustMethod = "fdr",  #多重假设检验校正方式
                        minGSSize = 1,   #注释的最小基因集，默认为10
                        maxGSSize = 500,  #注释的最大基因集，默认为500
                        qvalueCutoff = 0.05,  #p值阈值
                        readable = TRUE)  #基因ID转换为基因名
    dim(ego_ALL)
    head(ego_ALL@result)
    ego_ALL@result$ONTOLOGY %>% table
    ego_ALL <- setReadable(ego_ALL, OrgDb = org.Hs.eg.db)
    write.csv(summary(ego_ALL),paste0(out.dir,"/GO-all_enrich_",datasource,".csv"),row.names =FALSE)

  #   ##Visualization####
  #   library(stringr)
  #   library(cowplot)
  #   library(ggplot2)
  #   p = dotplot(ego_ALL,
  #               showCategory=nshow,
  #               font.size=10,
  #               split="ONTOLOGY",
  #               color ="pvalue",
  #               title="EnrichmentGO_all_dot")+ 
  #     facet_grid(ONTOLOGY~.,scale="free")+
  #     scale_y_discrete(labels=function(x) str_wrap(x, width=80)) #点图，按富集的数从大到小的 
  #   ggsave(plot=p,path = paste0(out.dir), filename = paste0("GO-all_dotplot_",datasource,".pdf"),dpi = 300, width = 25, height = 20, units = "cm")
  #   
    ##可视化--条形图
    p=barplot(ego_ALL,showCategory=nshow,
              font.size=5,
              split="ONTOLOGY",
              color ="pvalue",title="EnrichmentGO_ALL")+
      facet_grid(ONTOLOGY~.,scale="free")+
      scale_y_discrete(labels=function(x) str_wrap(x, width=60))#条状图，按p从小到大排，绘制前10个Term
    ggsave(plot=p,path = paste0(out.dir), filename = paste0("GO-all_barplot_",datasource,".pdf"),dpi = 300, width = 30, height = 20, units = "cm")

  #   ##GO term tree 
  #   BP <- enrichGO(gene = test$ENTREZID,  #基因列表(转换的ID)
  #                  keyType = "ENTREZID",  #指定的基因ID类型，默认为ENTREZID
  #                  OrgDb=org.Hs.eg.db,  #物种对应的org包
  #                  # universe = names(geneList), #背景基因集  If missing, the all genes listed in the database (eg TERM2GENE table) will be used as background.
  #                  ont = "BP",   #CC细胞组件，MF分子功能，BP生物学过程
  #                  pvalueCutoff = 0.01,  #p值阈值
  #                  pAdjustMethod = "fdr",  #多重假设检验校正方式
  #                  minGSSize = 1,   #注释的最小基因集，默认为10
  #                  maxGSSize = 500,  #注释的最大基因集，默认为500
  #                  qvalueCutoff = 0.05,  #p值阈值
  #                  readable = TRUE)  #基因ID转换为基因名
  #   # BiocManager::install("enrichplot")
  #   library(enrichplot)
  #   p=goplot(BP, showCategory = nshow)
  #   ggsave(plot=p,path = paste0(out.dir), filename = paste0("GO-BP-Tree_",datasource,".pdf"),dpi = 300, width = 30, height = 20, units = "cm")
  #   
  #   ##Disease Ontology (DO) 分析和富集分析，enrichDO函数，主要用于鉴定疾病相关的目标基因，gseDO函数主要用于DO的基因富集分析DOSE也可以用于癌症基因的网络(NCG）以及疾病基因网络(DGN)的富集分析
  #   #install.packages("ggnewscale")
  #   library(enrichplot)
  #   library(ggnewscale)
  #   ed = enrichDO(test$ENTREZID, pvalueCutoff=0.05)
  #   met <- pairwise_termsim(ed)
  #   p=emapplot(met,showCategory = nshow)
  #   ggsave(plot=p,path = paste0(out.dir), filename = paste0("DO_",datasource,".pdf"),dpi = 300, width = 30, height = 20, units = "cm")
  }
  
  out.dir = paste0(outdir_subset_final,'/GOKEGG')
  for (cluster in unique(all.markers.ranked$cluster)){
    if (cluster == 3){
      gene_of_interest = all.markers.ranked[(all.markers.ranked$cluster == cluster & all.markers.ranked$avg_log2FC > 0.5),"gene"]
    }else
    {gene_of_interest = all.markers.ranked[(all.markers.ranked$cluster == cluster & all.markers.ranked$avg_log2FC > 1.25),"gene"]}
    
    GO_Kegg_Enrichment_Vis(gene_of_interest,datasource=paste0("Fib_markers_cluster",cluster),nshow=5,out.dir = out.dir)
    print(paste0("cluster done: ", cluster))
  }
  
  GO_Kegg_Enrichment_Vis(c("FN1","COL1A1","MMP11", "CTHRC1", "COL1A2", 'COL3A1', 'SPARC', 'COL5A2', 'POSTN', 'COL11A1'),
                         datasource = 'top10C0',
                         nshow=10,fromtype="SYMBOL",totype=c("ENSEMBL","ENTREZID"),orgdb="org.Hs.eg.db", out.dir)
  
  #################################################################################
  # ##2.4 find DS of PLAU in epithelial and in Fibroblast ####
  # #ligand-target-receivercell
  # ligand_target_merged_allDS = read.csv(paste0(wd,"/cancer/datafile/outputs_5cancers/new_outputs/MCDM/Ligand_target_allDS_of_allUR/ligand_target_merged_withFC.csv"), sep=",",header=T)
  # head(ligand_target_merged_allDS)
  # ligand_target_merged_allDS_of_scUR = ligand_target_merged_allDS[ligand_target_merged_allDS$ligand %in% scUR$ligand,]
  # ligand_target_merged_scDS_of_scUR = ligand_target_merged_allDS[ligand_target_merged_allDS$ligand %in% scUR,] %>% .[.$target %in% scDS ,]
  # 
  # ligand_target_PLAU_to_Epi = ligand_target_merged_allDS[ligand_target_merged_allDS$ligand=="PLAU",] %>%
  #   .[.$receiver =="Epithelial",]
  # ligand_target_PLAU_to_Epi.wide = dcast(ligand_target_PLAU_to_Epi,target~cancer,value.var="avg_log2FC")
  # write.table(ligand_target_PLAU_to_Epi,file = paste0(outdir_subset_final,"/ligand_target_PLAU_to_Epi.txt"),sep="\t", col.names = T, row.names = F,quote=F)
  # 
  # ligand_target_PLAU_to_Fib = ligand_target_merged_allDS[ligand_target_merged_allDS$ligand=="PLAU",] %>%
  #   .[.$receiver =="Fibroblast",]
  # ligand_target_PLAU_to_Fib.wide = dcast(ligand_target_PLAU_to_Fib,target~cancer,value.var="avg_log2FC")
  # write.table(ligand_target_PLAU_to_Fib,file = paste0(outdir_subset_final,"/ligand_target_PLAU_to_Fib.txt"),sep="\t", col.names = T, row.names = F,quote=F)
  # 
  # ###2.4.1 GO and KEGG for DS ####
  # out.dir = outdir_subset_final
  # gene_of_interest = ligand_target_PLAU_to_Epi$target %>% unique
  # GO_Kegg_Enrichment_Vis(gene_of_interest,datasource="DS_of_PLAU_inEpi",nshow=10)
  # 
  # 
  # gene_of_interest = ligand_target_PLAU_to_Fib$target %>% unique
  # GO_Kegg_Enrichment_Vis(gene_of_interest,datasource="DS_of_PLAU_inFib",nshow=10)
  # 
  # 
  # #check the weight of DS of PLAU in IL-17 pathway
  # g = c('CEBPB' ,'JUN' ,'JUND' ,'FOS' ,'FOSB' ,'FOSL1' ,'HSP90AB1' ,'IFNG' ,'MMP3')
  # ligand_target_PLAU_to_Epi[ligand_target_PLAU_to_Epi$target %in% g,] %>% arrange(-weight)
  # ligand_target_PLAU_to_Fib[ligand_target_PLAU_to_Fib$target %in% g,] %>% arrange(-weight)
  # 
  # intersect(rownames(all.markers.top),ligand_target_PLAU_to_Fib$target %>% unique)
  # intersect(rownames(all.markers.top),ligand_target_PLAU_to_Epi$target %>% unique)
  # 
  # ##2.5 DEG of each cluster in each cancer between tumor vs. normal####
  #   cell="Fibroblast"
  #   load(paste0(outdir_subset_final,"/all_integrated_",cell,".RData"))
  #   all.integrated@meta.data <- all.integrated@meta.data[colnames(all.integrated@assays$RNA@counts),]
  #   DefaultAssay(all.integrated) <- "RNA"
  #   DEG_out.dir <- paste0(outdir_subset_final,"/DEGs_of_",cell)
  #   if (dir.exists(DEG_out.dir) != 1) {dir.create(DEG_out.dir) 
  #     print("create DEG_out.dir") } else { print("DEG_out.dir exsist")}
  #   DEG_method <- "MAST"  #wilcox or MAST和DESeq2
  #   
  #   ###2.5.1 calculate DEGs for each Fib cluster in each cancer ####
  #   for (i in unique(all.integrated$cancer_name)) {
  #     print(i)
  #     Idents(all.integrated) <- "cancer_name"
  #     temp <- subset(all.integrated, idents = i)
  #     for (y in unique(temp$seurat_clusters)){
  #       Idents(temp) <- "seurat_clusters"
  #       temp1 <- subset(temp, idents = y)
  #       Idents(temp1) <- "group"
  #     
  #       if(dim(table(temp1$seurat_clusters,temp1$group))[2] ==1){next} #skip if cell type only exists in one group
  #       if(sum(table(temp1$seurat_clusters,temp1$group)[,1]) < 3){next}
  #       if(sum(table(temp1$seurat_clusters,temp1$group)[,2]) < 3){next}
  #       
  #       iDEGs <- FindMarkers(temp1, ident.1 = "T", ident.2 = "N", verbose = FALSE,
  #                            min.pct = 0.1, #only select genes that expressed in more than 10% of cells 
  #                            test.use = DEG_method #"wilcox"(default), can change to wilcox or MAST和DESeq2
  #       )
  #       head(iDEGs, n = 15)
  #       # iDEGs <- iDEGs %>% filter(p_val_adj <= 0.05)
  #       write.csv(iDEGs,paste0(DEG_out.dir,"/",cell,"_cluster",y,"_",i,"_",DEG_method,".csv"),row.names = T, col.names = T, sep = ',',quote = F)
  #     }
  #   }
  #   
  #  ###2.5.2 calculate DEG for all Fibroblast in each cancer####
  #   load(paste0(outdir_subset_final,"/all_integrated_",cell,".RData"))
  #   all.integrated@meta.data <- all.integrated@meta.data[colnames(all.integrated@assays$RNA@counts),]
  #   DefaultAssay(all.integrated) <- "RNA"
  #   for (i in unique(all.integrated$cancer_name)) {
  #     print(i)
  #     Idents(all.integrated) <- "cancer_name"
  #     temp <- subset(all.integrated, idents = i)
  #     Idents(temp) <- "group"
  #       
  #     if(dim(table(temp$cancer_name,temp$group))[2] ==1){next} #skip if cell type only exists in one group
  #       
  #     iDEGs <- FindMarkers(temp, ident.1 = "T", ident.2 = "N", verbose = FALSE,
  #                            min.pct = 0.1, #only select genes that expressed in more than 10% of cells 
  #                            test.use = DEG_method #"wilcox"(default), can change to wilcox or MAST和DESeq2
  #       )
  #       head(iDEGs, n = 15)
  #       # iDEGs <- iDEGs %>% filter(p_val_adj <= 0.05)
  #       write.csv(iDEGs,paste0(DEG_out.dir,"/",cell,"_",i,"_all_",DEG_method,".csv"),row.names = T, col.names = T, sep = ',',quote = F)
  #     }
    
  
  
#################################
  #3. re-sub merging for Epithelial ####
  ##3.1. subset celltype from each cancers ####
  Epithelial_sublist = list()
  
  outdir = "/home/yelzh67/Projects/cancer/datafile/outputs_5cancers/new_outputs/Seurat_subset_celltyping/2022Feb12"
  file_list= list.files(path=outdir, pattern="all_integrated_Celltype_final_")
  file_list=file_list[file_list!="all_integrated_Celltype_final_nopericyteLung.RData" ] %>% .[. != 'all_integrated_Celltype_final_withpericyteLung.RData']
  cancer_name= sapply(strsplit(file_list,".RData"),'[[',1)
  cancer_name= sapply(strsplit(cancer_name,"_"),'[[',5)
  
  for (i in 1:length(file_list)) {
    load(paste0(outdir,"/",file_list[i]))
    print(i)
    DefaultAssay(all.integrated) <- "RNA"
    seuratObj = all.integrated
    seuratObj@meta.data %>% head()
    
    seuratObj@meta.data$Celltype_final %>% table()
    seuratObj@meta.data$cancer_name = cancer_name[i]
    
    Epithelial_sub <- subset(x = seuratObj, 
                             subset = Celltype_final == "Epithelial")
    
    Epithelial_sublist = c(Epithelial_sublist,Epithelial_sub)
    
    print(paste0(cancer_name[i], " done!"))
    rm(all.integrated,seuratObj)
  }
  
  save(Epithelial_sublist,file = paste0(outdir_subset_final,"/Epithelial_sublist.RData"))
  
  Epithelial_sublist[[1]]$Celltype_final %>% table
  Epithelial_sublist[[1]]@meta.data %>% head()
  
  
  ##3.2. merge the same celltype across cancers - seurat integration####
  library(SeuratDisk)
  library(SeuratWrappers)
  ###3.2.1 select integration features####
  Epithelial_sublist <- lapply(X = Epithelial_sublist, FUN = function(x) {
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })
  # select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = Epithelial_sublist)
  save(features,file = paste0(outdir_subset_final,"/features_Epithelial.RData"))
  
  ###3.2.2 integration  #### 
  load(paste0(outdir_subset_final,"/features_Epithelial.RData"))
  anchors <- FindIntegrationAnchors(object.list = Epithelial_sublist, anchor.features = features)
  all.integrated <- IntegrateData(anchorset = anchors)
  save(all.integrated,file = paste0(outdir_subset_final,"/all_integrated_Epithelial_beforePCA.RData"))
  rm(features)
  
  ###3.2.3 Clustering after integration#### till here
  # specify that we will perform downstream analysis on the corrected data note that the
  # original unmodified data still resides in the 'RNA' assay
  for (i in c("Epithelial")){
    i="Epithelial"
    load(paste0(outdir_subset_final,"/all_integrated_",i,"_beforePCA.RData"))
    all.integrated <- subset(x = all.integrated, 
                             subset = (cancer_name != "nopericyteLung" & cancer_name !="withpericyteLung"))
    all.integrated@meta.data <- all.integrated@meta.data[colnames(all.integrated@assays$RNA@counts),]
    DefaultAssay(all.integrated) <- "integrated"
    seed = 42
    nPCs = 20
    res= 0.2
    reduction= "umap"
    
    all.integrated <- ScaleData(all.integrated, verbose = T)
    all.integrated <- RunPCA(all.integrated,seed.use = seed, verbose = T)
    p <- ElbowPlot(all.integrated,ndims = 30)
    ggsave(plot=p,  path = outdir_subset_final, dpi = 300, width = 25, height = 20, units = "cm" ,filename=paste0(i,"_PCA-ElbowPlot.pdf"))
    all.integrated <- RunUMAP(all.integrated, seed.use = seed, reduction = "pca", dims = seq(nPCs))
    all.integrated <- RunTSNE(all.integrated, seed.use = seed, reduction = "pca", dims=seq(nPCs))
    all.integrated <- FindNeighbors(all.integrated, reduction = "pca", dims = seq(nPCs))
    all.integrated <- FindClusters(all.integrated, random.seed = seed, resolution = res)
    table(all.integrated$seurat_clusters)
    
    DefaultAssay(all.integrated) <- "RNA"
    all.integrated <- NormalizeData(all.integrated)
    all.integrated <- ScaleData(all.integrated, verbose = T)
    all.integrated@assays$RNA@counts[1:70,1:2]
    all.integrated@assays$RNA@data[1:70,1:2]
    all.integrated@assays$RNA@scale.data[1:70,1:2]
    
    all.integrated@assays$integrated@counts[1:70,1:2]
    all.integrated@assays$integrated@data[1:70,1:2]
    all.integrated@assays$integrated@scale.data[1:70,1:2]
    
    save(all.integrated, file = paste0(outdir_subset_final,"/all_integrated_",i,".RData"))
  }
  
  for (i in c("Epithelial")){
    i="Epithelial"
    load(paste0(outdir_subset_final,"/all_integrated_",i,".RData"))
    all.integrated@meta.data <- all.integrated@meta.data[colnames(all.integrated@assays$RNA@counts),]
    DefaultAssay(all.integrated) <- "integrated"
    seed = 42
    nPCs = 20
    res= 0.2
    reduction= "umap"
    
    p1 <- DimPlot(all.integrated, reduction = reduction, group.by = "group")
    p2 <- DimPlot(all.integrated, reduction = reduction, label = TRUE)
    # p3 <- DimPlot(all.integrated, reduction = reduction, split.by = "group") #To visualize the two conditions side-by-side, we can use the split.by argument to show each condition colored by cluster.
    p3 <- DimPlot(all.integrated, reduction = reduction, group.by = "cancer_name")
    p4 <- DimPlot(all.integrated, reduction = reduction, split.by = "group" )
    p =  (p1 + p2 + p3)
    ggsave(plot=p,  path = outdir_subset_final, width = 9, height = 3,filename=paste0("/DimPlot_celltypefinal_",i,"_",reduction,"_seed",seed,"_res",res,".pdf"))
    
    p=DotPlot(object = all.integrated, features = c("PLAU",'PLAUR','IL6','ACTA2','COL11A1','FN1','COL1A1','SPP1', 'COL4A1', 'COL18A1',  'CLEC11A', 'MDK','MMP3','IL17A','FOS','JUN','MMP11','CTHRC1'), assay="RNA",split.by = "group")  
    ggsave(plot=p,  path = outdir_subset_final, device = "pdf", width = 12, height = 10 ,filename=paste0("DotPlot","_",i,"_bygroup_",reduction,"_seed",seed,"_res",res,".pdf"))
    p=DotPlot(object = all.integrated, features = c("PLAU",'PLAUR','IL6','ACTA2','COL11A1','FN1','COL1A1',  'SPP1', 'COL4A1', 'COL18A1',  'CLEC11A', 'MDK','MMP3','IL17A','FOS','JUN','MMP11','CTHRC1'), assay="RNA")+RotatedAxis()
    ggsave(plot=p,  path = outdir_subset_final, device = "pdf", width = 12, height = 9 ,filename=paste0("DotPlot","_",i,"_",reduction,"_seed",seed,"_res",res,".pdf"))
    
    
    DefaultAssay(all.integrated) <- "RNA"
    p <- FeaturePlot(all.integrated, features = c("PLAU",'IL6','ACTA2'),split.by = "group", max.cutoff = 3, label = TRUE, cols = c("grey", "red"))
    ggsave(plot=p,  path = outdir_subset_final, device = "pdf", width = 6, height = 9 ,filename=paste0("Feature_Plot","_",i,"_bygroup_",reduction,"_seed",seed,"_res",res,".pdf"))
    
    p <- FeaturePlot(all.integrated, features = c("PLAU",'IL6','ACTA2'),split.by = "cancer_name", max.cutoff = 3, label = TRUE)
    ggsave(plot=p,  path = outdir_subset_final, device = "pdf",width = 15, height = 9 ,filename=paste0("Feature_Plot","_",i,"_bycancer_",reduction,"_seed",seed,"_res",res,".pdf"))
  }
  
  
  
  
  
  
  
  
  ###############END###########
  
  
  
  
  
  
  
  
  
  
  
  
    
  
  
  
  
  
  
  
  

##2.3 identify PLAU high ####
i="Fibroblast"
# load(paste0(outdir_subset_final,"/all_integrated_",i,".RData"))
library(ggpubr)
all.integrated$hPLAU = all.integrated@assays$RNA@data["PLAU",]
all.integrated$hPLAU[all.integrated$hPLAU == 0] <- NA
q = all.integrated$hPLAU %>% quantile(.,probs=c(0.5),na.rm=T)
all.integrated$hPLAU[all.integrated$hPLAU >= q] <- "high"
all.integrated$hPLAU[all.integrated$hPLAU < q] <- "low"
table(all.integrated$hPLAU,all.integrated$group)

p <- FeaturePlot(all.integrated, features = c("PLAU",'IL6','ACTA2','MDK','SPP1'),split.by = "hPLAU", max.cutoff = 3, label = TRUE, cols = c("grey", "red"))
ggsave(plot=p,  path = outdir_subset_final, device = "pdf",width = 6, height = 12 ,filename=paste0("Feature_Plot","_",i,"_byhPLAU_",reduction,"_seed",seed,"_res",res,".pdf"))

# p <- DotPlot(all.integrated, features = unique(genes_to_check),
#              assay='RNA' , ##can use this to check orignial or normalized expression
#              group.by = 'Celltype_final' ##this can change to other group method like 'seurat_clusters'
# )  + scale_color_viridis_c()  + coord_flip() 
# ggsave(plot=p,  path = outdir_subset_PCA, dpi = 300, width = 30, height = 45, units = "cm" ,filename=paste0("check_markers_celltypefinal_",i,"_seed",seed,"_res",res,".pdf"))
# 

all.integrated$PLAU = all.integrated@assays$RNA@data["PLAU",]
all.integrated$PLAU[all.integrated$PLAU == 0] <- NA
all.integrated$IL6 = all.integrated@assays$RNA@data["IL6",]
all.integrated$IL6[all.integrated$IL6 == 0] <- NA
all.integrated$ACTA2 = all.integrated@assays$RNA@data["ACTA2",]
all.integrated$ACTA2[all.integrated$ACTA2 == 0] <- NA

p1 = ggboxplot(all.integrated@meta.data,x = "cancer_name", y = "IL6", color = "hPLAU",
          short.panel.labs = F,palette = "jco") +
  stat_compare_means(label="p.format",aes(group=hPLAU)) #t.test
p2 = ggboxplot(all.integrated@meta.data,x = "cancer_name", y = "ACTA2", color = "hPLAU",
               short.panel.labs = F,palette = "jco") +
  stat_compare_means(label="p.format",aes(group=hPLAU)) #t.test
p1+p2
pdf(file=paste0(outdir_subset_final,"/Boxplot_bycancer_plau_IL6_ACTA2.pdf"), width = 10, height = 5)
print(p1+p2)
dev.off()

p1 = ggboxplot(all.integrated@meta.data,x = "hPLAU", y = "IL6", color = "hPLAU",
               short.panel.labs = F,palette = "jco") +
  stat_compare_means(label="p.format",aes(group=hPLAU)) #t.test
p2 = ggboxplot(all.integrated@meta.data,x = "hPLAU", y = "ACTA2", color = "hPLAU",
               short.panel.labs = F,palette = "jco") +
  stat_compare_means(label="p.format",aes(group=hPLAU)) #t.test
p1+p2
pdf(file=paste0(outdir_subset_final,"/Boxplot_plau_IL6_ACTA2.pdf"), width = 8, height = 5)
print(p1+p2)
dev.off()


data = all.integrated@meta.data
p1=ggplot(data, aes(x=PLAU, y=IL6)) +
  geom_point(alpha=0.7,aes(color=cancer_name,shape=group)) +
  # geom_cor(method = "spearman",use = "complete.obs") +
  stat_smooth(method = "lm", #glm,gam,loess: default value,rim
              col = "gray",
              se = FALSE,
              size = 0.5)  +
  #scale_y_reverse(limits = c(12, 0))+
  # geom_text_repel(aes(label = cancer_name),size=3 ) +
  # annotate("text",x=2,y=1,label="atop(Spearman_r==0.508,p==`0.0007`)",parse=T) +
  labs(x='PLAU',y='IL6') +
  theme_classic()
p2=ggplot(data, aes(x=PLAU, y=ACTA2)) +
  geom_point(alpha=0.7,aes(color=cancer_name,shape=group)) +
  # geom_cor(method = "spearman",use = "complete.obs") +
  stat_smooth(method = "lm", #glm,gam,loess: default value,rim
              col = "gray",
              se = FALSE,
              size = 0.5)  +
  #scale_y_reverse(limits = c(12, 0))+
  # geom_text_repel(aes(label = cancer_name),size=3 ) +
  # annotate("text",x=2,y=1,label="atop(Spearman_r==0.508,p==`0.0007`)",parse=T) +
  labs(x='PLAU',y='ACTA2') +
  theme_classic()
pdf(file=paste0(outdir_subset_final,"/Corr_plau_IL6_ACTA2.pdf"), width = 8, height =4)
print(p1+p2)
dev.off()


## 2.4 find marker genes ####
for (i in c("Fibroblast")){ #findallmarkers
  load(paste0(outdir_subset_final,"/all_integrated_",i,".RData"))
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










