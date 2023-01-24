#InferCNV
# module load R/4.0.4
# library("devtools")
# devtools::install_github("broadinstitute/infercnv")
rm(list=ls())
myPaths <- .libPaths() 
# myPaths <- c(myPaths, "/home/yelzh67/R/x86_64-pc-linux-gnu-library/4.0") # add new path
myPaths <- c(myPaths[2], myPaths[1])  # switch them
.libPaths(myPaths)  

# remotes::install_version("Seurat", version = "4.0.4")
library(Seurat)
library(clusterProfiler)  
library(DOSE)
library(topGO)
library(pathview)
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(infercnv)
# install.packages('AnnoProbe')
library(AnnoProbe)


location = "omika"
cancer_name = c("Lung","Liver","Colon","Ovary","Breast")

if (location == "mac"){  wd="/Users/cynthia_ye/Documents/OneDrive - Linköpings universitet"}
if (location == "omika"){ wd="/home/yelzh67/Projects"}

inputdir = paste0(wd,"/cancer/datafile/outputs_5cancers/new_outputs/Seurat_subset_celltyping/2022Feb12")
outdir= paste0(wd,"/cancer/datafile/outputs_5cancers/new_outputs/InferCNV/20220803-onlyepi")
if (dir.exists(outdir)==F){dir.create(outdir)
  print("outdir created")} else {print("outdir existed")}

################################################################
################################################################
i=1 done
i=2 
i=3 done
i=4 done
4:length(cancer_name)

for (i in c(2,5)){
  out.dir= paste0(outdir,"/",cancer_name[i])
  if (dir.exists(out.dir)==F){dir.create(out.dir)
    print("out.dir created")} else {print("out.dir existed")}
  
  #1. Prepare files####
  load(paste0(inputdir,"/all_integrated_Celltype_final_",cancer_name[i],".RData"))
  all.integrated@meta.data %>% head
  table(all.integrated@meta.data$Celltype_final)
  all.integrated@meta.data$celltypeforInferCNV = paste0(all.integrated@meta.data$group,"_",all.integrated@meta.data$Celltype_final)
  table(all.integrated@meta.data$celltypeforInferCNV)
  
  all.integrated =  subset(x = all.integrated, 
                           subset = celltypeforInferCNV %in% c("T_Epithelial","N_Epithelial")) #,"N_Epithelial","N_Fibroblast","T_Fibroblast","T_Endothelial"
  ## random subset 5000 cells if nrow(all.integrated@meta.data) >5000
  if (nrow(all.integrated@meta.data) >5000) {
    print("N cells >5000, randomly subset 5000 cells")
    all.integrated = subset(all.integrated, downsample = 5000)
  }
  ## get the matrix file
  matrix_file = as.data.frame(GetAssayData(all.integrated , slot='counts',assay='RNA'))
  matrix_file[1:5,1:5]
  
  groupinfo = data.frame(v1=colnames(all.integrated),
                         v2=all.integrated@meta.data$celltypeforInferCNV)
  table(groupinfo$v2,useNA =  "always")
  groupinfo[is.na(groupinfo$v2),]
  identical(colnames(matrix_file),groupinfo$v1)
  
  geneInfor=annoGene(rownames(all.integrated),"SYMBOL",'human')
  colnames(geneInfor)
  geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
  geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
  length(unique(geneInfor[,1]))
  head(geneInfor)
  dim(geneInfor)
  which(is.na(geneInfor$chr))
  geneInfor = geneInfor[is.na(geneInfor$chr) == FALSE,] #remove gene which the chromosomes is not properly assigned
  geneInfor = geneInfor[(geneInfor$chr == "chrM") == FALSE,] #remove gene which the chromosomesM
  geneInfor = geneInfor[(geneInfor$chr == "chrX") == FALSE,] #remove gene which the chromosomesX
  geneInfor = geneInfor[(geneInfor$chr == "chrY") == FALSE,] #remove gene which the chromosomesY
  
  # Make sure that chromosomes are ordered correctly 
  geneInfor$chr = factor(geneInfor$chr, 
                       levels = c("chr1", "chr2","chr3","chr4", "chr5", "chr6","chr7", "chr8", "chr9","chr10", "chr11", "chr12","chr13", "chr14", "chr15","chr16", "chr17", "chr18","chr19", "chr20", "chr21","chr22"))
  geneInfor = geneInfor[order(geneInfor$chr),]
  
  expFile=paste0(out.dir,'/expFile.txt')
  write.table(matrix_file,file = expFile,sep = '\t',quote = F,col.names = T,row.names = T)
  groupFiles=paste0(out.dir,'/groupFiles.txt')
  write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
  geneFile=paste0(out.dir,'/geneFile.txt')
  write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)
  
  ################################################################
  #2. Run InferCNV####
  options(stringsAsFactors = F)
  expFile=paste0(out.dir,'/expFile.txt')
  groupFiles=paste0(out.dir,'/groupFiles.txt')
  geneFile=paste0(out.dir,'/geneFile.txt')
  
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,
                                      annotations_file=groupFiles,
                                      delim="\t",
                                      gene_order_file= geneFile,
                                      ref_group_names=c("N_Epithelial")) #"N_Epithelial(non-malignant)",'T_Endothelial','T_Fibroblast',,"T_Fibroblast","T_Endothelial"
  unique(infercnv_obj@gene_order$chr)
  table(infercnv_obj@gene_order$chr)

  # Make sure that chrmosomes are ordered correctly 
  slot(infercnv_obj, "gene_order")[,"chr"] <- factor(slot(infercnv_obj, "gene_order")[,"chr"], 
                                                     levels = c("chr1", "chr2","chr3","chr4", "chr5", "chr6","chr7", "chr8", "chr9","chr10", "chr11", "chr12","chr13", "chr14", "chr15","chr16", "chr17", "chr18","chr19", "chr20", "chr21","chr22","chrX"))
   
  # perform infercnv operations to reveal cnv signal
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                               out_dir=out.dir,  # dir is auto-created for storing outputs
                               cluster_by_groups=T,   # perform separate clustering for the tumor cells according to the patient type, as defined in the cell annotations file.
                               denoise=T,
                               cluster_references = T,
                               analysis_mode = "subclusters",
                               HMM=T,
                               k_obs_groups=5 ,
                               num_threads=20,
                               output_format="pdf")
  
  # infercnv_obj = readRDS('preliminary.infercnv_obj')
  # infercnv::plot_cnv(infercnv_obj, #上两步得到的infercnv对象
  #                    plot_chr_scale = T, #画染色体全长，默认只画出（分析用到的）基因
  #                    out_dir = out.dir,
  #                    output_filename = "better_plot",output_format = "pdf" ) #保存为pdf文件
  #                    # custom_color_pal =  color.palette(c("#8DD3C7","white","#BC80BD"), c(2, 2)) #改颜色
  # 
}




################################################################