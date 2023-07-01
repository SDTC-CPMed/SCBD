# Identify unique URs for different tissue of origin
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
library(patchwork)
library(RColorBrewer)
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
library(patchwork)
library(RColorBrewer)
wd = '/Users/yelin.zhao/Library/CloudStorage/OneDrive-Personal/Github/SCBD'

input.dir=paste0(wd,"/Data")
outdir= paste0(wd,"/Output/MCTM/uniqMCTM")
if (dir.exists(outdir)==F){dir.create(outdir)
  print("outdir created")}
cancer_name = c('Lung','Liver','Ovary','Colon','Breast')
#########################
# find unique UR####
ligand_list = list()
for (cancer in cancer_name){
  ligands <- read.csv(paste0(wd,"/Output/MCTM/",cancer,"/PCC_positive_ligand.txt"), sep="\t",header= T)
  ligand_list = c(ligand_list,list(ligands$test_ligand %>% unique ) )
}
names(ligand_list) = cancer_name
ligand_list = c(ligand_list['Breast'],ligand_list['Colon'],ligand_list['Liver'],ligand_list['Lung'],ligand_list['Ovary'])

Breast_unique <- setdiff(ligand_list$Breast, unlist(ligand_list[c("Colon", "Liver", "Lung", "Ovary")]))
Colon_unique <- setdiff(ligand_list$Colon, unlist(ligand_list[c("Breast", "Liver", "Lung", "Ovary")]))
Liver_unique <- setdiff(ligand_list$Liver, unlist(ligand_list[c("Breast", "Colon", "Lung", "Ovary")]))
Lung_unique <- setdiff(ligand_list$Lung, unlist(ligand_list[c("Breast", "Liver", "Colon", "Ovary")]))
Ovary_unique <- setdiff(ligand_list$Ovary, unlist(ligand_list[c("Breast", "Liver", "Lung", "Colon")]))

merged_unique = data.frame()
for (cancer in cancer_name){
  i_cancer = get(paste0(cancer,'_unique')) %>% as.data.frame()
  i_cancer$cancer = cancer
  merged_unique = rbind(merged_unique,i_cancer)
}
colnames(merged_unique) = c('uniqueUR','cancer')

# check whether they encode secreted proteins
secre_HPA = read.csv(file = '/Users/yelin.zhao/Library/CloudStorage/OneDrive-Personal/Projects/cancer/datafile/outputs_5cancers/new_outputs/Biomarker_candidate _criterial/proteinatlas_download_2022Feb16_secreated.txt',sep='\t')
merged_unique = merged_unique %>% mutate(UR_secreated = ifelse(uniqueUR %in% secre_HPA$Gene,'Yes','No'))
table(merged_unique$cancer,merged_unique$UR_secreated)

# check which cell type contains the unique ligands
ligand_target_merged_FC = read.csv('/Users/yelin.zhao/Library/CloudStorage/OneDrive-Personal/Github/SCBD/Output/MCTM/shMCTM/ligand_target_merged_FC.csv',row.names = 1)
merged_unique_UR_FC = merge(merged_unique,ligand_target_merged_FC,by.x=c('uniqueUR','cancer'),by.y=c('ligand','cancer'),all.x=T)
merged_unique_UR_FC = (merged_unique_UR_FC %>% select(-weight) )%>% unique
write.csv(merged_unique_UR_FC,file=paste0(outdir,'/merged_unique_UR_FC.csv'))
merged_unique_UR_celltype = (merged_unique_UR_FC %>% select(-receiver, -target, -avg_log2FC.target, -avg_log2FC.ligand) )%>% unique
write.csv(merged_unique_UR_celltype,file=paste0(outdir,'/merged_unique_UR.csv'))

# check the FC of unique secreated URs
merged_unique_secreated = merged_unique_UR_FC[merged_unique_UR_FC$UR_secreated== 'Yes',] %>% select(-UR_secreated,-receiver, -target, -avg_log2FC.target) %>% unique
write.csv(merged_unique_secreated,file=paste0(outdir,'/merged_unique_secreated.csv'))
p = ggplot(merged_unique_secreated, aes(x = cancer, y = avg_log2FC.ligand, color = sender)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = uniqueUR), size = 2, hjust = 1.2)
  labs(title="log2FC of unique URs", y="log2FC")
pdf(file=paste0(outdir,"/log2FC_of_unique URs.pdf"), width = 7, height = 7)
print(p)
dev.off()

ggplot(merged_unique_secreated, aes(x = cancer, y = avg_log2FC.ligand, color = sender)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = uniqueUR), size = 2, hjust = 1.2)














