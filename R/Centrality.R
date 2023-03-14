library(reshape2)

input.dir=paste0(wd,"/Data")
outdir= paste0(wd,"/Output/MCTM")
if (dir.exists(outdir)==F){dir.create(outdir)
  print("outdir created")}
cancer_name = c("Lung","Liver","Colon","Ovary","Breast")

# Centrality heatmap  
centrality_merge <- matrix(NA, ncol=10) %>% as.data.frame()
colnames(centrality_merge) = c("node_degree_all","node_degree_in","node_degree_out","closeness","eigenvector centralities" , "Kleinberg's hub centrality scores", "Laplacian Centrality","Leverage Centrality","cancer","celltype")
for (cancer in cancer_name){
  a = read.csv(paste0(outdir,"/",cancer,"/Cell-cell_centrality_summary_PCC_positive.txt"),sep = "\t",row.names = 1) %>% as.data.frame()
  colnames(a) = lapply(colnames(a), function(x) {strsplit(x,"_") %>% unlist %>% .[2]})  %>% unlist  
  a = t(a)  %>% as.data.frame()
  a$cancer = cancer
  a$celltype = rownames(a)
  centrality_merge = rbind(centrality_merge,a)
}
centrality_merge = centrality_merge[2:dim(centrality_merge)[1],]


df = centrality_merge
head(df)
df.wide = dcast(df,celltype~cancer, value.var ="eigenvector centralities")
rowname = df.wide$celltype
df.wide = df.wide[,2:dim(df.wide)[2]]
df.wide[1:5,1:4]
df.wide[is.na(df.wide)==T] <- 0
df.wide = apply(df.wide,2,as.numeric)
rownames(df.wide) = rowname

p = Heatmap(t(df.wide),
            name = "Eigenvector centrality", #title of legend
            column_title = "Cell type",
            column_title_side = "bottom",   
            column_title_gp = gpar(fontsize = 15, fontface = "bold"),  
            column_names_gp = gpar(fontsize = 15),
            row_title = "Cancers",
            row_title_side ="right",   
            row_title_gp = gpar(fontsize = 15) ,
            row_names_gp = gpar(fontsize = 15), # Text size for row names
            cluster_columns = T,
            cluster_rows = T,
            col = colorRamp2(c(0,1), c("white", "red")),
            border = T
) 
pdf(file=paste0(outdir,"/heatmap_eigenvector_centrality.pdf"), width = 10, height = 5)
print(p)
dev.off()

