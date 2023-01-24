

##########################################################################################################
# Calculation of cell-cell centralities
#
# Assuming NicheNet_analysis.R has been run
##########################################################################################################
#install.packages("igraph")
#install.packages("CINNA")

starttime<-Sys.time()

library(igraph)
library(CINNA)
 


#########################################################
# Input
#########################################################
# LOAD NicheNet ligand interactions
all_ligand_activity_Yelin <- PCC_positive_ligand
# all_ligand_activity_Yelin <-  read.table(file = "/Users/cynthia_ye/Documents/OneDrive - Linköpings universitet/cancer/datafile/outputs_5cancers/new_outputs/Nichenet/Liver/all_ligand_activity.csv", header= T, sep=",",stringsAsFactor = F)

ligands_interaction <- all_ligand_activity_Yelin # cell-cell interactions + ligand information
ligands_interaction <- as.matrix(ligands_interaction)
head(ligands_interaction)
# Create graph
g <- graph_from_edgelist(el = ligands_interaction[,5:6], directed = T)
                 
out <- matrix(NA, nrow = 4, ncol = length(V(g)$name))

rownames(out) <- c("node_degree_all", "node_degree_in", "node_degree_out", "closeness")
colnames(out) <- paste("Cluster_", V(g)$name, sep="")

# Calculate centrality degree
out[1,] <- centr_degree(g, mode = "all")$res # all
out[2,] <- centr_degree(g, mode = "in")$res # in degree
out[3,] <- centr_degree(g, mode = "out")$res # out degree

# Calculate closeness
out[4,] <- closeness(g, mode = "all")

# Calculate centralities
#pr_cent<-proper_centralities(g)
#calculate_centralities(g, include = pr_cent[1:20])  %>% pca_centralities(scale.unit = TRUE)
proper_centralities(g)
centrality_matrix <- calculate_centralities(g, include = c("eigenvector centralities", "Kleinberg's hub centrality scores", 
                                                           "Laplacian Centrality", "Leverage Centrality", "Group Centrality"))
p = pca_centralities( centrality_matrix  ) #figure out the order of most important centrality types based on your graph structure
ggsave(plot=p,  path = outputdir, dpi = 300, width = 25, height = 20, units = "cm" ,filename="PCA_centrality_matrix.pdf")

tsne_centralities( centrality_matrix, dims = 2, perplexity = 1, scale=TRUE) #Another method for distinguishing which centrality measure has more information or in another words has more costs is using (t-SNE) t-Distributed Stochastic Neighbor Embedding analysis(Van Der Maaten 2014). 
p = visualize_heatmap( centrality_matrix , scale = TRUE  )
ggsave(plot=p,  path = outputdir, dpi = 300, width = 25, height = 20, units = "cm" ,filename="centrality_matrix.pdf")


centrality_matrix <- matrix(unlist(centrality_matrix), ncol = length(V(g)$name), byrow = T)
colnames(centrality_matrix) <- paste("Cluster_", V(g)$name, sep="")
rownames(centrality_matrix) <- c("eigenvector centralities", "Kleinberg's hub centrality scores", 
                                 "Laplacian Centrality", "Leverage Centrality", "Group Centrality")

out <- rbind(out, centrality_matrix)

write.table(out, file=paste0(outputdir,"/Cell-cell_centrality_summary_PCC_positive.txt"), sep="\t", col.names = NA, row.names = T)


endtime <-Sys.time()
print(starttime) 
print(endtime)


