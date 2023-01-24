#### looping for ligands activity ####
# cluster in all_expressed_genes which has DEGs as sender
# cluster in DEG_human as receiver
# genes in each cluster in DEG_human as gene_oi 
# background genes are all expressed genes in the cluster which is defined as receiver in the current loop
starttime<-Sys.time()

library(nichenetr) 
library(Matrix)
library(dplyr)

all_ligand_activity_Yelin <- data.frame(matrix(ncol=6,nrow=0))
colnames(all_ligand_activity_Yelin) <- c("test_ligand","auroc","aupr","pearson", "Sender", "Target")

best_upstream_ligands_all_interaction <- data.frame(matrix(ncol=3,nrow=0))
colnames(best_upstream_ligands_all_interaction) <- c("test_ligand", "Sender", "Target")

for(i in 1:dim(DEGs_human)[2]){ 

  for(j in 1:dim(DEGs_human)[2]){
    
    # if (j==7) { #One time use, because hv10000_singleR_fine cluster7 doesn't work as receiver cluster
    #   next
    # }
    sender_cluster = colnames(DEGs_human)[i]
    receiver_cluster = colnames(DEGs_human)[j]   
   
    if (receiver_cluster %in% colnames(DEGs_human) ==F){next}
    if (receiver_cluster %in% colnames(all_expressed_genes) ==F){next}
    
    expressed_genes_sender = as.vector(DEGs_human[,i])
    expressed_genes_receiver = as.vector(all_expressed_genes[,receiver_cluster])
    
    # gene set of interest - set DEG as gene set of interest ####
    geneset_oi <- as.vector(DEGs_human[, j] %>% .[. %in% rownames(ligand_target_matrix)] )
    
    # set background expressed genes  ####
    background_expressed_genes = as.vector(expressed_genes_receiver %>%
                                             .[. %in% rownames(ligand_target_matrix)])
    
    # Define a set of potential ligands ####
    # Putative ligand-receptor links were gathered from NicheNet’s ligand-receptor data sources.
    # lr_network = readRDS("lr_network.rds")
    # If wanted, users can remove ligand-receptor interactions that were predicted based on protein-protein interactions 
    # and only keep ligand-receptor interactions that are described in curated databases. 
    # To do this: uncomment following line of code:
    # lr_network = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
    
    ligands = lr_network %>% pull(from) %>% unique()        
    expressed_ligands = intersect(ligands,expressed_genes_sender)  
    
    if (length(expressed_ligands)==0) {
      next
    }
    
    receptors = lr_network %>% pull(to) %>% unique()
    expressed_receptors = intersect(receptors,expressed_genes_receiver)
    if (length(expressed_receptors)==0) {
      next
    }
    
    lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors)      #选取ligand 和receptor都要有表达的
    head(lr_network_expressed)
    
    potential_ligands = lr_network_expressed %>% pull(from) %>% unique()  #只有lr_network在当前数据中有表达的ligands才考虑作为potentially active ligands for the NicheNet analysis
    head(potential_ligands)
    
    if (length(potential_ligands)[1]==0){
      next}
    
    ligand_activities = predict_ligand_activities(geneset = geneset_oi, 
                                                  background_expressed_genes = background_expressed_genes, 
                                                  ligand_target_matrix = ligand_target_matrix, 
                                                  potential_ligands = potential_ligands)
    
    ligand_activities %>% arrange(-pearson)
    best_upstream_ligands_interaction <- ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
    
    best_upstream_ligands_interaction <- cbind(best_upstream_ligands_interaction, Sender=sender_cluster, Target=receiver_cluster)
    best_upstream_ligands_all_interaction <- rbind(best_upstream_ligands_all_interaction, best_upstream_ligands_interaction)
    
    ligand_activities <- cbind(ligand_activities,Sender=sender_cluster, Target=receiver_cluster)
    all_ligand_activity_Yelin <- rbind(all_ligand_activity_Yelin,ligand_activities)
    
  }
}
 
endtime <-Sys.time()
print(starttime)
print(endtime)



