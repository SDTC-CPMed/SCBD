
# A scalable, single-cell based strategy for prioritization of multi-omics pan-cancer biomarkers : **S**hared **C**ancer **B**iomarker **D**iscovery (SCBD)

**S**wedish **D**igital **T**win **C**onsortium **-** **C**enter for **P**ersonalized **M**edicine

## General information

The background problem is early cancer diagnosis improves survival, but is complicated by lack of symptoms until late stages. Multicancer early detection (MCED) tests may greatly improve cancer prognosis. Existing and emerging technologies for multi-omics analyses have shown that prioritization of biomarkers for an MCED is complicated by the involvement of thousands of genes, transcripts, proteins, and other types of molecules, which may vary not only between malignant cells, but also among other cell types in each cancer. Development of a clinically tractable MCED therefore calls for a systematic, and explainable, strategy to prioritize biomarkers that are shared across multiple cell types in different cancers. Here, we propose such a strategy, which we call Shared Cancer Biomarker Discovery (SCBD).

We first constructed multicellular tumor models of each cancer, as well as their shared characteristics (MCTMs, and shMCTM, respectively). These showed predicted interactions between upstream regulators in each cell type and their downstream genes in other cell types. The interactions were multi-directional, and showed no evident hierarchical organization in which tumor cells had a superior role. Ranking of cell types based on network tools and GWAS enrichment showed that tumor cells and fibroblasts had the relatively most important pathogenic roles. These analyzes were used to systematically prioritize genetic, mRNA and protein biomarkers based on functional understanding of underlying mechanisms, and their relative importance. Machine learning was applied to validate, optimal combinations of these biomarkers for MCED in independent GWAS, mRNA or protein profiling data from multiple cancers. 

We propose SCBD as a scalable, and explainable, strategy to prioritize biomarkers from existing and emerging multi-omics data, which may be used for MCED and prognosis. 


## Overview
![Figure 1.pdf](https://github.com/SDTC-CPMed/SCBD/files/11925299/Figure.1.pdf)

**Overview figure of the study**. A) Multicellular tumor models (MCTMs) were constructed based on scRNA-seq data of five common cancers.  MCTMs showed predicted Upstream Regulators (URs) in any cell type and their Downstream Genes in any other cell types (lines and dots, respectively). A shared MCTM (shMCTM) was extracted from the five MCTMs; Each dot represents one cell type, which is connected by predicted molecular interactions, the directions of which are indicated by arrows. B) Prioritization of shared upstream regulator genes (shURs), cell types and subtypes in the shMCTM led to the identification of potential genetic, mRNA or protein MCED biomarkers. C) The diagnostic and prognostic ability of those were analyzed using GWAS (12 cancers), transcriptomic, and proteomic data (9 cancers).

## Multicellular tumor models (MCTMs) and shared MCTM (shMCTM)

After scRNA-seq data had been denoised, clustered and differentially expressed genes (DEGs) had been calculated we applied NicheNet (https://github.com/saeyslab/nichenetr) to select ligand-target interactions between cell types that were predictive of the transcriptomic perturbation observed in the downstream cell type. This allowed creation of a directed multicellular disease model which reflected the altered information flow in disease. Using centrality in the MCDM we were able to rank cell types by their relative importance, which correlated well with the significance of GWAS enrichment among the DEGs of a cell type and the prediction precision for disease-relevant drugs. 

In this project, we created MCTM for each individual cancer dataset from five different cancers. Next, we compared these five MCTMs to identify shMCTM. 

### Prepare data file for creating MCTM, shMCTM, prioritize cell type

Before we start to create MCTMs, several files were needed as input files:

1. Differentially expressed genes (DEGs) between tumor vs. normal within each cell type. 

The DEG files prepared as below:

- A dataframe that contains gene name of DEGs within each cell type (one file per cancer).
<img width="898" alt="image" src="https://user-images.githubusercontent.com/98571115/206129210-e804f115-d3fc-4176-8140-70c860eb6ec9.png">

- A dataframe that contains logFC for all DEGs with each cell type/cancer (one file for all cancers).
<img width="425" alt="image" src="https://user-images.githubusercontent.com/98571115/206133430-a8fdbdef-f5ac-4079-8ac2-7e2f38554722.png">

2. Background genes within each cell type.

In this case study, background genes were defined as genes expressed in more than 1% of cells in this cell type. However, the percentage of cells can be change to other value based on each project's need. Background gene list organized as below (one file per cancer):
<img width="966" alt="image" src="https://user-images.githubusercontent.com/98571115/206129419-0eb174aa-5028-48b7-8136-9f9554952ca2.png">

3. Marker genes calculated for each cluster against all other clusters.
<img width="394" alt="image" src="https://user-images.githubusercontent.com/98571115/206136431-eef41185-6bf6-4dba-b37a-2cb1b44541d1.png">



## Environemental set-up

We have used R v4.0.4 and Python v3.8.3 for the development of the codes.

***R packages:***

Seurat: 4.1.0

SingleR: 1.4.1

nichenetr: 1.1.0

CINNA: 1.2.0

MAST: 1.16.0 

clusterProfiler: 3.18.1 

dplyr: 1.0.8  

tidyverse_1.3.1

igraph: 1.3.0

circlize: 0.4.14

ComplexHeatmap: 2.11.1

ggplot: 2_3.3.6 

patchwork: 1.1.2             
 

***Python:***

sklearn: 0.23.1

scipy: 1.9.3

pandas: 1.5.1

numpy: 1.23.4

matplotlib: 3.2.2








