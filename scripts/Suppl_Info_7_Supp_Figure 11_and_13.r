
#This script has the code for Supp
#Supplementary Information 7.
#Supplementary Figure 11 and #Supplementary Figure 13
#-----------------------------------
#Load libraries
# Load libraries
suppressMessages(require(Seurat))
suppressMessages(require(tidyverse))
suppressMessages(require(scales))
suppressMessages(require(cowplot))
suppressMessages(require(RCurl))
suppressMessages(require(Matrix))
suppressMessages(require(scCustomize))
suppressMessages(require(viridis))
suppressMessages(require("MAST"))
#-------------------------------


#Supplementary Figure 11

#Single cell RNA-seq clusters of cell types in murine caecal tumours
#----Plot Markers for the Main cell type-------
#plot the scaled expression across all cell types
# epithelial cells (Epcam), cancer associated fibroblast (CAF) (DCN, Col1a1), 
#T cells (Cd3e, Cd3g), Macrophages (Cd74, H2-Ab1, lyz2),
# Neutrophils (S100a8 and S100a9).

seurat_integrated <- readRDS("results/integrated_seurat_Resolution.rds")

DefaultAssay(seurat_integrated) <- "RNA"
Idents(object = seurat_integrated) <- "integrated_snn_res.1.2"
#Another marker for marcrophages is , "lyz2"
gene.list  <- c("Epcam","Dcn", "Col1a1", "Cd3e", "Cd3g", "Cd74", "H2-Ab1", "S100a8", "S100a9")

FeaturePlot_scCustom(seurat_integrated, features = gene.list, 
        label = FALSE, order = TRUE, reduction = "umap", figure_plot = TRUE, min.cutoff = 'q10')


 

