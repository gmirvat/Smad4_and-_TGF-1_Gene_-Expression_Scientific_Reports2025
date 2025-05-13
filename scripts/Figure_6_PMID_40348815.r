
#In this script I will show the code used in 
# PMID: 40348815  DOI: 10.1038/s41598-025-00908-4
#6a. Cluster tree
#6b. Heatmap
#6c. dotplots for different cell types
#6d. stack plot for the celll type See if the proportion of the cell type differ between the genotype stack plot for the celll type
#6e. UMAP plot for selected genes 
#6f. Violin plot for pseudobulk RNA expression
#author: "Mirvat"

setwd("~/Documents/Analysis/") # set it to your preferred directory

# Load libraries
#---------------------------
suppressMessages(require(Seurat))
suppressMessages(require(tidyverse))
suppressMessages(require(scales))
suppressMessages(require(cowplot))
suppressMessages(require(RCurl))
library('irlba')
library(Matrix)
library(scCustomize)
library(readxl)
#devtools::install_github("rpolicastro/scProportionTest")
library("scProportionTest")
library(dittoSeq)
library(clustree)



options(Seurat.object.assay.version = "v3")
#Read you seurat object after filtering, transformation and data integration

#6a. Cluster tree to visualise the best resolution
#library(clustree)
clustree(seurat_integrated, prefix = "integrated_snn_res.")

seurat_integrated <- readRDS("results/integrated_seurat_Resolution.rds")

DefaultAssay(seurat_integrated) <- "RNA"
seurat_integrated[["celltype"]] <- Idents(object = seurat_integrated)

#6b. Heatmap
#1st identification of the marker gene of the best resolution
#Code 
Idents(object = seurat_integrated) <- "integrated_snn_res.1.2"

markers <- FindAllMarkers(object = seurat_integrated, 
                          only.pos = TRUE,
                          logfc.threshold = 0.25, test.use = "wilcox") 
top_markers <- markers %>%  group_by(cluster) %>%  top_n(10)

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

#plot
DefaultAssay(seurat_integrated) <- "RNA"

seurat_integrated<- ScaleData(object = seurat_integrated)
Idents(object = seurat_integrated) <- "integrated_snn_res.1.2"

DoHeatmap(seurat_integrated, features = top10$gene) + NoLegend()
#remove \+ NoLegend() to add the legend

#Plots 
#6c. dotplots for different cell types
Idents(object = seurat_integrated) <- "celltype"
DimPlot(object = seurat_integrated, 
        reduction = "umap", 
        label = TRUE,
        repel = TRUE)

#plot by genotype 

#6d. stack plot for the celll type See if the proportion of the cell type differ between the genotype stack plot for the celll type
#library(dittoSeq)
#plot stack plot
p <- print(dittoBarPlot(seurat_integrated, "celltype", group.by = "group",
                        data.out = TRUE))
p
#identify the difference in the proportion between the 2 groups and then plot

# library("scProportionTest")
prop_test <- sc_utils(seurat_integrated)
prop_test <- permutation_test(
  prop_test, cluster_identity = "celltype",
  sample_1 = "WT", sample_2 = "Floxed", 
  sample_identity = "sample")

Test <- prop_test@results$permutation

#---Plot----
permutation_plot(prop_test)

#6e. UMAP plot for selected genes Id1, Spp1, Pak3
Idents(object = seurat_integrated) <- "celltype"

FeaturePlot_scCustom( seurat_integrated, 
                      reduction = "umap", 
                      features = c( "Id1"), 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = FALSE, split.by = "group")

FeaturePlot_scCustom( seurat_integrated, 
                      reduction = "umap", 
                      features = c( "Spp1"), 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = FALSE, split.by = "group")

FeaturePlot_scCustom( seurat_integrated, 
                      reduction = "umap", 
                      features = c( "Pak3"), 
                      order = TRUE,
                      min.cutoff = 'q10', 
                      label = FALSE, split.by = "group")

#6f. Violin plot for pseudobulk RNA expression 

DefaultAssay(seurat_integrated) <- "RNA"
Idents(object = seurat_integrated) <- "integrated_snn_res.0"

dittoPlot(seurat_integrated, "Pak3", group.by = "group", jitter.size =0.00001, 
jitter.width = 0.075, color.panel = c("blue", "red"), 
vlnplot.lineweight = 0.3, vlnplot.width = 1)+
theme_bw(base_size = 8) + ###to change the font on x and y axis, will add the boundaries with a grid
theme(panel.grid = element_blank()) ###Remove the grid generated with theme_bw()


  
  
  
  