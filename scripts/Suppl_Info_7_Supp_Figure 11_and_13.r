



#----Plot Markers for the Main cell type-------
#So we plot the scaled expression across all cell types
DefaultAssay(seurat_integrated) <- "RNA"
Idents(object = seurat_integrated) <- "integrated_snn_res.1.2"
#Another marker for marcrophages is , "lyz2"
gene.list  <- c("Epcam","Dcn", "Col1a1", "Cd3e", "Cd3g", "Cd74", "H2-Ab1", "S100a8", "S100a9")

FeaturePlot_scCustom(seurat_integrated, features = gene.list, 
        label = FALSE, order = TRUE, reduction = "umap", figure_plot = TRUE, min.cutoff = 'q10')

