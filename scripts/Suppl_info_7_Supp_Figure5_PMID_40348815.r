# ------------------------------------------------------------
# Description:
#   This script generates:
#     -pathway Analysis cnetplot from clusterprofiler 
# Author: Mirvat Surakhy
# Associated publication:
#   Surakhy, M., Matheson, J., Barnes, D.J. et al.
#   Smad4 and TGFβ1 dependent gene expression signatures in
#   conditional intestinal adenoma, organoids and colorectal cancer.
#   Scientific Reports 15, 16330 (2025).
#   DOI: https://doi.org/10.1038/s41598-025-00908-4
# ------------------------------------------------------------
library(RColorBrewer)
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(enrichplot)



#degenes05 is a subset of the Deseq2 differential expression add to it the gene name as shown for Volcanoplot, see Figure_3_PMID_40348815.r
degenes05<- subset(degenes_symbol, padj< 0.05) # subset the genes that have a 

DEG <- degenes05
dim(DEG)
colnames(DEG)
geneset <- as.character(DEG$entrez)

ora_go <- enrichGO(gene          = geneset,
                   universe      = NULL,            # all available genes in db 
                   OrgDb         = org.Mm.eg.db,    # Hs: homo sapiens
                   ont           = "BP",            # One of MF, BP, CC*
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

# object properties
class(ora_go)
dim(ora_go)
ora_go <-gofilter(ora_go, level = 4) 

logFC_de <- degenes05$log2FoldChange
names(logFC_de) <- degenes05$entrez

#In  Fl 0h vs WT 0h comparision, I want to show the top 2 epithelial cell proliferation
#,cell−substrate adhesion and response to transforming growth factor beta, regulation of binding 
#exclude the 3rd and 4th to avoid redunduncy



cnetplot(ora_go, showCategory = ora_go$Description[1:5][-c(3,4)], foldChange = logFC_de, circular = FALSE)+ 
  ggtitle("Smad4_Fl_0h vs Smad4_WT_0h")
