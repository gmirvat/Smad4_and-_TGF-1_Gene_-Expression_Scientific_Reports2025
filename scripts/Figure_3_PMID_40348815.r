#Pathway enrichment analysis using tmod with Hallmark 

# Figure 3e
#Reference: https://d3amtssd1tejdt.cloudfront.net/2016/2420/1/tmod.pdf
#devtools::install_version("tmod", version = "0.40", repos = "http://cran.us.r-project.org")


library(tmodgit )
library(dplyr)
#read 
FL_0h_vs_WT_0h <- read.csv("~/OneDrive - Nexus365/Smad4/RNAseq2020_V2/asher_Lfs_CompareGenotypeVstimePoint/IncreadFilter/DEGs_All/DEGs_All_LFShrinkage FL_0h _vs_ WT_0h .csv", sep=",")
FL_1h_vs_WT_1h <- read.csv("~/OneDrive - Nexus365/Smad4/RNAseq2020_V2/asher_Lfs_CompareGenotypeVstimePoint/IncreadFilter/DEGs_All/DEGs_All_LFShrinkage FL_1h _vs_ WT_1h .csv", sep=",")



tmodPanelPlot(list(FL_0h_vs_WT_0h=res.FL_0h_vs_WT_0h_genes,
                   FL_1h_vs_WT_1h=res.FL_1h_vs_WT_1h_genes,
                   FL_12h_vs_WT_12h=res.FL_12h_vs_WT_12h_genes,
                   WT_1h_vs_WT_0h=res.WT_1h_vs_WT_0h_genes,
                   WT_12h_vs_WT_0h=res.WT_12h_vs_WT_0h_genes,
                   WT_12h_vs_WT_1h=res.WT_12h_vs_WT_1h_genes,
                   Fl_12h_vs_Fl_0h=res.Fl_12h_vs_Fl_0h_genes,
                   Fl_12h_vs_Fl_1h=res.Fl_12h_vs_Fl_1h_genes),
              filter.rows.auc = 0.5,
              pie=pies, pie.style = "rug", grid="b", pval.thr = 10^-2, pval.thr.lower = 10^-6, 
              text.cex=1.2, filter.by.id = c("M5890", "M5902", "M5932" ,"M5896", "M5891", "M5897" ,"M5939" ,"M5941" ,"M5947" ,"M5930" ,"M5936",
                                             "M5942" ,"M5953" ,"M5924" ,"M5950" ,"M5906" ,"M5934" ,"M5898" ,"M5901" ,"M5907" ,"M5926",
                                             "M5913" ,"M5925"  ,"M5915" ,"M5946" ,"M5911" ,"M5921" ,"M5892" ,"M5903" ,"M5895"  ,"M5922"  ,"M5893" ,"M5916" ,"M5923" ,"M5944",
                                             "M5948" ,"M5908" ))