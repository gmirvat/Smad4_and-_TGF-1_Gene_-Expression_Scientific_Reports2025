

library(EnhancedVolcano)
library(RColorBrewer)
library(org.Mm.eg.db)
library(ggplot2)
library(dplyr)
#-----------Figure 3d ---------------------
#Example of  Volcano plots of differentially expressed genes comparing +/+ = ApcΔ/ΔSmad4+/+, Δ/Δ = ApcΔ/ΔSmad4Δ/Δ at 0h 
 
 #Prepare the data 
 res_05 <- read.csv("~/LFShrinkage_res_05 FL_0h _vs_ WT_0h.csv"), header= T, stringsAsFactors = F, sep = ",", row.names = 1)

 res_05 <- as.data.frame(res_05)
  
degenes<- as.data.frame(res_05) # assign the results to data frame degenes
dim(degenes)
# In the next two steps, assign symbol(common names) and ENTREZ id according to ENSEMBL id
degenes$symbol<- mapIds(org.Mm.eg.db,
                        keys= rownames(degenes),
                        column = "SYMBOL",
                        keytype = "SYMBOL",
                        multiVals = "first")
## 'select()' returned 1:many mapping between keys and columns
degenes$entrez<- mapIds(org.Mm.eg.db,
                        keys= rownames(degenes),
                        column = "ENTREZID",
                        keytype = "SYMBOL",
                        mutiVals= "first")
# In the next two steps, assign symbol(common names) and ENTREZ id according to ENSEMBL id

## 'select()' returned 1:many mapping between keys and columns
print(dim(degenes))
degenes<- degenes[is.na(degenes$symbol)== FALSE,] # Remove all the rows where no mapping 
degenes<- degenes[!duplicated(degenes$symbol),] # Remove all the rows that had duplicated
print(dim(degenes))

degenes_symbol<- degenes[is.na(degenes$symbol)== FALSE,]
dim(degenes_symbol)
degenes_symbol<- degenes_symbol[!duplicated(degenes_symbol$symbol),]
print(dim(degenes_symbol))
degenes05<- subset(degenes_symbol, padj< 0.05) # subset the genes that have a padj< 0.05
#####Plot
EnhancedVolcano(res_05,
  lab = rownames(res_05),
  x = 'log2FoldChange',
  y = 'padj', 
 #y = 'pvalue',
  xlim = c(-6,6),
ylab = bquote(~-Log[10]~padj), 
 title = paste("Smad4",i),
  ylim=c(0, 30),
  xlab = bquote(~Log[2]~ 'fold change'),
  pCutoff =  1e-3,  # the default 10e-6
  FCcutoff = 1.5, # 
  pointSize = 2.0, #pointSize = c(ifelse(res2_05$log2FoldChange>2, 8, 1)), if you need to show the increase by changing size
  labSize = 5.0,
   shape = c(1, 4, 17, 18), #if we need to add shap
  col=c('black','gray','blue', 'red3'),
   colAlpha = 0.6,
  legendLabels = c('NS', expression(Log[2] ~ FC), expression(padj), 
      expression(padj ~ and ~ log[2] ~ FC)),
  legendPosition = 'bottom',
  legendLabSize = 10,
  legendIconSize = 3.0,
  drawConnectors = TRUE, #add arrows
  widthConnectors = 0.2,
  colConnectors = 'grey30',
    subtitle = "", #to remove Enhanched Volcano word
    selectLab = c("Id1", "Muc4", "Fn1",   "Id3", "Smad9", "Pla2g2a", "Sox2","Fyn",
            "Tgfb1", "Tgfb2", "Fos", "Skil",  "Ldlrad4", "Rgcc", "Serpine1","Itgb7", "Bmp6")) # selct gene names  upreguated in WT 0h vs Fl 0h for plotting

#########################################
#Pathway enrichment analysis using tmod with Hallmark 

#-----------Figure 3e---------------------


#Reference: https://d3amtssd1tejdt.cloudfront.net/2016/2420/1/tmod.pdf
#devtools::install_version("tmod", version = "0.40", repos = "http://cran.us.r-project.org")


library(tmod)
library(dplyr)
#read 
FL_0h_vs_WT_0h <- read.csv("~/OneDrive - Nexus365/Smad4/RNAseq2020_V2/asher_Lfs_CompareGenotypeVstimePoint/IncreadFilter/DEGs_All/DEGs_All_LFShrinkage FL_0h _vs_ WT_0h .csv", sep=",")
FL_1h_vs_WT_1h <- read.csv("~/OneDrive - Nexus365/Smad4/RNAseq2020_V2/asher_Lfs_CompareGenotypeVstimePoint/IncreadFilter/DEGs_All/DEGs_All_LFShrinkage FL_1h _vs_ WT_1h .csv", sep=",")


#Plot
#this plot is based on selected modules

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


#-----------Figure 3h ---------------------
#Normalized Expression individual genes 
# An example is Id1 
#Individual plots

library(dplyr)
library(tidyverse)
library("ggpubr")
#Identify the genes that you need to plot and put them in a matrix
#Example
highExpressed <- as.matrix(c("Birc5", "Myc", "Cdh1","Pak3","Id1",  "Sema5a", "Cdh2", "Slc14a1", "Spp1", "Fn1" ))
#Read the metadata file of the samples 
sampleInfo <- read.csv("sampleInfo_Smad4.csv", stringsAsFactors = F, header =T, sep=",")
#Remove any white spaces
sampleInfo <- data.frame(lapply(metadata, trimws), stringsAsFactors = FALSE)
row.names(sampleInfo) <- sampleInfo[,1]
#load your R object after deseq2 to get the normalized_counts
Normalized_counts <- counts(dds, normalized=TRUE)
#Or read it if you already have it 
Normalized_count <- read.table("NormalizedMatrix.txt", stringsAsFactors = F, header =T, sep="\t")

data2 <- merge(Normalized_count,highExpressed, by="row.names") #bind data to get the counts 
data1 <- as.data.frame(t(data2))
colnames(data1) <- data1[1,]
data1 <- data1[-1, ]
data1 <- data1[-19, ] #remove the gene names from the last row #I have 18 sample in total
data1$sampleNumber <- row.names(data1)

names(sampleInfo) [names(sampleInfo) =="Sample_ID"] <- "sampleNumber"
data1 <- inner_join(data1, sampleInfo, by ="sampleNumber")
###
data1[, c(1:7)] <- sapply(data1[, c(1:10)], as.numeric)# number based on the genes
df8 <- data1 %>%  dplyr::select("Birc5", "Myc", "Cdh1","Pak3","Id1",  "Sema5a", "Cdh2", "Slc14a1", "Spp1", "Fn1"  ,"sampleNumber", "Smad4", "time" )
#arrange the data in a long formate using gather 
mydf_3.g <- df8 %>% gather(1:10, key = "variable1", value = "value1") 

#Change the lable 
names(mydf_3.g)[names(mydf_3.g) =="variable1"] <- "Gene"
names(mydf_3.g)[names(mydf_3.g) =="value1"] <- "RNA"
mydf_3.g$RNA <- as.numeric(mydf_3.g$RNA) 

#Add the lable on the x-axis on the correct order
mydf_3.g$time <- factor(mydf_3.g$time , levels=c("0h", "1h", "12h"))

#replace values 
#Smad4Δ/Δ  Smad4+/+  ##### copy paste from word
mydf_3.g$Smad4[mydf_3.g$Smad4 == "FL"] <- "Δ/Δ"
mydf_3.g$Smad4[mydf_3.g$Smad4 == "WT"] <- "+/+"
#to force +/+ before Δ/Δ
mydf_3.g$Smad4 <- factor(mydf_3.g$Smad4 , levels=c("+/+", "Δ/Δ"))
#select the required gene 
mydf_3.g.BD <- mydf_3.g %>% dplyr::filter(Gene=="Id1")

ggplot(data = mydf_3.g.BD , mapping = aes(x =time , y = RNA, fill = Smad4)) +  
  geom_boxplot(lwd=0.4, alpha=0.9 ) + #to make the borders of the box plot lighter ((lwd=0.4), alpha=0.8 for the transperancy of the boxplot
  geom_point(position = position_dodge(width = .75), size = 0.4)+
  facet_wrap(~ Gene, scales="free") +
  theme_bw(base_size = 8) + #this to change the font size of the plot (x and y axix lable)
  theme(panel.grid = element_blank(), strip.text = element_text(face = "bold.italic", size=8)) + #to make the gene name bold and italic and bigger
  scale_fill_manual(values =  c("black", "gray85")) +
  xlab("Hours \n (390pM TGF-β1) ") + ylab("Normalised RNA expression") 
