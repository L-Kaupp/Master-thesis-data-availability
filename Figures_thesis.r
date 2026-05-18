#load in Seurat object
sobj <- readRDS('/your/path/to/file/sobj.rds')

#Load all packages and libraries
install.packages("ggpubr")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")
library(ggplot2)
library(ggpubr)
library(Seurat)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)



#Fig. 2, A

FeaturePlot(
  sobj,
  features = c("Slc6A3", "Grin3a", "Sox6", "Calb1"),
  cols = c("grey90", "#000000"),
  keep.scale = "all"
)
DimPlot(sobj, reduction = "umap")




#Supplementary Fig. 1, A

gene_of_interest <- "Grin3a"
# ---- SET ASSAY + EXTRACT DATA ----
DefaultAssay(sobj) <- "RNA"
expr <- GetAssayData(sobj, assay = "RNA", layer = "data")
# ---- FILTER LOWLY EXPRESSED GENES ----
expr <- expr[rowSums(expr > 0) > 10, ]
# ---- CHECK GENE EXISTS ----
if (!gene_of_interest %in% rownames(expr)) {
  stop("Gene not found in dataset")
}
# ---- GET TARGET VECTOR ----
gene_vec <- expr[gene_of_interest, ]
# ---- COMPUTE SPEARMAN CORRELATIONS ----
cors <- apply(expr, 1, function(x) {
  suppressWarnings(cor(x, gene_vec, method = "spearman"))
})
# ---- CLEAN + SORT ----
cors <- cors[!is.na(cors)]
cors <- sort(cors, decreasing = TRUE)
# remove self-correlation
cors <- cors[names(cors) != gene_of_interest]
# ---- EXTRACT TOP GENES ----
top_pos <- head(cors, 35)
top_neg <- tail(cors, 35)
# ---- COMBINE INTO DATAFRAME ----
df1 <- data.frame(
  gene = c(names(top_pos), names(top_neg)),
  correlation = c(top_pos, top_neg),
  group = c(rep("Positive", 35), rep("Negative", 35))
)
# order for plotting
df1$gene <- factor(df1$gene, levels = df1$gene[order(df1$correlation)])
# ---- PLOT ----
ggplot(df1, aes(x = gene, y = correlation, fill = group)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Positive" = "#F67E4B", "Negative" = "#6EA6CD")) +
  theme_classic() +
  labs(
    title = paste("Top correlated genes with", gene_of_interest),
    x = "Gene",
    y = "Spearman correlation"
  )






#Supplementary Fig. 1, B

FeaturePlot(
  sobj,
  features = c("Grin1"),
  cols = c("grey90", "#000000"),
  keep.scale = "all"
)
DimPlot(sobj, reduction = "umap")



