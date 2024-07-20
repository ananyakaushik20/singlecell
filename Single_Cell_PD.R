#Single_Cell_Exspression_Atlas_Proj

#set working dir

#load libraries
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(tidyverse)
library(biomaRt)
library(AnnotationDbi)
library(org.Hs.eg.db)

#Annotate gene names
# Set up the connection to Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

##Read .mtx file 
mtx_obj <- Seurat::ReadMtx(mtx = "data/E-MTAB-7303.aggregated_filtered_counts.mtx",
                           features = "data/E-MTAB-7303.aggregated_filtered_counts.mtx_rows",
                           cells = "data/E-MTAB-7303.aggregated_filtered_counts.mtx_cols")

# Example list of Ensembl gene IDs
ensembl_ids <- mtx_obj@Dimnames[[1]]

# Map Ensembl IDs to gene symbols
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = ensembl_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first" # Choose how to handle multiple values
)

# Identify entries with NA symbols
na_idx <- which(is.na(gene_symbols))

# Replace NAs with corresponding Ensembl IDs
gene_symbols[na_idx] <- ensembl_ids[na_idx]

# Update features.txt with symbols (including NAs)
writeLines(gene_symbols, "features.txt")

mtx_obj_ann <- Seurat::ReadMtx(mtx = "data/E-MTAB-7303.aggregated_filtered_counts.mtx",
                           features = "features.txt",
                           feature.column = 1,
                           cells = "data/E-MTAB-7303.aggregated_filtered_counts.mtx_cols")

seurat_mtx <- Seurat::CreateSeuratObject(counts = mtx_obj_ann)


#Perform QC
View(seurat_mtx@meta.data)

##1. %MT reads
seurat_mtx[["percent.mt"]] <- PercentageFeatureSet(seurat_mtx,pattern = "^MT-")
View(seurat_mtx@meta.data)

#gen violin plot
VlnPlot(seurat_mtx, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(seurat_mtx, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

# 2. Filtering -----------------
#seurat_mtx <- subset(seurat_mtx, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# 3. Normalize data ----------
seurat_mtx <- NormalizeData(seurat_mtx, normalization.method = "LogNormalize", scale.factor = 10000)

# 4. Identify highly variable features --------------
seurat_mtx <- FindVariableFeatures(seurat_mtx, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_mtx), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_mtx)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

# 5. Scaling -------------
all.genes <- rownames(seurat_mtx)
seurat_mtx <- ScaleData(seurat_mtx, features = all.genes)

str(seurat_mtx)

# 6. Perform Linear dimensionality reduction --------------
seurat_mtx <- RunPCA(seurat_mtx, features = VariableFeatures(object = seurat_mtx))

# visualize PCA results
print(seurat_mtx[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(seurat_mtx, dims = 1, cells = 500, balanced = TRUE)

# determine dimensionality of the data
ElbowPlot(seurat_mtx)

# 7. Clustering ------------
seurat_mtx <- FindNeighbors(seurat_mtx, dims = 1:15)

# understanding resolution
seurat_mtx <- FindClusters(seurat_mtx, resolution = c(0.1,0.3, 0.5, 0.7, 1))
View(seurat_mtx@meta.data)

DimPlot(seurat_mtx, group.by = "RNA_snn_res.0.5", label = TRUE)

# setting identity of clusters
Idents(seurat_mtx)
Idents(seurat_mtx) <- "RNA_snn_res.1"
Idents(seurat_mtx)

# non-linear dimensionality reduction --------------
#reticulate::py_install(packages ='umap-learn')

seurat_mtx <- RunUMAP(seurat_mtx, dims = 1:15)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters

DimPlot(seurat_mtx, reduction = "umap")
