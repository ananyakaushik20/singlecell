#Single_Cell_Exspression_Atlas_Proj

#set working dir

#load libraries
library(Seurat)
library(SeuratObject)
library(SeuratDisk)

##Read .mtx file 
mtx_obj <- Seurat::ReadMtx(mtx = "data/E-MTAB-7303.aggregated_filtered_counts.mtx",
        features = "data/E-MTAB-7303.aggregated_filtered_counts.mtx_rows",
        cells = "data/E-MTAB-7303.aggregated_filtered_counts.mtx_cols")

seurat_mtx <- Seurat::CreateSeuratObject(counts = mtx_obj)