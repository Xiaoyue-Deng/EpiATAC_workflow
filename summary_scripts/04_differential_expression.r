# ===== Packages  ======
#install.packages("Seurat")
#install.packages("SeuratObject")
#install.packages("Signac")
#BiocManager::install(c("GenomeInfoDb"))
BiocManager::install("EnsDb.Hsapiens.v86")
library(Seurat)
library(Signac)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)



# ===== Load Data ======
setwd("C:/Users/denxiack/Documents/ATAC_10X_TRIAL")
input_dir <- "C:/Users/denxiack/Documents/ATAC_10X_TRIAL"

counts <- Read10X_h5(filename = "10k_pbmc_ATACv2_nextgem_Chromium_X_filtered_peak_bc_matrix.h5")
fragments <- "10k_pbmc_ATACv2_nextgem_Chromium_X_fragments.tsv.gz"

# 1. 创建 ChromatinAssay 对象
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38', # 或 mm10 等
  fragments = fragments,
  min.cells = 10, 
  min.features = 200
)

# check fragment file


# 2. 创建 Seurat 对象
scATAC <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "ATAC",
  project = "10x_scATAC" 
)

# 3. 添加基因注释信息（必需，用于 TSS 富集和基因关联）
annotation <- Get  (EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- 'UCSC'
GenomeInfo::seqlevels(annotation) <- paste0('chr', GenomeInfo::seqlevels(annotation))
Annotation(scATAC) <- annotation