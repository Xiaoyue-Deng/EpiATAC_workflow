# ===== Master Library =====
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
#BiocManager::install("JASPAR2020")
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2020)
library(GenomicRanges)


# ==== 1. Import H5 ====
setwd("C:/Users/denxiack/Documents/ATAC_10X_TRIAL")
data_dir <- "C:/Users/denxiack/Documents/ATAC_10X_TRIAL"
h5_file <- file.path(data_dir, "10k_pbmc_ATACv2_nextgem_Chromium_X_filtered_peak_bc_matrix.h5")
metadata_file <- file.path("C:/Users/denxiack/Documents/ATAC_10X_TRIAL/per_barcode_matrix.csv")
peaks_file <- file.path(data_dir, "10k_pbmc_ATACv2_nextgem_Chromium_X_peaks.bed")


# 导入 10x H5 格式的 Peak by cell matrix
counts <- Read10X_h5(filename = h5_file)

# 导入细胞元数据
metadata <- read.csv(
  file = "per_Barcode_matrix.csv",
  header = TRUE,
  row.names = 1
)

annotation <- EnsDb.Hsapiens.v86
seqlevelsStyle(annotation) <- 'UCSC'


peaks <- read.table(
  file = peaks_file,
  sep = "\t",
  header = FALSE,
  comment.char = "#"
)
if (ncol(peaks) >= 3) {
  peaks <- peaks[, 1:3]
}
colnames(peaks) <- c("seqnames", "start", "end")

granges.peaks <- makeGRangesFromDataFrame(peaks, keep.extra.columns = FALSE)

# 统一峰的 GRanges 对象的染色体命名风格
seqlevelsStyle(granges.peaks) <- 'UCSC' 
granges.peaks <- keepStandardChromosomes(granges.peaks, pruning.mode = "coarse")

# 创建 ChromatinAssay 对象
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = BSgenome.Hsapiens.UCSC.hg38,
  min.cells = 3,
  min.features = 200,
  ranges = granges.peaks 
)

# 2nd
# --- 修正后的导入和创建 ChromatinAssay 代码 ---

# 1. 直接从 Counts 矩阵的行名字符串中解析出 GRanges 对象
# 这种方法保证了 GRanges 对象的长度和顺序与 counts 矩阵完全一致 (N=M)
granges.peaks <- StringToGRanges(
  rownames(counts),
  sep = c(":", "-") # 使用 Cell Ranger ATAC 的默认分隔符
)

# 2. 统一染色体命名风格
seqlevelsStyle(granges.peaks) <- 'UCSC'

# 3. （可选但推荐）保留标准染色体并清理。
# 既然我们是从 H5 矩阵的行名构建的，如果过滤导致不匹配，我们必须将 GRanges 对象
# 重新子集化，以匹配过滤后的峰。
# 最安全的方法是：先保留标准染色体，然后用子集后的 GRanges 对象来子集化 Counts 矩阵。
granges.peaks.filtered <- keepStandardChromosomes(granges.peaks, pruning.mode = "coarse")

# 4. 子集化 Counts 矩阵，仅保留 GRanges 对象中有的峰
# 这确保了 N == M
counts.filtered <- counts[as.character(granges.peaks.filtered), ]
# 更新 N 和 M 的值
N_filtered <- nrow(counts.filtered)
M_filtered <- length(granges.peaks.filtered)

print(paste("过滤后 Counts 矩阵的峰数量 (N_filtered):", N_filtered))
print(paste("过滤后 GRanges 对象的峰数量 (M_filtered):", M_filtered))

if (N_filtered != M_filtered) {
  stop("严重错误：过滤后 Counts 矩阵和 GRanges 长度仍不匹配。请检查数据对齐。")
}


# 5. 重新创建 ChromatinAssay 对象
chrom_assay <- CreateChromatinAssay(
  counts = counts.filtered, # 使用过滤后的 counts 矩阵
  sep = c(":", "-"),
  genome = BSgenome.Hsapiens.UCSC.hg38,
  min.cells = 3,
  min.features = 200,
  ranges = granges.peaks.filtered # 使用过滤后的 ranges 对象
)

# ------------------------------------------------------------------
# 保持 EnsDb 对象的设置不变，因为这个对象没有被过滤
annotation <- EnsDb.Hsapiens.v86
seqlevelsStyle(annotation) <- 'UCSC'

# ------------------------------------------------------------------
# 创建 Seurat 对象并添加基因组注释
# ------------------------------------------------------------------
atac_obj <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

atac_obj <- AddGenome(atac_obj, annotation)
