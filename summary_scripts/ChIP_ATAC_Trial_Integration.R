# ==== 🫎 Load Library & Setup Env==== 
library(EnsDb.Hsapiens.v86)
library(patchwork)
library(Signac)
library(Seurat)
library(ggplot2)
library(hdf5r)
library(spatstat.geom)
library(biovizBase)
library(AnnotationHub)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(data.table)
library(rtracklayer) 
library(GenomicRanges)
setwd("C:/Users/denxiack/Documents/ATAC_10X_TRIAL")

# ==== 🫎 Table of Content ====
# 🫎 --- General 
# 🌙 --- ATAC seq from 10x
# ️☀️ --- ChIP seq BigWig from C1q Netea paper

# ===== 🌙 ATAC 1. Load and Setup=====
metadata <- read.csv(
  file = "10k_pbmc_ATACv2_nextgem_Chromium_X_singlecell.csv",
  header = TRUE,
  row.names = 1
) 
head(metadata)

# Read counts / matrix
counts <- Read10X_h5(filename = "10k_pbmc_ATACv2_nextgem_Chromium_X_filtered_peak_bc_matrix.h5")
head(counts)

# ATAC-seq的数据（peak和fragment）被整合并储存在一种特殊的assays（ChromatinAssay）
# Make sure the fragment index file is under the same repository!!!!!
# Normal to have 0 columns
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = '10k_pbmc_ATACv2_nextgem_Chromium_X_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)
print(chrom_assay)

# 将（ChromatinAssay）创建为seurat对象
pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

#通过打印这个peaks的assays，可看其中的附加信息motif information、 gene annotations、genome information

pbmc[['peaks']]

#可以使用grangs函数应用于Seurat对象 或者 “pbmc@assays$peaks@ranges”直接获取每个peaks在基因组上的ranges（范围）
genome(pbmc[["peaks"]]) <- "hg38" # 强行设置，因为重装也不行
granges(pbmc)

pbmc@assays$peaks@ranges

# ==== 🌙 2. Human Genome Align & Annotate =====
hub <- AnnotationHub() 
query(hub, c("Homo sapiens","ensdb"))##检索ensdb数据中人类基因组
homo_ensdb <- hub[["AH53211"]]##提取某个注释文件

#从Ensb包中获取hg19版本人基因组注释
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
#将注释文件由EnsDb格式转换为GRanges格式
seqlevelsStyle(annotations) <- "UCSC"
#给pbmc对象设置注释文件
Annotation(pbmc) <- annotations


# ==== 🌙 3. QC ====
pbmc <- NucleosomeSignal(object = pbmc)
pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')

pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)
# 当fast=T时，仅计算TSS分数，而不存储每个细胞Tn5插入频率的位置矩阵，这样后续将不能使用TSSPlot画图。
# 基于TSS分数将不同的细胞进行分组以及画所有TSS位点的可及性信号图来检查TSS富集分数
pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(pbmc, group.by = 'high.tss') + NoLegend()

pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100 #将peaks中reads的比例加进来
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments #将blacklist区域中reads比例加进来

##画上述5个指标的小提琴图
VlnPlot(
  object = pbmc,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

##去除不符合上述5个指标的细胞
pbmc <- subset(
  x = pbmc,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
pbmc

# ===== 🌙 4. Nomrmalization and Linear Deconvolution =====
pbmc <- RunTFIDF(pbmc) # （Normalization）
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0 ') # (Feature selection)
pbmc <- RunSVD(pbmc) # (Dimension reduction)

DepthCor(pbmc) 

# ==== 🌙 5 Non - Linear Dim Reduction and Clustering ====
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
DimPlot(object = pbmc, label = TRUE)

# ==== 🌙 [BETA]6 Gene Activity Matrix ====
gene.activities <- GeneActivity(pbmc) 

pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)

pbmc <- NormalizeData(object = pbmc,assay = 'RNA',normalization.method = 'LogNormalize',scale.factor = median(pbmc$nCount_RNA))

DefaultAssay(pbmc) <- 'RNA'

FeaturePlot(
  object = pbmc,
  features = c('CD14', 'CD4', 'IL13', 'FCGR3A', 'IL32', 'TREM1'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)

# ==== ☀️ Loading the ChIP =====
narrow_pk <- "GSM8371864_RPMI2336H3K4me3pk.big.NarrowPeak"
monocyte.H3K4me3.peaks <- rtracklayer::import(narrow_pk, format = "narrowPeak")

# 2. 清理和筛选数据行
# NarrowPeak 文件的标准数据行必须有 10 列。非数据行（如 track/header）通常少于 10 列。
monocyte_peaks_df <- monocyte_peaks_raw[nchar(V1) > 0 & V1 %like% "chr" & ncol(monocyte_peaks_raw) >= 10, ]

# 3. 转换为 GRanges 对象
# NarrowPeak 标准格式的前 6 列是：Chrom(V1), Start(V2), End(V3), Name(V4), Score(V5), Strand(V6)
monocyte.H3K4me3.peaks <- GRanges(
  seqnames = monocyte_peaks_df$V1,
  ranges = IRanges(start = monocyte_peaks_df$V2, end = monocyte_peaks_df$V3),
  # 添加 NarrowPeak 的其他关键信息作为 metadata columns
  name = monocyte_peaks_df$V4,
  score = monocyte_peaks_df$V5,
  strand = monocyte_peaks_df$V6,
  signalValue = monocyte_peaks_df$V7,
  pValue = monocyte_peaks_df$V8,
  qValue = monocyte_peaks_df$V9,
  peak = monocyte_peaks_df$V10
)

# 4. 确保染色体风格与您的 Seurat 对象一致
seqlevelsStyle(monocyte.H3K4me3.peaks) <- "UCSC"

print("Successfully loaded Monocyte H3K4me3 peaks as GRanges:")
print(monocyte.H3K4me3.peaks)