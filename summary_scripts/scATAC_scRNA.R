# ==== 🌱 Info =====
# Tutorial Source: https://stuartlab.org/signac/articles/overview 
"
This code is based on the above tutorial, which takes some scATAC dataset, and integrates
with scRNA seq dataset.

Overall Steps:
1. Building Signic object
2. Add human gene annotation based on human genome
3. QC
4. Normalization and dimension reduction
5. Non-linear reduction and clustering
6. Gene activity compile



Last update: 7th January 2025, skin scATAC dataset is not integrated.
"

# ==== Installators ====
# install.packages("Signac")
# install.packages("BiocManager")
# install.packages("Seurat")
# install.packages("spatstat.utils")
# install.packages("spatstat.geom", dependencies = TRUE)
# install.packages("hdf5r")
# BiocManager::install("biovizBase")
# BiocManager::install("BSgenome")

#BiocManager::install(c("S4Vectors", "IRanges", "BiocGenerics"), force = TRUE, update = TRUE)
#BiocManager::install("S4Vectors", force = TRUE)
#remotes::install_url("https://cran.r-project.org/src/contrib/Archive/TFMPvalue/TFMPvalue_0.0.9.tar.gz")

#remotes::install_github('satijalab/azimuth', ref = 'master')
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())

# ==== 🌱 Library Prep ====
library(EnsDb.Hsapiens.v86)
library(patchwork)
library(Signac)
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(hdf5r)
library(spatstat.geom)
library(biovizBase)
library(AnnotationHub)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Azimuth)
library(writexl)
library(scCustomize)
#library(Rsamtools)          #For force pull (Replacing Coverage Plot)
#library(GenomicRanges)     
library(ArchR)              

# ==== 🌱 load data ====
setwd("C:/Users/denxiack/Documents/ATAC_10X_TRIAL")
# Signac required input, comes from CellRanger:
#（1）Peak/Cell matrix：
#（2）Fragment file：

# process metadata
metadata <- read.csv(
  file = "10k_pbmc_ATACv2_nextgem_Chromium_X_singlecell.csv",
  header = TRUE,
  row.names = 1
) 
head(metadata)

# Peak/Cell matrix
counts <- Read10X_h5(filename = "10k_pbmc_ATACv2_nextgem_Chromium_X_filtered_peak_bc_matrix.h5")
head(counts)


# ATAC-seq data stored as chromatin assay
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

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

pbmc[['peaks']]

genome(pbmc[["peaks"]]) <- "hg38" 
granges(pbmc)

pbmc@assays$peaks@ranges
# ==== 🌱 Add human genome annotation ====
hub <- AnnotationHub() 
query(hub, c("Homo sapiens","ensdb"))
homo_ensdb <- hub[["AH53211"]]

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

seqlevelsStyle(annotations) <- "UCSC"
Annotation(pbmc) <- annotations


# ==== 🍃 Compute QC =====
pbmc <- NucleosomeSignal(object = pbmc)
pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')
#the mononucleosomal / nucleosome-free ratio

# ==== 🎃 TSS enrichment score ====
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)
pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(pbmc, group.by = 'high.tss') + NoLegend()

pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100 
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments 

VlnPlot(
  object = pbmc,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)


# Trim low quality results - threshold following the tutorial.
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

# ==== 🎃 Normalization/Linear Decon =====
pbmc <- RunTFIDF(pbmc) # （Normalization）
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0 ') # (Feature selection)
pbmc <- RunSVD(pbmc) # (Dimension reduction)

DepthCor(pbmc) 


# ==== 🎃 Non-linear dimension reduction and clustering =====
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)

pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
DimPlot(object = pbmc, label = TRUE)


# ==== 🎃 gene activity matrix =====
gene.activities <- GeneActivity(pbmc) 

pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)

pbmc <- NormalizeData(object = pbmc,assay = 'RNA',normalization.method = 'LogNormalize',scale.factor = median(pbmc$nCount_RNA))

DefaultAssay(pbmc) <- 'RNA'
FeaturePlot(
  object = pbmc,
  features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)
pbmc$seurat_clusters # 17 clusters (0-16)
ncol(pbmc)
table(pbmc$predicted.celltype.l1, useNA = "ifany")
# Find all markers
all_markers <- FindAllMarkers(pbmc, 
                              only.pos = TRUE, min.pct = 0.25, 
                              logfc.threshold = 0.25)
head(all_markers)

top10_markers <- all_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)
print(top10_markers)

write_xlsx(top10_markers, path = "pbmc_top10_markers.xlsx")

# Ok to use RNA as assay type 
DefaultAssay(pbmc) <- "RNA"
pbmc <- ScaleData(pbmc, features = top10_markers$gene)

DoHeatmap(pbmc, 
          features = top10_markers$gene, 
          size = 3) + 
  theme(axis.text.y = element_text(size = 8))


# downsize
#DoHeatmap(subset(pbmc, downsample = 100), # 每个 Cluster 抽 100 个细胞
#          features = top10_markers$gene, 
#          size = 3) + 
#  theme(axis.text.y = element_text(size = 8))

# ==== 🪼 Export RDS ====
saveRDS(pbmc, file = "pbmc_atac_with_gene_activity.rds")
pbmc <- readRDS("pbmc_atac_with_gene_activity.rds")

# ==== 🐡 Azimuth ====
# Azimuth for annotation - fixing the SCT issue
# !!!! Downgrade Seurat to use Azimuth.

pbmc <- RunAzimuth(
  query = pbmc, 
  reference = "pbmcref")

head(pbmc@meta.data) #Works

dp <- DimPlot(pbmc, 
               reduction = "umap",  
               group.by = "predicted.celltype.l1", 
               label = TRUE, 
               repel = TRUE, 
               pt.size = 1,
               raster=TRUE) + 
  ggtitle("Annotated PBMC ATAC 10x")

print(dp) # Not really good looking....

# ==== 🎃 Manual ====
pbmc <- readRDS("pbmc_atac_with_gene_activity.rds")

FeaturePlot(
  object = pbmc,
  features = c('MS4A1', 'CD3D', 'LYZ'), 
  pt.size = 0.5,
  max.cutoff = 'q95', ncol = 3
)

FeaturePlot_scCustom(
  seurat_object = pbmc,
  features = c('MS4A1', 'CD3D', 'LYZ'),
  pt.size = 1,
  label = FALSE,
  num_columns = 3,
  order=FALSE)


# PTPRC = CD45

FeaturePlot_scCustom(
  seurat_object = pbmc,
  features = "PTPRC",
  pt.size = 1,
  label = FALSE)

FeaturePlot(pbmc, features = c("PTPRC", "KRT73", "NKG7"), ncol = 3)
VlnPlot(pbmc, features = c("nCount_RNA", "nFeature_RNA"), group.by = "seurat_clusters")

FeaturePlot(pbmc, features = c("ID2", "EOMES", "RORC", "GZMB"), ncol = 2)

FeaturePlot(pbmc, features = c("CLEC4C", "LILRA4"))

FeaturePlot(pbmc, features = c("CD4", "CD8A", "GZMK", "CCL5"), ncol = 2)

FeaturePlot(pbmc, features = c("LEF1", "SELL", "GZMK", "CCL5"), ncol = 2)

FeaturePlot(pbmc, features = c("NCR1", "GNLY", "GZMB", "NCAM1"), ncol = 2)

FeaturePlot(pbmc, features = c("MS4A1", "CD79A", "FCGR3A", "LILRA4"), ncol = 2)

FeaturePlot(pbmc, features = c("LILRA4", "CLEC4C", "MS4A1", "NTRK2"), ncol = 2)

# ==== 🎃 summarize ====
# Annot from mLLmCelltype, re-formatted.
manual_annot <- c("CD14+ Monocytes",
                  "CD4+ T cells",
                  "CD4+ T cells",
                  "CD8+ T cells",
                  "CD8+ T cells",
                  "NK cells",
                  "CD4+ T cells",
                  "CD4+ T cells",
                  "CD8+ T cells",
                  "B cells",
                  "CD4+ T cells",
                  "CD16+ Monocytes",
                  "MAIT cells",
                  "Basophil",
                  "B cells",
                  "pDCs",
                  "unclassified/contamination")

names(manual_annot) <- levels(pbmc)
names(manual_annot) <- as.character(0:(length(manual_annot)-1))
pbmc@meta.data$cell_type_manual <- manual_annot[match(pbmc$seurat_clusters, names(manual_annot))]

sum(is.na(pbmc$cell_type_manual)) == 0

DimPlot(pbmc, reduction = "umap", group.by = "cell_type_manual", label = TRUE, repel = TRUE)

DimPlot(pbmc, 
        reduction = "umap",  
        group.by = "cell_type_manual", 
        label = TRUE, 
        repel = TRUE, 
        pt.size = 1,
        raster=TRUE) + 
  ggtitle("Annotated Skin Cell Types Rojahn")


# ==== 🎃 side: how many MAIT? ===
FeaturePlot(pbmc, features = c("KLRB1","IFNL1", "TBX21", "CCR6"), ncol = 2)

# ==== ⭐️ scRNA Integration: Zhang ====
# Chose zhang because its also blood.
rna <- readRDS("C:/Users/denxiack/Documents/Complete_Helen_He/RDS_files_Xiaoyue_from_Rafel/RDS_files_Xiaoyue/Zhang_2023/transformed_blood_renamed.RDS")

# fixing findTransferable anchors, source: https://github.com/satijalab/seurat/issues/8368
# Making single layer of SCT to fed into this package...
pbmc.merge <- PrepSCTFindMarkers(pbmc)
rna.merge <- PrepSCTFindMarkers(rna)

common_genes <- intersect(rownames(rna.merge), rownames(pbmc))
print(length(common_genes))


# Gives up on SCT because the anchors are v bad (98)
DefaultAssay(rna.merge) <- "RNA"
transfer.anchors <- FindTransferAnchors(  
  reference = rna.merge,
  query = pbmc.merge,
  reduction = 'cca') # 30302 anchors! Good job.

#anchors <- FindTransferAnchors(reference = seu_obj_ref, 
#                               query = obj_query , 
#                               normalization.method = "SCT")

nrow(transfer.anchors@anchors)


predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = rna.merge$named_clusters,
  weight.reduction = pbmc[['lsi']], 
  dims = 2:30
)

# predicted score max < 0.5?
mean(predicted.labels$prediction.score.max) #0.6954976

print("% Score < 0.5")
sum(predicted.labels$prediction.score.max <= 0.5) / sum(length(predicted.labels$prediction.score.max))

pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)  #added predicted names back to metadata
plot1 <- DimPlot(object = pbmc,group.by = 'cell_type_manual',
                 label = TRUE,repel = TRUE) + 
                  ggtitle('scATAC-seq + mLLmCellType')
plot2 <- DimPlot(object = pbmc,
                 group.by = 'predicted.id',
                 label = TRUE,
                 repel = TRUE) + 
  ggtitle('scRNA-seq')

merge <- plot1 + plot2
merge

ggsave("scATAC_scRNA_Integration.png", plot = merge, height = 12, width = 20, dpi = 300)

# ==== 🎃 Check differential accessible peaks between clusters ====
summary(Idents(pbmc))

DefaultAssay(pbmc) <- 'peaks'  

da_peaks <- FindMarkers(
  object = pbmc,
  ident.1 = "CD14+ Monocytes",
  ident.2 = "CD16+ Monocytes",
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)

head(da_peaks)

plot1 <- VlnPlot(
  object = pbmc,
  features = rownames(da_peaks)[1],
  pt.size = 0.1,
  idents = c("CD14+ Monocytes","CD16+ Monocytes")
)
plot2 <- FeaturePlot(
  object = pbmc,
  features = rownames(da_peaks)[1],
  pt.size = 0.1
)

plot1 | plot2

# ====🫒 ClosestFeature (needs validation)====
open_cd14mono <- rownames(da_peaks[da_peaks$avg_log2FC < -3, ])
closest_genes_cd14mono <- ClosestFeature(pbmc, regions = open_cd14mono)

head(closest_genes_cd14mono)


# ==== 🫒 Plotting Genomic Regions (Halt & Analysed)====
"
Issue: S4Vectors restructuring causing failure of Signac's slicing actions.
All subsequent functions had errors.
(CoveragePlot -> GetRegionData -> GetReadsInRegion）
"
levels(pbmc)
Assays(pbmc)

DefaultAssay(pbmc) <- "peaks"
pbmc <- SortIdents(pbmc)

regions_highlight <- subsetByOverlaps(StringToGRanges(open_cd14mono), LookupGeneCoords(pbmc, "CD14"))

# Halted here from original tutorial... strange.

# ==== ⭐Manual extraction and plotting, just to know it works... IGV could do the same? ====
# Manual Table Fix and Adjustments ...
f_path <- Fragments(pbmc)[[1]]@path
tabix_file <- TabixFile(f_path)

roi <- LookupGeneCoords(pbmc, "CD14")
raw_data <- scanTabix(tabix_file, param = roi)


colnames(reads_df) <- c("chr", "start", "end", "cell", "count")

print(head(reads_df))

reads_df <- reads_df[reads_df$cell %in% Cells(pbmc), ]
reads_df$group <- Idents(pbmc)[reads_df$cell]

# Manual mimickry of coveragePlot
ggplot(reads_df, aes(x = (start + end)/2, fill = group)) +
  geom_histogram(bins = 200) + 
  facet_wrap(~group, ncol = 1, scales = "free_y") +
  theme_classic() +
  labs(x = "Genome Position", y = "Count", title = "CD14 Coverage (Force Pull)")

#==== 🦪 AD ATAC dataset: Skin ====
# Note: Here I thought I found an AD ATAC dataset. However after closer investigation,
# turned out this dataset does not come with necessary arrow files to analyse. 
# Discontinued.
ad_atac <- readRDS("C:/Users/denxiack/Documents/ATAC_10X_TRIAL/GSM6190777_atacseq-lesional.rds/GSM6190777_atacseq-lesional.rds")
ad_atac

Assays(ad_atac) #? You are no seurat obj
class(ad_atac) # You are ArchR ojb...

# No arrow files...
getAvailableMatrices(ad_atac)
getArrowFiles(ad_atac)
