library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(biovizBase)
library(ggplot2)
library(cowplot)


# Reference: https://stuartlab.org/signac/articles/pbmc_vignette

setwd("D:/OneDrive/SciLifeLab/Endothelial_cells/20230920scATAC")


sample<-"H_GW10_14a"


## Data Pre-processing and Object Making
raw_data_dir<-"G:/Endothelial_cells/Batch07_2R_temp-20230915111246_ATAC/data/"
raw_file_counts<-paste0(raw_data_dir, "h5/2306022_GW10_H2_ATAC/filtered_peak_bc_matrix.h5")
raw_file_fragments<-paste0(raw_data_dir, "fragments/2306022_GW10_H2_ATAC/fragments.tsv.gz")
raw_file_metadata<-paste0(raw_data_dir, "fragments/2306022_GW10_H2_ATAC/singlecell.csv")

raw_H_GW10_14a_counts <- Read10X_h5(filename = raw_file_counts)
dim(raw_H_GW10_14a_counts)
# 

raw_H_GW10_14a_chrom_assay <- CreateChromatinAssay(
  counts = raw_H_GW10_14a_counts,
  sep = c(":", "-"),
  fragments = raw_file_fragments
)
# ChromatinAssay data with 


# call peaks using MACS3
macs3_peaks <- CallPeaks(
  object = raw_H_GW10_14a_chrom_assay,
  macs2.path = "/path/to/macs3"
)
# GRanges object with

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
macs3_peaks <- keepStandardChromosomes(macs3_peaks, pruning.mode = "coarse")
macs3_peaks <- subsetByOverlaps(x = macs3_peaks, ranges = blacklist_hg38_unified, invert = TRUE)
# GRanges object with


# quantify counts in each peak
macs3_counts <- FeatureMatrix(
  fragments = Fragments(raw_H_GW10_14a_chrom_assay),
  features = macs3_peaks
)
dim(macs3_counts)
# 

# create a new assay using the MACS3 peak set
macs3_chrom_assay <- CreateChromatinAssay(
  counts = macs3_counts,
  sep = c(":", "-"),
  fragments = Fragments(raw_H_GW10_14a_chrom_assay),
  min.cells = 10,
  min.features = 200
)
# ChromatinAssay data with 


# create the Seurat object
raw_H_GW10_14a_metadata <- read.csv(
  file = raw_file_metadata,
  header = TRUE,
  row.names = 1
)
str(raw_H_GW10_14a_metadata)
# 'data.frame':	

raw_H_GW10_14a <- CreateSeuratObject(
  counts = macs3_chrom_assay,
  project = sample,
  assay = "peaks",
  meta.data = raw_H_GW10_14a_metadata[-1,]
)
# An object of class Seurat 



# Get gene annotations
annotations_hg38 <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# add "chr" to the chromosome name so that consistent to input file
seqlevels(annotations_hg38) <- paste0('chr', seqlevels(annotations_hg38))

# add the gene information to the object
Annotation(raw_H_GW10_14a) <- annotations_hg38



## Computing QC Metrics
# compute nucleosome signal score per cell
raw_H_GW10_14a <- NucleosomeSignal(object = raw_H_GW10_14a)

# compute TSS enrichment score per cell
raw_H_GW10_14a <- TSSEnrichment(object = raw_H_GW10_14a, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
raw_H_GW10_14a$pct_reads_in_peaks <- raw_H_GW10_14a$peak_region_fragments / raw_H_GW10_14a$passed_filters * 100
raw_H_GW10_14a$blacklist_ratio <- raw_H_GW10_14a$blacklist_region_fragments / raw_H_GW10_14a$peak_region_fragments


DensityScatter(raw_H_GW10_14a, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
ggsave(paste0("QC_", sample, "_densityScatter.png"), width=6, height=4)

raw_H_GW10_14a$high.tss <- ifelse(raw_H_GW10_14a$TSS.enrichment > 3, 'High', 'Low')
TSSPlot(raw_H_GW10_14a, group.by = 'high.tss') + NoLegend()
ggsave(paste0("QC_", sample, "_TSSPlot.png"), width=6, height=4)


raw_H_GW10_14a$nucleosome_group <- ifelse(raw_H_GW10_14a$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = raw_H_GW10_14a, group.by = 'nucleosome_group')
ggsave(paste0("QC_", sample, "_FragmentHistogram.png"), width=6, height=3)


# plot the distribution of each QC metric separately using a violin plot
VlnPlot(
  object = raw_H_GW10_14a,
  features = c('nCount_peaks', 'peak_region_fragments', 'pct_reads_in_peaks', 'nucleosome_signal', 'TSS.enrichment', 'blacklist_ratio'),
  pt.size = 0.1,
  ncol = 3
)
ggsave(paste0("QC_metrics_", sample, "_1_all.png"), width=8, height=8)


# remove cells that are outliers for these QC metrics.
H_GW10_14a <- subset(
  x = raw_H_GW10_14a,
  subset = nCount_peaks > 2000 &
    nCount_peaks < 75000 &
    peak_region_fragments > 2000 &
    peak_region_fragments < 50000 &
    pct_reads_in_peaks > 20 &
    blacklist_ratio < 0.025 &
    nucleosome_signal < 4 &
    TSS.enrichment > 3
)
H_GW10_14a

ncol(H_GW10_14a)/ncol(raw_H_GW10_14a)
# 

VlnPlot(
  object = H_GW10_14a,
  features = c('nCount_peaks', 'peak_region_fragments', 'pct_reads_in_peaks', 'nucleosome_signal', 'TSS.enrichment', 'blacklist_ratio'),
  pt.size = 0.1,
  ncol = 3
)
ggsave(paste0("QC_metrics_", sample, "_2_filter.png"), width=8, height=8)


saveRDS(H_GW10_14a, paste0(sample, "_QC.rds"))

