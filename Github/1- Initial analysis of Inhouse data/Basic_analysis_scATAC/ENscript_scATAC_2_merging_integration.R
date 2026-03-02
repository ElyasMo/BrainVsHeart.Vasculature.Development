library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(biovizBase)
library(ggplot2)
library(cowplot)




################################## Merging

raw_file_fragments_B_GW10_14a<-"~/Bioc/EndothelialCellAtlas/Raw_data/scATAC/Batch07/fragments/2306022_GW10_B2_ATAC/fragments.tsv.gz"
raw_file_fragments_H_GW10_14a<-"~/Bioc/EndothelialCellAtlas/Raw_data/scATAC/Batch07/fragments/2306022_GW10_H2_ATAC/fragments.tsv.gz"
raw_file_fragments_B_GW12_15a<-"~/Bioc/EndothelialCellAtlas/Raw_data/scATAC/Batch09/fragments/2306022_gw12Bra2_ATAC/fragments.tsv.gz"
raw_file_fragments_H_GW12_15a<-"~/Bioc/EndothelialCellAtlas/Raw_data/scATAC/Batch09/fragments/2306022_gw12Hea2_ATAC/fragments.tsv.gz"
raw_file_fragments_B_GW08_16a<-"~/Bioc/EndothelialCellAtlas/Raw_data/scATAC/Batch10/fragments/2306022_GW8_Bra_3_ATAC/fragments.tsv.gz"
raw_file_fragments_H_GW08_16a<-"~/Bioc/EndothelialCellAtlas/Raw_data/scATAC/Batch10/fragments/2306022_GW8_Hea_3_ATAC/fragments.tsv.gz"


## Creating a common peak set

# read in QCed object
B_GW10_14a <- readRDS("../20231012Batch7_scATAC/B_GW10_14a_QC.rds")
H_GW10_14a <- readRDS("../20231012Batch7_scATAC/H_GW10_14a_QC.rds")
B_GW12_15a <- readRDS("../20241125scATAC_2batches/B_GW12_15a_QC.rds")
H_GW12_15a <- readRDS("../20241125scATAC_2batches/H_GW12_15a_QC.rds")
B_GW08_16a <- readRDS("../20241125scATAC_2batches/B_GW08_16a_QC.rds")
H_GW08_16a <- readRDS("../20241125scATAC_2batches/H_GW08_16a_QC.rds")

# convert to genomic ranges
gr_B_GW10_14a <- granges(B_GW10_14a)
gr_H_GW10_14a <- granges(H_GW10_14a)
gr_B_GW12_15a <- granges(B_GW12_15a)
gr_H_GW12_15a <- granges(H_GW12_15a)
gr_B_GW08_16a <- granges(B_GW08_16a)
gr_H_GW08_16a <- granges(H_GW08_16a)


# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(gr_B_GW10_14a, gr_H_GW10_14a, gr_B_GW12_15a, gr_H_GW12_15a, gr_B_GW08_16a, gr_H_GW08_16a))
  

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)


combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]



## Create Fragment objects

# load metadata
md_B_GW10_14a<-B_GW10_14a@meta.data
md_H_GW10_14a<-H_GW10_14a@meta.data
md_B_GW12_15a<-B_GW12_15a@meta.data
md_H_GW12_15a<-H_GW12_15a@meta.data
md_B_GW08_16a<-B_GW08_16a@meta.data
md_H_GW08_16a<-H_GW08_16a@meta.data


# perform an initial filtering of low count cells
md_B_GW10_14a <- md_B_GW10_14a[md_B_GW10_14a$passed_filters > 500, ]
md_H_GW10_14a <- md_H_GW10_14a[md_H_GW10_14a$passed_filters > 500, ]
md_B_GW12_15a <- md_B_GW12_15a[md_B_GW12_15a$passed_filters > 500, ]
md_H_GW12_15a <- md_H_GW12_15a[md_H_GW12_15a$passed_filters > 500, ]
md_B_GW08_16a <- md_B_GW08_16a[md_B_GW08_16a$passed_filters > 500, ]
md_H_GW08_16a <- md_H_GW08_16a[md_H_GW08_16a$passed_filters > 500, ]

# create fragment objects
frags_B_GW10_14a <- CreateFragmentObject(
  path = raw_file_fragments_B_GW10_14a,
  cells = rownames(md_B_GW10_14a)
)

frags_H_GW10_14a <- CreateFragmentObject(
  path = raw_file_fragments_H_GW10_14a,
  cells = rownames(md_H_GW10_14a)
)

frags_B_GW12_15a <- CreateFragmentObject(
  path = raw_file_fragments_B_GW12_15a,
  cells = rownames(md_B_GW12_15a)
)

frags_H_GW12_15a <- CreateFragmentObject(
  path = raw_file_fragments_H_GW12_15a,
  cells = rownames(md_H_GW12_15a)
)

frags_B_GW08_16a <- CreateFragmentObject(
  path = raw_file_fragments_B_GW08_16a,
  cells = rownames(md_B_GW08_16a)
)

frags_H_GW08_16a <- CreateFragmentObject(
  path = raw_file_fragments_H_GW08_16a,
  cells = rownames(md_H_GW08_16a)
)



## Quantify peaks in each dataset
counts_B_GW10_14a <- FeatureMatrix(
  fragments = frags_B_GW10_14a,
  features = combined.peaks,
  cells = rownames(md_B_GW10_14a)
)

counts_H_GW10_14a <- FeatureMatrix(
  fragments = frags_H_GW10_14a,
  features = combined.peaks,
  cells = rownames(md_H_GW10_14a)
)

counts_B_GW12_15a <- FeatureMatrix(
  fragments = frags_B_GW12_15a,
  features = combined.peaks,
  cells = rownames(md_B_GW12_15a)
)

counts_H_GW12_15a <- FeatureMatrix(
  fragments = frags_H_GW12_15a,
  features = combined.peaks,
  cells = rownames(md_H_GW12_15a)
)

counts_B_GW08_16a <- FeatureMatrix(
  fragments = frags_B_GW08_16a,
  features = combined.peaks,
  cells = rownames(md_B_GW08_16a)
)

counts_H_GW08_16a <- FeatureMatrix(
  fragments = frags_H_GW08_16a,
  features = combined.peaks,
  cells = rownames(md_H_GW08_16a)
)


## Create the objects
assay_B_GW10_14a <- CreateChromatinAssay(counts_B_GW10_14a, fragments = frags_B_GW10_14a)
obj_B_GW10_14a <- CreateSeuratObject(assay_B_GW10_14a, assay = "ATAC", meta.data=md_B_GW10_14a)

assay_H_GW10_14a <- CreateChromatinAssay(counts_H_GW10_14a, fragments = frags_H_GW10_14a)
obj_H_GW10_14a <- CreateSeuratObject(assay_H_GW10_14a, assay = "ATAC", meta.data=md_H_GW10_14a)

assay_B_GW12_15a <- CreateChromatinAssay(counts_B_GW12_15a, fragments = frags_B_GW12_15a)
obj_B_GW12_15a <- CreateSeuratObject(assay_B_GW12_15a, assay = "ATAC", meta.data=md_B_GW12_15a)

assay_H_GW12_15a <- CreateChromatinAssay(counts_H_GW12_15a, fragments = frags_H_GW12_15a)
obj_H_GW12_15a <- CreateSeuratObject(assay_H_GW12_15a, assay = "ATAC", meta.data=md_H_GW12_15a)

assay_B_GW08_16a <- CreateChromatinAssay(counts_B_GW08_16a, fragments = frags_B_GW08_16a)
obj_B_GW08_16a <- CreateSeuratObject(assay_B_GW08_16a, assay = "ATAC", meta.data=md_B_GW08_16a)

assay_H_GW08_16a <- CreateChromatinAssay(counts_H_GW08_16a, fragments = frags_H_GW08_16a)
obj_H_GW08_16a <- CreateSeuratObject(assay_H_GW08_16a, assay = "ATAC", meta.data=md_H_GW08_16a)



## Merge objects

# add information to identify dataset of origin
obj_B_GW10_14a$dataset <- 'B_GW10_14a'
obj_H_GW10_14a$dataset <- 'H_GW10_14a'
obj_B_GW12_15a$dataset <- 'B_GW12_15a'
obj_H_GW12_15a$dataset <- 'H_GW12_15a'
obj_B_GW08_16a$dataset <- 'B_GW08_16a'
obj_H_GW08_16a$dataset <- 'H_GW08_16a'


# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
  x = obj_B_GW10_14a,
  y = list(obj_H_GW10_14a, obj_B_GW12_15a, obj_H_GW12_15a, obj_B_GW08_16a, obj_H_GW08_16a),
  add.cell.ids = c("B_GW10_14a", "H_GW10_14a", "B_GW12_15a", "H_GW12_15a", "B_GW08_16a", "H_GW08_16a")
)



# normalization
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')

DimPlot(combined, group.by = 'dataset', pt.size = 0.1)

saveRDS(combined, "scATAC_all6samples_merged.rds")




################################## Integration
p1 <-DimPlot(combined, group.by = 'dataset', pt.size = 0.1, shuffle = TRUE)


## Preprocessing

# compute LSI
obj_B_GW10_14a <- RunTFIDF(obj_B_GW10_14a)
obj_B_GW10_14a <- FindTopFeatures(obj_B_GW10_14a, min.cutoff = 10)
obj_B_GW10_14a <- RunSVD(obj_B_GW10_14a)
obj_B_GW10_14a<-RenameCells(obj_B_GW10_14a, add.cell.id = "B_GW10_14a")

obj_H_GW10_14a <- RunTFIDF(obj_H_GW10_14a)
obj_H_GW10_14a <- FindTopFeatures(obj_H_GW10_14a, min.cutoff = 10)
obj_H_GW10_14a <- RunSVD(obj_H_GW10_14a)
obj_H_GW10_14a<-RenameCells(obj_H_GW10_14a, add.cell.id = "H_GW10_14a")

obj_B_GW12_15a <- RunTFIDF(obj_B_GW12_15a)
obj_B_GW12_15a <- FindTopFeatures(obj_B_GW12_15a, min.cutoff = 10)
obj_B_GW12_15a <- RunSVD(obj_B_GW12_15a)
obj_B_GW12_15a<-RenameCells(obj_B_GW12_15a, add.cell.id = "B_GW12_15a")

obj_H_GW12_15a <- RunTFIDF(obj_H_GW12_15a)
obj_H_GW12_15a <- FindTopFeatures(obj_H_GW12_15a, min.cutoff = 10)
obj_H_GW12_15a <- RunSVD(obj_H_GW12_15a)
obj_H_GW12_15a<-RenameCells(obj_H_GW12_15a, add.cell.id = "H_GW12_15a")

obj_B_GW08_16a <- RunTFIDF(obj_B_GW08_16a)
obj_B_GW08_16a <- FindTopFeatures(obj_B_GW08_16a, min.cutoff = 10)
obj_B_GW08_16a <- RunSVD(obj_B_GW08_16a)
obj_B_GW08_16a<-RenameCells(obj_B_GW08_16a, add.cell.id = "B_GW08_16a")

obj_H_GW08_16a <- RunTFIDF(obj_H_GW08_16a)
obj_H_GW08_16a <- FindTopFeatures(obj_H_GW08_16a, min.cutoff = 10)
obj_H_GW08_16a <- RunSVD(obj_H_GW08_16a)
obj_H_GW08_16a<-RenameCells(obj_H_GW08_16a, add.cell.id = "H_GW08_16a")



## Integration

# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = list(obj_B_GW10_14a, obj_H_GW10_14a, obj_B_GW12_15a, obj_H_GW12_15a, obj_B_GW08_16a, obj_H_GW08_16a),
  anchor.features = rownames(combined),
  reduction = "rlsi",
  dims = 2:30
)


# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = combined[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)
p2 <- DimPlot(integrated, group.by = "dataset", pt.size = 0.1, shuffle = TRUE)


(p1 + ggtitle("Merged")) | (p2 + ggtitle("rLSI integrated"))
ggsave("scATAC_all6samples_integration_rLSI.png", width=10, height=4)



saveRDS(integrated, "scATAC_all6samples_integrated_rLSI.rds")



## Harmony integration
library(harmony)

combined.integrated <- RunHarmony(
  object = combined,
  group.by.vars = 'dataset',
  reduction = 'lsi',
  project.dim = FALSE
)


# re-compute the UMAP using corrected LSI embeddings
combined.integrated <- RunUMAP(combined.integrated, dims = 2:30, reduction = 'harmony')

saveRDS(combined.integrated, "scATAC_all6samples_integrated_Harmony.rds")

