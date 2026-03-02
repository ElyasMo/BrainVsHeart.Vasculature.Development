library(Seurat)
library(ggplot2)
library(cowplot)
library(SeuratDisk)
library(harmony)



## Load raw dataset
raw_data_dir<-"/home/liuyang/Bioc/EndothelialCellAtlas/Analysis_Yubo/QC matrix/"
raw_data_files <- list.files(raw_data_dir, ".rds", recursive = TRUE)

raw_data_files<-raw_data_files[!grepl("force", raw_data_files)]

raw_data_rds <- list()

for (i in seq_along(raw_data_files)) {
  raw_data_rds[[i]] <- readRDS(file.path(raw_data_dir, raw_data_files[i]))
}
names(raw_data_rds)<-sapply(raw_data_files, function(x) strsplit(x, "/")[[1]][2])


raw_data_sampleID<-names(raw_data_rds)



for (i in seq_along(raw_data_rds)) {
  raw_data_rds[[i]]$orig.ident <- raw_data_sampleID[i]
}


## Combine samples
EN_Batch01to08 <- merge(raw_data_rds[[1]], raw_data_rds[-1], 
                    add.cell.ids = raw_data_sampleID, project="ENatlas_dev_TT")


rm(raw_data_rds)
gc()

Idents(EN_Batch01to08)<-EN_Batch01to08$orig.ident



## Sample information
VlnPlot(EN_Batch01to08, features = c("XIST", "percent.chrY"), group.by = "orig.ident", pt.size = 0, ncol = 1)
ggsave("QC_Batch01to08_sex_identification.png", width=6, height=6)


sample_sex<-as.data.frame.matrix(table(EN_Batch01to08$orig.ident,  EN_Batch01to08$Inferred_sex))
sample_sex<-transform(sample_sex, Sex=ifelse(sample_sex[,1]>0, "Female", "Male"))


sample_phase<-as.data.frame.matrix(table(EN_Batch01to08$orig.ident,  EN_Batch01to08$Phase))
colnames(sample_phase)<-paste0("Count_", colnames(sample_phase))


raw_data_batch<-substr(sapply(raw_data_files, function(x) strsplit(x, "/")[[1]][1]), 1, 7)
names(raw_data_batch)<-raw_data_sampleID
names(raw_data_files)<-raw_data_sampleID


sample_ID<-names(table(EN_Batch01to08$orig.ident))

sample_info<-table(EN_Batch01to08$orig.ident)
sample_info<-data.frame(Original_ID=names(sample_info), Cell_count=as.integer(sample_info), Inferred_sex=sample_sex$Sex)
sample_info$Gestational_week<-sapply(sample_ID, function(x) strsplit(x, "_")[[1]][2])
sample_info$Organ<-ifelse(sapply(sample_ID, function(x) strsplit(x, "_")[[1]][1])=="B", "Brain", "Heart")
sample_info$Sample_ID<-sapply(sample_ID, function(x) strsplit(x, "_")[[1]][3])
sample_info$Batch_ID<-raw_data_batch[sample_ID]
sample_info<-cbind(sample_info, sample_phase)
sample_info$QC_file<-raw_data_files[sample_ID]


## Quality check
feats <- c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo", "percent.hb", "percent.plat")

VlnPlot(EN_Batch01to08, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 2) +
  NoLegend()

VlnPlot(EN_Batch01to08, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) +
  NoLegend()



## Filter genes
genes_all<-rownames(EN_Batch01to08)  # 36601
genes_selected<-genes_all[grepl("MALAT1", genes_all) | grepl("^RP[SL]", genes_all)
                          | grepl("^HB[^(P)]", genes_all) | grepl("^MT-", genes_all)]  # 129
genes_filter<-setdiff(genes_all, genes_selected)  # 36472

EN_Batch01to08_filter <- subset(EN_Batch01to08, features = genes_filter)


# detection based gene filtering
counts <- GetAssayData(object = EN_Batch01to08_filter, slot = "counts")
nonzero <- counts > 0


selected_f <- rownames(EN_Batch01to08_filter)[Matrix::rowSums(nonzero) >= 30]
EN_Batch01to08_filter <- subset(EN_Batch01to08_filter, features = selected_f)


VlnPlot(EN_Batch01to08_filter, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 2) +
  NoLegend()

VlnPlot(EN_Batch01to08_filter, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) +
  NoLegend()

rm(nonzero, counts)
gc()


## Quality check again
EN_Batch01to08_filter <- NormalizeData(EN_Batch01to08_filter)

EN_Batch01to08_filter <- PercentageFeatureSet(EN_Batch01to08_filter, "^MT-", col.name = "percent.mito")
EN_Batch01to08_filter <- PercentageFeatureSet(EN_Batch01to08_filter, "^HB[^(P)]", col.name = "percent.hb")
EN_Batch01to08_filter <- PercentageFeatureSet(EN_Batch01to08_filter, "^RP[SL]", col.name = "percent.ribo")
EN_Batch01to08_filter <- PercentageFeatureSet(EN_Batch01to08_filter, "PECAM1|PF4", col.name = "percent.plat")

VlnPlot(EN_Batch01to08_filter, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 2) +
  NoLegend()

VlnPlot(EN_Batch01to08_filter, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) +
  NoLegend()




## Analysis
EN_Batch01to08_analysis <- FindVariableFeatures(EN_Batch01to08_filter, selection.method = "vst", nfeatures = 5000)
VariableFeaturePlot(EN_Batch01to08_analysis) + scale_y_continuous(trans="log2")

ggplot(EN_Batch01to08_analysis@assays$RNA@meta.features, aes(x=vst.variable, y=vst.variance.standardized)) + geom_violin()

rm(EN_Batch01to08_filter)
gc()


EN_Batch01to08_analysis <- ScaleData(EN_Batch01to08_analysis, vars.to.regress = c("nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score"))


EN_Batch01to08_analysis <- RunPCA(EN_Batch01to08_analysis, npcs = 100)

VizDimLoadings(EN_Batch01to08_analysis, dims = 1:6, reduction = "pca", ncol = 3)

ElbowPlot(EN_Batch01to08_analysis, reduction = "pca", ndims = 100)


EN_Batch01to08_analysis <- RunUMAP(EN_Batch01to08_analysis, reduction = "pca", dims = 1:40)



## convert to scanpy object
EN_Batch01to08_analysis_diet<-DietSeurat(EN_Batch01to08_analysis, dimreducs = names(EN_Batch01to08_analysis@reductions))

SaveH5Seurat(EN_Batch01to08_analysis_diet, filename = "EN_Batch01to08_3_analysis.h5Seurat")

Convert("EN_Batch01to08_3_analysis.h5Seurat", dest = "h5ad")



## Harmony Integration
EN_Batch01to08_integrated <- RunHarmony(EN_Batch01to08_analysis, group.by.vars = "orig.ident", reduction = "pca",
                                    dims.use = 1:100, assay.use = "RNA", reduction.save="harmony", 
                                    max.iter.harmony=25, max.iter.cluster = 25, plot_convergence = TRUE)
# Harmony converged after 5 iterations
ElbowPlot(EN_Batch01to08_integrated, reduction = "harmony", ndims = 100)


EN_Batch01to08_integrated <- RunUMAP(EN_Batch01to08_integrated, dims = 1:40, n.neighbors = 20, reduction = "harmony", reduction.name = "harmony_umap_d40k20")



saveRDS(EN_Batch01to08_integrated, "EN_Batch01to08_4_integration_Harmony_dims100_UMAP_d40k20.rds")

