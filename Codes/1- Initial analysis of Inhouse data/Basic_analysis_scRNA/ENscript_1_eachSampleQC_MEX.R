library(Seurat)
library(DoubletFinder)
library(ggplot2)
library(cowplot)
library(Matrix)


## Read raw dataset
matrix_dir<-"G:/Endothelial_cells/Batch01_4R_P01S2206012-220747BCD-sf-1201/data/matrix/GW9_B/"
# !!! Change

barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
barcode.names <- read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
feature.names <- read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
sum(duplicated(feature.names$V2))
# should be 10

feature.names$V2[feature.names$V1=="ENSG00000285053"]<-"GGPS1-TBCE"
feature.names$V2[feature.names$V1=="ENSG00000261186"]<-"LINC01238.2"
feature.names$V2[feature.names$V1=="ENSG00000271858"]<-"CACNA2D2-AS1"
feature.names$V2[feature.names$V1=="ENSG00000280987"]<-"MATR3.2"
feature.names$V2[feature.names$V1=="ENSG00000234229"]<-"ENSG00000234229"
feature.names$V2[feature.names$V1=="ENSG00000284024"]<-"MSANTD7"
feature.names$V2[feature.names$V1=="ENSG00000261480"]<-"GOLGA8M.2"
feature.names$V2[feature.names$V1=="ENSG00000286070"]<-"ENSG00000286070"
feature.names$V2[feature.names$V1=="ENSG00000286237"]<-"ARMCX5-GPRASP2.2"
feature.names$V2[feature.names$V1=="ENSG00000269226"]<-"TMSB15C"

sum(duplicated(feature.names$V2))
# should be 0

colnames(mat) <- barcode.names$V1
rownames(mat) <- feature.names$V2



## Make seurat dataset
sample_id<-"B_GW09_02"
# !!! Change

ncol(mat)
# 17058
DF_nExpRate<- 0.08
# !!! Change according to ncol(mat)
# Reference: https://cdn.10xgenomics.com/image/upload/v1668017706/support-documents/CG000315_ChromiumNextGEMSingleCell3-_GeneExpression_v3.1_DualIndex__RevE.pdf
  

object_raw <- CreateSeuratObject(counts = mat, project = sample_id)
# 36601 features across 17058 samples within 1 assay 


object_raw$Age<-"GW09"
object_raw$Organ<-"Brain"
object_raw$Batch<-"Batch01"
# !!! Change


feats <- c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo", "percent.hb", "percent.plat")
genes.table <- read.csv("D:/OneDrive/KI/Materials/Resources/genes.table.csv", row.names = 1)




## Calculate QC
object_raw <- PercentageFeatureSet(object_raw, "^MT-", col.name = "percent.mito")
object_raw <- PercentageFeatureSet(object_raw, "^RP[SL]", col.name = "percent.ribo")
object_raw <- PercentageFeatureSet(object_raw, "^HB[^(P)]", col.name = "percent.hb")
object_raw <- PercentageFeatureSet(object_raw, "PECAM1|PF4", col.name = "percent.plat")

VlnPlot(object_raw, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) + NoLegend()
ggsave(paste0("QC_metrics_", sample_id, "_1_all.png"), width=8, height=6.5)



object_raw <- NormalizeData(object_raw, verbose = F)

## Sex identification
Sex_genes.table <- genes.table[genes.table$external_gene_name %in% rownames(object_raw),]
Sex_chrY.gene <- Sex_genes.table$external_gene_name[Sex_genes.table$chromosome_name == "Y"]
object_raw$percent.chrY <- colSums(object_raw@assays$RNA@data[Sex_chrY.gene, ])/colSums(object_raw@assays$RNA@data)*100

VlnPlot(object_raw, features = c("XIST", "percent.chrY"), group.by = "orig.ident")
ggsave(paste0("QC_", sample_id, "_sex_identification.png"), width=6, height=4)


object_raw$Inferred_sex<-"Male"
# !!! Change according to the sex identification plot



## Cell cycle
object_raw <- CellCycleScoring(object = object_raw, 
                               g2m.features = cc.genes.updated.2019$g2m.genes, s.features = cc.genes.updated.2019$s.genes)
VlnPlot(object_raw, features = c("S.Score", "G2M.Score"), group.by = "orig.ident")
ggsave(paste0("QC_", sample_id, "_cell_cycle.png"), width=6, height=4)



## Quality check
selected_c <- WhichCells(object_raw, expression = nFeature_RNA > 500 & nFeature_RNA < 10000)
length(selected_c)/ncol(object_raw)
# 0.9913823
# Warn if less than 0.95
object_filter <- subset(object_raw, cells = selected_c)

selected_ribo <- WhichCells(object_filter, expression = percent.ribo > 5)
length(selected_ribo)/ncol(object_filter)
# 0.9329431
# Warn if less than 0.95
object_filter <- subset(object_filter, cells = selected_ribo)

selected_mito <- WhichCells(object_filter, expression = percent.mito < 25)
length(selected_mito)/ncol(object_filter)
# 0.9958167
object_filter <- subset(object_filter, cells = selected_mito)


object_filter
# 36601 features across 15711 samples within 1 assay 
ncol(object_filter)/ncol(object_raw)
# 0.9210341
# Warn if less than 0.9


VlnPlot(object_filter, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) + NoLegend()
ggsave(paste0("QC_metrics_", sample_id, "_2_filter_cells.png"), width=8, height=6.5)




## Predict doublets
object_filter <- FindVariableFeatures(object_filter, verbose = F)
object_filter <- ScaleData(object_filter, vars.to.regress = c("nFeature_RNA", "percent.mito"), verbose = F)
object_filter <- RunPCA(object_filter, verbose = F, npcs = 30)
object_filter <- RunUMAP(object_filter, dims = 1:20, verbose = F)


DF_sweep.res <- paramSweep_v3(object_filter)
DF_sweep.stats <- summarizeSweep(DF_sweep.res, GT = FALSE)
DF_bcmvn <- find.pK(DF_sweep.stats)

png(paste0("DoubetFinder_", sample_id, "_BCmetric.png"), width=600)
barplot(DF_bcmvn$BCmetric, names.arg = DF_bcmvn$pK, las=2)
dev.off()

DF_pK_bcmvn <- as.numeric(as.character(DF_bcmvn$pK[which.max(DF_bcmvn$BCmetric)]))

DF_nExp <- round(ncol(object_filter) * DF_nExpRate)
# 1257

object_filter <- doubletFinder_v3(object_filter, pN = 0.25, pK = DF_pK_bcmvn, nExp = DF_nExp, PCs=1:20)
DF_name <- colnames(object_filter@meta.data)[grepl("DF.classification", colnames(object_filter@meta.data))]

plot_grid(ncol = 2, DimPlot(object_filter, group.by = "Phase") + NoAxes(),
          DimPlot(object_filter, group.by = DF_name) + NoAxes())
ggsave(paste0("DoubetFinder_", sample_id, "_UMAP.png"), width=8, height=4)

VlnPlot(object_filter, features = "nFeature_RNA", group.by = DF_name, pt.size = 0.1)
ggsave(paste0("DoubetFinder_", sample_id, "_nFeature.png"), width=5, height=4)


object_DF <- object_filter[, object_filter@meta.data[, DF_name] == "Singlet"]
# 36601 features across 14454 samples within 1 assay 

ncol(object_DF)/ncol(object_filter)
# 0.9199924
ncol(object_DF)/ncol(object_raw)
# 0.8473444
# Warn if less than 0.8


VlnPlot(object_DF, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) + NoLegend()
ggsave(paste0("QC_metrics_", sample_id, "_3_remove_doublets.png"), width=8, height=6.5)



## Diet & save
object_diet<-DietSeurat(object_DF)
object_diet@meta.data<-object_diet@meta.data[,-c(16,17)]

saveRDS(object_diet, paste0(sample_id, ".rds"))
write.csv(object_diet@meta.data, file=paste0(sample_id, "_metadata.csv"), quote = FALSE)
