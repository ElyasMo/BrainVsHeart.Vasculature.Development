library(Seurat)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(dplyr)
library(dendsort)



## Load data
EN_Batch01to08_Vasculature <- readRDS("EN_Batch01to08_annotated_L2_Vasculature.rds")


############## L3: ECs
Key(EN_Batch01to08_Vasculature@reductions$harmony_umap_L1)<-"harmonyUmapL1_"
Key(EN_Batch01to08_Vasculature@reductions$harmony_umap_L2)<-"harmonyUmapL2_"

EN_Batch01to08_ECs<-subset(EN_Batch01to08_Vasculature, L2_Cell_type=="ECs")
EN_Batch01to08_ECs$L2_RNA_snn_res.0.5<-droplevels(EN_Batch01to08_ECs$L2_RNA_snn_res.0.5)


EN_Batch01to08_ECs <- RunUMAP(EN_Batch01to08_ECs, reduction = "harmony", dims = 1:40, n.neighbors = 20, reduction.name = "harmony_umap_L3", reduction.key = "harmonyUmapL3_")



## Clustering
EN_Batch01to08_ECs <- FindNeighbors(EN_Batch01to08_ECs, dims = 1:40, k.param = 20, reduction = "harmony")

EN_Batch01to08_ECs <- FindClusters(EN_Batch01to08_ECs, graph.name = "RNA_snn", resolution = 0.7, algorithm = 1)

EN_Batch01to08_ECs$RNA_snn_res.0.7<-factor(EN_Batch01to08_ECs$RNA_snn_res.0.7, levels = seq(0, length(unique(EN_Batch01to08_ECs$RNA_snn_res.0.7))-1))
colnames(EN_Batch01to08_ECs@meta.data)[20]<-"L3_RNA_snn_res.0.7"


############## annotation
annotation_L3_ECs_re<-setNames(rep("", length(table(EN_Batch01to08_ECs$L3_RNA_snn_res.0.7))), names(table(EN_Batch01to08_ECs$L3_RNA_snn_res.0.7)))
annotation_L3_ECs_re<-annotate_by_clusterN(annotation_L3_ECs_re, "aEC_1", c(1))
annotation_L3_ECs_re<-annotate_by_clusterN(annotation_L3_ECs_re, "aEC_2", c(11))
annotation_L3_ECs_re<-annotate_by_clusterN(annotation_L3_ECs_re, "vEC_1", c(6))
annotation_L3_ECs_re<-annotate_by_clusterN(annotation_L3_ECs_re, "vEC_2", c(17))
annotation_L3_ECs_re<-annotate_by_clusterN(annotation_L3_ECs_re, "cEC_1", c(4))
annotation_L3_ECs_re<-annotate_by_clusterN(annotation_L3_ECs_re, "cEC_2", c(5))
annotation_L3_ECs_re<-annotate_by_clusterN(annotation_L3_ECs_re, "cEC_3", c(8))
annotation_L3_ECs_re<-annotate_by_clusterN(annotation_L3_ECs_re, "cEC_4", c(18))
annotation_L3_ECs_re<-annotate_by_clusterN(annotation_L3_ECs_re, "tip_1", c(7))
annotation_L3_ECs_re<-annotate_by_clusterN(annotation_L3_ECs_re, "tip_2", c(14))
annotation_L3_ECs_re<-annotate_by_clusterN(annotation_L3_ECs_re, "endo_like", c(2))
annotation_L3_ECs_re<-annotate_by_clusterN(annotation_L3_ECs_re, "precursor", c(16))
annotation_L3_ECs_re<-annotate_by_clusterN(annotation_L3_ECs_re, "EC_mito", c(13))
annotation_L3_ECs_re<-annotate_by_clusterN(annotation_L3_ECs_re, "endo_1", c(15))
annotation_L3_ECs_re<-annotate_by_clusterN(annotation_L3_ECs_re, "endo_2", c(20))
annotation_L3_ECs_re<-annotate_by_clusterN(annotation_L3_ECs_re, "EC_lym", c(12))
annotation_L3_ECs_re<-annotate_by_clusterN(annotation_L3_ECs_re, "debri", c(0, 3, 9, 19))
annotation_L3_ECs_re<-annotate_by_clusterN(annotation_L3_ECs_re, "doublet", c(10))


annotation_L3_ECs_re_merge <- annotation_L3_ECs_re
annotation_L3_ECs_re_merge[grep("aEC", annotation_L3_ECs_re_merge)] <- "Arterial"
annotation_L3_ECs_re_merge[grep("vEC", annotation_L3_ECs_re_merge)] <- "Venous"
annotation_L3_ECs_re_merge[grep("cEC", annotation_L3_ECs_re_merge)] <- "Capillary"
annotation_L3_ECs_re_merge[grep("tip", annotation_L3_ECs_re_merge)] <- "Tip cells"
annotation_L3_ECs_re_merge[grep("mito", annotation_L3_ECs_re_merge)] <- "Mitotic"
annotation_L3_ECs_re_merge[grep("pre", annotation_L3_ECs_re_merge)] <- "Precursors"
annotation_L3_ECs_re_merge[grep("endo", annotation_L3_ECs_re_merge)] <- "Endocardial"
annotation_L3_ECs_re_merge[grep("lym", annotation_L3_ECs_re_merge)] <- "Lymphatic"

data.frame(Cell_types = annotation_L3_ECs_re_merge, Subtypes = annotation_L3_ECs_re)

EN_Batch01to08_ECs$L3_ECs_Cell_types<-annotation_L3_ECs_re_merge[EN_Batch01to08_ECs$L3_RNA_snn_res.0.7]

EN_Batch01to08_ECs$L3_ECs_Subtypes<-annotation_L3_ECs_re[EN_Batch01to08_ECs$L3_RNA_snn_res.0.7]

saveRDS(EN_Batch01to08_ECs, "EN_Batch01to08_5_Vasculature_ECs_reanno.rds")



############################# Function #############################

annotate_by_clusterN<-function(annotation_vector, cell_type, cluster_names, ignore=FALSE) {
  cluster_names<-as.character(cluster_names)
  tmp_test<-annotation_vector[cluster_names]
  if(all(unique(tmp_test)=="")) {
    annotation_vector[cluster_names]<-cell_type
  } else {
    warning("Please_check:", paste(names(tmp_test[tmp_test!=""]), collapse = ", "))
    if (ignore) {
      annotation_vector[names(tmp_test[tmp_test==""])]<-cell_type
    }
  }
  return(annotation_vector)
}
