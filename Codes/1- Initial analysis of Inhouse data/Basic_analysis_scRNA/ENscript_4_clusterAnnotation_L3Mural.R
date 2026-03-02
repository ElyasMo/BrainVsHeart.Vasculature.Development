library(Seurat)
library(ggplot2)
library(cowplot)



## Load data
EN_Batch01to08_Vasculature <- readRDS("EN_Batch01to08_annotated_L2_Vasculature.rds")
Key(EN_Batch01to08_Vasculature@reductions$harmony_umap_L1)<-"harmonyUmapL1_"
Key(EN_Batch01to08_Vasculature@reductions$harmony_umap_L2)<-"harmonyUmapL2_"

EN_Batch01to08_Mural<-subset(EN_Batch01to08_Vasculature, L2_Cell_type=="Mural")

rm(EN_Batch01to08_Vasculature)
gc()

EN_Batch01to08_Mural$L2_RNA_snn_res.0.5<-droplevels(EN_Batch01to08_Mural$L2_RNA_snn_res.0.5)


## Analysis

EN_Batch01to08_analysis <- RunPCA(EN_Batch01to08_Mural, npcs = 100)

EN_Batch01to08_analysis <- RunUMAP(EN_Batch01to08_analysis, reduction = "pca", dims = 1:40, reduction.name = "pca_umap", reduction.key = "pcaUmap_")



############################## CCA Integration
## Load data
EN_Batch01to08_analysis_diet <- DietSeurat(EN_Batch01to08_analysis)


alldata.list <- SplitObject(EN_Batch01to08_analysis_diet, split.by = "orig.ident")

for (i in 1:length(alldata.list)) {
  alldata.list[[i]] <- NormalizeData(alldata.list[[i]], verbose = FALSE)
  alldata.list[[i]] <- FindVariableFeatures(alldata.list[[i]], selection.method = "vst",
                                            nfeatures = 2000, verbose = FALSE)
}

alldata.anchors <- FindIntegrationAnchors(object.list = alldata.list, dims = 1:50, reduction = "cca")

alldata.int <- IntegrateData(anchorset = alldata.anchors, dims = 1:50, new.assay.name = "CCA", k.weight = 80)

alldata.int <- ScaleData(alldata.int, verbose = FALSE)
alldata.int <- RunPCA(alldata.int, npcs = 100, verbose = FALSE)
ElbowPlot(alldata.int, reduction = "pca", ndims = 100)
alldata.int <- RunUMAP(alldata.int, dims = 1:50)


saveRDS(alldata.int, "EN_Batch01to08_Mural_integration_CCA.rds")



############################## Annotation
EN_Batch01to08_Mural<-readRDS("EN_Batch01to08_Mural_integration_CCA.rds")



## Clustering
EN_Batch01to08_Mural <- FindNeighbors(EN_Batch01to08_Mural, dims = 1:40, k.param = 20, reduction = "pca")
names(EN_Batch01to08_Mural@graphs)
# "CCA_nn"  "CCA_snn"

EN_Batch01to08_Mural <- FindClusters(EN_Batch01to08_Mural, graph.name = "CCA_snn", resolution = 0.6, algorithm = 1)

EN_Batch01to08_Mural$CCA_snn_res.0.6<-factor(EN_Batch01to08_Mural$CCA_snn_res.0.6, levels = seq(0, length(unique(EN_Batch01to08_Mural$CCA_snn_res.0.6))-1))

EN_Batch01to08_Mural <- SetIdent(EN_Batch01to08_Mural, value = "CCA_snn_res.0.6")


EN_Batch01to08_Mural <- RunUMAP(EN_Batch01to08_Mural, reduction = "pca", dims = 1:40, n.neighbors = 20, reduction.name = "umap")


DefaultAssay(EN_Batch01to08_Mural) <- "RNA"



## FindAllMarkers
L3_Mural_markers <- FindAllMarkers(EN_Batch01to08_Mural, only.pos = TRUE, assay = "RNA")
write.table(L3_Mural_markers, "DEGs_L3_Mural_cluster.txt", sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)



# Annotation based on DEGs & GO
annotation_L3_Mural_DEGs<-setNames(rep("", length(table(EN_Batch01to08_Mural$CCA_snn_res.0.6))), names(table(EN_Batch01to08_Mural$CCA_snn_res.0.6)))
annotation_L3_Mural_DEGs<-annotate_by_clusterN(annotation_L3_Mural_DEGs, "Debri", c(11, 12))
annotation_L3_Mural_DEGs<-annotate_by_clusterN(annotation_L3_Mural_DEGs, "SMCs", c(9))
annotation_L3_Mural_DEGs<-annotate_by_clusterN(annotation_L3_Mural_DEGs, "PCs", c(3))
annotation_L3_Mural_DEGs<-annotate_by_clusterN(annotation_L3_Mural_DEGs, "Fibs", c(0, 2, 8))
annotation_L3_Mural_DEGs<-annotate_by_clusterN(annotation_L3_Mural_DEGs, "Mitotic", c(4, 5, 6))
annotation_L3_Mural_DEGs<-annotate_by_clusterN(annotation_L3_Mural_DEGs, "Precursor", c(7))
annotation_L3_Mural_DEGs<-annotate_by_clusterN(annotation_L3_Mural_DEGs, "Responsive", c(1))
annotation_L3_Mural_DEGs<-annotate_by_clusterN(annotation_L3_Mural_DEGs, "Mural", c(10))


EN_Batch01to08_Mural$L3_Mural_Cell_type<-annotation_L3_Mural_DEGs[EN_Batch01to08_Mural$CCA_snn_res.0.6]



## Remove debri 
EN_Batch01to08_Mural@meta.data[,21]<-NULL
colnames(EN_Batch01to08_Mural@meta.data)[20]<-"L3_CCA_snn_res.0.6"

EN_Batch01to08_Mural_pure <- subset(EN_Batch01to08_Mural, L3_Mural_Cell_type != "Debri")



EN_Batch01to08_Mural_pure$L3_CCA_snn_res.0.6<-droplevels(EN_Batch01to08_Mural_pure$L3_CCA_snn_res.0.6)


EN_Batch01to08_Mural_pure <- RunUMAP(EN_Batch01to08_Mural_pure, reduction = "pca", dims = 1:40, n.neighbors = 20, reduction.name = "umap_pure", reduction.key = "harmonyUmapL3pure_")


saveRDS(EN_Batch01to08_Mural_pure, "EN_Batch01to08_5_Vasculature_Mural_pure.rds")




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
