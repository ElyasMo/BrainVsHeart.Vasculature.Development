library(Seurat)
library(ggplot2)
library(cowplot)
library(clustree)
library(RColorBrewer)
library(dplyr)
library(clusterProfiler)


## Load data
EN_Batch01to08_integrated <- readRDS("EN_Batch01to08_4_integration_Harmony_dims100_UMAP_d40k20.rds")



############### L1: CNS major cell types
## Clustering
ElbowPlot(EN_Batch01to08_integrated, reduction = "harmony", ndims = 100)

EN_Batch01to08_final <- FindNeighbors(EN_Batch01to08_integrated, dims = 1:40, k.param = 15, reduction = "harmony")

EN_Batch01to08_final <- FindClusters(EN_Batch01to08_final, graph.name = "RNA_snn", resolution = 0.35, algorithm = 1)


## resolution = 0.35
EN_Batch01to08_final$RNA_snn_res.0.35<-factor(EN_Batch01to08_final$RNA_snn_res.0.35, levels = seq(0, length(unique(EN_Batch01to08_final$RNA_snn_res.0.35))-1))

EN_Batch01to08_final <- SetIdent(EN_Batch01to08_final, value = "RNA_snn_res.0.35")



## FindAllMarkers
L1_markers <- FindAllMarkers(EN_Batch01to08_final, only.pos = TRUE)


L1_markers %>%
  group_by(cluster) %>%
  top_n(-100, p_val_adj) %>%
  top_n(50, pct.1-pct.2) %>%
  top_n(10, avg_log2FC) -> L1_top10
table(L1_top10$cluster)


annotation_L1_DEGs<-setNames(rep("", length(table(Idents(EN_Batch01to08_final)))), names(table(Idents(EN_Batch01to08_final))))
annotation_L1_DEGs<-annotate_by_clusterN(annotation_L1_DEGs, "Debri", c(2, 21))
annotation_L1_DEGs<-annotate_by_clusterN(annotation_L1_DEGs, "ECs", c(0, 4, 6, 10, 17))
annotation_L1_DEGs<-annotate_by_clusterN(annotation_L1_DEGs, "Mural", c(1, 12, 20, 22))
annotation_L1_DEGs<-annotate_by_clusterN(annotation_L1_DEGs, "Vasculature", c(13, 14))
annotation_L1_DEGs<-annotate_by_clusterN(annotation_L1_DEGs, "Neurons", c(3))
annotation_L1_DEGs<-annotate_by_clusterN(annotation_L1_DEGs, "NPCs", c(5, 18))
annotation_L1_DEGs<-annotate_by_clusterN(annotation_L1_DEGs, "EPCs", c(8))
annotation_L1_DEGs<-annotate_by_clusterN(annotation_L1_DEGs, "Immune", c(7, 19))
annotation_L1_DEGs<-annotate_by_clusterN(annotation_L1_DEGs, "HSCs", c(9, 16))
annotation_L1_DEGs<-annotate_by_clusterN(annotation_L1_DEGs, "CMCs", c(11, 15))


annotation_L1_DEGs[annotation_L1_DEGs=="ECs"]<-"Vasculature"
annotation_L1_DEGs[annotation_L1_DEGs=="Mural"]<-"Vasculature"


EN_Batch01to08_final$L1_Cell_type<-annotation_L1_DEGs[EN_Batch01to08_final$RNA_snn_res.0.35]
saveRDS(EN_Batch01to08_final, "EN_Batch01to08_annotated_L1.rds")





############### L2: vasculature cell types
# Make L2 seurat object
EN_Batch01to08_Vasculature<-subset(EN_Batch01to08_final, L1_Cell_type=="Vasculature")

colnames(EN_Batch01to08_Vasculature@meta.data)[16]<-"L1_RNA_snn_res.0.35"
EN_Batch01to08_Vasculature$L1_RNA_snn_res.0.35<-droplevels(EN_Batch01to08_Vasculature$L1_RNA_snn_res.0.35)
names(EN_Batch01to08_Vasculature@reductions)[2]<-"harmony_umap_L1"


# L2 UMAP
EN_Batch01to08_Vasculature <- RunUMAP(EN_Batch01to08_Vasculature, dims = 1:40, n.neighbors = 20, reduction = "harmony", reduction.name = "harmony_umap_L2_d40k20")


## Clustering
EN_Batch01to08_Vasculature <- FindNeighbors(EN_Batch01to08_Vasculature, dims = 1:40, k.param = 15, reduction = "harmony")

EN_Batch01to08_Vasculature <- FindClusters(EN_Batch01to08_Vasculature, graph.name = "RNA_snn", resolution = 0.5, algorithm = 1)



## resolution = 0.5
EN_Batch01to08_Vasculature <- SetIdent(EN_Batch01to08_Vasculature, value = "RNA_snn_res.0.5")



## FindAllMarkers
L2_markers <- FindAllMarkers(EN_Batch01to08_Vasculature, only.pos = TRUE)


## GO analysis
GO_database <- 'org.Hs.eg.db'

annotation_L2_markers_gene_list<-list()
annotation_L2_markers_id_list<-list()
annotation_L2_markers_GO_list<-list()
annotation_L2_markers_GOs_list<-list()

for (tmp in levels(EN_Batch01to08_Vasculature$RNA_snn_res.0.5)) {
  annotation_L2_markers_gene_list[[tmp]]<-rownames(L2_markers)[L2_markers$cluster==tmp]
  annotation_L2_markers_id_list[[tmp]]<- bitr(annotation_L2_markers_gene_list[[tmp]],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  annotation_L2_markers_GO_list[[tmp]]<-enrichGO(annotation_L2_markers_id_list[[tmp]]$ENTREZID, OrgDb = GO_database, keyType = "ENTREZID",
                                                 ont = "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = T)
}

for (tmp in levels(EN_Batch01to08_Vasculature$RNA_snn_res.0.5)) {
  re = try({annotation_L2_markers_GOs_list[[tmp]] = simplify(annotation_L2_markers_GO_list[[tmp]])}, silent = TRUE)
  if (inherits(re, 'try-error')) {
    annotation_L2_markers_GOs_list[[tmp]] = 'Troubles here. Possibly no significant GO term enriched.'
  } else {
    barplot(annotation_L2_markers_GOs_list[[tmp]], split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")
    ggsave(paste0("GO_annotation_L2_markers_", tmp, ".png"), width=7, height=8)
  }
}


# Annotation based on DEGs
annotation_L2_DEGs<-setNames(rep("", length(table(Idents(EN_Batch01to08_Vasculature)))), names(table(Idents(EN_Batch01to08_Vasculature))))
annotation_L2_DEGs<-annotate_by_clusterN(annotation_L2_DEGs, "Debri", c(8, 13, 14, 20))
annotation_L2_DEGs<-annotate_by_clusterN(annotation_L2_DEGs, "Doublet", c(9, 11))
annotation_L2_DEGs<-annotate_by_clusterN(annotation_L2_DEGs, "EC_arterial", c(5))
annotation_L2_DEGs<-annotate_by_clusterN(annotation_L2_DEGs, "EC_ven_cap", c(2, 3, 17))
annotation_L2_DEGs<-annotate_by_clusterN(annotation_L2_DEGs, "EC_capillary", c(1, 4, 18))
annotation_L2_DEGs<-annotate_by_clusterN(annotation_L2_DEGs, "EC_lym", c(12))
annotation_L2_DEGs<-annotate_by_clusterN(annotation_L2_DEGs, "Fibs", c(0))
annotation_L2_DEGs<-annotate_by_clusterN(annotation_L2_DEGs, "SMCs", c(19))
annotation_L2_DEGs<-annotate_by_clusterN(annotation_L2_DEGs, "Mural", c(7, 16))
annotation_L2_DEGs<-annotate_by_clusterN(annotation_L2_DEGs, "Brain_unknown", c(21, 23))
annotation_L2_DEGs<-annotate_by_clusterN(annotation_L2_DEGs, "Heart_unknown", c(15, 22))
annotation_L2_DEGs<-annotate_by_clusterN(annotation_L2_DEGs, "Heart_endocardial", c(6))
annotation_L2_DEGs<-annotate_by_clusterN(annotation_L2_DEGs, "Heart_epi_pro", c(10))
# Heart_epi_pro: epicardial progenitor cells


annotation_L2_DEGs[grep("EC_", annotation_L2_DEGs)]<-"ECs"
annotation_L2_DEGs[annotation_L2_DEGs %in% c("Fibs", "SMCs")]<-"Mural"
annotation_L2_DEGs[annotation_L2_DEGs == "Heart_endocardial"]<-"ECs"
annotation_L2_DEGs[annotation_L2_DEGs == "Heart_epi_pro"]<-"Mural"
annotation_L2_DEGs[c("21", "23")]<-"Mural"
annotation_L2_DEGs[c("15")]<-"Mural"
annotation_L2_DEGs[c("22")]<-"ECs"



EN_Batch01to08_Vasculature$L2_Cell_type<-annotation_L2_DEGs[EN_Batch01to08_Vasculature$RNA_snn_res.0.5]
colnames(EN_Batch01to08_Vasculature@meta.data)[18]<-"L2_RNA_snn_res.0.5"
names(EN_Batch01to08_Vasculature@reductions)[3]<-"harmony_umap_L2"


saveRDS(EN_Batch01to08_Vasculature, "EN_Batch01to08_annotated_L2_Vasculature.rds")



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
