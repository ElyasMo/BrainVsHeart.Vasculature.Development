library(Seurat)
library(ggplot2)
library(cowplot)



############################## scRNA-seq Combine annotation
## Load data
EN_Batch01to08_final <- readRDS("EN_Batch01to08_annotated_L1.rds")

############### L1: Major cell types
colnames(EN_Batch01to08_final@meta.data)[16]<-"L1_RNA_snn_res.0.35"

EN_Batch01to08_final@meta.data[EN_Batch01to08_final$L1_Cell_type == "CMCs", "L1_Cell_type"] <- "CMs"
EN_Batch01to08_final@meta.data[EN_Batch01to08_final$L1_Cell_type == "EPCs", "L1_Cell_type"] <- "NPCs"

names(EN_Batch01to08_final@reductions)[2]<-"harmony_umap_L1"



############### L2: vasculature cell types
EN_Batch01to08_Vasculature <- readRDS("EN_Batch01to08_annotated_L2_Vasculature.rds")

MD_L2_cell_id <- colnames(EN_Batch01to08_Vasculature)
MD_L2_cell_type <- EN_Batch01to08_Vasculature$L2_Cell_type
MD_L2_cluster_number <- EN_Batch01to08_Vasculature$L2_RNA_snn_res.0.5
names(MD_L2_cluster_number) <- MD_L2_cell_id

RD_L2_harmony_umap <- EN_Batch01to08_Vasculature@reductions$harmony_umap_L2


EN_Batch01to08_final$L2_Cell_type <- EN_Batch01to08_final$L1_Cell_type
EN_Batch01to08_final@meta.data[MD_L2_cell_id, "L2_Cell_type"] <- MD_L2_cell_type


EN_Batch01to08_final$L2_RNA_snn_res.0.5 <- NA
EN_Batch01to08_final@meta.data[MD_L2_cell_id, "L2_RNA_snn_res.0.5"] <- MD_L2_cluster_number


rm(EN_Batch01to08_Vasculature)
gc()


############### L3_ECs
EN_Batch01to08_ECs <- readRDS("EN_Batch01to08_5_Vasculature_ECs_reanno.rds")



EN_Batch01to08_ECs$L3ECs_Subtypes <- EN_Batch01to08_ECs$L3_ECs_Subtypes
EN_Batch01to08_ECs@meta.data[EN_Batch01to08_ECs$L3_ECs_Subtypes=="debri", "L3ECs_Subtypes"] <- "Debri"
EN_Batch01to08_ECs@meta.data[EN_Batch01to08_ECs$L3_ECs_Subtypes=="doublet", "L3ECs_Subtypes"] <- "Doublet"
EN_Batch01to08_ECs@meta.data[EN_Batch01to08_ECs$L3_ECs_Subtypes=="cEC_4", "L3ECs_Subtypes"] <- "corEC_late"
EN_Batch01to08_ECs@meta.data[EN_Batch01to08_ECs$L3_ECs_Subtypes=="precursor", "L3ECs_Subtypes"] <- "SV"
EN_Batch01to08_ECs@meta.data[EN_Batch01to08_ECs$L3_ECs_Subtypes=="endo_2", "L3ECs_Subtypes"] <- "endo_cushion"
EN_Batch01to08_ECs@meta.data[EN_Batch01to08_ECs$L3_ECs_Subtypes=="endo_1", "L3ECs_Subtypes"] <- "endo_EndoMT"



EN_Batch01to08_ECs$L3ECs_Cell_types <- EN_Batch01to08_ECs$L3_ECs_Cell_types
EN_Batch01to08_ECs@meta.data[EN_Batch01to08_ECs$L3_ECs_Cell_types=="Arterial", "L3ECs_Cell_types"] <- "ECs_Arterial"
EN_Batch01to08_ECs@meta.data[EN_Batch01to08_ECs$L3_ECs_Cell_types=="Capillary", "L3ECs_Cell_types"] <- "ECs_Cpillary"
EN_Batch01to08_ECs@meta.data[EN_Batch01to08_ECs$L3_ECs_Cell_types=="debri", "L3ECs_Cell_types"] <- "Debri"
EN_Batch01to08_ECs@meta.data[EN_Batch01to08_ECs$L3_ECs_Cell_types=="doublet", "L3ECs_Cell_types"] <- "Doublet"
EN_Batch01to08_ECs@meta.data[EN_Batch01to08_ECs$L3_ECs_Cell_types=="Endocardial", "L3ECs_Cell_types"] <- "Endocardial-like"
EN_Batch01to08_ECs@meta.data[EN_Batch01to08_ECs$L3_ECs_Cell_types=="Lymphatic", "L3ECs_Cell_types"] <- "ECs_Lymphatic"
EN_Batch01to08_ECs@meta.data[EN_Batch01to08_ECs$L3_ECs_Cell_types=="Mitotic", "L3ECs_Cell_types"] <- "ECs_Mitotic"
EN_Batch01to08_ECs@meta.data[EN_Batch01to08_ECs$L3_ECs_Cell_types=="Precursors", "L3ECs_Cell_types"] <- "Sinus Venosus"
EN_Batch01to08_ECs@meta.data[EN_Batch01to08_ECs$L3_ECs_Cell_types=="Venous", "L3ECs_Cell_types"] <- "ECs_Venous"




MD_L3_ECs_cell_id <- colnames(EN_Batch01to08_ECs)
MD_L3_ECs_cell_type <- EN_Batch01to08_ECs$L3ECs_Cell_types
MD_L3_ECs_cell_subtype <- EN_Batch01to08_ECs$L3ECs_Subtypes
MD_L3_ECs_cluster_number <- EN_Batch01to08_ECs$L3_RNA_snn_res.0.7
names(MD_L3_ECs_cluster_number) <- MD_L3_ECs_cell_id

RD_L3_ECs_harmony_umap <- EN_Batch01to08_ECs@reductions$harmony_umap_L3


EN_Batch01to08_final$L3ECs_Cell_type <- EN_Batch01to08_final$L2_Cell_type
EN_Batch01to08_final@meta.data[MD_L3_ECs_cell_id, "L3ECs_Cell_type"] <- MD_L3_ECs_cell_type

EN_Batch01to08_final$L3ECs_Subtype <- EN_Batch01to08_final$L3ECs_Cell_type
EN_Batch01to08_final@meta.data[MD_L3_ECs_cell_id, "L3ECs_Subtype"] <- MD_L3_ECs_cell_subtype


EN_Batch01to08_final$L3ECs_RNA_snn_res.0.7 <- NA
EN_Batch01to08_final@meta.data[MD_L3_ECs_cell_id, "L3ECs_RNA_snn_res.0.7"] <- MD_L3_ECs_cluster_number


rm(EN_Batch01to08_ECs)
gc()




############### L3 Mural
EN_Batch01to08_Mural<-readRDS("EN_Batch01to08_5_Vasculature_Mural_pure.rds")


EN_Batch01to08_Mural$L3Mural_Subtypes <- EN_Batch01to08_Mural$L3_Mural_Subtypes
EN_Batch01to08_Mural@meta.data[EN_Batch01to08_Mural$L3_Mural_Subtypes=="Pre_1", "L3Mural_Subtypes"] <- "pre_transient"
EN_Batch01to08_Mural@meta.data[EN_Batch01to08_Mural$L3_Mural_Subtypes=="Pre_2", "L3Mural_Subtypes"] <- "pre_EndMT"
EN_Batch01to08_Mural@meta.data[EN_Batch01to08_Mural$L3_Mural_Subtypes=="Nervous", "L3Mural_Subtypes"] <- "pre_supporting"
EN_Batch01to08_Mural@meta.data[EN_Batch01to08_Mural$L3_Mural_Subtypes=="Mito_1", "L3Mural_Subtypes"] <- "Mural_mito_1"
EN_Batch01to08_Mural@meta.data[EN_Batch01to08_Mural$L3_Mural_Subtypes=="Mito_2", "L3Mural_Subtypes"] <- "Mural_mito_2"
EN_Batch01to08_Mural@meta.data[EN_Batch01to08_Mural$L3_Mural_Subtypes=="Mito_3", "L3Mural_Subtypes"] <- "Mural_mito_3"



EN_Batch01to08_Mural$L3Mural_Cell_types <- EN_Batch01to08_Mural$L3_Mural_Cell_types
EN_Batch01to08_Mural@meta.data[EN_Batch01to08_Mural$L3_Mural_Cell_types=="Fibs", "L3Mural_Cell_types"] <- "Mural_Fibroblasts"
EN_Batch01to08_Mural@meta.data[EN_Batch01to08_Mural$L3_Mural_Cell_types=="Nervous", "L3Mural_Cell_types"] <- "Mural_Precursors"
EN_Batch01to08_Mural@meta.data[EN_Batch01to08_Mural$L3_Mural_Cell_types=="Precursors", "L3Mural_Cell_types"] <- "Mural_Precursors"
EN_Batch01to08_Mural@meta.data[EN_Batch01to08_Mural$L3_Mural_Cell_types=="Mitotic", "L3Mural_Cell_types"] <- "Mural_Mitotic"
EN_Batch01to08_Mural@meta.data[EN_Batch01to08_Mural$L3_Mural_Cell_types=="PCs", "L3Mural_Cell_types"] <- "Mural_Pericytes"
EN_Batch01to08_Mural@meta.data[EN_Batch01to08_Mural$L3_Mural_Cell_types=="SMCs", "L3Mural_Cell_types"] <- "Mural_SMCs"





MD_L3_Mural_cell_id <- colnames(EN_Batch01to08_Mural)
MD_L3_Mural_cell_type <- EN_Batch01to08_Mural$L3Mural_Cell_types
MD_L3_Mural_cell_subtype <- EN_Batch01to08_Mural$L3Mural_Subtypes
MD_L3_Mural_cluster_number <- EN_Batch01to08_Mural$L3_CCA_snn_res.0.6
names(MD_L3_Mural_cluster_number) <- MD_L3_Mural_cell_id

RD_L3_Mural_CCA_umap <- EN_Batch01to08_Mural@reductions$umap_pure


EN_Batch01to08_final$L3Mural_Cell_type <- EN_Batch01to08_final$L2_Cell_type
EN_Batch01to08_final@meta.data[MD_L3_Mural_cell_id, "L3Mural_Cell_type"] <- MD_L3_Mural_cell_type

EN_Batch01to08_final@meta.data[EN_Batch01to08_final$L3Mural_Cell_type=="Mural", "L3Mural_Cell_type"] <- "Debri"

EN_Batch01to08_final$L3Mural_Subtype <- EN_Batch01to08_final$L3Mural_Cell_type
EN_Batch01to08_final@meta.data[MD_L3_Mural_cell_id, "L3Mural_Subtype"] <- MD_L3_Mural_cell_subtype

EN_Batch01to08_final@meta.data[EN_Batch01to08_final$L3Mural_Subtype=="Mural", "L3Mural_Subtype"] <- "Debri"

EN_Batch01to08_final$L3Mural_CCA_snn_res.0.6 <- NA
EN_Batch01to08_final@meta.data[MD_L3_Mural_cell_id, "L3Mural_CCA_snn_res.0.6"] <- MD_L3_Mural_cluster_number


rm(EN_Batch01to08_Mural)
gc()








############### L3 combine
EN_Batch01to08_final$L3_Cell_type <- EN_Batch01to08_final$L2_Cell_type
EN_Batch01to08_final@meta.data[MD_L3_ECs_cell_id, "L3_Cell_type"] <- MD_L3_ECs_cell_type
EN_Batch01to08_final@meta.data[MD_L3_Mural_cell_id, "L3_Cell_type"] <- MD_L3_Mural_cell_type
EN_Batch01to08_final@meta.data[EN_Batch01to08_final$L3_Cell_type=="Mural", "L3_Cell_type"] <- "Debri"

EN_Batch01to08_final$L3_Subtype <- EN_Batch01to08_final$L2_Cell_type
EN_Batch01to08_final@meta.data[MD_L3_ECs_cell_id, "L3_Subtype"] <- MD_L3_ECs_cell_subtype
EN_Batch01to08_final@meta.data[MD_L3_Mural_cell_id, "L3_Subtype"] <- MD_L3_Mural_cell_subtype
EN_Batch01to08_final@meta.data[EN_Batch01to08_final$L3_Subtype=="Mural", "L3_Subtype"] <- "Debri"




saveRDS(EN_Batch01to08_final, "EN_Batch01to08_6_All.rds")





############################## scRNA-seq Update annotation
EN_all <- readRDS("EN_Batch01to08_6_All.rds")

EN_all$newL3_Subtype <- EN_all$L3_Subtype
EN_all$newL3_Subtype[EN_all$L3_Subtype=="endo_EndoMT"] <- "CMs"
EN_all$newL3_Subtype[EN_all$L3_Subtype=="SV"] <- "SV_related"
EN_all$newL3_Subtype[EN_all$L3_Subtype=="endo_cushion"] <- "endoMT"


EN_all$newL3_Cell_type <- EN_all$L3_Cell_type
EN_all$newL3_Cell_type[EN_all$L3_Subtype=="endo_EndoMT"] <- "CMs"
EN_all$newL3_Cell_type[EN_all$L3_Subtype=="SV"] <- "SinusVenosus_related"
EN_all$newL3_Cell_type[EN_all$L3_Cell_type=="ECs_Cpillary"] <- "ECs_Capillary"

EN_all$newL3_Cell_type[EN_all$L3_Cell_type == "Mural_Pericytes"] <- "Pericytes"
EN_all$newL3_Cell_type[EN_all$L3_Cell_type == "Mural_SMCs"] <- "SMCs"
EN_all$newL3_Cell_type[EN_all$L3_Cell_type == "Mural_Fibroblasts"] <- "Fibroblasts"


EN_all$newL2_Cell_type <- EN_all$L2_Cell_type
EN_all$newL2_Cell_type[EN_all$L3_Subtype=="endo_EndoMT"] <- "CMs"


EN_all$newL1_Cell_type <- EN_all$L1_Cell_type
EN_all$newL1_Cell_type[EN_all$L3_Subtype=="endo_EndoMT"] <- "CMs"






################### Subset vasculature
CT_L3_Cell_type <- unique(EN_all$newL3_Cell_type)
CT_ECs_Cell_type <- c(CT_L3_Cell_type[grepl("ECs", CT_L3_Cell_type)], "Endocardial-like", "SinusVenosus_related", "Tip cells")
CT_Mural_Cell_type <- c(CT_L3_Cell_type[grepl("Mural", CT_L3_Cell_type)], "Fibroblasts", "Pericytes", "SMCs")
CT_Vasculature_Cell_type <- c(CT_ECs_Cell_type, CT_Mural_Cell_type)


EN_vasculature <- subset(EN_all, newL3_Cell_type %in% CT_Vasculature_Cell_type)


EN_vasculature <- RunUMAP(EN_vasculature, reduction = "harmony", dims = 1:40, n.neighbors = 20, reduction.name = "harmony_umap_vasculature", reduction.key = "harmonyUmapVasculature_")






################### Further Subset

EN_brain_ECs <- subset(EN_vasculature, Organ == "Brain" & newL3_Cell_type %in% setdiff(CT_ECs_Cell_type, c("SinusVenosus_related", "ECs_Lymphatic")) & newL3_Subtype != "corEC_late")
saveRDS(EN_brain_ECs, "EN_Batch01to08_7_Brain_ECs.rds")


EN_heart_ECs <- subset(EN_vasculature, Organ == "Heart" & newL3_Cell_type %in% setdiff(CT_ECs_Cell_type, "ECs_Lymphatic") & newL3_Subtype != "vEC_2")
saveRDS(EN_heart_ECs, "EN_Batch01to08_7_Heart_ECs.rds")


EN_brain_Mural <- subset(EN_vasculature, Organ == "Brain" & newL3_Cell_type %in% CT_Mural_Cell_type)
saveRDS(EN_brain_Mural, "EN_Batch01to08_7_Brain_Mural.rds")


EN_heart_Mural <- subset(EN_vasculature, Organ == "Heart" & newL3_Cell_type %in% CT_Mural_Cell_type)
saveRDS(EN_heart_Mural, "EN_Batch01to08_7_Heart_Mural.rds")

