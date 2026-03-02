library(Seurat)
library(clustree)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(SeuratDisk)
library(slingshot)
library(tradeSeq)
library(clusterProfiler)
library(ComplexHeatmap)


################### Load data
EN_brain_ECs_pure <- readRDS("EN_Batch01to08_7_Brain_ECs.rds")

pal <- c(RColorBrewer::brewer.pal(9,"Set1"),RColorBrewer::brewer.pal(9,"Pastel1"),RColorBrewer::brewer.pal(8,"Accent"),RColorBrewer::brewer.pal(8,"Set2"),RColorBrewer::brewer.pal(8,"Pastel2") ,scales::hue_pal()(8))

EN_brain_ECs_nonGW24 <- subset(EN_brain_ECs_pure, Age != "GW24")




################### retain only venous, capillary, and artery ECs
EN_brain_ECs_nonGW24_onlyVCA <- subset(EN_brain_ECs_nonGW24, newL3_Subtype %in% c("aEC_1", "cEC_3", "cEC_1", "vEC_1"))

EN_brain_ECs_nonGW24_onlyVCA <- RunUMAP(EN_brain_ECs_nonGW24_onlyVCA, dims = 1:40, #50
                                        reduction = "harmony",
                                        reduction.name = "umap_Trajectory", reduction.key = "umapTrajectory_",
                                        min.dist = 0.001,
                                        metric = "euclidean",
                                        n.neighbors = 10)


### Clustering
EN_brain_ECs_nonGW24_onlyVCA <- FindNeighbors(EN_brain_ECs_nonGW24_onlyVCA, dims = 1:40, k.param = 20, reduction = "harmony")

EN_brain_ECs_nonGW24_onlyVCA <- FindClusters(EN_brain_ECs_nonGW24_onlyVCA, graph.name = "RNA_snn", resolution = 0.7, algorithm = 1)

colnames(EN_brain_ECs_nonGW24_onlyVCA@meta.data)[which(colnames(EN_brain_ECs_nonGW24_onlyVCA@meta.data)=="RNA_snn_res.0.7")] <- "Trajectory_RNA_snn_res.0.7"

clusters_all <- names(table(EN_brain_ECs_nonGW24_onlyVCA$Trajectory_RNA_snn_res.0.7))





################### Subset cluster
clusters_to_remove <- c(10, 8, 7)
EN_brain_ECs_selected <- subset(EN_brain_ECs_nonGW24_onlyVCA, Trajectory_RNA_snn_res.0.7 %in% setdiff(clusters_all, clusters_to_remove))


EN_brain_ECs_selected$Trajectory_RNA_snn_res.0.7 <- droplevels(EN_brain_ECs_selected$Trajectory_RNA_snn_res.0.7)


EN_brain_ECs_selected <- RunUMAP(EN_brain_ECs_selected, dims = 1:40, #50
                               reduction = "harmony",
                               reduction.name = "umap_Trajectory", reduction.key = "umapTrajectory_",
                               min.dist = 0.001,
                               metric = "euclidean",
                               n.neighbors = 10)


################### Run Slingshot
Brain_ECs_sce <- as.SingleCellExperiment(EN_brain_ECs_selected)
reducedDims(Brain_ECs_sce) <- list(umap_Trajectory = EN_brain_ECs_selected@reductions$umap_Trajectory@cell.embeddings)


Brain_ECs_clustering <- EN_brain_ECs_selected$Trajectory_RNA_snn_res.0.7
pal_clustering <- pal[seq(length(table(EN_brain_ECs_selected$Trajectory_RNA_snn_res.0.7)))]
names(pal_clustering)<-sort(unique(Brain_ECs_clustering))




Brain_ECs_sce_sW_C9 <- slingshot(Brain_ECs_sce, clusterLabels = 'Trajectory_RNA_snn_res.0.7', reducedDim = 'umap_Trajectory', start.clus = "9",
                                       extend = "n", stretch=1, thresh=0.1)
Slingshot_sW_C9<-SlingshotDataSet(Brain_ECs_sce_sW_C9)





### Convert for velocity analysis
EN_brain_ECs_selected_diet<-DietSeurat(EN_brain_ECs_selected, scale.data = TRUE, dimreducs = c("umap_Trajectory"))

SaveH5Seurat(EN_brain_ECs_selected_diet, filename = "EN_brain_ECs_nonGW24_onlyVCA_selected_trajectory.h5Seurat")

Convert("EN_brain_ECs_nonGW24_onlyVCA_selected_trajectory.h5Seurat", dest = "h5ad")




pdf("Trajecotry_Brain_ECs_nonGW24_onlyVCA_selected_UMAP_clustering_sW_C9_lineages_curves.pdf", width=6, height=6)
plot(reducedDims(Brain_ECs_sce_sW_C9)$umap_Trajectory, col = pal_clustering[Brain_ECs_clustering], pch=16, asp = 1, cex=0.5)
lines(SlingshotDataSet(Brain_ECs_sce_sW_C9), lwd=3)
dev.off()






################### TradeSeq for C9 trajectory

## TradeSeq analysis
Comb1_counts <- as.matrix(EN_brain_ECs_selected@assays$RNA@counts[EN_brain_ECs_selected@assays$RNA@var.features, ])

Comb1_pseudotime <- slingPseudotime(Slingshot_sW_C9, na = FALSE)
Comb1_cellWeights <- slingCurveWeights(Slingshot_sW_C9)

Comb1_samples <- sample( colnames(EN_brain_ECs_selected) , size = 1000)
Comb1_counts_sampling <- Comb1_counts[rowSums(Comb1_counts > 3) > ncol(Comb1_counts)/100, Comb1_samples]



set.seed(6)
Comb1_icMat_seed6 <- evaluateK(counts = Comb1_counts_sampling, k=3:10, pseudotime = Comb1_pseudotime[Comb1_samples, ], cellWeights=Comb1_cellWeights[Comb1_samples, ], verbose = T, parallel=T)
plot_evalutateK_results(Comb1_icMat_seed6)

set.seed(4)
Comb1_icMat_seed4 <- evaluateK(counts = Comb1_counts_sampling, k=3:10, pseudotime = Comb1_pseudotime[Comb1_samples, ], cellWeights=Comb1_cellWeights[Comb1_samples, ], verbose = T, parallel=T)
plot_evalutateK_results(Comb1_icMat_seed4)


## fit GAM
TradeSeq_Comb1 <- fitGAM(counts = Comb1_counts, pseudotime = Comb1_pseudotime, cellWeights=Comb1_cellWeights, 
                                    nknots=7, verbose = T, parallel=T)
# Note! Modify nknots





## pattern test
Comb1_patternRes <- patternTest(TradeSeq_Comb1)
Comb1_patternRes_filter <- filter_TradeSeq_results(Comb1_patternRes)

Comb1_patternRes_top200 <- rownames(Comb1_patternRes_filter)[1:min(200,nrow(Comb1_patternRes_filter))]
# Check the top 200 genes due to three lineages

Comb1_smoothed <- generate_lineage_smoothed_matrix(Slingshot_sW_C9, EN_brain_ECs_selected, Comb1_patternRes_top200)
Comb1_smoothed_sorted <- plot_lineage_heatmap(Comb1_smoothed, Slingshot_sW_C9, order="hclust", outPrefix = "brain_ECs_nonGW24_onlyVCA_selected_top200", height=22)


## Identify lineage associated genes
Comb1_genes_L1_candidate <- Comb1_smoothed_sorted[c(1:56, 190:196)]
Comb1_genes_L2_candidate <- Comb1_smoothed_sorted[c(9:20, 57:113)]
Comb1_genes_L3_candidate <- Comb1_smoothed_sorted[c(114:200)]

for (g in Comb1_genes_L1_candidate) {
  plot_gene_expression(Slingshot_sW_C9, EN_brain_ECs_selected, g, outPrefix = "brain_ECs_nonGW24_onlyVCA_selected_vECs")
}

for (g in Comb1_genes_L2_candidate) {
  plot_gene_expression(Slingshot_sW_C9, EN_brain_ECs_selected, g, outPrefix = "brain_ECs_nonGW24_onlyVCA_selected_aECs")
}

for (g in Comb1_genes_L3_candidate) {
  plot_gene_expression(Slingshot_sW_C9, EN_brain_ECs_selected, g, outPrefix = "brain_ECs_nonGW24_onlyVCA_selected_cECs")
}



Comb1_genes_L1_selected <- c("SELE", "ADGRG6", "CCL14", "CEBPD", "FNIP2", "ICAM1", "MMRN1", "PLAT", "PLVAP", "S100A10", "SEMA3A", "TMEM176B", "VWF")
Comb1_genes_L2_selected <- c("GJA4", "ARL15", "UNC5B", "C12orf75", "HEY1", "CD9", "DIAPH2", "MCTP1", "DLL4", "DST", "KCTD12", "EGFL8", "EYS", "FILIP1", "JAG2", "PIK3R3", "SEMA3G", "SLC45A4", "STOM")
Comb1_genes_L3_selected <- c("ADGRL3", "APCDD1", "BSG", "CD320", "CDH6", "CETP", "CHGA", "CSRP2", "DAB2", "GMFG", "HOPX", "HPGD", "PRDX2", "RGCC", "SLC7A8", "SLC16A1", "SLC38A5")

dir.create("nonGW24_onlyVCA_selected_TradeSeq/Selected", recursive=TRUE)
sapply(Comb1_genes_L1_selected, function(gene) {
  src <- file.path("nonGW24_onlyVCA_selected_TradeSeq", paste0("Gene_brain_ECs_nonGW24_onlyVCA_selected_vECs_", gene, ".png"))
  if (file.exists(src)) file.copy(src, "nonGW24_onlyVCA_selected_TradeSeq/Selected") else warning(paste("Not exists:", src))
})

sapply(Comb1_genes_L2_selected, function(gene) {
  src <- file.path("nonGW24_onlyVCA_selected_TradeSeq", paste0("Gene_brain_ECs_nonGW24_onlyVCA_selected_aECs_", gene, ".png"))
  if (file.exists(src)) file.copy(src, "nonGW24_onlyVCA_selected_TradeSeq/Selected") else warning(paste("Not exists:", src))
})

sapply(Comb1_genes_L3_selected, function(gene) {
  src <- file.path("nonGW24_onlyVCA_selected_TradeSeq", paste0("Gene_brain_ECs_nonGW24_onlyVCA_selected_cECs_", gene, ".png"))
  if (file.exists(src)) file.copy(src, "nonGW24_onlyVCA_selected_TradeSeq/Selected") else warning(paste("Not exists:", src))
})


## select genes
DotPlot(EN_brain_ECs_selected, features = Comb1_genes_L1_selected, group.by = "Age") +
  coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust=1))
ggsave("DotPlot_brain_ECs_nonGW24_onlyVCA_selected_vECs_selected.pdf", width=6, height=6)

DotPlot(EN_brain_ECs_selected, features = Comb1_genes_L2_selected, group.by = "Age") +
  coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust=1))
ggsave("DotPlot_brain_ECs_nonGW24_onlyVCA_selected_aECs_selected.pdf", width=6, height=6)

DotPlot(EN_brain_ECs_selected, features = Comb1_genes_L3_selected, group.by = "Age") +
  coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust=1))
ggsave("DotPlot_brain_ECs_nonGW24_onlyVCA_selected_cECs_selected.pdf", width=6, height=6)


Comb1_selected_smoothed <- generate_lineage_smoothed_matrix(Slingshot_sW_C9, EN_brain_ECs_selected, c(Comb1_genes_L1_selected, Comb1_genes_L2_selected, Comb1_genes_L3_selected))
Comb1_selected_smoothed_sorted <- plot_lineage_heatmap(Comb1_selected_smoothed, Slingshot_sW_C9, order="max", outPrefix = "brain_ECs_nonGW24_onlyVCA_selected_selected", height = 6)


dir.create("nonGW24_onlyVCA_selected_final/Selected", recursive=TRUE)
for (gene in Comb1_genes_L1_selected) {
  FeaturePlot(EN_brain_ECs_selected, gene, reduction = "umap_Trajectory")
  ggsave(paste0("nonGW24_onlyVCA_selected_final/Selected/FeaturePlot_brain_ECs_nonGW24_onlyVCA_selected_vECs_", gene, ".png"), width=5, height=4)
}
for (gene in Comb1_genes_L2_selected) {
  FeaturePlot(EN_brain_ECs_selected, gene, reduction = "umap_Trajectory")
  ggsave(paste0("nonGW24_onlyVCA_selected_final/Selected/FeaturePlot_brain_ECs_nonGW24_onlyVCA_selected_aECs_", gene, ".png"), width=5, height=4)
}
for (gene in Comb1_genes_L3_selected) {
  FeaturePlot(EN_brain_ECs_selected, gene, reduction = "umap_Trajectory")
  ggsave(paste0("nonGW24_onlyVCA_selected_final/Selected/FeaturePlot_brain_ECs_nonGW24_onlyVCA_selected_cECs_", gene, ".png"), width=5, height=4)
}






############################# Function #############################

filter_TradeSeq_results<-function(TradeSeq_results, filter=TRUE, thr_FDR=0.05) {
  TradeSeq_results$pvalue[is.na(TradeSeq_results$pvalue)] <- 1
  TradeSeq_results <- TradeSeq_results[order(TradeSeq_results$waldStat,decreasing = T),]
  TradeSeq_results$FDR <- p.adjust(TradeSeq_results$pvalue,method = "BH")
  if (filter) {
    TradeSeq_results <- TradeSeq_results[ TradeSeq_results$FDR < thr_FDR , ]
  }
  return(TradeSeq_results)
}

generate_lineage_smoothed_matrix <- function(Slingshot_object, Seurat_object, gene_list, nknots = 20, spar = 0.8, weight_threshold = 1) {
  pseudotime <- slingPseudotime(Slingshot_object, na = FALSE)
  cellWeights <- slingCurveWeights(Slingshot_object)
  lims <- quantile(pseudotime,c(0.02, 0.98), na.rm = TRUE)
  gene_expression <- Seurat_object@assays$RNA@data[gene_list,]
  n_lineages <- ncol(pseudotime)
  
  build_column_structure <- function(n_lineages, nknots) {
    if (n_lineages == 2) {
      c(rev(paste0("Lineage1_", seq_len(nknots))), 
        "Separator",
        paste0("Lineage2_", seq_len(nknots)))
    } else {
      unlist(lapply(seq_len(n_lineages), function(l) {
        c(paste0("Lineage", l, "_", seq_len(nknots)), 
          if (l < n_lineages) "Separator")
      }))
    }
  }
  col_structure <- build_column_structure(n_lineages, nknots)
  
  res <- matrix(
    nrow = length(gene_list),
    ncol = length(col_structure),
    dimnames = list(gene_list, col_structure)
  )
  
  
  for (g in seq_along(gene_list)) {
    gene <- gene_list[g]
    lineage_data <- vector("list", n_lineages)
    
    for (l in seq_len(n_lineages)) {
      valid_idx <- which(
        cellWeights[, l] >= weight_threshold &
          pseudotime[, l] > lims[1] &
          pseudotime[, l] < lims[2]
      )
      
      if (length(valid_idx) < 4) {
        lineage_data[[l]] <- rep(NA, nknots)
        next
      }
      
      tryCatch({
        fit <- smooth.spline(
          x = pseudotime[valid_idx, l],
          y = gene_expression[gene, valid_idx],
          nknots = nknots,
          spar = spar
        )
        sm <- spline(fit, n = nknots)
        lineage_data[[l]] <- sm$y
      }, error = function(e) {
        message("Smoothing failed for ", gene, " in Lineage", l)
        lineage_data[[l]] <- rep(NA, nknots)
      })
    }
    
    if (n_lineages == 2) {
      merged <- c(rev(lineage_data[[1]]), NA, lineage_data[[2]])
    } else {
      merged <- unlist(lapply(seq_len(n_lineages), function(l) {
        c(lineage_data[[l]], if (l < n_lineages) NA)
      }))
    }
    
    res[g, ] <- merged
  }
  
  return(res)
}

plot_lineage_heatmap <- function(lineage_sm_results, Slingshot_object, order="hclust", outPrefix=NULL, height=12) {
  sm_row_scale <- t(apply(lineage_sm_results,1,function(x){scale(x,T,T)}))
  
  separator_pos <- which(colnames(lineage_sm_results) == "Separator")
  n_lineages <- length(separator_pos) + 1
  starts <- c(1, separator_pos + 1)
  ends <- c(separator_pos - 1, ncol(lineage_sm_results))
  
  lineage_labels <- Slingshot_object@lineages[1:n_lineages]
  
  annot_labels <- rep("", ncol(lineage_sm_results))
  for (l in seq_len(n_lineages)) {
    lbls <- lineage_labels[[l]]
    
    seq_range <- if (n_lineages == 2) {
      if (l == 1) {
        seq(ends[l]+1, starts[l], length.out = length(lbls))
      } else {
        seq(starts[l]-1, ends[l], length.out = length(lbls))
      }
    } else {
      seq(starts[l], ends[l], length.out = length(lbls))
    }
    
    annot_labels[round(seq_range)] <- lbls
  }
  
  top_anno <- HeatmapAnnotation(
    NodeLabels = anno_text(
      annot_labels,
      rot = 0,
      location = unit(2, "mm"),
      just = "center",
      gp = gpar(fontsize = 10)
    ),
    annotation_height = unit(0.5, "cm")
  )
  
  
  if (! is.null(outPrefix)) {
    pdf(paste0("Heatmap_", outPrefix, ".pdf"), width = 4 + n_lineages*2, height=height, useDingbats = F)
  }
  
  if (order=="hclust") {
    hclust_res <- hclust(as.dist((1 - cor(t(lineage_sm_results), use = "pairwise.complete.obs"))/2), method = "complete")
    gene_order <- hclust_res$order
    draw(Heatmap(
      matrix = sm_row_scale,
      name = "Z-score",
      na_col = "white",
      col = colorRampPalette(c("grey95", "grey70", "firebrick", "firebrick4"))(90),
      cluster_rows = as.dendrogram(hclust_res),
      cluster_columns = FALSE,
      top_annotation = top_anno,
      show_row_names = TRUE,
      row_names_side = "right",
      row_names_gp = gpar(fontsize = 10)
    ))
  } else if (order=="max") {
    sm_row_max <- apply(sm_row_scale,1,function(x){ which.max(x) })
    gene_order <- order(sm_row_max)
    draw(Heatmap(
      matrix = sm_row_scale[gene_order,],
      name = "Z-score",
      na_col = "white",
      col = colorRampPalette(c("grey95", "grey70", "firebrick", "firebrick4"))(90),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      top_annotation = top_anno,
      show_row_names = TRUE,
      row_names_side = "right",
      row_names_gp = gpar(fontsize = 10)
    ))
  }
  
  if (! is.null(outPrefix)) {
    dev.off()
  }
  return(rownames(lineage_sm_results)[gene_order])
}

plot_gene_expression <- function(Slingshot_object, Seurat_object, gene, weight_threshold = 1, legend=TRUE, lineage_legend=NULL, outPrefix=NULL, outFormat=c("png", "pdf")) {
  pseudotime <- slingPseudotime(Slingshot_object, na = FALSE)
  cellWeights <- slingCurveWeights(Slingshot_object)
  lims <- quantile(pseudotime,c(0.02, 0.98), na.rm = TRUE)
  gene_expression <- Seurat_object@assays$RNA@data[gene, rownames(pseudotime)]
  
  n_lineages <- ncol(pseudotime)
  
  colors <- viridis::viridis(n_lineages)
  
  
  if (! is.null(outPrefix)) {
    if (outFormat=="png") {
      png(paste0("Gene_", outPrefix, "_", gene, ".png"), width = 2000, height = 1600, res = 300)
    } else if (outFormat=="pdf") {
      pdf(paste0("Gene_", outPrefix, "_", gene, ".pdf"), width = 5, height = 4)
    }
  }

  par(mar = c(5, 5, 4, 2) + 0.1)
  plot(0, type = "n", 
       xlim = lims, ylim = range(gene_expression),
       main = gene, cex.main = 2,
       xlab = "Pseudotime", ylab = "Expression Level", 
       cex.lab = 1.5, cex.axis = 1.2, las=1)
  
  for(i in seq_len(n_lineages)){
    valid_cells <- (cellWeights[,i] >= weight_threshold) & (pseudotime[,i] > lims[1]) & (pseudotime[,i] < lims[2])
    cell_ids <- colnames(Seurat_object)[valid_cells]
    points(pseudotime[cell_ids, i], gene_expression[cell_ids]+runif(n = length(cell_ids),0.0,0.02),
           cex=.5, pch=16, col=adjustcolor(colors[i], alpha.f = 0.3))
    
    tryCatch({
      spline_fit <- smooth.spline(pseudotime[cell_ids, i], 
                                  gene_expression[cell_ids],
                                  nknots = 20, spar = 0.8)
      
      lines(spline_fit$x, spline_fit$y, lwd = 6, col = colors[i])
      lines(spline_fit$x, spline_fit$y, lwd = 1, col = "black")
    }, error = function(e) message("Skipped smoothing for lineage ", i))
  }
  
  if (legend) {
    if(is.null(lineage_legend)) {
      lineage_legend <- colnames(pseudotime)
    }
    legend("topright", 
           legend = lineage_legend,
           col = colors, lwd = 4,
           cex = 1.2, bg = "white")
  }
  
  if (! is.null(outPrefix)) {
    dev.off()
  }
}
