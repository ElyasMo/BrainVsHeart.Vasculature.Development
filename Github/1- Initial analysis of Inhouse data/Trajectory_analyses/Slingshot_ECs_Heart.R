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
EN_heart_ECs_pure <- readRDS("EN_Batch01to08_7_Heart_ECs.rds")

pal <- c(RColorBrewer::brewer.pal(9,"Set1"),RColorBrewer::brewer.pal(9,"Pastel1"),RColorBrewer::brewer.pal(8,"Accent"),RColorBrewer::brewer.pal(8,"Set2"),RColorBrewer::brewer.pal(8,"Pastel2") ,scales::hue_pal()(8))

EN_heart_ECs_nonGW24 <- subset(EN_heart_ECs_pure, Age != "GW24" & newL3_Subtype != "corEC_late")



################### Arteriovenous specialization
EN_heart_ECs_nonGW24_AVS <- subset(EN_heart_ECs_nonGW24, newL3_Subtype %in% c("aEC_1", "cEC_2", "vEC_1", "tip_1"))

EN_heart_ECs_nonGW24_AVS <- RunUMAP(EN_heart_ECs_nonGW24_AVS, dims = 1:40, #50
                                        reduction = "harmony",
                                        reduction.name = "umap_Trajectory", reduction.key = "umapTrajectory_",
                                        min.dist = 0.001,
                                        metric = "euclidean",
                                        n.neighbors = 10)




### Clustering
EN_heart_ECs_nonGW24_AVS <- FindNeighbors(EN_heart_ECs_nonGW24_AVS, dims = 1:40, k.param = 20, reduction = "harmony")

EN_heart_ECs_nonGW24_AVS <- FindClusters(EN_heart_ECs_nonGW24_AVS, graph.name = "RNA_snn", resolution = 1, algorithm = 1)

colnames(EN_heart_ECs_nonGW24_AVS@meta.data)[which(colnames(EN_heart_ECs_nonGW24_AVS@meta.data)=="RNA_snn_res.1")] <- "Trajectory_RNA_snn_res.1"

clusters_AVS_all <- names(table(EN_heart_ECs_nonGW24_AVS$Trajectory_RNA_snn_res.1))




### Subset cluster
clusters_AVS_to_remove <- c(12, 13, 10, 7)
EN_heart_ECs_nonGW24_AVS_selected <- subset(EN_heart_ECs_nonGW24_AVS, Trajectory_RNA_snn_res.1 %in% setdiff(clusters_AVS_all, clusters_AVS_to_remove))

EN_heart_ECs_nonGW24_AVS_selected$Trajectory_RNA_snn_res.1 <- droplevels(EN_heart_ECs_nonGW24_AVS_selected$Trajectory_RNA_snn_res.1)


EN_heart_ECs_nonGW24_AVS_selected <- RunUMAP(EN_heart_ECs_nonGW24_AVS_selected, dims = 1:40, #50
                                 reduction = "harmony",
                                 reduction.name = "umap_Trajectory", reduction.key = "umapTrajectory_",
                                 min.dist = 0.001,
                                 metric = "euclidean",
                                 n.neighbors = 10)


### Run Slingshot
Heart_ECs_nonGW24_AVS_selected_sce <- as.SingleCellExperiment(EN_heart_ECs_nonGW24_AVS_selected)
reducedDims(Heart_ECs_nonGW24_AVS_selected_sce) <- list(umap_Trajectory = EN_heart_ECs_nonGW24_AVS_selected@reductions$umap_Trajectory@cell.embeddings)


Heart_ECs_nonGW24_AVS_selected_clustering <- EN_heart_ECs_nonGW24_AVS_selected$Trajectory_RNA_snn_res.1
pal_clustering_nonGW24_AVS_selected <- pal[seq(length(table(Heart_ECs_nonGW24_AVS_selected_clustering)))]
names(pal_clustering_nonGW24_AVS_selected)<-sort(unique(Heart_ECs_nonGW24_AVS_selected_clustering))



Heart_ECs_nonGW24_AVS_selected_sce <- slingshot(Heart_ECs_nonGW24_AVS_selected_sce, clusterLabels = 'Trajectory_RNA_snn_res.1', reducedDim = 'umap_Trajectory', start.clus='6',
                                 extend = "n", stretch=1, thresh=0.1)
Slingshot_nonGW24_AVS_selected_sW_C6<-SlingshotDataSet(Heart_ECs_nonGW24_AVS_selected_sce)




### Convert for velocity analysis
EN_heart_ECs_nonGW24_AVS_selected_diet<-DietSeurat(EN_heart_ECs_nonGW24_AVS_selected, scale.data = TRUE, dimreducs = c("umap_Trajectory"))

SaveH5Seurat(EN_heart_ECs_nonGW24_AVS_selected_diet, filename = "EN_heart_ECs_nonGW24_AVS_selected_trajectory.h5Seurat")

Convert("EN_heart_ECs_nonGW24_AVS_selected_trajectory.h5Seurat", dest = "h5ad")




pdf("Trajecotry_Heart_ECs_nonGW24_AVS_selected_UMAP_clustering_sW_C6_lineages_curves.pdf", width=6, height=6)
plot(reducedDims(Heart_ECs_nonGW24_AVS_selected_sce)$umap_Trajectory, col = pal_clustering_nonGW24_AVS_selected[Heart_ECs_nonGW24_AVS_selected_clustering], pch=16, asp = 1)
lines(SlingshotDataSet(Heart_ECs_nonGW24_AVS_selected_sce), lwd=3)
dev.off()



## TradeSeq analysis
Comb1_counts <- as.matrix(EN_heart_ECs_nonGW24_AVS_selected@assays$RNA@counts[EN_heart_ECs_nonGW24_AVS_selected@assays$RNA@var.features, ])

Comb1_pseudotime <- slingPseudotime(Slingshot_nonGW24_AVS_selected_sW_C6, na = FALSE)
Comb1_cellWeights <- slingCurveWeights(Slingshot_nonGW24_AVS_selected_sW_C6)

Comb1_samples <- sample( colnames(EN_heart_ECs_nonGW24_AVS_selected) , size = 1000)
Comb1_counts_sampling <- Comb1_counts[rowSums(Comb1_counts > 3) > ncol(Comb1_counts)/100, Comb1_samples]



set.seed(6)
Comb1_icMat_seed6 <- evaluateK(counts = Comb1_counts_sampling, k=3:10, pseudotime = Comb1_pseudotime[Comb1_samples, ], cellWeights=Comb1_cellWeights[Comb1_samples, ], verbose = T, parallel=T)
plot_evalutateK_results(Comb1_icMat_seed6)


set.seed(4)
Comb1_icMat_seed4 <- evaluateK(counts = Comb1_counts_sampling, k=3:10, pseudotime = Comb1_pseudotime[Comb1_samples, ], cellWeights=Comb1_cellWeights[Comb1_samples, ], verbose = T, parallel=T)
plot_evalutateK_results(Comb1_icMat_seed4)



## fit GAM
set.seed(6)
TradeSeq_Comb1 <- fitGAM(counts = Comb1_counts, pseudotime = Comb1_pseudotime, cellWeights=Comb1_cellWeights, 
                                    nknots=9, verbose = T, parallel=T)
# Note! Modify nknots




Slingshot_Comb1 <- Slingshot_nonGW24_AVS_selected_sW_C6
Seurat_Comb1 <- EN_heart_ECs_nonGW24_AVS_selected


## pattern test
Comb1_patternRes <- patternTest(TradeSeq_Comb1)
Comb1_patternRes_filter <- filter_TradeSeq_results(Comb1_patternRes)

Comb1_patternRes_top200 <- rownames(Comb1_patternRes_filter)[1:min(200,nrow(Comb1_patternRes_filter))]
# Check the top 200 genes due to three lineages

Comb1_smoothed <- generate_lineage_smoothed_matrix(Slingshot_Comb1, Seurat_Comb1, Comb1_patternRes_top200)
Comb1_smoothed_sorted <- plot_lineage_heatmap(Comb1_smoothed, Slingshot_Comb1, order="hclust", outPrefix = "heart_ECs_nonGW24_AVS_selected_top200", height=22)


## Identify lineage associated genes
Comb1_genes_L1_candidate <- Comb1_smoothed_sorted[c(1:63)]
Comb1_genes_L2_candidate <- Comb1_smoothed_sorted[c(64:115, 188:200)]
Comb1_genes_L3_candidate <- Comb1_smoothed_sorted[c(116:200)]

for (g in Comb1_genes_L1_candidate) {
  plot_gene_expression(Slingshot_Comb1, Seurat_Comb1, g, outPrefix = "heart_ECs_nonGW24_AVS_selected_vECs")
}

for (g in Comb1_genes_L2_candidate) {
  plot_gene_expression(Slingshot_Comb1, Seurat_Comb1, g, outPrefix = "heart_ECs_nonGW24_AVS_selected_aECs")
}

for (g in Comb1_genes_L3_candidate) {
  plot_gene_expression(Slingshot_Comb1, Seurat_Comb1, g, outPrefix = "heart_ECs_nonGW24_AVS_selected_cECs")
}


Comb1_genes_L1_selected <- c("CCL2", "C11orf96", "CEBPD", "COL3A1", "CTHRC1", "IER3", "IRF1", "KRT18", "NCOA7", "PLVAP", "POSTN", "TFPI")
Comb1_genes_L2_selected <- c("ADM", "APOA1", "ARL15", "BTNL9", "CLDN5", "CLEC14A", "CXCL12", "EGFL8", "GJA4", "GRIA2", "JAG2", "MECOM", "PDE4D", "PIK3R3", "PREX2", "PRND", "RBP7", "SERPINE2", "SOX5", "SPECC1", "SRGN", "TSPAN2", "UNC5B", "VEGFC")
Comb1_genes_L3_selected <- c("ABCC4", "ABLIM1", "ADAMTS9", "ADGRF5", "AKAP12", "ANO2", "ARHGEF28", "CACNA1C", "CD36", "CLIC5", "DACH1", "EBF1", "ESYT2", "FAM155A", "FLT1", "HECW2", "ITGA6", "KCNB2", "KITLG", "MEOX1", "NRP1", "PRKCH", "QRFPR", "SPTBN1", "TCIM", "TNNT3", "VAV3")

dir.create("nonGW24_AVS_selected_TradeSeq/Selected", recursive=TRUE)
sapply(Comb1_genes_L1_selected, function(gene) {
  src <- file.path("nonGW24_AVS_selected_TradeSeq", paste0("Gene_heart_ECs_nonGW24_AVS_selected_vECs_", gene, ".png"))
  if (file.exists(src)) file.copy(src, "nonGW24_AVS_selected_TradeSeq/Selected") else warning(paste("Not exists:", src))
})

sapply(Comb1_genes_L2_selected, function(gene) {
  src <- file.path("nonGW24_AVS_selected_TradeSeq", paste0("Gene_heart_ECs_nonGW24_AVS_selected_aECs_", gene, ".png"))
  if (file.exists(src)) file.copy(src, "nonGW24_AVS_selected_TradeSeq/Selected") else warning(paste("Not exists:", src))
})

sapply(Comb1_genes_L3_selected, function(gene) {
  src <- file.path("nonGW24_AVS_selected_TradeSeq", paste0("Gene_heart_ECs_nonGW24_AVS_selected_cECs_", gene, ".png"))
  if (file.exists(src)) file.copy(src, "nonGW24_AVS_selected_TradeSeq/Selected") else warning(paste("Not exists:", src))
})



## select genes
DotPlot(Seurat_Comb1, features = Comb1_genes_L1_selected, group.by = "Age") +
  coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust=1))
ggsave("DotPlot_heart_ECs_nonGW24_AVS_selected_vECs_selected.pdf", width=6, height=6)

DotPlot(Seurat_Comb1, features = Comb1_genes_L2_selected, group.by = "Age") +
  coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust=1))
ggsave("DotPlot_heart_ECs_nonGW24_AVS_selected_aECs_selected.pdf", width=6, height=6)

DotPlot(Seurat_Comb1, features = Comb1_genes_L3_selected, group.by = "Age") +
  coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust=1))
ggsave("DotPlot_heart_ECs_nonGW24_AVS_selected_cECs_selected.pdf", width=6, height=6)

Comb1_selected_smoothed <- generate_lineage_smoothed_matrix(Slingshot_Comb1, Seurat_Comb1, c(Comb1_genes_L1_selected, Comb1_genes_L2_selected, Comb1_genes_L3_selected))
Comb1_selected_smoothed_sorted <- plot_lineage_heatmap(Comb1_selected_smoothed, Slingshot_Comb1, order="max", outPrefix = "heart_ECs_nonGW24_AVS_selected_selected", height = 8)


rm(Slingshot_Comb1)
rm(Seurat_Comb1)
gc()


dir.create("nonGW24_AVS_selected/Selected", recursive=TRUE)
for (gene in Comb1_genes_L1_selected) {
  FeaturePlot(EN_heart_ECs_nonGW24_AVS_selected, gene, reduction = "umap_Trajectory")
  ggsave(paste0("nonGW24_AVS_selected/Selected/FeaturePlot_heart_ECs_nonGW24_AVS_selected_vECs_", gene, ".png"), width=5, height=4)
}
for (gene in Comb1_genes_L2_selected) {
  FeaturePlot(EN_heart_ECs_nonGW24_AVS_selected, gene, reduction = "umap_Trajectory")
  ggsave(paste0("nonGW24_AVS_selected/Selected/FeaturePlot_heart_ECs_nonGW24_AVS_selected_aECs_", gene, ".png"), width=5, height=4)
}
for (gene in Comb1_genes_L3_selected) {
  FeaturePlot(EN_heart_ECs_nonGW24_AVS_selected, gene, reduction = "umap_Trajectory")
  ggsave(paste0("nonGW24_AVS_selected/Selected/FeaturePlot_heart_ECs_nonGW24_AVS_selected_cECs_", gene, ".png"), width=5, height=4)
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
