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
EN_Batch01to08_Mural <- readRDS("EN_Batch01to08_5_Vasculature_Mural_pure.rds")
EN_brain_Mural_pure <- subset(EN_Batch01to08_Mural, Organ == "Brain")


pal <- c(RColorBrewer::brewer.pal(9,"Set1"),RColorBrewer::brewer.pal(9,"Pastel1"),RColorBrewer::brewer.pal(8,"Accent"),RColorBrewer::brewer.pal(8,"Set2"),RColorBrewer::brewer.pal(8,"Pastel2") ,scales::hue_pal()(8))


EN_brain_Mural_nonGW24 <- subset(EN_brain_Mural_pure, Age != "GW24")




################### EndMT
EN_brain_Mural_nonGW24_EMT <- subset(EN_brain_Mural_nonGW24, newL3_Subtype %in% c("pre_EndMT", "pre_transient", "SMCs", "PCs"))

EN_brain_Mural_nonGW24_EMT <- RunUMAP(EN_brain_Mural_nonGW24_EMT, dims = 1:40, #50
                                      reduction = "pca",
                                      reduction.name = "umap_Trajectory", reduction.key = "umapTrajectory_",
                                      min.dist = 0.001,
                                      metric = "euclidean",
                                      n.neighbors = 10)



### Clustering
EN_brain_Mural_nonGW24_EMT <- FindNeighbors(EN_brain_Mural_nonGW24_EMT, dims = 1:40, k.param = 20, reduction = "pca")

EN_brain_Mural_nonGW24_EMT <- FindClusters(EN_brain_Mural_nonGW24_EMT, graph.name = "CCA_snn", resolution = 1.5, algorithm = 1)

colnames(EN_brain_Mural_nonGW24_EMT@meta.data)[which(colnames(EN_brain_Mural_nonGW24_EMT@meta.data)=="CCA_snn_res.1.5")] <- "Trajectory_CCA_snn_res.1.5"


### C6 subclusering
EN_brain_Mural_nonGW24_EMT_C6 <- subset(EN_brain_Mural_nonGW24_EMT, Trajectory_CCA_snn_res.1.5 == 6)

EN_brain_Mural_nonGW24_EMT_C6 <- RunUMAP(EN_brain_Mural_nonGW24_EMT_C6, dims = 1:40, #50
                                         reduction = "pca",
                                         reduction.name = "umap_Trajectory", reduction.key = "umapTrajectory_",
                                         min.dist = 0.001,
                                         metric = "euclidean",
                                         n.neighbors = 10)


EN_brain_Mural_nonGW24_EMT_C6 <- FindNeighbors(EN_brain_Mural_nonGW24_EMT_C6, dims = 1:40, k.param = 20, reduction = "pca")

EN_brain_Mural_nonGW24_EMT_C6 <- FindClusters(EN_brain_Mural_nonGW24_EMT_C6, graph.name = "CCA_snn", resolution = 0.6, algorithm = 1)

EN_brain_Mural_nonGW24_EMT_C6$Trajectory_CCA_snn_res.0.6 <- paste0("6.", EN_brain_Mural_nonGW24_EMT_C6$CCA_snn_res.0.6)




### write-back
EN_brain_Mural_nonGW24_EMT$Trajectory_CCA_snn_res.1.5 <- as.character(EN_brain_Mural_nonGW24_EMT$Trajectory_CCA_snn_res.1.5)
EN_brain_Mural_nonGW24_EMT@meta.data[rownames(EN_brain_Mural_nonGW24_EMT_C6@meta.data), "Trajectory_CCA_snn_res.1.5"] <- as.character(EN_brain_Mural_nonGW24_EMT_C6$Trajectory_CCA_snn_res.0.6)






### Subset cluster
clusters_EMT_to_include <- c("6.1", as.character(c(4, 1, 3, 5, 6, 7)))
EN_brain_Mural_nonGW24_EMT_selected <- subset(EN_brain_Mural_nonGW24_EMT, Trajectory_CCA_snn_res.1.5 %in% clusters_EMT_to_include)


EN_brain_Mural_nonGW24_EMT_selected <- RunUMAP(EN_brain_Mural_nonGW24_EMT_selected, dims = 1:40, #50
                                               reduction = "pca",
                                               reduction.name = "umap_Trajectory", reduction.key = "umapTrajectory_",
                                               min.dist = 0.001,
                                               metric = "euclidean",
                                               n.neighbors = 10)



### Run Slingshot
Brain_Mural_nonGW24_EMT_selected_sce <- as.SingleCellExperiment(EN_brain_Mural_nonGW24_EMT_selected)
reducedDims(Brain_Mural_nonGW24_EMT_selected_sce) <- list(umap_Trajectory = EN_brain_Mural_nonGW24_EMT_selected@reductions$umap_Trajectory@cell.embeddings)


Brain_Mural_nonGW24_EMT_selected_clustering <- EN_brain_Mural_nonGW24_EMT_selected$Trajectory_CCA_snn_res.1.5
pal_clustering_nonGW24_EMT_selected <- pal[seq(length(table(Brain_Mural_nonGW24_EMT_selected_clustering)))]
names(pal_clustering_nonGW24_EMT_selected)<-sort(unique(Brain_Mural_nonGW24_EMT_selected_clustering))



Brain_Mural_nonGW24_EMT_selected_sce <- slingshot(Brain_Mural_nonGW24_EMT_selected_sce, clusterLabels = 'Trajectory_CCA_snn_res.1.5', reducedDim = 'umap_Trajectory', 
                                                  extend = "n", stretch=1, thresh=0.1)
Slingshot_nonGW24_EMT_selected<-SlingshotDataSet(Brain_Mural_nonGW24_EMT_selected_sce)



### Convert for velocity analysis
EN_brain_Mural_nonGW24_EMT_selected_diet<-DietSeurat(EN_brain_Mural_nonGW24_EMT_selected, scale.data = TRUE, dimreducs = c("umap_Trajectory"), assays = "CCA")

SaveH5Seurat(EN_brain_Mural_nonGW24_EMT_selected_diet, filename = "EN_brain_Mural_nonGW24_EMT_selected_trajectory.h5Seurat")

Convert("EN_brain_Mural_nonGW24_EMT_selected_trajectory.h5Seurat", dest = "h5ad")





pdf("Trajecotry_Brain_Mural_nonGW24_EMT_selected_UMAP_clustering_lineages_curves.pdf", width=6, height=6)
plot(reducedDims(Brain_Mural_nonGW24_EMT_selected_sce)$umap_Trajectory, col = pal_clustering_nonGW24_EMT_selected[Brain_Mural_nonGW24_EMT_selected_clustering], pch=16, asp = 1)
lines(SlingshotDataSet(Brain_Mural_nonGW24_EMT_selected_sce), lwd=3)
dev.off()






################### TradeSeq for trajectory
Brain_Seurat_Traj1 <- EN_brain_Mural_nonGW24_EMT_selected
Slingshot_Brain_Traj1 <- Slingshot_nonGW24_EMT_selected

## TradeSeq analysis
Brain_Traj1_counts <- as.matrix(Brain_Seurat_Traj1@assays$RNA@counts[Brain_Seurat_Traj1@assays$CCA@var.features, ])
Brain_Traj1_pseudotime <- slingPseudotime(Slingshot_Brain_Traj1, na = FALSE)
Brain_Traj1_cellWeights <- slingCurveWeights(Slingshot_Brain_Traj1)

Brain_Traj1_samples <- sample( colnames(Brain_Seurat_Traj1) , size = 1000)
Brain_Traj1_counts_sampling <- Brain_Traj1_counts[rowSums(Brain_Traj1_counts > 3) > ncol(Brain_Traj1_counts)/100, Brain_Traj1_samples]



set.seed(6)
Brain_Traj1_icMat_seed6 <- evaluateK(counts = Brain_Traj1_counts_sampling, k=3:10, pseudotime = Brain_Traj1_pseudotime[Brain_Traj1_samples, ], cellWeights=Brain_Traj1_cellWeights[Brain_Traj1_samples, ], verbose = T, parallel=T)
plot_evalutateK_results(Brain_Traj1_icMat_seed6)


set.seed(4)
Brain_Traj1_icMat_seed4 <- evaluateK(counts = Brain_Traj1_counts_sampling, k=3:10, pseudotime = Brain_Traj1_pseudotime[Brain_Traj1_samples, ], cellWeights=Brain_Traj1_cellWeights[Brain_Traj1_samples, ], verbose = T, parallel=T)
plot_evalutateK_results(Brain_Traj1_icMat_seed4)



## fit GAM
TradeSeq_Brain_Traj1 <- fitGAM(counts = Brain_Traj1_counts, pseudotime = Brain_Traj1_pseudotime, cellWeights=Brain_Traj1_cellWeights, 
                                   nknots=7, verbose = T, parallel=T)
# Note! Modify nknots


## pattern test
Brain_Traj1_patternRes <- patternTest(TradeSeq_Brain_Traj1)
Brain_Traj1_patternRes_filter <- filter_TradeSeq_results(Brain_Traj1_patternRes)

Brain_Traj1_patternRes_top100 <- rownames(Brain_Traj1_patternRes_filter)[1:min(100,nrow(Brain_Traj1_patternRes_filter))]

Brain_Traj1_smoothed <- generate_lineage_smoothed_matrix(Slingshot_Brain_Traj1, Brain_Seurat_Traj1, Brain_Traj1_patternRes_top100)
Brain_Traj1_smoothed_sorted <- plot_lineage_heatmap(Brain_Traj1_smoothed, Slingshot_Brain_Traj1, order="hclust", outPrefix = "brain_Mural_nonGW24_EMT_selected_top100", height=22)


## Identify lineage associated genes
Brain_Traj1_genes_L1_candidate <- Brain_Traj1_smoothed_sorted[c(1:3, 11:43)]
Brain_Traj1_genes_L2_candidate <- Brain_Traj1_smoothed_sorted[c(4:10, 16:100)]

for (g in Brain_Traj1_genes_L1_candidate) {
  plot_gene_expression(Slingshot_Brain_Traj1, Brain_Seurat_Traj1, g, outPrefix = "brain_Mural_nonGW24_EMT_selected_PCs")
}

for (g in Brain_Traj1_genes_L2_candidate) {
  plot_gene_expression(Slingshot_Brain_Traj1, Brain_Seurat_Traj1, g, outPrefix = "brain_Mural_nonGW24_EMT_selected_SMCs")
}


Brain_Traj1_genes_L1_selected <- c("JUNB", "NFKBIA", "KCNJ8", "TESC", "CCL2", "CEBPD")
Brain_Traj1_genes_L2_selected <- c("ACTA2", "BCL2", "CACNA1C", "DLC1", "DMD", "EBF1", "ENOX1", "EPS8", "FBXL7", "INPP4B", "MEF2C", "PARD3B", "PDE3A", "PDE4D", "PDLIM5", "PEAK1", "PRKG1", "SLC8A1", "SLIT3", "SOX5", "SYNPO2", "SYTL2", "TAGLN", "TEX41", "TRPC6", "UACA", "UTRN", "ZFHX3")

dir.create("Brain_nonGW24_EMT_selected_TradeSeq/Selected", recursive=TRUE)
sapply(Brain_Traj1_genes_L1_selected, function(gene) {
  src <- file.path("Brain_nonGW24_EMT_selected_TradeSeq", paste0("Gene_brain_Mural_nonGW24_EMT_selected_PCs_", gene, ".png"))
  if (file.exists(src)) file.copy(src, "Brain_nonGW24_EMT_selected_TradeSeq/Selected") else warning(paste("Not exists:", src))
})

sapply(Brain_Traj1_genes_L2_selected, function(gene) {
  src <- file.path("Brain_nonGW24_EMT_selected_TradeSeq", paste0("Gene_brain_Mural_nonGW24_EMT_selected_SMCs_", gene, ".png"))
  if (file.exists(src)) file.copy(src, "Brain_nonGW24_EMT_selected_TradeSeq/Selected") else warning(paste("Not exists:", src))
})


### select genes
DotPlot(Brain_Seurat_Traj1, features = Brain_Traj1_genes_L1_selected, group.by = "Age") +
  coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust=1))
ggsave("DotPlot_brain_Mural_nonGW24_EMT_selected_PCs_selected.pdf", width=6, height=6)

DotPlot(Brain_Seurat_Traj1, features = Brain_Traj1_genes_L2_selected, group.by = "Age") +
  coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust=1))
ggsave("DotPlot_brain_Mural_nonGW24_EMT_selected_SMCs_selected.pdf", width=6, height=6)


Brain_Traj1_selected_smoothed <- generate_lineage_smoothed_matrix(Slingshot_Brain_Traj1, Brain_Seurat_Traj1, c(Brain_Traj1_genes_L1_selected, Brain_Traj1_genes_L2_selected))
Brain_Traj1_selected_smoothed_sorted <- plot_lineage_heatmap(Brain_Traj1_selected_smoothed, Slingshot_Brain_Traj1, order="max", outPrefix = "brain_Mural_nonGW24_EMT_selected_selected", height=6)



rm(Brain_Seurat_Traj1)
rm(Slingshot_Brain_Traj1)
gc()


dir.create("Brain_nonGW24_EMT_selected/Selected", recursive=TRUE)
DefaultAssay(EN_brain_Mural_nonGW24_EMT_selected) <- "RNA"

for (gene in Brain_Traj1_genes_L1_selected) {
  FeaturePlot(EN_brain_Mural_nonGW24_EMT_selected, gene, reduction = "umap_Trajectory")
  ggsave(paste0("Brain_nonGW24_EMT_selected/Selected/FeaturePlot_brain_mural_nonGW24_EMT_selected_PCs_", gene, ".png"), width=5, height=4)
}
for (gene in Brain_Traj1_genes_L2_selected) {
  FeaturePlot(EN_brain_Mural_nonGW24_EMT_selected, gene, reduction = "umap_Trajectory")
  ggsave(paste0("Brain_nonGW24_EMT_selected/Selected/FeaturePlot_brain_mural_nonGW24_EMT_selected_SMCs_", gene, ".png"), width=5, height=4)
}

write.table(EN_brain_Mural_nonGW24_EMT_selected@meta.data[, c(28:32)], "EN_brain_mural_ECs_nonGW24_EMT_selected_metaData.txt", sep="\t", row.names = TRUE, col.names = TRUE, quote=FALSE)




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
