library(Seurat)
library(ggplot2)
library(EWCE)
library(MAGMA.Celltyping)
library(dplyr)
library(tidyr)
library("clusterProfiler")

dir_software <- "D:/OneDrive/SciLifeLab/Materials/Softwares/MAGMA"
dir_working <- "EWCE_CTD"

adata <- readRDS("total_CD31filt_double_filtered_iNMF_wclusters_BHproj_final2.rds")


expData <- adata@assays$RNA@data
metaData <- adata@meta.data
metaData[metaData$Organ=="lung", "Organ"] <- "Lung"

genes_adata <- rownames(adata)
rm(adata)
gc()


annotLevels <- list(CellType = metaData$Annotation_res0.5)

# Generate EWCE-CTD file
fNames <- EWCE::generate_celltype_data(
  exp = expData,
  annotLevels = annotLevels,
  groupName = paste0("EN_organ_", "all"),
  savePath = dir_working,
  input_species = "human"
)





############ MAGMA
organs <- c("Brain", "Heart", "Eye", "Kidney", "Liver", "Lung", "Systemic", "Neuropsychiatric")

# ====================================== Generate EWCE-CTD files by organ ======================================
# Iterate through each organ to generate organ-specific CTD files
for (tmp_organ in organs) {
  # Step 1: Filter all cells of the current organ (keep only samples from target organ)
  meta_subset_organ <- metaData %>%
    filter(Organ == tmp_organ)
  
  # Step 2: Extract expression data for current organ
  expData_subset_organ <- expData[, rownames(meta_subset_organ), drop = FALSE]
  
  # Step 3: Construct annotation levels
  annotLevels_organ <- list(CellType = meta_subset_organ$Annotation_res0.5)
  
  # Step 4: Generate EWCE-CTD file for the current organ
  fNames_organ <- EWCE::generate_celltype_data(
    exp = expData_subset_organ,
    annotLevels = annotLevels_organ,
    groupName = paste0("EN_organ_", tmp_organ),
    savePath = dir_working,
    input_species = "human"
  )
  
  rm(meta_subset_organ, expData_subset_organ, annotLevels_organ, fNames_organ)
  gc()
}

rm(expData)
gc()

# 2. Load MAGMA path files for all organs
for (tmp_organ in organs) {
  load(paste0("MAGMA_files/20251121Batch01to08_MAGMA_paths_", tmp_organ, ".RData"))
}


# ====================================== Outer loop: iterate through all organs (calculate all cell types within each organ) ======================================
for (tmp_organ in organs) {
  # Skip organs without CTD files
  tmp_ctd_path <- paste0(dir_working, "/ctd_EN_organ_", tmp_organ, ".rda")
  if (!file.exists(tmp_ctd_path)) {
    print(paste("Skipping organ", tmp_organ, ": no corresponding CTD file"))
    next
  }
  
  # Create result directory for current organ
  tmp_dir_working <- paste0("./MAGMA_results_", tmp_organ, "/")
  if (!dir.exists(tmp_dir_working)) dir.create(tmp_dir_working, recursive = TRUE)

  # Load CTD file for current organ
  load(tmp_ctd_path)
  print(paste("Successfully loaded CTD file for organ", tmp_organ, ", containing cell types:", paste(unique(ctd[[1]]$annot), collapse = ", ")))

  ##################################### Calculation
  merged_results_celltypes_current <- list()
  merged_results_celltypes_current_Linear <- list()
  merged_results_celltypes_current_Top <- list()

  # Get MAGMA paths for current organ
  tmp_magma_paths <- get(paste0("magma_paths_", tmp_organ))

  # Iterate through all GWAS traits for current organ
  for (trait in names(tmp_magma_paths)) {
    CTD_name <- paste0("byOrgan_", tmp_organ, "_", trait)

    # Run MAGMA cell type enrichment analysis
    tmp_assoc <- MAGMA.Celltyping::celltype_associations_pipeline(
      magma_dirs = dirname(file.path(tmp_magma_paths[[trait]])),
      ctd = ctd,
      ctd_species = "human",
      ctd_name = CTD_name,
      run_linear = TRUE,
      run_top10 = TRUE,
      save_dir = paste(tmp_dir_working, sep='/')
    )

    # Extract Linear results
    if (length(tmp_assoc[[1]][["ctAssocsLinear"]]) > 0) {
      tmp_linear <- MAGMA.Celltyping::merge_results(
        MAGMA_results = tmp_assoc,
        filetype = "ctAssocsLinear",
        level = 1,
        save_dir = paste(tmp_dir_working, "Linear", sep='/')
      )
      merged_results_celltypes_current_Linear[[trait]] <- tmp_linear
    }

    # Extract Top results
    if (length(tmp_assoc[[1]][["ctAssocsTop"]]) > 0) {
      tmp_top <- MAGMA.Celltyping::merge_results(
        MAGMA_results = tmp_assoc,
        filetype = "ctAssocsTop",
        level = 1,
        save_dir = paste(tmp_dir_working, "Top", sep='/')
      )
      merged_results_celltypes_current_Top[[trait]] <- tmp_top
    }

    merged_results_celltypes_current[[trait]] <- tmp_assoc
  }

  # Assign calculation results for current organ to global variables
  assign(paste0("merged_results_", tmp_organ), merged_results_celltypes_current)
  assign(paste0("merged_results_", tmp_organ, "_Linear"), merged_results_celltypes_current_Linear)
  assign(paste0("merged_results_", tmp_organ, "_Top"), merged_results_celltypes_current_Top)


  ############################## Combination
  combined_results_celltypes <- list()
  selected_results_celltypes <- list()
  
  # Merge Linear and Top type results
  for (suffix in c("Linear", "Top")) {
    tmp_merged_all <- get(paste0("merged_results_", tmp_organ, "_", suffix))
    tmp_combined <- lapply(names(tmp_merged_all), function(trait) {
      data <- tmp_merged_all[[trait]]
      if (is.null(data) || nrow(data) == 0) return(NULL)
      
      data %>%
        mutate(
          GWAS = trait,
          Organ = tmp_organ,
          mlog10_P = -log10(P),
          mlog10_FDR = -log10(FDR)
        ) %>%
        select(GWAS, Organ, Celltype_id, P, FDR, EnrichmentMode, OBS_GENES, BETA, BETA_STD, SE, mlog10_P, mlog10_FDR) %>%
        rename(Celltype = Celltype_id)
    }) %>%
      bind_rows()
    
    combined_results_celltypes[[suffix]] <- tmp_combined
  }
  
  # Filter optimal result for each (GWAS-cell type) pair (minimum FDR)
  tmp_all_combined <- bind_rows(combined_results_celltypes[["Linear"]], combined_results_celltypes[["Top"]])
  if (nrow(tmp_all_combined) > 0) {
    tmp_selected <- tmp_all_combined %>%
      group_by(GWAS, Celltype) %>%
      slice_min(FDR, n = 1, with_ties = FALSE) %>%
      ungroup()
  } else {
    tmp_selected <- NULL
    print(paste("Warning: organ", tmp_organ, "has no valid enrichment results"))
  }
  
  selected_results_celltypes <- tmp_selected
  
  # Assign merged results for current organ to global variables
  assign(paste0("combined_results_", tmp_organ), combined_results_celltypes)
  assign(paste0("selected_results_", tmp_organ), selected_results_celltypes)
  
  
  ############################## Save all results for current organ
  save(
    merged_results_celltypes_current, merged_results_celltypes_current_Linear, merged_results_celltypes_current_Top,
    combined_results_celltypes, selected_results_celltypes,
    file = paste0(tmp_dir_working, "MAGMA_results_", tmp_organ, ".RData")
  )
  
  # Clean up temporary variables for current organ (including CTD) to free memory
  rm(
    ctd, merged_results_celltypes_current, merged_results_celltypes_current_Linear, merged_results_celltypes_current_Top,
    combined_results_celltypes, selected_results_celltypes, tmp_selected, tmp_all_combined, tmp_magma_paths
  )
  gc()
}

# Clean up global temporary variables
rm(list = ls(pattern = "tmp_"))
gc()





## additional analyses
# 1. Brain, Eye and Neuropsychiatric disease
tmp_organ <- "Brain"
tmp_magma_paths <- magma_paths_Eye
tmp_magma_paths <- magma_paths_Neuropsychiatric
for (trait in names(tmp_magma_paths)) {
  # Run MAGMA
}

merged_results_Brain <- c(merged_results_Brain, merged_results_celltypes_current)
merged_results_Brain_Linear <- c(merged_results_Brain_Linear, merged_results_celltypes_current_Linear)
merged_results_Brain_Top <- c(merged_results_Brain_Top, merged_results_celltypes_current_Top)

save(
  merged_results_Brain, merged_results_Brain_Linear, merged_results_Brain_Top,
  combined_results_Brain, selected_results_Brain,
  file = paste0(tmp_dir_working, "MAGMA_results_Brain_additional.RData")
)




# 2. All organs, Systemic disease
tmp_organ <- "all"
tmp_magma_paths <- magma_paths_Systemic


# 3. Combine all organs


## disease-organ corresponding
list_organ_disease <- list()
for (tmp_organ in organs) {
  tmp_paths <- load(paste0("../20251121MAGMA/MAGMA_files/20251121Batch01to08_MAGMA_paths_", tmp_organ, ".RData"))
  list_organ_disease[[tmp_organ]] <- names(get(tmp_paths))
}
matrix_organ_disease <- do.call(rbind, sapply(names(list_organ_disease), function(x) cbind(x, list_organ_disease[[x]])))
colnames(matrix_organ_disease) <- c("Organ", "Disease")
rownames(matrix_organ_disease) <- matrix_organ_disease[,"Disease"]


byCellType_combine <- rbind(selected_results_Brain, selected_results_Heart, selected_results_Kidney, selected_results_Liver, selected_results_Lung, selected_results_all)
byCellType_combine$Organ_disease <- matrix_organ_disease[byCellType_combine$GWAS, "Organ"]
byCellType_combine <- byCellType_combine[order(byCellType_combine$Organ_disease),]
byCellType_combine <- byCellType_combine[! (byCellType_combine$Celltype %in% c("Immune_like_1", "Immune_like_2", "Neural_Prog_like")), ]
byCellType_combine <- byCellType_combine[byCellType_combine$Organ_disease != "Eye", ]


results2heatmap(byCellType_combine, label="byOrgan_combine", value_type="BETA", BETA_symmetric = TRUE, BETA_limit = 2, output_format = "pdf", 
                fontsize = 18, height=800, width=650)





########### Specific organs
byCellType_combine[byCellType_combine$FDR<0.1,1:6]


# 1. Brain
load(paste0(dir_working, "/ctd_EN_organ_Brain.rda"))
ctd_Brain <- ctd

specificity_Brain_aEC <- get_gene_specificity_from_ctd(ctd_Brain, "EC_Arterial")


genesOut_Brain_StrokeLargeArtery <- read.table("D:/OneDrive/SciLifeLab/Endothelial_cells/20251121MAGMA/MAGMA_files/Brain/MAGMA_Files/Large_artery_stroke_20231025.tsv.35UP.10DOWN/Large_artery_stroke_20231025.tsv.35UP.10DOWN.genes.out", header = TRUE, row.names = 1)
genesOut_Brain_StrokeLargeArtery_df <- bitr_dataframe(genesOut_Brain_StrokeLargeArtery, species="Human", fromType = "ENTREZID", toType="SYMBOL")

genes_Brain_StrokeLargeArtery <- intersect(rownames(genesOut_Brain_StrokeLargeArtery_df), genes_adata)
genesOut_Brain_StrokeLargeArtery_combine <- data.frame(genesOut_Brain_StrokeLargeArtery_df[genes_Brain_StrokeLargeArtery, "ZSTAT", drop=FALSE], Specificity=specificity_Brain_aEC[genes_Brain_StrokeLargeArtery,])




# Custom gene list
genes_LAA <- c("FURIN", "SLCO1B1", "LIPA", "FLT1", "TSPAN2")

# Call function (label custom genes)
result_lm_label <- plot_magma_regression(
  df = genesOut_Brain_StrokeLargeArtery_combine,
  cell_type = "EC_Arterial",
  disease = "Large Arterial Stroke",
  organ = "Brain",
  reg_type = "linear",
  custom_gene_list = genes_LAA,
  show_reg_label = FALSE,
  use_density_color = FALSE,
  point_color = "#2E86AB",
  save_path = "./Brain_LAS_Linear_Regression_label.pdf"
)





# 2. Heart
load(paste0(dir_working, "/ctd_EN_organ_Heart.rda"))
ctd_Heart <- ctd

specificity_Heart_cEC_CNS <- get_gene_specificity_from_ctd(ctd_Heart, "EC_Capillary_CNS")


genesOut_Heart_CardiovascularDisease <- read.table("D:/OneDrive/SciLifeLab/Endothelial_cells/20251121MAGMA/MAGMA_files/Heart/MAGMA_Files/Cardiovascular_disease_20240301.tsv.35UP.10DOWN/Cardiovascular_disease_20240301.tsv.35UP.10DOWN.level1.byOrgan_Heart_Cardiovascular_disease_top10.gsa.genes.out.txt", header = TRUE, row.names = 1)
genesOut_Heart_CardiovascularDisease_df <- bitr_dataframe(genesOut_Heart_CardiovascularDisease, species="Human", fromType = "ENTREZID", toType="SYMBOL")

genes_Heart_CardiovascularDisease <- intersect(rownames(genesOut_Heart_CardiovascularDisease_df), genes_adata)
genesOut_Heart_CardiovascularDisease_combine <- data.frame(genesOut_Heart_CardiovascularDisease_df[genes_Heart_CardiovascularDisease,7, drop=FALSE], Specificity=specificity_Heart_cEC_CNS[genes_Heart_CardiovascularDisease,])




# Custom gene list
genes_CAD <- c("LEPR", "PCSK1", "CDH11", "ADA2", "ITLN1", "GSTZ1", "MTHFS")

# Call function (label custom genes)
result_top_label <- plot_magma_top10(
  df = genesOut_Heart_CardiovascularDisease_combine,
  cell_type = "EC_Capillary_CNS",
  disease = "Cardiovascular Disease",
  organ = "Heart",
  custom_gene_list = genes_CAD,
  save_path = "./Heart_CAD_Top10_label.pdf"
)





# 3. Lung
load(paste0(dir_working, "/ctd_EN_organ_Lung.rda"))
ctd_Lung <- ctd
rm(ctd)

specificity_Lung_cEC_CNS <- get_gene_specificity_from_ctd(ctd_Lung, "Mural_Pericytes")


genesOut_Lung_PulmonaryFibrosis <- read.table("D:/OneDrive/SciLifeLab/Endothelial_cells/20251121MAGMA/MAGMA_files/Lung/MAGMA_Files/Pulmonary_fibrosis_20251127.tsv.35UP.10DOWN/Pulmonary_fibrosis_20251127.tsv.35UP.10DOWN.level1.byOrgan_Lung_Pulmonary_fibrosis_top10.gsa.genes.out.txt", header = TRUE, row.names = 1)
genesOut_Lung_PulmonaryFibrosis_df <- bitr_dataframe(genesOut_Lung_PulmonaryFibrosis, species="Human", fromType = "ENTREZID", toType="SYMBOL")

genes_Lung_PulmonaryFibrosis <- intersect(rownames(genesOut_Lung_PulmonaryFibrosis_df), genes_adata)
genesOut_Lung_PulmonaryFibrosis_combine <- data.frame(genesOut_Lung_PulmonaryFibrosis_df[genes_Lung_PulmonaryFibrosis,7, drop=FALSE], Specificity=specificity_Lung_cEC_CNS[genes_Lung_PulmonaryFibrosis,])



# Custom gene list
genes_PF <- c("MAPT", "FUT3", "TTYH1")

# Call function (label custom genes)
result_Lung_top_label <- plot_magma_top10(
  df = genesOut_Lung_PulmonaryFibrosis_combine,
  cell_type = "Mural_Pericytes",
  disease = "Pulmonary Fibrosis",
  organ = "Lung",
  use_density_color = FALSE,
  custom_gene_list = genes_PF,
  save_path = "./Lung_PF_Top10_label.pdf"
)



### 4. Compare
to_compare_brain <- byCellType_combine[byCellType_combine$GWAS %in% c("Stroke_large_artery", "Brain_injury"),]


p <- ggplot(to_compare_brain, aes(x = Celltype, y = mlog10_FDR, fill = GWAS)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8, alpha = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 1) +
  scale_fill_manual(
    values = c("Stroke_large_artery" = "#2E86AB", 
               "Brain_injury" = "#F1A299"),
    name = "Disease Type"
  ) +
  labs(
    x = "Cell Subtype",
    y = "-log10(FDR)",
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray50"),
    plot.caption = element_text(size = 9, hjust = 0, color = "gray40"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 9, angle = 60, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

ggsave(
  filename = "./Brain_FDR_Comparison_Barplot.pdf",
  plot = p,
  width = 8, height = 5, dpi = 300, bg = "white"
)




########################### Function

results2heatmap <- function(results_combine, 
                            label = "MAGMA", 
                            height = 650, 
                            width = 480,
                            cluster_rows = FALSE,
                            cluster_cols = FALSE,
                            show_rownames = TRUE,
                            show_colnames = TRUE,
                            main_title = NULL,
                            color_palette = NULL,
                            fontsize = 8,
                            fontsize_row = 10,
                            fontsize_col = 10,
                            output_format = "png",
                            return_plot = FALSE,
                            value_type = "FDR",  # "FDR" or "BETA"
                            BETA_symmetric = TRUE,  # whether BETA values should have symmetric colors
                            BETA_limit = NULL,  # manually set BETA value color range
                            breaks = NULL) {  # manually set color breaks
  
  # Check required packages
  require("pheatmap")
  require("dplyr")
  require("tidyr")
  require("tibble")
  require("RColorBrewer")
  
  # Check value_type parameter
  value_type <- toupper(value_type)
  if(!value_type %in% c("FDR", "BETA")) {
    stop("value_type parameter must be 'FDR' or 'BETA'")
  }
  
  # Check input data format
  required_cols <- c("Celltype", "GWAS", "mlog10_FDR", "FDR")
  if(value_type == "BETA") {
    required_cols <- c(required_cols, "BETA")
    if(!"BETA" %in% colnames(results_combine)) {
      stop("When value_type='BETA', input data must contain 'BETA' column")
    }
  }
  
  if(!all(required_cols %in% colnames(results_combine))) {
    stop(paste("Input data must contain the following columns:", paste(required_cols, collapse = ", ")))
  }
  
  # select value to plot based on value_type
  value_col <- ifelse(value_type == "FDR", "mlog10_FDR", "BETA")
  
  # convert to heatmap matrix
  heatmap_matrix <- results_combine %>%
    pivot_wider(
      id_cols = GWAS,
      names_from = Celltype,
      values_from = .data[[value_col]]
    ) %>%
    column_to_rownames("GWAS")
  
  # Create significance marker matrix
  sig_matrix <- results_combine %>%
    mutate(
      sig = case_when(
        FDR < 0.001 ~ "***",
        FDR < 0.01  ~ "**",
        FDR < 0.05  ~ "*",
        FDR < 0.1   ~ "·",
        TRUE ~ ""
      )
    ) %>%
    pivot_wider(
      id_cols = GWAS,
      names_from = Celltype,
      values_from = sig
    ) %>%
    column_to_rownames("GWAS")
  
  # Replace NA in sig_matrix with empty string
  sig_matrix[is.na(sig_matrix)] <- ""
  
  # Set default color gradient
  if(is.null(color_palette)) {
    if(value_type == "FDR") {
      # FDR heatmap uses low-saturation blue gradient
      color_palette <- colorRampPalette(c("#F7FBFF", "#C6DBEF", "#9ECAE1", 
                                          "#6BAED6", "#4292C6", "#2171B5", 
                                          "#08519C", "#08306B"))(100)
    } else if(value_type == "BETA") {
      # BETA heatmap uses divergent color scheme
      if(BETA_symmetric) {
        color_palette <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", 
                                            "#D1E5F0", "#F7F7F7", "#FDDBC7", 
                                            "#F4A582", "#D6604D", "#B2182B"))(100)
      } else {
        color_palette <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", 
                                            "#D1E5F0", "#F7F7F7", "#FDDBC7", 
                                            "#F4A582", "#D6604D", "#B2182B"))(100)
      }
    }
  }
  
  # Set default title
  if(is.null(main_title)) {
    if(value_type == "FDR") {
      main_title <- paste("GWAS vs Celltype - mlog10(FDR) Heatmap", 
                          if(label != "MAGMA") paste("-", label))
    } else if(value_type == "BETA") {
      main_title <- paste("GWAS vs Celltype - BETA Value Heatmap", 
                          if(label != "MAGMA") paste("-", label))
    }
  }
  
  # Initialize annotation-related variables
  annotation_row <- NULL
  annotation_colors <- NULL
  
  # Create row annotations only when Organ_disease column exists
  if("Organ_disease" %in% colnames(results_combine)) {
    # Build low-saturation color list (for Organ_disease annotation)
    low_sat_palettes <- list(
      brewer.pal(9, "Pastel1"),
      brewer.pal(8, "Pastel2"),
      brewer.pal(12, "Set3"),
      brewer.pal(8, "Set2"),
      brewer.pal(8, "Accent"),
      c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
        "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
    )
    all_low_sat_colors <- unique(unlist(low_sat_palettes))
    if(length(all_low_sat_colors) < 100) {
      custom_low_sat_colors <- c(
        "#B3CDE3", "#8C96C6", "#8856A7", "#810F7C",
        "#FED9A6", "#FBB4AE", "#DECBE4", "#CCEBC5",
        "#CCEBC5", "#A8DDB5", "#7BCCC4", "#4EB3D3",
        "#2B8CBE", "#0868AC", "#084081", "#F0F9E8",
        "#FBB4AE", "#FDCDAC", "#FEE0B6", "#F1B6DA",
        "#E6F5D0", "#B8E186", "#7FBC41", "#4D9221",
        "#D8DAEB", "#B2ABD2", "#8073AC", "#542788",
        "#E6F5C9", "#B8E186", "#7FBC41", "#4D9221"
      )
      all_low_sat_colors <- unique(c(all_low_sat_colors, custom_low_sat_colors))
    }
    
    # Create row annotations
    annotation_row <- results_combine %>%
      select(GWAS, Organ_disease) %>%
      distinct() %>%
      column_to_rownames("GWAS")
    
    # Check if there are duplicate GWAS corresponding to different Organ_disease
    if(nrow(annotation_row) != length(unique(results_combine$GWAS))) {
      warning("Some GWAS correspond to multiple Organ_disease values, taking the first match")
    }
    
    # Assign colors to Organ_disease
    organ_disease_levels <- unique(annotation_row$Organ_disease)
    organ_disease_colors <- all_low_sat_colors[1:length(organ_disease_levels)]
    names(organ_disease_colors) <- organ_disease_levels
    annotation_colors <- list(Organ_disease = organ_disease_colors)
  }
  
  # Handle BETA symmetry issue
  if(value_type == "BETA" && BETA_symmetric && is.null(breaks)) {
    BETA_values <- na.omit(as.vector(as.matrix(heatmap_matrix)))
    if(length(BETA_values) > 0) {
      max_abs <- max(abs(BETA_values))
      if(is.null(BETA_limit)) {
        BETA_limit <- max_abs
      }
      if(BETA_limit == 0) {
        BETA_limit <- 1
      }
      breaks <- seq(-BETA_limit, BETA_limit, length.out = length(color_palette) + 1)
    }
  }
  if(value_type == "BETA" && !is.null(BETA_limit) && is.null(breaks)) {
    breaks <- seq(-BETA_limit, BETA_limit, length.out = length(color_palette) + 1)
  }
  
  # Create heatmap parameters
  heatmap_params <- list(
    mat = heatmap_matrix,
    display_numbers = sig_matrix,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    number_color = "black",
    main = main_title,
    fontsize = fontsize,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    color = color_palette,
    na_col = "gray90",
    border_color = ifelse(nrow(heatmap_matrix) * ncol(heatmap_matrix) > 10000, 
                          NA, "grey60"),
    scale = "none",
    silent = return_plot,
    na_number = ""
  )
  
  # Add to parameters only when annotations exist
  if(!is.null(annotation_row)) {
    heatmap_params$annotation_row <- annotation_row
  }
  if(!is.null(annotation_colors)) {
    heatmap_params$annotation_colors <- annotation_colors
  }
  
  # If breaks are set, add to parameters
  if(!is.null(breaks)) {
    heatmap_params$breaks <- breaks
  }
  
  # Create heatmap
  heatmap_plot <- do.call(pheatmap, heatmap_params)
  
  # Save image
  if(!return_plot) {
    filename_base <- paste0("Heatmap_", label, "_", value_type)
    
    if(output_format == "png") {
      filename <- paste0(filename_base, ".png")
      png(filename, height = height, width = width, res = NA)
      print(heatmap_plot)
      dev.off()
      message(paste("Heatmap saved as:", filename))
    } else if(output_format == "pdf") {
      filename <- paste0(filename_base, ".pdf")
      pdf(filename, height = height/72, width = width/72)
      print(heatmap_plot)
      dev.off()
      message(paste("Heatmap saved as:", filename))
    }
  }
  
  # Return graphics object (optional)
  if(return_plot) {
    return(heatmap_plot)
  }
}

bitr_dataframe <- function(raw_df, species="Mouse", fromType = "ENSEMBL", toType="SYMBOL") {
  require("clusterProfiler")
  if (species=="Mouse") {
    df_gene_conversion <- bitr(rownames(raw_df), fromType = fromType, toType = toType, OrgDb = 'org.Mm.eg.db')
  } else if (species=="Human") {
    df_gene_conversion <- bitr(rownames(raw_df), fromType = fromType, toType = toType, OrgDb = 'org.Hs.eg.db')
  }
  
  df_gene_conversion_unique <- df_gene_conversion[!duplicated(df_gene_conversion[,1]) & !duplicated(df_gene_conversion[,2]), ]
  
  new_df_count <- raw_df[df_gene_conversion_unique[,1],]
  rownames(new_df_count) <- df_gene_conversion_unique[,2]
  
  return(new_df_count)
}

get_gene_specificity_from_ctd <- function(ctd, col_names, gene_list = NULL) {
  # Load required dependency packages
  if (!require(Matrix)) {
    stop("Please install and load Matrix package first: install.packages('Matrix'), library(Matrix)")
  }
  if (!require(dplyr)) {
    stop("Please install and load dplyr package first: install.packages('dplyr'), library(dplyr)")
  }
  
  # ===================== Step 1: Input validation =====================
  # Validate CTD object
  if (missing(ctd) || is.null(ctd) || is.null(ctd[[1]]) || is.null(ctd[[1]][["specificity"]])) {
    stop("Input ctd invalid: please provide a ctd object containing [[1]][['specificity']] matrix")
  }
  
  # Extract specificity matrix
  spec_matrix <- ctd[[1]][["specificity"]]
  
  # Validate if column names exist
  if (missing(col_names)) {
    stop("Please specify column names to extract (col_names parameter)!")
  }
  col_names <- as.character(col_names)
  valid_cols <- intersect(col_names, colnames(spec_matrix))
  invalid_cols <- setdiff(col_names, colnames(spec_matrix))
  
  # Warn about invalid column names and keep only valid columns
  if (length(invalid_cols) > 0) {
    warning(paste(
      "The following", length(invalid_cols), "column names do not exist! Automatically filtered:",
      paste(invalid_cols, collapse = ", ")
    ))
  }
  if (length(valid_cols) == 0) {
    stop("No valid column names! Please check col_names parameter, valid column names:", paste(colnames(spec_matrix), collapse = ", "))
  }
  
  # ===================== Step 2: Extract all data for specified columns and convert to data frame =====================
  # Extract valid columns and convert to data frame (rows = genes, columns = cell types, values = specificity values)
  spec_matrix_sub <- as.matrix(spec_matrix[, valid_cols, drop = FALSE])
  spec_matrix_sub <- spec_matrix_sub[!apply(is.na(spec_matrix_sub), 1, all), , drop = FALSE]
  spec_df <- as.data.frame(spec_matrix_sub)
  rownames(spec_df) <- rownames(spec_matrix_sub)
  
  # ===================== Step 3: Filter based on gene list =====================
  if (!is.null(gene_list)) {
    # Convert gene list to character type
    gene_list <- as.character(gene_list)
    
    # Keep only genes present in spec_df
    valid_genes <- intersect(gene_list, rownames(spec_df))
    invalid_genes <- setdiff(gene_list, rownames(spec_df))
    
    # Warn about invalid genes
    if (length(invalid_genes) > 0) {
      warning(paste(
        "The following", length(invalid_genes), "genes were not found in CTD and have been automatically deleted:",
        paste(head(invalid_genes, 10), collapse = ", "),
        if (length(invalid_genes) > 10) paste("... (total", length(invalid_genes), "genes)") else ""
      ))
    }
    
    # Filter valid genes
    spec_df <- spec_df[valid_genes, , drop = FALSE]
    spec_df <- spec_df[match(valid_genes, rownames(spec_df)), , drop = FALSE]
  } else {
    spec_df <- spec_df[order(rownames(spec_df)), , drop = FALSE]
    message(paste("No gene list specified, returning specificity values for", paste(valid_cols, collapse = ", "), "columns, total", nrow(spec_df), "genes"))
  }
  
  # ===================== Step 4: Format output (ensure correct numeric types) =====================
  spec_df <- dplyr::mutate_all(spec_df, as.numeric)
  
  # ===================== Step 5: Return results =====================
  return(spec_df)
}

plot_magma_regression <- function(df, 
                                  cell_type, 
                                  disease, 
                                  organ,
                                  custom_gene_list = NULL,
                                  top_pct = 0.1,
                                  reg_type = "linear",
                                  show_reg_label = TRUE,
                                  use_density_color = TRUE,
                                  point_color = "#2E86AB",
                                  point_alpha = 0.3,
                                  ci_color = "gray80",
                                  line_color = "#2874A6",
                                  width = 9, 
                                  height = 7,
                                  save_path = "./MAGMA_ZSTAT_Specificity_Regression.png") {
  if (!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
  if (!require(dplyr)) {install.packages("dplyr"); library(dplyr)}
  if (!require(ggrepel)) {install.packages("ggrepel"); library(ggrepel)}
  if (use_density_color) {
    if (!require(ggpointdensity)) {install.packages("ggpointdensity"); library(ggpointdensity)}
  }
  if (reg_type == "lm_robust") {
    if (!require(robustbase)) {install.packages("robustbase"); library(robustbase)}
  }
  
  # ===================== Custom utility functions =====================
  cap_first <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }
  
  # ===================== Step 1: Data preprocessing =====================
  # Calculate top percentage threshold (only consider positive values)
  zstat_positive <- df$ZSTAT[df$ZSTAT > 0]
  specificity_positive <- df$Specificity[df$Specificity > 0]
  
  # Check if sufficient positive value data exists
  if (length(zstat_positive) == 0) {
    stop("ZSTAT has no positive value data, cannot calculate top percentage threshold")
  }
  if (length(specificity_positive) == 0) {
    stop("Specificity has no positive value data, cannot calculate top percentage threshold")
  }
  
  zstat_top_threshold <- quantile(zstat_positive, probs = 1 - top_pct, na.rm = TRUE)
  specificity_top_threshold <- quantile(specificity_positive, probs = 1 - top_pct, na.rm = TRUE)
  
  
  plot_data <- df %>%
    tibble::rownames_to_column("Gene") %>%
    filter(!is.na(ZSTAT) & !is.na(Specificity)) %>%
    mutate(
      is_custom_gene = ifelse(!is.null(custom_gene_list) & Gene %in% custom_gene_list, "Custom", "Other"),
      is_label = ifelse(is_custom_gene == "Custom", "Labeled", "Unlabeled"),
      is_zstat_top = ifelse(ZSTAT >= zstat_top_threshold, TRUE, FALSE),
      is_specificity_top = ifelse(Specificity >= specificity_top_threshold, TRUE, FALSE),
      is_top_gene = case_when(
        is_zstat_top & is_specificity_top ~ "Both",
        is_zstat_top & !is_specificity_top ~ "ZSTAT_top",
        !is_zstat_top & is_specificity_top ~ "Specificity_top",
        TRUE ~ "Neither"
      )
    ) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("Gene")
  
  reg_data <- plot_data %>%
    tibble::rownames_to_column("Gene") %>%
    select(Gene, ZSTAT, Specificity)
  
  # ===================== Step 2: Regression analysis =====================
  reg_formula <- y ~ x
  
  if (reg_type == "linear") {
    reg_model <- lm(ZSTAT ~ Specificity, data = reg_data)
    r_squared <- summary(reg_model)$r.squared
    beta <- coef(reg_model)["Specificity"]
    se <- summary(reg_model)$coefficients["Specificity", "Std. Error"]
    p_val <- summary(reg_model)$coefficients["Specificity", "Pr(>|t|)"]
    reg_label <- paste0(
      "Linear Regression\n",
      "y = ", round(beta, 4), "x + ", round(coef(reg_model)[1], 4), "\n",
      "R² = ", round(r_squared, 4), ", P = ", format.pval(p_val, digits = 3)
    )
  } else if (reg_type == "loess") {
    reg_model <- loess(ZSTAT ~ Specificity, data = reg_data, span = 0.75)
    pred_loess <- predict(reg_model)
    r_squared <- cor(reg_data$ZSTAT, pred_loess)^2
    beta <- NA
    se <- NA
    p_val <- NA
    reg_label <- paste0(
      "LOESS Regression (span=0.75)\n",
      "R² ≈ ", round(r_squared, 4)
    )
  } else if (reg_type == "lm_robust") {
    reg_model <- lmrob(ZSTAT ~ Specificity, data = reg_data, 
                       control = lmrob.control(max.it = 1000, tuning.chi = 1.345))
    pred_robust <- predict(reg_model)
    r_squared <- cor(reg_data$ZSTAT, pred_robust)^2
    beta <- coef(reg_model)["Specificity"]
    se <- summary(reg_model)$coefficients["Specificity", "Std. Error"]
    p_val <- summary(reg_model)$coefficients["Specificity", "Pr(>|t|)"]
    reg_label <- paste0(
      "Robust Linear Regression\n",
      "y = ", round(beta, 4), "x + ", round(coef(reg_model)[1], 4), "\n",
      "R² = ", round(r_squared, 4), ", P = ", format.pval(p_val, digits = 3)
    )
  } else {
    stop("unsupported reg_type!")
  }
  
  reg_result_summary <- list(
    type = reg_type,
    formula = as.character(reg_formula),
    r_squared = r_squared,
    coefficients = if (!is.na(beta)) c(intercept = coef(reg_model)[1], slope = beta) else NULL,
    se = se,
    p_value = p_val,
    top_thresholds = list(
      zstat_top_threshold = zstat_top_threshold,
      specificity_top_threshold = specificity_top_threshold,
      top_pct = top_pct
    )
  )
  
  # ===================== Step 3: Plotting =====================
  p <- ggplot(reg_data, aes(x = Specificity, y = ZSTAT))
  
  if (use_density_color) {
    p <- p +
      geom_pointdensity(
        aes(color = after_stat(density)), 
        size = 1.8, alpha = 0.8
      ) +
      scale_color_viridis_c(option = "plasma", name = "Point Density")
  } else {
    p <- p +
      geom_point(
        color = point_color, 
        alpha = point_alpha,
        size = 1.8
      )
  }
  
  if (reg_type == "linear") {
    p <- p +
      geom_smooth(method = "lm", formula = reg_formula, 
                  color = line_color, linewidth = 1.2, alpha = 0.7, 
                  se = TRUE, fill = ci_color)
  } else if (reg_type == "loess") {
    p <- p +
      geom_smooth(method = "loess", formula = reg_formula, span = 0.75,
                  color = line_color, linewidth = 1.2, alpha = 0.7, 
                  se = TRUE, fill = ci_color)
  } else if (reg_type == "lm_robust") {
    robust_pred <- data.frame(
      Specificity = seq(min(reg_data$Specificity), max(reg_data$Specificity), length.out = 100)
    )
    robust_pred$ZSTAT <- predict(reg_model, newdata = robust_pred)
    
    robust_pred_ci <- cbind(
      robust_pred,
      predict(reg_model, newdata = robust_pred, se.fit = TRUE)[c("fit", "se.fit")]
    )
    robust_pred_ci$ymin <- robust_pred_ci$fit - 1.96*robust_pred_ci$se.fit
    robust_pred_ci$ymax <- robust_pred_ci$fit + 1.96*robust_pred_ci$se.fit
    
    p <- p +
      geom_ribbon(
        data = robust_pred_ci,
        aes(x = Specificity, ymin = ymin, ymax = ymax),
        alpha = 0.7, fill = ci_color
      ) +
      geom_line(data = robust_pred, aes(x = Specificity, y = ZSTAT),
                color = line_color, linewidth = 1.2, alpha = 0.7)
  }
  
  if (!is.null(custom_gene_list)) {
    p <- p +
      geom_point(
        data = reg_data %>% filter(Gene %in% custom_gene_list),
        aes(x = Specificity, y = ZSTAT),
        color = "#E64B35", size = 3, shape = 21, fill = "#E64B3580", stroke = 1,
        alpha = 1
      ) +
      geom_text_repel(
        data = reg_data %>% filter(Gene %in% custom_gene_list),
        aes(label = Gene),
        size = 3.8, color = "black", fontface = "bold",
        box.padding = 0.6, point.padding = 0.4,
        max.overlaps = 30, segment.color = "gray50"
      )
  }
  
  if (show_reg_label && reg_type != "loess") {
    p <- p +
      geom_label(
        x = max(reg_data$Specificity)*0.7, 
        y = max(reg_data$ZSTAT)*0.9,
        label = reg_label, 
        hjust = 0, vjust = 1,
        size = 4, color = line_color, fontface = "bold",
        fill = "white", alpha = 0.8, 
        label.size = 0.5
      )
  }
  
  p <- p +
    labs(
      x = paste0("Cell Type Specificity (", cell_type, ")"),
      y = "ZSTAT (MAGMA Enrichment)",
      title = paste0("MAGMA ZSTAT vs. Cell Specificity - Linear Mode"),
      subtitle = paste0(organ, " - ", disease),
      caption = paste0("Total genes: ", nrow(reg_data), 
                       if (!is.null(custom_gene_list)) paste0(" | Custom genes labeled: ", sum(reg_data$Gene %in% custom_gene_list)) else "")
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray50"),
      plot.caption = element_text(size = 9, hjust = 0, color = "gray40"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      legend.position = if (use_density_color) "right" else "none",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      panel.grid.minor = element_blank()
    )
  
  # ===================== Step 4: Save image =====================
  ggsave(
    filename = save_path,
    plot = p,
    width = width, height = height, dpi = 300,
    bg = "white"
  )
  
  # ===================== Step 5: Return results =====================
  message(paste("Image saved to:", save_path))
  message(paste("R² =", round(reg_result_summary$r_squared, 4)))
  message(paste("ZSTAT top ", top_pct*100, "% threshold: ", round(zstat_top_threshold, 4), sep = ""))
  message(paste("Specificity top ", top_pct*100, "% threshold: ", round(specificity_top_threshold, 4), sep = ""))
  
  # Count top gene distribution
  top_gene_counts <- table(plot_data$is_top_gene)
  for (category in names(top_gene_counts)) {
    message(paste("  ", category, ": ", top_gene_counts[category], " genes", sep = ""))
  }
  
  if (reg_type %in% c("linear", "lm_robust")) {
    message(paste("P =", format.pval(reg_result_summary$p_value, digits = 3)))
  }
  
  return(list(
    plot = p,
    plot_data = plot_data %>% select(ZSTAT, Specificity, is_top_gene, is_label),
    reg_result = reg_result_summary
  ))
}

plot_magma_top10 <- function(df, 
                             cell_type, 
                             disease, 
                             organ,
                             top_pct = 0.1,
                             custom_gene_list = NULL,
                             use_density_color = TRUE,
                             point_color = "#3C5488",
                             width = 9, 
                             height = 7,
                             save_path = "./MAGMA_ZSTAT_Specificity_Top10.png") {
  if (!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
  if (!require(dplyr)) {install.packages("dplyr"); library(dplyr)}
  if (!require(ggrepel)) {install.packages("ggrepel"); library(ggrepel)}
  if (use_density_color) {
    if (!require(ggpointdensity)) {install.packages("ggpointdensity"); library(ggpointdensity)}
  }
  
  # ===================== Step 1: Data preprocessing =====================
  # Calculate top percentage threshold (only consider positive values)
  zstat_positive <- df$ZSTAT[df$ZSTAT > 0]
  specificity_positive <- df$Specificity[df$Specificity > 0]
  
  # Check if sufficient positive value data exists
  if (length(zstat_positive) == 0) {
    stop("ZSTAT has no positive value data, cannot calculate top percentage threshold")
  }
  if (length(specificity_positive) == 0) {
    stop("Specificity has no positive value data, cannot calculate top percentage threshold")
  }
  
  zstat_top_threshold <- quantile(zstat_positive, probs = 1 - top_pct, na.rm = TRUE)
  specificity_top_threshold <- quantile(specificity_positive, probs = 1 - top_pct, na.rm = TRUE)
  
  plot_data <- df %>%
    tibble::rownames_to_column("Gene") %>%
    filter(!is.na(ZSTAT) & !is.na(Specificity)) %>%
    mutate(
      # column 1: is_top_gene
      is_zstat_top = ifelse(ZSTAT >= zstat_top_threshold, TRUE, FALSE),
      is_specificity_top = ifelse(Specificity >= specificity_top_threshold, TRUE, FALSE),
      is_top_gene = case_when(
        is_zstat_top & is_specificity_top ~ "Both",
        is_zstat_top & !is_specificity_top ~ "ZSTAT_top",
        !is_zstat_top & is_specificity_top ~ "Specificity_top",
        TRUE ~ "Neither"
      ),
      # column 2: is_label
      is_label = ifelse(!is.null(custom_gene_list) & Gene %in% custom_gene_list, "Labeled", "Unlabeled"),
      is_top_simple = case_when(
        !is.null(custom_gene_list) & Gene %in% custom_gene_list ~ "Top",
        is.null(custom_gene_list) & is_zstat_top & is_specificity_top ~ "Top",
        TRUE ~ "Other"
      )
    ) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("Gene")
  
  # ===================== Step 2: Plotting =====================
  p <- ggplot(plot_data %>% tibble::rownames_to_column("Gene"), 
              aes(x = Specificity, y = ZSTAT))
  
  if (use_density_color) {
    p <- p +
      geom_pointdensity(
        aes(color = after_stat(density)), 
        size = 1.8, alpha = 0.8
      ) +
      scale_color_viridis_c(option = "plasma", name = "Point Density")
  } else {
    p <- p +
      geom_point(
        color = point_color, 
        size = 1.8, alpha = 0.8
      )
  }
  
  # Depending on whether there is a custom gene list, choose different labeling methods
  if (!is.null(custom_gene_list)) {
    # With custom gene list: only mark custom genes
    top_genes_to_plot <- plot_data %>% 
      tibble::rownames_to_column("Gene") %>% 
      filter(is_label == "Labeled")
  } else {
    # No custom gene list: mark all top genes (both ZSTAT and Specificity top)
    top_genes_to_plot <- plot_data %>% 
      tibble::rownames_to_column("Gene") %>% 
      filter(is_top_gene == "Both")
  }
  
  # Plot highlighted points
  if (nrow(top_genes_to_plot) > 0) {
    p <- p +
      geom_point(
        data = top_genes_to_plot,
        aes(x = Specificity, y = ZSTAT),
        color = "#E64B35", size = 3, shape = 21, fill = "#E64B3580", stroke = 1
      )
  }
  
  # Always draw topN% threshold line (regardless of whether there is a custom gene list)
  p <- p +
    # geom_hline(yintercept = zstat_top_threshold, 
    #            linetype = "dashed", color = "#E64B35", linewidth = 0.8) +
    geom_vline(xintercept = specificity_top_threshold, 
               linetype = "dashed", color = "#E64B35", linewidth = 0.8)
  
  if (nrow(top_genes_to_plot) > 0) {
    p <- p +
      geom_text_repel(
        data = top_genes_to_plot,
        aes(label = Gene),
        size = 3.8, color = "black", fontface = "bold",
        box.padding = 0.6, point.padding = 0.4,
        max.overlaps = 30, segment.color = "gray50"
      )
  }
  
  caption_text <- ifelse(
    is.null(custom_gene_list),
    paste0("ZSTAT top ", top_pct*100, "% threshold (positive values): ", round(zstat_top_threshold, 4), 
           "\nSpecificity top ", top_pct*100, "% threshold (positive values): ", round(specificity_top_threshold, 4)),
    paste0("Custom genes labeled: ", nrow(top_genes_to_plot), "/", length(custom_gene_list))
  )
  
  p <- p +
    labs(
      x = paste0("Cell Type Specificity (", cell_type, ")"),
      y = "ZSTAT (MAGMA Enrichment)",
      title = paste0("MAGMA Enrichment ZSTAT vs. Cell Specificity - Top 10% Mode"),
      subtitle = paste0(organ, " - ", disease),
      caption = caption_text
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray50"),
      plot.caption = element_text(size = 9, hjust = 0, color = "gray40"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      legend.position = if (use_density_color) "right" else "none",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9)
    )
  
  # ===================== Step 3: Save image =====================
  ggsave(
    filename = save_path,
    plot = p,
    width = width, height = height, dpi = 300,
    bg = "white"
  )
  
  # ===================== Step 4: Return results =====================
  # Organize plot_data (keep core columns)
  plot_data_final <- plot_data %>%
    select(ZSTAT, Specificity, is_top_gene, is_label)
  
  # Count top gene distribution
  top_gene_counts <- table(plot_data_final$is_top_gene)
  
  # Message output
  message(paste("Image saved to:", save_path))
  message(paste("ZSTAT top ", top_pct*100, "% threshold (positive values): ", round(zstat_top_threshold, 4), sep = ""))
  message(paste("Specificity top ", top_pct*100, "% threshold (positive values): ", round(specificity_top_threshold, 4), sep = ""))
  
  for (category in names(top_gene_counts)) {
    message(paste("  ", category, ": ", top_gene_counts[category], " genes", sep = ""))
  }
  
  # Return ggplot object + processed data frame
  return(list(
    plot = p,
    plot_data = plot_data_final,
    thresholds = list(
      zstat_top_threshold = zstat_top_threshold,
      specificity_top_threshold = specificity_top_threshold,
      top_pct = top_pct
    )
  ))
}
