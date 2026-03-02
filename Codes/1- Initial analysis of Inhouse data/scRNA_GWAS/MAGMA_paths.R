library(MAGMA.Celltyping)


dir_software <- "D:/OneDrive/SciLifeLab/Materials/Softwares/MAGMA"
dir_working <- "./MAGMA_files/"

# Tutorial: https://neurogenomics.github.io/MAGMA_Celltyping/articles/MAGMA.Celltyping.html




###### Heart
Trait_category <- "Heart"

files_Heart <- read.table("Prepare_GWAS_Heart.txt", sep="\t", header = TRUE)

magma_paths_Heart <- list()
for (i in seq(nrow(files_Heart))) {
  Trait_name_full <- files_Heart[i,2]
  Trait_name <- files_Heart[i,5]
  Trait_N<-files_Heart[i,6]
  Trait_population<-files_Heart[i,7]
  
  
  Trait_path <- paste0(dir_working, Trait_category, '/', Trait_name_full)
  
  magma_paths_Heart[[Trait_name]] <- MAGMA.Celltyping::map_snps_to_genes(
    path_formatted = Trait_path,
    genome_build = 'GRCh38', storage_dir = dir_software, N=Trait_N, population=Trait_population)
  
  
  file.rename(paste0(magma_paths_Heart[[Trait_name]], ".txt"), magma_paths_Heart[[Trait_name]])
}

save(magma_paths_Heart, file="MAGMA_files/20251121Batch01to08_MAGMA_paths_Heart.RData")




###### Brain
Trait_category <- "Brain"

files_Brain <- read.table("Prepare_GWAS_Brain.txt", sep="\t", header = TRUE, quote="")

magma_paths_Brain <- list()
for (i in seq(nrow(files_Brain))) {
  Trait_name_full <- gsub(".tsv", "", files_Brain[i,2])
  Trait_name <- files_Brain[i,5]
  Trait_N<-files_Brain[i,6]
  Trait_population<-files_Brain[i,7]
  
  
  Trait_path <- paste0(dir_working, Trait_category, '/', Trait_name_full, ".tsv")
  
  magma_paths_Brain[[Trait_name]] <- MAGMA.Celltyping::map_snps_to_genes(
    path_formatted = Trait_path,
    genome_build = 'GRCh38', storage_dir = dir_software, N=Trait_N, population=Trait_population)
  
  
  file.rename(paste0(magma_paths_Brain[[Trait_name]], ".txt"), magma_paths_Brain[[Trait_name]])
}

save(magma_paths_Brain, file="MAGMA_files/20251121Batch01to08_MAGMA_paths_Brain.RData")




###### Lung
Trait_category <- "Lung"

files_Lung <- read.table("Prepare_GWAS_Lung.txt", sep="\t", header = TRUE, quote="")

magma_paths_Lung <- list()
for (i in seq(nrow(files_Lung))) {
  Trait_name_full <- gsub(".tsv", "", files_Lung[i,2])
  Trait_name <- files_Lung[i,5]
  Trait_N<-files_Lung[i,6]
  Trait_population<-files_Lung[i,7]
  
  
  Trait_path <- paste0(dir_working, Trait_category, '/', Trait_name_full, ".tsv")
  
  magma_paths_Lung[[Trait_name]] <- MAGMA.Celltyping::map_snps_to_genes(
    path_formatted = Trait_path,
    genome_build = 'GRCh38', storage_dir = dir_software, N=Trait_N, population=Trait_population)
  
  
  file.rename(paste0(magma_paths_Lung[[Trait_name]], ".txt"), magma_paths_Lung[[Trait_name]])
}

save(magma_paths_Lung, file="MAGMA_files/20251121Batch01to08_MAGMA_paths_Lung.RData")




###### Liver
Trait_category <- "Liver"

files_Liver <- read.table("Prepare_GWAS_Liver.txt", sep="\t", header = TRUE, quote="")

magma_paths_Liver <- list()
for (i in seq(nrow(files_Liver))) {
  Trait_name_full <- gsub(".tsv", "", files_Liver[i,2])
  Trait_name <- files_Liver[i,5]
  Trait_N<-files_Liver[i,6]
  Trait_population<-files_Liver[i,7]
  
  
  Trait_path <- paste0(dir_working, Trait_category, '/', Trait_name_full, ".tsv")
  
  magma_paths_Liver[[Trait_name]] <- MAGMA.Celltyping::map_snps_to_genes(
    path_formatted = Trait_path,
    genome_build = 'GRCh38', storage_dir = dir_software, N=Trait_N, population=Trait_population)
  
  
  file.rename(paste0(magma_paths_Liver[[Trait_name]], ".txt"), magma_paths_Liver[[Trait_name]])
}

save(magma_paths_Liver, file="MAGMA_files/20251121Batch01to08_MAGMA_paths_Liver.RData")




###### Kidney
Trait_category <- "Kidney"

files_Kidney <- read.table("Prepare_GWAS_Kidney.txt", sep="\t", header = TRUE, quote="")

magma_paths_Kidney <- list()
for (i in seq(nrow(files_Kidney))) {
  Trait_name_full <- gsub(".tsv", "", files_Kidney[i,2])
  Trait_name <- files_Kidney[i,5]
  Trait_N<-files_Kidney[i,6]
  Trait_population<-files_Kidney[i,7]
  
  
  Trait_path <- paste0(dir_working, Trait_category, '/', Trait_name_full, ".tsv")
  
  magma_paths_Kidney[[Trait_name]] <- MAGMA.Celltyping::map_snps_to_genes(
    path_formatted = Trait_path,
    genome_build = 'GRCh38', storage_dir = dir_software, N=Trait_N, population=Trait_population)
  
  
  file.rename(paste0(magma_paths_Kidney[[Trait_name]], ".txt"), magma_paths_Kidney[[Trait_name]])
}

save(magma_paths_Kidney, file="MAGMA_files/20251121Batch01to08_MAGMA_paths_Kidney.RData")




###### Eye
Trait_category <- "Eye"

files_Eye <- read.table("Prepare_GWAS_Eye.txt", sep="\t", header = TRUE, quote="")

magma_paths_Eye <- list()
for (i in seq(nrow(files_Eye))) {
  Trait_name_full <- gsub(".tsv", "", files_Eye[i,2])
  Trait_name <- files_Eye[i,5]
  Trait_N<-files_Eye[i,6]
  Trait_population<-files_Eye[i,7]
  
  
  Trait_path <- paste0(dir_working, Trait_category, '/', Trait_name_full, ".tsv")
  
  magma_paths_Eye[[Trait_name]] <- MAGMA.Celltyping::map_snps_to_genes(
    path_formatted = Trait_path,
    genome_build = 'GRCh38', storage_dir = dir_software, N=Trait_N, population=Trait_population)
  
  
  file.rename(paste0(magma_paths_Eye[[Trait_name]], ".txt"), magma_paths_Eye[[Trait_name]])
}

save(magma_paths_Eye, file="MAGMA_files/20251121Batch01to08_MAGMA_paths_Eye.RData")




###### Systemic
Trait_category <- "Systemic"

files_Systemic <- read.table("Prepare_GWAS_Systemic.txt", sep="\t", header = TRUE, quote="")

magma_paths_Systemic <- list()
for (i in seq(nrow(files_Systemic))) {
  Trait_name_full <- gsub(".tsv", "", files_Systemic[i,2])
  Trait_name <- files_Systemic[i,5]
  Trait_N<-files_Systemic[i,6]
  Trait_population<-files_Systemic[i,7]
  
  
  Trait_path <- paste0(dir_working, Trait_category, '/', Trait_name_full, ".tsv")
  
  magma_paths_Systemic[[Trait_name]] <- MAGMA.Celltyping::map_snps_to_genes(
    path_formatted = Trait_path,
    genome_build = 'GRCh38', storage_dir = dir_software, N=Trait_N, population=Trait_population)
  
  
  file.rename(paste0(magma_paths_Systemic[[Trait_name]], ".txt"), magma_paths_Systemic[[Trait_name]])
}

for (i in seq(nrow(files_Systemic))) {
  Trait_name_full <- files_Systemic[i,2]
  Trait_name <- files_Systemic[i,5]
  magma_paths_Systemic[[Trait_name]] <- 
    paste0("D:\\OneDrive\\SciLifeLab\\Endothelial_cells\\20251121MAGMA\\MAGMA_files\\", 
           Trait_category, "\\MAGMA_Files\\", Trait_name_full, ".35UP.10DOWN\\", Trait_name_full, ".35UP.10DOWN.genes.out")
}
save(magma_paths_Systemic, file="MAGMA_files/20251121Batch01to08_MAGMA_paths_Systemic.RData")




###### Neuropsychiatric
Trait_category <- "Neuropsychiatric"

files_Neuropsychiatric <- read.table("Prepare_GWAS_Neuropsychiatric.txt", sep="\t", header = TRUE, quote="")

magma_paths_Neuropsychiatric <- list()
for (i in seq(nrow(files_Neuropsychiatric))) {
  Trait_name_full <- gsub(".tsv", "", files_Neuropsychiatric[i,2])
  Trait_name <- files_Neuropsychiatric[i,5]
  Trait_N<-files_Neuropsychiatric[i,6]
  Trait_population<-files_Neuropsychiatric[i,7]
  
  
  Trait_path <- paste0(dir_working, Trait_category, '/', Trait_name_full, ".tsv")
  
  magma_paths_Neuropsychiatric[[Trait_name]] <- MAGMA.Celltyping::map_snps_to_genes(
    path_formatted = Trait_path,
    genome_build = 'GRCh38', storage_dir = dir_software, N=Trait_N, population=Trait_population)
  
  
  file.rename(paste0(magma_paths_Neuropsychiatric[[Trait_name]], ".txt"), magma_paths_Neuropsychiatric[[Trait_name]])
}

for (i in seq(nrow(files_Neuropsychiatric))) {
  Trait_name_full <- files_Neuropsychiatric[i,2]
  Trait_name <- files_Neuropsychiatric[i,5]
  magma_paths_Neuropsychiatric[[Trait_name]] <- 
    paste0("D:\\OneDrive\\SciLifeLab\\Endothelial_cells\\20251121MAGMA\\MAGMA_files\\", 
           Trait_category, "\\MAGMA_Files\\", Trait_name_full, ".35UP.10DOWN\\", Trait_name_full, ".35UP.10DOWN.genes.out")
}

save(magma_paths_Neuropsychiatric, file="MAGMA_files/20251121Batch01to08_MAGMA_paths_Neuropsychiatric.RData")

