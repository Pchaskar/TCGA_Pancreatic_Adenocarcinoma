# HeatMap of DEGs

# Set options for future globals and memory limit
options(future.globals.maxSize = 10000 * 1024 ^ 2)  # Set maximum size for future globals
library(BiocParallel)

# Register parallel backend
register(MulticoreParam(12))  # Use 12 cores

# Path to the scripts
library("rstudioapi")

# the following line is for getting the path of your current open file
script_path <- getActiveDocumentContext()$path

dirpath <- dirname(script_path)

# The next line set the working directory to the relevant one:
setwd(dirname(dirpath))

# tidyverse core packages
library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(purrr)
library(nclust)

source("./scripts/Marker.R")

# Restore the object
TCGA_PAAD <- readRDS("../00_RAW_DATA/RDS/TCGA_PAAD.rds")

# Get Cell annotations from scRNAseq
meta_info_cells <- readRDS("../../Prasad_PublicData_scRNAseq_PDAC_2024/02_ANALYZED_DATA/RDS/meta_info_cells.rds")

# filtered counts
vst_counts <- TCGA_PAAD$vst_counts # If it's a dense matrix

# Create colData and rowData
meta_data <- as.data.frame(TCGA_PAAD$meta_data)

# update the sample type
# Recode sample types
meta_data$sample_type <- recode(meta_data$sample_type,
                                "Primary Tumor" = "Tumor",
                                "Solid Tissue Normal" = "Normal")

# remove TCGA_PAAD
rm(TCGA_PAAD)

##--
# All counts
# Subset the matrix based on matching row names
deg_counts <- vst_counts
dim(deg_counts)

# Clustering
scaled_counts <- t(scale(t(deg_counts)))

# Add cell annotations
gene_meta_data <- as.data.frame(scaled_counts[,1])
colnames(gene_meta_data) <- "TMP"

# Ensure the row names are preserved as a column for matching
meta_info_cells <- meta_info_cells %>% rownames_to_column(var = "Gene")
gene_meta_data <- gene_meta_data %>% rownames_to_column(var = "Gene")

# Join the data frames based on the 'Gene' column
updated_gene_meta_data <- gene_meta_data %>%
  left_join(meta_info_cells %>% select(Gene, Cell), by = "Gene")

# Optionally, restore row names
updated_gene_meta_data <- updated_gene_meta_data %>% column_to_rownames(var = "Gene")

hist <- coldmap(scaled_counts, method = "ward")

saveRDS(hist, "./RDS/hist_TCGA.rds")

coldmap(
  scaled_counts,
  clust = hist,
  saturation = TRUE,
  ctag = make_tag(
    meta_data,
    varnames = c(
      "tissue_type",
      "ajcc_pathologic_stage",
      "tissue_or_organ_of_origin",
      "primary_diagnosis"
    ),
    cols = c("green4", "brown", "purple", "black")
  ),
  ctag.space = 8,
  rmarg = 1,
  rtag = make_tag(
    updated_gene_meta_data,
    varnames = c("Cell"),
    cols = c("black")
  ),
  rtag.space = 2,
  rlab = list(
    list(lit_panc = lit_panc),
    list(basal = basal),
    list(classical = classical),
    list(emt = emt),
    list(neuro_endo = neuro_endo)
  ),
  col.rlab  = col.blca,
  rdend.col = list(
list(c(1, 7000), "brown"),
list(c(9000, 10000), "purple"),
list(c(11900, 13100), "violet"),
list(c(14000, 15000), "green")
)
)

