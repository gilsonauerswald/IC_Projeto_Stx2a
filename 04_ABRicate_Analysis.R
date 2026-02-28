#' @title Pangenomic Profile Analysis: Virulome & Resistome
#' @description This script processes ABRicate output tables to analyze 
#' Virulence factors, Antimicrobial Resistance (AMR) genes, and Plasmids.
#' It performs data cleaning, normalization, Spearman correlation analysis,
#' and generates a unified multi-panel figure for publication.
#' 
#' @details 
#' Input requirements:
#' 1. Virulome summary (.tabular)
#' 2. Resistome/Plasmidome summary (.tabular)
#' 3. Metadata (.tsv) with columns 'ABRicate' and 'Nome_Real'
#'
#' @author [Seu Nome] (ORCID: 0000-XXXX-XXXX-XXXX)
#' @date 2026-02-17
#' @version 1.1.5 (Fix Overlapping Points)
#' @license MIT License
#' @doi 10.5281/zenodo.XXXXX (Placeholder)
#'
#' @dependencies
#' - R version 4.5.2
#' - Core Packages: tidyverse, ggplot2, patchwork, openxlsx, viridis, ggrepel

# ==============================================================================
# 1. SETUP AND LIBRARIES
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(openxlsx)
  library(viridis)
  library(ggrepel)
  library(scales)
})

# ==============================================================================
# 2. DATA PROCESSING FUNCTIONS
# ==============================================================================

#' Load Metadata
#' @description Loads and validates the metadata file.
#' @param file_path Path to the TSV file.
#' @return A tibble with metadata.
load_metadata <- function(file_path) {
  if (!file.exists(file_path)) stop("Metadata file not found: ", file_path)
  
  meta <- readr::read_tsv(file_path, show_col_types = FALSE)
  
  required_cols <- c("ABRicate", "Nome_Real")
  if (!all(required_cols %in% names(meta))) {
    stop("Metadata must contain columns: ", paste(required_cols, collapse = ", "))
  }
  return(meta)
}

#' Map Sample Names
#' @description Translates filename IDs to real sample names.
#' @param file_names Vector of filenames from ABRicate output.
#' @param metadata Metadata tibble.
#' @return Vector of mapped names.
map_sample_names <- function(file_names, metadata) {
  lookup_map <- setNames(metadata$Nome_Real, metadata$ABRicate)
  new_names <- lookup_map[file_names]
  # Keep original if no match found
  new_names[is.na(new_names)] <- file_names[is.na(new_names)]
  return(new_names)
}

#' Process ABRicate Data
#' @description Reads raw Galaxy/ABRicate tabular output and normalizes format.
#' @param file_path Path to .tabular file.
#' @param type String describing the data type (e.g., "Virulome").
#' @return Tidy tibble with gene presence/absence.
process_abricate_data <- function(file_path, type) {
  if (!file.exists(file_path)) stop("Input file not found: ", file_path)
  
  raw_data <- readr::read_tsv(file_path, show_col_types = FALSE, comment = "")
  # Rename first column generic 'FILE' for consistency
  names(raw_data)[1] <- "FILE"
  
  clean_data <- raw_data %>%
    dplyr::select(-dplyr::any_of("NUM_FOUND")) %>%
    dplyr::mutate(dplyr::across(-FILE, as.character)) %>%
    tidyr::pivot_longer(cols = -FILE, names_to = "GENE", values_to = "IDENTITY") %>%
    dplyr::mutate(
      IDENTITY = suppressWarnings(as.numeric(IDENTITY)),
      # Binary presence: Identity >= 95%
      IS_PRESENT = ifelse(!is.na(IDENTITY) & IDENTITY >= 95, 1, 0),
      DATA_TYPE = type
    )
  
  return(clean_data)
}

#' Save Excel Table
#' @description Exports processed data to Excel.
save_excel_table <- function(data, file_name, dir_out) {
  path <- file.path(dir_out, paste0(file_name, ".xlsx"))
  openxlsx::write.xlsx(list(Data = data), path)
}

# ==============================================================================
# 3. VISUALIZATION FUNCTIONS
# ==============================================================================

#' Generate Scientific Theme
theme_nature <- function() {
  theme_bw(base_size = 12) +
    theme(
      text = element_text(family = "sans"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(color = "black"),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      legend.position = "bottom"
    )
}

#' Plot Virulome Heatmap
plot_virulome <- function(data, metadata) {
  target_genes <- c("stx1A", "stx1B", "stx2A", "stx2B", "eae", "aggR", "ehxA")
  
  plot_data <- data %>%
    dplyr::filter(GENE %in% target_genes) %>%
    dplyr::mutate(
      SAMPLE_ID = map_sample_names(FILE, metadata),
      IDENTITY = tidyr::replace_na(IDENTITY, 0)
    )
  
  ggplot(plot_data, aes(x = GENE, y = SAMPLE_ID, fill = IDENTITY)) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(IDENTITY > 0, round(IDENTITY, 0), "")), size = 3) +
    scale_fill_gradientn(colors = c("grey95", "#fee0d2", "#de2d26"), limits = c(0, 100)) +
    labs(title = "Virulence Gene Profile", x = "Gene Target", y = NULL, fill = "Identity %") +
    theme_nature() +
    theme(axis.text.x = element_text(face = "italic"))
}

#' Plot Resistome Burden
plot_resistome_burden <- function(data, metadata) {
  plot_data <- data %>%
    dplyr::group_by(FILE) %>%
    dplyr::summarise(Count = sum(IS_PRESENT, na.rm = TRUE)) %>%
    dplyr::mutate(SAMPLE_ID = map_sample_names(FILE, metadata))
  
  ggplot(plot_data, aes(x = reorder(SAMPLE_ID, -Count), y = Count)) +
    geom_bar(stat = "identity", fill = "#2c7bb6", color = "black", width = 0.7) +
    geom_text(aes(label = Count), vjust = -0.5, size = 3.5) +
    labs(title = "Antimicrobial Resistance Burden", x = NULL, y = "Gene Count (n)") +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#' Plot Plasmid Profile
plot_plasmidome <- function(data, metadata) {
  plot_data <- data %>%
    dplyr::filter(IS_PRESENT == 1) %>%
    dplyr::mutate(
      SAMPLE_ID = map_sample_names(FILE, metadata),
      # Extract replicon family (regex)
      Replicon_Family = stringr::str_extract(GENE, "^[A-Za-z0-9]+")
    ) 
  
  ggplot(plot_data, aes(x = SAMPLE_ID, fill = Replicon_Family)) +
    geom_bar(position = "stack", color = "black", linewidth = 0.2) +
    scale_fill_viridis_d(option = "turbo") +
    labs(title = "Plasmid Replicon Content", x = NULL, y = "Count", fill = "Family") +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#' Plot Correlation (Resistome vs Plasmidome)
#' @description Updates: Uses JITTER to prevent point overlap. Full Join for safety.
plot_correlation <- function(res_data, plasmid_data, metadata) {
  # Aggregation
  r_sum <- res_data %>% dplyr::group_by(FILE) %>% dplyr::summarise(R = sum(IS_PRESENT, na.rm=TRUE))
  p_sum <- plasmid_data %>% dplyr::group_by(FILE) %>% dplyr::summarise(P = sum(IS_PRESENT, na.rm=TRUE))
  
  # Join logic update: Full join to ensure strains with 0 hits in one dataset are kept
  combined <- dplyr::full_join(r_sum, p_sum, by = "FILE") %>%
    dplyr::mutate(
      R = tidyr::replace_na(R, 0),
      P = tidyr::replace_na(P, 0),
      SAMPLE_ID = map_sample_names(FILE, metadata)
    )
  
  # Variance Check
  if(sd(combined$R) == 0 || sd(combined$P) == 0) {
    warning("Zero variance detected. Correlation plot might look flat.")
  }
  
  # Statistics
  cor_test <- cor.test(combined$R, combined$P, method = "spearman", exact = FALSE)
  rho_val <- round(cor_test$estimate, 3)
  p_val   <- format.pval(cor_test$p.value, digits = 3, eps = 0.001)
  
  # PLOT UPDATE: Jitter added to separate overlapping points
  ggplot(combined, aes(x = P, y = R)) +
    # Trend line
    geom_smooth(method = "lm", formula = y ~ x, color = "grey50", fill = "grey90", alpha = 0.5, linetype = "dashed") +
    
    # Points: Jitter added here!
    geom_point(aes(fill = SAMPLE_ID), 
               size = 5, 
               shape = 21, 
               color = "black", 
               stroke = 0.8,
               position = position_jitter(width = 0.3, height = 0.3), # ESPALHA OS PONTOS
               alpha = 0.9) +
    
    # Stats Annotation
    annotate("text", x = min(combined$P), y = max(combined$R), 
             label = paste0("Spearman rho = ", rho_val, "\np-value = ", p_val), 
             hjust = 0, vjust = 1, size = 4, fontface = "italic") +
    
    # Colors & Legend
    scale_fill_viridis_d(option = "turbo", name = "Strain Reference") +
    
    labs(title = "Plasmids vs. Resistome", 
         x = "Plasmid Replicons (n)", y = "AMR Genes (n)") +
    
    theme_nature() +
    theme(legend.position = "right")
}

# ==============================================================================
# 4. MAIN EXECUTION BLOCK
# ==============================================================================

# Configuration
CONFIG <- list(
  FILE_VIR   = "Galaxy41-[ABRicate Summary on dataset 11-20].tabular",
  FILE_RES   = "Galaxy42-[ABRicate Summary on dataset 21-30].tabular",
  FILE_PLASM = "Galaxy42-[ABRicate Summary on dataset 21-30].tabular", # Usually same file as resistome in this pipeline
  FILE_META  = "Metadados.tsv",
  DIR_OUT    = "outputs"
)

main <- function() {
  
  # 1. Environment Setup
  current_wd <- getwd()
  cat("\n[INFO] Environment: R v4.5.2\n")
  cat("[INFO] Working Directory:", current_wd, "\n")
  
  dir_output <- file.path(current_wd, CONFIG$DIR_OUT)
  if (!dir.exists(dir_output)) {
    dir.create(dir_output, recursive = TRUE)
    cat("[INFO] Output directory created:", dir_output, "\n")
  }
  
  # 2. Load Data
  cat("[1/4] Loading and Processing Data...\n")
  meta <- load_metadata(file.path(current_wd, CONFIG$FILE_META))
  
  data_vir   <- process_abricate_data(file.path(current_wd, CONFIG$FILE_VIR), "Virulome")
  data_res   <- process_abricate_data(file.path(current_wd, CONFIG$FILE_RES), "Resistome")
  data_plasm <- process_abricate_data(file.path(current_wd, CONFIG$FILE_PLASM), "Plasmidome")
  
  # 3. Generate Individual Plots
  cat("[2/4] Generating Figures...\n")
  p1 <- plot_virulome(data_vir, meta)
  p2 <- plot_resistome_burden(data_res, meta)
  p3 <- plot_plasmidome(data_plasm, meta)
  p4 <- plot_correlation(data_res, data_plasm, meta)
  
  # 4. Compose Unified Figure
  cat("[3/4] Composing Panel...\n")
  layout_design <- (p1 | p2) / (p3 | p4) +
    patchwork::plot_annotation(
      tag_levels = 'A',
      title = "Integrated Genomic Characterization: Virulome, Resistome and Mobilome",
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )
  
  # Save
  ggsave(
    filename = file.path(dir_output, "Figure_Unified_Genomic_Profile.png"),
    plot = layout_design,
    width = 18, height = 14, dpi = 300, bg = "white"
  )
  
  # 5. Export Data Tables
  cat("[4/4] Exporting Tables...\n")
  
  # Prepare for export (Add real names)
  export_vir <- data_vir %>% dplyr::mutate(Strain = map_sample_names(FILE, meta)) %>% dplyr::relocate(Strain)
  export_res <- data_res %>% dplyr::mutate(Strain = map_sample_names(FILE, meta)) %>% dplyr::relocate(Strain)
  
  save_excel_table(export_vir, "Table_Virulome_Processed", dir_output)
  save_excel_table(export_res, "Table_Resistome_Processed", dir_output)
  
  # Save Session Info
  writeLines(c(
    paste("Pipeline Run Date:", Sys.time()),
    paste("R Version:", R.version.string),
    "Bioconductor: v3.22 (Expected)",
    capture.output(sessionInfo())
  ), file.path(dir_output, "session_info_ABRicate.txt"))
  
  cat("\n✓ PIPELINE COMPLETED SUCCESSFULLY.\n")
  cat("✓ Results saved in:", dir_output, "\n")
}

# Execute Main
if (!interactive()) {
  main()
} else {
  cat("Script loaded. Ensure input files are in root directory and run 'main()'.\n")
}