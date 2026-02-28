#' @title Genomic Completeness Analysis Pipeline: BUSCO
#' @description This script processes BUSCO tabular reports to assess genomic 
#' completeness. It integrates metadata, performs comparative statistical 
#' analysis (Pathogenic vs. Commensal), and generates publication-ready 
#' visualizations (Barplots, Boxplots, Heatmaps).
#' 
#' @details 
#' Execution requires input files (Galaxy *.tabular and Metadados.tsv) 
#' to be located in the same directory as this script.
#'
#' @author [Seu Nome] (ORCID: 0000-XXXX-XXXX-XXXX)
#' @date 2026-02-17
#' @version 1.2.1 (Wilcoxon Ties Fix)
#' @license MIT License
#' @doi 10.5281/zenodo.XXXXX (Placeholder)
#'
#' @dependencies
#' - R version 4.5.2
#' - Bioconductor version 3.22
#' - Core Packages: tidyverse (>= 2.0.0), ggplot2, patchwork, openxlsx, viridis, ggpubr

# ==============================================================================
# 1. SETUP AND LIBRARIES
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(openxlsx)
  library(viridis)
  library(ggpubr)
  library(scales)
})

# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

#' Parse Single BUSCO File
#' @description Reads a raw Galaxy/BUSCO tabular file, calculates percentages,
#' and attempts to map the dataset ID to a real sample name using metadata.
#' 
#' @param file_path String. Path to the .tabular file.
#' @param metadata Tibble. The loaded metadata table.
#' @return A tibble with summarized BUSCO stats for one sample.
parse_busco_file <- function(file_path, metadata) {
  
  # 1. Extract Dataset ID from filename (e.g., "dataset 1")
  # Uses regex to find the number after "dataset"
  dataset_num <- stringr::str_extract(basename(file_path), "dataset\\s+(\\d+)", group = 1)
  
  if (is.na(dataset_num)) {
    warning("Could not extract dataset number from filename: ", basename(file_path))
    dataset_num <- "Unknown"
  }
  
  # 2. Define BUSCO columns specs
  cols_spec <- c("Busco_id", "Status", "Sequence", "Gene_Start", "Gene_End", 
                 "Strand", "Score", "Length", "OrthoDB_url", "Description")
  
  # 3. Read Data (Skip comments)
  raw_data <- readr::read_tsv(file_path, comment = "#", col_names = cols_spec, 
                              show_col_types = FALSE)
  
  # 4. Calculate Stats
  summary_stats <- raw_data %>%
    dplyr::count(Status) %>%
    tidyr::complete(Status = c("Complete", "Duplicated", "Fragmented", "Missing"), 
                    fill = list(n = 0)) %>%
    dplyr::mutate(Percentage = n / sum(n) * 100)
  
  # 5. Metadata Mapping Logic (Heuristic matching)
  # Construct potential keys found in the "BUSCO" column of metadata
  keys <- c(
    paste0("Busco on dataset ", dataset_num, " Full table - Specific lineage"),
    paste0("Busco on dataset ", dataset_num, "_ Short summary - Specific lineage"),
    paste0("Busco on dataset ", dataset_num)
  )
  
  # Attempt match
  matched_name <- metadata$Nome_Real[metadata$BUSCO %in% keys]
  
  # Fallback if no match found
  final_name <- if (length(matched_name) > 0) matched_name[1] else paste0("Dataset_", dataset_num)
  
  return(summary_stats %>% dplyr::mutate(Sample = final_name))
}

#' Perform Statistical Analysis
#' @description Compares 'Pathogenic' vs 'Commensal' groups using T-test or Wilcoxon.
#' @param data Wide format tibble containing 'Group' and 'Total_Complete'.
#' @return A list containing the test results and a tibble for reporting.
run_statistics <- function(data) {
  
  # Define groups
  pathogenic <- data %>% dplyr::filter(Group == "Pathogenic") %>% dplyr::pull(Total_Complete)
  reference  <- data %>% dplyr::filter(Group == "Commensal (Ref)") %>% dplyr::pull(Total_Complete)
  
  # Check data sufficiency (Need 1 ref and >= 2 pathogens)
  if (length(reference) != 1 || length(pathogenic) < 2) {
    warning("Insufficient data for statistics (Need 1 Reference and >2 Pathogens).")
    return(NULL)
  }
  
  # Normality Test (Shapiro-Wilk)
  shapiro <- shapiro.test(pathogenic)
  is_normal <- shapiro$p.value > 0.05
  
  # Hypothesis Test (One-sample vs Reference Value)
  if (is_normal) {
    test_res <- t.test(pathogenic, mu = reference)
    method <- "One-sample t-test"
    stat_str <- sprintf("t=%.2f", test_res$statistic)
  } else {
    # FIX: Added 'exact = FALSE' to handle ties and zeroes (identical values to ref)
    # This suppresses warnings and uses the normal approximation, which is standard.
    test_res <- wilcox.test(pathogenic, mu = reference, exact = FALSE)
    method <- "One-sample Wilcoxon test"
    stat_str <- sprintf("V=%.1f", test_res$statistic)
  }
  
  # Result Tibble
  stats_table <- tibble::tibble(
    Analysis = "Genomic Completeness Comparison",
    Group_Tested = paste0("Pathogenic (n=", length(pathogenic), ")"),
    Reference_Strain = "K-12 MG1655",
    Reference_Value = reference,
    Pathogenic_Mean = mean(pathogenic),
    Normality_P = shapiro$p.value,
    Test_Method = method,
    Test_Statistic = stat_str,
    P_Value = test_res$p.value,
    Significance = dplyr::if_else(test_res$p.value < 0.05, "Significant (*)", "ns")
  )
  
  return(list(test_obj = test_res, table = stats_table, method = method, ref_val = reference))
}

#' Generate Scientific Theme
#' @description Returns a standardized ggplot2 theme.
theme_nature <- function() {
  theme_minimal(base_family = "sans", base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
      axis.title = element_text(face = "bold", size = 10),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    )
}

# ==============================================================================
# 3. MAIN EXECUTION BLOCK (RELATIVE PATHS)
# ==============================================================================

# Configuration: File names expected in the root directory
CONFIG <- list(
  FILE_META_NAME = "Metadados.tsv",
  PATTERN_BUSCO  = "Galaxy.*Busco.*Full table.*\\.tabular$", 
  DIR_OUT        = "outputs"
)

SAMPLE_ORDER <- c(
  "K-12 MG1655 (negative control)", "O157:H7 Sakai (reference)", "O157:H7 EDL933",
  "O145:H28 RM12581", "O121:H19 FWSEC0006", "O111:H8 7-58 72A",
  "O104:H4 2011C 3493", "O103:H2 12009", "O45:H2 FWSEC0003", "O26:H11 11368"
)

main <- function() {
  
  # 1. Environment Setup (Local Root)
  current_wd <- getwd()
  cat("\n[INFO] Environment: R v4.5.2 / Bioconductor v3.22\n")
  cat("[INFO] Working Directory:", current_wd, "\n")
  
  dir_output <- file.path(current_wd, CONFIG$DIR_OUT)
  if (!dir.exists(dir_output)) {
    dir.create(dir_output, recursive = TRUE)
    cat("[INFO] Output directory created:", dir_output, "\n")
  }
  
  # 2. Load Metadata (Check current dir only)
  cat("[1/6] Loading Metadata...\n")
  path_meta <- file.path(current_wd, CONFIG$FILE_META_NAME)
  
  if (!file.exists(path_meta)) {
    stop("\n[CRITICAL ERROR] 'Metadados.tsv' not found in the current directory.\n",
         "  -> Expected at: ", path_meta, "\n",
         "  -> ACTION: Move the file to this folder or setwd() correctly.")
  }
  metadados <- readr::read_tsv(path_meta, show_col_types = FALSE)
  
  # 3. Find and Process BUSCO Files (Check current dir only)
  cat("[2/6] Processing BUSCO Files...\n")
  busco_files <- list.files(path = current_wd, pattern = CONFIG$PATTERN_BUSCO, 
                            full.names = TRUE, ignore.case = TRUE)
  
  if (length(busco_files) == 0) {
    stop("\n[CRITICAL ERROR] No BUSCO files found in the current directory.\n",
         "  -> Pattern used: ", CONFIG$PATTERN_BUSCO)
  }
  
  # Iterate and Bind
  busco_results <- purrr::map_dfr(busco_files, ~parse_busco_file(.x, metadados))
  
  # 4. Prepare Data for Analysis (Wide Format)
  cat("[3/6] Consolidating Data...\n")
  table_wide <- busco_results %>%
    dplyr::select(Sample, Status, Percentage) %>%
    tidyr::pivot_wider(names_from = Status, values_from = Percentage, values_fill = 0) %>%
    dplyr::mutate(
      Total_Complete = Complete + Duplicated,
      Group = dplyr::case_when(
        stringr::str_detect(Sample, "K-12|MG1655") ~ "Commensal (Ref)",
        TRUE ~ "Pathogenic"
      )
    )
  
  # Export Data
  openxlsx::write.xlsx(table_wide, file.path(dir_output, "Table_BUSCO_Completeness.xlsx"))
  
  # 5. Statistical Analysis
  cat("[4/6] Running Statistics...\n")
  stats_res <- run_statistics(table_wide)
  
  if (!is.null(stats_res)) {
    readr::write_csv(stats_res$table, file.path(dir_output, "Statistical_Analysis_Results.csv"))
    cat(" > Statistics computed:", stats_res$method, "\n")
  }
  
  # 6. Visualization
  cat("[5/6] Generating Plots...\n")
  
  # A. Stacked Bar Plot
  p1 <- busco_results %>%
    dplyr::mutate(Sample = factor(Sample, levels = SAMPLE_ORDER),
                  Status = factor(Status, levels = c("Complete", "Duplicated", "Fragmented", "Missing"))) %>%
    ggplot(aes(x = Sample, y = Percentage, fill = Status)) +
    geom_col(color = "black", width = 0.8, linewidth = 0.2) +
    scale_fill_manual(values = c("Complete"="#5CB85C", "Duplicated"="#5BC0DE", 
                                 "Fragmented"="#F0AD4E", "Missing"="#D9534F")) +
    geom_hline(yintercept = 95, linetype = "dotted", alpha = 0.5) +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(title = "A. BUSCO Genome Assessment Profile", x = NULL, y = "Percentage (%)")
  
  # B. Statistical Boxplot (Only if stats available)
  if (!is.null(stats_res)) {
    p2 <- table_wide %>%
      dplyr::filter(Group == "Pathogenic") %>%
      ggplot(aes(x = Group, y = Total_Complete)) +
      geom_boxplot(fill = "gray95", outlier.shape = NA, width = 0.5) +
      geom_jitter(width = 0.1, size = 3, aes(color = Total_Complete)) +
      scale_color_viridis(option = "D", guide = "none") +
      geom_hline(aes(yintercept = stats_res$ref_val, linetype = "Ref (K-12)"), color = "#D9534F") +
      scale_linetype_manual(name = "", values = "dashed") +
      annotate("text", x = 1, y = max(table_wide$Total_Complete) + 0.2, 
               label = paste0("p = ", format.pval(stats_res$table$P_Value, digits=3)), 
               size = 3.5, fontface = "italic") +
      theme_nature() +
      labs(title = "B. Statistical Comparison", 
           subtitle = paste0("Pathogenic vs K-12 (", stats_res$method, ")"),
           y = "Total Completeness (%)", x = NULL)
  } else {
    p2 <- ggplot() + theme_void() # Empty placeholder
  }
  
  # C. Heatmap
  p3 <- busco_results %>%
    dplyr::filter(Status %in% c("Complete", "Duplicated")) %>%
    dplyr::mutate(Sample = factor(Sample, levels = SAMPLE_ORDER)) %>%
    ggplot(aes(x = Status, y = Sample, fill = Percentage)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.1f", Percentage)), color = "white", size = 3, fontface = "bold") +
    scale_fill_viridis(option = "D", begin = 0.2, end = 0.8, name = "%") +
    theme_nature() +
    theme(legend.position = "right") +
    labs(title = "C. Completeness Heatmap", x = NULL, y = NULL)
  
  # Compose
  fig_final <- (p1 / (p2 | p3)) + 
    patchwork::plot_layout(heights = c(1.5, 1)) +
    patchwork::plot_annotation(
      title = "Genomic Completeness Analysis: E. coli Pathogenic vs Commensal Strains",
      theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
    )
  
  # Save
  cat("[6/6] Saving Figures...\n")
  ggsave(file.path(dir_output, "Figure_BUSCO_Statistical_Analysis.png"), 
         fig_final, width = 12, height = 10, dpi = 600, bg = "white")
  
  # Save Session Info
  writeLines(c(
    paste("Pipeline Run Date:", Sys.time()),
    paste("R Version:", R.version.string),
    "Bioconductor: v3.22 (Expected)",
    capture.output(sessionInfo())
  ), file.path(dir_output, "session_info_BUSCO.txt"))
  
  cat("\n✓ PIPELINE COMPLETED SUCCESSFULLY.\n")
  cat("✓ Results saved in:", dir_output, "\n")
}

# Execute Main
if (!interactive()) {
  main()
} else {
  cat("Script loaded. Ensure BUSCO files and Metadata are in working directory and run 'main()' to execute.\n")
}