# ==============================================================================
# Title:       Pan-genome Analysis of Shiga Toxin-Producing Escherichia coli
#              (STEC) Using Roary and Prokka Output Data
#
# Description: Constructs the pan-genome from Roary gene presence/absence
#              matrices to isolate core genes conserved across 100% of
#              pathogenic STEC strains and absent in the negative control
#              (K-12). Validates the ubiquity of the stx2a target through
#              Heaps' Law modelling, diversity indices, hierarchical clustering
#              heatmaps, and UpSet intersection plots.
#
# Module:      3 – Comparative Pan-genomics and Core Genome Delimitation
#
# Input files:
#   Galaxy41_Roary_Summary_statistics.csv    – Roary summary table
#   Galaxy42_Roary_Core_Gene_Alignment.fasta – Core-gene alignment
#   Galaxy43_Roary_Gene_Presence_Absence.csv – Binary presence/absence matrix
#   Galaxy44_Roary_Accessory_Gene_Table.csv  – Accessory gene table
#   Galaxy45_Roary_Core_Accessory_Table.csv  – Core + accessory merged table
#   Metadados.csv / Metadados.tsv            – Strain name mapping (optional)
#
# Output:
#   Graficos/      – High-resolution figures (PNG, 600 DPI)
#   Tabelas/       – Exported result tables (CSV)
#   Documentacao/  – Figure legends, session info, and README
#
# Usage:
#   Rscript 05_Prokka_Roary_Analysis_revised.R
#   Or source interactively in RStudio (File > Source with Echo)
#
# Dependencies:  See Section 1 (auto-installed if absent).
#   R         >= 4.3.0
#   Bioconductor >= 3.18
#
# Author(s):   [AUTHOR NAME]
# Affiliation: [Department], [Institution], [Country]
# Contact:     [email@institution.edu]
# ORCID:       [0000-0000-0000-0000]
#
# Created:      2026-02-13
# Last modified: 2026-02-13
# Version:      1.0.0
#
# Citation:
#   [Author] et al. (2026) [Article title]. [Journal]. DOI: [10.xxxx/xxxxx]
#
#   Software dependencies:
#   Seemann T (2014) Prokka: rapid prokaryotic genome annotation.
#     Bioinformatics 30(14):2068-2069. doi:10.1093/bioinformatics/btu153
#   Page AJ et al. (2015) Roary: the fast pan-genome pipeline.
#     Bioinformatics 31(22):3691-3693. doi:10.1093/bioinformatics/btv421
#
# License:     MIT License
#   Copyright (c) 2026 [Author Name]
#   Permission is hereby granted, free of charge, to any person obtaining
#   a copy of this software to use, copy, modify, merge, publish, and
#   distribute copies of the Software without restriction, subject to the
#   conditions of the MIT License.
#
# Repository:  https://github.com/[user]/[repo]
# DOI (repo):  https://doi.org/10.5281/zenodo.[id]
# ==============================================================================


# ==============================================================================
# SECTION 1 – INITIAL SETUP AND PACKAGE MANAGEMENT
# ==============================================================================

cat("\n")
cat("==============================================================================\n")
cat("     PAN-GENOME ANALYSIS – ROARY / PROKKA DATA PROCESSING                   \n")
cat("==============================================================================\n\n")

execution_start <- Sys.time()

cat(sprintf("Start time : %s\n", format(execution_start, "%Y-%m-%d %H:%M:%S")))
cat(sprintf("R version  : %s\n\n", R.version.string))

# Suppress scientific notation for readability in printed output
options(scipen = 999)

# Global seed for all stochastic operations (permutations, sampling)
# Value 42 is arbitrary; fixing it guarantees reproducible output across runs
set.seed(42)

# --- 1.1  BiocManager --------------------------------------------------------
cat("1. Checking BiocManager...\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cran.r-project.org")
}

# --- 1.2  Package lists -------------------------------------------------------
bioc_packages <- c(
  "ComplexHeatmap",  # Publication-quality heatmaps (Gu et al., 2016)
  "Biostrings",      # Biological sequence containers and algorithms
  "GenomicRanges"    # Genomic interval operations
)

cran_packages <- c(
  "tidyverse",    # Data manipulation and ggplot2 visualisation ecosystem
  "data.table",   # High-performance tabular data I/O and manipulation
  "ggpubr",       # ggplot2 publication-ready helpers
  "patchwork",    # Multi-panel figure composition
  "cowplot",      # Publication-grade plot layouts
  "viridis",      # Perceptually uniform, colour-blind-safe palettes
  "RColorBrewer", # Qualitative and sequential colour palettes
  "ggsci",        # Scientific journal colour palettes
  "ggrepel",      # Non-overlapping text labels in ggplot2
  "scales",       # Axis scale formatting utilities
  "vegan",        # Community ecology: diversity indices (Oksanen et al.)
  "Cairo",        # High-quality raster graphics device
  "UpSetR",       # UpSet intersection plots (Conway et al., 2017)
  "pheatmap",     # Clustered heatmaps with annotations
  "dendextend"    # Dendrogram manipulation and styling
)

# --- 1.3  Installation helper ------------------------------------------------

#' Install a package if it is not already available
#'
#' @description
#' Checks whether \code{package} can be loaded via
#' \code{\link[base]{requireNamespace}} and installs it silently if not.
#' Supports both CRAN and Bioconductor repositories.
#'
#' @param package Character scalar. Name of the R package to check/install.
#' @param source  Character scalar. Either \code{"CRAN"} (default) or
#'   \code{"Bioconductor"}.
#'
#' @return Invisible \code{NULL}. Called for its side-effect.
#'
#' @examples
#' \dontrun{
#' install_if_missing("ggplot2", source = "CRAN")
#' install_if_missing("ComplexHeatmap", source = "Bioconductor")
#' }
install_if_missing <- function(package, source = "CRAN") {
  if (!requireNamespace(package, quietly = TRUE)) {
    message(sprintf("  Installing %s from %s...", package, source))
    if (source == "Bioconductor") {
      BiocManager::install(package, update = FALSE, ask = FALSE)
    } else {
      install.packages(
        package,
        repos        = "https://cran.r-project.org",
        dependencies = TRUE,
        quiet        = TRUE
      )
    }
  }
  invisible(NULL)
}

cat("2. Verifying and installing required packages...\n")
invisible(lapply(bioc_packages,  install_if_missing, source = "Bioconductor"))
invisible(lapply(cran_packages,  install_if_missing, source = "CRAN"))
cat("   All required packages verified.\n\n")


# ==============================================================================
# SECTION 2 – LIBRARY LOADING
# ==============================================================================

cat("3. Loading libraries...\n")

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ComplexHeatmap)
  library(Biostrings)
  library(ggpubr)
  library(patchwork)
  library(cowplot)
  library(viridis)
  library(RColorBrewer)
  library(ggsci)
  library(ggrepel)
  library(scales)
  library(vegan)
  library(Cairo)
  library(UpSetR)
  library(pheatmap)
  library(dendextend)
})

cat("   Libraries loaded successfully.\n\n")


# ==============================================================================
# SECTION 3 – GLOBAL PARAMETERS AND DIRECTORY STRUCTURE
# ==============================================================================

cat("4. Configuring directories and global parameters...\n")

# Resolve the script's own directory, falling back to the current working
# directory when sourced non-interactively (e.g., via Rscript on the command
# line). This ensures relative paths remain valid regardless of where the
# script is invoked.
script_dir <- tryCatch(
  {
    if (requireNamespace("rstudioapi", quietly = TRUE) &&
        rstudioapi::isAvailable()) {
      dirname(rstudioapi::getSourceEditorContext()$path)
    } else {
      # Fallback for Rscript execution
      dirname(sys.frames()[[1]]$ofile)
    }
  },
  error = function(e) getwd()
)

if (is.null(script_dir) || !nzchar(script_dir)) {
  script_dir <- getwd()
}

setwd(script_dir)
cat(sprintf("  Working directory : %s\n", script_dir))

input_dir  <- script_dir
output_dir <- script_dir

# Create output subdirectories (silent if already present)
output_dirs <- c("Graficos", "Tabelas", "Documentacao")
invisible(lapply(
  file.path(output_dir, output_dirs),
  dir.create,
  showWarnings = FALSE,
  recursive    = TRUE
))

# --- 3.1  Graphical constants -------------------------------------------------
# Resolution and dimensions follow Nature / Nature Methods figure guidelines:
# minimum 300 DPI for colour images; 600 DPI recommended for line art.
DPI             <- 600L   # dots per inch for all saved figures
WIDTH_INCHES    <- 10     # default figure width  (inches)
HEIGHT_INCHES   <- 8      # default figure height (inches)
FONT_SIZE_BASE  <- 12     # base font size (points)
FONT_SIZE_TITLE <- 14     # title font size (points)
FONT_SIZE_AXIS  <- 10     # axis-label font size (points)

# --- 3.2  Colour palettes -----------------------------------------------------
# Viridis option D is perceptually uniform and safe for the most common forms
# of colour-vision deficiency (deuteranopia, protanopia).
PALETTE_CATEGORIES <- viridis(4L, option = "D")
names(PALETTE_CATEGORIES) <- c("Core", "Soft Core", "Shell", "Cloud")

# Extended qualitative palette for individual strain labels
PALETTE_STRAINS <- c(brewer.pal(9L, "Set1"), brewer.pal(8L, "Set2"))

cat("   Parameters configured.\n\n")


# ==============================================================================
# SECTION 4 – HELPER FUNCTIONS
# ==============================================================================

cat("5. Defining helper functions...\n")

#' Save a ggplot object as a high-resolution PNG via the Cairo device
#'
#' @description
#' Writes \code{plot} to \file{Graficos/<filename>} in the global
#' \code{output_dir}, using the Cairo PNG device for anti-aliased rendering
#' suitable for print publication.
#'
#' @param plot     A \code{ggplot} or \code{patchwork} object.
#' @param filename Character scalar. File name including the \file{.png}
#'   extension (e.g., \code{"Fig1_Distribution.png"}).
#' @param width    Numeric scalar. Figure width in inches (default:
#'   \code{WIDTH_INCHES}).
#' @param height   Numeric scalar. Figure height in inches (default:
#'   \code{HEIGHT_INCHES}).
#' @param dpi      Integer scalar. Output resolution in dots per inch
#'   (default: \code{DPI}).
#'
#' @return Invisible character scalar: the full file path of the saved figure.
#'
#' @examples
#' \dontrun{
#' p <- ggplot(mtcars, aes(wt, mpg)) + geom_point()
#' save_figure(p, "test_figure.png", width = 8, height = 6)
#' }
save_figure <- function(plot,
                        filename,
                        width  = WIDTH_INCHES,
                        height = HEIGHT_INCHES,
                        dpi    = DPI) {
  filepath <- file.path(output_dir, "Graficos", filename)

  CairoPNG(
    filename = filepath,
    width    = width,
    height   = height,
    units    = "in",
    dpi      = dpi,
    bg       = "white"
  )
  print(plot)
  dev.off()

  message(sprintf("  Figure saved: %s", filename))
  invisible(filepath)
}


#' Export a data frame to CSV using data.table::fwrite
#'
#' @description
#' Writes \code{df} to \file{Tabelas/<filename>} in \code{output_dir}.
#' Uses \code{\link[data.table]{fwrite}} for fast, locale-independent output.
#'
#' @param df       A \code{data.frame} or \code{data.table}.
#' @param filename Character scalar. File name including the \file{.csv}
#'   extension (e.g., \code{"Table_Core_Genes.csv"}).
#'
#' @return Invisible character scalar: the full file path of the saved table.
save_table <- function(df, filename) {
  filepath <- file.path(output_dir, "Tabelas", filename)
  fwrite(df, filepath, sep = ",", row.names = FALSE)
  message(sprintf("  Table saved: %s", filename))
  invisible(filepath)
}


#' Apply a publication-ready ggplot2 theme
#'
#' @description
#' Extends \code{\link[ggplot2]{theme_bw}} with formatting conventions
#' appropriate for print journals: bold axis titles, minimal grid lines,
#' and a clean white background. Font sizes are parameterised through the
#' global constants \code{FONT_SIZE_AXIS}, \code{FONT_SIZE_BASE}, and
#' \code{FONT_SIZE_TITLE}.
#'
#' @param base_size Numeric scalar. Base font size in points
#'   (default: \code{FONT_SIZE_BASE}).
#'
#' @return A \code{\link[ggplot2]{theme}} object.
#'
#' @examples
#' \dontrun{
#' ggplot(mtcars, aes(wt, mpg)) +
#'   geom_point() +
#'   theme_publication()
#' }
theme_publication <- function(base_size = FONT_SIZE_BASE) {
  theme_bw(base_size = base_size) +
    theme(
      # Grid – retain only major lines for readability
      panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"),
      panel.grid.minor = element_blank(),

      # Axis text and titles
      axis.text    = element_text(size = FONT_SIZE_AXIS, colour = "black"),
      axis.title   = element_text(size = FONT_SIZE_BASE, face = "bold"),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),

      # Plot title and subtitle
      plot.title    = element_text(
        size   = FONT_SIZE_TITLE, face = "bold",
        hjust  = 0.5, margin = margin(b = 10)
      ),
      plot.subtitle = element_text(
        size   = FONT_SIZE_BASE, hjust = 0.5,
        margin = margin(b = 10)
      ),

      # Legend
      legend.text       = element_text(size = FONT_SIZE_AXIS),
      legend.title      = element_text(size = FONT_SIZE_BASE, face = "bold"),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.key        = element_rect(fill = "white", colour = NA),

      # Outer margins (t, r, b, l)
      plot.margin = margin(10, 10, 10, 10)
    )
}


#' Load and return strain metadata from CSV or TSV
#'
#' @description
#' Searches \code{input_dir} for \file{Metadados.csv} or
#' \file{Metadados.tsv} (in that order) and loads the first file found.
#' Returns \code{NULL} silently when no file is present, allowing
#' downstream code to fall back to the original Roary column names.
#'
#' @return A \code{data.table} with at least the columns \code{Roary}
#'   (Roary column identifier) and \code{Nome_Real} (display name), or
#'   \code{NULL} if no metadata file is found.
load_metadata <- function() {
  csv_path <- file.path(input_dir, "Metadados.csv")
  tsv_path <- file.path(input_dir, "Metadados.tsv")

  if (file.exists(csv_path)) {
    meta <- fread(csv_path, sep = ",", header = TRUE)
    message(sprintf("  Metadata loaded: %d mappings (Metadados.csv)", nrow(meta)))
    return(meta)
  }

  if (file.exists(tsv_path)) {
    meta <- fread(tsv_path, sep = "\t", header = TRUE)
    message(sprintf("  Metadata loaded: %d mappings (Metadados.tsv)", nrow(meta)))
    return(meta)
  }

  message("  No metadata file found; using original strain identifiers.")
  NULL
}


#' Map Roary column identifiers to display strain names
#'
#' @description
#' Performs an exact lookup of each element of \code{original_names} in the
#' \code{Roary} column of \code{metadata}, returning the corresponding
#' \code{Nome_Real} value. Unmatched identifiers are returned unchanged with
#' a warning.
#'
#' @param original_names Character vector of Roary column names to remap.
#' @param metadata A \code{data.table} as returned by \code{load_metadata},
#'   with columns \code{Roary} and \code{Nome_Real}. If \code{NULL}, the
#'   function returns \code{original_names} unmodified.
#'
#' @return Character vector of the same length as \code{original_names}.
map_strain_names <- function(original_names, metadata) {
  if (is.null(metadata)) return(original_names)

  required_cols <- c("Roary", "Nome_Real")
  if (!all(required_cols %in% colnames(metadata))) {
    warning(
      "Columns 'Roary' and 'Nome_Real' not found in metadata; ",
      "keeping original names."
    )
    return(original_names)
  }

  mapped <- vapply(original_names, function(name) {
    idx <- which(metadata[["Roary"]] == name)
    if (length(idx) > 0L) {
      new_name <- metadata[["Nome_Real"]][idx[[1L]]]
      message(sprintf("  Mapping: %s -> %s", name, new_name))
      new_name
    } else {
      warning(sprintf("No mapping found for strain: %s", name))
      name
    }
  }, character(1L))

  unname(mapped)
}


cat("   Helper functions defined.\n\n")


# ==============================================================================
# SECTION 5 – DATA IMPORT
# ==============================================================================

cat("6. Importing datasets...\n")

metadata <- load_metadata()

# -- Galaxy41: Roary summary statistics ---------------------------------------
cat("  Loading Galaxy41_Roary_Summary_statistics.csv...\n")
summary_stats <- fread(
  file.path(input_dir, "Galaxy41_Roary_Summary_statistics.csv"),
  header    = FALSE,
  col.names = c("Category", "Definition", "Count")
)

# -- Galaxy42: Core-gene alignment (FASTA) ------------------------------------
cat("  Loading Galaxy42_Roary_Core_Gene_Alignment.fasta...\n")
core_alignment_file <- file.path(input_dir, "Galaxy42_Roary_Core_Gene_Alignment.fasta")

if (file.exists(core_alignment_file)) {
  core_alignment <- readDNAStringSet(core_alignment_file)
  cat(sprintf("  Core alignment loaded: %d sequences\n", length(core_alignment)))
} else {
  core_alignment <- NULL
  warning("Galaxy42_Roary_Core_Gene_Alignment.fasta not found; continuing without it.")
  cat("  Core alignment unavailable; continuing analysis.\n")
}

# -- Galaxy43: Gene presence/absence matrix -----------------------------------
cat("  Loading Galaxy43_Roary_Gene_Presence_Absence.csv...\n")
gene_pa <- fread(
  file.path(input_dir, "Galaxy43_Roary_Gene_Presence_Absence.csv"),
  header          = TRUE,
  stringsAsFactors = FALSE
)

# -- Galaxy44: Accessory gene table -------------------------------------------
cat("  Loading Galaxy44_Roary_Accessory_Gene_Table.csv...\n")
accessory_genes <- fread(
  file.path(input_dir, "Galaxy44_Roary_Accessory_Gene_Table.csv"),
  header          = TRUE,
  stringsAsFactors = FALSE,
  fill            = TRUE   # rows with variable field counts are allowed
)

# -- Galaxy45: Core + accessory merged table ----------------------------------
cat("  Loading Galaxy45_Roary_Core_Accessory_Table.csv...\n")
core_accessory <- fread(
  file.path(input_dir, "Galaxy45_Roary_Core_Accessory_Table.csv"),
  header          = TRUE,
  stringsAsFactors = FALSE,
  fill            = TRUE
)

cat(sprintf("\n  Datasets loaded: 5\n"))
cat(sprintf("  Genes in presence/absence matrix: %d\n\n", nrow(gene_pa)))


# ==============================================================================
# SECTION 6 – PRE-PROCESSING AND VALIDATION
# ==============================================================================

cat("7. Pre-processing and data validation...\n")

# -- 6.1  Dimension check ------------------------------------------------------
# Expected totals are hard-coded from the study design (10 STEC strains,
# 11,453 orthologous gene clusters from Roary v3.13.0 run in Galaxy).
EXPECTED_GENES  <- 11453L
EXPECTED_STRAINS <- 10L

total_genes_obs <- nrow(gene_pa)
if (abs(total_genes_obs - EXPECTED_GENES) < 50L) {
  cat(sprintf("  [OK] Gene count: %d (expected ~%d)\n",
              total_genes_obs, EXPECTED_GENES))
} else {
  warning(sprintf(
    "Gene count differs from expectation: %d vs %d",
    total_genes_obs, EXPECTED_GENES
  ))
}

# Roary encodes strain columns with the pattern "Prokka_on_dataset_*__gff"
strain_col_idx <- grep("^Prokka_on_dataset_.*__gff$", colnames(gene_pa))
n_strains      <- length(strain_col_idx)
cat(sprintf("  Strains identified: %d\n", n_strains))

# -- 6.2  Strain name mapping --------------------------------------------------
original_col_names <- colnames(gene_pa)

if (!is.null(metadata)) {
  cat("  Applying strain name mapping from metadata...\n")
  strain_names <- map_strain_names(
    original_col_names[strain_col_idx],
    metadata
  )
  colnames(gene_pa)[strain_col_idx] <- strain_names
} else {
  # Simplify default Roary column identifiers for readability
  strain_names <- str_replace_all(
    original_col_names[strain_col_idx],
    c("^Prokka_on_dataset_" = "Strain_", "__gff$" = "")
  )
  colnames(gene_pa)[strain_col_idx] <- strain_names
}

cat("\n  Strain names after processing:\n")
walk(seq_along(strain_names), ~ cat(sprintf("    %d. %s\n", .x, strain_names[.x])))

# Guard: duplicated strain names will break matrix indexing downstream
if (anyDuplicated(strain_names) > 0L) {
  stop(
    "Duplicate strain names detected: ",
    paste(strain_names[duplicated(strain_names)], collapse = ", "),
    "\nResolve duplicates in Metadados.csv before proceeding."
  )
}
cat("\n  [OK] All strain names are unique.\n")

# -- 6.3  Binary presence/absence matrix ---------------------------------------
cat("\n  Building binary presence/absence matrix...\n")

# Extract strain columns; empty strings and NA encode absence (0)
pa_matrix_raw <- gene_pa[, ..strain_names]
pa_binary <- pa_matrix_raw[, lapply(.SD, function(x) {
  as.integer(nzchar(x) & !is.na(x))
})]

gene_names <- gene_pa[["Gene"]]

cat(sprintf("  Binary matrix: %d genes × %d strains\n",
            nrow(pa_binary), ncol(pa_binary)))

# -- 6.4  Gene category classification ----------------------------------------
# Thresholds follow the Roary convention:
#   Core     : present in >= 99% of strains  (all 10 strains = 100%)
#   Soft core: present in 95-99% of strains
#   Shell    : present in 15-95% of strains
#   Cloud    : present in < 15% of strains   (rare / strain-specific)
#
# NOTE on Soft Core with n = 10 strains:
#   95% of 10 = 9.5  and  99% of 10 = 9.9 — both non-integers.
#   Since presence counts are whole numbers, no gene can have a percentage
#   strictly between 95% and 99% with 10 strains (9 strains = 90%, 10 = 100%).
#   Soft Core count = 0 is therefore mathematically expected, not a data error.
cat("  Classifying genes into pan-genome categories...\n")

gene_pa[, n_strains_present := rowSums(pa_binary)]
gene_pa[, pct_strains        := (n_strains_present / n_strains) * 100]

gene_pa[, category := cut(
  pct_strains,
  breaks = c(-Inf, 15, 95, 99, Inf),
  labels = c("Cloud", "Shell", "Soft Core", "Core"),
  right  = FALSE
)]

category_tbl <- gene_pa[, .(count = .N,
                             percent = round(.N / nrow(gene_pa) * 100, 1)),
                         by = category][order(-count)]

cat("  Gene category distribution:\n")
print(as.data.frame(category_tbl), row.names = FALSE)

cat("\n  Pre-processing complete.\n\n")


# ==============================================================================
# SECTION 7 – STATISTICAL ANALYSES
# ==============================================================================

cat("8. Running statistical analyses...\n")

# -- 7.1  Descriptive statistics -----------------------------------------------
cat("  7.1  Descriptive statistics...\n")

stats_pangenome <- data.frame(
  Metric = c(
    "Total genes (pan-genome)",
    "Core genome (99-100%)",
    "Soft core genome (95-99%)",
    "Shell genes (15-95%)",
    "Cloud genes (0-15%)",
    "Mean singleton genes per strain",
    "Mean genome size (genes)"
  ),
  Value = c(
    nrow(gene_pa),
    sum(gene_pa[["category"]] == "Core",      na.rm = TRUE),
    sum(gene_pa[["category"]] == "Soft Core", na.rm = TRUE),
    sum(gene_pa[["category"]] == "Shell",     na.rm = TRUE),
    sum(gene_pa[["category"]] == "Cloud",     na.rm = TRUE),
    # Mean number of genes present exclusively in one strain (singletons)
    mean(vapply(strain_names, function(s) {
      sum(pa_binary[[s]] == 1L & rowSums(pa_binary) == 1L)
    }, integer(1L))),
    # Mean total genes annotated per genome (genome size proxy)
    mean(colSums(pa_binary))
  )
)

# -- 7.2  Strain-specific (singleton) genes ------------------------------------
cat("  7.2  Identifying strain-specific (singleton) genes...\n")

# A singleton gene is present in exactly one strain
singleton_per_strain <- data.frame(
  Strain       = strain_names,
  Unique_genes = vapply(strain_names, function(s) {
    sum(pa_binary[[s]] == 1L & rowSums(pa_binary) == 1L)
  }, integer(1L))
) |>
  dplyr::arrange(dplyr::desc(Unique_genes))

cat(sprintf("  Total singleton genes: %d\n",
            sum(singleton_per_strain[["Unique_genes"]])))

# -- 7.3  Pairwise gene-sharing matrix -----------------------------------------
cat("  7.3  Computing pairwise gene-sharing matrix...\n")

# Element [i,j] = number of genes present in both strain i and strain j.
# tcrossprod(t(M)) is equivalent to M %*% t(M) and yields this directly in one
# vectorised operation — no loop needed. dimnames are assigned separately for
# clarity.
pa_matrix_num <- as.matrix(pa_binary)

sharing_matrix <- tcrossprod(t(pa_matrix_num))
dimnames(sharing_matrix) <- list(strain_names, strain_names)

# -- 7.4  Alpha-diversity indices (Shannon and Simpson) ------------------------
cat("  7.4  Computing alpha-diversity indices...\n")

# vegan::diversity expects samples as rows; transpose accordingly
pa_transposed  <- t(pa_matrix_num)
shannon_index  <- vegan::diversity(pa_transposed, index = "shannon")
simpson_index  <- vegan::diversity(pa_transposed, index = "simpson")

diversity_df <- data.frame(
  Strain  = strain_names,
  Shannon = shannon_index,
  Simpson = simpson_index
) |> dplyr::arrange(dplyr::desc(Shannon))

cat(sprintf("  Mean Shannon index : %.3f\n", mean(shannon_index)))
cat(sprintf("  Mean Simpson index : %.3f\n", mean(simpson_index)))

# -- 7.5  Heaps' Law pan-genome model ------------------------------------------
# Heaps' Law: N(g) = κ × g^γ
#   γ < 1 → open pan-genome (new genes accumulate indefinitely)
#   γ ≥ 1 → closed pan-genome (gene pool is finite)
# Parameters estimated by OLS regression on log-transformed data:
#   log(N) = log(κ) + γ × log(g)
cat("  7.5  Fitting Heaps' Law pan-genome model...\n")

set.seed(42)  # Fix randomisation for strain addition order
strain_order   <- sample(seq_len(n_strains))

pangenome_accum <- integer(n_strains)
core_accum      <- integer(n_strains)

for (i in seq_len(n_strains)) {
  strains_so_far   <- strain_order[seq_len(i)]
  sub_matrix       <- pa_matrix_num[, strains_so_far, drop = FALSE]

  # Pan-genome: union of all genes seen in ≥ 1 strain
  pangenome_accum[i] <- sum(rowSums(sub_matrix) > 0L)

  # Core genome: intersection – genes in ALL strains considered so far
  core_accum[i]      <- sum(rowSums(sub_matrix) == i)
}

heaps_model <- lm(log(pangenome_accum) ~ log(seq_len(n_strains)))
kappa_est   <- exp(coef(heaps_model)[[1L]])
gamma_est   <- coef(heaps_model)[[2L]]

cat(sprintf("  Heaps' Law parameters:\n"))
cat(sprintf("    κ (kappa) = %.2f\n",  kappa_est))
cat(sprintf("    γ (gamma) = %.4f\n",  gamma_est))
cat(sprintf("    Pan-genome type: %s (γ %s 1)\n",
            ifelse(gamma_est < 1, "OPEN", "CLOSED"),
            ifelse(gamma_est < 1, "<",    "≥")))

accumulation_df <- data.frame(
  n_genomes       = seq_len(n_strains),
  pangenome       = pangenome_accum,
  core            = core_accum,
  pangenome_model = kappa_est * seq_len(n_strains) ^ gamma_est
)

cat("\n  Statistical analyses complete.\n\n")


# ==============================================================================
# SECTION 8 – INDIVIDUAL FIGURES
# ==============================================================================

cat("9. Generating individual figures...\n\n")

# Prepare labelled distribution data frame
dist_df <- category_tbl |>
  dplyr::mutate(
    label = sprintf("%s\n%d genes (%.1f%%)", category, count, percent)
  )

# -- Figure 1: Pan-genome category distribution --------------------------------
cat("  Fig 1 – Pan-genome category distribution...\n")

fig1_distribution <- ggplot(dist_df,
                             aes(x = category, y = count, fill = category)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.5) +
  geom_text(
    aes(label = sprintf("%d\n(%.1f%%)", count, percent)),
    vjust = -0.5, size = 4, fontface = "bold"
  ) +
  scale_fill_manual(values = PALETTE_CATEGORIES) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.15)),
    labels = scales::comma
  ) +
  labs(
    title    = "Pan-genome Structure Distribution",
    subtitle = sprintf("Total: %s genes across %d strains",
                       scales::comma(nrow(gene_pa)), n_strains),
    x        = "Gene Category",
    y        = "Number of Genes",
    fill     = "Category"
  ) +
  theme_publication() +
  theme(legend.position = "none")

save_figure(fig1_distribution, "Fig01_Pangenome_Category_Distribution.png",
            width = 10, height = 8)

# -- Figure 2: Pan-genome accumulation curve with Heaps' Law fit ---------------
cat("  Fig 2 – Pan-genome accumulation curve...\n")

fig2_accumulation <- ggplot(accumulation_df, aes(x = n_genomes)) +
  geom_point(aes(y = pangenome), colour = "#E74C3C", size = 3) +
  geom_line(aes(y = pangenome), colour = "#E74C3C", linewidth = 1) +
  geom_line(
    aes(y = pangenome_model), colour = "#3498DB",
    linetype = "dashed", linewidth = 1
  ) +
  annotate(
    "text",
    x     = n_strains * 0.05,
    # Position annotation within the visible y-range:
    # the y-axis does NOT start at 0 (expand = lower 0%), so max * 0.3 falls
    # below the axis minimum (≈ 5 700). Using min + 15% of range keeps the
    # text inside the plot area regardless of the data spread.
    y     = min(accumulation_df[["pangenome"]]) +
              diff(range(accumulation_df[["pangenome"]])) * 0.15,
    label = sprintf("Heaps' Law:\nN = %.1f \u00d7 n^%.3f\n\u03b3 = %.3f (%s)",
                    kappa_est, gamma_est, gamma_est,
                    ifelse(gamma_est < 1, "Open", "Closed")),
    hjust = 0, size = 4, colour = "#3498DB"
  ) +
  scale_x_continuous(breaks = seq_len(n_strains)) +
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0.05, 0.1))  # 5% lower margin keeps annotation visible
  ) +
  labs(
    title    = "Pan-genome Accumulation Curve",
    subtitle = "Cumulative gene count with sequential genome addition",
    x        = "Number of Genomes",
    y        = "Cumulative Number of Genes"
  ) +
  theme_publication()

save_figure(fig2_accumulation, "Fig02_Pangenome_Accumulation_Curve.png",
            width = 10, height = 8)

# -- Figure 3: Core genome reduction curve ------------------------------------
cat("  Fig 3 – Core genome reduction curve...\n")

fig3_core_reduction <- ggplot(accumulation_df,
                               aes(x = n_genomes, y = core)) +
  geom_point(colour = "#27AE60", size = 3) +
  geom_line(colour  = "#27AE60", linewidth = 1) +
  geom_hline(
    yintercept = dplyr::last(accumulation_df[["core"]]),
    linetype   = "dashed", colour = "grey40"
  ) +
  annotate(
    "text",
    x     = n_strains * 0.6,
    y     = max(accumulation_df[["core"]]) * 0.8,
    label = sprintf("Core genome plateau:\n%d genes",
                    dplyr::last(accumulation_df[["core"]])),
    hjust = 0, size = 4, colour = "#27AE60"
  ) +
  scale_x_continuous(breaks = seq_len(n_strains)) +
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.1))
  ) +
  labs(
    title    = "Core Genome Reduction Curve",
    subtitle = "Core genome size with sequential genome addition",
    x        = "Number of Genomes",
    y        = "Number of Core Genes"
  ) +
  theme_publication()

save_figure(fig3_core_reduction, "Fig03_Core_Genome_Reduction_Curve.png",
            width = 10, height = 8)

# -- Figure 4: Gene presence/absence heatmap (representative sample) -----------
cat("  Fig 4 – Presence/absence heatmap...\n")

heatmap_matrix <- as.matrix(pa_binary)
rownames(heatmap_matrix) <- gene_names

gene_annotation_df <- data.frame(Category = gene_pa[["category"]])
rownames(gene_annotation_df) <- gene_names

# Downsample to a representative subset (max 500 genes per category) because
# the full matrix (11,453 × 10) would produce an unreadable figure. The random
# sample is seeded for reproducibility.
set.seed(42)
# slice_sample(prop = 1) shuffles all rows within each category group;
# slice_head(n = 500L) then caps at 500 per group, returning fewer rows
# for categories that have < 500 genes without raising an error.
sampled_genes <- gene_pa |>
  dplyr::group_by(category) |>
  dplyr::slice_sample(prop = 1) |>
  dplyr::slice_head(n = 500L) |>
  dplyr::pull(Gene)

valid_genes <- sampled_genes[sampled_genes %in% rownames(heatmap_matrix)]
n_sampled   <- length(valid_genes)

if (n_sampled < 10L) {
  warning("Fewer than 10 valid genes after sampling; heatmap skipped.")
} else {
  cat(sprintf("  Heatmap: %d genes × %d strains\n", n_sampled, n_strains))

  heatmap_sample      <- heatmap_matrix[valid_genes, ]
  annotation_sample   <- gene_annotation_df[valid_genes, , drop = FALSE]
  annotation_colours  <- list(
    Category = setNames(PALETTE_CATEGORIES, c("Core", "Soft Core", "Shell", "Cloud"))
  )

  png(
    filename = file.path(output_dir, "Graficos", "Fig04_Gene_Presence_Absence_Heatmap.png"),
    width    = WIDTH_INCHES,
    height   = HEIGHT_INCHES,
    units    = "in",
    res      = DPI,
    type     = "cairo"
  )

  pheatmap(
    heatmap_sample,
    # NOTE: pheatmap uses 'color', not 'colour'. Using 'colour' silently falls
    # back to the default blue-red gradient — that is why the heatmap rendered
    # with colour instead of the intended binary white / dark-navy palette.
    color                    = c("white", "#2C3E50"),
    border_color             = NA,
    cluster_rows             = TRUE,
    cluster_cols             = TRUE,
    clustering_distance_rows = "binary",
    clustering_distance_cols = "binary",
    clustering_method        = "average",
    annotation_row           = annotation_sample,
    annotation_colors        = annotation_colours,
    show_rownames            = FALSE,
    show_colnames            = TRUE,
    fontsize                 = 10,
    fontsize_col             = 10,
    main                     = "Gene Presence/Absence Matrix (Representative Sample)",
    legend_breaks            = c(0, 1),
    legend_labels            = c("Absent", "Present")
  )

  dev.off()
  cat("  Figure saved: Fig04_Gene_Presence_Absence_Heatmap.png\n")
}

# -- Figure 5: Strain-specific (singleton) genes -------------------------------
cat("  Fig 5 – Strain-specific genes...\n")

fig5_singletons <- ggplot(
  singleton_per_strain,
  aes(x = reorder(Strain, Unique_genes), y = Unique_genes)
) +
  geom_bar(
    stat      = "identity",
    fill      = "#9B59B6",
    colour    = "black",
    linewidth = 0.5
  ) +
  geom_text(aes(label = Unique_genes), hjust = -0.2, size = 3.5) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title    = "Strain-Specific Genes (Singletons)",
    subtitle = "Genes present exclusively in one strain",
    x        = "Strain",
    y        = "Number of Unique Genes"
  ) +
  theme_publication()

save_figure(fig5_singletons, "Fig05_Singleton_Genes_Per_Strain.png",
            width = 10, height = 8)

# -- Figure 6: UpSet intersection plot -----------------------------------------
cat("  Fig 6 – UpSet intersection plot...\n")

# Build a named list: each element is the set of gene names present in that
# strain. UpSetR::fromList converts this to the required binary matrix.
gene_sets_by_strain <- setNames(
  lapply(strain_names, function(s) gene_names[pa_binary[[s]] == 1L]),
  strain_names
)

png(
  filename = file.path(output_dir, "Graficos", "Fig06_UpSet_Gene_Intersections.png"),
  width    = 12,
  height   = 8,
  units    = "in",
  res      = DPI,
  type     = "cairo"
)

# suppressWarnings() silences three deprecation warnings emitted by UpSetR's
# own internals (aes_string(), size aesthetic, element_line(size=)):
# these originate inside the package source — not in our code — and have been
# reported upstream. They do not affect plot correctness or reproducibility.
suppressWarnings(
  UpSetR::upset(
    UpSetR::fromList(gene_sets_by_strain),
    nsets        = n_strains,
    nintersects  = 30,       # display top 30 intersections by frequency
    order.by     = "freq",
    decreasing   = TRUE,
    mb.ratio     = c(0.6, 0.4),
    text.scale   = c(1.5, 1.3, 1.2, 1.2, 1.5, 1.3),
    main.bar.color = "#E67E22",
    matrix.color   = "#E67E22",
    sets.bar.color = "#3498DB"
  )
)

dev.off()
cat("  Figure saved: Fig06_UpSet_Gene_Intersections.png\n")

# -- Figure 7: Alpha-diversity indices -----------------------------------------
cat("  Fig 7 – Alpha-diversity indices...\n")

fig7a_shannon <- ggplot(
  diversity_df,
  aes(x = reorder(Strain, Shannon), y = Shannon)
) +
  geom_bar(
    stat = "identity", fill = "#1ABC9C",
    colour = "black", linewidth = 0.5
  ) +
  geom_text(aes(label = sprintf("%.3f", Shannon)), hjust = -0.2, size = 3) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Shannon Diversity Index", x = "Strain", y = "Shannon Index") +
  theme_publication()

fig7b_simpson <- ggplot(
  diversity_df,
  aes(x = reorder(Strain, Simpson), y = Simpson)
) +
  geom_bar(
    stat = "identity", fill = "#E91E63",
    colour = "black", linewidth = 0.5
  ) +
  geom_text(aes(label = sprintf("%.3f", Simpson)), hjust = -0.2, size = 3) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Simpson Diversity Index",
    # Simpson saturates at 1 for binary data with large N:
    # 1 - sum(p_i^2) = 1 - N*(1/N)^2 = 1 - 1/N → 1 as N → ∞.
    # All strains show 1.000; Shannon (panel A) is the informative metric.
    subtitle = "Note: saturates at 1 with large binary gene sets (see legend)",
    x = "Strain",
    y = "Simpson Index"
  ) +
  theme_publication()

fig7_diversity <- (fig7a_shannon / fig7b_simpson) +
  patchwork::plot_annotation(
    title    = "Genetic Diversity Indices",
    subtitle = "Alpha-diversity comparison across STEC strains",
    theme    = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  )

save_figure(fig7_diversity, "Fig07_Alpha_Diversity_Indices.png",
            width = 10, height = 12)

cat("\n  All individual figures generated.\n\n")


# ==============================================================================
# SECTION 9 – COMPOSITE MULTI-PANEL FIGURE
# ==============================================================================

cat("10. Assembling composite multi-panel figure...\n")

# Compact versions reduce font clutter in small sub-panels
# IMPORTANT: fig1_compact needs extra top margin (0.25 vs 0.15) because geom_text
# labels with vjust = -0.5 extend above bars. In standalone Fig01 at 10×8 inches,
# 15% is sufficient, but in the 2×2 composite at 12×10 inches with 4 panels,
# the reduced subplot height causes clipping of the "4281 (37.4%)" label on the
# Shell bar. Increasing to 25% ensures visibility in all layouts.
fig1_compact <- fig1_distribution +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.25)),  # was 0.15 in standalone version
    labels = scales::comma
  ) +
  theme(plot.subtitle = element_blank(), plot.title = element_text(size = 12))
fig2_compact <- fig2_accumulation +
  theme(plot.subtitle = element_blank(), plot.title = element_text(size = 12))
fig3_compact <- fig3_core_reduction +
  theme(plot.subtitle = element_blank(), plot.title = element_text(size = 12))
fig5_compact <- fig5_singletons +
  theme(plot.subtitle = element_blank(), plot.title = element_text(size = 12))
fig7a_compact <- fig7a_shannon +
  theme(plot.title = element_text(size = 12))

# 3 × 2 layout using patchwork design string
design_layout <- "
AB
CD
EF
"

fig_composite <- patchwork::wrap_plots(
  fig1_compact, fig2_compact,
  fig3_compact, fig5_compact,
  fig7a_compact, patchwork::plot_spacer(),
  design = design_layout
) +
  patchwork::plot_annotation(
    title    = "Pan-genome Analysis – E. coli STEC",
    subtitle = sprintf(
      "%d strains | %s genes total | Open pan-genome (\u03b3 = %.4f)",
      n_strains,
      scales::comma(nrow(gene_pa)),
      gamma_est
    ),
    tag_levels = "A",
    theme = theme(
      plot.title    = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      plot.caption  = element_text(size = 10, hjust = 1)
    )
  )

save_figure(fig_composite, "FigS01_Composite_All_Analyses.png",
            width = 14, height = 16)

# -- 9.1  Article figure for manuscript Section 3.3 ----------------------------
cat("  Generating 2×2 figure for manuscript Section 3.3...\n")

fig_article_s33 <- (fig1_compact + fig2_compact) /
  (fig3_compact + fig5_compact) +
  patchwork::plot_annotation(
    title    = "Pan-genome Structure and Core Genome Analysis",
    subtitle = sprintf(
      "E. coli STEC strains (n = %d) | Total genes: %s | Open pan-genome (\u03b3 = %.4f)",
      n_strains,
      format(nrow(gene_pa), big.mark = " "),
      gamma_est
    ),
    caption = sprintf(
      "Date: %s | Core genome: %s genes (%.1f%%) | Accessory: %.1f%% | Analysis: Roary v3.13.0",
      Sys.Date(),
      format(dplyr::last(accumulation_df[["core"]]), big.mark = " "),
      sum(gene_pa[["category"]] == "Core",    na.rm = TRUE) / nrow(gene_pa) * 100,
      sum(gene_pa[["category"]] %in% c("Shell", "Cloud"), na.rm = TRUE) /
        nrow(gene_pa) * 100
    ),
    tag_levels = "A",
    theme = theme(
      plot.title    = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, margin = margin(b = 10)),
      plot.caption  = element_text(size = 9,  hjust = 0.5, margin = margin(t = 10))
    )
  )

save_figure(fig_article_s33,
            "FigS02_Manuscript_Section3.3_Pangenome_CoreGenome.png",
            width = 12, height = 10)

cat("  Manuscript composite figure (Section 3.3) saved as FigS02.\n\n")


# ==============================================================================
# SECTION 10 – RESULT EXPORT
# ==============================================================================

cat("11. Exporting result tables...\n")

# Pan-genome summary statistics
save_table(stats_pangenome, "Table_Pangenome_Statistics.csv")

# Core genome gene list
core_genes_tbl <- gene_pa |>
  dplyr::filter(category == "Core") |>
  dplyr::select(Gene, Annotation, n_strains_present, pct_strains) |>
  dplyr::arrange(Gene)

save_table(core_genes_tbl, "Table_Core_Genes.csv")

# Singleton genes per strain
save_table(singleton_per_strain, "Table_Singleton_Genes_Per_Strain.csv")

# Pairwise gene-sharing matrix
sharing_df <- as.data.frame(sharing_matrix)
sharing_df <- cbind(Strain = rownames(sharing_df), sharing_df)
save_table(sharing_df, "Table_Gene_Sharing_Matrix.csv")

# Alpha-diversity indices
save_table(diversity_df, "Table_Diversity_Indices.csv")

cat("\n")


# ==============================================================================
# SECTION 11 – FIGURE LEGENDS
# ==============================================================================

cat("12. Writing figure legends...\n")

core_count_val  <- sum(gene_pa[["category"]] == "Core",      na.rm = TRUE)
shell_count_val <- sum(gene_pa[["category"]] == "Shell",     na.rm = TRUE)
cloud_count_val <- sum(gene_pa[["category"]] == "Cloud",     na.rm = TRUE)
core_pct_val    <- round(core_count_val  / nrow(gene_pa) * 100, 1)
shell_pct_val   <- round(shell_count_val / nrow(gene_pa) * 100, 1)
cloud_pct_val   <- round(cloud_count_val / nrow(gene_pa) * 100, 1)
core_final_val  <- dplyr::last(accumulation_df[["core"]])
total_singles   <- sum(singleton_per_strain[["Unique_genes"]])

legends <- sprintf(
  "
================================================================================
FIGURE LEGENDS – PAN-GENOME ANALYSIS
================================================================================
Project  : Pan-genome Analysis of E. coli Shiga Toxin-Producing (STEC)
Date     : %s
================================================================================

FIGURE 1: PAN-GENOME CATEGORY DISTRIBUTION
Bar chart showing the distribution of orthologous gene clusters across four
categories: Core (99–100%% of strains), Soft Core (95–99%%), Shell (15–95%%),
and Cloud (< 15%%). The pan-genome of %d STEC strains encompasses %s genes.
Core genome: %s genes (%.1f%%); Shell: %s genes (%.1f%%);
Cloud: %s genes (%.1f%%).

FIGURE 2: PAN-GENOME ACCUMULATION CURVE
Cumulative gene count as genomes are added sequentially in random order.
Red solid line: observed data; blue dashed line: Heaps' Law fit
(N = %.2f × n^%.4f, γ = %.4f). γ < 1 confirms an open pan-genome.

FIGURE 3: CORE GENOME REDUCTION CURVE
Core genome size as a function of the number of strains included. The
plateau at %s genes (%d strains) defines the conserved gene repertoire.

FIGURE 4: GENE PRESENCE/ABSENCE HEATMAP (REPRESENTATIVE SAMPLE)
Binary heatmap (black = present, white = absent) for a stratified random
sample of genes (≤ 500 per category). Rows and columns are ordered by
hierarchical clustering (binary distance, average linkage).

FIGURE 5: STRAIN-SPECIFIC GENES (SINGLETONS)
Horizontal bar chart of genes present exclusively in each strain
(total singletons: %d; range: %d – %d).

FIGURE 6: UPSET INTERSECTION PLOT
Top 30 pairwise intersections of gene sets across strains, ordered by
frequency. Vertical bars: intersection size; horizontal bars: strain set size.

FIGURE 7: ALPHA-DIVERSITY INDICES
(A) Shannon and (B) Simpson diversity indices per strain.
Mean Shannon: %.3f; Mean Simpson: %.3f.

================================================================================
",
  Sys.Date(),
  n_strains,
  format(nrow(gene_pa), big.mark = " "),
  format(core_count_val,  big.mark = " "), core_pct_val,
  format(shell_count_val, big.mark = " "), shell_pct_val,
  format(cloud_count_val, big.mark = " "), cloud_pct_val,
  kappa_est, gamma_est, gamma_est,
  format(core_final_val, big.mark = " "),
  n_strains,
  total_singles,
  min(singleton_per_strain[["Unique_genes"]]),
  max(singleton_per_strain[["Unique_genes"]]),
  mean(diversity_df[["Shannon"]]),
  mean(diversity_df[["Simpson"]])
)

writeLines(legends, file.path(output_dir, "Documentacao", "Figure_Legends.txt"))
cat("  Figure legends saved: Documentacao/Figure_Legends.txt\n\n")


# ==============================================================================
# SECTION 12 – REPRODUCIBILITY DOCUMENTATION
# ==============================================================================

cat("13. Writing reproducibility documentation...\n")

# -- 12.1  session_info --------------------------------------------------------
sink(file.path(output_dir, "Documentacao", "SessionInfo.txt"))
cat(sprintf("Analysis date     : %s\n",  Sys.time()))
cat(sprintf("Working directory : %s\n\n", getwd()))
print(sessionInfo())
sink()

# -- 12.2  README  -------------------------------------------------------------
readme <- sprintf(
  "
================================================================================
PAN-GENOME ANALYSIS – E. coli STEC
================================================================================
Script          : 05_Prokka_Roary_Analysis_revised.R
Analysis date   : %s
R version       : %s
Bioconductor    : %s
================================================================================

SUMMARY
-------
Total genes (pan-genome) : %s
Strains analysed         : %d

Pan-genome composition:
  Core (99–100%%)   : %s genes (%.1f%%)
  Soft Core (95–99%%): %s genes (%.1f%%)
  Shell (15–95%%)   : %s genes (%.1f%%)
  Cloud (< 15%%)    : %s genes (%.1f%%)

Heaps' Law model:
  kappa (κ) : %.2f
  gamma (γ) : %.4f
  Type      : %s pan-genome

================================================================================
OUTPUT STRUCTURE
================================================================================

Graficos/
  Fig01_Pangenome_Category_Distribution.png
  Fig02_Pangenome_Accumulation_Curve.png
  Fig03_Core_Genome_Reduction_Curve.png
  Fig04_Gene_Presence_Absence_Heatmap.png
  Fig05_Singleton_Genes_Per_Strain.png
  Fig06_UpSet_Gene_Intersections.png
  Fig07_Alpha_Diversity_Indices.png
  FigS01_Composite_All_Analyses.png        (supplementary – all panels)
  FigS02_Manuscript_Section3.3_Pangenome_CoreGenome.png  (2x2 article figure)

Tabelas/
  Table_Pangenome_Statistics.csv
  Table_Core_Genes.csv
  Table_Singleton_Genes_Per_Strain.csv
  Table_Gene_Sharing_Matrix.csv
  Table_Diversity_Indices.csv

Documentacao/
  Figure_Legends.txt
  SessionInfo.txt
  README.txt  (this file)

================================================================================
CITATION
================================================================================
[Author] et al. (2026) [Title]. [Journal]. DOI: [10.xxxx/xxxxx]

Software:
  Seemann T (2014) Prokka. Bioinformatics 30(14):2068.
    doi:10.1093/bioinformatics/btu153
  Page AJ et al. (2015) Roary. Bioinformatics 31(22):3691.
    doi:10.1093/bioinformatics/btv421

================================================================================
",
  Sys.Date(),
  R.version.string,
  as.character(BiocManager::version()),
  format(nrow(gene_pa), big.mark = " "),
  n_strains,
  format(core_count_val,  big.mark = " "), core_pct_val,
  format(sum(gene_pa[["category"]] == "Soft Core", na.rm = TRUE), big.mark = " "),
  round(sum(gene_pa[["category"]] == "Soft Core", na.rm = TRUE) / nrow(gene_pa) * 100, 1),
  format(shell_count_val, big.mark = " "), shell_pct_val,
  format(cloud_count_val, big.mark = " "), cloud_pct_val,
  kappa_est, gamma_est,
  ifelse(gamma_est < 1, "OPEN", "CLOSED")
)

writeLines(readme, file.path(output_dir, "Documentacao", "README.txt"))
cat("  README saved: Documentacao/README.txt\n\n")


# ==============================================================================
# SECTION 13 – FINALISATION
# ==============================================================================

execution_end      <- Sys.time()
execution_duration <- as.numeric(difftime(execution_end, execution_start,
                                          units = "mins"))

cat("\n")
cat("==============================================================================\n")
cat("                       ANALYSIS COMPLETE                                     \n")
cat("==============================================================================\n\n")
cat("Output summary:\n")
cat("  9 figures        (PNG, 600 DPI) in Graficos/\n")
cat("  5 result tables  (CSV)          in Tabelas/\n")
cat("  3 documentation files           in Documentacao/\n\n")
cat(sprintf("All outputs saved to : %s\n", output_dir))
cat(sprintf("Total execution time : %.2f minutes\n\n", execution_duration))
cat("==============================================================================\n")
cat("  SessionInfo and README are available in Documentacao/ for reproducibility. \n")
cat("==============================================================================\n\n")

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
