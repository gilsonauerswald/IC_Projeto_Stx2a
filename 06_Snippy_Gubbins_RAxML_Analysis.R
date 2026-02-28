#!/usr/bin/env Rscript
# ==============================================================================
# Script: 06_snippy_gubbins_raxml_analysis.R
#
# Purpose
#   Generate publication-quality phylogenomics figures from a Snippy + Gubbins +
#   RAxML workflow, including:
#     (i) Maximum-likelihood tree (rectangular; “preview-matched” styling)
#     (ii) Maximum-likelihood tree (circular; compact)
#     (iii) Pairwise SNP distance heatmap (post-Gubbins)
#     (iv) Quantification of recombination filtering impact (pre vs post Gubbins)
#     (v) Composite multi-panel figure (A–D)
#
# Manuscript mapping
#   Update the mapping below to match your manuscript before submission:
#     - Fig_Phylogenetic_Rectangular.png  -> Figure ?? (panel/section ??)
#     - Fig_Phylogenetic_Tree.png         -> Figure ?? (panel/section ??)
#     - Fig_SNP_Distance_Heatmap.png      -> Figure ?? (panel/section ??)
#     - Fig_Recombination_Impact.png      -> Figure ?? (panel/section ??)
#     - Fig_Composite_Panel.png           -> Figure ?? (panel/section ??)
#
# Inputs (expected in the working directory unless overridden below)
#   - Core alignment FASTA (Snippy: “clean_full_aln”)
#   - Filtered polymorphic sites FASTA (Gubbins output)
#   - RAxML bipartitions tree file
#   - Metadata CSV (required columns: Patogeno, Nome_Real)
#
# Outputs (written to out_dir; PNG at 600 dpi by default)
#   - Fig_Phylogenetic_Rectangular.png
#   - Fig_Phylogenetic_Tree.png
#   - Fig_SNP_Distance_Heatmap.png
#   - Fig_Recombination_Impact.png
#   - Fig_Composite_Panel.png
#   - run_metadata.txt         (runtime parameters + file paths)
#   - session_info.txt         (full sessionInfo() for provenance)
#
# Computational environment (authoritative; update only if you truly changed it)
#   - R:           4.5.2
#   - Bioconductor: 3.22
#
# Reproducibility notes
#   - The script is deterministic under fixed inputs. A seed is set defensively.
#   - Package auto-install is disabled by default to avoid uncontrolled versions.
#     If you want auto-install, set INSTALL_MISSING_PACKAGES=1 (see below).
#   - The script records sessionInfo() and key parameters to output files.
#
# Licensing / attribution
#   - Define a license before publication (MIT/BSD/Apache are typical for code).
#   - Add authorship/contact and a preferred citation once the repository is public.
#
# ==============================================================================

options(stringsAsFactors = FALSE, warn = 1)
set.seed(1)

# ----------------------------
# 0) CONFIGURATION
# ----------------------------

# You may override a few fields via environment variables:
#   OUT_DIR                    (default ".")
#   PNG_DPI                    (default 600)
#   INSTALL_MISSING_PACKAGES   (0/1; default 0)
#
# Rationale: environment variables are CI-friendly and avoid extra dependencies
# (e.g., optparse) for a “supplementary script” use case.

config <- list(
  # Input files (Galaxy pipeline outputs)
  file_core_aln = "Galaxy22-[snippy-clean_full_aln on dataset 21 cleaned core alignment].fasta",
  file_gubbins_aln = "Galaxy23-[Gubbins on dataset 22 Filtered Polymorphic Sites fasta].fasta",
  file_bipart = "Galaxy28-[RAxML on dataset 23_ Bipartitions].txt",
  file_metadata = "Metadados.csv",

  # Output
  out_dir = Sys.getenv("OUT_DIR", unset = "."),
  export_png = TRUE,
  png_dpi = as.integer(Sys.getenv("PNG_DPI", unset = "600")),

  # Figure sizes (inches). Defaults target a typical double-column width.
  fig_size = list(
    phylo_rect = c(width = 12.0, height = 7.1),
    phylo_circ = c(width = 10.0, height = 8.0),
    heatmap = c(width = 10.0, height = 9.0),
    recomb = c(width = 15.0, height = 5.5),
    panel = c(width = 12.0, height = 14.5)
  ),

  # Typography
  base_family = "sans",
  base_size = 10,

  # Thresholds
  bootstrap_min_label = 70,

  # Behavior
  install_missing_packages = identical(Sys.getenv("INSTALL_MISSING_PACKAGES", unset = "0"), "1")
)

# Curated, colorblind-safe palette (Okabe–Ito inspired)
palette <- list(
  group = c(
    "Commensal (K-12)" = "#999999",
    "O157:H7" = "#E69F00",
    "Big Six" = "#56B4E9",
    "Hybrid (O104:H4)" = "#CC79A7",
    "Unknown" = "#666666"
  ),
  stx2a = c(
    "Present" = "#D55E00",
    "Absent" = "#0072B2"
  )
)

# stx2a carrier status (from ABRicate results / manuscript)
stx2a_carriers <- c(
  "O157_H7_Sakai_fna",     # reference (may appear as "Reference" in some tree files)
  "O145_H28_RM12581_fna",
  "O121_H19_FWSEC0006_fna",
  "O104_H4_2011C_3493_fna",
  "O103_H2_12009_fna"
)

# Serogroup classification used for annotation
serogroup_info <- data.frame(
  tip_id = c(
    "K_12_MG1655_fna", "Reference", "O157_H7_EDL933_fna",
    "O145_H28_RM12581_fna", "O121_H19_FWSEC0006_fna",
    "O111_H8_7_58_72A_fna", "O104_H4_2011C_3493_fna",
    "O103_H2_12009_fna", "O45_H2_FWSEC0003_fna", "O26_H11_11368_fna"
  ),
  group = c(
    "Commensal (K-12)", "O157:H7", "O157:H7",
    "Big Six", "Big Six",
    "Big Six", "Hybrid (O104:H4)",
    "Big Six", "Big Six", "Big Six"
  )
)

# ----------------------------
# 1) DEPENDENCIES (controlled)
# ----------------------------

#' Check and load a package with controlled installation behavior.
#'
#' This helper avoids silent version drift by default. In a manuscript workflow,
#' the preferred approach is to capture an environment (e.g., renv.lock or a
#' container) rather than auto-install at runtime.
#'
#' @param pkg Package name.
#' @param source One of "CRAN" or "Bioconductor".
#' @param install_missing If TRUE, attempt installation when missing.
#' @return Invisibly returns TRUE when the package is available and loaded.
#' @noRd
require_pkg <- function(pkg, source = c("CRAN", "Bioconductor"), install_missing = FALSE) {
  source <- match.arg(source)

  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!isTRUE(install_missing)) {
      stop(
        "Required package not installed: ", pkg, "\n",
        "Install it and re-run. ",
        if (source == "Bioconductor") {
          "Example: BiocManager::install(\"" %+% pkg %+% "\")"
        } else {
          "Example: install.packages(\"" %+% pkg %+% "\")"
        },
        call. = FALSE
      )
    }

    message("[INSTALL] ", pkg, " (", source, ")")
    if (source == "Bioconductor") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)
      }
      BiocManager::install(pkg, update = FALSE, ask = FALSE, quiet = TRUE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
    }
  }

  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  invisible(TRUE)
}

#' Load all required packages.
#'
#' @param install_missing_packages If TRUE, missing packages are installed.
#' @return NULL (called for side effects).
#' @noRd
load_packages <- function(install_missing_packages = FALSE) {
  # Bioconductor
  require_pkg("ggtree", "Bioconductor", install_missing_packages)
  require_pkg("treeio", "Bioconductor", install_missing_packages)
  require_pkg("Biostrings", "Bioconductor", install_missing_packages)

  # CRAN
  require_pkg("ape", "CRAN", install_missing_packages)
  require_pkg("ggplot2", "CRAN", install_missing_packages)
  require_pkg("ggnewscale", "CRAN", install_missing_packages)
  require_pkg("viridis", "CRAN", install_missing_packages)
  require_pkg("cowplot", "CRAN", install_missing_packages)
  require_pkg("scales", "CRAN", install_missing_packages)
  require_pkg("ragg", "CRAN", install_missing_packages)

  # Tidyverse-adjacent utilities (modern replacements for reshape2/read.csv patterns)
  require_pkg("readr", "CRAN", install_missing_packages)
  require_pkg("dplyr", "CRAN", install_missing_packages)
  require_pkg("tibble", "CRAN", install_missing_packages)
}

# Small infix to avoid paste0 noise in error strings.
`%+%` <- function(x, y) paste0(x, y)

# ----------------------------
# 2) I/O + VALIDATION HELPERS
# ----------------------------

#' Stop with a clear message if any required inputs are missing.
#'
#' @param paths Character vector of file paths.
#' @return NULL (called for side effects).
#' @noRd
stop_if_missing_files <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    stop(
      "Missing input files:\n  ",
      paste(missing, collapse = "\n  "),
      call. = FALSE
    )
  }
}

#' Read metadata CSV robustly (UTF-8 preferred; fallback to latin1).
#'
#' @param path Path to CSV.
#' @return A data.frame with strings preserved (no factors).
#' @noRd
read_metadata_csv <- function(path) {
  out <- tryCatch(
    readr::read_csv(path, show_col_types = FALSE, locale = readr::locale(encoding = "UTF-8")) |>
      as.data.frame(),
    error = function(e) NULL
  )

  if (is.null(out)) {
    out <- readr::read_csv(path, show_col_types = FALSE, locale = readr::locale(encoding = "latin1")) |>
      as.data.frame()
  }

  out
}

#' Rename row/column names of a matrix using a named mapping.
#'
#' @param mat A numeric matrix.
#' @param mapping Named character vector; names are old IDs, values are display labels.
#' @return The same matrix with updated dimnames.
#' @noRd
rename_matrix_dimnames <- function(mat, mapping) {
  rn <- rownames(mat)
  cn <- colnames(mat)

  if (!is.null(rn)) {
    rownames(mat) <- ifelse(rn %in% names(mapping), mapping[rn], rn)
  }
  if (!is.null(cn)) {
    colnames(mat) <- ifelse(cn %in% names(mapping), mapping[cn], cn)
  }

  mat
}

#' Extract bootstrap values from a rooted phylo object into a data.frame.
#'
#' @param tree An ape::phylo object.
#' @return data.frame with columns: node (integer), bootstrap (numeric).
#' @noRd
get_bootstrap_df <- function(tree) {
  n_tips <- ape::Ntip(tree)
  n_nodes <- tree$Nnode
  node_ids <- (n_tips + 1):(n_tips + n_nodes)

  raw_labels <- tree$node.label
  if (is.null(raw_labels)) {
    raw_labels <- rep(NA_character_, n_nodes)
  }

  # root(resolve.root = TRUE) can add an unlabeled node; pad/trim to match Nnode
  if (length(raw_labels) < n_nodes) {
    raw_labels <- c(raw_labels, rep(NA_character_, n_nodes - length(raw_labels)))
  }
  if (length(raw_labels) > n_nodes) {
    raw_labels <- raw_labels[seq_len(n_nodes)]
  }

  tibble::tibble(
    node = node_ids,
    bootstrap = suppressWarnings(as.numeric(raw_labels))
  ) |>
    as.data.frame()
}

#' Save a ggplot/cowplot object to PNG using ragg at fixed DPI.
#'
#' @param plot A plot object.
#' @param filename_base Output base name (without extension).
#' @param width Width in inches.
#' @param height Height in inches.
#' @param out_dir Output directory.
#' @param dpi Resolution (dots per inch).
#' @return Path to the saved file.
#' @noRd
save_figure_png <- function(plot, filename_base, width, height, out_dir, dpi) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  out_path <- file.path(out_dir, paste0(filename_base, ".png"))
  ragg::agg_png(
    filename = out_path,
    width = width,
    height = height,
    units = "in",
    res = dpi,
    background = "white"
  )
  print(plot)
  dev.off()

  out_path
}

#' Best-effort retrieval of the current script path (when executed via Rscript).
#'
#' @return Character scalar path, or NULL if not available.
#' @noRd
get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    sub("^--file=", "", file_arg)
  } else {
    NULL
  }
}

#' Write run metadata and session provenance to disk.
#'
#' @param cfg Configuration list.
#' @return NULL.
#' @noRd
write_provenance_files <- function(cfg) {
  if (!dir.exists(cfg$out_dir)) {
    dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  meta_path <- file.path(cfg$out_dir, "run_metadata.txt")
  meta_lines <- c(
    paste0("started_at\t", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("script\t", normalizePath(get_script_path() %||% "unknown", winslash = "/", mustWork = FALSE)),
    paste0("file_core_aln\t", cfg$file_core_aln),
    paste0("file_gubbins_aln\t", cfg$file_gubbins_aln),
    paste0("file_bipart\t", cfg$file_bipart),
    paste0("file_metadata\t", cfg$file_metadata),
    paste0("out_dir\t", cfg$out_dir),
    paste0("png_dpi\t", cfg$png_dpi),
    paste0("bootstrap_min_label\t", cfg$bootstrap_min_label),
    paste0("install_missing_packages\t", cfg$install_missing_packages)
  )
  writeLines(meta_lines, con = meta_path)

  sess_path <- file.path(cfg$out_dir, "session_info.txt")
  writeLines(capture.output(sessionInfo()), con = sess_path)

  invisible(NULL)
}

# Provide `%||%` to avoid extra dependencies.
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

# ----------------------------
# 3) PLOTTING THEME
# ----------------------------

#' Publication theme for non-tree plots.
#'
#' @param base_size Base font size.
#' @param base_family Font family.
#' @return A ggplot2 theme object.
#' @noRd
theme_publication <- function(base_size = config$base_size, base_family = config$base_family) {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      text = ggplot2::element_text(color = "black"),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = base_size, color = "gray35", hjust = 0.5),
      axis.title = ggplot2::element_text(face = "bold"),
      panel.grid.major = ggplot2::element_line(color = "gray92", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", linewidth = 0.6),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.background = ggplot2::element_rect(fill = "white", color = "gray70", linewidth = 0.3)
    )
}

# ----------------------------
# 4) DATA LOADING + ANNOTATION
# ----------------------------

#' Read an RAxML tree and root it on an outgroup.
#'
#' @param file_bipart Path to RAxML bipartitions file.
#' @param outgroup Tip label to root on.
#' @return Rooted ape::phylo object.
#' @noRd
read_tree_rooted <- function(file_bipart, outgroup = "K_12_MG1655_fna") {
  tree <- ape::read.tree(file_bipart)
  ape::root(tree, outgroup = outgroup, resolve.root = TRUE)
}

#' Build a tip label mapping from metadata.
#'
#' @param metadata Metadata data.frame with columns Patogeno and Nome_Real.
#' @return Named character vector (old -> display).
#' @noRd
build_name_map <- function(metadata) {
  required_cols <- c("Patogeno", "Nome_Real")
  if (!all(required_cols %in% names(metadata))) {
    stop(
      "Metadados.csv must contain columns: ",
      paste(required_cols, collapse = ", "),
      call. = FALSE
    )
  }

  nm <- setNames(metadata$Nome_Real, metadata$Patogeno)

  # Some tree files use "Reference" instead of the original Sakai ID.
  if ("O157_H7_Sakai_fna" %in% names(nm)) {
    nm["Reference"] <- nm[["O157_H7_Sakai_fna"]]
  }

  nm
}

#' Rename tree tip labels using a mapping.
#'
#' @param tree A phylo object.
#' @param mapping Named character vector (old -> display).
#' @return Updated phylo object.
#' @noRd
rename_tree_tips <- function(tree, mapping) {
  tree$tip.label <- ifelse(
    tree$tip.label %in% names(mapping),
    mapping[tree$tip.label],
    tree$tip.label
  )
  tree
}

#' Build per-tip annotation used by ggtree.
#'
#' @param tree_display Rooted tree with display labels.
#' @param name_map Named character vector mapping original IDs -> display labels.
#' @return data.frame with columns: label, original_id, group, stx2a.
#' @noRd
build_tip_annotation <- function(tree_display, name_map) {
  tip_df <- tibble::tibble(label = tree_display$tip.label)

  # Reverse mapping: display -> original.
  rev_map <- setNames(names(name_map), name_map)

  tip_df <- tip_df |>
    dplyr::mutate(
      original_id = vapply(
        label,
        function(x) if (x %in% names(rev_map)) rev_map[[x]] else x,
        character(1)
      )
    )

  infer_group <- function(orig_id, disp_label) {
    if (grepl("K[-_ ]?12|MG1655", orig_id, ignore.case = TRUE) ||
        grepl("K[-_ ]?12|MG1655", disp_label, ignore.case = TRUE)) {
      return("Commensal (K-12)")
    }

    if (grepl("^O157", orig_id) || grepl("^O157", disp_label) || grepl("O157:H7", disp_label, fixed = TRUE)) {
      return("O157:H7")
    }

    if (grepl("^O104", orig_id) || grepl("O104", disp_label)) {
      return("Hybrid (O104:H4)")
    }

    if (grepl("^O(26|45|103|111|121|145)", orig_id) || grepl("O(26|45|103|111|121|145)", disp_label)) {
      return("Big Six")
    }

    "Unknown"
  }

  tip_df <- tip_df |>
    dplyr::rowwise() |>
    dplyr::mutate(
      group = {
        idx <- which(serogroup_info$tip_id == original_id)
        if (length(idx) > 0) {
          serogroup_info$group[idx][1]
        } else {
          infer_group(original_id, label)
        }
      }
    ) |>
    dplyr::ungroup()

  stx2a_tree_ids <- unique(c("Reference", stx2a_carriers))
  tip_df <- tip_df |>
    dplyr::mutate(stx2a = ifelse(original_id %in% stx2a_tree_ids, "Present", "Absent"))

  as.data.frame(tip_df)
}

#' Read alignments and compute SNP-distance matrices (pre and post Gubbins).
#'
#' Implementation details:
#'   - Reads FASTA as DNAbin matrices (ape::read.dna with as.matrix=TRUE).
#'   - Attempts polymorphic-site extraction via ape::seg.sites; falls back to a
#'     Biostrings consensus-matrix based extractor when needed.
#'
#' @param file_core_aln Path to core alignment FASTA.
#' @param file_gubbins_aln Path to Gubbins-filtered polymorphic sites FASTA.
#' @param name_map Mapping of IDs -> display labels (applied to matrix dimnames).
#' @return List containing site counts and distance matrices.
#' @noRd
read_alignments_and_distances <- function(file_core_aln, file_gubbins_aln, name_map) {
  dna_core <- ape::read.dna(file_core_aln, format = "fasta", as.matrix = TRUE)
  dna_gub <- ape::read.dna(file_gubbins_aln, format = "fasta", as.matrix = TRUE)

  if (!is.matrix(dna_core)) dna_core <- as.matrix(dna_core)
  if (!is.matrix(dna_gub)) dna_gub <- as.matrix(dna_gub)

  safe_ncol <- function(x) {
    d <- dim(x)
    if (is.null(d) || length(d) < 2) return(0L)
    as.integer(d[2])
  }

  extract_poly_from_fasta <- function(fasta_path, dna_matrix) {
    dss <- Biostrings::readDNAStringSet(fasta_path)

    # Use only A/C/G/T for polymorphism; treat N/gaps/ambiguous as missing.
    cm <- Biostrings::consensusMatrix(dss, as.prob = FALSE)
    bases <- intersect(c("A", "C", "G", "T"), rownames(cm))
    if (length(bases) < 4) {
      stop(
        "Consensus matrix missing expected base rows (A/C/G/T). Check FASTA encoding: ",
        fasta_path,
        call. = FALSE
      )
    }

    cm4 <- cm[bases, , drop = FALSE]
    n_present <- colSums(cm4 > 0)
    poly_pos <- which(n_present > 1)
    if (length(poly_pos) == 0L) return(NULL)

    dna_matrix[, poly_pos, drop = FALSE]
  }

  dna_core_poly <- tryCatch(
    ape::seg.sites(dna_core, strict = FALSE),
    error = function(e) NULL
  )

  if (is.null(dna_core_poly) || safe_ncol(dna_core_poly) == 0L) {
    dna_core_poly <- extract_poly_from_fasta(file_core_aln, dna_core)
  }

  if (is.null(dna_core_poly) || safe_ncol(dna_core_poly) == 0L) {
    stop(
      "Unable to extract polymorphic sites from core alignment. ",
      "Check that the core alignment contains variation and is a proper alignment matrix.",
      call. = FALSE
    )
  }

  if (safe_ncol(dna_gub) == 0L) {
    stop("Gubbins FASTA appears empty or unparsable: ", file_gubbins_aln, call. = FALSE)
  }

  n_poly_input <- safe_ncol(dna_core_poly)
  n_poly_output <- safe_ncol(dna_gub)

  sites_removed <- max(0L, n_poly_input - n_poly_output)
  pct_removed <- if (n_poly_input > 0L) round(100 * sites_removed / n_poly_input, 2) else NA_real_

  dist_core <- as.matrix(ape::dist.dna(dna_core_poly, model = "N", pairwise.deletion = TRUE))
  dist_gub <- as.matrix(ape::dist.dna(dna_gub, model = "N", pairwise.deletion = TRUE))

  dist_core <- rename_matrix_dimnames(dist_core, name_map)
  dist_gub <- rename_matrix_dimnames(dist_gub, name_map)

  # Ensure consistent strain order between matrices for comparisons.
  common <- rownames(dist_gub)[rownames(dist_gub) %in% rownames(dist_core)]
  if (length(common) < 2) {
    stop(
      "Recombination impact: strain sets differ between core alignment and Gubbins output (insufficient overlap).",
      call. = FALSE
    )
  }

  dist_core <- dist_core[common, common, drop = FALSE]
  dist_gub <- dist_gub[common, common, drop = FALSE]

  list(
    dna_core_poly = dna_core_poly,
    dna_gubbins = dna_gub,
    n_poly_input = n_poly_input,
    n_poly_output = n_poly_output,
    sites_removed = sites_removed,
    pct_removed = pct_removed,
    dist_core = dist_core,
    dist_gubbins = dist_gub
  )
}

# ----------------------------
# 5) PLOTTING FUNCTIONS
# ----------------------------

#' Rectangular ML tree (preview-matched styling).
#'
#' @param tree_display Rooted tree with display labels.
#' @param tip_data Tip annotation data.frame.
#' @param boot_df Bootstrap data.frame from get_bootstrap_df().
#' @param cfg Configuration list (for bootstrap threshold and typography).
#' @return ggplot object.
#' @noRd
plot_tree_rectangular_preview <- function(tree_display, tip_data, boot_df, cfg = config) {
  # Precompute bootstrap labels outside aes() to avoid evaluation-scope issues.
  boot_df2 <- boot_df |>
    dplyr::mutate(
      bootstrap_label = dplyr::if_else(
        !is.na(bootstrap) & bootstrap >= cfg$bootstrap_min_label,
        bootstrap,
        NA_real_
      )
    )

  p <- ggtree::ggtree(
    tree_display,
    layout = "rectangular",
    ladderize = TRUE,
    linewidth = 1.0,
    color = "black"
  ) %<+% tip_data %<+% boot_df2

  p <- p +
    ggtree::theme_tree2() +
    ggplot2::labs(
      title = "Core Genome SNP-Based ML Phylogeny of STEC Strains",
      subtitle = "RAxML (GTR+GAMMA) | 1,000 Bootstrap Replicates | Rooted on K-12 MG1655",
      x = "Substitutions per site",
      y = "taxa"
    ) +
    ggtree::geom_tiplab(
      ggplot2::aes(color = group),
      size = 4.2,
      offset = 0.010,
      fontface = "bold",
      show.legend = FALSE
    ) +
    ggtree::geom_tippoint(
      ggplot2::aes(color = group, shape = stx2a),
      size = 3.2,
      stroke = 0.4
    ) +
    ggtree::geom_nodelab(
      ggplot2::aes(
        subset = !isTip & !is.na(bootstrap_label),
        label = bootstrap_label
      ),
      size = 3.2,
      color = "black",
      fontface = "plain",
      hjust = -0.25,
      vjust = -0.55
    ) +
    ggplot2::scale_color_manual(
      name = "Serogroup / stx2a",
      values = palette$group,
      breaks = c("O157:H7", "Big Six", "Hybrid (O104:H4)", "Commensal (K-12)")
    ) +
    ggplot2::scale_shape_manual(
      name = NULL,
      values = c("Present" = 17, "Absent" = 16),
      breaks = c("Present", "Absent"),
      labels = c("stx2a Present", "stx2a Absent")
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(order = 1, override.aes = list(shape = 16, size = 3.5)),
      shape = ggplot2::guide_legend(
        order = 2,
        override.aes = list(
          color = c(palette$stx2a[["Present"]], palette$stx2a[["Absent"]]),
          size = 3.5
        )
      )
    ) +
    ggplot2::theme(
      text = ggplot2::element_text(family = cfg$base_family),
      plot.title = ggplot2::element_text(face = "bold", size = 18, hjust = 0.5, margin = ggplot2::margin(b = 6)),
      plot.subtitle = ggplot2::element_text(face = "bold", size = 14, hjust = 0.5, margin = ggplot2::margin(b = 12)),
      axis.title = ggplot2::element_text(face = "plain", size = 13),
      axis.text = ggplot2::element_text(size = 12, color = "black"),
      panel.grid = ggplot2::element_blank(),
      legend.position = c(0.83, 0.40),
      legend.justification = c(0, 0),
      legend.background = ggplot2::element_rect(fill = "white", color = "gray70", linewidth = 0.3),
      legend.title = ggplot2::element_text(face = "bold", size = 12),
      legend.text = ggplot2::element_text(size = 11),
      plot.margin = ggplot2::margin(10, 30, 10, 10)
    )

  xmax <- max(ape::node.depth.edgelength(tree_display), na.rm = TRUE) * 1.20
  p + ggplot2::coord_cartesian(xlim = c(NA, xmax), clip = "off")
}

#' Circular ML tree (compact).
#'
#' @param tree_display Rooted tree with display labels.
#' @param tip_data Tip annotation data.frame.
#' @param boot_df Bootstrap data.frame.
#' @param n_sites Number of recombination-free SNP sites.
#' @param cfg Configuration list.
#' @return ggplot object.
#' @noRd
plot_tree_circular_compact <- function(tree_display, tip_data, boot_df, n_sites, cfg = config) {
  # Precompute bootstrap labels outside aes() to avoid evaluation-scope issues.
  boot_df2 <- boot_df |>
    dplyr::mutate(
      bootstrap_label = dplyr::if_else(
        !is.na(bootstrap) & bootstrap >= cfg$bootstrap_min_label,
        bootstrap,
        NA_real_
      )
    )

  p <- ggtree::ggtree(tree_display, layout = "circular", linewidth = 0.85, color = "gray25") %<+% tip_data %<+% boot_df2

  p <- p +
    ggtree::geom_tiplab(
      ggplot2::aes(color = group),
      size = 3.0,
      offset = 0.015,
      fontface = "bold",
      show.legend = FALSE
    ) +
    ggtree::geom_tippoint(ggplot2::aes(color = group, shape = stx2a), size = 2.6, stroke = 0.35) +
    ggtree::geom_nodelab(
      ggplot2::aes(
        subset = !isTip & !is.na(bootstrap_label),
        label = bootstrap_label
      ),
      size = 2.3,
      color = "gray25",
      hjust = -0.25,
      vjust = -0.45
    ) +
    ggplot2::scale_color_manual(
      name = "Serogroup / stx2a",
      values = palette$group,
      breaks = c("O157:H7", "Big Six", "Hybrid (O104:H4)", "Commensal (K-12)")
    ) +
    ggplot2::scale_shape_manual(
      name = NULL,
      values = c("Present" = 17, "Absent" = 16),
      breaks = c("Present", "Absent"),
      labels = c("stx2a Present", "stx2a Absent")
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(order = 1, override.aes = list(shape = 16, size = 3.5)),
      shape = ggplot2::guide_legend(
        order = 2,
        override.aes = list(
          color = c(palette$stx2a[["Present"]], palette$stx2a[["Absent"]]),
          size = 3.5
        )
      )
    ) +
    ggplot2::labs(
      title = "Maximum Likelihood Phylogenetic Tree of STEC Strains",
      subtitle = paste0(
        "GTR+GAMMA | 1,000 bootstraps | ",
        scales::comma(n_sites),
        " recombination-free SNP sites"
      )
    ) +
    ggplot2::theme(
      text = ggplot2::element_text(family = cfg$base_family),
      plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 9, color = "gray40", hjust = 0.5),
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold", size = 12),
      legend.text = ggplot2::element_text(size = 11),
      legend.background = ggplot2::element_rect(fill = "white", color = "gray70", linewidth = 0.3),
      plot.margin = ggplot2::margin(14, 14, 14, 14)
    )

  max_x <- max(p$data$x, na.rm = TRUE)
  if (is.finite(max_x) && max_x > 0) {
    p <- p + ggplot2::xlim(0, max_x * 1.22)
  }

  p
}

#' SNP distance heatmap (post-Gubbins).
#'
#' @param dist_mat Pairwise SNP distance matrix (numeric).
#' @param tree_display Rooted tree with display labels (used for ordering).
#' @param subtitle_nsites Number of sites for subtitle text.
#' @return ggplot object.
#' @noRd
plot_snp_heatmap <- function(dist_mat, tree_display, subtitle_nsites) {
  dist_df <- as.data.frame(as.table(dist_mat), stringsAsFactors = FALSE) |>
    dplyr::rename(strain1 = Var1, strain2 = Var2, snps = Freq)

  tree_order <- rev(ggtree::get_taxa_name(ggtree::ggtree(tree_display, ladderize = TRUE)))
  dist_df <- dist_df |>
    dplyr::mutate(
      strain1 = factor(strain1, levels = tree_order),
      strain2 = factor(strain2, levels = tree_order)
    )

  y_limits <- rev(tree_order)
  maxv <- max(dist_df$snps, na.rm = TRUE)
  dist_df <- dist_df |>
    dplyr::mutate(label_col = ifelse(snps > maxv * 0.62, "white", "black"))

  ggplot2::ggplot(dist_df, ggplot2::aes(x = strain1, y = strain2, fill = snps)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = scales::comma(snps), color = label_col), size = 2.6) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_viridis_c(
      name = "Pairwise\nSNP Count",
      option = "inferno",
      direction = -1,
      begin = 0.05,
      end = 0.95
    ) +
    ggplot2::scale_y_discrete(limits = y_limits) +
    ggplot2::labs(
      title = "Pairwise SNP Distance Matrix (Recombination-Free Core Genome)",
      subtitle = paste0(scales::comma(subtitle_nsites), " polymorphic sites after Gubbins filtering"),
      x = NULL,
      y = NULL
    ) +
    theme_publication(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 8),
      axis.text.y = ggplot2::element_text(face = "bold", size = 8),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right",
      plot.title = ggplot2::element_text(size = 12),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    ) +
    ggplot2::coord_fixed()
}

#' Multi-panel figure quantifying recombination filtering impact.
#'
#' Panels:
#'   A) Density of pairwise SNP distances before vs after Gubbins
#'   B) Scatter plot: before vs after pairwise SNP distances
#'   C) Bar chart: sites (input, removed, output)
#'
#' @param dist_core Pairwise SNP distance matrix (core polymorphic sites).
#' @param dist_gubbins Pairwise SNP distance matrix (post-Gubbins).
#' @param n_poly_input Count of polymorphic sites in the core alignment.
#' @param n_poly_output Count of polymorphic sites after Gubbins.
#' @param sites_removed Recombinant sites removed by Gubbins.
#' @param pct_removed Percent removed.
#' @param internal_labels If TRUE, panels are labeled A/B/C internally.
#' @param include_main_title If TRUE, adds a main title/subtitle (for standalone figure).
#' @return cowplot object.
#' @noRd
plot_recombination_impact <- function(dist_core,
                                     dist_gubbins,
                                     n_poly_input,
                                     n_poly_output,
                                     sites_removed,
                                     pct_removed,
                                     internal_labels = TRUE,
                                     include_main_title = TRUE) {
  core_df <- as.data.frame(as.table(dist_core), stringsAsFactors = FALSE) |>
    dplyr::rename(strain1 = Var1, strain2 = Var2, snps = Freq) |>
    dplyr::mutate(dataset = "Before Gubbins\n(Core polymorphic sites)")

  gub_df <- as.data.frame(as.table(dist_gubbins), stringsAsFactors = FALSE) |>
    dplyr::rename(strain1 = Var1, strain2 = Var2, snps = Freq) |>
    dplyr::mutate(dataset = "After Gubbins\n(Recombination filtered)")

  combined <- dplyr::bind_rows(core_df, gub_df) |>
    dplyr::filter(strain1 != strain2, is.finite(snps))

  if (length(unique(combined$dataset)) < 2) {
    stop(
      "Recombination impact: missing 'Before Gubbins' or 'After Gubbins' distances. Verify inputs.",
      call. = FALSE
    )
  }

  p_a <- ggplot2::ggplot(combined, ggplot2::aes(x = snps, fill = dataset)) +
    ggplot2::geom_density(alpha = 0.55, color = "gray35", linewidth = 0.5) +
    ggplot2::scale_fill_manual(values = c(
      "Before Gubbins\n(Core polymorphic sites)" = "#E69F00",
      "After Gubbins\n(Recombination filtered)" = "#0072B2"
    )) +
    ggplot2::scale_x_continuous(labels = scales::comma) +
    ggplot2::labs(title = if (isTRUE(internal_labels)) "A" else NULL, x = "Pairwise SNP distance", y = "Density") +
    theme_publication(base_size = 10) +
    ggplot2::theme(
      legend.position = c(0.63, 0.86),
      legend.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0),
      panel.grid.minor = ggplot2::element_blank()
    )

  pairs_df <- tibble::tibble(
    before = dist_core[lower.tri(dist_core)],
    after = dist_gubbins[lower.tri(dist_gubbins)]
  ) |>
    dplyr::filter(is.finite(before), is.finite(after), before > 0)

  if (nrow(pairs_df) == 0) {
    stop(
      "Recombination impact: no valid pairwise comparisons for Panel B. Check polymorphic sites extraction and strain overlap.",
      call. = FALSE
    )
  }

  pairs_df <- pairs_df |>
    dplyr::mutate(reduction_pct = round(100 * (before - after) / before, 1))

  p_b <- ggplot2::ggplot(pairs_df, ggplot2::aes(x = before, y = after)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(color = reduction_pct), size = 3, alpha = 0.9, stroke = 0.25) +
    ggplot2::scale_color_viridis_c(name = "SNP reduction\n(%)", option = "plasma", begin = 0.1, end = 0.9) +
    ggplot2::scale_x_continuous(labels = scales::comma) +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::labs(
      title = if (isTRUE(internal_labels)) "B" else NULL,
      x = "SNP distance before Gubbins",
      y = "SNP distance after Gubbins"
    ) +
    theme_publication(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0))

  summary_stats <- tibble::tibble(
    metric = factor(
      c("Polymorphic\nSites (Input)", "Sites Removed\n(Recombination)", "Clean Sites\n(Output)"),
      levels = c("Polymorphic\nSites (Input)", "Sites Removed\n(Recombination)", "Clean Sites\n(Output)")
    ),
    value = c(n_poly_input, sites_removed, n_poly_output),
    category = c("Input", "Removed", "Output")
  )

  p_c <- ggplot2::ggplot(summary_stats, ggplot2::aes(x = metric, y = value, fill = category)) +
    ggplot2::geom_col(width = 0.65, color = "gray35", linewidth = 0.3) +
    ggplot2::geom_text(
      ggplot2::aes(label = scales::comma(value)),
      vjust = -0.5,
      size = 3.5,
      fontface = "bold"
    ) +
    ggplot2::scale_fill_manual(values = c("Input" = "#56B4E9", "Removed" = "#D55E00", "Output" = "#009E73")) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15)), labels = scales::comma) +
    ggplot2::labs(title = if (isTRUE(internal_labels)) "C" else NULL, x = NULL, y = "Number of sites") +
    theme_publication(base_size = 10) +
    ggplot2::theme(legend.position = "none", plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0))

  if (!isTRUE(internal_labels)) {
    p_a <- p_a + ggplot2::theme(plot.title = ggplot2::element_blank())
    p_b <- p_b + ggplot2::theme(plot.title = ggplot2::element_blank())
    p_c <- p_c + ggplot2::theme(plot.title = ggplot2::element_blank())
  }

  grid <- cowplot::plot_grid(p_a, p_b, p_c, ncol = 3, rel_widths = c(1, 1, 0.95), align = "h", axis = "tb")

  title_main <- "Impact of Recombination Filtering on Phylogenetic Signal"
  title_sub <- paste0(
    "Gubbins removed ",
    scales::comma(sites_removed),
    " recombinant sites (",
    pct_removed,
    "%) from polymorphic core alignment"
  )

  if (isTRUE(include_main_title)) {
    cowplot::ggdraw() +
      cowplot::draw_plot(grid, 0, 0, 1, 0.92) +
      cowplot::draw_label(title_main, x = 0.5, y = 0.97, fontface = "bold", size = 13) +
      cowplot::draw_label(title_sub, x = 0.5, y = 0.935, size = 9, color = "gray40")
  } else {
    grid
  }
}

#' Build a unified publication-style panel from the four module figures.
#'
#' @param p_rect Rectangular tree plot.
#' @param p_circ Circular tree plot.
#' @param p_heat Heatmap plot.
#' @param p_recomb Recombination panel (titles removed variant).
#' @return cowplot object.
#' @noRd
plot_unified_panel <- function(p_rect, p_circ, p_heat, p_recomb) {
  strip_titles <- function(p) p + ggplot2::labs(title = NULL, subtitle = NULL)

  p_rect2 <- strip_titles(p_rect) +
    ggplot2::labs(y = NULL) +
    ggplot2::theme(legend.position = c(0.83, 0.52), plot.margin = ggplot2::margin(6, 10, 6, 6))

  p_circ2 <- strip_titles(p_circ) +
    ggplot2::theme(legend.position = "none", plot.margin = ggplot2::margin(6, 6, 6, 6))

  p_heat2 <- strip_titles(p_heat) +
    ggplot2::theme(plot.margin = ggplot2::margin(6, 6, 6, 6))

  top <- cowplot::plot_grid(
    p_rect2,
    ncol = 1,
    labels = "A",
    label_size = 16,
    label_fontface = "bold",
    label_x = 0.01,
    label_y = 0.99,
    hjust = 0,
    vjust = 1
  )

  mid <- cowplot::plot_grid(
    p_circ2,
    p_heat2,
    ncol = 2,
    rel_widths = c(1.00, 1.15),
    labels = c("B", "C"),
    label_size = 16,
    label_fontface = "bold",
    label_x = 0.01,
    label_y = 0.99,
    hjust = 0,
    vjust = 1
  )

  bot <- cowplot::plot_grid(
    p_recomb,
    ncol = 1,
    labels = "D",
    label_size = 16,
    label_fontface = "bold",
    label_x = 0.01,
    label_y = 0.99,
    hjust = 0,
    vjust = 1
  )

  cowplot::plot_grid(top, mid, bot, ncol = 1, rel_heights = c(1.05, 1.00, 0.88), align = "v", axis = "lr")
}

# ----------------------------
# 6) MAIN
# ----------------------------

#' Main entry point.
#'
#' @param cfg Configuration list.
#' @return Invisibly returns a list of generated file paths.
#' @noRd
main <- function(cfg = config) {
  message("============================================================")
  message("STEC Phylogenomics - Module 4 (Submission-ready; PNG export)")
  message("Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  message("Output dir: ", cfg$out_dir)
  message("============================================================")

  load_packages(install_missing_packages = cfg$install_missing_packages)

  stop_if_missing_files(c(cfg$file_core_aln, cfg$file_gubbins_aln, cfg$file_bipart, cfg$file_metadata))

  # Record provenance early (even if something fails later).
  write_provenance_files(cfg)

  metadata <- read_metadata_csv(cfg$file_metadata)
  name_map <- build_name_map(metadata)

  # Tree
  tree_rooted <- read_tree_rooted(cfg$file_bipart, outgroup = "K_12_MG1655_fna")
  tree_display <- rename_tree_tips(tree_rooted, name_map)
  boot_df <- get_bootstrap_df(tree_display)
  tip_data <- build_tip_annotation(tree_display, name_map)

  # Alignments + distances
  aln <- read_alignments_and_distances(cfg$file_core_aln, cfg$file_gubbins_aln, name_map)

  # Figures
  p_rect <- plot_tree_rectangular_preview(tree_display, tip_data, boot_df, cfg)
  fig_rect <- save_figure_png(
    p_rect,
    "Fig_Phylogenetic_Rectangular",
    width = cfg$fig_size$phylo_rect[["width"]],
    height = cfg$fig_size$phylo_rect[["height"]],
    out_dir = cfg$out_dir,
    dpi = cfg$png_dpi
  )

  p_circ <- plot_tree_circular_compact(tree_display, tip_data, boot_df, n_sites = aln$n_poly_output, cfg = cfg)
  fig_circ <- save_figure_png(
    p_circ,
    "Fig_Phylogenetic_Tree",
    width = cfg$fig_size$phylo_circ[["width"]],
    height = cfg$fig_size$phylo_circ[["height"]],
    out_dir = cfg$out_dir,
    dpi = cfg$png_dpi
  )

  p_heat <- plot_snp_heatmap(aln$dist_gubbins, tree_display, subtitle_nsites = aln$n_poly_output)
  fig_heat <- save_figure_png(
    p_heat,
    "Fig_SNP_Distance_Heatmap",
    width = cfg$fig_size$heatmap[["width"]],
    height = cfg$fig_size$heatmap[["height"]],
    out_dir = cfg$out_dir,
    dpi = cfg$png_dpi
  )

  p_recomb <- plot_recombination_impact(
    dist_core = aln$dist_core,
    dist_gubbins = aln$dist_gubbins,
    n_poly_input = aln$n_poly_input,
    n_poly_output = aln$n_poly_output,
    sites_removed = aln$sites_removed,
    pct_removed = aln$pct_removed
  )
  fig_recomb <- save_figure_png(
    p_recomb,
    "Fig_Recombination_Impact",
    width = cfg$fig_size$recomb[["width"]],
    height = cfg$fig_size$recomb[["height"]],
    out_dir = cfg$out_dir,
    dpi = cfg$png_dpi
  )

  p_recomb_panel <- plot_recombination_impact(
    dist_core = aln$dist_core,
    dist_gubbins = aln$dist_gubbins,
    n_poly_input = aln$n_poly_input,
    n_poly_output = aln$n_poly_output,
    sites_removed = aln$sites_removed,
    pct_removed = aln$pct_removed,
    internal_labels = FALSE,
    include_main_title = FALSE
  )

  p_panel <- plot_unified_panel(p_rect, p_circ, p_heat, p_recomb_panel)
  fig_panel <- save_figure_png(
    p_panel,
    "Fig_Composite_Panel",
    width = cfg$fig_size$panel[["width"]],
    height = cfg$fig_size$panel[["height"]],
    out_dir = cfg$out_dir,
    dpi = cfg$png_dpi
  )

  # Console summary
  message("\n============================================================")
  message("ANALYSIS SUMMARY")
  message("============================================================")
  message("Strains (tips): ", ape::Ntip(tree_display))
  message("Polymorphic sites (input; core alignment): ", scales::comma(aln$n_poly_input))
  message("Polymorphic sites (output; Gubbins): ", scales::comma(aln$n_poly_output))
  message("Recombination sites removed: ", scales::comma(aln$sites_removed), " (", aln$pct_removed, "%)")
  message("============================================================")
  message("Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

  invisible(list(
    fig_rect = fig_rect,
    fig_circ = fig_circ,
    fig_heat = fig_heat,
    fig_recomb = fig_recomb,
    fig_panel = fig_panel
  ))
}

if (sys.nframe() == 0) {
  main(config)
}