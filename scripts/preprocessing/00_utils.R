# ---
# jupyter:
#   jupytext:
#     formats: R:percent
#   kernelspec:
#     display_name: R
#     language: R
#     name: ir
# ---

# %% [markdown]
# # Shared Utility Functions

# %% [markdown]
# ## Library loading and configuration

# %% tags=["cell-1"]
library(tidyverse)
library(svglite)
library(Seurat)
library(patchwork)
library(viridis)
library(ggforce)
theme_set( theme_bw() )
root.dir <- here::here()
"%ni%" <- Negate("%in%")
opt.recalc <- FALSE
data.dir <- "../../R4.4.0/"
intermediate.prefix <- "20250604_"

# Chick-specific sex chromosome gene lists (used by 01_chick_preprocessing.R)
W.genes <- tryCatch(
    read_tsv(glue::glue("../../data/chick_W_genes.tsv"), col_names = "W", show_col_types = FALSE) %>% pull(W),
    error = function(e) { message("Note: chick_W_genes.tsv not found (ok for non-chick scripts)"); character(0) }
)
Z.genes <- tryCatch(
    read_tsv(glue::glue("../../data/chick_Z_genes.tsv"), col_names = "Z", show_col_types = FALSE) %>% pull(Z),
    error = function(e) { message("Note: chick_Z_genes.tsv not found (ok for non-chick scripts)"); character(0) }
)

# %% [markdown]
# ## export_seurat(): Seurat to h5ad/MEX export

# %% tags=["cell-3"]
export_seurat <- function(s, dir_path = "export", assay = "RNA", prefix = "", parallel = FALSE) {
    timestamp <- function(msg) {
        message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
    }
    timestamp("Starting export process...")
    dir.create(dir_path, showWarnings = FALSE)
    timestamp("Joining layer matrices...")
    joined_assay <- JoinLayers(s[[assay]])
    counts <- LayerData(joined_assay, layer = "counts")
    data   <- LayerData(joined_assay, layer = "data")
    # File paths
    raw_counts_temp <- file.path(dir_path, glue::glue("{prefix}raw_counts.mtx"))
    norm_data_temp  <- file.path(dir_path, glue::glue("{prefix}normalized_data.mtx"))
    # Matrix writing function
    write_matrix <- function(matrix, filepath, type) {
        tryCatch({
            Matrix::writeMM(matrix, file = filepath)
            timestamp(glue::glue("{type} matrix written successfully. File size: {file.size(filepath)} bytes"))
            TRUE
        }, error = function(e) {
            timestamp(glue::glue("Error writing {type} matrix: {e$message}"))
            FALSE
        })
    }
    # Compression function
    compress_file <- function(filepath, type) {
        tryCatch({
            system2("gzip", args = c("-f", filepath), stdout = TRUE, stderr = TRUE)
            if(file.exists(paste0(filepath, ".gz"))) {
                timestamp(glue::glue("{type} compression successful. Compressed size: {file.size(paste0(filepath, '.gz'))} bytes"))
                TRUE
            } else {
                timestamp(glue::glue("{type} compression failed - output file not found"))
                FALSE
            }
        }, error = function(e) {
            timestamp(glue::glue("Error compressing {type}: {e$message}"))
            FALSE
        })
    }
    # Handle matrix writing and compression based on parallel parameter
    if (parallel) {
        timestamp("Starting parallel matrix export...")
        # Write matrices in parallel
        p1 <- parallel::mcparallel(write_matrix(counts, raw_counts_temp, "Raw counts"))
        p2 <- parallel::mcparallel(write_matrix(data, norm_data_temp, "Normalized data"))
        # Wait for matrix writing to complete
        write_results <- parallel::mccollect()
        if (!all(unlist(write_results))) {
            stop("Error in parallel matrix writing")
        }
        # Parallel compression
        timestamp("Starting parallel compression...")
        p3 <- parallel::mcparallel(compress_file(raw_counts_temp, "Raw counts"))
        p4 <- parallel::mcparallel(compress_file(norm_data_temp, "Normalized data"))
        compress_results <- parallel::mccollect()
        if (!all(unlist(compress_results))) {
            stop("Error in parallel compression")
        }
    } else {
        timestamp("Starting sequential matrix export...")
        # Sequential matrix writing
        if (!write_matrix(counts, raw_counts_temp, "Raw counts")) {
            stop("Error writing raw counts matrix")
        }
        if (!write_matrix(data, norm_data_temp, "Normalized data")) {
            stop("Error writing normalized data matrix")
        }
        # Sequential compression
        timestamp("Starting sequential compression...")
        if (!compress_file(raw_counts_temp, "Raw counts")) {
            stop("Error compressing raw counts")
        }
        if (!compress_file(norm_data_temp, "Normalized data")) {
            stop("Error compressing normalized data")
        }
    }
    # Export features
    features <- rownames(counts)
    timestamp(glue::glue("Writing {length(features)} features to {dir_path}/{prefix}features.tsv..."))
    write.table(
        features,
        file = file.path(dir_path, glue::glue("{prefix}features.tsv")),
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE
    )
    # Export cell barcodes
    cells <- colnames(counts)
    timestamp(glue::glue("Writing {length(cells)} cell barcodes to {dir_path}/{prefix}barcodes.tsv..."))
    write.table(
        cells,
        file = file.path(dir_path, glue::glue("{prefix}barcodes.tsv")),
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE
    )
    # Metadata export
    timestamp("Exporting metadata...")
    meta <- s@meta.data
    write.table(
        meta,
        file = file.path(dir_path, glue::glue("{prefix}metadata.tsv")),
        sep = "\t",
        quote = FALSE,
        row.names = TRUE
    )
    timestamp(glue::glue("Metadata with {ncol(meta)} columns exported"))
    # Dimensional reductions export
    timestamp("Processing dimensional reductions...")
    red_names <- Reductions(s)
    for (red in red_names) {
        timestamp(glue::glue("Exporting {red} reduction..."))
        embedding <- Embeddings(s[[red]])
        # Add prefix to column names for easier parsing in Python
        colnames(embedding) <- paste0(red, "_", colnames(embedding))
        write.table(
            embedding,
            file = file.path(dir_path, glue::glue("{prefix}{red}_embeddings.tsv")),
            sep = "\t",
            quote = FALSE,
            row.names = TRUE
        )
        timestamp(glue::glue("{red} reduction with {ncol(embedding)} dimensions exported"))
    }
    timestamp("Export complete!")
}

# %% [markdown]
# ## plot_axial_expression(): axial expression plots with metadata support

# %% tags=["cell-5"]
#' Plot axial gene expression along DV or NT axis with metadata support
#'
#' Extended version supporting:
#' - Gene expression (original functionality)
#' - Continuous metadata (e.g., nCount_RNA, percent.mt)
#' - Categorical metadata as fractions (e.g., Phase_G1, Phase_S, Phase_G2M)
#' - Stacked bar histogram for categorical metadata
#'
#' @param seurat_obj Seurat object
#' @param genes Character vector of gene names (default NULL)
#' @param metadata_columns Character vector of continuous metadata column names (default NULL)
#' @param metadata_categorical Character vector of categorical metadata column names (default NULL)
#' @param metadata_colors Named vector of colors for metadata (optional)
#' @param metadata_linetype Line type for metadata ("dashed", "dotted", "solid")
#' @param add_scale_note Add caption about different scales (default TRUE)
#' @param axis_score Axis to plot along ("DV.Score" or "NT.Score")
#' @param normalize Normalize all curves to 0-1 (default FALSE)
#' @param ... All other original parameters (conditional_lines, add_histogram, etc.)
#'
#' @return patchwork object with main plot and optional histogram
#'
#' @examples
#' # Backward compatible - genes only
#' plot_axial_expression(seurat_obj, genes = c("FGF8", "PAX6"))
#'
#' # Genes + continuous metadata
#' plot_axial_expression(seurat_obj, genes = c("FGF8"),
#'                       metadata_columns = c("nCount_RNA"))
#'
#' # Genes + categorical metadata (with stacked histogram!)
#' plot_axial_expression(seurat_obj, genes = c("FGF8"),
#'                       metadata_categorical = c("Phase"))
#'
plot_axial_expression <- function(
  seurat_obj,
  genes = NULL,
  axis_score = "DV.Score",
  # NEW PARAMETERS for metadata
  metadata_columns = NULL,
  metadata_categorical = NULL,
  metadata_colors = NULL,
  metadata_linetype = "dashed",
  add_scale_note = TRUE,
  # ORIGINAL PARAMETERS
  normalize = FALSE,
  conditional_lines = FALSE,
  density_threshold = 50,
  n_bins = 50,
  window_size = 3,
  add_histogram = TRUE,
  histogram_bins = 20,
  smooth_method = "loess",
  smooth_span = 0.5,
  gam_k = 20,
  add_labels = FALSE,
  label_fontsize = 10,
  label_at = "max",
  label_method = "force",
  label_nudge_y = 0.1,
  label_force = 5,
  line_alpha = 1,
  reference_genes = NULL,
  reference_color = "grey70",
  reference_alpha = 0.05,
  text_size = 20,
  plot_height = 5,
  histogram_height = 1,
  slot = "data"
) {
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  require(patchwork)
  # ========================================================================
  # STEP 1: Data Fetching and Preparation
  # ========================================================================
  # Initialize variables to fetch
  vars_to_fetch <- c(axis_score)
  # Combine main genes and reference genes if provided
  if (!is.null(genes) || !is.null(reference_genes)) {
    all_genes <- unique(c(genes, reference_genes))
    vars_to_fetch <- c(vars_to_fetch, all_genes)
  } else {
    all_genes <- NULL
  }
  # Add continuous metadata columns
  if (!is.null(metadata_columns)) {
    vars_to_fetch <- c(vars_to_fetch, metadata_columns)
  }
  # Add categorical metadata columns (will expand to fractions)
  categorical_fraction_cols <- NULL
  if (!is.null(metadata_categorical)) {
    vars_to_fetch <- c(vars_to_fetch, metadata_categorical)
  }
  # Fetch data from Seurat object
  raw_data <- FetchData(seurat_obj, vars = unique(vars_to_fetch), slot = slot)
  # ========================================================================
  # STEP 2: Expand Categorical Metadata to Fractions
  # ========================================================================
  categorical_mapping <- list()  # Store category -> column mappings for histogram
  if (!is.null(metadata_categorical)) {
    for (meta_col in metadata_categorical) {
      # Get unique categories
      categories <- sort(unique(raw_data[[meta_col]]))
      categories <- categories[!is.na(categories)]
      category_cols <- c()
      # Create binary indicator for each category
      for (cat in categories) {
        frac_col <- paste0(meta_col, "_", cat)
        raw_data[[frac_col]] <- as.numeric(raw_data[[meta_col]] == cat)
        category_cols <- c(category_cols, frac_col)
      }
      categorical_fraction_cols <- c(categorical_fraction_cols, category_cols)
      categorical_mapping[[meta_col]] <- category_cols
    }
  }
  # ========================================================================
  # STEP 3: Pivot to Long Format and Label Data Types
  # ========================================================================
  # Columns to plot (exclude original categorical and axis)
  cols_to_plot <- setdiff(names(raw_data), c(axis_score, metadata_categorical))
  plot_data <- raw_data %>%
    select(all_of(c(axis_score, cols_to_plot))) %>%
    pivot_longer(cols = !any_of(axis_score), names_to = "variable", values_to = "expression") %>%
    mutate(
      data_type = case_when(
        variable %in% reference_genes ~ "reference_gene",
        variable %in% genes ~ "gene",
        variable %in% metadata_columns ~ "metadata_continuous",
        variable %in% categorical_fraction_cols ~ "metadata_categorical_fraction",
        TRUE ~ "unknown"
      ),
      is_reference = variable %in% reference_genes,
      is_categorical_fraction = data_type == "metadata_categorical_fraction"
    )
  # ========================================================================
  # STEP 4: Apply Smoothing
  # ========================================================================
  # Create base plot with smoothing
  if (smooth_method == "loess") {
    p_base <- plot_data %>%
      ggplot(aes_string(x = axis_score, y = "expression", colour = "variable")) +
      geom_smooth(linewidth = 1, method = "loess", span = smooth_span, se = FALSE)
  } else {
    p_base <- plot_data %>%
      ggplot(aes_string(x = axis_score, y = "expression", colour = "variable")) +
      geom_smooth(linewidth = 1, method = "gam", formula = y ~ s(x, k = gam_k), se = FALSE)
  }
  # Build plot to extract smooth data
  built_plot <- ggplot_build(p_base)
  smooth_data <- built_plot$data[[1]]
  # Get variable mapping and colors
  variable_mapping <- built_plot$plot$scales$get_scales("colour")$get_limits()
  if (is.null(variable_mapping)) {
    variable_mapping <- sort(unique(built_plot$plot$data$variable))
  }
  # Generate colors
  n_vars <- length(variable_mapping)
  default_colors <- scales::hue_pal()(n_vars)
  names(default_colors) <- variable_mapping
  # Apply custom colors for metadata if provided
  if (!is.null(metadata_colors)) {
    for (var_name in names(metadata_colors)) {
      if (var_name %in% names(default_colors)) {
        default_colors[var_name] <- metadata_colors[var_name]
      }
    }
  }
  # Add variable info to smooth data
  smooth_data <- smooth_data %>%
    mutate(
      variable = variable_mapping[group],
      data_type = case_when(
        variable %in% reference_genes ~ "reference_gene",
        variable %in% genes ~ "gene",
        variable %in% metadata_columns ~ "metadata_continuous",
        variable %in% categorical_fraction_cols ~ "metadata_categorical_fraction",
        TRUE ~ "unknown"
      ),
      is_reference = variable %in% reference_genes,
      is_categorical_fraction = data_type == "metadata_categorical_fraction",
      base_color = colour,
      plot_color = ifelse(is_reference, reference_color, colour)
    )
  # ========================================================================
  # STEP 5: Normalization
  # ========================================================================
  if (normalize) {
    smooth_data <- smooth_data %>%
      group_by(variable) %>%
      mutate(
        y_norm = if (is_categorical_fraction[1]) {
          # Categorical fractions: already 0-1, just ensure max=1
          pmin(y, 1.0)
        } else {
          # Genes and continuous metadata: normalize to 0-1
          y_range <- max(y, na.rm = TRUE) - min(y, na.rm = TRUE)
          if (y_range > 0) {
            (y - min(y, na.rm = TRUE)) / y_range
          } else {
            y
          }
        }
      ) %>%
      ungroup()
    y_var <- "y_norm"
    y_limits <- c(0, 1)
    # BACKWARD COMPATIBILITY: Use exact original label when no metadata
    has_metadata <- !is.null(metadata_columns) || !is.null(metadata_categorical)
    if (has_metadata) {
      y_label <- "Normalized values (0-1)"
      # Add scale note if mixing data types
      if (add_scale_note && !is.null(metadata_categorical)) {
        y_label <- paste0(y_label, "\n(categorical fractions: max=1, others: data-normalized)")
      }
    } else {
      # ORIGINAL: exact label from original function
      y_label <- "Normalized expression (0-1)"
    }
  } else {
    # BACKWARD COMPATIBILITY: Only create y_norm if we have metadata that needs special handling
    has_metadata <- !is.null(metadata_columns) || !is.null(metadata_categorical)
    if (has_metadata) {
      # Need to clip categorical fractions at 1.0
      smooth_data <- smooth_data %>%
        mutate(
          y_norm = if_else(is_categorical_fraction,
                          pmin(y, 1.0),  # Clip categorical at 1.0
                          y)
        )
      y_var <- "y_norm"
    } else {
      # ORIGINAL BEHAVIOR: Use raw 'y' column
      y_var <- "y"
    }
    y_label <- "Average log-normalized expression"
    y_limits <- NULL
    # Add scale note if mixing data types
    if (add_scale_note && has_metadata) {
      notes <- c()
      if (!is.null(metadata_columns)) {
        notes <- c(notes, "continuous metadata may have different scales")
      }
      if (!is.null(metadata_categorical)) {
        notes <- c(notes, "categorical shown as fractions (0-1)")
      }
      if (length(notes) > 0) {
        y_label <- paste0(y_label, "\n(", paste(notes, collapse = "; "), ")")
      }
    }
  }
  # ========================================================================
  # STEP 6: Handle Conditional Lines (if requested)
  # ========================================================================
  if (conditional_lines) {
    # Calculate density
    density_data <- plot_data %>%
      group_by(variable) %>%
      mutate(score_bin = cut(.data[[axis_score]], breaks = n_bins, include.lowest = TRUE)) %>%
      group_by(variable, score_bin) %>%
      summarise(
        n_points = n(),
        bin_center = mean(range(.data[[axis_score]], na.rm = TRUE), na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      group_by(variable) %>%
      arrange(bin_center) %>%
      mutate(bin_id = row_number())
    # Determine line type based on density
    density_with_linetype <- density_data %>%
      group_by(variable) %>%
      mutate(
        rolling_mean = zoo::rollapply(
          n_points, width = 2*window_size + 1,
          FUN = mean, partial = TRUE, fill = NA, align = "center"
        ),
        is_dense = rolling_mean >= density_threshold
      ) %>%
      ungroup()
    # Join with smooth data
    smooth_data <- smooth_data %>%
      group_by(variable) %>%
      mutate(
        x_bin = cut(x, breaks = n_bins, include.lowest = TRUE),
        bin_id = as.numeric(x_bin)
      ) %>%
      left_join(
        density_with_linetype %>% select(variable, bin_id, is_dense),
        by = c("variable", "bin_id")
      ) %>%
      mutate(
        line_type_density = ifelse(is_dense %in% TRUE, "solid", "dashed")
      ) %>%
      ungroup()
    # Override line type for metadata
    smooth_data <- smooth_data %>%
      mutate(
        final_line_type = case_when(
          data_type %in% c("metadata_continuous", "metadata_categorical_fraction") ~ metadata_linetype,
          TRUE ~ line_type_density
        ),
        # BACKWARD COMPATIBILITY: Keep 'line_type' field for label filtering
        line_type = line_type_density
      )
    # Create segments
    smooth_data_segmented <- smooth_data %>%
      arrange(variable, x) %>%
      group_by(variable) %>%
      mutate(
        linetype_change = final_line_type != lag(final_line_type, default = first(final_line_type)),
        segment_id = cumsum(linetype_change),
        unique_group = paste0(variable, "_", segment_id)
      ) %>%
      ungroup()
    # Create plot with conditional lines
    p_main <- ggplot()
    for (lt in unique(smooth_data_segmented$final_line_type)) {
      lt_data <- smooth_data_segmented %>% filter(final_line_type == lt)
      if (nrow(lt_data) > 0) {
        p_main <- p_main + geom_line(
          data = lt_data,
          aes(x = x, y = .data[[y_var]], colour = plot_color, group = unique_group),
          linewidth = 1, linetype = lt, alpha = line_alpha
        )
      }
    }
    p_main <- p_main + scale_colour_identity()
  } else {
    # Simple plot with different line types by data type
    p_main <- ggplot(smooth_data) +
      geom_line(
        data = smooth_data %>% filter(data_type == "gene"),
        aes(x = x, y = .data[[y_var]], colour = plot_color, group = variable),
        linewidth = 1, linetype = "solid", alpha = line_alpha
      )
    if (!is.null(reference_genes)) {
      p_main <- p_main + geom_line(
        data = smooth_data %>% filter(data_type == "reference_gene"),
        aes(x = x, y = .data[[y_var]], colour = plot_color, group = variable),
        linewidth = 1, linetype = "solid", alpha = reference_alpha
      )
    }
    if (!is.null(metadata_columns)) {
      p_main <- p_main + geom_line(
        data = smooth_data %>% filter(data_type == "metadata_continuous"),
        aes(x = x, y = .data[[y_var]], colour = plot_color, group = variable),
        linewidth = 1, linetype = metadata_linetype, alpha = line_alpha
      )
    }
    if (!is.null(metadata_categorical)) {
      p_main <- p_main + geom_line(
        data = smooth_data %>% filter(data_type == "metadata_categorical_fraction"),
        aes(x = x, y = .data[[y_var]], colour = plot_color, group = variable),
        linewidth = 1, linetype = "dotted", alpha = line_alpha
      )
    }
    p_main <- p_main + scale_colour_identity()
  }
  # ========================================================================
  # STEP 7: Label Handling (existing logic, applied to all variables)
  # ========================================================================
  if (add_labels) {
    # Prepare label data based on label_at parameter
    if (label_at == "max") {
      if (conditional_lines && exists("smooth_data_segmented")) {
        # ORIGINAL BEHAVIOR: Filter for solid lines only when conditional_lines is TRUE
        label_data <- smooth_data_segmented %>%
          filter(line_type == "solid") %>%
          group_by(variable, is_reference, plot_color) %>%
          filter(!is.na(.data[[y_var]])) %>%
          filter(.data[[y_var]] == max(.data[[y_var]], na.rm = TRUE)) %>%
          slice(1) %>%
          select(variable, x, all_of(y_var), is_reference, plot_color) %>%
          dplyr::rename(label_x = x, label_y = all_of(y_var))
      } else {
        label_data <- smooth_data %>%
          group_by(variable, is_reference, plot_color) %>%
          filter(!is.na(.data[[y_var]])) %>%
          filter(.data[[y_var]] == max(.data[[y_var]], na.rm = TRUE)) %>%
          slice(1) %>%
          select(variable, x, all_of(y_var), is_reference, plot_color) %>%
          dplyr::rename(label_x = x, label_y = all_of(y_var))
      }
    } else if (label_at == "endpoints") {
      # ORIGINAL BEHAVIOR: Use smooth_data and summarise to get endpoint with higher value
      label_data <- smooth_data %>%
        group_by(variable, is_reference, plot_color) %>%
        filter(!is.na(.data[[y_var]])) %>%
        summarise(
          left_x = min(x),
          left_y = .data[[y_var]][which.min(x)],
          right_x = max(x),
          right_y = .data[[y_var]][which.max(x)],
          label_x = ifelse(left_y > right_y, left_x, right_x),
          label_y = ifelse(left_y > right_y, left_y, right_y),
          .groups = "drop"
        )
    }
    # Apply labeling method
    if (label_method == "circle_repel") {
      # Combine circles with repelled labels
      require(ggforce)
      require(ggrepel)
      # Add circles (with line colors)
      for(i in 1:nrow(label_data)) {
        row <- label_data[i,]
        label_alpha <- ifelse(row$is_reference, reference_alpha, 0.1)
        # Add circle with line color
        p_main <- p_main +
          ggforce::geom_mark_circle(
            data = row,
            aes(x = label_x, y = label_y),
            fill = row$plot_color,
            expand = unit(2, "mm"),
            label.fontsize = 0,  # Hide the built-in label
            con.colour = NA,      # No connector line
            alpha = label_alpha
          )
      }
      # Determine nudge direction based on y position
      # If labels are near the top (>0.8 normalized), nudge down; otherwise nudge up
      label_data <- label_data %>%
        mutate(
          nudge_direction = ifelse(label_y > 0.8, -abs(label_nudge_y), abs(label_nudge_y))
        )
      # Add repelled labels on top
      p_main <- p_main +
        geom_label_repel(
          data = label_data,
          aes(x = label_x, y = label_y, label = variable),
          color = label_data$plot_color,  # Use line colors for text
          alpha = ifelse(label_data$is_reference,
                        pmin(reference_alpha * 5, 0.8),  # Cap reference alpha at 0.8
                        1),  # Full opacity for main genes
          fill = alpha("white", 0.8),
          size = label_fontsize/3,
          nudge_y = label_data$nudge_direction,  # Use calculated nudge direction
          force = label_force,
          box.padding = 0.5,
          point.padding = 0.3,
          segment.color = "grey50",
          segment.alpha = 0.4,
          label.padding = unit(0.15, "lines"),
          max.overlaps = Inf,
          direction = "both",  # Allow movement in both x and y
          show.legend = FALSE
        )
    } else {
      # DEFAULT: Original force method (backward compatible)
      require(ggforce)
      for(i in 1:nrow(label_data)) {
        row <- label_data[i,]
        label_alpha <- ifelse(row$is_reference, reference_alpha, 0.1)
        label_color <- ifelse(row$is_reference, reference_color, "gray60")
        p_main <- p_main +
          ggforce::geom_mark_circle(
            data = row,
            aes(x = label_x, y = label_y, label = variable),
            fill = row$plot_color,
            expand = unit(2, "mm"),
            label.fontsize = label_fontsize,
            label.buffer = unit(3, "mm"),
            con.colour = label_color,
            con.type = "straight",
            label.colour = ifelse(row$is_reference, reference_color, "black"),
            con.cap = unit(0, "mm"),
            alpha = label_alpha
          )
      }
    }
  }
  # ========================================================================
  # STEP 8: Finalize Main Plot
  # ========================================================================
  # Y-axis expansion matching original function
  y_expansion <- if (add_labels && label_method == "circle_repel") {
    expansion(add = 0, mult = c(0.02, 0.1))  # More space at top for labels
  } else {
    expansion(add = 0, mult = c(0, 0.05))  # Original expansion: no bottom, 5% top
  }
  p_main <- p_main +
    expand_limits(y = 0) +
    scale_y_continuous(limits = y_limits, expand = y_expansion) +
    labs(x = axis_score, y = y_label) +
    theme(
      text = element_text(size = text_size),
      legend.position = "none"
    )
  # ========================================================================
  # STEP 9: Histogram (with Stacked Bars for Categorical Metadata!)
  # ========================================================================
  if (add_histogram) {
    # Check if we have categorical metadata
    has_categorical <- !is.null(metadata_categorical) && length(categorical_mapping) > 0
    if (has_categorical) {
      # ===== STACKED BAR CHART FOR CATEGORICAL METADATA =====
      # For now, support one categorical column at a time for stacking
      # (multiple would be too complex)
      primary_categorical <- names(categorical_mapping)[1]
      category_columns <- categorical_mapping[[primary_categorical]]
      # Prepare histogram data with category proportions
      hist_data <- raw_data %>%
        mutate(bin = cut(.data[[axis_score]], breaks = histogram_bins, include.lowest = TRUE)) %>%
        group_by(bin) %>%
        summarise(
          across(all_of(category_columns), mean, na.rm = TRUE),
          bin_center = mean(.data[[axis_score]], na.rm = TRUE),
          n = n(),
          .groups = 'drop'
        ) %>%
        pivot_longer(cols = all_of(category_columns), names_to = "category", values_to = "proportion")
      # Get colors for categories (match main plot colors)
      category_colors <- default_colors[category_columns]
      names(category_colors) <- category_columns
      # Create stacked bar chart (flipped downward with negative values)
      p_hist <- ggplot(hist_data, aes(x = bin_center, y = -proportion, fill = category)) +
        geom_col(position = "stack", width = diff(range(raw_data[[axis_score]])) / histogram_bins) +
        scale_fill_manual(values = category_colors,
                         labels = gsub(paste0(primary_categorical, "_"), "", names(category_colors))) +
        scale_y_continuous(limits = c(-1, 0), expand = c(0, 0),
                          labels = function(x) abs(x),  # Show positive labels
                          name = "Fraction") +
        labs(x = axis_score, fill = primary_categorical) +
        theme_minimal(base_size = text_size * 0.75) +
        theme(
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          legend.position = "right",
          legend.key.size = unit(0.4, "cm"),
          legend.text = element_text(size = text_size * 0.6),
          plot.margin = margin(0, 5, 5, 5)
        )
    } else {
      # ===== REGULAR HISTOGRAM (original behavior - flipped downward) =====
      p_hist <- raw_data %>%
        ggplot(aes(x = .data[[axis_score]])) +
        geom_histogram(aes(y = -after_stat(count)), bins = histogram_bins,
                      fill = "gray60", alpha = 0.7) +
        scale_y_continuous(expand = c(0, 0),
                          labels = function(x) abs(x),
                          name = "# of cells") +
        theme_minimal(base_size = text_size * 0.75) +
        theme(
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          plot.margin = margin(0, 5, 5, 5)
        )
    }
    # Combine plots: MAIN plot on TOP, histogram on BOTTOM
    combined_plot <- p_main / p_hist +
      plot_layout(heights = c(plot_height, histogram_height))
    return(combined_plot)
  } else {
    return(p_main)
  }
}
