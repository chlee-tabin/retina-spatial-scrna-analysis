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
# # Figure 6H + Figure S23: Area-specific Differential Gene Expression

# %%
source("../preprocessing/00_utils.R")

# %% [markdown]
# ## Output directory setup

# %%
FIGURES_BASE <- file.path(dirname(sys.frame(1)$ofile %||% "."), "..", "..", "figures")
dir.create(file.path(FIGURES_BASE, "Figure6"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(FIGURES_BASE, "Figure_SF23"), recursive = TRUE, showWarnings = FALSE)
# The manuscript Supplementary Table 1 (chick, Fig 6H) / Table 2 (human, Fig
# S23) are emitted by the dedicated companion script
# scripts/analysis/area_significant_deg.R (writes to
# docs/supplementary_tables/SuppTable{1,2}_*.csv). This figure script no
# longer writes its own tables; the volcanos below drive from `deg_results`
# in memory.

# %% [markdown]
# ## Utility: plot.retina3() region selection function
#
# This helper is shared between the chick (F6H) and human (SF23) blocks below
# -- one definition, two consumers.

# %% tags=["cell-236"]
plot.retina3 <- function(obj, gene, bin_size = 50,
                         percentile = 1.00, highlight_quadrants = NULL,
                         n_grid = 10, highlight_color = "blue", highlight_alpha = 0.3) {
  # Extract cell metadata including barcodes
  cell_data <- FetchData(
    obj,
    vars = c("DV.Score", "NT.Score", gene)
  )
  # Add cell barcodes as a column
  cell_data$barcode <- rownames(cell_data)
  # Get data range for creating partition lines
  dv_range <- range(cell_data$DV.Score)
  nt_range <- range(cell_data$NT.Score)
  # Calculate grid positions (0 to n_grid-1)
  dv_partitions <- dv_range[1] + (dv_range[2] - dv_range[1]) * (0:n_grid)/n_grid
  nt_partitions <- nt_range[1] + (nt_range[2] - nt_range[1]) * (0:n_grid)/n_grid
  # Initialize vector to store selected cell barcodes
  selected_cells <- c()
  # Extract cells from highlighted quadrants if specified
  if (!is.null(highlight_quadrants)) {
    for (quad in highlight_quadrants) {
      nt_idx <- quad[1]  # x coordinate
      dv_idx <- quad[2]  # y coordinate
      # Define quadrant boundaries
      nt_min <- nt_partitions[nt_idx + 1]
      nt_max <- nt_partitions[nt_idx + 2]
      dv_min <- dv_partitions[dv_idx + 1]
      dv_max <- dv_partitions[dv_idx + 2]
      # Find cells in this quadrant
      quadrant_cells <- cell_data$barcode[
        cell_data$NT.Score >= nt_min &
        cell_data$NT.Score < nt_max &
        cell_data$DV.Score >= dv_min &
        cell_data$DV.Score < dv_max
      ]
      selected_cells <- c(selected_cells, quadrant_cells)
    }
  }
  # Create grid label positions (centers of grid cells)
  label_positions <- expand.grid(
    x = nt_partitions[-length(nt_partitions)] + diff(nt_partitions)/2,
    y = dv_partitions[-length(dv_partitions)] + diff(dv_partitions)/2
  )
  # Add grid indices to label positions
  label_positions$label <- sprintf("(%d,%d)",
                                 rep(0:(n_grid-1), times = n_grid), # NT index (x)
                                 rep(0:(n_grid-1), each = n_grid))  # DV index (y)
  # Create highlight rectangles
  highlight_rects <- NULL
  if (!is.null(highlight_quadrants)) {
    highlight_rects <- do.call(rbind, lapply(highlight_quadrants, function(quad) {
      nt_idx <- quad[1]  # x coordinate
      dv_idx <- quad[2]  # y coordinate
      data.frame(
        xmin = nt_partitions[nt_idx + 1],
        xmax = nt_partitions[nt_idx + 2],
        ymin = dv_partitions[dv_idx + 1],
        ymax = dv_partitions[dv_idx + 2]
      )
    }))
  }
  # Create the plot
  p <- ggplot(cell_data, aes(x = NT.Score, y = DV.Score, z = !!sym(gene))) +
    stat_summary_2d(fun = mean, bins = bin_size)
  # Add highlight rectangles if specified
  if (!is.null(highlight_rects)) {
    p <- p + geom_rect(data = highlight_rects,
                      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                      fill = highlight_color, alpha = highlight_alpha,
                      inherit.aes = FALSE)
  }
  # Add grid lines and labels
  p <- p +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "salmon") +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "salmon") +
    geom_vline(xintercept = nt_partitions, linetype = "dotted",
               colour = "navy", linewidth = 0.5, alpha = 0.3) +
    geom_hline(yintercept = dv_partitions, linetype = "dotted",
               colour = "navy", linewidth = 0.5, alpha = 0.3) +
    geom_text(data = label_positions,
              aes(x = x, y = y, label = label),
              inherit.aes = FALSE,
              size = 3, alpha = 0.7) +
    theme_bw() +
    theme(text = element_text(size = 20)) +
    labs(title = gene)
  # Extract and process values for color scale
  pb <- ggplot_build(p)
  values <- pb$data[[1]]$value
  pval <- quantile(values, probs = percentile)
  # Adjust percentile if needed
  if (pval == min(values)) {
    while (pval == min(values) && percentile < 1.0) {
      percentile <- percentile + 0.01
      pval <- quantile(values, probs = percentile)
    }
    adjusted <- "*"
  } else if (pval == max(values)) {
    while (pval == max(values) && percentile > 0.0) {
      percentile <- percentile - 0.01
      pval <- quantile(values, probs = percentile)
    }
    adjusted <- "*"
  } else {
    adjusted <- ""
  }
  # Finalize and print the plot
  final_plot <- p +
    scale_fill_viridis(
      option = "viridis",
      limits = c(min(values), pval),
      oob = scales::squish
    ) +
    labs(
      fill = "expression",
      caption = glue::glue("Percentile cut-off: {percentile}{adjusted}")
    )
  print(final_plot)
  # Return the selected cell barcodes
  return(unique(selected_cells))
}

# %% [markdown]
# ## F6H: Chick area selection (HAA, Temporal, Nasal, Dorsal, Ventral, DVcentral, NTcentral)

# %% [markdown]
# ### Select area: HAA

# %% tags=["cell-238"]
options( repr.plot.width = 10.5, repr.plot.height = 10.5 )
# `haa.area` (not `fovea.area`) -- chick uses the HAA (high-acuity area)
# label; the variable name should match the manuscript / SuppTable1 naming.
haa.area <-
plot.retina3(
    fabp7,
    "CYP26C1",
    highlight_quadrants = list(
        c(7,5),c(6,5),
        c(7,4),c(6,4)
    ),
    highlight_color = "blue",
    highlight_alpha = 0.5,
    percentile = 0.95,
    n_grid = 10
)
# Gate-drift forcing function: the chick HAA gate must produce exactly 5,971
# cells (= the published Fig 6H HAA arm; 20.6% of 29,025 chick RPCs). The
# companion analysis script scripts/analysis/area_significant_deg.R encodes
# the same gate in a range-tuple parameterization with a matching stopifnot,
# so if anyone edits the quadrant list above without updating the
# analysis-script counterpart, this assertion trips before silently-wrong
# Fig 6H or Supp Table 1 outputs. Note: no `is.null` guard -- if plot.retina3
# ever returns NULL on a degenerate gate, `length(NULL) == 5971L` is FALSE
# and stopifnot trips. That is exactly the failure mode we want to catch.
stopifnot(length(haa.area) == 5971L)
ggsave( file.path(FIGURES_BASE, "Figure6", "F6H_select_HAA.png"), width = 10.5, height = 10.5 )

# %% [markdown]
# ### Select area: Temporal

# %% tags=["cell-240"]
options( repr.plot.width = 10.5, repr.plot.height = 10.5 )
temporal.area <-
plot.retina3(
    fabp7, 
    "FOXD1", 
    highlight_quadrants = list(
        c(0,0),c(1,0),c(2,0),
        c(0,1),c(1,1),c(2,1),
        c(0,2),c(1,2),c(2,2),
        c(0,3),c(1,3),c(2,3),
        c(0,4),c(1,4),c(2,4),
        c(0,5),c(1,5),c(2,5),
        c(0,6),c(1,6),c(2,6)
    ),
    highlight_color = "blue",
    highlight_alpha = 0.5,
    percentile = 0.95,
    n_grid = 7
)
ggsave( file.path(FIGURES_BASE, "Figure6", "F6H_select_Temporal.png"), width = 10.5, height = 10.5 )

# %% [markdown]
# ### Select area: Nasal

# %% tags=["cell-242"]
options( repr.plot.width = 10.5, repr.plot.height = 10.5 )
nasal.area <-
plot.retina3(
    fabp7, 
    "FOXG1", 
    highlight_quadrants = list(
        c(5,0),c(6,0),
        c(5,1),c(6,1),
        c(5,2),c(6,2),
        c(5,3),c(6,3),
        c(5,4),c(6,4),
        c(5,5),c(6,5),
        c(5,6),c(6,6)
    ),
    highlight_color = "blue",
    highlight_alpha = 0.5,
    percentile = 0.95,
    n_grid = 7
)
ggsave( file.path(FIGURES_BASE, "Figure6", "F6H_select_Nasal.png"), width = 10.5, height = 10.5 )

# %% [markdown]
# ### Select area: Dorsal

# %% tags=["cell-244"]
options( repr.plot.width = 10.5, repr.plot.height = 10.5 )
dorsal.area <-
plot.retina3(
    fabp7, 
    "ALDH1A1", 
    highlight_quadrants = list(
        c(0,6),c(1,6),c(2,6),c(3,6),c(4,6),c(5,6),c(6,6),
        c(0,5),c(1,5),c(2,5),c(3,5),c(4,5),c(5,5),c(6,5),
        c(0,4),c(1,4),c(2,4),c(3,4),c(4,4),c(5,4),c(6,4)
    ),
    highlight_color = "blue",
    highlight_alpha = 0.5,
    percentile = 0.95,
    n_grid = 7
)
ggsave( file.path(FIGURES_BASE, "Figure6", "F6H_select_Dorsal.png"), width = 10.5, height = 10.5 )

# %% [markdown]
# ### Select area: Ventral

# %% tags=["cell-246"]
options( repr.plot.width = 10.5, repr.plot.height = 10.5 )
ventral.area <-
plot.retina3(
    fabp7, 
    "VAX1", 
    highlight_quadrants = list(
        c(0,0),c(1,0),c(2,0),c(3,0),c(4,0),c(5,0),c(6,0),
        c(0,1),c(1,1),c(2,1),c(3,1),c(4,1),c(5,1),c(6,1),
        c(0,2),c(1,2),c(2,2),c(3,2),c(4,2),c(5,2),c(6,2)
    ),
    highlight_color = "blue",
    highlight_alpha = 0.5,
    percentile = 0.95,
    n_grid = 7
)
ggsave( file.path(FIGURES_BASE, "Figure6", "F6H_select_Ventral.png"), width = 10.5, height = 10.5 )

# %% [markdown]
# ### Select area: DVcentral

# %% tags=["cell-248"]
options( repr.plot.width = 10.5, repr.plot.height = 10.5 )
DVcentral.area <-
plot.retina3(
    fabp7, 
    "BMP2", 
    highlight_quadrants = list(
        c(0,3),c(1,3),c(2,3),c(3,3),c(4,3),c(5,3),c(6,3)
    ),
    highlight_color = "blue",
    highlight_alpha = 0.5,
    percentile = 0.95,
    n_grid = 7
)
ggsave( file.path(FIGURES_BASE, "Figure6", "F6H_select_DVcentral.png"), width = 10.5, height = 10.5 )

# %% [markdown]
# ### Select area: NTcentral

# %% tags=["cell-250"]
options( repr.plot.width = 10.5, repr.plot.height = 10.5 )
NTcentral.area <-
plot.retina3(
    fabp7, 
    "CYP1B1", 
    highlight_quadrants = list(
        c(5,0),c(5,1),c(5,2),c(5,3),c(5,4),c(5,5),c(5,6),c(5,7),c(5,8),
        c(6,0),c(6,1),c(6,2),c(6,3),c(6,4),c(6,5),c(6,6),c(6,7),c(6,8)
    ),
    highlight_color = "blue",
    highlight_alpha = 0.5,
    percentile = 0.95,
    n_grid = 9
)
ggsave( file.path(FIGURES_BASE, "Figure6", "F6H_select_NTcentral.png"), width = 10.5, height = 10.5 )

# %% [markdown]
# ## Pseudobulk DEG functions

# %% [markdown]
# ## Pseudobulk DEG

# %% tags=["cell-254"]
# Pseudobulk DEG runner -- iterates over the area.* columns added by
# add.area.columns() and fits a glmGamPoi `~ library + area` per region. Used
# by both the chick (F6H) and human (SF23) blocks below; donor_col differs
# (genotype vs. sample) and is passed explicitly at the call site.
run_all_pseudobulk_deg <- function(seurat_obj, donor_col, library_col = "library", min_cells = 50) {
    require(glmGamPoi)
    require(SingleCellExperiment)
    require(Seurat)
    # Get all area columns
    area_cols <- grep("^area\\.", colnames(seurat_obj@meta.data), value = TRUE)
    # Store results for each area
    all_results <- list()
    # Track per-region failures so we can surface them loudly at the end --
    # silently dropping a failed territory from a 7-panel volcano produces a
    # 6-panel figure that reads as "no significant DEGs" rather than "the DEG
    # pipeline crashed", which is the worst silent-failure mode here.
    region_failures <- list()
    # Process each area
    for(area_col in area_cols) {
        message("\nProcessing ", area_col)
        tryCatch({
            # Create minimal Seurat object with required metadata
            minimal_obj <- CreateSeuratObject(
                counts = LayerData(seurat_obj, layer = "counts"),
                meta.data = data.frame(
                    library = seurat_obj@meta.data[[library_col]],
                    area = seurat_obj@meta.data[[area_col]],
                    donor = seurat_obj@meta.data[[donor_col]],
                    row.names = colnames(seurat_obj)
                )
            )
            # Normalize data
            minimal_obj <- NormalizeData(minimal_obj, verbose = FALSE)
            # Calculate mean expression from the normalized data used in analysis
            normalized_data <- LayerData(minimal_obj, layer = "data")
            gene_mean_exp <- rowMeans(normalized_data)
            names(gene_mean_exp) <- rownames(normalized_data)  # Ensure proper names
            # Convert to SingleCellExperiment
            sce <- as.SingleCellExperiment(minimal_obj)
            # Create pseudobulk data
            pseudobulk_data <- glmGamPoi::pseudobulk(
                sce,
                group_by = vars(
                    library,
                    area,
                    donor
                )
            )
            # Gate out pseudobulk samples below min_cells. (Previously min_cells was
            # accepted but never applied; tiny library x area x donor groups -- down to
            # single cells -- inflate dispersion for moderate-expression genes. min_cells
            # = 50 is a standard pseudobulk floor; samples below it are dropped.)
            .grp_n <- minimal_obj@meta.data %>%
                dplyr::count(library, area, donor, name = "ncell") %>%
                dplyr::mutate(dplyr::across(c(library, area, donor), as.character))
            .cd <- as.data.frame(colData(pseudobulk_data))
            .cd$.ord <- seq_len(nrow(.cd))
            .nvec <- dplyr::left_join(
                dplyr::mutate(.cd, dplyr::across(c(library, area, donor), as.character)),
                .grp_n, by = c("library", "area", "donor")
            )
            .nvec <- .nvec$ncell[order(.nvec$.ord)]
            pseudobulk_data <- pseudobulk_data[, .nvec >= min_cells]
            # Create design matrix (drop levels emptied by the min_cells gate)
            colData(pseudobulk_data)$area <- droplevels(factor(colData(pseudobulk_data)$area))
            colData(pseudobulk_data)$library <- droplevels(factor(colData(pseudobulk_data)$library))
            # Fit model using formula notation
            fit <- glm_gp(
                pseudobulk_data,
                design = ~ library + area
            )
            # Test for differential expression
            de_results <- test_de(fit, contrast = "area1")
            # Fix gene names in DE results - they should match the pseudobulk
            # row names exactly. A mismatch would mean glmGamPoi dropped or
            # reordered rows internally; we abort rather than silently relabel
            # with rownames(fit), because the gene-name -> log2FC mapping is
            # the entire point of the table.
            actual_gene_names <- rownames(pseudobulk_data)
            if (nrow(de_results) != length(actual_gene_names)) {
                stop("DE row count (", nrow(de_results),
                     ") does not match pseudobulk gene count (",
                     length(actual_gene_names), ") for ", area_col)
            }
            rownames(de_results) <- actual_gene_names
            # Add area information, mean expression, and clean up results
            gene_names <- rownames(de_results)
            de_results <- de_results %>%
                as.data.frame() %>%
                dplyr::mutate(
                    area = sub("^area\\.", "", area_col),
                    gene = gene_names,
                    meanExp = gene_mean_exp[gene_names]
                ) %>%
                dplyr::select(gene, area, meanExp, everything())
            # Store results
            all_results[[area_col]] <- de_results
        }, error = function(e) {
            message("Error processing ", area_col, ": ", e$message)
            region_failures[[area_col]] <<- e$message
        })
    }
    # If anything failed, surface it as a warning so it shows up in
    # `warnings()` -- `message()` to stderr alone gets buried under the
    # glmGamPoi progress chatter. The full per-region error messages are
    # attached as an attribute on the returned data.frame for post-mortem
    # inspection (`attr(deg_results, "region_failures")`).
    if (length(region_failures) > 0L) {
        warning(immediate. = TRUE,
                sprintf("Pseudobulk DEG failed for %d of %d territories: %s. ",
                        length(region_failures), length(area_cols),
                        paste(names(region_failures), collapse = ", ")),
                "The volcano below will be missing these panels.")
    }
    # Remove any NULL results from errors
    all_results <- all_results[!sapply(all_results, is.null)]
    # Combine all results
    if (length(all_results) > 0) {
        final_results <- do.call(rbind, all_results)
        rownames(final_results) <- NULL
    } else {
        final_results <- NULL
    }
    # Attach the per-region failure messages (if any) so they survive the
    # function return and are inspectable from the REPL.
    if (length(region_failures) > 0L && !is.null(final_results))
        attr(final_results, "region_failures") <- region_failures
    return(final_results)
}

# %% tags=["cell-255"]
areas <-
list(
    "HAA"       = haa.area,
    "Temporal"  = temporal.area,
    "Nasal"     = nasal.area,
    "Dorsal"    = dorsal.area,
    "Ventral"   = ventral.area,
    "DVcentral" = DVcentral.area,
    "NTcentral" = NTcentral.area
)

# %% tags=["cell-256"]
# Function to add area columns to Seurat object
add.area.columns <- function(seurat_obj, areas) {
    # Create a data frame with cell barcodes as rownames
    area_matrix <- matrix(0, 
                         nrow = ncol(seurat_obj), 
                         ncol = length(areas),
                         dimnames = list(colnames(seurat_obj), 
                                       paste0("area.", names(areas))))
    # Fill in 1s for cells in each area
    for (area_name in names(areas)) {
        area_matrix[areas[[area_name]], paste0("area.", area_name)] <- 1
    }
    # Convert to data frame
    area_df <- as.data.frame(area_matrix)
    # Add columns to Seurat object
    seurat_obj <- AddMetaData(seurat_obj, area_df)
    return(seurat_obj)
}
fabp7 <- add.area.columns(fabp7, areas)

# %% [markdown]
# ## F6H: Run DEG analysis (chick)

# %% tags=["cell-257"]
# Run the analysis
deg_results <- run_all_pseudobulk_deg(fabp7, donor_col = "genotype", library_col = "library", min_cells = 50)
# deg_results

# %% [markdown]
# ### Supplementary Table 1 (chick): area DEG -- emitted by area_significant_deg.R
#
# The full statistically-significant area-DEG table (adj-p < 0.05, no 2-fold
# headline cut, all 7 territories, `min_cells = 50`) is emitted by the dedicated
# script `scripts/analysis/area_significant_deg.R`, which produces
# `docs/supplementary_tables/SuppTable1_chick_area_significant.csv` and the
# matching human SuppTable2 from the same gates. That table is the corrected
# manuscript Supp Table 1 (chick) / Table 2 (human). The volcano below uses
# `deg_results` in-memory.

# %% [markdown]
# ## F6H: Volcano plot (chick)

# %% tags=["cell-260"]
options( repr.plot.width = 10.5, repr.plot.height = 7 )
p <-
deg_results %>%
mutate(
    flag = abs(lfc) > 1 & adj_pval < 0.05,
#     lfc = ifelse( abs(lfc) > 5, sign(lfc)*5, lfc )
) %>% {
    ggplot( ., aes( x = lfc, y = -log10( adj_pval ) ) ) +
    geom_point( 
        aes( 
            colour = flag, 
            size = sqrt(meanExp)
        ),
        alpha = 0.5
    ) +
    geom_vline( xintercept = c(-1, 1), linetype = "dashed", colour = "salmon" ) +
    geom_hline( yintercept = -log10( 0.05 ), linetype = "dashed", colour = "salmon" ) +
    # annotate(
    #     "text",
    #     x = -10,
    #     y = 18,
    #     label = "Fovea",
    #     hjust = 0,
    #     size = 5
    # ) +
    ggforce::geom_mark_circle(
        data = . %>% 
                dplyr::arrange( desc( lfc ) ) %>%
                group_by( area ) %>%
                dplyr::filter(
                    adj_pval < 0.05,
                    abs(lfc) > 1,
                    # Drop LOC* (uncharacterized) chick gene-symbol labels
                    # from the volcano callout set.
                    !grepl( "^LOC", name )
                ) %>%
                dplyr::filter( 
#                    gene %in% c( "CYP26C1", "MYOF", "HHEX", "AKR1D1" ) ),
                    row_number() <= 5 # Top 5
                ),
        aes(
            label = gene,
            group = gene,
#            fill = after_scale(alpha(colour, 0.1))
        ),
        alpha = 1,
        expand = 0.02,
        con.cap = unit(0, "mm"),
        con.arrow = arrow(
            length = unit(1, "mm"),
            ends = "last",
            type = "closed"
        ),
        label.buffer = unit(1, 'mm'),
        label.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
        label.fontsize = 15,  # geom_mark_circle uses points
        label.fill = "white"
    ) +
    scale_y_continuous( expand = expansion( add = 0, mult = c(0, 0.05) ) ) +
    scale_colour_manual( 
        values = c( 
            "TRUE" = "salmon", 
            "FALSE" = "grey" 
        )
    ) +
    theme(
    #    legend.position = "none",
        text = element_text( size = 20 )
    ) +
    labs(
        x = "log2 fold change"
    ) + 
    guides(
        flag = "none",
        colour = "none",
        size = guide_legend(title = "expression\nlevel")
    ) +
    facet_wrap( area ~ ., ncol = 7 )
}
options( repr.plot.width = 49, repr.plot.height = 7 )
p
ggsave( p, file=file.path(FIGURES_BASE, "Figure6", "F6H_volcano.png"), width = 49, height = 7 )

# %% [markdown]
# ## SF23: Human area selection

# %% [markdown]
# ### Human RPC sample composition

# %% tags=["cell-293"]
human@meta.data %>%
dplyr::count( sample )

# (plot.retina3 defined once near the top of this script and reused here for
# the human SF23 selections.)

# %% [markdown]
# #### Select area: Fovea

# %% tags=["cell-297"]
options( repr.plot.width = 10.5, repr.plot.height = 10.5 )
fovea.area <-
plot.retina3(
    human,
    "CYP26C1",
    highlight_quadrants = list(
        c(3,5),c(4,5),
        c(3,4),c(4,4)
    ),
    highlight_color = "blue",
    highlight_alpha = 0.5,
    percentile = 0.95,
    n_grid = 10
)
# Gate-drift forcing function: the human fovea gate must produce exactly 2,236
# cells (= the marker-guided localization). Mirror of the chick HAA assertion
# in cell-238 and the matching stopifnot in
# scripts/analysis/area_significant_deg.R. No `is.null` guard, by design --
# `length(NULL) == 2236L` is FALSE and stopifnot trips, which is what we want.
stopifnot(length(fovea.area) == 2236L)
ggsave( file.path(FIGURES_BASE, "Figure_SF23", "SF23_select_Fovea.pdf"), width = 10.5, height = 10.5 )

# %% [markdown]
# #### Select area: Temporal

# %% tags=["cell-299"]
options( repr.plot.width = 10.5, repr.plot.height = 10.5 )
temporal.area <-
plot.retina3(
    human, 
    "FOXD1", 
    highlight_quadrants = list(
        c(0,0),c(1,0),c(2,0),
        c(0,1),c(1,1),c(2,1),
        c(0,2),c(1,2),c(2,2),
        c(0,3),c(1,3),c(2,3),
        c(0,4),c(1,4),c(2,4),
        c(0,5),c(1,5),c(2,5),
        c(0,6),c(1,6),c(2,6)
    ),
    highlight_color = "blue",
    highlight_alpha = 0.5,
    percentile = 0.95,
    n_grid = 7
)
ggsave( file.path(FIGURES_BASE, "Figure_SF23", "SF23_select_Temporal.pdf"), width = 10.5, height = 10.5 )

# %% [markdown]
# #### Select area: Nasal

# %% tags=["cell-301"]
options( repr.plot.width = 10.5, repr.plot.height = 10.5 )
nasal.area <-
plot.retina3(
    human, 
    "FOXG1", 
    highlight_quadrants = list(
        c(5,0),c(6,0),c(4,0),
        c(5,1),c(6,1),c(4,1),
        c(5,2),c(6,2),c(4,2),
        c(5,3),c(6,3),c(4,3),
        c(5,4),c(6,4),c(4,4),
        c(5,5),c(6,5),c(4,5),
        c(5,6),c(6,6),c(4,6)
    ),
    highlight_color = "blue",
    highlight_alpha = 0.5,
    percentile = 0.95,
    n_grid = 7
)
ggsave( file.path(FIGURES_BASE, "Figure_SF23", "SF23_select_Nasal.pdf"), width = 10.5, height = 10.5 )

# %% [markdown]
# #### Select area: Dorsal

# %% tags=["cell-303"]
options( repr.plot.width = 10.5, repr.plot.height = 10.5 )
dorsal.area <-
plot.retina3(
    human, 
    "ALDH1A1", 
    highlight_quadrants = list(
        c(0,6),c(1,6),c(2,6),c(3,6),c(4,6),c(5,6),c(6,6),
        c(0,5),c(1,5),c(2,5),c(3,5),c(4,5),c(5,5),c(6,5),
        c(0,4),c(1,4),c(2,4),c(3,4),c(4,4),c(5,4),c(6,4)
    ),
    highlight_color = "blue",
    highlight_alpha = 0.5,
    percentile = 0.95,
    n_grid = 7
)
ggsave( file.path(FIGURES_BASE, "Figure_SF23", "SF23_select_Dorsal.pdf"), width = 10.5, height = 10.5 )

# %% [markdown]
# #### Select area: Ventral

# %% tags=["cell-305"]
options( repr.plot.width = 10.5, repr.plot.height = 10.5 )
ventral.area <-
plot.retina3(
    human, 
    "VAX1", 
    highlight_quadrants = list(
        c(0,0),c(1,0),c(2,0),c(3,0),c(4,0),c(5,0),c(6,0),
        c(0,1),c(1,1),c(2,1),c(3,1),c(4,1),c(5,1),c(6,1),
        c(0,2),c(1,2),c(2,2),c(3,2),c(4,2),c(5,2),c(6,2)
    ),
    highlight_color = "blue",
    highlight_alpha = 0.5,
    percentile = 0.95,
    n_grid = 7
)
ggsave( file.path(FIGURES_BASE, "Figure_SF23", "SF23_select_Ventral.pdf"), width = 10.5, height = 10.5 )

# %% [markdown]
# #### Select area: DVcentral

# %% tags=["cell-307"]
options( repr.plot.width = 10.5, repr.plot.height = 10.5 )
DVcentral.area <-
plot.retina3(
    human, 
    "BMP2", 
    highlight_quadrants = list(
        c(0,3),c(1,3),c(2,3),c(3,3),c(4,3),c(5,3),c(6,3)
    ),
    highlight_color = "blue",
    highlight_alpha = 0.5,
    percentile = 0.95,
    n_grid = 7
)
ggsave( file.path(FIGURES_BASE, "Figure_SF23", "SF23_select_DVcentral.pdf"), width = 10.5, height = 10.5 )

# %% [markdown]
# #### Select area: NTcentral

# %% tags=["cell-309"]
options( repr.plot.width = 10.5, repr.plot.height = 10.5 )
NTcentral.area <-
plot.retina3(
    human, 
    "CYP1B1", 
    highlight_quadrants = list(
        c(5,0),c(5,1),c(5,2),c(5,3),c(5,4),c(5,5),c(5,6),c(5,7),c(5,8),
        c(4,0),c(4,1),c(4,2),c(4,3),c(4,4),c(4,5),c(4,6),c(4,7),c(4,8)
    ),
    highlight_color = "blue",
    highlight_alpha = 0.5,
    percentile = 0.95,
    n_grid = 9
)
ggsave( file.path(FIGURES_BASE, "Figure_SF23", "SF23_select_NTcentral.pdf"), width = 10.5, height = 10.5 )

# %% [markdown]
# ## SF23: Human pseudobulk DEG + volcano

# %% [markdown]
# #### Pseudobulk DEG

# (run_all_pseudobulk_deg and add.area.columns defined once near the top of
# the chick block and reused here for the human SF23 pipeline.)

# %% tags=["cell-314"]
# Human uses "Fovea" (matching SuppTable2 + the manuscript) rather than the
# chick "HAA" label -- the volcano facet strip then matches the table region
# name without cross-document inconsistency.
areas <-
list(
    "Fovea"     = fovea.area,
    "Temporal"  = temporal.area,
    "Nasal"     = nasal.area,
    "Dorsal"    = dorsal.area,
    "Ventral"   = ventral.area,
    "DVcentral" = DVcentral.area,
    "NTcentral" = NTcentral.area
)

# %% tags=["cell-315"]
human <- add.area.columns(human, areas)

# %% tags=["cell-316"]
# Run the analysis
deg_results <- run_all_pseudobulk_deg(human, donor_col = "sample", library_col = "library", min_cells = 50)
# deg_results

# %% tags=["cell-317"]
# Inspect the |log2FC|>1 & adj-p<0.05 human DEGs (volcano-flagged set). The
# corrected manuscript Supplementary Table 2 (full adj-p<0.05, no fold cap)
# is emitted separately by scripts/analysis/area_significant_deg.R.
deg_results %>%
dplyr::filter( abs(lfc) > 1, adj_pval < 0.05 ) %>%
dplyr::arrange( area, desc(lfc) )

# %% tags=["cell-319"]
options( repr.plot.width = 10.5, repr.plot.height = 7 )
p <-
deg_results %>%
mutate(
    flag = abs(lfc) > 1 & adj_pval < 0.05,
#     lfc = ifelse( abs(lfc) > 5, sign(lfc)*5, lfc )
) %>% {
    ggplot( ., aes( x = lfc, y = -log10( adj_pval ) ) ) +
    geom_point( 
        aes( 
            colour = flag, 
            size = sqrt(meanExp)
        ),
        alpha = 0.5
    ) +
    geom_vline( xintercept = c(-1, 1), linetype = "dashed", colour = "salmon" ) +
    geom_hline( yintercept = -log10( 0.05 ), linetype = "dashed", colour = "salmon" ) +
    # annotate(
    #     "text",
    #     x = -10,
    #     y = 18,
    #     label = "Fovea",
    #     hjust = 0,
    #     size = 5
    # ) +
    ggforce::geom_mark_circle(
        data = . %>% 
                dplyr::arrange( desc( lfc ) ) %>%
                group_by( area ) %>%
                dplyr::filter(
                    adj_pval < 0.05,
                    abs(lfc) > 1,
                    meanExp > 0.2,
                    # Drop human AC* contigs and HBA/HBB hemoglobin-cluster
                    # genes (RBC contamination) from the volcano callout set.
                    !grepl( "^(AC|HBA|HBB)", name )
                ) %>%
                dplyr::filter( 
#                    gene %in% c( "CYP26C1", "MYOF", "HHEX", "AKR1D1" ) ),
                    row_number() <= 5 # Top 5
                ),
        aes(
            label = gene,
            group = gene,
#            fill = after_scale(alpha(colour, 0.1))
        ),
        alpha = 1,
        expand = 0.02,
        con.cap = unit(0, "mm"),
        con.arrow = arrow(
            length = unit(1, "mm"),
            ends = "last",
            type = "closed"
        ),
        label.buffer = unit(1, 'mm'),
        label.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
        label.fontsize = 15,  # geom_mark_circle uses points
        label.fill = "white"
    ) +
    scale_y_continuous( expand = expansion( add = 0, mult = c(0, 0.05) ) ) +
    scale_colour_manual( 
        values = c( 
            "TRUE" = "salmon", 
            "FALSE" = "grey" 
        )
    ) +
    theme(
    #    legend.position = "none",
        text = element_text( size = 20 )
    ) +
    labs(
        x = "log2 fold change"
    ) + 
    guides(
        flag = "none",
        colour = "none",
        size = guide_legend(title = "expression\nlevel")
    ) +
    facet_wrap( area ~ ., ncol = 7 )
}
options( repr.plot.width = 49, repr.plot.height = 7 )
p
ggsave( p, file=file.path(FIGURES_BASE, "Figure_SF23", "SF23_volcano.pdf"), width = 49, height = 7 )
