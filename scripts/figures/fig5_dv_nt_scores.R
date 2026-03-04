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
# # Figure 5B-E: DV and NT Score Generation from Chicken scRNA-seq

# %%
source("../preprocessing/00_utils.R")

# %% [markdown]
# ## Output directory setup

# %%
FIGURES_BASE <- file.path(dirname(sys.frame(1)$ofile %||% "."), "..", "..", "figures")
dir.create(file.path(FIGURES_BASE, "Figure5"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(FIGURES_BASE, "Figure5", "variants"), recursive = TRUE, showWarnings = FALSE)

# %% [markdown]
# ## F5B: DV marker spatial expression (9 anchor genes)

# %% tags=["cell-177"]
# Example usage:
p <- plot_axial_expression(
    fabp7, 
    genes = c("TBX5", "TBX2", "TBX3", "ALDH1A1", "EFNB2", "EFNB1", "VAX1", "CHRDL1", "ALDH1A3"),
    axis_score = "DV.Score",
    normalize = FALSE,
    conditional_lines = TRUE,
    density_threshold = 50,
    line_alpha = 0.5,
    n_bins = 50,             
    smooth_span = 0.5,
    add_labels = TRUE,
    label_at = "max",        # Labels at maximum points (or "endpoints")
    add_histogram = TRUE
)
options(repr.plot.width = 14, repr.plot.height = 9)
p
# ggsave(p, filename = "figures2/Fig51.png", width = 14, height = 9)
# ggsave(p, filename = "figures2/Fig51.pdf", width = 14, height = 9)

# %% tags=["cell-179"]
# Example usage:
p <- plot_axial_expression(
    fabp7, 
    genes = c("TBX5", "TBX2", "TBX3", "ALDH1A1", "EFNB2", "EFNB1", "VAX1", "CHRDL1", "ALDH1A3"),
    axis_score = "DV.Score",
    normalize = TRUE,
    conditional_lines = TRUE,
    density_threshold = 50,
    line_alpha = 0.8,
    n_bins = 50,             
    smooth_span = 0.5,
    add_labels = TRUE,
    label_at = "max",        # Labels at maximum points (or "endpoints")
    label_method = "circle_repel",
    label_fontsize = 20,
    add_histogram = TRUE
)
options(repr.plot.width = 14, repr.plot.height = 9)
p
ggsave(p, filename = file.path(FIGURES_BASE, "Figure5", "F5B_dv_markers.png"), width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure5", "F5B_dv_markers.svg"), width = 14, height = 9)
# ggsave(p, filename = "figures2/Fig51_normalized.png", width = 14, height = 9)
# ggsave(p, filename = "figures2/Fig51_normalized.pdf", width = 14, height = 9)

# %% [markdown]
# ## F5D: DV validation with independent genes + reference overlay

# %% tags=["cell-182"]
# Example usage:
p <- plot_axial_expression(
    fabp7,
    genes = c("TBX5", "TBX2", "TBX3", "ALDH1A1", "EFNB2", "EFNB1", "VAX1", "CHRDL1", "ALDH1A3", "FGF8", "CYP26C1", "BMP2"),
    axis_score = "DV.Score",
    normalize = FALSE,
    conditional_lines = TRUE,
    density_threshold = 50,
    line_alpha = 0.5,
    n_bins = 50,             
    smooth_span = 0.5,
    add_labels = TRUE,
    label_at = "max",        # Labels at maximum points (or "endpoints")
    add_histogram = TRUE
)
options(repr.plot.width = 14, repr.plot.height = 9)
p
# ggsave(p, filename = "figures2/Fig52.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure5", "variants", "F5D_raw.pdf"), width = 14, height = 9)

# %% tags=["cell-183"]
p <- plot_axial_expression(
    fabp7, 
    genes = c( "FGF8", "CYP26C1", "BMP2" ),
    reference_genes = c( "TBX5", "TBX2", "TBX3", "ALDH1A1", "EFNB2", "EFNB1", "VAX1", "CHRDL1", "ALDH1A3" ),
    reference_color = "grey70",
    reference_alpha = 0.3,
    axis_score = "DV.Score",
    normalize = FALSE,
    conditional_lines = TRUE,
    density_threshold = 50,
    line_alpha = 0.5,
    n_bins = 50,             
    smooth_span = 0.5,
    add_labels = TRUE,
    label_at = "max",        # Labels at maximum points (or "endpoints")
    add_histogram = TRUE
)
options(repr.plot.width = 14, repr.plot.height = 9)
p
# ggsave(p, filename = "figures2/Fig52_reference.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure5", "variants", "F5D_reference.pdf"), width = 14, height = 9)

# %% [markdown]
# ### Figure 52 (normalized)

# %% tags=["cell-185"]
# Example usage:
p <- plot_axial_expression(
    fabp7,
    genes = c("TBX5", "TBX2", "TBX3", "ALDH1A1", "EFNB2", "EFNB1", "VAX1", "CHRDL1", "ALDH1A3", "FGF8", "CYP26C1", "BMP2"),
    axis_score = "DV.Score",
    normalize = TRUE,
    conditional_lines = TRUE,
    density_threshold = 50,
    line_alpha = 0.8,
    n_bins = 50,             
    smooth_span = 0.5,
    add_labels = TRUE,
    label_at = "max",        # Labels at maximum points (or "endpoints")
    add_histogram = TRUE
)
options(repr.plot.width = 14, repr.plot.height = 9)
p
# ggsave(p, filename = "figures2/Fig52_normalized.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure5", "variants", "F5D_normalized.pdf"), width = 14, height = 9)

# %% tags=["cell-186"]
p_circle_repel <- plot_axial_expression(
    fabp7, 
    genes = c("FGF8", "CYP26C1", "BMP2"),
    reference_genes = c("TBX5", "TBX2", "TBX3", "ALDH1A1", "EFNB2", 
                       "EFNB1", "VAX1", "CHRDL1", "ALDH1A3"),
    reference_color = "grey70",
    reference_alpha = 0.3,
    axis_score = "DV.Score",
    normalize = TRUE,
    conditional_lines = TRUE,
    density_threshold = 50,
    line_alpha = 0.8,
    n_bins = 50,
    smooth_span = 0.5,
    add_labels = TRUE,
    label_at = "max",
    label_method = "circle_repel",  # NEW: Circles + repelled labels
    label_fontsize = 20,
    label_nudge_y = 0.08,           # Magnitude of nudge (auto-adjusts direction)
    label_force = 10,               # Increase for more label separation
    add_histogram = TRUE
)
options(repr.plot.width = 14, repr.plot.height = 9)
p_circle_repel
# ggsave(p_circle_repel, filename = "figures2/Fig52_normalized_reference.png", width = 14, height = 9)
# ggsave(p_circle_repel, filename = "figures2/Fig52_normalized_reference.pdf", width = 14, height = 9)
ggsave(p_circle_repel, filename = file.path(FIGURES_BASE, "Figure5", "F5D_dv_validation.png"), width = 14, height = 9)
ggsave(p_circle_repel, filename = file.path(FIGURES_BASE, "Figure5", "F5D_dv_validation.svg"), width = 14, height = 9)

# %% tags=["cell-187"]
p_circle_repel <- plot_axial_expression(
    fabp7, 
    genes = c("FGF8", "CYP26A1", "CYP26C1", "BMP2"),
    reference_genes = c("TBX5", "TBX2", "TBX3", "ALDH1A1", "EFNB2", 
                       "EFNB1", "VAX1", "CHRDL1", "ALDH1A3"),
    reference_color = "grey70",
    reference_alpha = 0.3,
    axis_score = "DV.Score",
    normalize = TRUE,
    conditional_lines = TRUE,
    density_threshold = 50,
    line_alpha = 0.5,
    n_bins = 50,
    smooth_span = 0.5,
    add_labels = TRUE,
    label_at = "max",
    label_method = "circle_repel",  # NEW: Circles + repelled labels
    label_nudge_y = 0.08,           # Magnitude of nudge (auto-adjusts direction)
    label_force = 10,               # Increase for more label separation
    add_histogram = TRUE
)
options(repr.plot.width = 14, repr.plot.height = 9)
p_circle_repel
# ggsave(p_circle_repel, filename = "figures2/Fig52_normalized_reference2.png", width = 14, height = 9)
ggsave(p_circle_repel, filename = file.path(FIGURES_BASE, "Figure5", "variants", "F5D_normalized_reference2.pdf"), width = 14, height = 9)

# %% [markdown]
# ## F5C: NT marker spatial expression (7 anchor genes)

# %% tags=["cell-211"]
# Example usage:
p <- plot_axial_expression(
    fabp7,
    genes = c("FOXG1", "SOHO-1", "HMX1", "EFNA5", "EFNA2", "FOXD1", "EPHA3"),
    axis_score = "NT.Score",
    normalize = FALSE,
    conditional_lines = TRUE,
    density_threshold = 50,
    line_alpha = 0.5,
    n_bins = 50,             
    smooth_span = 0.5,
    add_labels = TRUE,
    label_at = "max",        # Labels at maximum points (or "endpoints")
    add_histogram = TRUE
)
options(repr.plot.width = 14, repr.plot.height = 9)
p
# ggsave(p, filename = "figures2/Fig53.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure5", "variants", "F5C_raw.pdf"), width = 14, height = 9)

# %% tags=["cell-213"]
# Example usage:
p <- plot_axial_expression(
    fabp7,
    genes = c("FOXG1", "SOHO-1", "HMX1", "EFNA5", "EFNA2", "FOXD1", "EPHA3"),
    axis_score = "NT.Score",
    normalize = TRUE,
    conditional_lines = TRUE,
    density_threshold = 50,
    line_alpha = 0.8,
    n_bins = 50,             
    smooth_span = 0.5,
    add_labels = TRUE,
    label_at = "max",        # Labels at maximum points (or "endpoints")
    label_method = "circle_repel",
    label_fontsize = 20,
    add_histogram = TRUE
)
options(repr.plot.width = 14, repr.plot.height = 9)
p
# ggsave(p, filename = "figures2/Fig53_normalized.png", width = 14, height = 9)
# ggsave(p, filename = "figures2/Fig53_normalized.pdf", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure5", "F5C_nt_markers.svg"), width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure5", "F5C_nt_markers.png"), width = 14, height = 9)

# %% [markdown]
# ## F5E: NT validation with independent genes + reference overlay

# %% tags=["cell-215"]
# Example usage:
p <- plot_axial_expression(
    fabp7,
    genes = c("FOXG1", "SOHO-1", "HMX1", "EFNA5", "EFNA2", "FOXD1", "EPHA3",
        "FGF8", "CYP26C1", "BMP2", "CYP1B1"),
    axis_score = "NT.Score",
    normalize = FALSE,
    conditional_lines = TRUE,
    density_threshold = 50,
    line_alpha = 0.5,
    n_bins = 50,             
    smooth_span = 0.5,
    add_labels = TRUE,
    label_at = "max",        # Labels at maximum points (or "endpoints")
    add_histogram = TRUE
)
options(repr.plot.width = 14, repr.plot.height = 9)
p
# ggsave(p, filename = "figures2/Fig54.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure5", "variants", "F5E_raw.pdf"), width = 14, height = 9)

# %% tags=["cell-216"]
p <- plot_axial_expression(
    fabp7, 
    genes = c( "FGF8", "CYP26C1", "BMP2", "CYP1B1" ),
    reference_genes = c( "FOXG1", "SOHO-1", "HMX1", "EFNA5", "EFNA2", "FOXD1", "EPHA3" ),
    reference_color = "grey70",
    reference_alpha = 0.3,
    axis_score = "NT.Score",
    normalize = FALSE,
    conditional_lines = TRUE,
    density_threshold = 50,
    line_alpha = 0.5,
    n_bins = 50,             
    smooth_span = 0.5,
    add_labels = TRUE,
    label_at = "max",        # Labels at maximum points (or "endpoints")
    add_histogram = TRUE
)
options(repr.plot.width = 14, repr.plot.height = 9)
p
ggsave(p, filename = file.path(FIGURES_BASE, "Figure5", "variants", "F5E_reference.png"), width = 14, height = 9)

# %% [markdown]
# ### Figure 54 (normalized)

# %% tags=["cell-218"]
# Example usage:
p <- plot_axial_expression(
    fabp7,
    genes = c("FOXG1", "SOHO-1", "HMX1", "EFNA5", "EFNA2", "FOXD1", "EPHA3",
        "FGF8", "CYP26C1", "BMP2", "CYP1B1"),
    axis_score = "NT.Score",
    normalize = TRUE,
    conditional_lines = TRUE,
    density_threshold = 50,
    line_alpha = 0.5,
    n_bins = 50,             
    smooth_span = 0.5,
    add_labels = TRUE,
    label_at = "max",        # Labels at maximum points (or "endpoints")
    add_histogram = TRUE
)
options(repr.plot.width = 14, repr.plot.height = 9)
p
# ggsave(p, filename = "figures2/Fig54_normalized.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure5", "variants", "F5E_normalized.pdf"), width = 14, height = 9)

# %% tags=["cell-219"]
p <- plot_axial_expression(
    fabp7, 
    genes = c( "FGF8", "CYP26C1", "BMP2", "CYP1B1" ),
    reference_genes = c( "FOXG1", "SOHO-1", "HMX1", "EFNA5", "EFNA2", "FOXD1", "EPHA3" ),
    reference_color = "grey70",
    reference_alpha = 0.3,
    axis_score = "NT.Score",
    normalize = TRUE,
    conditional_lines = TRUE,
    density_threshold = 50,
    line_alpha = 0.5,
    n_bins = 50,             
    smooth_span = 0.5,
    add_labels = TRUE,
    label_at = "max",        # Labels at maximum points (or "endpoints")
    label_method = "circle_repel",  # NEW: Circles + repelled labels
    label_fontsize = 20,
    add_histogram = TRUE
)
options(repr.plot.width = 14, repr.plot.height = 9)
p
# ggsave(p, filename = "figures2/Fig54_normalized_reference.png", width = 14, height = 9)
# ggsave(p, filename = "figures2/Fig54_normalized_reference.pdf", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure5", "F5E_nt_validation.png"), width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure5", "F5E_nt_validation.svg"), width = 14, height = 9)

# %% [markdown]
# ### Figure 54 (CYP1B1)

# %% tags=["cell-221"]
show.genes <- c( "CYP1B1", "FGF8" )
DV.region.min <- -1
DV.region.max <- 1
p6.new <-
FetchData(
    fabp7,
    vars = c(
        "DV.Score",
        "NT.Score",
        show.genes
    ),
    slot = "data"
) %>%
pivot_longer(cols = !any_of(c("DV.Score", "NT.Score")), names_to = "gene", values_to = "expression") %>%
mutate(
    DV.category = case_when(
        DV.Score > DV.region.min & DV.Score < DV.region.max ~ "DV.HAA",
        TRUE ~ "DV.rest"
    )
) %>%
dplyr::select( -DV.Score ) %>%
mutate(
    gene.group = case_when(
        gene %in% genes.nasal ~ "Nasal",
        gene %in% genes.temporal ~ "Temporal",
        gene %in% c("CYP1B1") ~ "Stripe",
        TRUE ~ "other"
    ),
    gene.label = paste0( gene, "(", DV.category, ")" )
) %>%
ggplot(aes(x = NT.Score, y = expression, colour = gene.label )) +
geom_smooth(linewidth = 1, linetype = "dashed", se = FALSE) +
expand_limits(y = 0) +
scale_y_continuous(expand = expansion(add = 0, mult = c(0, 0.05))) +
labs(y = "Average log-normalized expression") +
theme(text = element_text(size = 20))
# Build the plot to extract data AND the correct gene mapping
built_plot <- ggplot_build(p6.new)
smooth_data <- built_plot$data[[1]]
# CORRECT WAY: Extract the gene names from the plot's scale
# This gives us the mapping between group numbers and gene names
gene_mapping <- built_plot$plot$scales$get_scales("colour")$get_limits()
# If the above doesn't work, try this alternative:
if(is.null(gene_mapping)) {
  gene_mapping <- sort(unique(built_plot$plot$data$gene.label))  # ggplot2 uses alphabetical order for characters
}
print("Gene mapping (in order of group numbers):")
print(gene_mapping)
# For each gene, find the MAXIMUM point (could be anywhere on the curve)
max_points <- smooth_data %>%
  group_by(group) %>%
  filter(!is.na(y)) %>%
  # Find the point with maximum y value for each group
  slice_max(y, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  # CORRECT MAPPING: Use the gene_mapping vector which is in the correct order
  mutate(gene.label = gene_mapping[group]) %>%
  select(gene.label, NT.Score = x, expression = y, color_hex = colour)
# Print to verify the maximum points
print("Maximum points for each gene:")
print(max_points)
# Optional: Get more detailed information about where maxima occur
max_points_detailed <- smooth_data %>%
  group_by(group) %>%
  filter(!is.na(y)) %>%
  slice_max(y, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    gene.label = gene_mapping[group],
    position_relative = (x - min(smooth_data$x)) / (max(smooth_data$x) - min(smooth_data$x)),
    position_type = case_when(
      position_relative < 0.2 ~ "left side",
      position_relative > 0.8 ~ "right side",
      TRUE ~ "middle"
    )
  ) %>%
  select(gene.label, NT.Score = x, max_expression = y, position_type)
print("Detailed maximum points:")
print(max_points_detailed)
# Now recreate the plot with markers at the maximum points
p6.new_with_labels <- FetchData(
    fabp7,
    vars = c(
        "DV.Score",
        "NT.Score",
        show.genes
    ),
    slot = "data"
) %>%
pivot_longer(cols = !any_of(c("DV.Score", "NT.Score")), names_to = "gene", values_to = "expression") %>%
mutate(
    DV.category = case_when(
        DV.Score > DV.region.min & DV.Score < DV.region.max ~ "DV.HAA",
        TRUE ~ "DV.rest"
    )
) %>%
dplyr::select( -DV.Score ) %>%
mutate(
    gene.group = case_when(
        gene %in% genes.nasal ~ "Nasal",
        gene %in% genes.temporal ~ "Temporal",
        gene %in% c("CYP1B1") ~ "Stripe",
        TRUE ~ "other"
    ),
    gene.label = paste0( gene, "(", DV.category, ")" )
) %>%
ggplot(aes(x = NT.Score, y = expression, colour = gene.label)) +
geom_smooth(linewidth = 1, linetype = "dashed", se = FALSE) +
geom_vline( xintercept =1, linetype = "dashed", colour = "salmon" ) +
# Add circles with labels at the MAXIMUM points
ggforce::geom_mark_circle(
    data = max_points,
    aes(x = NT.Score, y = expression, 
        label = gene.label,
        fill = gene.label),  # Map fill to gene for color matching
    expand = unit(2, "mm"),
    label.fontsize = 10,
    label.buffer = unit(3, "mm"),
    con.colour = "gray60",
    con.type = "straight",
    label.colour = "black",
    con.cap = unit(0, "mm"),
    alpha = 0.1  # Add transparency to the fill
) +
expand_limits(y = 0) +
scale_y_continuous(expand = expansion(add = 0, mult = c(0, 0.05))) +
scale_fill_discrete(guide = "none") +  # Hide the fill legend
labs(
    y = "Average log-normalized expression"
) +
theme(
    text = element_text(size = 20),
    legend.position = "none"  # Hide legend since we have labels
)
# Display the plot
options(repr.plot.width = 14, repr.plot.height = 7)
p6.new_with_labels
# Save the plot
ggsave(
    p6.new_with_labels, 
    filename=file.path(FIGURES_BASE, "Figure5", "variants", "F5E_CYP1B1.png"),
    width = 14, height = 7
)
Fig54_CYP1B1_legend <- glue::glue(
"
# Figure54 (CYB1B1) legend: Supplemental figure related to Figure 3
Expression of TBX2 related to FGF8 expression in HAA region. 
The aggregated log-normalized expression curve along the DV-axis is separated by the naso-temporal region aligning to HAA.
Note that TBX2 expression tapers at the point of FGF8 expression increase (dotted red line)
"
)
# Save to markdown file
writeLines(Fig54_CYP1B1_legend, file.path(FIGURES_BASE, "Figure5", "variants", "F5E_CYP1B1_legend.md"))
