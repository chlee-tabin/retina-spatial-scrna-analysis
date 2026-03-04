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
# # Figure S18: DV and NT Scores from Mouse and Human Retina scRNA-seq

# %%
source("../preprocessing/00_utils.R")

# %% [markdown]
# ## Output directory setup

# %%
FIGURES_BASE <- file.path(dirname(sys.frame(1)$ofile %||% "."), "..", "..", "figures")
dir.create(file.path(FIGURES_BASE, "Figure_SF18"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(FIGURES_BASE, "Figure_SF18", "variants"), recursive = TRUE, showWarnings = FALSE)

# %% [markdown]
# ## SF18C: Human DV axial expression

# %% [markdown]
# ### Figure 81

# %% tags=["cell-272"]
p <- plot_axial_expression(
    human, 
    genes = c( "TBX5", "TBX2", "TBX3", "ALDH1A1", "EFNB2", "EFNB1", "VAX1", "CHRDL1", "ALDH1A3" ),
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
# ggsave(p, filename = "figures2/Fig81.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18C_raw.pdf"), width = 14, height = 9)

# %% [markdown]
# ### Figure 81 (normalized)

# %% tags=["cell-274"]
p <- plot_axial_expression(
    human, 
    genes = c( "TBX5", "TBX2", "TBX3", "ALDH1A1", "EFNB2", "EFNB1", "VAX1", "CHRDL1", "ALDH1A3" ),
    axis_score = "DV.Score",
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
# ggsave(p, filename = "figures2/Fig81_normalized.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "SF18C_human_dv_markers.pdf"), width = 14, height = 9)

# %% [markdown]
# ### Figure 82

# %% tags=["cell-276"]
p <- plot_axial_expression(
    human, 
    genes = c( "TBX5", "TBX2", "TBX3", "ALDH1A1", "EFNB2", "EFNB1", "VAX1", "CHRDL1", "ALDH1A3", "CYP26C1", "CYP26A1" ),
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
# ggsave(p, filename = "figures2/Fig82.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18C_validation_raw.pdf"), width = 14, height = 9)

# %% tags=["cell-277"]
p <- plot_axial_expression(
    human, 
    genes = c(  "CYP26C1", "CYP26A1", "FGF8" ),
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
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18C_reference.png"), width = 14, height = 9)

# %% [markdown]
# ### Figure 82 (normalized)

# %% tags=["cell-279"]
p <- plot_axial_expression(
    human, 
    genes = c( "TBX5", "TBX2", "TBX3", "ALDH1A1", "EFNB2", "EFNB1", "VAX1", "CHRDL1", "ALDH1A3", "CYP26C1", "CYP26A1" ),
    axis_score = "DV.Score",
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
# ggsave(p, filename = "figures2/Fig82_normalized.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18C_normalized.pdf"), width = 14, height = 9)

# %% tags=["cell-280"]
p <- plot_axial_expression(
    human, 
    genes = c(  "CYP26C1", "CYP26A1" ),
    reference_genes = c( "TBX5", "TBX2", "TBX3", "ALDH1A1", "EFNB2", "EFNB1", "VAX1", "CHRDL1", "ALDH1A3" ),
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
    label_at = "max",        # Labels at maximum points (or "endpoints")
    add_histogram = TRUE
)
options(repr.plot.width = 14, repr.plot.height = 9)
p
# ggsave(p, filename = "figures2/Fig82_normalized_reference.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18C_normalized_reference.pdf"), width = 14, height = 9)

# %% [markdown]
# ## SF18D: Human NT axial expression

# %% [markdown]
# ### Figure 83

# %% tags=["cell-282"]
p <- plot_axial_expression(
    human, 
    genes = c( "FOXG1", "HMX1", "EFNA5", "EFNA2", "FOXD1", "EPHA3"
 ),
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
# ggsave(p, filename = "figures2/Fig83.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18D_raw.pdf"), width = 14, height = 9)

# %% [markdown]
# ### Figure 83 (normalized)

# %% tags=["cell-285"]
p <- plot_axial_expression(
    human, 
    genes = c( "FOXG1", "HMX1", "EFNA5", "EFNA2", "FOXD1", "EPHA3"
 ),
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
# ggsave(p, filename = "figures2/Fig83_normalized.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "SF18D_human_nt_markers.pdf"), width = 14, height = 9)

# %% [markdown]
# ### Figure 84

# %% tags=["cell-287"]
p <- plot_axial_expression(
    human, 
    genes = c( "FOXG1", "HMX1", "EFNA5", "EFNA2", "FOXD1", "EPHA3", "CYP26A1", "CYP26C1", "CYP1B1" ),
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
# ggsave(p, filename = "figures2/Fig84.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18D_validation_raw.pdf"), width = 14, height = 9)

# %% tags=["cell-288"]
p <- plot_axial_expression(
    human, 
    genes = c( "CYP26A1", "CYP26C1", "CYP1B1" ),
    reference_genes = c( "FOXG1", "HMX1", "EFNA5", "EFNA2", "FOXD1", "EPHA3" ),
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
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18D_reference.png"), width = 14, height = 9)

# %% [markdown]
# ### Figure 84 (normalized)

# %% tags=["cell-290"]
p <- plot_axial_expression(
    human, 
    genes = c( "FOXG1", "HMX1", "EFNA5", "EFNA2", "FOXD1", "EPHA3", "CYP26A1", "CYP26C1", "CYP1B1" ),
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
    add_histogram = TRUE
)
options(repr.plot.width = 14, repr.plot.height = 9)
p
# ggsave(p, filename = "figures2/Fig84_normalized.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18D_normalized.pdf"), width = 14, height = 9)

# %% tags=["cell-291"]
p <- plot_axial_expression(
    human, 
    genes = c( "CYP26A1", "CYP26C1", "CYP1B1" ),
    reference_genes = c( "FOXG1", "HMX1", "EFNA5", "EFNA2", "FOXD1", "EPHA3" ),
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
    add_histogram = TRUE
)
options(repr.plot.width = 14, repr.plot.height = 9)
p
# ggsave(p, filename = "figures2/Fig84_normalized_reference.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18D_normalized_reference.pdf"), width = 14, height = 9)

# %% [markdown]
# ## SF18A: Mouse DV axial expression

# %% [markdown]
# ### Figure 85

# %% tags=["cell-330"]
p <- plot_axial_expression(
    mouse, 
    genes = c( "Tbx5", "Tbx2", "Tbx3", "Aldh1a1", "Vax1", "Chrdl1", "Aldh1a3" ),
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
# ggsave(p, filename = "figures2/Fig85.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18A_raw.pdf"), width = 14, height = 9)

# %% [markdown]
# ### Figure 85 (normalized)

# %% tags=["cell-332"]
p <- plot_axial_expression(
    mouse, 
    genes = c( "Tbx5", "Tbx2", "Tbx3", "Aldh1a1", "Vax1", "Chrdl1", "Aldh1a3" ),
    axis_score = "DV.Score",
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
# ggsave(p, filename = "figures2/Fig85_normalized.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "SF18A_mouse_dv_markers.pdf"), width = 14, height = 9)

# %% [markdown]
# ### Figure 86

# %% tags=["cell-334"]
p <- plot_axial_expression(
    mouse, 
    genes = c( "Tbx5", "Tbx2", "Tbx3", "Aldh1a1", "Efnb2", "Efnb1", "Vax1", "Chrdl1", "Aldh1a3", "Cyp26c1", "Cyp26a1", "Bmp2" ),
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
# ggsave(p, filename = "figures2/Fig86.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18A_validation_raw.pdf"), width = 14, height = 9)

# %% tags=["cell-335"]
p <- plot_axial_expression(
    mouse, 
    genes = c(  "Cyp26c1", "Cyp26a1", "Bmp2" ),
    reference_genes = c( "Tbx5", "Tbx2", "Tbx3", "Aldh1a1", "Efnb2", "Efnb1", "Vax1", "Chrdl1", "Aldh1a3" ),
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
# ggsave(p, filename = "figures2/Fig86_reference.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18A_reference.pdf"), width = 14, height = 9)

# %% [markdown]
# ### Figure 86 (normalized)

# %% tags=["cell-337"]
p <- plot_axial_expression(
    mouse, 
    genes = c(  "Cyp26c1", "Cyp26a1", "Bmp2" ),
    reference_genes = c( "Tbx5", "Tbx2", "Tbx3", "Aldh1a1", "Efnb2", "Efnb1", "Vax1", "Chrdl1", "Aldh1a3" ),
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
    label_at = "max",        # Labels at maximum points (or "endpoints")
    label_method = "circle_repel",  # NEW: Circles + repelled labels
    label_nudge_y = 0.08,           # Magnitude of nudge (auto-adjusts direction)
    label_force = 10,               # Increase for more label separation
    add_histogram = TRUE
)
options(repr.plot.width = 14, repr.plot.height = 9)
p
# ggsave(p, filename = "figures2/Fig86_normalized_reference.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18A_normalized_reference.pdf"), width = 14, height = 9)

# %% tags=["cell-338"]
p <- plot_axial_expression(
    mouse, 
    genes = c( "Tbx5", "Tbx2", "Tbx3", "Aldh1a1", "Efnb2", "Efnb1", "Vax1", "Chrdl1", "Aldh1a3", "Cyp26c1", "Cyp26a1", "Bmp2" ),
    axis_score = "DV.Score",
    normalize = TRUE,
    conditional_lines = TRUE,
    density_threshold = 50,
    line_alpha = 0.5,
    n_bins = 50,             
    smooth_span = 0.5,
    add_labels = TRUE,
    label_at = "max",        # Labels at maximum points (or "endpoints")
    label_method = "circle_repel",  # NEW: Circles + repelled labels
    label_nudge_y = 0.08,           # Magnitude of nudge (auto-adjusts direction)
    label_force = 10,               # Increase for more label separation
    add_histogram = TRUE
)
options(repr.plot.width = 14, repr.plot.height = 9)
p
# ggsave(p, filename = "figures2/Fig86_normalized.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18A_normalized.pdf"), width = 14, height = 9)

# %% [markdown]
# ## SF18B: Mouse NT axial expression

# %% [markdown]
# ### Figure 87

# %% tags=["cell-340"]
p <- plot_axial_expression(
    mouse, 
    genes = c( "Foxg1", "Hmx1", "Efna5", "Efna2", "Foxd1", "Epha3" ),
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
# ggsave(p, filename = "figures2/Fig87.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18B_raw.pdf"), width = 14, height = 9)

# %% [markdown]
# ### Figure 87 (normalized)

# %% tags=["cell-342"]
p <- plot_axial_expression(
    mouse, 
    genes = c( "Foxg1", "Hmx1", "Efna5", "Efna2", "Foxd1", "Epha3" ),
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
# ggsave(p, filename = "figures2/Fig87_normalized.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "SF18B_mouse_nt_markers.pdf"), width = 14, height = 9)

# %% [markdown]
# ### Figure 88

# %% tags=["cell-344"]
p <- plot_axial_expression(
    mouse, 
    genes = c( "Foxg1", "Hmx1", "Efna5", "Efna2", "Foxd1", "Epha3", "Cyp26a1", "Cyp26c1", "Cyp1b1" ),
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
# ggsave(p, filename = "figures2/Fig88.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18B_validation_raw.pdf"), width = 14, height = 9)

# %% tags=["cell-345"]
p <- plot_axial_expression(
    mouse, 
    genes = c( "Cyp26a1", "Cyp26c1", "Cyp1b1" ),
    reference_genes = c( "Foxg1", "Hmx1", "Efna5", "Efna2", "Foxd1", "Epha3" ),
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
# ggsave(p, filename = "figures2/Fig88_reference.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18B_reference.pdf"), width = 14, height = 9)

# %% [markdown]
# ### Figure 88 (normalized)

# %% tags=["cell-347"]
p <- plot_axial_expression(
    mouse, 
    genes = c( "Foxg1", "Hmx1", "Efna5", "Efna2", "Foxd1", "Epha3", "Cyp26a1", "Cyp26c1", "Cyp1b1" ),
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
    label_nudge_y = 0.08,           # Magnitude of nudge (auto-adjusts direction)
    label_force = 10,               # Increase for more label separation
    add_histogram = TRUE
)
options(repr.plot.width = 14, repr.plot.height = 9)
p
# ggsave(p, filename = "figures2/Fig88_normalized.
