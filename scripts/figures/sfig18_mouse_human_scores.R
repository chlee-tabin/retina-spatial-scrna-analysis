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
source(file.path(here::here(), "scripts/preprocessing/00_utils.R"))

# %%
# Requires in-session Seurat objects (see README, "Reproducibility scope"):
# - `human`: produced by source()'ing scripts/preprocessing/03_human_preprocessing.R
#   (its input is a pre-publication intermediate, not distributed).
# - `mouse`: must be the Cell Ranger 9.0.1 / GRCm39 re-aligned RPC Seurat object
#   (rebuilt from scripts/realign_mouse/ or the deposited E13.5-E16 h5ad,
#   carrying DV.Score / NT.Score in its metadata).
#   scripts/preprocessing/04_mouse_preprocessing.R yields the SUPERSEDED pre-CR9
#   4-library object — do not use it for this figure.
if (!exists("human"))
    stop("Object 'human' not found: source scripts/preprocessing/03_human_preprocessing.R in this R session first (see README).")
if (!exists("mouse"))
    stop("Object 'mouse' not found. It must be the CR9/GRCm39 re-aligned mouse RPC Seurat object ",
         "(produced from the scripts/realign_mouse pipeline; 04_mouse_preprocessing.R yields the superseded pre-CR9 object). See README.")
# exists() alone would accept 04's superseded object (which is also named `mouse`);
# 25,202 cells is that object's deterministic count — reject it outright.
# This guard checks cell count only.
if (ncol(mouse) == 25202)
    stop("`mouse` has 25,202 cells — this is 04_mouse_preprocessing.R's superseded pre-CR9 object, ",
         "not the CR9/GRCm39 re-aligned RPC object (26,505 cells, E13.5-E16). See README.")

# %% [markdown]
# ## Output directory setup

# %%
# here::here() anchors to the repo root (via the committed .here sentinel, so it
# also works in ZIP/Zenodo archives without .git); the script works under
# Rscript or source() from any working directory inside the repository.
FIGURES_BASE <- file.path(here::here(), "figures")
dir.create(file.path(FIGURES_BASE, "Figure_SF18"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(FIGURES_BASE, "Figure_SF18", "variants"), recursive = TRUE, showWarnings = FALSE)

# %% [markdown]
# ## SF18A: Human DV axial expression

# %% [markdown]
# ### SF18A: Marker expression

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
# ggsave(p, filename = "figures2/SF18A.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18A_raw.pdf"), width = 14, height = 9)

# %% [markdown]
# ### SF18A: Marker expression (normalized)

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
# ggsave(p, filename = "figures2/SF18A_normalized.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "SF18A_human_dv_markers.pdf"), width = 14, height = 9)

# %% [markdown]
# ### SF18A': Validation expression

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
# ggsave(p, filename = "figures2/SF18A_prime.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18A_prime_validation_raw.pdf"), width = 14, height = 9)

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
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18A_prime_reference.png"), width = 14, height = 9)

# %% [markdown]
# ### SF18A': Validation expression (normalized)

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
# ggsave(p, filename = "figures2/SF18A_prime_normalized.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18A_prime_normalized.pdf"), width = 14, height = 9)

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
# ggsave(p, filename = "figures2/SF18A_prime_normalized_reference.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18A_prime_human_dv_validation.pdf"), width = 14, height = 9)

# %% [markdown]
# ## SF18B: Human NT axial expression

# %% [markdown]
# ### SF18B: Marker expression

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
# ggsave(p, filename = "figures2/SF18B.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18B_raw.pdf"), width = 14, height = 9)

# %% [markdown]
# ### SF18B: Marker expression (normalized)

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
# ggsave(p, filename = "figures2/SF18B_normalized.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "SF18B_human_nt_markers.pdf"), width = 14, height = 9)

# %% [markdown]
# ### SF18B': Validation expression

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
# ggsave(p, filename = "figures2/SF18B_prime.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18B_prime_validation_raw.pdf"), width = 14, height = 9)

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
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18B_prime_reference.png"), width = 14, height = 9)

# %% [markdown]
# ### SF18B': Validation expression (normalized)

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
# ggsave(p, filename = "figures2/SF18B_prime_normalized.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18B_prime_normalized.pdf"), width = 14, height = 9)

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
# ggsave(p, filename = "figures2/SF18B_prime_normalized_reference.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18B_prime_human_nt_validation.pdf"), width = 14, height = 9)

# %% [markdown]
# ## SF18C: Mouse DV axial expression

# %% [markdown]
# ### SF18C: Marker expression

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
# ggsave(p, filename = "figures2/SF18C_mouse_dv_markers.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18C_mouse_dv_raw.pdf"), width = 14, height = 9)

# %% [markdown]
# ### SF18C: Marker expression (normalized)

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
# ggsave(p, filename = "figures2/SF18C_mouse_dv_markers_normalized.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "SF18C_mouse_dv_markers.pdf"), width = 14, height = 9)

# %% [markdown]
# ### SF18C: Validation expression

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
# ggsave(p, filename = "figures2/SF18C_mouse_dv_validation.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18C_mouse_dv_validation_raw.pdf"), width = 14, height = 9)

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
# ggsave(p, filename = "figures2/SF18C_mouse_dv_validation_reference.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18C_mouse_dv_reference.pdf"), width = 14, height = 9)

# %% [markdown]
# ### SF18C: Validation expression (normalized)

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
# ggsave(p, filename = "figures2/SF18C_mouse_dv_validation_normalized_reference.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18C_mouse_dv_normalized_reference.pdf"), width = 14, height = 9)

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
# ggsave(p, filename = "figures2/SF18C_mouse_dv_validation_normalized.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18C_mouse_dv_normalized.pdf"), width = 14, height = 9)

# %% [markdown]
# ## SF18C: Mouse NT axial expression

# %% [markdown]
# ### SF18C: Marker expression

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
# ggsave(p, filename = "figures2/SF18C_mouse_nt_markers.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18C_mouse_nt_raw.pdf"), width = 14, height = 9)

# %% [markdown]
# ### SF18C: Marker expression (normalized)

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
# ggsave(p, filename = "figures2/SF18C_mouse_nt_markers_normalized.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "SF18C_mouse_nt_markers.pdf"), width = 14, height = 9)

# %% [markdown]
# ### SF18C: Validation expression

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
# ggsave(p, filename = "figures2/SF18C_mouse_nt_validation.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18C_mouse_nt_validation_raw.pdf"), width = 14, height = 9)

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
# ggsave(p, filename = "figures2/SF18C_mouse_nt_validation_reference.png", width = 14, height = 9)
ggsave(p, filename = file.path(FIGURES_BASE, "Figure_SF18", "variants", "SF18C_mouse_nt_reference.pdf"), width = 14, height = 9)

# %% [markdown]
# ### SF18C: Validation expression (normalized)

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
# ggsave(p, filename = "figures2/SF18C_mouse_nt_validation_normalized.
