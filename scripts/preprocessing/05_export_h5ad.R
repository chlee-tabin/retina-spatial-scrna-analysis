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
# # Export Seurat Objects to h5ad/MEX Format
#
# This script documents the h5ad exports that serve as the interface
# between R preprocessing and Python figure scripts.
#
# **Note**: The export calls are embedded in the preprocessing scripts
# that produce each object. This script can be used to re-export
# without re-running the full pipeline.

# %%
source("00_utils.R")

DATA_DIR <- Sys.getenv("RETINA_DATA_DIR", unset = "../../data")

# %% [markdown]
# ## Export full chick dataset
# Produced by `01_chick_preprocessing.R` (checkpoint 2)

# %%
retina <- readRDS(file.path(DATA_DIR, "20250604_01_retina.rds"))
export_seurat(retina, dir_path = file.path(DATA_DIR, "20250604chick"), assay = "RNA")

# %% [markdown]
# ## Export chick RPC subset (with DV/NT scores)
# Produced by `02_chick_dv_nt_scoring.R` (checkpoint 3)

# %%
fabp7 <- readRDS(file.path(DATA_DIR, "20250604_02_fabp7.rds"))
export_seurat(fabp7, dir_path = file.path(DATA_DIR, "20250604chick.RPC"), assay = "RNA")

# %% [markdown]
# ## Export mouse RPC subset
# CR9 revision: the mouse RPC Seurat object is produced by the Cell Ranger 9.0.1 /
# GRCm39 realignment pipeline (7 libraries / 4 GEO series; see `scripts/realign/` in
# the analysis repo), superseding the prior `04_mouse_preprocessing.R` 4-library object.

# %%
# mouse <- readRDS(file.path(DATA_DIR, "mouse_cr9_RPC.rds"))   # CR9 RPC Seurat (DV/NT scored)
# export_seurat(mouse, dir_path = file.path(DATA_DIR, "20260528mouse_cr9.RPC"), assay = "RNA")

# %% [markdown]
# ## Output files consumed by Python scripts
#
# | Export directory | Python input | Consumer |
# |-----------------|-------------|----------|
# | `20250604chick/` | Not directly used by figure scripts | GEO submission |
# | `20250604chick.RPC/` | `20250604_chick_RPC.h5ad` | fig6, fig7, sfig15-17 |
# | `20250604_human_RPC.h5ad` | Same | fig7, fig8, sfig19-22 |
# | `20260528_mouse_RPC_cr9_e13e16.h5ad` | Same | fig7, sfig18, sfig19 |

# %%
sessionInfo()
