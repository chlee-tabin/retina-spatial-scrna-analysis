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
# # Mouse Retina Preprocessing: Load Public Data, Subset RPCs, Score

# %%
source("00_utils.R")

# %% [markdown]
# ## Load mouse Seurat object

# %% tags=["cell-321"]
mouse <- readRDS( "../../data/20240815_02_mouse_integrated_RPCs.rds" )

# %% [markdown]
# ## Subset to E14/E16 stages

# %% tags=["cell-325"]
mouse <- subset(
    mouse,
    subset = age %in% c("E14", "E16") | library == "Balasubramanian"
)
mouse@meta.data %>%
dplyr::count( library, age )

# %% [markdown]
# ## Compute DV/NT scores

# %% tags=["cell-326"]
mouse <- 
  mouse %>%
  JoinLayers() %>%
  AddModuleScore(
      features = list(
            c(
                # Dorsal genes
                "Chrdl1",  # Chrdl1
                "Aldh1a3",  # Aldh1a3
                "Vax1"   # Vax1
            )
      ),
      ctrl = 3,
      name = 'Ventral.Score'
  ) %>%
  AddModuleScore(
      features = list(
          c(
                # Ventral genes
                "Tbx2", # Tbx2
                "Tbx3", # Tbx3
                "Tbx5", # Tbx5
                "Aldh1a1"  # Aldh1a1
            )
      ),
      ctrl = 4,
      name = 'Dorsal.Score'
  ) %>%
  AddModuleScore(
      features = list(
            c(
                # Nasal genes
                "Foxg1", # Foxg1
                "Hmx1", # Hmx1
                "Efna5", # Efna5
                "Efna2"  # Efna2
            )
      ),
      ctrl = 4,
      name = 'Nasal.Score'
  ) %>%
  AddModuleScore(
      features = list(
          c(
                # Temporal genes
                "Foxd1",
                "Epha3"
              # "FOXD1", "EPHA3"
            )
      ),
      ctrl = 2,
      name = 'Temporal.Score'
  )
mouse$DV.Score <- mouse$Dorsal.Score1 - mouse$Ventral.Score1
mouse$NT.Score <- mouse$Nasal.Score1 - mouse$Temporal.Score1

# %% [markdown]
# ## Export mouse RPC h5ad

# %% [markdown]
# #### Re-export correct age range

# %% tags=["cell-328"]
export_seurat( mouse, dir_path = "20251007mouse.RPC", assay = "RNA", prefix = "20251007mouse.RPC_" )
