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
# # Human Retina Preprocessing: Load Public Data, Subset RPCs, Score

# %%
source("00_utils.R")

# %% [markdown]
# ## Load human Seurat object

# %% tags=["cell-263"]
human <- readRDS( "../../data/20240815_01_human_retina_extended.rds" )

# %% [markdown]
# ## Subset to RPC clusters

# %% tags=["cell-266"]
human <- subset( human, subset = harmony_clusters.leiden.0.4 %in% c(1,5,7) )

# %% [markdown]
# ## Compute DV/NT scores

# %% tags=["cell-269"]
human <- 
  human %>%
  JoinLayers() %>%
  AddModuleScore(
      features = list(
            c(
                # Ventral genes
                "CHRDL1",  # Chrdl1
                "ALDH1A3",  # Aldh1a3
                "VAX1"   # Vax1
#                 "VAX1", "CHRDL1", "ALDH1A3"# , "EPHB2", "EPHB3", "PAX2" ,
# #                "SMOC1", "BMPR1B", "SHISA2", "RDH10", "FRAS1", "SORCS2", "IRX3"
            )
      ),
      ctrl = 5,
      name = 'Ventral.Score'
  ) %>%
  AddModuleScore(
      features = list(
          c(
                # Dorsal genes
                "EFNB1", # Efnb1
                "TBX3", # Tbx3
                "TBX5", # Tbx5
                "EFNB2", # Efnb2
                "TBX2", # Tbx2
                "ALDH1A1"  # Aldh1a1
                # "TBX5",  "TBX2",   "TBX3", "ALDH1A1",  "EFNB2", "EFNB1"# , "GDF6", "NR2F2" ,
                # # "BAMBI", "CRABP1", "ID2",  "ID3",     "CCND2", "ENOX1", "ANOS1", "UNC5B", "LYPD6"
            )
      ),
      ctrl = 6,
      name = 'Dorsal.Score'
  ) %>%
  AddModuleScore(
      features = list(
            c(
                # Nasal genes
                "FOXG1",
                "HMX1",
                "EFNA5",
                "EFNA2"
                # "FOXG1", "HMX1", "EFNA5", "EFNA2" # listed by Heer
            )
      ),
      ctrl = 5,
      name = 'Nasal.Score'
  ) %>%
  AddModuleScore(
      features = list(
          c(
        # Temporal genes
                "FOXD1",
                "EPHA3"
              # "FOXD1", "EPHA3"
            )
      ),
      ctrl = 5,
      name = 'Temporal.Score'
  )

# %% tags=["cell-270"]
human$DV.Score <- human$Dorsal.Score1 - human$Ventral.Score1
human$NT.Score <- human$Nasal.Score1 - human$Temporal.Score1

# %% [markdown]
# ## Export human RPC h5ad

# %%
export_seurat( human, dir_path = "../../data/20250604human.RPC", assay = "RNA", prefix = "20250604human.RPC_" )
