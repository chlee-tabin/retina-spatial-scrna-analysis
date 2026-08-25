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
# # Chick DV/NT Scoring: Spatial Axis Score Computation for RPCs

# %%
source("00_utils.R")

# %%
# Output of 01_chick_preprocessing.R; also deposited in GEO GSE322831
# (download GSE322831_20250604_01_retina.rds and strip the GSE322831_ prefix; see README).
if (exists("retina")) {
    message("retina: using in-session object (", ncol(retina), " cells); GEO object of record is data/20250604_01_retina.rds")
} else {
    retina <- readRDS(file.path(here::here(), "data", "20250604_01_retina.rds"))
}

# %% [markdown]
# ## Subset RPCs and compute DV score

# %% tags=["cell-142"]
fabp7 <- subset( retina, subset = annotation == "RPC" )
fabp7

# %% tags=["cell-143"]
fabp7 <- 
  fabp7 %>%
  JoinLayers() %>%
  AddModuleScore(
      features = list(
            c(
                "VAX1", "CHRDL1", "ALDH1A3"# , "EPHB2", "EPHB3", "PAX2" ,
#                "SMOC1", "BMPR1B", "SHISA2", "RDH10", "FRAS1", "SORCS2", "IRX3"
            )
      ),
      ctrl = 5,
      name = 'Ventral.Score'
  ) %>%
  AddModuleScore(
      features = list(
          c(
                "TBX5",  "TBX2",   "TBX3", "ALDH1A1",  "EFNB2", "EFNB1"# , "GDF6", "NR2F2" ,
                # "BAMBI", "CRABP1", "ID2",  "ID3",     "CCND2", "ENOX1", "ANOS1", "UNC5B", "LYPD6"
            )
      ),
      ctrl = 6,
      name = 'Dorsal.Score'
  )

# %% tags=["cell-144"]
fabp7$DV.Score <- fabp7$Dorsal.Score1 - fabp7$Ventral.Score1

# %% [markdown]
# ## Compute NT score

# %% tags=["cell-148"]
fabp7 <- 
  fabp7 %>%
  AddModuleScore(
      features = list(
            c(
                "FOXG1", "SOHO-1", "HMX1", "EFNA5", "EFNA2" # nasal marker set curated by H.N.V. Joisher
            )
      ),
      ctrl = 5,
      name = 'Nasal.Score'
  ) %>%
  AddModuleScore(
      features = list(
          c(
              "FOXD1", "EPHA3"
            )
      ),
      ctrl = 5,
      name = 'Temporal.Score'
  )

# %% tags=["cell-149"]
fabp7$NT.Score <- fabp7$Nasal.Score1 - fabp7$Temporal.Score1

# %% [markdown]
# ## Checkpoint 3: save scored RPC object

# %% tags=["cell-230"]
tictoc::tic()
saveRDS( fabp7, file=glue::glue( "../../data/{intermediate.prefix}02_fabp7.rds" ) )
tictoc::toc()

# %% [markdown]
# ## Export RPC MEX for Python (h5ad assembled Python-side)
# Uses `export_seurat()` from `00_utils.R`

# %% tags=["cell-234"]
export_seurat( fabp7, "../../data/20250604chick.RPC", prefix="20250604chickRPC_" )
