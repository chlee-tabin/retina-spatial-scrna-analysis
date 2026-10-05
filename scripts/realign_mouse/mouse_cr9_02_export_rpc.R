#!/usr/bin/env Rscript
# =============================================================================
# mouse_cr9_02_export_rpc.R  --  Stage 2 of the CR9 mouse retina atlas
#
# Load the Stage-1 integrated+annotated object, extract the RPC population
# (unified annotation; RPC class), confirm marker gating, and export both the
# RPC subset and the full atlas as MEX for the Python h5ad assembler.
# Also writes RPC-restricted Cyp26c1 / Gm32342 detection tables (apples-to-
# apples vs chick CYP26C1 0.92%) and an E13-E16 cross-species-aligned subset.
#
# DV/NT scores are fixed-gene module scores computed on the full atlas in
# Stage 1 (embedding-independent) -> inherited unchanged by the RPC subset.
# =============================================================================

suppressMessages({
  library(SingleCellExperiment); library(Matrix); library(tidyverse)
  library(Seurat); library(presto)
})
set.seed(1234)
options(future.globals.maxSize = 12e9)

SCR <- Sys.getenv("WORKDIR")
if (!nzchar(SCR)) stop("Set WORKDIR to the mouse re-alignment working directory (see README).")
OUT <- file.path(SCR, "atlas")
log <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")

export_mex <- function(s, dir_path, assay = "RNA", prefix = "") {
  dir.create(dir_path, showWarnings = FALSE, recursive = TRUE)
  ja <- JoinLayers(s[[assay]])
  counts <- LayerData(ja, layer = "counts"); data <- LayerData(ja, layer = "data")
  rc_gz <- file.path(dir_path, paste0(prefix, "raw_counts.mtx.gz"))
  nd_gz <- file.path(dir_path, paste0(prefix, "normalized_data.mtx.gz"))
  if (!file.exists(rc_gz)) { Matrix::writeMM(counts, sub("\\.gz$", "", rc_gz)); system2("gzip", c("-f", sub("\\.gz$", "", rc_gz))) } else log("skip existing", rc_gz)
  if (!file.exists(nd_gz)) { Matrix::writeMM(data,   sub("\\.gz$", "", nd_gz)); system2("gzip", c("-f", sub("\\.gz$", "", nd_gz))) } else log("skip existing", nd_gz)
  write.table(rownames(counts), file.path(dir_path, paste0(prefix, "features.tsv")),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(colnames(counts), file.path(dir_path, paste0(prefix, "barcodes.tsv")),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(s@meta.data, file.path(dir_path, paste0(prefix, "metadata.tsv")),
              sep = "\t", quote = FALSE, row.names = TRUE)
  for (red in Reductions(s)) {
    emb <- Embeddings(s[[red]]); colnames(emb) <- paste0(red, "_", colnames(emb))
    write.table(emb, file.path(dir_path, paste0(prefix, red, "_embeddings.tsv")),
                sep = "\t", quote = FALSE, row.names = TRUE)
  }
  log("MEX export:", dir_path, "(", ncol(counts), "cells )")
}

retina <- readRDS(file.path(OUT, "mouse_cr9_integrated.rds"))
log("loaded atlas:", ncol(retina), "cells")
retina <- JoinLayers(retina)

# ---- RPC extraction (unified annotation) ----
rpc <- subset(retina, subset = annotation_class == "RPC")
log("RPC cells:", ncol(rpc), sprintf("(%.1f%% of atlas)", 100 * ncol(rpc) / ncol(retina)))

# marker-gate sanity (Pax6/Sox2/Vsx2/Rax positive; Atoh7/Neurod4 low)
gate_genes <- intersect(c("Pax6","Sox2","Vsx2","Rax","Hes1","Atoh7","Neurod4","Neurod1"), rownames(rpc))
gate <- FetchData(rpc, vars = gate_genes, layer = "counts")
gate_tbl <- tibble(gene = gate_genes,
                   frac_pos = map_dbl(gate_genes, ~ mean(gate[[.x]] > 0)),
                   mean_logn = map_dbl(gate_genes, ~ mean(FetchData(rpc, .x, layer = "data")[[1]])))
write_tsv(gate_tbl, file.path(OUT, "rpc_marker_gate.tsv"))
print(gate_tbl)

# ---- RPC subset composition + DV/NT coverage by stage ----
write_tsv(as_tibble(rpc@meta.data) %>% count(stage, annotation),
          file.path(OUT, "rpc_composition_stage.tsv"))
write_tsv(as_tibble(rpc@meta.data) %>% group_by(stage) %>%
            summarise(n = n(),
                      DV_mean = mean(DV.Score), DV_sd = sd(DV.Score),
                      DV_iqr  = IQR(DV.Score),
                      NT_mean = mean(NT.Score), NT_sd = sd(NT.Score),
                      NT_iqr  = IQR(NT.Score), .groups = "drop"),
          file.path(OUT, "rpc_dvnt_by_stage.tsv"))

# ---- RPC-restricted Cyp26 / Gm32342 (apples-to-apples vs chick 0.92%) ----
foi <- intersect(c("Cyp26c1","Cyp26a1","Cyp26b1","Gm32342","Fgf8","Aldh1a3"), rownames(rpc))
cn <- FetchData(rpc, vars = foi, layer = "counts")
dn <- FetchData(rpc, vars = foi, layer = "data")
det <- as_tibble(cn > 0) %>% setNames(paste0("det_", foi))
mnl <- as_tibble(dn)     %>% setNames(paste0("mean_", foi))
gx  <- bind_cols(as_tibble(rpc@meta.data) %>% select(stage, annotation, library), det, mnl)
summ <- function(df, ...) df %>% group_by(...) %>%
  summarise(n_cells = n(),
            across(starts_with("det_"),  ~ mean(.x), .names = "frac_{.col}"),
            across(starts_with("mean_"), ~ mean(.x)), .groups = "drop")
write_tsv(summ(gx),        file.path(OUT, "rpc_cyp26_overall.tsv"))
write_tsv(summ(gx, stage), file.path(OUT, "rpc_cyp26_by_stage.tsv"))
write_tsv(summ(gx, stage, annotation), file.path(OUT, "rpc_cyp26_by_stage_anno.tsv"))

# ---- export RPC MEX (full stage range) ----
export_mex(rpc, file.path(OUT, "export_rpc"), assay = "RNA", prefix = "20260528_mouse_RPC_cr9_")

# ---- E13-E16 cross-species-aligned RPC subset used for the manuscript figures ----
rpc_e13e16 <- subset(rpc, subset = stage %in% c("E13.5", "E14", "E15.5", "E16"))
log("RPC E13-E16:", ncol(rpc_e13e16), "cells")
export_mex(rpc_e13e16, file.path(OUT, "export_rpc_e13e16"), assay = "RNA", prefix = "20260528_mouse_RPC_cr9_e13e16_")

cat("\n================ STAGE 2 SUMMARY ================\n")
cat("atlas", ncol(retina), "| RPC(all-stage)", ncol(rpc), "| RPC(E13-E16)", ncol(rpc_e13e16), "\n")
print(read_tsv(file.path(OUT, "rpc_cyp26_by_stage.tsv"), show_col_types = FALSE) %>%
        select(stage, n_cells, starts_with("frac_det_Cyp26c1"), starts_with("frac_det_Gm32342")))
log("STAGE 2 COMPLETE")
