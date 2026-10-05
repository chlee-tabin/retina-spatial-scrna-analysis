#!/usr/bin/env Rscript
# =============================================================================
# mouse_cr9_01_build.R  --  Stage 1 of the CR9 mouse retina atlas
#
# Build + integrate + comprehensively annotate the 10-library CR9.0.1 / GRCm39
# mouse retina atlas, score DV/NT, and cross-check de-novo calls against the
# Clark (GSE118614) umap2_CellType labels and the Balasubramanian graphclust.
#
# Mirrors scripts/preprocessing/01_chick_preprocessing.R. Built with R 4.4.3:
# load SingleCellExperiment (+glmGamPoi) BEFORE Seurat.
#
# Outputs (all under <SCR>/atlas/):
#   mouse_cr9_integrated.rds          full annotated integrated object
#   bala_filter_report.tsv            SCP-control overlap
#   doublet_summary.tsv               scDblFinder per library
#   qc_filter_report.tsv              cells kept per library + reason
#   markers_leiden_{0.2,0.4}.tsv      presto wilcoxauc
#   cluster_signature_zscore.tsv      cluster x cell-type-signature z-matrix
#   annotation_map.tsv                leiden_0.2 cluster -> annotation (argmax)
#   composition_{cluster,stage,lib}.tsv
#   clark_crosstab_cluster.tsv        leiden_0.2 x Clark umap2_CellType (clark libs)
#   clark_crosstab_annotation.tsv     annotation x Clark umap2_CellType
#   clark_join_report.tsv             per-library Clark match rate
#   bala_graphclust_crosstab.tsv
#   cyp26_gm32342_by_{cluster,annotation,stage,annostage}.tsv
#   full atlas MEX export -> <SCR>/atlas/export_full/
# =============================================================================

suppressMessages({
  library(SingleCellExperiment)   # env note: before Seurat
  library(glmGamPoi)
  library(Matrix)
  library(tidyverse)
  library(Seurat)
  library(harmony)
  library(presto)
  library(scDblFinder)
  library(scater)
})

set.seed(1234)
options(future.globals.maxSize = 12e9)
options(Seurat.object.assay.version = "v5")

SCR        <- Sys.getenv("WORKDIR")
if (!nzchar(SCR)) stop("Set WORKDIR to the mouse re-alignment working directory (see README).")
COUNTS     <- file.path(SCR, "counts")
OUT        <- file.path(SCR, "atlas")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
SCP_CTRL   <- Sys.getenv("SCP_CTRL", unset = file.path(SCR, "SCP1618", "SCP1618_control_bcs.txt"))
CLARK_META <- Sys.getenv("CLARK_META", unset = file.path(SCR, "GSE118614_meta", "GSE118614_barcodes.tsv.gz"))
BALA_GRAPH <- Sys.getenv("BALA_GRAPH", unset = file.path(SCR, "GSE139904", "analysis", "clustering", "graphclust", "clusters.csv"))

log <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")

# ---- compact MEX exporter (self-contained copy of 00_utils.R::export_seurat) ----
export_mex <- function(s, dir_path, assay = "RNA", prefix = "") {
  dir.create(dir_path, showWarnings = FALSE, recursive = TRUE)
  ja     <- JoinLayers(s[[assay]])
  counts <- LayerData(ja, layer = "counts")
  data   <- LayerData(ja, layer = "data")
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
    emb <- Embeddings(s[[red]])
    colnames(emb) <- paste0(red, "_", colnames(emb))
    write.table(emb, file.path(dir_path, paste0(prefix, red, "_embeddings.tsv")),
                sep = "\t", quote = FALSE, row.names = TRUE)
  }
  log("MEX export complete:", dir_path)
}

# ---- leiden availability probe (fail fast, before the expensive load) ----
# RETICULATE_PYTHON must select a Python environment with leidenalg.
# The batch wrapper exports it from the active Python interpreter.
rp <- Sys.getenv("RETICULATE_PYTHON", unset = "")
if (nzchar(rp)) reticulate::use_python(rp, required = TRUE)
leiden_ok <- tryCatch({ reticulate::import("leidenalg"); TRUE }, error = function(e) FALSE)
log("reticulate python:", tryCatch(reticulate::py_config()$python, error = function(e) "??"))
log("leidenalg importable:", leiden_ok)
if (!leiden_ok) stop("leidenalg not importable -- Seurat algorithm=4 (Leiden) will fail. Aborting before load.")

# ---- checkpoint / resume guard ----
# The heavy block (load -> QC -> integrate -> cluster -> markers) saves a
# pre-annotation checkpoint. Downstream annotation bugs then re-run cheaply.
CKPT <- file.path(OUT, "mouse_cr9_integrated_preanno.rds")
if (file.exists(CKPT)) {
  log("RESUME: loading pre-annotation checkpoint", CKPT)
  retina <- readRDS(CKPT)
} else {

# ---- manifest ----
manifest <- tribble(
  ~crid,                  ~stage,  ~source,
  "clark_E14_rep1",       "E14",   "clark",
  "clark_E14_rep2",       "E14",   "clark",
  "clark_E16",            "E16",   "clark",
  "clark_E18_rep2",       "E18",   "clark",
  "clark_E18_rep3",       "E18",   "clark",
  "clark_P0",             "P0",    "clark",
  "bala_E13p5_WT",        "E13.5", "bala",
  "wu_E13p5_atoh7het",    "E13.5", "wu",
  "logiudice_E15p5_rep1", "E15.5", "logiudice",
  "logiudice_E15p5_rep2", "E15.5", "logiudice"
)

# ---- load ----
obj_list <- list()
for (i in seq_len(nrow(manifest))) {
  crid    <- manifest$crid[i]
  mtx_dir <- file.path(COUNTS, crid, "outs", "filtered_feature_bc_matrix")
  stopifnot(dir.exists(mtx_dir))
  m <- Read10X(mtx_dir)
  s <- CreateSeuratObject(counts = m, project = crid, min.cells = 0, min.features = 0)
  s$cellranger_id <- crid
  s$library       <- crid
  s$stage         <- manifest$stage[i]
  s$source        <- manifest$source[i]
  s$barcode       <- colnames(s)      # "<16bp>-1"
  s$provenance    <- "whole_retina"
  log(sprintf("loaded %-22s %6d cells x %d genes", crid, ncol(s), nrow(s)))
  obj_list[[crid]] <- s
}

# ---- Balasubramanian WT filter to SCP1618 Control (4,642) ----
ctrl_bcs <- trimws(readLines(SCP_CTRL)); ctrl_bcs <- ctrl_bcs[nzchar(ctrl_bcs)]
bala  <- obj_list[["bala_E13p5_WT"]]
n_b   <- ncol(bala)
keep  <- intersect(colnames(bala), ctrl_bcs)
ov_recount <- length(keep) / n_b
ov_ctrl    <- length(keep) / length(ctrl_bcs)
log(sprintf("bala: recount %d | SCP-control %d | overlap %d (%.1f%% recount, %.1f%% control)",
            n_b, length(ctrl_bcs), length(keep), 100 * ov_recount, 100 * ov_ctrl))
stopifnot(length(keep) > 1000)        # sanity gate -- expect a few thousand
bala <- subset(bala, cells = keep)
bala$library    <- "bala_WT"
bala$provenance <- "peripheral_sorted_v3_SCPcontrol"
obj_list[["bala_E13p5_WT"]] <- bala
write_tsv(tibble(metric = c("recount_cells", "scp_control_list", "overlap",
                            "frac_of_recount", "frac_of_control"),
                 value  = c(n_b, length(ctrl_bcs), length(keep), ov_recount, ov_ctrl)),
          file.path(OUT, "bala_filter_report.tsv"))

# ---- gene sets ----
all_genes  <- rownames(obj_list[[1]])
mito.genes <- grep("^mt-", all_genes, value = TRUE)
rbc.genes  <- intersect(c("Hba-a1","Hba-a2","Hba-x","Hbb-bh1","Hbb-bh2","Hbb-bs","Hbb-bt","Hbb-y"), all_genes)
log("mito genes:", length(mito.genes), "| rbc genes:", length(rbc.genes))

data("cc.genes.updated.2019", package = "Seurat")
s.genes   <- intersect(str_to_title(cc.genes.updated.2019$s.genes),   all_genes)
g2m.genes <- intersect(str_to_title(cc.genes.updated.2019$g2m.genes), all_genes)
log(sprintf("cc matched: S %d/%d  G2M %d/%d", length(s.genes),
            length(cc.genes.updated.2019$s.genes), length(g2m.genes),
            length(cc.genes.updated.2019$g2m.genes)))

# ---- QC per library (scDblFinder seeded; mito/rbc fractions; cell cycle) ----
qc_one <- function(s) {
  s   <- NormalizeData(s, verbose = FALSE)
  cnt <- LayerData(s, layer = "counts")
  tot <- Matrix::colSums(cnt)
  s$percent.mito <- if (length(mito.genes)) Matrix::colSums(cnt[mito.genes, , drop = FALSE]) / tot else 0
  s$percent.rbc  <- if (length(rbc.genes))  Matrix::colSums(cnt[rbc.genes,  , drop = FALSE]) / tot else 0
  set.seed(1234)                                  # scDblFinder is stochastic
  sce <- SingleCellExperiment(assays = list(counts = cnt,
                                            logcounts = LayerData(s, layer = "data")))
  sce <- scDblFinder(sce)
  s$scDblFinder.score <- sce$scDblFinder.score
  s$scDblFinder.class <- sce$scDblFinder.class
  s <- CellCycleScoring(s, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
  s
}
obj_list <- lapply(obj_list, function(s) { log("QC", unique(s$library)); qc_one(s) })

write_tsv(map_dfr(obj_list, ~ as_tibble(.x@meta.data) %>% count(library, scDblFinder.class)),
          file.path(OUT, "doublet_summary.tsv"))

# ---- merge ----
retina <- merge(obj_list[[1]], y = obj_list[-1], add.cell.ids = manifest$crid)
log("merged:", ncol(retina), "cells across", length(unique(retina$library)), "libraries")

# ---- filter: singlet + mito<0.20 + nFeature>=200 (documented thresholds) ----
MITO_MAX <- 0.20; NFEAT_MIN <- 200
md <- retina@meta.data
retina$keep <- md$scDblFinder.class == "singlet" & md$percent.mito < MITO_MAX & md$nFeature_RNA >= NFEAT_MIN
write_tsv(as_tibble(retina@meta.data) %>%
            group_by(library) %>%
            summarise(n = n(),
                      n_singlet = sum(scDblFinder.class == "singlet"),
                      n_mito_ok = sum(percent.mito < MITO_MAX),
                      n_feat_ok = sum(nFeature_RNA >= NFEAT_MIN),
                      n_keep    = sum(keep), .groups = "drop"),
          file.path(OUT, "qc_filter_report.tsv"))
retina <- subset(retina, subset = keep)
log("after filter:", ncol(retina), "cells")

# ---- normalize / HVG / scale (regress cc + mito) / PCA ----
retina <- NormalizeData(retina, verbose = FALSE)
retina <- FindVariableFeatures(retina, verbose = FALSE)
retina <- ScaleData(retina, vars.to.regress = c("G2M.Score", "S.Score", "percent.mito"), verbose = FALSE)
retina <- RunPCA(retina, npcs = 50, verbose = FALSE)
log("PCA done")

# ---- Harmony integrate by library ----
retina <- IntegrateLayers(retina, method = HarmonyIntegration,
                          orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE)
retina <- FindNeighbors(retina, reduction = "harmony", dims = 1:30, verbose = FALSE)
for (res in c(0.1, 0.2, 0.4)) {
  retina <- FindClusters(retina, resolution = res, algorithm = 4, method = "igraph",
                         cluster.name = sprintf("leiden_%.1f", res), verbose = FALSE)
  log(sprintf("Leiden res %.1f -> %d clusters", res,
              length(unique(retina@meta.data[[sprintf("leiden_%.1f", res)]]))))
}
retina <- RunUMAP(retina, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", verbose = FALSE)
retina <- JoinLayers(retina)
log("integration + clustering + UMAP done")

# ---- markers (presto wilcoxauc) ----
for (res in c("leiden_0.2", "leiden_0.4")) {
  mk <- presto::wilcoxauc(retina, res)
  write_tsv(as_tibble(mk), file.path(OUT, paste0("markers_", res, ".tsv")))
}

  saveRDS(retina, CKPT)
  log("checkpoint saved (pre-annotation):", CKPT)
}   # end heavy block / resume guard

all_genes <- rownames(retina)

# ---- comprehensive cell-type signatures (developmental mouse retina) ----
sigs <- list(
  Early_RPC   = c("Sfrp2","Fgf15","Ccnd1","Sox2","Hes1","Hes5","Hmga2","Lin28b","Rax","Vsx2","Pax6","Nes","Sox9"),
  Late_RPC    = c("Nfia","Nfib","Nfix","Ascl1","Gadd45a","Vsx2","Pax6","Sox2","Ccnd1"),
  Neurogenic  = c("Atoh7","Neurog2","Neurod1","Neurod4","Dll1","Dll3","Hes6","Btg2","Otx2","Gadd45a","Dlx1","Dlx2"),
  RGC         = c("Pou4f1","Pou4f2","Pou4f3","Isl1","Rbpms","Sncg","Nefl","Nefm","Elavl4","Gap43","Nhlh2"),
  Amacrine    = c("Tfap2a","Tfap2b","Tfap2c","Gad1","Gad2","Slc6a9","Slc32a1","Pax6","Prox1"),
  Horizontal  = c("Onecut1","Onecut2","Onecut3","Lhx1","Prox1","Tfap2b"),
  Cone_PR     = c("Crx","Otx2","Prdm1","Thrb","Rxrg","Pde6h","Gnat2","Arr3","Opn1sw","Gngt2"),
  Rod         = c("Nrl","Nr2e3","Rho","Gnat1","Pde6b","Nt5e","Recoverin"),
  Bipolar     = c("Vsx1","Cabp5","Grik1","Grm6","Trpm1","Isl1","Otx2"),
  Muller      = c("Glul","Rlbp1","Apoe","Slc1a3","Clu","Aqp4","Crym","Sox9","Vim"),
  Microglia   = c("Aif1","C1qa","C1qb","Cx3cr1","Ptprc","Tmem119","Csf1r"),
  RBC         = c("Hba-a1","Hba-a2","Hba-x","Hbb-bs","Hbb-bt","Hbb-bh1"),
  Vascular    = c("Cldn5","Pecam1","Cdh5","Kdr","Flt1","Emcn"),
  RPE         = c("Mlana","Tyr","Tyrp1","Pmel","Rpe65","Mitf","Ttr"),
  Astrocyte   = c("Gfap","Pax2","Aqp4","S100b")
)
sigs <- lapply(sigs, function(g) intersect(unique(g), all_genes))
sig_report <- tibble(signature = names(sigs), n_genes = lengths(sigs),
                     genes = map_chr(sigs, ~ paste(.x, collapse = ",")))
write_tsv(sig_report, file.path(OUT, "signature_genes_used.tsv"))

set.seed(1234)
# AddModuleScore with a NAMED list names columns "<name><list-name>" in Seurat 5
# (not "<name>1..N"), so capture the columns it actually adds rather than assume.
before_cols <- colnames(retina@meta.data)
retina <- AddModuleScore(retina, features = unname(sigs), name = "sigscore_", nbin = 20, ctrl = 50)
score_cols <- setdiff(colnames(retina@meta.data), before_cols)            # append order = signature order
score_cols <- score_cols[order(as.integer(gsub("\\D", "", score_cols)))]   # sort by trailing integer
log("AddModuleScore added:", paste(score_cols, collapse = ", "))
stopifnot(length(score_cols) == length(sigs))
names(score_cols) <- names(sigs)

# per-cluster mean module score (leiden_0.2), z-scored per signature, argmax -> annotation
clmean <- as_tibble(retina@meta.data) %>%
  group_by(leiden_0.2) %>%
  summarise(across(all_of(unname(score_cols)), ~ mean(.x)), .groups = "drop")
mat <- as.matrix(clmean[, unname(score_cols)]); rownames(mat) <- as.character(clmean$leiden_0.2)
colnames(mat) <- names(sigs)
zmat <- scale(mat)                                   # z per signature (column) across clusters
zdf  <- as_tibble(zmat, rownames = "leiden_0.2")
write_tsv(zdf, file.path(OUT, "cluster_signature_zscore.tsv"))
anno_map <- tibble(leiden_0.2 = rownames(zmat),
                   annotation_argmax = colnames(zmat)[max.col(zmat, ties.method = "first")])

# ---- curated annotation: argmax PROPOSES, Clark umap2_CellType + markers CONFIRM ----
# Discordant small clusters corrected by convergence of Clark labels + presto markers
# (clark_crosstab_cluster.tsv + markers_leiden_0.2.tsv):
#   1,2 Early+Late RPC mix (argmax Early_RPC)                         -> RPC
#   4   Clark PR-precursors+Cones; Neurod1/Otx2/Prdm1/Cngb3 (Rod)     -> Photoreceptor
#   7   Clark Early RPCs 928/950; Aldh1a1/Trpm3/Mecom, no Tyr/Rpe65   -> RPC  (peripheral/CMZ)
#   8   Clark Early-RPC+RPE/Margin; Pax2/Nfia/Vim/Sparc  (Astrocyte)  -> Margin_Glia (optic stalk/margin)
#   9   Clark Amacrine 470; Tfap2b/Isl1/Megf10           (Bipolar)    -> Amacrine
#  11   Clark RPE/Margin catch-all but C1qa/b/c/Tyrobp/Ctss           -> Microglia (kept; more specific than Clark)
override <- c("1"="RPC", "2"="RPC", "4"="Photoreceptor", "7"="RPC", "8"="Margin_Glia", "9"="Amacrine")
class_of <- c(RPC="RPC", Early_RPC="RPC", Late_RPC="RPC", Neurogenic="Neurogenic",
              RGC="Neuron", Amacrine="Neuron", Horizontal="Neuron",
              Cone_PR="Neuron", Rod="Neuron", Bipolar="Neuron", Photoreceptor="Neuron",
              Muller="Glia", Astrocyte="Glia", Margin_Glia="Glia", Microglia="Immune",
              RBC="NonNeural", Vascular="NonNeural", RPE="NonNeural")
anno_map$annotation <- anno_map$annotation_argmax
mm <- match(names(override), as.character(anno_map$leiden_0.2))
stopifnot(!anyNA(mm))
# GUARD: `override` is keyed by leiden_0.2 cluster NUMBER, stable ONLY for this
# seeded run. A re-cluster (different inputs/versions) renumbers clusters, silently
# rebinding number -> biology. Assert each overridden cluster still carries the argmax
# we curated it against (listed above); a renumber changes that
# argmax -> loud abort instead of a silent mislabel. If this fires, re-derive the
# override from annotation_map.tsv + markers_leiden_0.2.tsv.
override_argmax <- c("1"="Early_RPC","2"="Early_RPC","4"="Rod","7"="RPE","8"="Astrocyte","9"="Bipolar")
arg_at <- setNames(anno_map$annotation_argmax, as.character(anno_map$leiden_0.2))[names(override_argmax)]
if (any(arg_at != override_argmax)) stop("annotation override keyed by cluster NUMBER is STALE -- cluster(s) ",
  paste(names(override_argmax)[arg_at != override_argmax], collapse=","),
  " no longer carry their curated argmax (clusters renumbered). Re-derive the override map.")
anno_map$annotation[mm] <- unname(override)
anno_map$annotation_class <- unname(class_of[anno_map$annotation])
stopifnot(!anyNA(anno_map$annotation_class))
write_tsv(anno_map, file.path(OUT, "annotation_map.tsv"))

# named lookup cluster -> label, indexed by each cell's cluster, then unname()
# so Seurat's $<- treats it as a positional vector (a named RHS makes Seurat try
# to match names to cell barcodes -> "No cell overlap" error).
cl2anno  <- setNames(anno_map$annotation,       as.character(anno_map$leiden_0.2))
cl2class <- setNames(anno_map$annotation_class, as.character(anno_map$leiden_0.2))
retina$annotation       <- unname(cl2anno [as.character(retina$leiden_0.2)])
retina$annotation_class <- unname(cl2class[as.character(retina$leiden_0.2)])
stopifnot(!anyNA(retina$annotation))

# ---- DV / NT scoring (VERBATIM from 04_mouse_preprocessing.R; Vax1 ventral) ----
dv_nt_genes <- c("Chrdl1","Aldh1a3","Vax1","Tbx2","Tbx3","Tbx5","Aldh1a1",
                 "Foxg1","Hmx1","Efna5","Efna2","Foxd1","Epha3")
missing_dvnt <- setdiff(dv_nt_genes, all_genes)
# DV/NT anchors are REQUIRED: a silently-dropped anchor would redefine the
# manuscript's topographic axis. Fail loudly rather than warn.
if (length(missing_dvnt)) stop("DV/NT anchor gene(s) missing -- would silently redefine the topographic axis: ", paste(missing_dvnt, collapse=","))
set.seed(1234)
retina <- retina %>%
  AddModuleScore(features = list(c("Chrdl1","Aldh1a3","Vax1")),        ctrl = 3, name = "Ventral.Score") %>%
  AddModuleScore(features = list(c("Tbx2","Tbx3","Tbx5","Aldh1a1")),   ctrl = 4, name = "Dorsal.Score") %>%
  AddModuleScore(features = list(c("Foxg1","Hmx1","Efna5","Efna2")),   ctrl = 4, name = "Nasal.Score") %>%
  AddModuleScore(features = list(c("Foxd1","Epha3")),                  ctrl = 2, name = "Temporal.Score")
retina$DV.Score <- retina$Dorsal.Score1 - retina$Ventral.Score1
retina$NT.Score <- retina$Nasal.Score1  - retina$Temporal.Score1
log("DV/NT scored")

# ---- composition tables ----
write_tsv(as_tibble(retina@meta.data) %>% count(leiden_0.2, annotation, annotation_class),
          file.path(OUT, "composition_cluster.tsv"))
write_tsv(as_tibble(retina@meta.data) %>% count(stage, annotation_class) %>%
            pivot_wider(names_from = annotation_class, values_from = n, values_fill = 0),
          file.path(OUT, "composition_stage.tsv"))
write_tsv(as_tibble(retina@meta.data) %>% count(library, stage, source),
          file.path(OUT, "composition_lib.tsv"))

# ---- Clark cross-check (umap2_CellType) ----
clark_cols <- c("rowname","barcode","sample","age","num_genes_expressed","Total_mRNAs",
                "umap_cluster","umap_coord1","umap_coord2","umap_coord3",
                "used_for_pseudotime","umap2_CellType")
clark <- read_tsv(CLARK_META, col_names = clark_cols, skip = 1, show_col_types = FALSE)
clark <- clark %>%
  mutate(bc16 = sub("^[^.]+\\.", "", barcode))   # "<16bp>-1" after stripping "<age>."
# assert Clark (age, bc16) uniqueness
dup_clark <- clark %>% count(age, bc16) %>% filter(n > 1) %>% nrow()
log("Clark (age,bc16) duplicate keys:", dup_clark)
clark_key <- clark %>% distinct(age, bc16, .keep_all = TRUE) %>%
  select(age, bc16, umap2_CellType, clark_umap_cluster = umap_cluster)

mc <- as_tibble(retina@meta.data, rownames = "cell") %>%
  filter(source == "clark") %>%
  mutate(age = stage) %>%
  left_join(clark_key, by = c("age", "barcode" = "bc16"))
join_rep <- mc %>% group_by(library) %>%
  summarise(n = n(), matched = sum(!is.na(umap2_CellType)),
            match_rate = matched / n, .groups = "drop")
write_tsv(join_rep, file.path(OUT, "clark_join_report.tsv"))
write_tsv(mc %>% filter(!is.na(umap2_CellType)) %>% count(leiden_0.2, umap2_CellType),
          file.path(OUT, "clark_crosstab_cluster.tsv"))
write_tsv(mc %>% filter(!is.na(umap2_CellType)) %>% count(annotation, umap2_CellType),
          file.path(OUT, "clark_crosstab_annotation.tsv"))

# ---- bala graphclust cross-check ----
bg <- read_csv(BALA_GRAPH, show_col_types = FALSE) %>% rename(barcode = Barcode, bala_graphclust = Cluster)
mb <- as_tibble(retina@meta.data) %>% filter(library == "bala_WT") %>%
  left_join(bg, by = "barcode")
write_tsv(mb %>% count(annotation, bala_graphclust), file.path(OUT, "bala_graphclust_crosstab.tsv"))
log("bala graphclust match:", sum(!is.na(mb$bala_graphclust)), "/", nrow(mb))

# ---- Cyp26 / Gm32342 detection tables ----
foi <- intersect(c("Cyp26c1","Cyp26a1","Cyp26b1","Gm32342","Fgf8","Pax6","Rax","Vsx2","Atoh7","Aldh1a3"), all_genes)
expr <- FetchData(retina, vars = foi, layer = "data")          # log-norm
cnts <- FetchData(retina, vars = foi, layer = "counts")        # raw counts
det  <- as_tibble(cnts > 0) %>% setNames(paste0("det_", foi))
meanl<- as_tibble(expr) %>% setNames(paste0("mean_", foi))
gx   <- bind_cols(as_tibble(retina@meta.data) %>%
                    select(leiden_0.2, annotation, annotation_class, stage, library), det, meanl)
summ <- function(df, ...) {
  df %>% group_by(...) %>%
    summarise(n_cells = n(),
              across(starts_with("det_"),  ~ mean(.x), .names = "frac_{.col}"),
              across(starts_with("mean_"), ~ mean(.x)),
              .groups = "drop")
}
write_tsv(summ(gx, leiden_0.2),                  file.path(OUT, "cyp26_gm32342_by_cluster.tsv"))
write_tsv(summ(gx, annotation),                  file.path(OUT, "cyp26_gm32342_by_annotation.tsv"))
write_tsv(summ(gx, stage),                       file.path(OUT, "cyp26_gm32342_by_stage.tsv"))
write_tsv(summ(gx, annotation, stage),           file.path(OUT, "cyp26_gm32342_by_annostage.tsv"))

# ---- save + full-atlas MEX export ----
saveRDS(retina, file.path(OUT, "mouse_cr9_integrated.rds"))
log("saved rds:", file.path(OUT, "mouse_cr9_integrated.rds"))
export_mex(retina, file.path(OUT, "export_full"), assay = "RNA", prefix = "mouse_cr9_full_")

# ---- console summary ----
cat("\n================ SUMMARY ================\n")
cat("cells:", ncol(retina), "| genes:", nrow(retina), "\n")
print(as_tibble(retina@meta.data) %>% count(annotation_class, annotation) %>% arrange(annotation_class, desc(n)))
cat("\nper-stage class composition:\n")
print(as_tibble(retina@meta.data) %>% count(stage, annotation_class) %>%
        pivot_wider(names_from = annotation_class, values_from = n, values_fill = 0))
cat("\nClark match rate:\n"); print(join_rep)
log("STAGE 1 COMPLETE")
