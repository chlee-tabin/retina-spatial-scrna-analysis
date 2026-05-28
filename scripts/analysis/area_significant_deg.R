#!/usr/bin/env Rscript
# Generates manuscript Supplementary Table 1 (chick, Fig. 6H) and Supplementary
# Table 2 (human, Fig. S23): the full statistically-significant area-DEG tables.
# These tables match the manuscript Supp Table 1 / 2 legend exactly: "genes
# enriched or de-enriched within spatial domains ... gene name, average
# expression, log2 fold-change between in-group and out-group cells, adjusted
# p-value, and significance." All adj-p < 0.05 genes are reported (no >= 2-fold
# headline cut -- that was the over-conservative published filter; relaxing it
# and applying the previously-unenforced `min_cells = 50` floor recovers the
# full list).
#
# Method mirrors Fig 6H: per region, region-vs-rest pseudobulk DEG (glmGamPoi
# `~ library + area`, contrast `area1`), with `min_cells = 50`. Region gates
# replicate the figure script's `plot.retina3` quadrant selections; the chick
# HAA and human fovea gates are validated against the published cell counts
# (5,971 / 2,236) with a `stopifnot` forcing function below.
#
# Run from the repo root: `Rscript scripts/analysis/area_significant_deg.R`.
# Inputs (override via the RETINA_DATA_DIR env var if your layout differs):
#   data/20250604_02_fabp7.rds         (chick RPC, from preprocessing
#                                       scripts/preprocessing/02_chick_dv_nt_scoring.R)
#   data/20250604human.RPC/            (human MEX export; barcodes/features/
#                                       metadata + raw counts mtx, with the
#                                       baked-in DV.Score / NT.Score used for
#                                       Fig. S23)
# Outputs:
#   docs/supplementary_tables/SuppTable1_chick_area_significant.csv
#   docs/supplementary_tables/SuppTable2_human_area_significant.csv

suppressWarnings(suppressMessages({
  library(SingleCellExperiment); library(glmGamPoi)
  library(Seurat); library(Matrix); library(dplyr); library(readr)
  library(here)
}))

# RETINA_DATA_DIR matches the project-wide convention (scripts/README.md and
# scripts/preprocessing/05_export_h5ad.R). here::here() anchors the default to
# the repo root regardless of working directory.
DATA_DIR <- Sys.getenv("RETINA_DATA_DIR", unset = here::here("data"))
OUT_DIR  <- here::here("docs", "supplementary_tables")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Region -> (defining marker, NT index range lo:hi, DV index range lo:hi,
# n_grid). The human high-acuity region is named Fovea to match anatomy
# (chick: HAA). Gates = rectangular quadrant ranges replicating the figure
# script's `plot.retina3` `highlight_quadrants` selections at the same
# `n_grid` for each region.
#
# `expected_n` is the published cell count for the gate (asserted as a forcing
# function below). Only HAA / Fovea are pinned in the published artefacts —
# the other 12 gates currently run with `expected_n = NA` (the per-region
# cell-count summary printed at end-of-run can be used to populate these once
# the figure-script counterpart gates are likewise pinned).
specs <- list(
  chick = list(
    HAA       = list(m="CYP26C1", ni=c(6,7), di=c(4,5), n=10, expected_n = 5971L),
    Temporal  = list(m="FOXD1",   ni=c(0,2), di=c(0,6), n=7,  expected_n = NA_integer_),
    Nasal     = list(m="FOXG1",   ni=c(5,6), di=c(0,6), n=7,  expected_n = NA_integer_),
    Dorsal    = list(m="ALDH1A1", ni=c(0,6), di=c(4,6), n=7,  expected_n = NA_integer_),
    Ventral   = list(m="VAX1",    ni=c(0,6), di=c(0,2), n=7,  expected_n = NA_integer_),
    DVcentral = list(m="BMP2",    ni=c(0,6), di=c(3,3), n=7,  expected_n = NA_integer_),
    NTcentral = list(m="CYP1B1",  ni=c(5,6), di=c(0,8), n=9,  expected_n = NA_integer_)),
  human = list(
    Fovea     = list(m="CYP26C1", ni=c(3,4), di=c(4,5), n=10, expected_n = 2236L),
    Temporal  = list(m="FOXD1",   ni=c(0,2), di=c(0,6), n=7,  expected_n = NA_integer_),
    Nasal     = list(m="FOXG1",   ni=c(4,6), di=c(0,6), n=7,  expected_n = NA_integer_),
    Dorsal    = list(m="ALDH1A1", ni=c(0,6), di=c(4,6), n=7,  expected_n = NA_integer_),
    Ventral   = list(m="VAX1",    ni=c(0,6), di=c(0,2), n=7,  expected_n = NA_integer_),
    DVcentral = list(m="BMP2",    ni=c(0,6), di=c(3,3), n=7,  expected_n = NA_integer_),
    NTcentral = list(m="CYP1B1",  ni=c(4,5), di=c(0,8), n=9,  expected_n = NA_integer_)))

gate <- function(dv, nt, spec){
  # Cells in the rectangular quadrant range (strict <, as validated for HAA / Fovea).
  dr <- range(dv); nr <- range(nt)
  dvp <- dr[1] + (dr[2]-dr[1]) * (0:spec$n) / spec$n
  ntp <- nr[1] + (nr[2]-nr[1]) * (0:spec$n) / spec$n
  (nt >= ntp[spec$ni[1]+1] & nt < ntp[spec$ni[2]+2]) &
  (dv >= dvp[spec$di[1]+1] & dv < dvp[spec$di[2]+2])
}

region_deg <- function(obj_counts, md, gene_mean, sel, donor_col, library_col, min_cells = 50){
  area <- factor(as.integer(sel), levels = c(0, 1))           # "1" = in-region -> coefficient area1
  minimal <- CreateSeuratObject(
    counts = obj_counts,
    meta.data = data.frame(library = as.character(md[[library_col]]),
                           area    = area,
                           donor   = as.character(md[[donor_col]]),
                           row.names = rownames(md)))
  sce <- as.SingleCellExperiment(minimal)
  pb  <- glmGamPoi::pseudobulk(sce, group_by = vars(library, area, donor))
  grp <- minimal@meta.data %>%
    dplyr::count(library, area, donor, name = "ncell") %>%
    dplyr::mutate(dplyr::across(c(library, area, donor), as.character))
  cd  <- as.data.frame(colData(pb)); cd$.ord <- seq_len(nrow(cd))
  nvec <- dplyr::left_join(
    dplyr::mutate(cd, dplyr::across(c(library, area, donor), as.character)),
    grp, by = c("library", "area", "donor"))
  nvec <- nvec$ncell[order(nvec$.ord)]
  pb <- pb[, nvec >= min_cells]
  colData(pb)$area    <- droplevels(factor(colData(pb)$area))
  colData(pb)$library <- droplevels(factor(colData(pb)$library))
  if (nlevels(colData(pb)$area) < 2) return(NULL)
  fit <- glm_gp(pb, design = ~ library + area)
  de  <- as.data.frame(test_de(fit, contrast = "area1"))
  de$gene    <- rownames(pb)
  de$meanExp <- gene_mean[de$gene]  # species-wide mean log-normalized expression
                                    # (matches figure-script volcano size aesthetic)
  de
}

do_species <- function(obj, sp_name, donor_col, library_col){
  # Validate the upstream object carries everything `region_deg` will reach
  # for. Failing here is far more informative than the downstream "subscript
  # out of bounds" / all-NA propagation a missing column would cause.
  required_meta <- c(library_col, donor_col, "DV.Score", "NT.Score")
  missing_meta  <- setdiff(required_meta, colnames(obj@meta.data))
  if (length(missing_meta) > 0L)
    stop(sprintf("[%s] missing required metadata columns: %s",
                 sp_name, paste(missing_meta, collapse = ", ")))
  # Value-level non-NA check: an all-NA DV.Score / NT.Score column would pass
  # the existence check above and then propagate NA selections through the
  # gate, producing degenerate empty / all-cells groupings that the HAA /
  # Fovea expected_n stopifnots would catch but the other 12 gates would not.
  if (anyNA(obj@meta.data$DV.Score) || anyNA(obj@meta.data$NT.Score))
    stop(sprintf(
      "[%s] DV.Score / NT.Score contains NA values (%d / %d non-NA of %d cells)",
      sp_name,
      sum(!is.na(obj@meta.data$DV.Score)),
      sum(!is.na(obj@meta.data$NT.Score)),
      nrow(obj@meta.data)))

  obj <- NormalizeData(obj, verbose = FALSE)
  md  <- obj@meta.data; dv <- md$DV.Score; nt <- md$NT.Score
  L   <- LayerData(obj, layer = "data"); gene_mean <- rowMeans(L)
  ct  <- LayerData(obj, layer = "counts")
  out         <- list()                          # successful per-region significant-DEG frames
  gate_counts <- integer(0)                      # per-region gate sizes (for end-of-run summary)
  skipped     <- character(0)                    # region names dropped by the min_cells gate
  for (rg in names(specs[[sp_name]])) {
    spec <- specs[[sp_name]][[rg]]
    sel  <- gate(dv, nt, spec)
    gate_counts[rg] <- sum(sel)
    # Forcing function against gate drift: HAA / Fovea cell counts are
    # validated against the published Fig 6H / S23. The figure script
    # (scripts/figures/fig6h_sfig23_area_deg.R) encodes the same quadrants
    # in a different parameterization (corner-pair vs. range-tuple); if
    # anyone edits one without the other, this stopifnot trips before a
    # silently wrong Supp Table gets written. Expected counts: chick HAA
    # = 5,971 (= the published Fig 6H HAA arm); human Fovea = 2,236.
    if (!is.na(spec$expected_n))
      stopifnot(sum(sel) == spec$expected_n)
    de <- region_deg(ct, md, gene_mean, sel, donor_col, library_col, min_cells = 50)
    if (is.null(de)) {
      skipped <- c(skipped, rg)
      cat(sprintf("[%s/%s] SKIPPED (a level emptied by min_cells)\n", sp_name, rg))
      next
    }
    sig <- de %>%
      dplyr::filter(adj_pval < 0.05) %>%
      dplyr::transmute(
        gene,
        region = rg,
        average_expression = round(meanExp, 4),
        log2FC = round(lfc, 3),
        pval     = signif(pval, 3),
        adj_pval = signif(adj_pval, 3),
        direction = ifelse(lfc > 0, "enriched", "de-enriched")) %>%
      dplyr::arrange(desc(log2FC))
    mk <- spec$m
    mk_lfc <- if (mk %in% de$gene) round(de$lfc[de$gene == mk], 2) else NA
    up <- sig %>% dplyr::filter(log2FC > 0)
    cat(sprintf(
      "[%s/%-9s] gate %5d cells | sig %3d up / %3d down | marker %s lfc=%s (up-rank %s)\n",
      sp_name, rg, sum(sel), sum(sig$log2FC > 0), sum(sig$log2FC < 0), mk, mk_lfc,
      ifelse(mk %in% up$gene, which(up$gene == mk), "NA")))
    out[[rg]] <- sig
  }
  # Forcing function: every spec'd region must have either produced a sig frame
  # (in `out`) or been recorded as a skip (in `skipped`). Anything else means
  # an exception bypassed the loop and we'd be writing a partial CSV.
  stopifnot(length(out) + length(skipped) == length(specs[[sp_name]]))
  res <- dplyr::bind_rows(out)
  tbl_num <- if (sp_name == "chick") "1" else "2"
  f <- file.path(OUT_DIR, sprintf("SuppTable%s_%s_area_significant.csv", tbl_num, sp_name))
  write_csv(res, f)
  cat(sprintf("[%s] wrote %s : %d significant region-gene rows (%d unique genes)\n",
              sp_name, basename(f), nrow(res), dplyr::n_distinct(res$gene)))
  cat(sprintf("[%s] per-region gate sizes: %s\n",
              sp_name,
              paste(names(gate_counts), gate_counts, sep="=", collapse=", ")))
  if (length(skipped) > 0L)
    cat(sprintf("[%s] regions skipped by min_cells gate: %s\n",
                sp_name, paste(skipped, collapse=", ")))
  cat("\n")
}

# ----- chick ----------------------------------------------------------------
cat("loading chick...\n")
chick <- readRDS(file.path(DATA_DIR, "20250604_02_fabp7.rds"))
# Fail-fast (not no-op) if the upstream object lacks an RNA assay -- a
# renamed assay (e.g. RNA_GRCh38, SCT) would silently flip LayerData() to a
# different matrix and produce a wrong gene_mean / wrong DEG fit.
stopifnot("RNA" %in% SeuratObject::Assays(chick))
SeuratObject::DefaultAssay(chick) <- "RNA"
# Seurat v5 splits counts/data into per-batch layers; `LayerData(..., layer =
# "counts")` returns one of them, not a joined matrix. We need the joined
# matrix downstream — let any JoinLayers failure surface rather than silently
# fall back to the un-joined object (which silently misreports gene_mean and
# the DEG fit).
chick <- JoinLayers(chick)
do_species(chick, "chick", donor_col = "genotype", library_col = "library")
rm(chick); gc()

# ----- human ----------------------------------------------------------------
cat("loading human...\n")
EXP <- file.path(DATA_DIR, "20250604human.RPC")
PFX <- "20250604human.RPC_"
counts <- ReadMtx(
  mtx      = file.path(EXP, paste0(PFX, "raw_counts.mtx.gz")),
  cells    = file.path(EXP, paste0(PFX, "barcodes.tsv")),
  features = file.path(EXP, paste0(PFX, "features.tsv")),
  feature.column = 1, cell.column = 1)
meta <- read.delim(file.path(EXP, paste0(PFX, "metadata.tsv")),
                   check.names = FALSE, stringsAsFactors = FALSE)
# Recover the barcode-keyed row names. read.delim auto-detects row.names only
# when the header has one fewer column than the data rows; if the metadata
# was written with explicit `row.names = TRUE` that's the common case, but
# fall back to a column whose values match the count-matrix barcodes if not.
rn <- if (!is.null(rownames(meta)) && all(colnames(counts) %in% rownames(meta))) {
  rownames(meta)
} else {
  match_frac <- sapply(meta, function(col) mean(colnames(counts) %in% as.character(col)))
  bc <- names(which(match_frac > 0.99))
  if (length(bc) == 0L)
    stop("Could not locate a barcode column in human metadata.tsv ",
         "(no column matches >99% of count-matrix colnames; best column ",
         "match was ", round(max(match_frac, na.rm = TRUE), 3), "). ",
         "Check the human MEX export at ", EXP)
  as.character(meta[[bc[1]]])
}
# Duplicate barcodes in `rn` would let `meta[colnames(counts), ]` silently
# keep only the first row per barcode -- so cells sharing a barcode would all
# inherit one cell's metadata (wrong library / donor / DV / NT). Catch that
# BEFORE assigning rownames.
stopifnot(!anyDuplicated(rn))
rownames(meta) <- rn
# After reorder, the metadata row order must align *exactly* with the count
# matrix columns. `identical(...)` catches both ordering mismatches and
# silent NA-rowname rows (the latter would produce `meta[NA, ]` -> NA row).
meta <- meta[colnames(counts), , drop = FALSE]
stopifnot(identical(rownames(meta), colnames(counts)))
human <- CreateSeuratObject(counts = counts, meta.data = meta)
do_species(human, "human", donor_col = "sample", library_col = "library")
cat("DONE\n")
