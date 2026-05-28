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
# Inputs (override via the DATA_DIR env var if your layout differs):
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
}))

DATA_DIR <- Sys.getenv("DATA_DIR", unset = "data")
OUT_DIR  <- "docs/supplementary_tables"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Region -> (defining marker, NT index range lo:hi, DV index range lo:hi,
# n_grid). The human high-acuity region is named Fovea to match anatomy
# (chick: HAA). Gates = rectangular quadrant ranges replicating the figure
# script's `plot.retina3` `highlight_quadrants` selections at the same
# `n_grid` for each region.
specs <- list(
  chick = list(
    HAA       = list(m="CYP26C1", ni=c(6,7), di=c(4,5), n=10),
    Temporal  = list(m="FOXD1",   ni=c(0,2), di=c(0,6), n=7),
    Nasal     = list(m="FOXG1",   ni=c(5,6), di=c(0,6), n=7),
    Dorsal    = list(m="ALDH1A1", ni=c(0,6), di=c(4,6), n=7),
    Ventral   = list(m="VAX1",    ni=c(0,6), di=c(0,2), n=7),
    DVcentral = list(m="BMP2",    ni=c(0,6), di=c(3,3), n=7),
    NTcentral = list(m="CYP1B1",  ni=c(5,6), di=c(0,8), n=9)),
  human = list(
    Fovea     = list(m="CYP26C1", ni=c(3,4), di=c(4,5), n=10),
    Temporal  = list(m="FOXD1",   ni=c(0,2), di=c(0,6), n=7),
    Nasal     = list(m="FOXG1",   ni=c(4,6), di=c(0,6), n=7),
    Dorsal    = list(m="ALDH1A1", ni=c(0,6), di=c(4,6), n=7),
    Ventral   = list(m="VAX1",    ni=c(0,6), di=c(0,2), n=7),
    DVcentral = list(m="BMP2",    ni=c(0,6), di=c(3,3), n=7),
    NTcentral = list(m="CYP1B1",  ni=c(4,5), di=c(0,8), n=9)))

gate <- function(dv, nt, sp){
  # Cells in the rectangular quadrant range (strict <, as validated for HAA / Fovea).
  dr <- range(dv); nr <- range(nt)
  dvp <- dr[1] + (dr[2]-dr[1]) * (0:sp$n) / sp$n
  ntp <- nr[1] + (nr[2]-nr[1]) * (0:sp$n) / sp$n
  (nt >= ntp[sp$ni[1]+1] & nt < ntp[sp$ni[2]+2]) &
  (dv >= dvp[sp$di[1]+1] & dv < dvp[sp$di[2]+2])
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
  de$meanExp <- gene_mean[de$gene]
  de
}

do_species <- function(obj, sp_name, donor_col, library_col){
  obj <- NormalizeData(obj, verbose = FALSE)
  md  <- obj@meta.data; dv <- md$DV.Score; nt <- md$NT.Score
  L   <- LayerData(obj, layer = "data"); gene_mean <- rowMeans(L)
  ct  <- LayerData(obj, layer = "counts")
  out <- list()
  for (rg in names(specs[[sp_name]])) {
    sp  <- specs[[sp_name]][[rg]]
    sel <- gate(dv, nt, sp)
    # Forcing function against gate drift: HAA / Fovea cell counts are
    # validated against the published Fig 6H / S23. The figure script
    # (scripts/figures/fig6h_sfig23_area_deg.R) encodes the same quadrants
    # in a different parameterization (corner-pair vs. range-tuple); if
    # anyone edits one without the other, this stopifnot trips before a
    # silently wrong Supp Table gets written. Expected counts: chick HAA
    # = 5,971 (= the published Fig 6H HAA arm); human Fovea = 2,236.
    if (sp_name == "chick" && rg == "HAA")   stopifnot(sum(sel) == 5971L)
    if (sp_name == "human" && rg == "Fovea") stopifnot(sum(sel) == 2236L)
    de <- region_deg(ct, md, gene_mean, sel, donor_col, library_col)
    if (is.null(de)) {
      cat(sprintf("[%s/%s] SKIPPED (a level emptied by min_cells)\n", sp_name, rg))
      next
    }
    sig <- de %>%
      dplyr::transmute(
        gene,
        region = rg,
        average_expression = round(meanExp, 4),
        log2FC = round(lfc, 3),
        pval     = signif(pval, 3),
        adj_pval = signif(adj_pval, 3),
        direction = ifelse(lfc > 0, "enriched", "de-enriched")) %>%
      dplyr::filter(de$adj_pval < 0.05) %>%
      dplyr::arrange(desc(log2FC))
    mk <- sp$m
    mk_lfc <- if (mk %in% de$gene) round(de$lfc[de$gene == mk], 2) else NA
    up <- sig %>% dplyr::filter(log2FC > 0)
    cat(sprintf(
      "[%s/%-9s] gate %5d cells | sig %3d up / %3d down | marker %s lfc=%s (up-rank %s)\n",
      sp_name, rg, sum(sel), sum(sig$log2FC > 0), sum(sig$log2FC < 0), mk, mk_lfc,
      ifelse(mk %in% up$gene, which(up$gene == mk), "NA")))
    out[[rg]] <- sig
  }
  res <- dplyr::bind_rows(out)
  tbl_num <- if (sp_name == "chick") "1" else "2"
  f <- file.path(OUT_DIR, sprintf("SuppTable%s_%s_area_significant.csv", tbl_num, sp_name))
  write_csv(res, f)
  cat(sprintf("[%s] wrote %s : %d significant region-gene rows (%d unique genes)\n\n",
              sp_name, basename(f), nrow(res), dplyr::n_distinct(res$gene)))
}

# ----- chick ----------------------------------------------------------------
cat("loading chick...\n")
chick <- readRDS(file.path(DATA_DIR, "20250604_02_fabp7.rds"))
if ("RNA" %in% SeuratObject::Assays(chick)) SeuratObject::DefaultAssay(chick) <- "RNA"
chick <- tryCatch(JoinLayers(chick), error = function(e) chick)
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
rn <- if (!is.null(rownames(meta)) && all(colnames(counts) %in% rownames(meta))) {
  rownames(meta)
} else {
  bc <- names(which(sapply(meta, function(c) mean(colnames(counts) %in% as.character(c)) > 0.99)))[1]
  as.character(meta[[bc]])
}
rownames(meta) <- rn
meta <- meta[colnames(counts), , drop = FALSE]
human <- CreateSeuratObject(counts = counts, meta.data = meta)
do_species(human, "human", donor_col = "sample", library_col = "library")
cat("DONE\n")
