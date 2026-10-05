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
# # Supplementary Tables 1 & 2: area-specific differentially expressed genes
#
# Supplementary Table 1 (chick) and Supplementary Table 2 (human): for each
# spatial region of Figure S24, every gene with adjusted p < 0.05 in the
# region-vs-rest pseudobulk test, with average expression (mean log-normalized
# expression over all cells), log2 fold change (in-region vs. rest), p value,
# adjusted p value and direction (enriched / de-enriched). No fold-change
# cut-off is applied; the 2-fold line is only the highlight in the S24
# volcanoes.
#
# Method (same as `sfig24_area_deg.R`): per region, pseudobulk counts by
# library x area x donor, groups with >= 50 cells kept, glmGamPoi fit
# `~ library + area`, contrast `area1`. Donor = `genotype` (chick) / `sample`
# (human).
#
# Regions are rectangular quadrant ranges on an n x n grid over the DV.Score
# and NT.Score ranges -- the same gates `sfig24_area_deg.R` draws, written
# here as (NT index range, DV index range, n). The HAA gate counts are
# asserted (5,971 chick / 2,236 human cells) in both scripts, so an edit to one
# copy of a gate without the other stops the run.
#
# Inputs: `data/20250604_02_fabp7.rds` (chick) and
# `data/20250604_human_RPC.h5ad` (human) from GEO GSE322831.
# Outputs: `figures/Tables/SuppTable1_chick_area_significant.csv` (6,383
# region-gene rows) and `figures/Tables/SuppTable2_human_area_significant.csv`
# (1,396 rows). Columns: gene, region, average_expression, log2FC, pval,
# adj_pval, direction.

# %%
source(file.path(here::here(), "scripts/preprocessing/00_utils.R"))
suppressWarnings(suppressMessages({ library(SingleCellExperiment); library(glmGamPoi); library(Matrix) }))

OUTDIR <- file.path(here::here(), "figures", "Tables")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
MIN_CELLS <- 50

# region -> marker used for the S24 selection map (reported only), NT index
# range, DV index range, grid size n. The human high-acuity gate
# (temporal-of-center) is labelled "HAA" in Supplementary Table 2.
specs <- list(
    chick = list(
        HAA       = list(m = "CYP26C1", ni = c(6, 7), di = c(4, 5), n = 10),
        Temporal  = list(m = "FOXD1",   ni = c(0, 2), di = c(0, 6), n = 7),
        Nasal     = list(m = "FOXG1",   ni = c(5, 6), di = c(0, 6), n = 7),
        Dorsal    = list(m = "ALDH1A1", ni = c(0, 6), di = c(4, 6), n = 7),
        Ventral   = list(m = "VAX1",    ni = c(0, 6), di = c(0, 2), n = 7),
        DVcentral = list(m = "BMP2",    ni = c(0, 6), di = c(3, 3), n = 7),
        NTcentral = list(m = "CYP1B1",  ni = c(5, 6), di = c(0, 8), n = 9)),
    human = list(
        HAA       = list(m = "CYP26C1", ni = c(3, 4), di = c(4, 5), n = 10),
        Temporal  = list(m = "FOXD1",   ni = c(0, 2), di = c(0, 6), n = 7),
        Nasal     = list(m = "FOXG1",   ni = c(4, 6), di = c(0, 6), n = 7),
        Dorsal    = list(m = "ALDH1A1", ni = c(0, 6), di = c(4, 6), n = 7),
        Ventral   = list(m = "VAX1",    ni = c(0, 6), di = c(0, 2), n = 7),
        DVcentral = list(m = "BMP2",    ni = c(0, 6), di = c(3, 3), n = 7),
        NTcentral = list(m = "CYP1B1",  ni = c(4, 5), di = c(0, 8), n = 9)))
expected_haa <- c(chick = 5971L, human = 2236L)

# Cells inside the quadrant range (lower bound inclusive, upper bound strict,
# as in plot.retina3).
gate <- function(dv, nt, sp) {
    dr <- range(dv); nr <- range(nt)
    dvp <- dr[1] + (dr[2] - dr[1]) * (0:sp$n) / sp$n
    ntp <- nr[1] + (nr[2] - nr[1]) * (0:sp$n) / sp$n
    (nt >= ntp[sp$ni[1] + 1] & nt < ntp[sp$ni[2] + 2]) & (dv >= dvp[sp$di[1] + 1] & dv < dvp[sp$di[2] + 2])
}

# %%
region_deg <- function(counts, md, gene_mean, sel, donor_col, library_col, min_cells = MIN_CELLS) {
    area <- factor(as.integer(sel), levels = c(0, 1))          # "1" = in region -> coefficient area1
    minimal <- CreateSeuratObject(counts = counts,
        meta.data = data.frame(library = as.character(md[[library_col]]), area = area,
                               donor = as.character(md[[donor_col]]), row.names = rownames(md)))
    pb <- glmGamPoi::pseudobulk(as.SingleCellExperiment(minimal), group_by = vars(library, area, donor))
    grp <- minimal@meta.data %>% dplyr::count(library, area, donor, name = "ncell") %>%
        dplyr::mutate(dplyr::across(c(library, area, donor), as.character))
    cd <- as.data.frame(colData(pb)); cd$.ord <- seq_len(nrow(cd))
    nvec <- dplyr::left_join(dplyr::mutate(cd, dplyr::across(c(library, area, donor), as.character)),
                             grp, by = c("library", "area", "donor"))
    nvec <- nvec$ncell[order(nvec$.ord)]
    pb <- pb[, nvec >= min_cells]
    colData(pb)$area <- droplevels(factor(colData(pb)$area))
    colData(pb)$library <- droplevels(factor(colData(pb)$library))
    if (nlevels(colData(pb)$area) < 2) return(NULL)
    fit <- glm_gp(pb, design = ~ library + area)
    de <- as.data.frame(test_de(fit, contrast = "area1")); de$gene <- rownames(pb)
    de$meanExp <- gene_mean[de$gene]
    de
}

area_tables <- function(obj, species, donor_col, library_col, table_no) {
    obj <- NormalizeData(obj, verbose = FALSE)
    md <- obj@meta.data; dv <- md$DV.Score; nt <- md$NT.Score
    gene_mean <- rowMeans(LayerData(obj, layer = "data")); ct <- LayerData(obj, layer = "counts")
    out <- list()
    for (rg in names(specs[[species]])) {
        sp <- specs[[species]][[rg]]; sel <- gate(dv, nt, sp)
        if (rg == "HAA") stopifnot(sum(sel) == expected_haa[[species]])
        de <- region_deg(ct, md, gene_mean, sel, donor_col, library_col)
        if (is.null(de)) { cat(sprintf("[%s/%s] skipped: a level emptied by min_cells\n", species, rg)); next }
        sig <- de %>%
            dplyr::filter(adj_pval < 0.05) %>%
            dplyr::transmute(gene, region = rg, average_expression = round(meanExp, 4),
                             log2FC = round(lfc, 3), pval = signif(pval, 3), adj_pval = signif(adj_pval, 3),
                             direction = ifelse(lfc > 0, "enriched", "de-enriched")) %>%
            dplyr::arrange(desc(log2FC))
        up <- sig %>% dplyr::filter(log2FC > 0)
        mk_lfc <- if (sp$m %in% de$gene) round(de$lfc[de$gene == sp$m], 2) else NA
        cat(sprintf("[%s/%-9s] gate %5d cells | adj-p<0.05: %4d enriched / %4d de-enriched | marker %s log2FC=%s (enriched rank %s)\n",
                    species, rg, sum(sel), sum(sig$log2FC > 0), sum(sig$log2FC < 0), sp$m, mk_lfc,
                    ifelse(sp$m %in% up$gene, which(up$gene == sp$m), "NA")))
        out[[rg]] <- sig
    }
    res <- dplyr::bind_rows(out)
    f <- file.path(OUTDIR, sprintf("SuppTable%s_%s_area_significant.csv", table_no, species))
    write_csv(res, f)
    cat(sprintf("[%s] wrote %s: %d region-gene rows (%d unique genes)\n\n",
                species, basename(f), nrow(res), dplyr::n_distinct(res$gene)))
}

# %% [markdown]
# ## Supplementary Table 1: chick

# %%
chick <- readRDS(file.path(here::here(), "data", "20250604_02_fabp7.rds"))
if ("RNA" %in% SeuratObject::Assays(chick)) SeuratObject::DefaultAssay(chick) <- "RNA"
chick <- tryCatch(JoinLayers(chick), error = function(e) chick)
area_tables(chick, "chick", donor_col = "genotype", library_col = "library", table_no = "1")
rm(chick); invisible(gc())

# %% [markdown]
# ## Supplementary Table 2: human

# %%
human <- read_h5ad_as_seurat(file.path(here::here(), "data", "20250604_human_RPC.h5ad"))
area_tables(human, "human", donor_col = "sample", library_col = "library", table_no = "2")
