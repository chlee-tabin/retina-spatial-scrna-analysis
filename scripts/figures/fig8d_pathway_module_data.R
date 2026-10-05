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
# # Figure 8D (data): per-bin pathway-module maps, chick vs. human
#
# Computes the binned DV x NT maps behind Figure 8D from the deposited GEO
# objects and writes them as small wide CSVs (rows = `dvidx`, `ntidx`, one row
# per module; columns = occupied spatial bins) to `data/fig8d_pathway_modules/`.
# The committed copies of those CSVs let `fig8d_pathway_module_topography.R`
# render the figure without the large objects; re-running this script
# overwrites them, so `git diff data/fig8d_pathway_modules/` is the
# reproduction check.
#
# Modules (Seurat `AddModuleScore`, seed 42, `ctrl = max(n_genes, 5)`):
#
# - **FGF ligand**: the species' acuity-area FGF ligand as a single gene
#   (per-bin mean log-normalized expression) -- chick FGF8, human FGF7. The
#   choice is data-driven: the ligand survey printed below ranks every FGF
#   ligand by enrichment at the HAA / fovea gate.
# - **FGF feedback**: SPRY1/2/4, DUSP4/5/6, ETV4/5, SPRED1/2, FGFBP3, IL17RD.
# - **RA degradation**: CYP26A1/B1/C1. **RA synthesis**: ALDH1A1/2/3, RDH10.
# - **BMP ligand**: BMP2/4/7. GDF6 (a dorsal GDF-subfamily ligand) and the
#   sparse BMP5/6 are left out because their dorsal signal swamps the
#   patterning-BMP topography.
# - **BMP target**: ID1/2/3, BAMBI (written, not plotted).
# - **BMP antagonist**: SOSTDC1, NOG, GREM1, GREM2, FST. CHRDL1 (Ventroptin),
#   the ventral DV-axis antagonist, is left out because its ventral peak
#   swamps the module. Every BMP gene, including the excluded ones, is still
#   reported individually at the gate below.
#
# Grid: 15 x 15 equal-width bins over the DV.Score and NT.Score ranges; bins
# with >= 12 cells are kept.

# %%
source(file.path(here::here(), "scripts/preprocessing/00_utils.R"))
suppressWarnings(suppressMessages(library(Matrix)))

OUTDIR <- file.path(here::here(), "data", "fig8d_pathway_modules")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

NB <- 15; MINCELL <- 12

modules <- list(
    FGFfeedback   = c("SPRY1","SPRY2","SPRY4","DUSP6","DUSP4","DUSP5","ETV4","ETV5",
                      "SPRED1","SPRED2","FGFBP3","IL17RD"),
    RAdegradation = c("CYP26A1","CYP26B1","CYP26C1"),
    RAsynthesis   = c("ALDH1A1","ALDH1A2","ALDH1A3","RDH10"),
    BMPligand     = c("BMP2","BMP4","BMP7"),
    BMPtarget     = c("ID1","ID2","ID3","BAMBI"),
    BMPantagonist = c("SOSTDC1","NOG","GREM1","GREM2","FST")
)
fgf_ligands <- paste0("FGF", c(1:14, 16:23))
bmp_report  <- c("BMP2","BMP4","BMP5","BMP6","BMP7","GDF6","ID1","ID2","ID3","MSX1","MSX2",
                 "SMAD6","SMAD7","BAMBI","CHRDL1","SOSTDC1","GREM1","NOG","FST","GREM2")

# HAA (chick, central-nasal) / fovea (human, temporal-of-center) gate: four
# (NT, DV) quadrants of a 10 x 10 grid -- the same gate as the area-DEG HAA arm
# (sfig24_area_deg.R and supptables1_2_area_deg.R).
gates <- list(chick = list(c(7,5), c(6,5), c(7,4), c(6,4)),
              human = list(c(3,5), c(4,5), c(3,4), c(4,4)))

# %%
gate_cells <- function(dv, nt, quads, n = 10) {
    dr <- range(dv); nr <- range(nt)
    dvp <- dr[1] + (dr[2] - dr[1]) * (0:n) / n; ntp <- nr[1] + (nr[2] - nr[1]) * (0:n) / n
    sel <- rep(FALSE, length(dv))
    for (q in quads) sel <- sel | (nt >= ntp[q[1] + 1] & nt < ntp[q[1] + 2] &
                                   dv >= dvp[q[2] + 1] & dv < dvp[q[2] + 2])
    sel
}

# Mean in-gate vs. rest (log-normalized), Wilcoxon p, and peak bin per gene.
gate_report <- function(L, genes, sel, peak_of) {
    bind_rows(lapply(intersect(genes, rownames(L)), function(g) {
        x <- L[g, ]; mi <- mean(x[sel]); mo <- mean(x[!sel])
        p <- tryCatch(suppressWarnings(wilcox.test(x[sel], x[!sel])$p.value), error = function(e) NA)
        data.frame(gene = g, mean_in = round(mi, 4), mean_out = round(mo, 4), diff = round(mi - mo, 4),
                   log2FC = round(log2((mi + 1e-4) / (mo + 1e-4)), 3), p = signif(p, 3), peak = peak_of(x))
    })) %>% arrange(desc(diff))
}

module_maps <- function(obj, species, ligand) {
    obj <- NormalizeData(obj, verbose = FALSE); genes <- rownames(obj)
    for (m in names(modules)) {
        g <- intersect(modules[[m]], genes)
        message(sprintf("[%s] %-14s genes present: %s", species, m, paste(g, collapse = ", ")))
        if (length(g) >= 2) { set.seed(42); obj <- AddModuleScore(obj, features = list(g), name = m, ctrl = max(length(g), 5)) }
    }
    md <- obj@meta.data; L <- LayerData(obj, layer = "data"); dv <- md$DV.Score; nt <- md$NT.Score
    dvb <- cut(dv, seq(min(dv), max(dv), length.out = NB + 1), labels = FALSE, include.lowest = TRUE)
    ntb <- cut(nt, seq(min(nt), max(nt), length.out = NB + 1), labels = FALSE, include.lowest = TRUE)
    binid <- (dvb - 1) * NB + ntb; tab <- table(binid)
    good <- sort(as.integer(names(tab)[tab >= MINCELL])); keep <- binid %in% good
    binf <- factor(binid[keep], levels = good)
    smap <- function(s) as.numeric(tapply(s[keep], binf, mean)[as.character(good)])
    mcol <- function(m) { col <- paste0(m, "1"); if (col %in% colnames(md)) smap(md[[col]]) else rep(NA, length(good)) }
    dvidx <- ((good - 1) %/% NB) + 1; ntidx <- ((good - 1) %% NB) + 1
    peak_of <- function(x) { b <- good[which.max(smap(x))]; sprintf("DV%d-NT%d", ((b - 1) %/% NB) + 1, ((b - 1) %% NB) + 1) }
    lig <- if (ligand %in% genes) smap(L[ligand, ]) else rep(NA, length(good))

    write_wide <- function(rows, labels, file) {
        M <- rbind(dvidx, ntidx, rows); dimnames(M) <- NULL       # columns -> V1..Vn
        df <- data.frame(rowlabel = labels, round(as.data.frame(M), 4), check.names = FALSE)
        write_csv(df, file.path(OUTDIR, file))
    }
    fgfra <- rbind(lig, mcol("FGFfeedback"), mcol("RAdegradation"), mcol("RAsynthesis"))
    write_wide(fgfra, c("dvidx", "ntidx", "FGFligand", "FGFfeedback", "RAdegradation", "RAsynthesis"),
               sprintf("figdata_%s.csv", species))
    bmp <- rbind(mcol("BMPligand"), mcol("BMPtarget"), mcol("BMPantagonist"))
    write_wide(bmp, c("dvidx", "ntidx", "BMPligand", "BMPtarget", "BMPantagonist"),
               sprintf("figdata_bmp_%s.csv", species))
    cat(sprintf("\n[%s] %d cells, %d bins kept (>= %d cells), FGF ligand = %s\n",
                species, ncol(obj), length(good), MINCELL, ligand))

    sel <- gate_cells(dv, nt, gates[[species]])
    cat(sprintf("=== %s: gate %d cells / rest %d cells ===\n", species, sum(sel), sum(!sel)))
    cat("--- FGF ligands at the gate ---\n");   print(gate_report(L, fgf_ligands, sel, peak_of), row.names = FALSE)
    cat("--- BMP genes at the gate ---\n");     print(gate_report(L, bmp_report, sel, peak_of), row.names = FALSE)
    invisible(NULL)
}

# %% [markdown]
# ## Chick (deposited RPC Seurat object, n = 29,025)

# %%
chick <- readRDS(file.path(here::here(), "data", "20250604_02_fabp7.rds"))
if ("RNA" %in% SeuratObject::Assays(chick)) SeuratObject::DefaultAssay(chick) <- "RNA"
chick <- tryCatch(JoinLayers(chick), error = function(e) chick)
module_maps(chick, "chick", ligand = "FGF8")
rm(chick); invisible(gc())

# %% [markdown]
# ## Human (R-export-stage RPCs; published DV/NT scores)

# %%
# R-export-stage human RPCs (23,031 cells) -- the object the published Fig. 8D used.
human <- read_human_rexport()
module_maps(human, "human", ligand = "FGF7")
cat("\nwrote", OUTDIR, "\n")
