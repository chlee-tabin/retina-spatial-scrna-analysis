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
# # Figure 8D: Pathway-module topography of the high-acuity area, chick vs. human
#
# Renders the six-column module map (FGF ligand | FGF feedback | BMP ligand |
# BMP antagonist | RA degradation | RA synthesis) for chick (top; HAA =
# central-nasal) and human (bottom; fovea = temporal-of-center). Each panel is
# z-scored across its own bins (module scores and single genes are not on a
# common absolute scale, so the spatial pattern is what is compared); colors
# are clipped at +/-2 SD. The green rectangle marks the HAA / fovea gate.
#
# Input: the per-bin CSVs in `data/fig8d_pathway_modules/`, committed so this
# renders with tidyverse alone. They are produced from the GEO objects by
# `fig8d_pathway_module_data.R` (same 15 x 15 grid, bins >= 12 cells).
# Titleless by convention; the caption is the manuscript legend.

# %%
suppressWarnings(suppressMessages(library(tidyverse)))

DD     <- file.path(here::here(), "data", "fig8d_pathway_modules")
OUTDIR <- file.path(here::here(), "figures", "Figure8")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# Read the requested module rows from a wide CSV (rowlabel + V1..Vn bin columns)
# and z-score each panel across its own bins.
load_items <- function(path, sp, items) {
    m <- suppressMessages(read_csv(path, show_col_types = FALSE))
    lab <- m$rowlabel; M <- as.matrix(m[, -1])
    dvidx <- as.integer(M[lab == "dvidx", ]); ntidx <- as.integer(M[lab == "ntidx", ])
    purrr::map_dfr(items, function(it) {
        v <- as.numeric(M[lab == it, ]); s <- sd(v, na.rm = TRUE)
        tibble(species = sp, item = it, ntidx = ntidx, dvidx = dvidx,
               z = if (is.finite(s) && s > 0) (v - mean(v, na.rm = TRUE)) / s else v * 0)
    })
}

fgfra <- c("FGFligand", "FGFfeedback", "RAdegradation", "RAsynthesis")
bmp   <- c("BMPligand", "BMPantagonist")
df <- bind_rows(
    load_items(file.path(DD, "figdata_chick.csv"),     "chick", fgfra),
    load_items(file.path(DD, "figdata_human.csv"),     "human", fgfra),
    load_items(file.path(DD, "figdata_bmp_chick.csv"), "chick", bmp),
    load_items(file.path(DD, "figdata_bmp_human.csv"), "human", bmp)
)

# %%
item_levels <- c("FGFligand", "FGFfeedback", "BMPligand", "BMPantagonist", "RAdegradation", "RAsynthesis")
item_labs <- c(
    "FGF ligand\n(chick: FGF8 / human: FGF7)",
    "FGF-feedback\n(SPRY/DUSP/ETV5)",
    "BMP ligand\n(BMP2/4/7)",
    "BMP antagonist\n(SOSTDC1/NOG/GREM/FST)",
    "RA-degradation\n(CYP26)",
    "RA-synthesis\n(ALDH1A)"
)
df$item <- factor(df$item, levels = item_levels, labels = item_labs)

sp_labs <- c("chick  (HAA = central-nasal)", "human  (fovea = temporal-center)")
df$species <- factor(df$species, levels = c("chick", "human"), labels = sp_labs)

# HAA / fovea gate boxes: the (NT, DV) quadrants of the 10 x 10 gate grid in
# fig8d_pathway_module_data.R (keep in sync), mapped onto the 15 x 15 display grid
# (bin i spans fractions (i-1)/15..i/15, so fraction f sits at f * 15 + 0.5).
gates <- list(chick = list(c(7,5), c(6,5), c(7,4), c(6,4)),
              human = list(c(3,5), c(4,5), c(3,4), c(4,4)))
gate_box <- function(q, n = 10, nb = 15) {
    nt <- sapply(q, `[`, 1); dv <- sapply(q, `[`, 2)
    c(xmin = min(nt), xmax = max(nt) + 1, ymin = min(dv), ymax = max(dv) + 1) / n * nb + 0.5
}
boxes <- bind_cols(species = factor(sp_labs, levels = levels(df$species)),
                   as_tibble(do.call(rbind, lapply(gates, gate_box))))
# Species-specific FGF-ligand label, top-left inside the FGF-ligand column.
ligtext <- tibble(
    species = factor(sp_labs, levels = levels(df$species)),
    item = factor(item_labs[1], levels = levels(df$item)),
    ntidx = c(2.7, 3.7), dvidx = 13, lab = c("FGF8", "FGF7")
)

p <- ggplot(df, aes(ntidx, dvidx, fill = z)) +
    geom_tile() +
    facet_grid(species ~ item) +
    geom_rect(data = boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE, fill = NA, color = "green3", linewidth = 0.9) +
    geom_text(data = ligtext, aes(ntidx, dvidx, label = lab),
              inherit.aes = FALSE, fontface = "bold", size = 3, hjust = 0, vjust = 1) +
    scale_fill_gradient2(low = "#2166ac", mid = "grey95", high = "#b2182b", midpoint = 0,
                         limits = c(-2, 2), oob = scales::squish, name = "z\n(per panel)") +
    coord_fixed(expand = FALSE) +
    labs(x = "NT:  temporal → nasal", y = "DV:  ventral → dorsal") +
    theme_minimal(base_size = 11) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(), strip.text = element_text(face = "bold"))

ggsave(file.path(OUTDIR, "F8D_pathway_module_topography.png"), p, width = 19, height = 7.2, dpi = 200)
ggsave(file.path(OUTDIR, "F8D_pathway_module_topography.pdf"), p, width = 19, height = 7.2, device = cairo_pdf)
cat("wrote", file.path(OUTDIR, "F8D_pathway_module_topography.{png,pdf}"), "\n")
