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
# # Figure S13: Generation of DV Score from Chicken Retina scRNA-seq

# %%
source("../preprocessing/00_utils.R")

# %% [markdown]
# ## Output directory setup

# %%
FIGURES_BASE <- file.path(dirname(sys.frame(1)$ofile %||% "."), "..", "..", "figures")
dir.create(file.path(FIGURES_BASE, "Figure_SF13"), recursive = TRUE, showWarnings = FALSE)

# %% [markdown]
# ## SF13B: Dorsal vs Ventral score scatter

# %% tags=["cell-167"]
p8 <-
FetchData(
    fabp7,
    vars = c("Ventral.Score1", "Dorsal.Score1"),
    layer = "data"
) %>%
ggplot( aes( x = Ventral.Score1, y = Dorsal.Score1 ) ) +
geom_point( alpha = 0.1 ) +
geom_density_2d() +
scale_x_continuous( expand = expansion( add = 0, mult = 0) ) +
scale_y_continuous( expand = expansion( add = 0, mult = 0) ) +
theme(
    text = element_text( size = 20 )
)
options( repr.plot.width = 7, repr.plot.height = 7 )
p8
ggsave(p8, filename=file.path(FIGURES_BASE, "Figure_SF13", "SF13B_dv_scatter.png"),  width = 7, height = 7 )

# %% [markdown]
# ## SF13A: DV score distribution

# %% tags=["cell-169"]
p8.1 <-
FetchData(
  fabp7,
  vars = c("DV.Score", "Dorsal.Score1", "Ventral.Score1"),
  slot = "data"
) %>%
pivot_longer( cols = everything(), names_to = "score.category", values_to = "score" ) %>%
ggplot(
    aes(
        x = score
    )
) +
geom_histogram() +
scale_y_continuous( expand = expansion( add = 0, mult = c(0, 0.05) ) ) +
facet_wrap( score.category ~ . ) +
theme(
    text = element_text( size = 20 )
)
options( repr.plot.width = 10, repr.plot.height = 5 )
p8.1
ggsave(p8.1, filename=file.path(FIGURES_BASE, "Figure_SF13", "SF13A_dv_distribution.png"),  width = 10, height = 5 )

# %% [markdown]
# ## SF13C: Gene expression along DV axis (greyed anchor genes)

# %% tags=["cell-171"]
genes.dorsal <- c(
        "TBX5",  "TBX3",   "TBX2", "ALDH1A1", "EFNB2", "EFNB1", "GDF6", "NR2F2" ,
        "BAMBI", "CRABP1", "ID2",  "ID3",     "CCND2", "ENOX1", "ANOS1", "UNC5B", "LYPD6"
)
genes.ventral <- c(
        "VAX1", "CHRDL1", "ALDH1A3", "EPHB2", "EPHB3", "PAX2" ,
        "SMOC1", "BMPR1B", "SHISA2", "RDH10", "FRAS1", "SORCS2", "IRX3"
)
p6 <-
FetchData(
    fabp7,
    vars = c(
        "DV.Score", 
        "TBX5", "TBX3", "TBX2", "FGF8", "VAX1", "CYP26A1"
        # "PTPRU", 
        # "VSNL1"
    ),
    slot = "data"
) %>%
pivot_longer( cols = !any_of("DV.Score"), names_to = "gene", values_to = "expression" ) %>%
mutate(
    gene.group = case_when(
        gene %in% genes.dorsal ~ "Dorsal",
        gene %in% genes.ventral ~ "Ventral",
        gene %in% c("FGF8", "CYP26A1", "CYP26C1", "BMP2") ~ "Stripe",
        TRUE ~ "other"
    )
) %>%
# dplyr::filter( expression > 0 ) %>%
ggplot( aes( x = DV.Score, y = expression, colour = gene ) ) +
geom_smooth( 
    linewidth = 2, 
    linetype = "dashed", 
    se = FALSE, 
    method = "loess",
    span = 0.2  # Smaller span = less smoothing (default is 0.75); 20% of the closest points are used for regression
) +
expand_limits( y = 0 ) +
scale_y_continuous( expand = expansion( add = 0, mult = c(0, 0.05) ) ) +
theme(
    text = element_text( size = 20 )
)
options( repr.plot.width = 14, repr.plot.height = 7 )
p6
ggsave(p6, filename=file.path(FIGURES_BASE, "Figure_SF13", "SF13C_dv_expression.png"),  width = 14, height = 7 )

# %% [markdown]
# ## SF13D: Binned DV expression with embryo points

# %% [markdown]
# ### Figure S4C2

# %% tags=["cell-196"]
genes.dorsal <- c(
        "TBX5",  "TBX3",   "TBX2", "ALDH1A1", "EFNB2", "EFNB1", "GDF6", "NR2F2" ,
        "BAMBI", "CRABP1", "ID2",  "ID3",     "CCND2", "ENOX1", "ANOS1", "UNC5B", "LYPD6"
)
genes.ventral <- c(
        "VAX1", "CHRDL1", "ALDH1A3", "EPHB2", "EPHB3", "PAX2" ,
        "SMOC1", "BMPR1B", "SHISA2", "RDH10", "FRAS1", "SORCS2", "IRX3"
)
p7 <-
FetchData(
    fabp7,
    vars = c(
        "DV.Score", 
        genes.dorsal,
        genes.ventral,
        "FGF8", "CYP26A1", "CYP26C1", "BMP2", "ID1", "PTPRU"
    ),
    slot = "data"
) %>%
pivot_longer( cols = !any_of("DV.Score"), names_to = "gene", values_to = "expression" ) %>%
mutate(
    gene = factor(
        gene,
        levels = c(
            genes.dorsal,
            genes.ventral,
            "FGF8", "CYP26A1", "CYP26C1", "BMP2", "ID1", "PTPRU"
        )
    ),
    gene.group = case_when(
        gene %in% genes.dorsal ~ "Dorsal",
        gene %in% genes.ventral ~ "Ventral",
        gene %in% c("FGF8", "CYP26A1", "CYP26C1", "BMP2") ~ "Stripe",
        TRUE ~ "other"
    ),
    strip.color = case_when(
        gene %in% c("VAX1", "CHRDL1", "ALDH1A3") ~ "grey", # ventral
        gene %in% c("TBX5", "TBX2", "TBX3", "ALDH1A1", "EFNB2", "EFNB1") ~ "grey", #dorsal
        gene %in% c("FOXG1", "SOHO-1", "HMX1", "EFNA5", "EFNA2") ~ "grey", #nasal
        gene %in% c("FOXD1", "EPHA3") ~ "grey", # temporal
        TRUE ~ "white"
    )
) %>% {
# dplyr::filter( expression > 0 ) %>%
    ggplot( ., aes( x = DV.Score, y = expression, colour = gene.group ) ) +
    # geom_point( alpha = 0.1 ) +
    expand_limits( y = 0 ) +
    geom_smooth( linewidth = 2, linetype = "dashed", se = FALSE ) +
    scale_y_continuous( expand = expansion( add = 0, mult = c(0, 0.05) ) ) +
    theme(
        text = element_text( size = 20 )
    ) +
    ggh4x::facet_wrap2( 
        gene ~ ., 
        scale = "free_y",
        strip = ggh4x::strip_themed( background_x = ggh4x::elem_list_rect( fill = .$strip.color ) )
    ) 
}
options( repr.plot.width = 28, repr.plot.height = 14 )
p7
ggsave(p7, filename=file.path(FIGURES_BASE, "Figure_SF13", "SF13D_dv_all_genes.png"),  width = 28, height = 14 )

# %% tags=["cell-197"]
require(ggforce)
require(ggh4x)
# Select genes to plot (you can modify this list)
select.genes <- c("TBX5", "TBX3", "TBX2", "FGF8", "VAX1", "CHRDL1", "CYP26A1", "CYP26C1")
p9 <-
FetchData(
  fabp7,
  vars = c("DV.Score", "library", "genotype", select.genes),
  slot = "data"
) %>%
pivot_longer(cols = all_of(select.genes), names_to = "gene", values_to = "expression") %>%
mutate(
    gene = factor(
        gene,
        levels = c(
            select.genes
        )
    ),
    gene.group = case_when(
      gene %in% genes.dorsal ~ "Dorsal",
      gene %in% genes.ventral ~ "Ventral",
      gene %in% c("FGF8", "CYP26A1", "CYP26C1", "BMP2") ~ "Stripe",
      TRUE ~ "other"
    ),
    strip.color = case_when(
        gene %in% c("VAX1", "CHRDL1", "ALDH1A3") ~ "grey", # ventral
        gene %in% c("TBX5", "TBX2", "TBX3", "ALDH1A1", "EFNB2", "EFNB1") ~ "grey", #dorsal
        gene %in% c("FOXG1", "SOHO-1", "HMX1", "EFNA5", "EFNA2") ~ "grey", #nasal
        gene %in% c("FOXD1", "EPHA3") ~ "grey", # temporal
        TRUE ~ "white"
    )
) %>%
mutate(
    DV_bin = {
        min_score <- min(DV.Score, na.rm = TRUE)
        max_score <- max(DV.Score, na.rm = TRUE)
        cut(
            DV.Score,
            breaks = seq(min_score, max_score, length.out = 21),  # Reduced to 20 bins for readability
            labels = 1:20  # Using integer labels for bins
        )
    }
) %>%
group_by(gene, library, genotype, DV_bin, gene.group, strip.color) %>%
dplyr::filter(n() > 10) %>%
summarize(
    mean_expression = mean(expression), 
    mean_DV.Score = mean(DV.Score),
    strip.color = first(strip.color),  # CRITICAL: Preserve this column!
    .groups = "drop" 
) %>% {
    # Extract unique colors in correct factor level order within the pipe
    plot_data <- .
    unique_colors <- plot_data %>% 
        distinct(gene, strip.color) %>% 
        arrange(factor(gene, levels = select.genes)) %>%  # Order by factor levels
        pull(strip.color)
    ggplot(
        plot_data,
        aes(
            x = DV_bin,  # Use DV_bin for x-axis
            y = mean_expression, 
            fill = gene
        )
    ) +
      geom_boxplot(outlier.shape = NA) +
      ggbeeswarm::geom_beeswarm() +
      ggh4x::facet_wrap2( 
          gene ~ ., 
          scale = "free_y",
          ncol = 3,
          strip = ggh4x::strip_themed( 
              background_x = ggh4x::elem_list_rect( fill = unique_colors )
          )
      ) +
      theme_bw() +
      theme(
        text = element_text(size = 30),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      ) +
      scale_y_continuous(
          expand = expansion( add = 0, mult = c(0, 0.05) )
      ) +
      labs(
        x = "DV Score Bins",
        y = "Mean Expression",
        title = "Gene Expression Distribution Along DV Axis",
        subtitle = "Aggregated by genotype and DV Score Bins (at least 10 cells)"
      )
}
# Set plot dimensions
options(repr.plot.width = 28, repr.plot.height = 14)
p9
ggsave(p9, filename=file.path(FIGURES_BASE, "Figure_SF13", "SF13D_dv_binned_8genes.png"),  width = 28, height = 14 )

# %% tags=["cell-198"]
require(ggforce)
require(ggh4x)
# Select genes to plot (you can modify this list)
select.genes <- c("TBX5", "TBX3", "TBX2", "FGF8", "VAX1", "CHRDL1", "CYP26A1", "CYP26C1", "BMP2")
p9 <-
FetchData(
  fabp7,
  vars = c("DV.Score", "library", "genotype", select.genes),
  slot = "data"
) %>%
pivot_longer(cols = all_of(select.genes), names_to = "gene", values_to = "expression") %>%
mutate(
    gene = factor(
        gene,
        levels = c(
            select.genes
        )
    ),
    gene.group = case_when(
      gene %in% genes.dorsal ~ "Dorsal",
      gene %in% genes.ventral ~ "Ventral",
      gene %in% c("FGF8", "CYP26A1", "CYP26C1", "BMP2") ~ "Stripe",
      TRUE ~ "other"
    ),
    strip.color = case_when(
        gene %in% c("VAX1", "CHRDL1", "ALDH1A3") ~ "grey", # ventral
        gene %in% c("TBX5", "TBX2", "TBX3", "ALDH1A1", "EFNB2", "EFNB1") ~ "grey", #dorsal
        gene %in% c("FOXG1", "SOHO-1", "HMX1", "EFNA5", "EFNA2") ~ "grey", #nasal
        gene %in% c("FOXD1", "EPHA3") ~ "grey", # temporal
        TRUE ~ "white"
    )
) %>%
mutate(
    DV_bin = {
        min_score <- min(DV.Score, na.rm = TRUE)
        max_score <- max(DV.Score, na.rm = TRUE)
        cut(
            DV.Score,
            breaks = seq(min_score, max_score, length.out = 21),  # Reduced to 20 bins for readability
            labels = 1:20  # Using integer labels for bins
        )
    }
) %>%
group_by(gene, library, genotype, DV_bin, gene.group, strip.color) %>%
dplyr::filter(n() > 10) %>%
summarize(
    mean_expression = mean(expression), 
    mean_DV.Score = mean(DV.Score),
    strip.color = first(strip.color),  # CRITICAL: Preserve this column!
    .groups = "drop" 
) %>% {
    # Extract unique colors in correct factor level order within the pipe
    plot_data <- .
    unique_colors <- plot_data %>% 
        distinct(gene, strip.color) %>% 
        arrange(factor(gene, levels = select.genes)) %>%  # Order by factor levels
        pull(strip.color)
    ggplot(
        plot_data,
        aes(
            x = DV_bin,  # Use DV_bin for x-axis
            y = mean_expression, 
            fill = gene
        )
    ) +
      geom_boxplot(outlier.shape = NA) +
      ggbeeswarm::geom_beeswarm() +
      ggh4x::facet_wrap2( 
          gene ~ ., 
          scale = "free_y",
          ncol = 3,
          strip = ggh4x::strip_themed( 
              background_x = ggh4x::elem_list_rect( fill = unique_colors )
          )
      ) +
      theme_bw() +
      theme(
        text = element_text(size = 30),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      ) +
      scale_y_continuous(
          expand = expansion( add = 0, mult = c(0, 0.05) )
      ) +
      labs(
        x = "DV Score Bins",
        y = "Mean Expression",
        title = "Gene Expression Distribution Along DV Axis",
        subtitle = "Aggregated by genotype and DV Score Bins (at least 10 cells)"
      )
}
# Set plot dimensions
options(repr.plot.width = 28, repr.plot.height = 14)
p9
ggsave(p9, filename=file.path(FIGURES_BASE, "Figure_SF13", "SF13D_dv_binned_9genes.png"),  width = 28, height = 14 )
