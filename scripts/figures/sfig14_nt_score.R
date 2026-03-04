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
# # Figure S14: Generation of NT Score from Chicken Retina scRNA-seq

# %%
source("../preprocessing/00_utils.R")

# %% [markdown]
# ## Output directory setup

# %%
FIGURES_BASE <- file.path(dirname(sys.frame(1)$ofile %||% "."), "..", "..", "figures")
dir.create(file.path(FIGURES_BASE, "Figure_SF14"), recursive = TRUE, showWarnings = FALSE)

# %% [markdown]
# ## SF14B: Nasal vs Temporal score scatter

# %% tags=["cell-201"]
options( repr.plot.width = 7, repr.plot.height = 7 )
p10 <-
FetchData(
    fabp7,
    vars = c("Nasal.Score1", "Temporal.Score1"),
    slot = "data"
) %>%
ggplot( aes( x = Nasal.Score1, y = Temporal.Score1 ) ) +
geom_point( alpha = 0.1 ) +
geom_density_2d() +
scale_x_continuous( expand = expansion( add = 0, mult = 0) ) +
scale_y_continuous( expand = expansion( add = 0, mult = 0) ) +
theme(
    text = element_text( size = 20 )
)
options( repr.plot.width = 7, repr.plot.height = 7 )
p10
ggsave(p10, filename=file.path(FIGURES_BASE, "Figure_SF14", "SF14B_nt_scatter.png"),  width = 7, height = 7 )

# %% [markdown]
# ## SF14A: NT score distribution

# %% tags=["cell-204"]
p10.1 <-
FetchData(
  fabp7,
  vars = c("NT.Score", "Nasal.Score1", "Temporal.Score1"),
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
p10.1
ggsave(p10.1, filename=file.path(FIGURES_BASE, "Figure_SF14", "SF14A_nt_distribution.png"),  width = 10, height = 5 )

# %% [markdown]
# ## SF14C: Gene expression along NT axis (greyed anchor genes)

# %% tags=["cell-206"]
genes.nasal <- 
            c(
                "FOXG1", "SOHO-1", "HMX1", "EFNA5", "EFNA2" # listed by Heer
            )
genes.temporal <- c(
              "FOXD1", "EPHA3"
)
p11 <-
FetchData(
    fabp7,
    vars = c(
        "NT.Score", 
        genes.nasal,
        genes.temporal,
        "FGF8", "CYP26A1", "CYP26C1", "BMP2", "ID1", "CYP1B1"
    ),
    layer = "data"
) %>%
pivot_longer( cols = !any_of("NT.Score"), names_to = "gene", values_to = "expression" ) %>%
mutate(
    gene = factor(
        gene,
        levels = c(
            genes.nasal,
            genes.temporal,
            "FGF8", "CYP26A1", "CYP26C1", "BMP2", "ID1", "CYP1B1"
        )
    ),
    gene.group = case_when(
        gene %in% genes.nasal ~ "Nasal",
        gene %in% genes.temporal ~ "Temporal",
        gene %in% c("CYP1B1") ~ "Stripe",
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
    ggplot( ., aes( x = NT.Score, y = expression, colour = gene.group ) ) +
    # geom_point( alpha = 0.1 ) +
    geom_smooth( linewidth = 2, linetype = "dashed", se = FALSE ) +
    expand_limits( y = 0 ) +
    scale_y_continuous( expand = expansion( add = 0, mult = c(0, 0.05) ) ) +
    theme(
        text = element_text( size = 20 )
    ) +
    expand_limits( y = 0 ) +
    ggh4x::facet_wrap2( 
        gene ~ ., 
        scale = "free_y",
        strip = ggh4x::strip_themed( background_x = ggh4x::elem_list_rect( fill = .$strip.color ) )
    ) 
}
options( repr.plot.width = 28, repr.plot.height = 14 )
p11
ggsave(p11, filename=file.path(FIGURES_BASE, "Figure_SF14", "SF14C_nt_expression.png"),  width = 28, height = 14 )

# %% tags=["cell-209"]
p12 <-
FetchData(
    fabp7,
    vars = c(
        "NT.Score", 
        "CYP1B1", "FOXG1", "FOXD1", "FGF8", "SOHO-1", "HMX1", "CYP26C1", "CYP26A1"
    ),
    slot = "data"
) %>%
pivot_longer( cols = !any_of("NT.Score"), names_to = "gene", values_to = "expression" ) %>%
mutate(
    gene = factor(
        gene,
        levels = c(
            genes.nasal,
            genes.temporal,
            "FGF8", "CYP26A1", "CYP26C1", "BMP2", "ID1", "CYP1B1"
        )
    ),
    gene.group = case_when(
        gene %in% genes.nasal ~ "Nasal",
        gene %in% genes.temporal ~ "Temporal",
        gene %in% c("CYP1B1") ~ "Stripe",
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
# dplyr::filter( expression > 0 ) %>%
ggplot( aes( x = NT.Score, y = expression, colour = gene ) ) +
# geom_point( alpha = 0.1 ) +
geom_smooth( linewidth = 2, linetype = "dashed", se = FALSE ) +
expand_limits( y = 0 ) +
scale_y_continuous( expand = expansion( add = 0, mult = c(0, 0.05) ) ) +
theme(
    text = element_text( size = 20 )
)
options( repr.plot.width = 14, repr.plot.height = 7 )
p12
ggsave(p12, filename=file.path(FIGURES_BASE, "Figure_SF14", "SF14C_nt_selected_genes.png"),  width = 14, height = 7 )

# %% [markdown]
# ## SF14D: Binned NT expression with embryo points

# %% [markdown]
# ### Figure S4D3

# %% tags=["cell-223"]
# Set plot dimensions
options(repr.plot.width = 28, repr.plot.height = 14)
require(ggforce)
# Select genes to plot (you can modify this list)
select.genes <- c("CYP1B1", "CYP26C1", "FGF8", "FOXD1", "FOXG1", "HMX1", "SOHO-1")
require(ggforce)
require(ggh4x)
p13 <-
FetchData(
  fabp7,
  vars = c("NT.Score", "library", "genotype", select.genes),
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
        gene %in% genes.nasal ~ "Nasal",
        gene %in% genes.temporal ~ "Temporal",
        gene %in% c("CYP1B1") ~ "Stripe",
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
    NT_bin = {
        min_score <- min(NT.Score, na.rm = TRUE)
        max_score <- max(NT.Score, na.rm = TRUE)
        cut(
            NT.Score,
            breaks = seq(min_score, max_score, length.out = 21),  # Reduced to 20 bins for readability
            labels = 1:20  # Using integer labels for bins
        )
    }
) %>%
group_by(gene, library, genotype, NT_bin, gene.group, strip.color) %>%
dplyr::filter(n() > 10) %>%
summarize(
    mean_expression = mean(expression), 
    mean_NT.Score = mean(NT.Score),
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
            x = NT_bin,  # Use DV_bin for x-axis
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
        x = "NT Score Bins",
        y = "Mean Expression",
        title = "Gene Expression Distribution Along NT Axis",
        subtitle = "Aggregated by genotype and NT Score Bins (at least 10 cells)"
      )
}
# Set plot dimensions
options(repr.plot.width = 28, repr.plot.height = 14)
p13
ggsave(p13, filename=file.path(FIGURES_BASE, "Figure_SF14", "SF14D_nt_binned.png"),  width = 28, height = 14 )
