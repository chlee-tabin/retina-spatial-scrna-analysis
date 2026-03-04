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
# # Figure S12: Processing of scRNA-seq Datasets from Developing Chicken Retina

# %%
source("../preprocessing/00_utils.R")

# %% [markdown]
# ## Output directory setup

# %%
FIGURES_BASE <- file.path(dirname(sys.frame(1)$ofile %||% "."), "..", "..", "figures")
dir.create(file.path(FIGURES_BASE, "Figure_SF12"), recursive = TRUE, showWarnings = FALSE)

# %% [markdown]
# ## SF12B: Violin plot of marker genes

# %% tags=["cell-125"]
p2 <-
FetchData(
    retina,
    vars = c("annotation", "FABP7", "OTX2", "NEUROD1", "THRB", "RXRG", "TFAP2A", "ISL1", "GYPC", "IFI6")
) %>%
pivot_longer( cols = !any_of("annotation"), names_to = "gene", values_to = "expression" ) %>%
mutate(
    annotation2 = factor( annotation, levels = c("RPC", "OTX2+ neurogenic", "THRB+RXRG+ cone", "Horizontal/Amacrine", "Retinal ganglion", "Vascular", "Infected") ),
    gene2 = factor( gene, levels = c("FABP7", "OTX2", "NEUROD1", "THRB", "RXRG", "TFAP2A", "ISL1", "GYPC", "IFI6") )
) %>%
ggplot( aes( x = annotation2, y = expression ) ) +
geom_violin( scale = "width" ) +
facet_wrap( gene2 ~ . ) +
scale_y_continuous( expand = expansion( add = 0, mult = c(0, 0.05) ) ) +
labs(
    x = "Graph-based clusters (leiden=0.1)"
) +
theme(
    text = element_text( size = 20 ),
    axis.text.x = element_text( angle = 45, hjust = 1 )
)
options( repr.plot.width = 21, repr.plot.height = 14 )
p2
ggsave(p2, filename=file.path(FIGURES_BASE, "Figure_SF12", "SF12B_marker_violins.png"),  width = 21, height = 14 )

# %% [markdown]
# ## SF12A: UMAP of cell clusters

# %% tags=["cell-126"]
Idents( retina ) <- "annotation"
p3 <-
DimPlot(
    retina,
    group.by = c("annotation"),
    reduction = "umap.harmony",
    label = T
) +
theme(
    axis.text = element_blank(),
    axis.ticks = element_blank()
)
options( repr.plot.width = 9, repr.plot.height = 7 )
p3
ggsave( p3, filename=file.path(FIGURES_BASE, "Figure_SF12", "SF12A_umap_clusters.png"),  width = 9, height = 7 )

# %% [markdown]
# ## SF12C: UMAP QC panels (phase, tech, mito, gene count)

# %% tags=["cell-127"]
Idents( retina ) <- "annotation"
# p4 <-
# DimPlot(
#     retina,
#     group.by = c("Phase", "technology"),
#     reduction = "umap.harmony",
#     label = F
# ) &
# theme(
#     axis.text = element_blank(),
#     axis.ticks = element_blank()
# )
p4.1 <-
FetchData(
    retina,
    vars = c("Phase", "technology", "annotation", "umapharmony_1", "umapharmony_2")
) %>%
dplyr::slice_sample(prop = 1) %>% {
    ggplot( ., aes( x = umapharmony_1, y = umapharmony_2, colour = Phase ) ) +
    geom_point(  alpha = 0.5 ) +
    geom_text( 
        data = . %>%
            dplyr::group_by( annotation ) %>%
            summarize(
                umapharmony_1 = mean( umapharmony_1 ),
                umapharmony_2 = mean( umapharmony_2 )
            ),
        aes( label = annotation ),
        colour = "black"
    ) +
    guides( colour = guide_legend( override.aes = list( alpha = 1 ) ) ) +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.line.x.top = element_blank(),
        axis.line.y.right = element_blank(),
        axis.line = element_line( color = "black" ),
        panel.border = element_blank()
    )
}
p4.2 <-
FetchData(
    retina,
    vars = c("Phase", "technology", "annotation", "umapharmony_1", "umapharmony_2")
) %>%
dplyr::slice_sample(prop = 1) %>% {
    ggplot( ., aes( x = umapharmony_1, y = umapharmony_2, colour = technology ) ) +
    geom_point(  alpha = 0.5 ) +
    geom_text( 
        data = . %>%
            dplyr::group_by( annotation ) %>%
            summarize(
                umapharmony_1 = mean( umapharmony_1 ),
                umapharmony_2 = mean( umapharmony_2 )
            ),
        aes( label = annotation ),
        colour = "black"
    ) +
    guides( colour = guide_legend( override.aes = list( alpha = 1 ) ) ) +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.line.x.top = element_blank(),
        axis.line.y.right = element_blank(),
        axis.line = element_line( color = "black" ),
        panel.border = element_blank()
    )
}
p4.3 <-
FetchData(
    retina,
    vars = c("Phase", "technology", "annotation", "percent.mito", "umapharmony_1", "umapharmony_2")
) %>%
dplyr::arrange( percent.mito ) %>% 
{
    ggplot( ., aes( x = umapharmony_1, y = umapharmony_2, colour = percent.mito ) ) +
    geom_point( alpha = 0.5 ) +
    geom_label( 
        data = . %>%
            dplyr::group_by( annotation ) %>%
            summarize(
                umapharmony_1 = mean( umapharmony_1 ),
                umapharmony_2 = mean( umapharmony_2 )
            ),
        aes( label = annotation ),
        colour = "black"
    ) +
    scale_colour_viridis_c() +
#    guides( colour = guide_legend( override.aes = list( alpha = 1 ) ) ) +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.line.x.top = element_blank(),
        axis.line.y.right = element_blank(),
        axis.line = element_line( color = "black" ),
        panel.border = element_blank()
    )
}
p4.4 <-
FetchData(
    retina,
    vars = c("Phase", "technology", "annotation", "nFeature_RNA", "umapharmony_1", "umapharmony_2")
) %>%
dplyr::arrange( nFeature_RNA ) %>% 
{
    ggplot( ., aes( x = umapharmony_1, y = umapharmony_2, colour = nFeature_RNA ) ) +
    geom_point( alpha = 0.5 ) +
    geom_label( 
        data = . %>%
            dplyr::group_by( annotation ) %>%
            summarize(
                umapharmony_1 = mean( umapharmony_1 ),
                umapharmony_2 = mean( umapharmony_2 )
            ),
        aes( label = annotation ),
        colour = "black"
    ) +
    scale_colour_viridis_c() +
#    guides( colour = guide_legend( override.aes = list( alpha = 1 ) ) ) +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.line.x.top = element_blank(),
        axis.line.y.right = element_blank(),
        axis.line = element_line( color = "black" ),
        panel.border = element_blank()
    )
}
options( repr.plot.width = 14, repr.plot.height = 14 )
p4 <-
(p4.1 + p4.2) / (p4.3 + p4.4)
p4
# options( repr.plot.width = 14, repr.plot.height = 7 )
# p4
ggsave( p4, filename=file.path(FIGURES_BASE, "Figure_SF12", "SF12C_umap_qc.png"),  width = 14, height = 7 )

# %% [markdown]
# ## SF12E: Sex chromosome fractions

# %% tags=["cell-128"]
p5 <-
retina@meta.data %>%
dplyr::filter( scDblFinder.class == "singlet", grepl( "(donor|singlet)", genotype ) ) %>%
dplyr::select( library, barcode, genotype, percent.W, percent.Z ) %>%
group_by( library, genotype ) %>%
dplyr::filter( n() > 100 ) %>%
mutate(
    genotype2 = paste0( genotype, "(", n(), ")" )
) %>%
pivot_longer( cols = starts_with( "percent" ), names_to = "metric", values_to = "value" ) %>%
ggplot( aes( x = genotype2, y = value ) ) +
geom_violin( scale = "width" ) +
scale_y_continuous( expand = expansion( add = 0, mult = c(0, 0.05) ) ) +
facet_grid( metric ~ library, scale = "free" ) +
theme(
    text = element_text( size = 20 ),
    axis.text.x = element_text( angle = 45, hjust = 1 )
) +
labs(
    x = "Inferred genotype",
    y = "fraction per cell"
)
options( repr.plot.width = 28, repr.plot.height = 14 )
p5
ggsave( p5, filename=file.path(FIGURES_BASE, "Figure_SF12", "SF12E_sex_chromosomes.png"),  width = 28, height = 14 )

# %% [markdown]
# ## SF12D: Cell type composition per embryo

# %% tags=["cell-129"]
p5.2 <-
retina@meta.data %>%
dplyr::filter( scDblFinder.class == "singlet", grepl( "(donor|singlet)", genotype ) ) %>%
dplyr::select( library, barcode, genotype, annotation ) %>%
group_by( library, genotype ) %>%
dplyr::filter( n() > 100 ) %>%
mutate(
    genotype2 = paste0( library, ".", genotype, "(", n(), ")" )
) %>%
dplyr::count( genotype2, annotation ) %>%
mutate(
    fraction = n / sum(n)
) %>%
ggplot( aes( x = genotype2, y = fraction, fill = annotation ) ) +
geom_bar( stat = "identity" ) +
scale_y_continuous( expand = expansion( add = 0, mult = c(0, 0) ) ) +
theme(
    text = element_text( size = 20 ),
    axis.text.x = element_text( angle = 90, hjust = 1 )
) +
labs(
    x = "Inferred genotype",
    y = "Cell type annotation composition"
)
options( repr.plot.width = 28, repr.plot.height = 14 )
p5.2
ggsave( p5.2, filename=file.path(FIGURES_BASE, "Figure_SF12", "SF12D_cell_composition.png"),  width = 28, height = 14 )
