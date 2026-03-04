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
# # Chick Retina Preprocessing: Load, QC, Integration, Clustering

# %%
source("00_utils.R")

# %% [markdown]
# ## Data directory setup

# %% tags=["cell-9"]
# ln -s loCRM1rep4 bloCRM1rep4
# ln -s loCRMmulti bloCRMmulti
raw.dir <- "/../../nlonfat/"
meta.dir2 <- tribble(
  ~dir, ~name,
  "/../../nlonfat", "retina1",
  "/../../nlonfat", "retina2",
  "/../../nlonfat", "retina3",
  "/../../nlonfat", "retina4",
  "/../../nlonfat", "CRM1FLP",
  "/../../nlonfat", "CRM1e6",
  "/../../nlonfat", "CRM1rep4",
  "/../../nlonfat", "mCherry",
  "/../../nlonfat", "CRMmulti",
  "/../../nlonfat", "Emerson",
  "/../../nlonfat", "Rep1_neg",
  "/../../nlonfat", "Rep1_pos",
  "/../../nlonfat", "Rep2_neg",
  "/../../nlonfat", "Rep2_pos",
#  "/../nlonfat", "Rep3_neg_wrong",
  "/../../nlonfat", "Rep3_neg",
  "/../../nlonfat", "Rep3_pos",
)
meta.dir2$name <- paste0( "blo", meta.dir2$name )

# %% tags=["cell-10"]
# Just check whether all files are present
meta.dir2 %>% 
  mutate(
    path       = map2_chr( dir, name, function(x, y) { paste0( root.dir, x, "/", y, "/outs/", "possorted_genome_bam_modsuffix2.bam" ) } ),
    library    = paste0( "", name )
  ) %>%
  dplyr::filter( !file.exists( path ) )

# %% [markdown]
# ## Cell Ranger quality metrics

# %% tags=["cell-12"]
meta.stats <- 
  meta.dir2 %>% 
  mutate(
    path       = map2_chr( dir, name, function(x, y) { paste0( root.dir, x, "/", y, "/outs/", "metrics_summary.csv" ) } ),
    library    = paste0( "", name )
  ) %>%
  dplyr::filter( file.exists( path ) ) %>%
  mutate(
    metrics = map( path, 
                   ~read_csv( .x, show_col_types = FALSE ) %>% as.data.frame()
                )
  ) %>%
  unnest( cols = "metrics")
meta.stats %>%
  mutate_if( is.numeric,
    scales::comma_format()
  )  %>%
  select( -path ) %>%
# DT::datatable()   # still not pretty
  mutate(across(everything(), as.character)) %>%
  pivot_longer(cols = -name,
               names_to = "Category",
               values_to = "Value") %>%
  pivot_wider(names_from = name,
              values_from = Value)

# %% tags=["cell-13"]
meta.stats <- 
  meta.dir2 %>% 
  mutate(
    path       = map2_chr( dir, name, function(x, y) { paste0( root.dir, x, "/", y, "/outs/", "metrics_summary.csv" ) } ),
    library    = paste0( "", name )
  ) %>%
#  dplyr::filter( file.exists( path ) ) %>%
  mutate(
    metrics = map( path, 
                   ~read_csv( .x, show_col_types = FALSE ) %>% as.data.frame()
                )
  ) %>%
  unnest( cols = "metrics")
meta.stats %>%
  mutate_if( is.numeric,
    scales::comma_format()
  )  %>%
  select( -path ) %>%
# DT::datatable()   # still not pretty
  mutate(across(everything(), as.character)) %>%
  pivot_longer(cols = -name,
               names_to = "Category",
               values_to = "Value") %>%
  pivot_wider(names_from = name,
              values_from = Value)

# %% [markdown]
# ## Load 10X data into Seurat

# %% tags=["cell-16"]
tictoc::tic()
gex <-
  meta.dir2 %>%
  mutate(
    path       = map2_chr( dir, name, function(x, y) { paste0( root.dir, x, "/", y, "/outs/", "metrics_summary.csv" ) } ),
    library    = paste0( "", name )
  ) %>%
  dplyr::filter( file.exists( path ) ) %>%
  mutate(
    seurat = map2( 
              dir, library, 
              function( d, l ) {
                print(glue::glue("{root.dir}{d}/{l}/outs/filtered_feature_bc_matrix/" ))  # log
                # Original should be lo
                m <- Read10X(
                    glue::glue("{root.dir}{d}/{l}/outs/filtered_feature_bc_matrix/" )
                )
                s <- CreateSeuratObject( 
                  counts = m, 
                  project = l
                )
                s$library <- l
                s$barcode <- colnames(s)
                return(s)
             }
          )
  )
tictoc::toc()

# %% [markdown]
# ## Barcode filtering and library QC

# %% tags=["cell-17"]
# Clean up some mess of new cellranger run by filtering CBs
# This is coming from original individual runs
retina <- readRDS( "../../data/20240815_01_retina_sctransform.rds" ) # This is a merged library, so keep an eye of any barcode changes

# %% tags=["cell-18"]
gex2 <-
gex %>%
mutate(
    seurat = map(
        seurat,
        function(s) {
            l <- gsub("^blo", "", s@project.name)
            l <- gsub("_wrong$", "", l )
            count.before <- ncol(s)
            s <-
            subset(
                s,
                cells = paste0(
                    retina@meta.data %>%
                    rownames_to_column( "bc" ) %>%
                    mutate(
                        library2 = ifelse( is.na( library ), "Emerson", library )
                    ) %>%
                    dplyr::filter( library2 == l ) %>%
                    pull( bc ) %>%
                    substring( 0, 16 ),
                    "-1"
                )
            )
            count.original <- 
            retina@meta.data %>%
            rownames_to_column( "bc" ) %>%
            mutate(
                library2 = ifelse( is.na( library ), "Emerson", library )
            ) %>%
            dplyr::filter( library2 == l ) %>%
            pull( bc ) %>% 
            length() 
            count.after <- ncol(s)
            print(glue::glue("{l} :: original {count.original}. {count.before} -> {count.after}") )
            return(s)
        }
    )
)

# %% tags=["cell-22"]
gex$seurat[[which(gex$name == "bloRep3_neg")]] <- gex2$seurat[[which(gex$name == "bloRep3_neg")]] 
gex$seurat[[which(gex$name == "bloRep3_pos")]] <- gex2$seurat[[which(gex$name == "bloRep3_pos")]] 
gex$seurat[[which(gex$name == "bloRep2_pos")]] <- gex2$seurat[[which(gex$name == "bloRep2_pos")]]

# %% tags=["cell-23"]
rm(gex2)

# %% [markdown]
# ## Barcode overlap analysis

# %% [markdown]
# ## Barcode overlap
# Want to make sure that there is no suspicious cellular barcode overlap between libraries (due to index hopping mainly).

# %% tags=["cell-25"]
barcode_lists <- map(
    gex$seurat,
    function(x) {
#        Cells( subset( x, subset = scDblFinder.class == "singlet" ) )
        Cells( x )
    }
)
names(barcode_lists) <- map_chr( gex$seurat, ~.x@project.name )
# Generate the upset plot directly from the list of sets
options( repr.plot.width = 14, repr.plot.height = 7 )
UpSetR::upset(
    UpSetR::fromList(barcode_lists),
    nsets = length(barcode_lists),
    nintersects = NA,
    keep.order = TRUE,
    sets = names(barcode_lists),
    sets.bar.color = "#56B4E9",
    matrix.color = "#56B4E9",
    point.size = 3,
    line.size = 1,
    mainbar.y.label = "Number of Shared Barcodes",
    sets.x.label = "Total Barcodes per Library",
    order.by = "freq",
    show.numbers = TRUE,
    text.scale = c(1.5, 1.5, 1.2, 1, 1.2, 1),
    number.angles = 0
)

# %% [markdown]
# It looks like that there is some significant overlap barcodes that indicate "ghost" libraries.

# %% [markdown]
# ### Fix Rep2_pos duplicate ghosts.

# %% tags=["cell-29"]
# Look for duplicated barcode names.
map_dfr(
    gex$seurat,
    ~.x@meta.data
) %>%
group_by( library ) %>%
mutate(
    total = n()
) %>%
group_by( barcode ) %>%
dplyr::filter( n() > 1 ) %>%
ungroup() %>%
dplyr::count( library, total )

# %% [markdown]
# Quite suspicious that `bloRep2_pos` contains information from `bloRep3_pos`. Also there is suspicions that `bloRep3_neg` also has problems.

# %% tags=["cell-31"]
problem.barcode <-
map_dfr(
    gex$seurat,
    ~.x@meta.data
) %>%
group_by( library ) %>%
mutate(
    total = n()
) %>%
group_by( barcode ) %>%
dplyr::filter( n() > 1 ) %>%
dplyr::filter( library == "bloRep2_pos" ) %>%
pull( barcode )

# %% tags=["cell-32"]
gex$seurat[[which(gex$name == "bloRep3_pos")]]@meta.data %>% 
dplyr::filter( barcode %in% problem.barcode ) %>%
head()

# %% tags=["cell-33"]
gex$seurat[[which(gex$name == "bloRep2_pos")]]@meta.data %>% 
dplyr::filter( barcode %in% problem.barcode ) %>%
head()

# %% tags=["cell-34"]
gex$seurat[[which(gex$name == "bloRep2_pos")]] <- 
subset( 
    gex$seurat[[which(gex$name == "bloRep2_pos")]],
    subset = barcode %ni% problem.barcode
)

# %% tags=["cell-35"]
# Look for duplicated barcode names.
map_dfr(
    gex$seurat,
    ~.x@meta.data
) %>%
group_by( library ) %>%
mutate(
    total = n()
) %>%
group_by( barcode ) %>%
dplyr::filter( n() > 1 ) %>%
ungroup() %>%
dplyr::count( library, total )

# %% [markdown]
# ## Doublet detection (scDblFinder)

# %% [markdown]
# ## Doublet identification: scDblFinder
# https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html
# The package is coming from Aron Lun, and also [benchmark](https://f1000research.com/articles/10-979) seem to suggest it supersedes DoubletFinder. Note this could be still be benefitted if we have genotyping reconstruction working.
# We could actually run the doublet finder with information from the genotyping call, but instead, I want to inform the genotyping call by cleaning up doublets here.
#
# NOTE: scDblFinder is stochastic (generates random synthetic doublets to
# train its classifier). Without a fixed set.seed(), re-runs will flag
# slightly different cells as doublets (~1-2% variation). This cascades
# into downstream Harmony/Leiden clustering and RPC subset sizes.
# Validated 2026-03-04: rerun yielded 83,915 cells (vs 85,135 original),
# a 1.4% difference attributable to this stochasticity.

# %% tags=["cell-38"]
require(scDblFinder)
gex <- gex %>%
mutate(
    seurat = map(
        seurat,
        function(s) {
            print(s@project.name)
            s <- NormalizeData( s )
            require(Seurat)
            require(SingleCellExperiment)
            require(scater)
            counts_matrix     <- LayerData( s, layer = "counts" )
            normalized_matrix <- LayerData( s, layer = "data" )
            sce <- SingleCellExperiment(
                assays = list(
                    counts = counts_matrix,
                    logcounts = normalized_matrix
                )
            )
            sce_colData  <- s@meta.data
            colData(sce) <- DataFrame(sce_colData)
            # If we have some doublet information from orthogonal approach we could do here
            # sce <- scDblFinder(
            #     sce,
            #     knownDoublets = "multiseq.doublet",
            #     knownUse = "discard"
            # )
            sce <- scDblFinder(
                sce
            )
            s$scDblFinder.score <- sce$scDblFinder.score
            s$scDblFinder.class <- sce$scDblFinder.class
            return(s)
        }
    )
)

# %% tags=["cell-39"]
map_dfr(
    gex$seurat,
    ~.x@meta.data %>%
    dplyr::count( library, scDblFinder.class )
) %>%
pivot_wider( names_from = scDblFinder.class, values_from = n )

# %% [markdown]
# ## Intron/exon ratio

# %% [markdown]
# ## intron/exon ratio
# Recently, too **low** intron/exon ratio has been found to be problematic (cell debris instead of real cell with nucleus)

# %% tags=["cell-42"]
tictoc::tic()
gex <-
  gex %>%
  mutate(
    seurat = pmap( 
              list( dir, name, seurat ),
              function( d, l, s ) {
                print(paste0(l, "\n"))  # log
                if( file.exists( glue::glue("{root.dir}{d}/{l}/outs/umi_type_metric_modsuffix2.tsv") ) ) {
                    umi_metric_df <- read_tsv( 
                        show_col_types = FALSE,
                        glue::glue("{root.dir}{d}/{l}/outs/umi_type_metric_modsuffix2.tsv" ),
                        comment = '#'
                    )
                    # print( head( umi_metric_df ) )
                    if (nrow(umi_metric_df) > 0) {
                        temp <-
                        s@meta.data %>%
                        dplyr::select( !any_of( c("Exonic", "Intronic", "Intergenic", "Total_UMIs") ) ) %>%
                        left_join(
                              umi_metric_df %>%
                              dplyr::rename( barcode = "Barcode" ) %>%
                              mutate( barcode = gsub("-[0-9A-Z]$", "-1", barcode ) ),
                              by = "barcode"
                        )
                        s <- AddMetaData(s, metadata = temp)
                    }
                } else {
                    print( "file not found" )
                }
                return(s)
             }
          )
  )
tictoc::toc()

# %% tags=["cell-43"]
gex$seurat[[1]]@meta.data %>% head()

# %% [markdown]
# ## Cell cycle scoring and metadata

# %% [markdown]
# ## Attach cell-specific meta information
# fraction of mitochondrial content, Cell phase etc.

# %% tags=["cell-45"]
gene.list <- rownames( gex$seurat[[1]] )
data('cc.genes')
setdiff( cc.genes$s.genes, gene.list ) # Mlf1ip is Cenpu
cc.genes$s.genes <- c( 
  intersect( gene.list, cc.genes$s.genes ), 
  "CENPU"                # MLF1IP
#  "POLR1B",               # 
)
setdiff( cc.genes$s.genes, gene.list )

# %% tags=["cell-46"]
setdiff( cc.genes$g2m.genes, gene.list )

# %% tags=["cell-47"]
# AURKB   no orthologues
# HN1     JPT1
# KIF20B  ENSGALG00000034232
# CDC25C  no orthologues
# CDCA2   ENSGALG00000045842
# CDCA8   ENSGALG00000033792
# CDCA8   ENSGALG00000047563
# CDCA8   ENSGALG00000053871
# PSRC1   not present
# CENPA   CenpA
# grep("ENSGALG00000004161", gene.list, value = T)
# grep("CENP", gene.list, value = T)
# Past
# CKAP2L  ENSGALG00000050615
# ANLN    ENSGALG00000043642
# LBR     TM7SF2
# CENPE   ENSGALG00000013208
# CTCF    CTCFL
cc.genes$g2m.genes <- c( 
  intersect( gene.list, cc.genes$g2m.genes ), 
  "JPT1"
)
setdiff( cc.genes$s.genes, gene.list )

# %% tags=["cell-48"]
# No mitochondrial gene annotation!
mito.genes <- c("ND1", "MT-ND2", "MT-CO1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "ND6", "CYTB" )
# cat /n/groups/tabin/chlee/singlecell/reference/star.reference.gg6a.nlonfat.retro/Gg6a_nlonfat_retro/genes/genes.gtf  | gawk '($1 == "MT" && $3 == "gene")' 
setdiff( mito.genes, gene.list )
mito.genes <- intersect( mito.genes, gene.list )

# %% tags=["cell-50"]
rbc.genes <- c("HBZ", "HBBR", "HBA1", "HBE", "HBE1", "HBBA")
# HBE1 shuold be present!
grep( "^HB", gene.list, value = T)
setdiff( rbc.genes, gene.list )
intersect( rbc.genes, gene.list )
rbc.genes <- intersect( rbc.genes, gene.list )

# %% tags=["cell-51"]
tictoc::tic()
gex <-
  gex %>%
  mutate(
    seurat = map( 
              seurat,
              function( s ) {
                s$percent.rbc <- Matrix::colSums(
                    LayerData(
                        object = s,
                        layer = "counts"
                    )[rbc.genes, ]
                ) / Matrix::colSums(
                    LayerData(
                        object = s,
                        layer = "counts"
                    )
                )
                s$percent.mito <- Matrix::colSums(
                    LayerData(
                        object = s,
                        layer = "counts"
                    )[mito.genes, ]
                ) / Matrix::colSums(
                    LayerData(
                        object = s,
                        layer = "counts"
                    )
                )
                s$percent.W <- Matrix::colSums(
                    LayerData(
                        object = s,
                        layer = "counts"
                    )[W.genes, ]
                ) / Matrix::colSums(
                    LayerData(
                        object = s,
                        layer = "counts"
                    )
                )
                s$percent.Z <- Matrix::colSums(
                    LayerData(
                        object = s,
                        layer = "counts"
                    )[Z.genes, ]
                ) / Matrix::colSums(
                    LayerData(
                        object = s,
                        layer = "counts"
                    )
                )
                s <- NormalizeData( s ) %>%
                CellCycleScoring(
                      s.features = cc.genes$s.genes,
                      g2m.features = cc.genes$g2m.genes,
                      set.ident = TRUE
                )
                return(s)
             }
          )
  )
tictoc::toc()

# %% [markdown]
# ## Vireo genotyping

# %% [markdown]
# # Genotyping
# This is following up the vireo genotyping pipeline again. See `nlonfat/` original directory for details.

# %% tags=["cell-60"]
map_dfr(
    gex$seurat,
    ~.x@meta.data
) %>%
dplyr::count( library, scDblFinder.class ) %>%
pivot_wider( names_from = scDblFinder.class, values_from = n ) %>%
mutate(
    Dbl.fraction = doublet / (singlet + doublet )
)

# %% [markdown]
# Note that `bloCRMmulti` has significant doublet fraction (which scrublet unfortunately did not detect). So it is understandable that these vireo results did not work out well before.

# %% [markdown]
# ## Export singlet barcodes

# %% tags=["cell-63"]
gex$seurat[[1]]@meta.data %>% head()

# %% tags=["cell-64"]
data.dir
data.dir <- "../../R4.4.0"
root.dir

# %% tags=["cell-65"]
pmap(
    list(
        gex$dir,
        gex$name,
        gex$seurat
    ),
    function( d, n, s ) {
        s@meta.data %>%
        dplyr::filter( scDblFinder.class == "singlet" ) %>%
        dplyr::select( barcode ) %>%
        write_tsv( glue::glue( "{root.dir}{d}/{n}/outs/barcodes.singlet.txt" ), col_names = FALSE ) 
        glue::glue( "{root.dir}{d}/{n}/outs/barcodes.singlet.txt" )
    }
)

# %% [markdown]
# ## Genotyping survey and optimization

# %% [markdown]
# ## Optimal genotyping result survey

# %% tags=["cell-67"]
# To understand the expected number of genotypes, I can use the `bam_config.csv` file.
# Read the bam_config.csv file
bam_config <- read_csv(glue::glue("{root.dir}{raw.dir}bam_config.csv"), col_types = cols(notes = col_character()))
# Parse notes column and count genotypes per group (exploratory, not used downstream)
tryCatch({
  genotype_counts <- bam_config %>%
    filter(!is.na(notes), notes != "") %>%
    separate_rows(notes, sep = ";") %>%
    mutate(notes = str_trim(notes)) %>%
    filter(notes != "") %>%
    distinct(group_id, notes) %>%
    count(group_id, name = "n_genotypes") %>%
    arrange(group_id)
  print(genotype_counts)
  genotype_details <- bam_config %>%
    filter(!is.na(notes), notes != "") %>%
    separate_rows(notes, sep = ";") %>%
    mutate(notes = str_trim(notes)) %>%
    filter(notes != "") %>%
    distinct(group_id, notes) %>%
    arrange(group_id, notes)
  genotype_summary <- genotype_details %>%
    group_by(group_id) %>%
    summarise(
      n_genotypes = n(),
      genotype_list = paste(notes, collapse = "; "),
      .groups = "drop"
    )
  genotype_summary
}, error = function(e) {
  cat("Note: genotype survey skipped (notes column empty):", conditionMessage(e), "\n")
})

# %% tags=["cell-68"]
bam_config

# %% tags=["cell-69"]
genotyping.df <-
map_dfr(
    1:16,
    function(group.id) {
        genotyping.group <- glue::glue("group{group.id}")
        cellranger.id <- bam_config$cellranger_id[which(bam_config$group_id == group.id)]
        genotype.dir <- glue::glue("{root.dir}{raw.dir}{cellranger.id}/outs/pseudobulk_{genotyping.group}_merged/vireo_output/")
        tibble(
            string = list.files( glue::glue("{genotype.dir}") )
        ) %>%
        mutate(
            vireo.output.dir = glue::glue("{genotype.dir}{string}/"),
        ) %>%
        separate(string, into = c("prefix", "type", "N", "cell_count"),
                 sep = "_", remove = FALSE) %>%
        mutate(
            across(
                c(N, cell_count), 
                as.numeric
            )
        ) %>%
        select(-prefix, -type) %>%
        mutate(
            group = group.id,
            chunks = grep("chunk", list.files(vireo.output.dir), value = T)
        ) %>%
        dplyr::select( group, everything() )
    }
)
genotyping.df

# %% tags=["cell-70"]
genotyping.df <-
genotyping.df %>%
mutate(
    qc = map2(
        chunks,
        vireo.output.dir,
        function(chunk,v) {
            seeds <- list.files(glue::glue("{v}{chunk}"))
            map_dfr(
                seeds,
                function(s) {
                    log_file <- glue::glue("{v}{chunk}/{s}/_log.txt")
                    # Read theta values from _log.txt
                    if(file.exists(log_file)) {
                        # Read all lines from the file
                        log_content <- readLines(log_file)
                        # Extract log-likelihood
                        loglik_line <- log_content[1]
                        loglik_value <- as.numeric(gsub("logLik: ", "", loglik_line))
                        # Find the theta matrix (usually starts after "thetas:" line)
                        theta_start <- which(grepl("thetas:", log_content)) + 1
                        if(length(theta_start) > 0 && theta_start <= length(log_content)) {
                            # Extract and clean the matrix rows
                            # First row contains [[ at the beginning
                            row1_line <- log_content[theta_start]
                            # Last row contains ]] at the end 
                            row2_line <- log_content[theta_start + 1]
                            # Clean up the bracket symbols and split
                            row1_vals <- gsub("\\[\\[|\\]", "", row1_line)
                            row2_vals <- gsub("\\[|\\]\\]", "", row2_line)
                            # Convert to numeric vectors, filtering out empty strings
                            row1_nums <- as.numeric(strsplit(trimws(row1_vals), "\\s+")[[1]])
                            row2_nums <- as.numeric(strsplit(trimws(row2_vals), "\\s+")[[1]])
                            # Make sure we have 3 values for each row
                            if(length(row1_nums) == 3 && length(row2_nums) == 3) {
                                # Calculate mean for each genotype
                                theta_means <- row1_nums / (row1_nums + row2_nums)
                                # Calculate concentration
                                theta_concs <- row1_nums + row2_nums
                                # Calculate deviation from expected values [0, 0.5, 1]
                                expected_means <- c(0, 0.5, 1)
                                theta_deviation <- sum(abs(theta_means - expected_means))
                                # Create tibble with calculated values
                                return(tibble(
                                    chunk = chunk,
                                    seed = s,
                                    loglik = loglik_value,
                                    theta_mean_ref = theta_means[1],
                                    theta_mean_het = theta_means[2], 
                                    theta_mean_alt = theta_means[3],
                                    theta_conc_ref = theta_concs[1],
                                    theta_conc_het = theta_concs[2],
                                    theta_conc_alt = theta_concs[3],
                                    theta_deviation = theta_deviation
                                ))
                            }
                        }
                        # If we couldn't parse the theta matrix correctly, return just the loglik
                        return(tibble(
                            chunk = chunk,
                            seed = s,
                            loglik = loglik_value,
                            theta_mean_ref = NA,
                            theta_mean_het = NA, 
                            theta_mean_alt = NA,
                            theta_conc_ref = NA,
                            theta_conc_het = NA,
                            theta_conc_alt = NA,
                            theta_deviation = NA
                        ))
                    }
                }
            )
        }
    )
) %>%
unnest( cols = qc )
# Check how many runs we have per chunk
genotyping.df %>%
dplyr::count( group, N, cell_count, chunks )
# dplyr::count( chunk )

# %% tags=["cell-71"]
genotyping.df <-
genotyping.df %>%
mutate(
    summary.stats = pmap(
        list( vireo.output.dir, chunk, seed ),
        function( v, c, s ) {
            donor_file <- glue::glue("{v}{c}/{s}/donor_ids.tsv")
            if(file.exists(donor_file)) {
                donor_data <- read_tsv(
                    donor_file,
                    show_col_types = FALSE
                ) %>%
                dplyr::select( cell, donor_id ) %>%
                dplyr::count( donor_id ) %>%
                arrange( desc(n) )
                # Calculate summary statistics here
                total_cells <- sum(donor_data$n)
                doublet_cells <- sum(donor_data$n[str_detect(donor_data$donor_id, "doublet")])
                unassigned_cells <- sum(donor_data$n[str_detect(donor_data$donor_id, "unassigned")])
                assigned_cells <- total_cells - doublet_cells - unassigned_cells
                # Count unique donors with at least 5% of cells
                # (excluding doublets and unassigned)
                real_donors <- donor_data %>%
                    dplyr::filter(!str_detect(donor_id, "doublet|unassigned"))
                donors_above_5pct <- sum(real_donors$n >= (total_cells * 0.05))
                # Return both the detailed data and summary stats
                list(
                    detailed_counts = donor_data,
                    total_cells = total_cells,
                    doublet_cells = doublet_cells,
                    unassigned_cells = unassigned_cells,
                    assigned_cells = assigned_cells,
                    doublet_rate = doublet_cells / total_cells,
                    unassigned_rate = unassigned_cells / total_cells,
                    assigned_rate = assigned_cells / total_cells,
                    donors_above_5pct = donors_above_5pct
                )
            } else {
                list(
                    detailed_counts = tibble(),
                    total_cells = 0,
                    doublet_cells = 0,
                    unassigned_cells = 0,
                    assigned_cells = 0,
                    doublet_rate = 0,
                    unassigned_rate = 0,
                    assigned_rate = 0,
                    donors_above_5pct = 0
                )
            }
        }
    )   
)

# %% tags=["cell-72"]
bam_config

# %% tags=["cell-73"]
genotyping.df %>%
    mutate(
        total_cells = map_dbl(summary.stats, "total_cells"),
        doublet_cells = map_dbl(summary.stats, "doublet_cells"),
        unassigned_cells = map_dbl(summary.stats, "unassigned_cells"),
        assigned_cells = map_dbl(summary.stats, "assigned_cells"),
        doublet_rate = map_dbl(summary.stats, "doublet_rate"),
        unassigned_rate = map_dbl(summary.stats, "unassigned_rate"),
        assigned_rate = map_dbl(summary.stats, "assigned_rate"),
        donors_above_5pct = map_dbl(summary.stats, "donors_above_5pct")
    ) %>%
dplyr::select( group, N, seed, loglik, assigned_rate, donors_above_5pct ) %>%
left_join(
    bam_config %>%
    dplyr::rename( group = group_id )
) %>%
ggplot( aes( y = assigned_rate, x = log10(-loglik), colour = factor(N) ) ) +
geom_point() +
facet_wrap( ~ paste0( group, "_", cellranger_id ) ) +
expand_limits( y = 0 ) +
scale_y_continuous( expand = expansion( add = 0, mult = c(0, 0.05) ) )

# %% tags=["cell-74"]
genotyping.df %>%
    mutate(
        total_cells = map_dbl(summary.stats, "total_cells"),
        doublet_cells = map_dbl(summary.stats, "doublet_cells"),
        unassigned_cells = map_dbl(summary.stats, "unassigned_cells"),
        assigned_cells = map_dbl(summary.stats, "assigned_cells"),
        doublet_rate = map_dbl(summary.stats, "doublet_rate"),
        unassigned_rate = map_dbl(summary.stats, "unassigned_rate"),
        assigned_rate = map_dbl(summary.stats, "assigned_rate"),
        donors_above_5pct = map_dbl(summary.stats, "donors_above_5pct")
    ) %>%
dplyr::select( group, N, seed, loglik, theta_deviation, assigned_rate, donors_above_5pct ) %>%
left_join(
    bam_config %>%
    dplyr::rename( group = group_id )
) %>%
ggplot( aes( y = assigned_rate, x = theta_deviation, colour = factor(N) ) ) +
geom_point() +
facet_wrap( ~ paste0( group, "_", cellranger_id ), scale = "free_x" ) +
expand_limits( y = 0 ) +
scale_y_continuous( expand = expansion( add = 0, mult = c(0, 0.05) ) )

# %% tags=["cell-75"]
genotyping.df %>%
    mutate(
        total_cells = map_dbl(summary.stats, "total_cells"),
        doublet_cells = map_dbl(summary.stats, "doublet_cells"),
        unassigned_cells = map_dbl(summary.stats, "unassigned_cells"),
        assigned_cells = map_dbl(summary.stats, "assigned_cells"),
        doublet_rate = map_dbl(summary.stats, "doublet_rate"),
        unassigned_rate = map_dbl(summary.stats, "unassigned_rate"),
        assigned_rate = map_dbl(summary.stats, "assigned_rate"),
        donors_above_5pct = map_dbl(summary.stats, "donors_above_5pct")
    ) %>%
dplyr::select( group, N, seed, loglik, doublet_rate, donors_above_5pct ) %>%
left_join(
    bam_config %>%
    dplyr::rename( group = group_id )
) %>%
ggplot( aes( x = donors_above_5pct, y = doublet_rate, colour = factor(N) ) ) +
geom_point() +
facet_wrap( ~ paste0( group, "_", cellranger_id ) ) +
expand_limits( y = 0 ) +
scale_y_continuous( expand = expansion( add = 0, mult = c(0, 0.05) ) )

# %% [markdown]
# * `Rep2` series seems to have issues. Thinking back, I realize that all the `RepX` series positive and negatives should have been combined together!
# * `CRMmulti` as well as `CRM1FLP` seems to have found the spot (2)
# * `CRM1e6`, `CRM1rep4` failed but may indicate that this is a single embryo.
# So in sum,
# * bloretina1 (3), bloretina2 (3), bloretina3 (3), bloretina (3), bloRep1 (4), bloRep2 (1?), bloRep3 (4), CRM1FLP (2), CRM1multi (2), blomCherry (2), bloEmerson (4) likely?

# %% tags=["cell-78"]
genotyping.df %>%
    mutate(
        total_cells = map_dbl(summary.stats, "total_cells"),
        doublet_cells = map_dbl(summary.stats, "doublet_cells"),
        unassigned_cells = map_dbl(summary.stats, "unassigned_cells"),
        assigned_cells = map_dbl(summary.stats, "assigned_cells"),
        doublet_rate = map_dbl(summary.stats, "doublet_rate"),
        unassigned_rate = map_dbl(summary.stats, "unassigned_rate"),
        assigned_rate = map_dbl(summary.stats, "assigned_rate"),
        donors_above_5pct = map_dbl(summary.stats, "donors_above_5pct")
    ) %>%
left_join(
    bam_config %>%
    dplyr::rename( group = group_id )
) %>%
dplyr::select( 
    cellranger_id,
    group, N, seed, chunk,
    loglik, theta_deviation, 
    donors_above_5pct,
    assigned_rate, 
    doublet_rate,
    unassigned_rate
) %>%
group_by( cellranger_id ) %>%
arrange( group, desc( assigned_rate ) ) %>%
dplyr::filter( row_number() <= 3 ) # top 3

# %% [markdown]
# ## Assign genotypes to cells

# %% [markdown]
# ## Assignment of provisional genotypes

# %% tags=["cell-82"]
lookup.table <-
tribble(
    ~group, ~N, ~seed,
    1, 5, "seed_746352",
    2, 5, "seed_314097",
    3, 5, "seed_121335",
    4, 5, "seed_468164",
    5, 3, "seed_114884",
#   6,
#   7,
    8, 5, "seed_340412",
    9, 3, "seed_254075",
   10, 5, "seed_737487",
   11, 5, "seed_639570",
   12, 5, "seed_319565",
#  13,
#  14,
   15, 5, "seed_680849",
   16, 5, "seed_11706"
) %>%
left_join(
    genotyping.df %>%
    dplyr::select( group, N, seed, chunk, vireo.output.dir )
) %>%
left_join(
    bam_config %>%
    dplyr::rename( group = group_id )
)
# mutate(
#     genotypes = pmap_dfr(
#         group, N, seed,
#         function( group, N, seed ) {
#         }
#     )
# )
lookup.table

# %% tags=["cell-83"]
gex <-
gex %>%
mutate(
    seurat = map(
        seurat,
        function(s) {
            print( s@project.name )
            library.name <- s$library %>% unique()
            # s$genotype <- NULL
            v <- lookup.table %>%
            dplyr::filter( cellranger_id == library.name ) %>%
            pull( vireo.output.dir )
            chunk <- lookup.table %>%
            dplyr::filter( cellranger_id == library.name ) %>%
            pull( chunk )
            seed <- lookup.table %>%
            dplyr::filter( cellranger_id == library.name ) %>%
            pull( seed )
            if (length(v) > 0) {
                vireo_result <- read_tsv( glue::glue("{v}{chunk}/{seed}/donor_ids.tsv") ) %>%
                dplyr::rename( barcode = cell ) %>%
                mutate(
                    barcode = gsub("-[0-9A-Z]$", "-1", barcode )
                )
                # print(donor_file %>% head())
                s$genotype <- NULL
                s$genotype <- with(
                    s@meta.data,
                    s@meta.data %>%
                    left_join(
                        vireo_result,
                        by = c("barcode")
                    ) %>%
                    pull( donor_id )
                )
                s$prob_max <- NULL
                s$prob_max <- with(
                    s@meta.data,
                    s@meta.data %>%
                    left_join(
                        vireo_result,
                        by = c("barcode")
                    ) %>%
                    pull( prob_max )
                )
                s$n_vars <- NULL
                s$n_vars <- with(
                    s@meta.data,
                    s@meta.data %>%
                    left_join(
                        vireo_result,
                        by = c("barcode")
                    ) %>%
                    pull( n_vars )
                )                
            } else {
                s$genotype <- "singlet"
                s$prob_max <- 1
                s$n_vars <- 0
            }
            s
        }
    )
)

# %% [markdown]
# ## Checkpoint 1: save raw Seurat

# %% tags=["cell-91"]
tictoc::tic()
saveRDS( gex, file=glue::glue( "../../data/{intermediate.prefix}01_gex.rds" ) )
tictoc::toc()

# %% [markdown]
# ## Merge and filter cells

# %% [markdown]
# # Clustering

# %% [markdown]
# ## Merge & Filter
# Given the continuous instability of genotype assignment we will drop any further filtering process and use the whole cell population as is for clustering. In such case, the critical issue is to control for the sex, as there is no dosage compensation at the transcriptional level in aves.

# %% tags=["cell-95"]
retina <- merge(
    x = gex$seurat[[1]],
    y = gex$seurat[c(2:nrow(gex))]
)

# %% tags=["cell-96"]
retina@meta.data %>% nrow()
temp1 <- retina@meta.data %>%
dplyr::count( library )
temp1

# %% tags=["cell-97"]
retina <- subset( 
    retina, 
    subset = scDblFinder.class == "singlet" & 
    grepl( "(donor|singlet)", genotype ) &
    Intronic / (Intronic + Exonic) > 0.1
)

# %% tags=["cell-98"]
retina@meta.data %>% 
nrow()
temp2 <- retina@meta.data %>%
dplyr::count( library )
temp2

# %% [markdown]
# ## Unintegrated clustering (PCA, UMAP)

# %% [markdown]
# ## Unintegrated clustering
# SCTransform seems to actually **introduce** batch effect for some reason.

# %% tags=["cell-100"]
options(future.globals.maxSize = 3e+09)
retina <-
retina %>%
NormalizeData() %>%
FindVariableFeatures() %>%
#SCTransform() %>%
#ScaleData( vars.to.regress = "CC.Difference" ) %>%
ScaleData( vars.to.regress = c("G2M.Score", "S.Score", "percent.mito", "percent.W", "percent.Z") ) %>% # control by sex
RunPCA( npcs = 50, verbose = F)

# %% tags=["cell-101"]
options( repr.plot.width = 14, repr.plot.height = 7 )
ElbowPlot( retina, ndims = 50 )

# %% tags=["cell-102"]
retina <- 
retina %>% 
FindNeighbors( dims = 1:30, reduction = "pca" ) %>%
FindClusters( resolution = 0.2, cluster.name = "unintegrated_clusters" )

# %% tags=["cell-103"]
retina <- RunUMAP( 
    retina, 
    dims = 1:30, 
    reduction = "pca", 
    reduction.name = "umap.unintegrated",
    verbose = F
)

# %% [markdown]
# ## Harmony batch integration

# %% [markdown]
# ## Harmony Integration

# %% tags=["cell-108"]
retina <- IntegrateLayers(
    object = retina,
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = "harmony",
    verbose = T
)

# %% tags=["cell-109"]
options( repr.plot.width = 14, repr.plot.height = 7 )
ElbowPlot( retina, ndims = 50 )

# %% tags=["cell-110"]
retina <- 
retina %>%
FindNeighbors( reduction = "harmony", dims = 1:30 ) %>%
FindClusters( resolution = 0.2, cluster.name = "harmony_clusters.leiden.0.2", algorithm = 4 ) %>%
RunUMAP( reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony" )

# %% tags=["cell-111"]
retina <- 
retina %>%
# FindNeighbors( reduction = "harmony", dims = 1:15 ) %>%
FindClusters(
    resolution = 0.1,
    cluster.name = "harmony_clusters.leiden.0.1",
    algorithm = 4
) %>%
FindClusters( 
    resolution = 0.4, 
    cluster.name = "harmony_clusters.leiden.0.4",
    algorithm = 4 
)

# %% [markdown]
# ## Cluster annotation

# %% [markdown]
# ## Clustering annotation

# %% tags=["cell-118"]
options( repr.plot.width = 14, repr.plot.height = 7 )
Idents( retina ) <- "harmony_clusters.leiden.0.1"
p1 <-
DimPlot(
    retina,
    group.by = c("harmony_clusters.leiden.0.4","harmony_clusters.leiden.0.1"),
    reduction = "umap.harmony",
    label = T
) +
theme(
    axis.text = element_blank()
)
p1

# %% tags=["cell-119"]
options( repr.plot.width = 14, repr.plot.height = 7 )
VlnPlot(
    retina,
    group.by = "harmony_clusters.leiden.0.4",
    features = c("percent.mito", "nFeature_RNA"),
    pt.size = 0
)

# %% tags=["cell-120"]
temp <- presto::wilcoxauc( JoinLayers(retina), 'harmony_clusters.leiden.0.4' )

# %% tags=["cell-121"]
# Antiviral response
temp %>%
dplyr::filter( group %in% c(12) ) %>%
arrange( group, pct_out - pct_in ) %>%
dplyr::group_by( group ) %>%
dplyr::filter( row_number() <= 20 ) %>%
dplyr::select( feature, group, pct_in, pct_out )

# %% tags=["cell-122"]
options( repr.matrix.max.rows = 1000 )
temp <- presto::wilcoxauc( JoinLayers(retina), 'harmony_clusters.leiden.0.1' )
temp %>%
arrange( group, pct_out - pct_in ) %>%
dplyr::group_by( group ) %>%
dplyr::filter( row_number() <= 5 ) %>%
dplyr::select( feature, group, pct_in, pct_out )

# %% tags=["cell-123"]
retina$annotation <- with(
    retina@meta.data,
    case_when(
        harmony_clusters.leiden.0.4 == 12 ~ "Infected",
        harmony_clusters.leiden.0.1 == 1 ~ "RPC",
        harmony_clusters.leiden.0.1 == 2 ~ "OTX2+ neurogenic",
        harmony_clusters.leiden.0.1 == 3 ~ "THRB+RXRG+ cone",
        harmony_clusters.leiden.0.1 == 4 ~ "Horizontal/Amacrine", # TFAP2A/B
        harmony_clusters.leiden.0.1 == 5 ~ "Retinal ganglion", # ISL1
        harmony_clusters.leiden.0.1 == 6 ~ "Vascular"
    )
)

# %% [markdown]
# ## Checkpoint 2: save integrated Seurat

# %% tags=["cell-134"]
tictoc::tic()
saveRDS( retina, file=glue::glue( "../../data/{intermediate.prefix}01_retina.rds" ) )
tictoc::toc()

# %% [markdown]
# ## Export to h5ad (full dataset)
# Uses `export_seurat()` from `00_utils.R`

# %% tags=["cell-138"]
export_seurat( retina, "../../data/20250604chick", prefix="20250604_" )
