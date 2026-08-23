# Standard processing of single-cell RNA datasets

## Processing of chick developing retina datasets

For chick dataset, 15 libraries from 10X Chromium v2, v3.1 Single cell (single, dual index) protocol, coming from 1-2 days of ex ovo culture from E5-6 chick retina has been processed in house according according to 10X Genomics user guidelines. Some of the samples had retroviral GFP barcode infection, or CRE-GFP plasmid transfection for different study purposes but such reads were ignored. One additional publicly available (GSE142244) was included, reprocessed starting from FASTQ files and undergoing re-processing with generation of count matrix (cellranger v9.0.0) with the identical chick reference.

From the generated count matrix, scDblFinder (Lun et al. 2022; Germain P, Lun A, Garcia Meixide C, Macnair W, Robinson M (2022). “Doublet identification in single-cell sequencing data using scDblFinder.” f1000research. doi:10.12688/f1000research.73600.2. ) identified doublets that were removed. The doublet-filtered data underwent SNP detection (cellsnp-lite v1.2.3), genotype-based demultiplexing (vireoSNP v0.5.8) for further doublet removal, and identification of unique genotypes. For genotype-based demultiplexing, at least 100 different random seeds were applied to obtain the best fit. The identity of unique genotypes (used for pseudobulk later) were validated with the sex chromosome genes (W/Z). Simple intro and exon ratio was computed to identify droplet transcriptomes that had low intron presence to further filter putative debris (Montserrat-Ayuso & Esteve-Codina, 2024). The amount of ambient RNA was assessed using SoupX (v1.6.2, Young & Behjati, 2020), but found to be similar across and at low level, so no further action was taken to correct for the ambient RNA effect. Overall, among initial 110,061 cells, 85,135 single-cell transcriptome were subject to standard processing in Seurat v5 framework (Hao et al., 2023, Hao et al., 2021, Stuart et al. 2019, Butler et al., 2018, Satija et al., 2015) for normalization, integration (Harmony, Korsunsky et al. 2019), and graph-based clustering (leiden, Traag et al. 2019). For clustering, cell cycle score (G2M/S Score, computed from each libraries separately), mitochondrial content, and W/Z gene content were regressed out. Leiden resolution=0.1 showed good demarcation of cell clusters, including the retinal progenitor (RPC) populations, and identification of cell transcriptomes with signatures of viral infection (due to retroviral barcode infection), which were filtered out.

The resulting RPC population amount to 23,706 single-cell transcriptomes, derived from 30 individual embryos (based on inferred genotypes by vireo). From the data, Dorsal, Ventral, Nasal and Temporal scores were computed with the combination of {TBX5, TBX2, TBX3, ALDH1A1, EFNB2, EFNB1}, {VAX1, CHRDL1, ALDH1A3}, {FOXG1, SOHO-1, HMX1, EFNA5, EFNA2}, and {FOXD1, EPHA3} respectively. Subsequently, the DV score and NT score was constructed by simply subtracting Dorsal score to Ventral score, Nasal score with Temporal score respectively. There was little difference of score distribution whether the score was computed in library or technology level (10X v2, v3.1) or altogether, so the integrated dataset as a whole was used to compute the score.

## Processing of human developing retina datasets

Publicly available datasets from three datasets were processed (GSE246169, GSE138002, GSE234963). Only whole retina samples between PCW7-11 were used for integration. The count matrices were merged, and subject to similar process as chick dataset of normalization, integration, and graph-based clustering. The only difference is not regressing for sex chromosome gene content since mammal have sex chromosome dosage compensation. After clustering, RPC population was identified based on the publicly available meta information overlap (retinal progenitors and mitotic cell population, which overlap with RPC cluster after cell-cycle regression). For the RPC population, the same set of genes were used to construct the DV score and NT score.

## Processing of mouse developing retina datasets

Publicly available dataset from GSE139904 (E13) and GSE118614 (E14,16) were processed accordingly. The count matrices were merged after harmonization of gene names. Mitochondrial content and cell-cycle scores were regressed out during the normalization and scaling process. After clustering, RPC population was identified based on the publicly available meta information overlap. The subsetted RPC population was used to construct the DV score and NT score.


# Spatial Expression Analysis

## Spatial Binning and Smoothing

Most of the spatial analyses were conducted in python (3.11), scanpy (v1.9.6), anndata (v0.10.3), and custom python module (github URL). 

To reconstruct spatial gene expression patterns, we discretized the continuous DV and NT coordinates into regular grids. Species-specific parameters were optimized based on cell density and spatial resolution requirements:

For chicken samples, we used 51×51 spatial bins (2,601 total), requiring a minimum of 3 cells per bin and 20 gene counts. For human samples, we employed 40×40 bins (1,600 total) with thresholds of 3 cells and 15 gene counts. Mouse samples utilized 51×51 bins with higher thresholds of 5 cells and 30 gene counts due to greater cell density.

Gaussian smoothing (σ = 1.0) was applied to the binned expression values to reduce noise while preserving spatial patterns. Expression values exceeding the 93rd-95th percentile (species-dependent) were clipped to minimize outlier effects in visualization

### Spatial Pattern Detection

We identified spatially variable genes using a custom greedy algorithm that maximizes spatial autocorrelation. For each gene, we calculated Moran's I statistic to quantify spatial clustering. Genes with significant spatial variance (p < 0.01 after Benjamini-Hochberg correction) were retained for pattern analysis.

Spatial similarity between genes was assessed using cosine similarity of their smoothed expression patterns. Genes with correlation coefficients > 0.9 were considered to have highly similar spatial distributions. Pre-defined spatial anchor genes were used to orient patterns along anatomical axes (Dorsal: *BMP4*, *TBX5*; Ventral: *VAX2*, *PAX2*; Temporal: *FOXG1*, *SOX2*; Nasal: *FOXD1*, *ALX1* for human/chicken with appropriate orthologs for mouse).

### Differentially expressed genes in regions

The differentially expressed gene analysis used glmGamPoi (Ahlmann-Eltze & Huber et al. 2020) to first pseudobulk the single cell transcriptome based on the defined region. Only pseudobulks containing more than 100 cells were used.


## Cross-Species Comparative Analysis

### Orthology Mapping

To enable cross-species comparisons, we mapped genes between species using a dual-database approach. Primary orthology information was retrieved from MyGene.info API (accessed via mygene Python client v3.2.2), with Ensembl BioMart serving as a fallback source. Species were identified using NCBI taxonomy IDs: 9606 (human), 10090 (mouse), and 9031 (chicken).

Gene naming conventions were strictly maintained with uppercase nomenclature for human and chicken genes (e.g., *FGF8*, *TBX5*, *PAX6*) and sentence case for mouse genes (e.g., *Fgf8*, *Tbx5*, *Pax6*). Orthology mapping achieved 40-80% coverage depending on the gene set and species pair.

### Conservation Analysis

We performed systematic comparison of spatial expression patterns across species. For each spatially variable gene cluster identified within a species, we identified orthologous genes in the other two species and compared their spatial distributions. Conservation was quantified using Pearson correlation of spatial expression patterns, with r > 0.7 indicating conserved spatial expression.

Gene sets showing conserved spatial patterns were subjected to Gene Ontology (GO) enrichment analysis using g:Profiler (g:GOSt) to identify biological processes under spatial regulatory control across vertebrate evolution.

## Statistical Analysis

### Correlation and Similarity Metrics

Pearson correlation coefficients were calculated for all pairwise gene comparisons within spatial bins. Spatial patterns were considered significantly correlated at r > 0.9 (p < 0.001). Cosine similarity was used as an additional metric for high-dimensional expression pattern comparisons.

### Multiple Testing Correction

All statistical tests involving multiple comparisons were adjusted using the Benjamini-Hochberg false discovery rate (FDR) method with a significance threshold of q < 0.05.

### Software and Reproducibility

Complete analysis code, including all parameters and random seeds, is available at [repository URL]. Processed data files in h5ad format containing expression matrices, spatial coordinates, and cell annotations are deposited at [data repository].

## Data Availability

Raw sequencing data and processed expression matrices are available through [accession numbers]. 