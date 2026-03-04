#!/usr/bin/env python3
"""
Spatial Gene Expression Analysis Module
======================================

This module provides functions for analyzing spatial gene expression patterns
across different species and performing cross-species comparisons.

Author: Generated for spatial transcriptomics analysis
"""

from __future__ import annotations
import sys
import os
import warnings
import yaml
from dataclasses import dataclass
from typing import List, Dict, Sequence, Optional, Tuple, Union
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.ndimage import gaussian_filter
import sklearn.metrics
from sklearn.preprocessing import normalize
import matplotlib.pyplot as plt
import anndata as ad

warnings.filterwarnings('ignore')

# Load control genes from YAML file
def load_control_genes(config_file: str = "control_genes.yaml") -> Dict[str, Dict[str, List[str]]]:
    """Load control genes configuration from YAML file"""
    # If config_file is just a filename, look for it in the script directory
    if not os.path.isabs(config_file) and os.path.dirname(config_file) == "":
        script_dir = os.path.dirname(os.path.abspath(__file__))
        config_file = os.path.join(script_dir, config_file)
    
    # Convert to absolute path for display
    config_file_abs = os.path.abspath(config_file)
    
    try:
        with open(config_file, 'r') as f:
            control_genes = yaml.safe_load(f)
            print(f"✓ Loaded control genes from: {config_file_abs}")
            return control_genes
    except FileNotFoundError:
        print(f"Warning: {config_file_abs} not found. Using minimal default control genes.")
        # Minimal fallback if file is missing
        return {
            'chick': {'Fovea': ['FGF8'], 'Dorsal': ['TBX5'], 'Ventral': ['EPHB2']},
            'human': {'Fovea': ['FGF8'], 'Dorsal': ['TBX5'], 'Ventral': ['EPHB2']},
            'mouse': {'Fovea': ['Fgf8'], 'Dorsal': ['Tbx5'], 'Ventral': ['Ephb2']}
        }

# Load the control genes (users can edit control_genes.yaml to customize)
CONTROL_GENES = load_control_genes()

@dataclass
class SpatialAnalysisParams:
    """Parameters for spatial expression analysis"""
    # Spatial binning parameters
    bin_size: int = 51
    adaptive_binning: bool = False
    target_cells_per_bin: int = 10
    
    # Expression filtering parameters
    min_gene_count: int = 100
    min_cells_per_pixel: int = 5
    percentile_clip: float = 0.95
    
    # Image processing parameters
    smooth_sigma: float = 2.0
    mask_count_threshold: int = 5
    
    # Similarity analysis parameters
    similarity_metric: str = "cosine"  # "cosine" or "linear"
    row_normalize: bool = False
    
    # Marker selection parameters
    n_markers: int = 20
    max_rank_fraction: float = -1.0
    n_neighbors: int = 10
    
    # Visualization parameters
    figsize: Tuple[int, int] = (5, 5)
    cmap: str = "viridis"
    
    @staticmethod
    def adaptive_bin_size(n_cells: int, target_cells_per_bin: int = 10) -> int:
        """
        Calculate adaptive bin size based on cell density.
        
        Logic: Use sqrt(n_cells / target_cells_per_bin) to maintain roughly 
        constant cell density per spatial bin across datasets.
        
        Args:
            n_cells: Total number of cells in dataset
            target_cells_per_bin: Desired average cells per spatial bin
            
        Returns:
            Recommended bin size
        """
        return int(np.sqrt(n_cells / target_cells_per_bin))
    
    @classmethod 
    def derive_parameters_from_data(cls, n_cells: int, species: str = None, 
                                   cell_density: str = "auto") -> 'SpatialAnalysisParams':
        """
        Derive parameters transparently based on data characteristics.
        
        This method shows the parameter selection logic that was originally
        in the notebook analysis.
        
        Args:
            n_cells: Number of cells in the dataset
            species: Species type for gene naming conventions
            cell_density: "high" (>50k cells), "medium" (10-50k), "low" (<10k), or "auto"
            
        Returns:
            SpatialAnalysisParams with explanatory logging
        """
        
        # Determine cell density category
        if cell_density == "auto":
            if n_cells > 50000:
                cell_density = "high"
            elif n_cells > 10000:
                cell_density = "medium"
            else:
                cell_density = "low"
        
        print(f"\n--- Parameter Derivation for {n_cells:,} cells ({cell_density} density) ---")
        
        # Adaptive bin size calculation
        suggested_bin_size = cls.adaptive_bin_size(n_cells, target_cells_per_bin=10)
        print(f"Adaptive bin size calculation: sqrt({n_cells}/10) = {suggested_bin_size}")
        
        # Apply density-based adjustments
        if cell_density == "high":
            # High cell density: can use larger bins, more smoothing
            bin_size = max(51, suggested_bin_size)
            smooth_sigma = 2.0
            min_cells_per_pixel = 5
            min_gene_count = 100
            percentile_clip = 0.95
            mask_count_threshold = 5
            print(f"High density adjustments: bin_size={bin_size}, smooth_sigma={smooth_sigma}")
            
        elif cell_density == "medium":
            # Medium density: moderate parameters
            bin_size = min(51, max(25, suggested_bin_size))
            smooth_sigma = 1.0
            min_cells_per_pixel = 3
            min_gene_count = 50
            percentile_clip = 0.93
            mask_count_threshold = 3
            print(f"Medium density adjustments: bin_size={bin_size}, smooth_sigma={smooth_sigma}")
            
        else:  # low density
            # Low cell density: smaller bins, less smoothing, permissive filtering
            bin_size = min(25, suggested_bin_size)
            smooth_sigma = 0.5
            min_cells_per_pixel = 1
            min_gene_count = 0  # Very permissive for sparse data
            percentile_clip = 0.90
            mask_count_threshold = 1
            print(f"Low density adjustments: bin_size={bin_size}, smooth_sigma={smooth_sigma}")
        
        # Species-specific gene naming adjustments
        species_note = ""
        if species and species.lower() == 'mouse':
            species_note = " (mouse gene naming: Sentence case)"
        elif species and species.lower() in ['human', 'chick']:
            species_note = f" ({species} gene naming: UPPERCASE)"
        
        print(f"Final parameters{species_note}:")
        print(f"  bin_size={bin_size} (spatial resolution)")
        print(f"  smooth_sigma={smooth_sigma} (Gaussian smoothing)")
        print(f"  min_cells_per_pixel={min_cells_per_pixel} (spatial filtering)")
        print(f"  min_gene_count={min_gene_count} (expression filtering)")
        print(f"  percentile_clip={percentile_clip} (outlier clipping)")
        
        return cls(
            bin_size=bin_size,
            smooth_sigma=smooth_sigma,
            min_cells_per_pixel=min_cells_per_pixel,
            min_gene_count=min_gene_count,
            percentile_clip=percentile_clip,
            mask_count_threshold=mask_count_threshold
        )
    
    @classmethod
    def get_species_defaults(cls, species: str, show_reasoning: bool = True) -> 'SpatialAnalysisParams':
        """
        Get default parameters optimized for specific species.
        
        Note: These are hard-coded defaults from the original analysis.
        For transparent parameter selection based on your data, use 
        derive_parameters_from_data() instead.
        
        Args:
            species: 'chick', 'human', or 'mouse'
            show_reasoning: Whether to print the reasoning behind parameter choices
        """
        defaults = {
            'chick': {
                'params': cls(
                    bin_size=51,
                    min_gene_count=100,
                    min_cells_per_pixel=5,
                    percentile_clip=0.95,
                    smooth_sigma=2.0,
                    mask_count_threshold=5
                ),
                'reasoning': [
                    "Chick retina dataset: ~32K cells (high density)",
                    "bin_size=51: Provides good spatial resolution for large dataset",
                    "smooth_sigma=2.0: Moderate smoothing to capture spatial patterns",
                    "min_gene_count=100: Standard filtering for well-powered dataset",
                    "percentile_clip=0.95: Conservative clipping to preserve dynamic range"
                ]
            },
            'human': {
                'params': cls(
                    bin_size=25,
                    min_gene_count=0,
                    min_cells_per_pixel=1,
                    percentile_clip=0.90,
                    smooth_sigma=0.5,
                    mask_count_threshold=2
                ),
                'reasoning': [
                    "Human fetal retina dataset: ~6.7K cells (low-medium density)",  
                    "bin_size=25: Smaller bins due to lower cell density",
                    "smooth_sigma=0.5: Minimal smoothing to preserve spatial detail",
                    "min_gene_count=0: Permissive filtering for sparse data",
                    "min_cells_per_pixel=1: Allow single-cell spatial bins",
                    "percentile_clip=0.90: More aggressive clipping for noisy data"
                ]
            },
            'mouse': {
                'params': cls(
                    bin_size=51,
                    min_gene_count=100,
                    min_cells_per_pixel=5,
                    percentile_clip=0.90,
                    smooth_sigma=2.0,
                    mask_count_threshold=5
                ),
                'reasoning': [
                    "Mouse retina dataset: ~61K cells (high density)",
                    "bin_size=51: Large bins suitable for high cell density",
                    "smooth_sigma=2.0: Moderate smoothing for pattern detection",
                    "min_gene_count=100: Standard filtering for well-powered dataset",
                    "Note: Mouse uses Sentence case gene naming (Fgf8, Tbx5, etc.)"
                ]
            }
        }
        
        species_lower = species.lower()
        if species_lower not in defaults:
            raise ValueError(f"Species '{species}' not recognized. Available: {list(defaults.keys())}")
        
        config = defaults[species_lower]
        
        if show_reasoning:
            print(f"\n--- {species.upper()} Default Parameters ---")
            for reason in config['reasoning']:
                print(f"  • {reason}")
            print()
        
        return config['params']
    
    def explain_parameters(self) -> None:
        """Print explanation of current parameter values"""
        print("Current Parameter Settings:")
        print(f"  bin_size={self.bin_size}")
        print(f"    → Creates {self.bin_size}×{self.bin_size} = {self.bin_size**2} spatial bins")
        print(f"  smooth_sigma={self.smooth_sigma}")
        print(f"    → Gaussian smoothing kernel width (0=no smoothing)")
        print(f"  min_gene_count={self.min_gene_count}")
        print(f"    → Exclude genes with < {self.min_gene_count} total counts")
        print(f"  min_cells_per_pixel={self.min_cells_per_pixel}")
        print(f"    → Mask spatial bins with < {self.min_cells_per_pixel} cells")
        print(f"  percentile_clip={self.percentile_clip}")
        print(f"    → Clip expression values above {self.percentile_clip*100}th percentile")

@dataclass
class NeighborInfo:
    """Information about neighboring genes"""
    name: str
    q_ij: float   # Q(i, j)
    q_ji: float   # Q(j, i)  
    r_ij: int     # rank of j in row-i (0-based)
    r_ji: int     # rank of i in row-j

@dataclass
class MarkerSetInfo:
    """Information about marker gene sets"""
    name: str
    neighbors: List[NeighborInfo]

class SpatialExpressionAnalyzer:
    """Main class for spatial expression analysis"""
    
    def __init__(self, params: SpatialAnalysisParams):
        self.params = params
        self.adata = None
        self.gene_names = None
        self.gene_info = None
        self.images = None
        self.similarity_matrix = None
        self.counts = None
        self.edges = None
        self.dv_mid = None
        self.nt_mid = None
        
    def load_data(self, adata_path: str) -> None:
        """Load AnnData object from file"""
        self.adata = ad.read_h5ad(adata_path)
        print(f"Loaded data: {self.adata.n_obs} cells × {self.adata.n_vars} genes")
        
    def adaptive_bin_size(self, n_cells: int, target_cells_per_bin: int = 10) -> int:
        """Calculate adaptive bin size based on cell density"""
        return int(np.sqrt(n_cells / target_cells_per_bin))
        
    def discretize_coordinates(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Discretize spatial coordinates into bins"""
        if self.adata is None:
            raise ValueError("Data not loaded. Call load_data() first.")
            
        # Get spatial coordinates
        dv = self.adata.obs["DV.Score"].to_numpy(float)
        nt = self.adata.obs["NT.Score"].to_numpy(float)
        
        # Determine bin size
        if self.params.adaptive_binning:
            bin_size = self.adaptive_bin_size(self.adata.n_obs, self.params.target_cells_per_bin)
            print(f"Using adaptive bin size: {bin_size}")
        else:
            bin_size = self.params.bin_size
            
        # Create bins
        eps = 1e-10
        dv_edges = np.linspace(dv.min() - eps, dv.max() + eps, bin_size + 1)
        nt_edges = np.linspace(nt.min() - eps, nt.max() + eps, bin_size + 1)
        
        # Digitize coordinates
        dv_idx = np.digitize(dv, dv_edges) - 1
        nt_idx = np.digitize(nt, nt_edges) - 1
        pix_idx = dv_idx * bin_size + nt_idx
        
        # Store for later use
        self.edges = {"DV": dv_edges, "NT": nt_edges}
        self.dv_mid = np.argmin(np.abs(dv_edges))
        self.nt_mid = np.argmin(np.abs(nt_edges))
        
        # Calculate cell counts per pixel
        n_pix = bin_size * bin_size
        pix_counts = np.bincount(pix_idx, minlength=n_pix).astype(np.float32)
        pix_counts[pix_counts == 0] = 1.0  # avoid division by zero
        self.counts = pix_counts.reshape(bin_size, bin_size)
        
        return pix_idx, pix_counts, bin_size
        
    def filter_genes(self, X: sp.csr_matrix) -> Tuple[sp.csr_matrix, np.ndarray, pd.DataFrame]:
        """Filter genes based on expression thresholds"""
        gene_tot = np.array(X.sum(axis=0)).flatten()
        
        # Apply minimum count filter
        gene_kept = np.where(gene_tot > self.params.min_gene_count)[0]
        X_filtered = X[:, gene_kept]
        gene_names = self.adata.var_names[gene_kept]
        
        # Create gene info dataframe
        gene_info = pd.DataFrame({
            "name": gene_names,
            "tot": gene_tot[gene_kept]
        })
        
        print(f"Filtered genes: {len(gene_names)} / {len(self.adata.var_names)} "
              f"(min_count >= {self.params.min_gene_count})")
        
        self.gene_names = gene_names
        self.gene_info = gene_info
        
        return X_filtered, gene_names, gene_info
        
    def create_spatial_images(self) -> np.ndarray:
        """Create spatial expression images"""
        if self.adata is None:
            raise ValueError("Data not loaded. Call load_data() first.")
            
        # Discretize coordinates
        pix_idx, pix_counts, bin_size = self.discretize_coordinates()
        
        # Build design matrix
        N = self.adata.n_obs
        design = sp.coo_matrix((
            np.ones(N, dtype=np.float32),
            (np.arange(N, dtype=np.int32), pix_idx)),
            shape=(N, bin_size * bin_size)
        ).tocsr()
        
        # Extract and filter expression matrix
        X = self.adata.X
        if not sp.issparse(X):
            X = sp.csr_matrix(X, dtype=np.float32)
        else:
            X = X.astype(np.float32)
            
        X_filtered, gene_names, gene_info = self.filter_genes(X)
        
        # Aggregate expression by pixels
        pixel_matrix = X_filtered.T @ design
        pixel_matrix = pixel_matrix.multiply(1.0 / pix_counts).tocsr()
        
        # Convert to dense images with clipping
        G = pixel_matrix.shape[0]
        images = np.empty((G, bin_size, bin_size), dtype=np.float32)
        
        chunk = 1024
        for start in range(0, G, chunk):
            end = min(start + chunk, G)
            block = pixel_matrix[start:end].toarray()
            
            # Apply clipping
            nz_mask = (block > 0) & (pix_counts > self.params.min_cells_per_pixel)
            for i in range(block.shape[0]):
                nz_vals = block[i, nz_mask[i]]
                if nz_vals.size:
                    thresh = np.quantile(nz_vals, self.params.percentile_clip)
                    block[i, block[i] > thresh] = thresh
                    
            images[start:end] = block.reshape(-1, bin_size, bin_size)
            
        # Apply smoothing
        if self.params.smooth_sigma > 0:
            for i in range(images.shape[0]):
                images[i] = gaussian_filter(images[i], self.params.smooth_sigma)
                
        self.images = images
        print(f"Created spatial images: {images.shape}")
        
        return images
        
    def calculate_similarity_matrix(self) -> np.ndarray:
        """Calculate gene similarity matrix"""
        if self.images is None:
            raise ValueError("Images not created. Call create_spatial_images() first.")
            
        # Create mask for valid pixels
        G, H, W = self.images.shape
        P = H * W
        mask = self.counts.ravel() > self.params.min_cells_per_pixel
        n_pix = mask.sum()
        
        # Extract expression matrix for valid pixels
        mtx = self.images.reshape(G, P)[:, mask]
        
        # Calculate similarity
        if self.params.similarity_metric == "cosine":
            sim_mtx = sklearn.metrics.pairwise.cosine_similarity(mtx)
        elif self.params.similarity_metric == "linear":
            sim_mtx = sklearn.metrics.pairwise.linear_kernel(
                normalize(mtx, axis=1, norm="l2")
            )
        else:
            raise ValueError(f"Unknown similarity metric: {self.params.similarity_metric}")
            
        self.similarity_matrix = sim_mtx
        print(f"Calculated similarity matrix: {sim_mtx.shape}")
        
        return sim_mtx
        
    def show_gene_image(self, gene_name: str, show_parameters: bool = True, **kwargs) -> None:
        """Display spatial expression image for a gene"""
        if self.images is None or self.gene_names is None:
            raise ValueError("Images not created. Call create_spatial_images() first.")
            
        if gene_name not in self.gene_names:
            raise ValueError(f"Gene '{gene_name}' not found in dataset")
            
        # Get gene index
        gene2idx = {gene: i for i, gene in enumerate(self.gene_names)}
        idx = gene2idx[gene_name]
        
        # Get spatial max intensity for this gene
        spatial_max = float(np.max(self.images[idx]))
        
        # Set default parameters
        plot_params = {
            'smooth_sigma': 0,  # Already smoothed
            'counts': self.counts,
            'mask_count': self.params.mask_count_threshold,
            'figsize': self.params.figsize,
            'cmap': self.params.cmap,
            'mask_zero': 1e-6,
            'log1p': False,
            'vmin': None,
            'vmax': None,
            'title': gene_name,
            'show_parameters': show_parameters,
            'spatial_max': spatial_max
        }
        plot_params.update(kwargs)
        
        self._plot_gene_image(gene_name, self.images[idx], **plot_params)
        
    def _plot_gene_image(self, gene: str, img: np.ndarray, 
                        smooth_sigma: float = 0,
                        counts: Optional[np.ndarray] = None,
                        mask_count: int = 2,
                        figsize: Tuple[int, int] = (5, 5),
                        cmap: str = "viridis",
                        mask_zero: float = 1e-6,
                        log1p: bool = False,
                        vmin: Optional[float] = None,
                        vmax: Optional[float] = None,
                        title: Optional[str] = None,
                        show_parameters: bool = False,
                        spatial_max: float = 0.0,
                        **imshow_kwargs) -> None:
        """Internal function to plot gene expression images"""
        
        if smooth_sigma > 0:
            img = gaussian_filter(img, smooth_sigma)
            
        if log1p:
            img = np.log1p(img)
            
        # Create mask
        if counts is not None and mask_count > 0:
            mask = (counts < mask_count)
        else:
            mask = np.zeros_like(img, dtype=bool)
            
        img = np.ma.masked_where(mask, img)
        w, h = img.shape
        
        fig, ax = plt.subplots(figsize=figsize)
        im = ax.imshow(
            img,
            origin="lower",
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            **imshow_kwargs
        )
        
        # Add reference lines
        if self.dv_mid is not None and self.nt_mid is not None:
            ax.axhline(y=self.dv_mid, color="black", lw=0.5)
            ax.axvline(x=self.nt_mid, color="black", lw=0.5)
            
        ax.set_xlim(0, w)
        ax.set_ylim(0, h)
        
        # Enhanced title with spatial max intensity
        if show_parameters and spatial_max > 0:
            enhanced_title = f"{title or gene} (max: {spatial_max:.2f})"
            ax.set_title(enhanced_title, fontweight="bold")
        else:
            ax.set_title(title or gene, fontweight="bold")
            
        ax.set_xlabel("Temporal  ← NT →  Nasal")
        ax.set_ylabel("Ventral  ← DV →  Dorsal")
        ax.set_aspect("equal")
        
        # Add parameter caption at bottom right
        if show_parameters:
            caption_text = (f"percentile={self.params.percentile_clip:.2f}, "
                          f"σ={self.params.smooth_sigma}")
            ax.text(0.98, 0.02, caption_text, 
                   transform=ax.transAxes,
                   fontsize=8, ha='right', va='bottom',
                   bbox=dict(boxstyle='round,pad=0.3', 
                            facecolor='white', alpha=0.8))
        
        cbar = fig.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label("Relative expression", rotation=270, labelpad=15)
        plt.tight_layout()
        plt.show()
        
    def get_gene_correlations(self, gene_name: str, top_n: int = 10) -> List[Tuple[str, float]]:
        """Get top correlated genes for a given gene"""
        if self.similarity_matrix is None:
            raise ValueError("Similarity matrix not calculated. Call calculate_similarity_matrix() first.")
            
        if gene_name not in self.gene_names:
            raise ValueError(f"Gene '{gene_name}' not found in dataset")
            
        gene2idx = {gene: i for i, gene in enumerate(self.gene_names)}
        idx = gene2idx[gene_name]
        
        # Get similarities and sort
        similarities = self.similarity_matrix[idx, :]
        top_indices = np.argsort(-similarities)[:top_n]
        
        correlations = [(self.gene_names[i], similarities[i]) for i in top_indices]
        return correlations
        
    def run_full_analysis(self, adata_path: str) -> Dict:
        """Run complete spatial analysis pipeline"""
        print("=== Starting Spatial Expression Analysis ===")
        
        # Load data
        self.load_data(adata_path)
        
        # Create spatial images
        images = self.create_spatial_images()
        
        # Calculate similarity matrix
        similarity_matrix = self.calculate_similarity_matrix()
        
        results = {
            'images': images,
            'similarity_matrix': similarity_matrix,
            'gene_names': self.gene_names,
            'gene_info': self.gene_info,
            'counts': self.counts,
            'analyzer': self
        }
        
        print("=== Analysis Complete ===")
        return results

    def show_gene_correlations_panel(self, gene_name: str, top_n: int = 8, 
                                   ncols: int = 3, figsize: Tuple[int, int] = (12, 8),
                                   cmap: str = "viridis", show_colorbar: bool = False,
                                   correlation_threshold: float = 0.0,
                                   show_parameters: bool = True) -> None:
        """
        Create a compact multi-panel figure showing a gene's expression pattern 
        alongside its top correlated genes.
        
        Args:
            gene_name: Name of the main gene to display
            top_n: Number of top correlated genes to show (including the main gene)
            ncols: Number of columns in the subplot grid
            figsize: Figure size (width, height)
            cmap: Colormap for the expression patterns
            show_colorbar: Whether to show individual colorbars for each subplot
            correlation_threshold: Minimum correlation to include genes
            show_parameters: Whether to show processing parameters as caption
            
        Layout:
            - Top left: Main gene expression pattern
            - Other panels: Top correlated genes ordered by correlation strength
            - Titles show: gene name, correlation (r), and max spatial intensity
            - Bottom right caption shows: percentile cutoff, smoothing, overall max intensity
        """
        if self.images is None or self.gene_names is None:
            raise ValueError("Images not created. Call create_spatial_images() first.")
            
        if gene_name not in self.gene_names:
            raise ValueError(f"Gene '{gene_name}' not found in dataset")
        
        # Get correlations for the main gene
        correlations = self.get_gene_correlations(gene_name, top_n=top_n)
        
        # Filter by correlation threshold
        filtered_correlations = [(g, c) for g, c in correlations if c >= correlation_threshold]
        
        if len(filtered_correlations) == 0:
            print(f"No genes found with correlation >= {correlation_threshold}")
            return
            
        # Limit to requested number
        genes_to_plot = filtered_correlations[:top_n]
        n_genes = len(genes_to_plot)
        
        # Calculate grid dimensions
        nrows = int(np.ceil(n_genes / ncols))
        
        # Create subplot grid
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
        if nrows == 1 and ncols == 1:
            axes = [axes]
        elif nrows == 1 or ncols == 1:
            axes = axes.flatten()
        else:
            axes = axes.flatten()
        
        # Gene to index mapping
        gene2idx = {gene: i for i, gene in enumerate(self.gene_names)}
        
        # Calculate spatial image max values for each gene
        spatial_max_values = {}
        for gene, _ in genes_to_plot:
            idx = gene2idx[gene]
            # Get max value from the processed spatial image
            spatial_max_values[gene] = float(np.max(self.images[idx]))
        
        # Calculate overall max across displayed genes for context
        overall_max = max(spatial_max_values.values())
        
        # Plot each gene
        for i, (gene, correlation) in enumerate(genes_to_plot):
            if i >= len(axes):
                break
                
            ax = axes[i]
            idx = gene2idx[gene]
            img = self.images[idx].copy()
            
            # Apply mask if counts are available
            if self.counts is not None:
                mask = (self.counts < self.params.mask_count_threshold)
                img = np.ma.masked_where(mask, img)
            
            # Plot the image
            im = ax.imshow(
                img,
                origin="lower",
                cmap=cmap,
                aspect="equal"
            )
            
            # Add reference lines if available
            if self.dv_mid is not None and self.nt_mid is not None:
                ax.axhline(y=self.dv_mid, color="black", lw=0.5, alpha=0.7)
                ax.axvline(x=self.nt_mid, color="black", lw=0.5, alpha=0.7)
            
            # Get spatial max value for this gene
            spatial_max = spatial_max_values[gene]
            
            # Set title with correlation and spatial max intensity
            if i == 0:  # Main gene
                title = f"{gene}\n(main gene, max: {spatial_max:.2f})"
                ax.set_title(title, fontweight="bold", fontsize=10)
            else:
                title = f"{gene}\n(r = {correlation:.3f}, max: {spatial_max:.2f})"
                ax.set_title(title, fontsize=9)
            
            # Remove axis labels for cleaner look
            ax.set_xticks([])
            ax.set_yticks([])
            
            # Add colorbar if requested
            if show_colorbar:
                cbar = plt.colorbar(im, ax=ax, shrink=0.6, aspect=10)
                cbar.ax.tick_params(labelsize=8)
        
        # Hide unused subplots
        for i in range(n_genes, len(axes)):
            axes[i].set_visible(False)
        
        # Add overall title
        fig.suptitle(f"Expression Patterns: {gene_name} and Top Correlated Genes", 
                    fontsize=14, fontweight="bold", y=1.0)
        
        # Add axis labels only to bottom-left subplot
        if len(axes) > 0:
            axes[0].set_xlabel("Temporal ← NT → Nasal", fontsize=8)
            axes[0].set_ylabel("Ventral ← DV → Dorsal", fontsize=8)
        
        # Add parameter caption at bottom right
        if show_parameters and len(axes) > 0:
            caption_text = (f"Parameters: percentile={self.params.percentile_clip:.2f}, "
                          f"σ={self.params.smooth_sigma}, max_intensity={overall_max:.2f}")
            
            # Place caption on the last visible subplot (bottom right)
            last_visible_idx = min(n_genes - 1, len(axes) - 1)
            if last_visible_idx >= 0:
                axes[last_visible_idx].text(0.98, 0.02, caption_text, 
                                          transform=axes[last_visible_idx].transAxes,
                                          fontsize=8, ha='right', va='bottom',
                                          bbox=dict(boxstyle='round,pad=0.3', 
                                                   facecolor='white', alpha=0.8))
        
        plt.tight_layout()
        plt.subplots_adjust(top=0.93)  # Make room for suptitle
        plt.show()
        
        # Print correlation and spatial intensity summary
        print(f"\nCorrelation & Spatial Intensity Summary for {gene_name}:")
        print("-" * 55)
        print(f"Processing: percentile_clip={self.params.percentile_clip:.2f}, "
              f"smooth_sigma={self.params.smooth_sigma}, max_intensity={overall_max:.2f}")
        print("-" * 55)
        for i, (gene, correlation) in enumerate(genes_to_plot):
            spatial_max = spatial_max_values[gene]
            if i == 0:
                print(f"  Main: {gene} (max: {spatial_max:.2f})")
            else:
                print(f"  {i:2d}. {gene}: r={correlation:.3f}, max={spatial_max:.2f}")

    def show_gene_correlations_compact(self, gene_name: str, top_n: int = 5, 
                                      figsize: Tuple[int, int] = (15, 3),
                                      cmap: str = "viridis", 
                                      show_parameters: bool = True) -> None:
        """
        Create a horizontal compact view of gene correlations.
        
        Args:
            gene_name: Name of the main gene to display
            top_n: Number of genes to show (including main gene)
            figsize: Figure size (width, height)
            cmap: Colormap for the expression patterns
            show_parameters: Whether to show processing parameters as caption
            
        Titles show: gene name, correlation (r), and max spatial intensity
        Caption shows: percentile cutoff, smoothing, max intensity
        """
        if self.images is None or self.gene_names is None:
            raise ValueError("Images not created. Call create_spatial_images() first.")
            
        if gene_name not in self.gene_names:
            raise ValueError(f"Gene '{gene_name}' not found in dataset")
        
        # Get correlations
        correlations = self.get_gene_correlations(gene_name, top_n=top_n)
        genes_to_plot = correlations[:top_n]
        
        # Create horizontal subplot grid
        fig, axes = plt.subplots(1, len(genes_to_plot), figsize=figsize)
        if len(genes_to_plot) == 1:
            axes = [axes]
        
        # Gene to index mapping
        gene2idx = {gene: i for i, gene in enumerate(self.gene_names)}
        
        # Calculate spatial image max values for each gene
        spatial_max_values = {}
        for gene, _ in genes_to_plot:
            idx = gene2idx[gene]
            # Get max value from the processed spatial image
            spatial_max_values[gene] = float(np.max(self.images[idx]))
        
        # Calculate overall max across displayed genes
        overall_max = max(spatial_max_values.values())
        
        # Plot each gene
        for i, (gene, correlation) in enumerate(genes_to_plot):
            ax = axes[i]
            idx = gene2idx[gene]
            img = self.images[idx].copy()
            
            # Apply mask if counts are available
            if self.counts is not None:
                mask = (self.counts < self.params.mask_count_threshold)
                img = np.ma.masked_where(mask, img)
            
            # Plot the image
            im = ax.imshow(
                img,
                origin="lower",
                cmap=cmap,
                aspect="equal"
            )
            
            # Add reference lines if available
            if self.dv_mid is not None and self.nt_mid is not None:
                ax.axhline(y=self.dv_mid, color="black", lw=0.5, alpha=0.7)
                ax.axvline(x=self.nt_mid, color="black", lw=0.5, alpha=0.7)
            
            # Get spatial max value for this gene
            spatial_max = spatial_max_values[gene]
            
            # Set title with correlation and spatial max intensity
            if i == 0:
                title = f"{gene}\n(main, {spatial_max:.2f})"
                ax.set_title(title, fontweight="bold", fontsize=10)
            else:
                title = f"{gene}\n(r={correlation:.3f}, {spatial_max:.2f})"
                ax.set_title(title, fontsize=9)
            
            # Remove axis ticks
            ax.set_xticks([])
            ax.set_yticks([])
            
            # Add labels only to first subplot
            if i == 0:
                ax.set_xlabel("T ← NT → N", fontsize=8)
                ax.set_ylabel("V ← DV → D", fontsize=8)
        
        # Add parameter caption at bottom right of last panel
        if show_parameters and len(axes) > 0:
            caption_text = (f"percentile={self.params.percentile_clip:.2f}, "
                          f"σ={self.params.smooth_sigma}, max_intensity={overall_max:.2f}")
            axes[-1].text(0.98, 0.02, caption_text, 
                         transform=axes[-1].transAxes,
                         fontsize=8, ha='right', va='bottom',
                         bbox=dict(boxstyle='round,pad=0.3', 
                                  facecolor='white', alpha=0.8))
        
        plt.suptitle(f"Gene Expression Correlations: {gene_name}", fontsize=12, fontweight="bold", y=1.0)
        plt.tight_layout()
        plt.subplots_adjust(top=0.85)
        plt.show()

    # ========================================================================
    # METADATA VISUALIZATION METHODS
    # ========================================================================

    def create_metadata_spatial_image(self, column_name: str,
                                      categorical: bool = False,
                                      compute_variance: bool = False) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray], Dict[str, np.ndarray]]:
        """
        Create spatial images from metadata columns.

        Args:
            column_name: Name of the column in adata.obs to visualize
            categorical: If True, treat as categorical and return dict of category proportions
            compute_variance: If True (and not categorical), also compute variance per bin

        Returns:
            - If categorical: Dict mapping category names to spatial proportion images
            - If continuous and compute_variance=False: Single mean image
            - If continuous and compute_variance=True: Tuple of (mean_image, variance_image)
        """
        if self.adata is None:
            raise ValueError("Data not loaded. Call load_data() first.")

        if column_name not in self.adata.obs.columns:
            raise ValueError(f"Column '{column_name}' not found in adata.obs")

        # Get discretization info (reuse from gene analysis)
        if not hasattr(self, 'edges') or self.edges is None:
            pix_idx, pix_counts, bin_size = self.discretize_coordinates()
        else:
            # Recompute to ensure consistency
            pix_idx, pix_counts, bin_size = self.discretize_coordinates()

        metadata_values = self.adata.obs[column_name]

        if categorical:
            # Handle categorical metadata
            return self._create_categorical_images(metadata_values, pix_idx, bin_size, pix_counts)
        else:
            # Handle continuous metadata
            return self._create_continuous_images(metadata_values, pix_idx, bin_size, pix_counts,
                                                  compute_variance=compute_variance)

    def _create_continuous_images(self, values: pd.Series, pix_idx: np.ndarray,
                                  bin_size: int, pix_counts: np.ndarray,
                                  compute_variance: bool = False) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
        """Create mean (and optionally variance) spatial images from continuous metadata."""
        n_pix = bin_size * bin_size

        # Convert to numeric, handling NaN
        numeric_values = pd.to_numeric(values, errors='coerce').to_numpy(dtype=np.float32)
        valid_mask = ~np.isnan(numeric_values)

        # Calculate mean per bin
        mean_per_bin = np.zeros(n_pix, dtype=np.float32)
        for i in range(n_pix):
            bin_mask = (pix_idx == i) & valid_mask
            if bin_mask.sum() > 0:
                mean_per_bin[i] = numeric_values[bin_mask].mean()
            else:
                mean_per_bin[i] = np.nan

        mean_image = mean_per_bin.reshape(bin_size, bin_size)

        # Apply smoothing
        if self.params.smooth_sigma > 0:
            # Handle NaN values in smoothing
            mask = ~np.isnan(mean_image)
            mean_image_filled = mean_image.copy()
            mean_image_filled[~mask] = 0
            mean_image = gaussian_filter(mean_image_filled, self.params.smooth_sigma)
            mean_image[~mask] = np.nan

        if not compute_variance:
            return mean_image

        # Calculate variance per bin
        var_per_bin = np.zeros(n_pix, dtype=np.float32)
        for i in range(n_pix):
            bin_mask = (pix_idx == i) & valid_mask
            if bin_mask.sum() > 1:  # Need at least 2 values for variance
                var_per_bin[i] = numeric_values[bin_mask].var()
            else:
                var_per_bin[i] = np.nan

        var_image = var_per_bin.reshape(bin_size, bin_size)

        # Apply smoothing to variance
        if self.params.smooth_sigma > 0:
            mask = ~np.isnan(var_image)
            var_image_filled = var_image.copy()
            var_image_filled[~mask] = 0
            var_image = gaussian_filter(var_image_filled, self.params.smooth_sigma)
            var_image[~mask] = np.nan

        return mean_image, var_image

    def _create_categorical_images(self, values: pd.Series, pix_idx: np.ndarray,
                                   bin_size: int, pix_counts: np.ndarray) -> Dict[str, np.ndarray]:
        """Create proportion images for each category in categorical metadata."""
        categories = values.dropna().unique()
        n_pix = bin_size * bin_size

        category_images = {}

        for category in categories:
            # Binary indicator: 1 if this category, 0 otherwise
            is_category = (values == category).to_numpy(dtype=np.float32)

            # Calculate proportion per bin
            proportion_per_bin = np.zeros(n_pix, dtype=np.float32)
            for i in range(n_pix):
                bin_mask = (pix_idx == i)
                n_cells_in_bin = bin_mask.sum()
                if n_cells_in_bin > 0:
                    proportion_per_bin[i] = is_category[bin_mask].mean()
                else:
                    proportion_per_bin[i] = np.nan

            proportion_image = proportion_per_bin.reshape(bin_size, bin_size)

            # Optional: apply minimal smoothing for categorical (usually you want sharp boundaries)
            # We'll skip smoothing for categorical by default to preserve category boundaries

            category_images[str(category)] = proportion_image

        return category_images

    def show_metadata_image(self, column_name: str, show_variance: bool = False,
                           variance_cmap: str = "plasma", show_parameters: bool = True,
                           figsize: Optional[Tuple[int, int]] = None, **kwargs) -> None:
        """
        Display spatial pattern of continuous metadata column.

        Args:
            column_name: Name of the continuous metadata column in adata.obs
            show_variance: If True, show both mean and variance as two panels
            variance_cmap: Colormap for variance panel (default: 'plasma')
            show_parameters: Whether to show processing parameters
            figsize: Figure size (auto-determined if None)
            **kwargs: Additional arguments passed to imshow
        """
        if self.adata is None:
            raise ValueError("Data not loaded. Call load_data() first.")

        if column_name not in self.adata.obs.columns:
            raise ValueError(f"Column '{column_name}' not found in adata.obs")

        # Check if column is numeric
        if not pd.api.types.is_numeric_dtype(self.adata.obs[column_name]):
            raise ValueError(f"Column '{column_name}' is not numeric. Use show_metadata_categorical() for categorical data.")

        # Create images
        if show_variance:
            mean_image, var_image = self.create_metadata_spatial_image(column_name,
                                                                        categorical=False,
                                                                        compute_variance=True)
            self._plot_mean_variance_panels(column_name, mean_image, var_image,
                                           variance_cmap=variance_cmap,
                                           show_parameters=show_parameters,
                                           figsize=figsize, **kwargs)
        else:
            mean_image = self.create_metadata_spatial_image(column_name,
                                                            categorical=False,
                                                            compute_variance=False)
            self._plot_single_metadata_image(column_name, mean_image,
                                            show_parameters=show_parameters,
                                            figsize=figsize, **kwargs)

    def _plot_single_metadata_image(self, column_name: str, image: np.ndarray,
                                    show_parameters: bool = True,
                                    figsize: Optional[Tuple[int, int]] = None,
                                    cmap: str = "viridis", **kwargs) -> None:
        """Plot a single metadata image (mean values)."""
        if figsize is None:
            figsize = self.params.figsize

        # Calculate statistics
        valid_vals = image[~np.isnan(image)]
        if len(valid_vals) > 0:
            spatial_max = float(np.nanmax(image))
            spatial_min = float(np.nanmin(image))
            spatial_mean = float(np.nanmean(image))
        else:
            spatial_max = spatial_min = spatial_mean = 0.0

        fig, ax = plt.subplots(figsize=figsize)

        # Mask invalid values
        masked_image = np.ma.masked_where(np.isnan(image), image)

        im = ax.imshow(
            masked_image,
            origin="lower",
            cmap=cmap,
            aspect='equal',
            **kwargs
        )

        # Add reference lines
        if self.dv_mid is not None and self.nt_mid is not None:
            ax.axhline(y=self.dv_mid, color="black", lw=0.5)
            ax.axvline(x=self.nt_mid, color="black", lw=0.5)

        # Title with statistics
        if show_parameters:
            title = f"{column_name}\n(mean: {spatial_mean:.2f}, max: {spatial_max:.2f}, min: {spatial_min:.2f})"
        else:
            title = column_name

        ax.set_title(title, fontweight="bold")
        ax.set_xlabel("Temporal  ← NT →  Nasal")
        ax.set_ylabel("Ventral  ← DV →  Dorsal")

        # Add parameter caption
        if show_parameters:
            caption_text = f"σ={self.params.smooth_sigma}, bin_size={self.params.bin_size}"
            ax.text(0.98, 0.02, caption_text,
                   transform=ax.transAxes,
                   fontsize=8, ha='right', va='bottom',
                   bbox=dict(boxstyle='round,pad=0.3',
                            facecolor='white', alpha=0.8))

        cbar = fig.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label(column_name, rotation=270, labelpad=15)
        plt.tight_layout()
        plt.show()

    def _plot_mean_variance_panels(self, column_name: str, mean_image: np.ndarray,
                                   var_image: np.ndarray, variance_cmap: str = "plasma",
                                   show_parameters: bool = True,
                                   figsize: Optional[Tuple[int, int]] = None, **kwargs) -> None:
        """Plot mean and variance as two side-by-side panels."""
        if figsize is None:
            figsize = (12, 5)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

        # Mask invalid values
        masked_mean = np.ma.masked_where(np.isnan(mean_image), mean_image)
        masked_var = np.ma.masked_where(np.isnan(var_image), var_image)

        # Calculate statistics
        mean_max = float(np.nanmax(mean_image)) if np.any(~np.isnan(mean_image)) else 0.0
        mean_avg = float(np.nanmean(mean_image)) if np.any(~np.isnan(mean_image)) else 0.0
        var_max = float(np.nanmax(var_image)) if np.any(~np.isnan(var_image)) else 0.0

        # Plot mean
        im1 = ax1.imshow(masked_mean, origin="lower", cmap="viridis", aspect='equal', **kwargs)
        if self.dv_mid is not None and self.nt_mid is not None:
            ax1.axhline(y=self.dv_mid, color="black", lw=0.5)
            ax1.axvline(x=self.nt_mid, color="black", lw=0.5)
        ax1.set_title(f"Mean: {column_name}\n(avg: {mean_avg:.2f}, max: {mean_max:.2f})", fontweight="bold")
        ax1.set_xlabel("Temporal  ← NT →  Nasal")
        ax1.set_ylabel("Ventral  ← DV →  Dorsal")
        cbar1 = fig.colorbar(im1, ax=ax1, shrink=0.8)
        cbar1.set_label("Mean", rotation=270, labelpad=15)

        # Plot variance
        im2 = ax2.imshow(masked_var, origin="lower", cmap=variance_cmap, aspect='equal', **kwargs)
        if self.dv_mid is not None and self.nt_mid is not None:
            ax2.axhline(y=self.dv_mid, color="black", lw=0.5)
            ax2.axvline(x=self.nt_mid, color="black", lw=0.5)
        ax2.set_title(f"Variance: {column_name}\n(max: {var_max:.2f})", fontweight="bold")
        ax2.set_xlabel("Temporal  ← NT →  Nasal")
        ax2.set_ylabel("Ventral  ← DV →  Dorsal")
        cbar2 = fig.colorbar(im2, ax=ax2, shrink=0.8)
        cbar2.set_label("Variance", rotation=270, labelpad=15)

        # Add parameter caption
        if show_parameters:
            caption_text = f"σ={self.params.smooth_sigma}, bin_size={self.params.bin_size}"
            ax2.text(0.98, 0.02, caption_text,
                   transform=ax2.transAxes,
                   fontsize=8, ha='right', va='bottom',
                   bbox=dict(boxstyle='round,pad=0.3',
                            facecolor='white', alpha=0.8))

        plt.tight_layout()
        plt.show()

    def show_metadata_categorical(self, column_name: str, ncol: int = 3,
                                  figsize: Optional[Tuple[int, int]] = None,
                                  cmap: str = "viridis", show_parameters: bool = True,
                                  vmin: float = 0.0, vmax: Optional[float] = None,
                                  percentile_clip: Optional[float] = None) -> None:
        """
        Display spatial composition of categorical metadata as faceted panels.

        Each panel shows the proportion of cells in each bin belonging to that category.

        Args:
            column_name: Name of the categorical metadata column in adata.obs
            ncol: Number of columns in the facet grid
            figsize: Figure size (auto-determined if None)
            cmap: Colormap for proportion heatmaps (default: 'viridis')
            show_parameters: Whether to show processing parameters
            vmin: Minimum value for colormap (default: 0.0 = 0%)
            vmax: Maximum value for colormap (default: None, auto-determined)
                  If None and percentile_clip is None, uses 1.0 (100%)
                  If None and percentile_clip is set, uses percentile value
            percentile_clip: Percentile for clipping max values (e.g., 0.95 for 95th percentile)
                           If None, uses self.params.percentile_clip
                           Set to 1.0 to disable clipping
        """
        if self.adata is None:
            raise ValueError("Data not loaded. Call load_data() first.")

        if column_name not in self.adata.obs.columns:
            raise ValueError(f"Column '{column_name}' not found in adata.obs")

        # Create categorical images
        category_images = self.create_metadata_spatial_image(column_name, categorical=True)

        categories = list(category_images.keys())
        n_categories = len(categories)

        if n_categories == 0:
            print(f"No categories found in column '{column_name}'")
            return

        # Determine vmax using percentile clipping
        if vmax is None:
            # Use percentile clipping if specified
            if percentile_clip is None:
                percentile_clip = self.params.percentile_clip

            if percentile_clip < 1.0:
                # Collect all non-NaN values across all categories
                all_values = []
                for img in category_images.values():
                    valid_vals = img[~np.isnan(img)]
                    if len(valid_vals) > 0:
                        all_values.extend(valid_vals)

                if len(all_values) > 0:
                    vmax = float(np.quantile(all_values, percentile_clip))
                    print(f"Using percentile clipping: {percentile_clip:.2%} → vmax={vmax:.2%}")
                else:
                    vmax = 1.0
            else:
                vmax = 1.0

        # Calculate grid dimensions
        nrow = int(np.ceil(n_categories / ncol))

        # Determine figure size
        if figsize is None:
            figsize = (ncol * 4, nrow * 3.5)

        fig, axes = plt.subplots(nrow, ncol, figsize=figsize)

        # Flatten axes for easier iteration
        if n_categories == 1:
            axes = [axes]
        else:
            axes = axes.flatten() if nrow > 1 else axes

        # Plot each category
        for i, category in enumerate(categories):
            ax = axes[i]
            image = category_images[category]

            # Mask invalid values
            masked_image = np.ma.masked_where(np.isnan(image), image)

            # Calculate statistics for this category
            max_prop = float(np.nanmax(image)) if np.any(~np.isnan(image)) else 0.0
            mean_prop = float(np.nanmean(image)) if np.any(~np.isnan(image)) else 0.0

            # Plot
            im = ax.imshow(masked_image, origin="lower", cmap=cmap, aspect='equal',
                          vmin=vmin, vmax=vmax)

            # Add reference lines
            if self.dv_mid is not None and self.nt_mid is not None:
                ax.axhline(y=self.dv_mid, color="white", lw=0.5, alpha=0.5)
                ax.axvline(x=self.nt_mid, color="white", lw=0.5, alpha=0.5)

            # Title with statistics
            ax.set_title(f"{category}\n(mean: {mean_prop:.2%}, max: {max_prop:.2%})",
                        fontsize=10, fontweight="bold")
            ax.set_xticks([])
            ax.set_yticks([])

            # Add labels only to leftmost and bottom subplots
            if i % ncol == 0:  # Leftmost
                ax.set_ylabel("V ← DV → D", fontsize=8)
            if i >= n_categories - ncol:  # Bottom row
                ax.set_xlabel("T ← NT → N", fontsize=8)

            # Colorbar for each panel
            cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            cbar.set_label("Proportion", rotation=270, labelpad=12, fontsize=8)

        # Hide unused subplots
        for i in range(n_categories, len(axes)):
            axes[i].set_visible(False)

        # Add parameter caption to figure
        if show_parameters:
            caption_text = f"{column_name}: bin_size={self.params.bin_size}, percentile_clip={percentile_clip:.2f}, vmax={vmax:.2%}"
            fig.text(0.5, 0.01, caption_text, ha='center', fontsize=8, style='italic', color='gray')

        # Use tight_layout with rect to leave room for caption
        plt.tight_layout(rect=[0, 0.03, 1, 1.0])
        plt.show()

    def show_cell_density(self, figsize: Optional[Tuple[int, int]] = None,
                          cmap: str = "viridis", show_parameters: bool = True,
                          log_scale: bool = False, **kwargs) -> None:
        """
        Display spatial cell density (number of cells per bin).

        Args:
            figsize: Figure size (auto-determined if None)
            cmap: Colormap for cell density (default: 'viridis')
            show_parameters: Whether to show processing parameters
            log_scale: If True, use log10 scale for cell counts (useful for large ranges)
            **kwargs: Additional arguments passed to imshow
        """
        if self.counts is None:
            raise ValueError("Cell counts not calculated. Call load_data() or discretize_coordinates() first.")

        if figsize is None:
            figsize = self.params.figsize

        # Get cell counts
        density = self.counts.copy()

        # Apply mask for low-count bins (same as gene expression visualization)
        mask = (self.counts < self.params.mask_count_threshold)

        # Apply log scale if requested
        if log_scale:
            # Add 1 to avoid log(0)
            density = np.log10(density + 1)
            scale_label = "log10(cells + 1)"
        else:
            scale_label = "cells"

        # Mask low-count bins
        density = np.ma.masked_where(mask, density)

        # Calculate statistics (only for valid bins)
        valid_counts = self.counts[~mask]
        total_cells = int(np.sum(self.counts))
        max_density = int(np.max(self.counts))
        mean_density = float(np.mean(valid_counts)) if len(valid_counts) > 0 else 0.0
        median_density = float(np.median(valid_counts)) if len(valid_counts) > 0 else 0.0
        bins_with_cells = int(np.sum(~mask))
        total_bins = self.counts.size

        fig, ax = plt.subplots(figsize=figsize)

        im = ax.imshow(
            density,
            origin="lower",
            cmap=cmap,
            aspect='equal',
            **kwargs
        )

        # Add reference lines
        if self.dv_mid is not None and self.nt_mid is not None:
            ax.axhline(y=self.dv_mid, color="black", lw=0.5)
            ax.axvline(x=self.nt_mid, color="black", lw=0.5)

        # Title with statistics
        if show_parameters:
            title = (f"Cell Density\n"
                    f"Total: {total_cells:,} cells, Max: {max_density} cells/bin, "
                    f"Mean: {mean_density:.1f}, Median: {median_density:.1f}")
        else:
            title = "Cell Density"

        ax.set_title(title, fontweight="bold")
        ax.set_xlabel("Temporal  ← NT →  Nasal")
        ax.set_ylabel("Ventral  ← DV →  Dorsal")

        # Add parameter caption
        if show_parameters:
            caption_text = (f"bin_size={self.params.bin_size}, "
                          f"bins with cells: {bins_with_cells}/{total_bins} ({100*bins_with_cells/total_bins:.1f}%)")
            ax.text(0.98, 0.02, caption_text,
                   transform=ax.transAxes,
                   fontsize=8, ha='right', va='bottom',
                   bbox=dict(boxstyle='round,pad=0.3',
                            facecolor='white', alpha=0.8))

        cbar = fig.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label(f"Cells per bin ({scale_label})", rotation=270, labelpad=15)
        plt.tight_layout()
        plt.show()

    def show_gene_image_with_variance(self, gene_name: str, variance_cmap: str = "plasma",
                                     show_parameters: bool = True,
                                     figsize: Optional[Tuple[int, int]] = None, **kwargs) -> None:
        """
        Display gene expression with both mean and variance panels.

        Args:
            gene_name: Name of the gene to display
            variance_cmap: Colormap for variance panel (default: 'plasma')
            show_parameters: Whether to show processing parameters
            figsize: Figure size (auto-determined if None)
            **kwargs: Additional arguments passed to imshow
        """
        if self.images is None or self.gene_names is None:
            raise ValueError("Images not created. Call create_spatial_images() first.")

        if gene_name not in self.gene_names:
            raise ValueError(f"Gene '{gene_name}' not found in dataset")

        # Get gene index and image
        gene2idx = {gene: i for i, gene in enumerate(self.gene_names)}
        idx = gene2idx[gene_name]
        mean_image = self.images[idx].copy()

        # Calculate variance image by binning raw expression values
        # Get raw expression for this gene
        if gene_name in self.adata.var_names:
            gene_idx_full = list(self.adata.var_names).index(gene_name)
            expression = self.adata.X[:, gene_idx_full]
            if sp.issparse(expression):
                expression = expression.toarray().flatten()
            else:
                expression = expression.flatten()

            # Get pixel indices
            pix_idx, pix_counts, bin_size = self.discretize_coordinates()

            # Calculate variance per bin
            n_pix = bin_size * bin_size
            var_per_bin = np.zeros(n_pix, dtype=np.float32)
            for i in range(n_pix):
                bin_mask = (pix_idx == i)
                if bin_mask.sum() > 1:
                    var_per_bin[i] = expression[bin_mask].var()
                else:
                    var_per_bin[i] = np.nan

            var_image = var_per_bin.reshape(bin_size, bin_size)

            # Apply smoothing to variance
            if self.params.smooth_sigma > 0:
                mask = ~np.isnan(var_image)
                var_image_filled = var_image.copy()
                var_image_filled[~mask] = 0
                var_image = gaussian_filter(var_image_filled, self.params.smooth_sigma)
                var_image[~mask] = np.nan

            # Plot both panels
            self._plot_mean_variance_panels(gene_name, mean_image, var_image,
                                           variance_cmap=variance_cmap,
                                           show_parameters=show_parameters,
                                           figsize=figsize, **kwargs)
        else:
            raise ValueError(f"Gene '{gene_name}' not found in original adata.var_names")

class MarkerSelectorPy:
    """Gene marker selection and clustering functionality"""
    
    def __init__(self, gene_names: Sequence[str], Q: np.ndarray,
                 valid_mask: Optional[np.ndarray] = None,
                 verbose: int = 0, rownorm: bool = True):
        self.gene_names = list(gene_names)
        self.name2idx = {g: i for i, g in enumerate(gene_names)}
        self.Q = normalize(np.asarray(Q), axis=1, norm="l1") if rownorm else np.asarray(Q)
        self.M = self.Q.shape[0]
        self.verbose = verbose
        
        if valid_mask is None:
            valid_mask = np.ones(self.M, dtype=bool)
        self.valid_mask = np.asarray(valid_mask, dtype=bool)
        
        self._rank_pos: Optional[np.ndarray] = None
        self.anchors_idx: List[int] = []
        
    def select_markers(self, K: int, fixed_anchors: Optional[Sequence[str]] = None,
                      eps: float = 1e-8) -> List[str]:
        """Select diverse marker genes using greedy selection"""
        if K > self.M:
            raise ValueError(f"K={K} exceeds number of features M={self.M}")
            
        fixed_anchors = fixed_anchors or []
        anchors: List[int] = []
        
        # Add fixed anchors
        fixed_set: set[int] = set()
        for name in fixed_anchors:
            idx = self.name2idx.get(name, None)
            if idx is None:
                if self.verbose:
                    print(f"[select_markers] fixed anchor {name} not found")
                continue
            if idx not in fixed_set:
                anchors.append(idx)
                fixed_set.add(idx)
                
        if len(fixed_set) > 0 and self.verbose:
            print(f"[select_markers] read {len(fixed_set)} fixed anchors")
            
        K0 = len(anchors)
        
        # Build orthonormal basis
        basis: List[np.ndarray] = []
        for idx in anchors:
            v = self.Q[idx].copy()
            for b in basis:
                v -= np.dot(b, v) * b
            if np.linalg.norm(v) > eps:
                basis.append(v / np.linalg.norm(v))
                
        if len(basis) > K:
            raise ValueError("More than K fixed anchors were provided")
            
        # Greedy selection
        while len(anchors) < K:
            best_idx, best_d2 = -1, -1.0
            for i in range(self.M):
                if not self.valid_mask[i] or i in fixed_set:
                    continue
                r = self.Q[i]
                proj = sum(np.dot(b, r) * b for b in basis) if basis else 0.0
                d2 = np.sum((r - proj) ** 2)
                if d2 > best_d2:
                    best_d2, best_idx = d2, i
                    
            if best_idx < 0:
                raise RuntimeError("Could not find a further valid anchor")
                
            anchors.append(best_idx)
            fixed_set.add(best_idx)
            
            v = self.Q[best_idx].copy()
            for b in basis:
                v -= np.dot(b, v) * b
            if np.linalg.norm(v) > eps:
                basis.append(v / np.linalg.norm(v))
            elif self.verbose:
                print(f"[select_markers] matrix may have rank < {K}")
                
            if self.verbose:
                k = len(anchors) - 1
                print(f"[select_markers] selected anchor {k}: {self.gene_names[best_idx]}")
                
        # Refinement step
        for t in range(K0, K):
            other_basis: List[np.ndarray] = []
            other_set: set[int] = {anchors[j] for j in range(K) if j != t}
            for jidx in other_set:
                v = self.Q[jidx].copy()
                for b in other_basis:
                    v -= np.dot(b, v) * b
                if np.linalg.norm(v) > eps:
                    other_basis.append(v / np.linalg.norm(v))
                else:
                    other_basis.append(np.zeros_like(v))
                    
            best_idx, best_d2 = anchors[t], -1.0
            for i in range(self.M):
                if not self.valid_mask[i] or i in other_set:
                    continue
                qi = self.Q[i]
                proj = sum(np.dot(b, qi) * b for b in other_basis) if other_basis else 0.0
                d2 = np.sum((qi - proj) ** 2)
                if d2 > best_d2:
                    best_d2, best_idx = d2, i
                    
            if best_idx != anchors[t] and self.verbose:
                print(f"[select_markers] replaced anchor {t}: "
                      f"{self.gene_names[anchors[t]]} -> {self.gene_names[best_idx]}")
            anchors[t] = best_idx
            
        self.anchors_idx = anchors
        return [self.gene_names[i] for i in anchors]
        
    def find_neighbors_to_anchors(self, m: int, max_rank_fraction: float = -1.0, 
                                 gene_to_label: Optional[Dict[str, str]] = None) -> List[MarkerSetInfo]:
        """Find neighboring genes for each anchor"""
        if not self.anchors_idx:
            raise RuntimeError("Call select_markers() first")
            
        # Build rank matrix
        if self._rank_pos is None:
            rank_pos = np.empty_like(self.Q, dtype=np.int32)
            for i in range(self.M):
                order = np.argsort(-self.Q[i])
                rank_pos[i, order] = np.arange(self.M)
                if self.verbose and i % 1000 == 0:
                    print(f"[find_neighbors] ranked row {i}")
            self._rank_pos = rank_pos
            if self.verbose:
                print("[find_neighbors] rank matrix built")
                
        rank_pos = self._rank_pos
        max_rank = self.M if max_rank_fraction <= 0 else int(max_rank_fraction * self.M)
        
        results: List[MarkerSetInfo] = []
        
        for k, anchor_idx in enumerate(self.anchors_idx):
            scores = np.maximum(rank_pos[anchor_idx], rank_pos[:, anchor_idx])
            scores[anchor_idx] = self.M
            
            cand_idx = np.argpartition(scores, m)[:m]
            cand_scores = scores[cand_idx]
            
            mask = cand_scores <= max_rank
            cand_idx = cand_idx[mask]
            cand_scores = cand_scores[mask]
            
            order = np.argsort(cand_scores)
            cand_idx = cand_idx[order]
            
            neighbors: List[NeighborInfo] = []
            for j in cand_idx:
                neighbors.append(NeighborInfo(
                    name=self.gene_names[j],
                    q_ij=float(self.Q[anchor_idx, j]),
                    q_ji=float(self.Q[j, anchor_idx]),
                    r_ij=int(rank_pos[anchor_idx, j]),
                    r_ji=int(rank_pos[j, anchor_idx]),
                ))
                
            if self.verbose:
                anchor_name = self.gene_names[anchor_idx]
                # Add label if available
                display_name = anchor_name
                if gene_to_label and anchor_name in gene_to_label:
                    display_name = f"{anchor_name} ({gene_to_label[anchor_name]})"
                    
                neigh_str = " ".join(f"{n.name}({n.r_ij}/{n.r_ji})" for n in neighbors)
                print(f"[find_neighbors] {k} {display_name} -> {neigh_str}")
                
            results.append(MarkerSetInfo(
                name=self.gene_names[anchor_idx],
                neighbors=neighbors,
            ))
            
        return results

def cluster_by_anchor(selector: MarkerSelectorPy, 
                     fixed_anchors: Optional[List[str]] = None,
                     null_q: Optional[np.ndarray] = None,
                     gene_to_label: Optional[Dict[str, str]] = None,
                     verbose: int = 0) -> Tuple[np.ndarray, Dict[str, List[str]], np.ndarray, Dict]:
    """
    Cluster genes by similarity to anchor genes
    
    Returns:
        cluster_idx: Array indicating which cluster each gene belongs to
        clusters: Dictionary mapping anchor names to their gene lists  
        D2: Similarity matrix between genes and anchors
        anchor_mapping: Dictionary with anchor indexing information for interrogation
    """
    
    if fixed_anchors is not None:
        fixed_anchors_idx = [selector.name2idx.get(g, None) for g in fixed_anchors]
        selector.anchors_idx = [a for a in fixed_anchors_idx if a is not None]
        if len(selector.anchors_idx) == 0:
            raise ValueError("No valid anchors were provided")
            
    if not selector.anchors_idx:
        raise RuntimeError("Call select_markers() first")
        
    K = len(selector.anchors_idx)
    M = selector.M
    Q = selector.Q
    
    # Create anchor matrix
    A = Q[selector.anchors_idx]
    if null_q is None:
        A = np.vstack((np.ones((1, M), dtype=Q.dtype) / M, A))
    else:
        A = np.vstack((null_q.reshape((1, -1)), A))
        
    # Calculate similarities
    D2 = sklearn.metrics.pairwise.cosine_similarity(Q, A)
    cluster_idx = np.argmax(D2, axis=1)
    
    # Build clusters dictionary
    clusters: Dict[str, List[str]] = {selector.gene_names[a]: [] for a in selector.anchors_idx}
    clusters["NULL"] = []
    
    for i, k in enumerate(cluster_idx):
        if k == 0:
            clusters["NULL"].append(selector.gene_names[i])
        else:
            anchor_name = selector.gene_names[selector.anchors_idx[k-1]]
            clusters[anchor_name].append(selector.gene_names[i])
    
    # Create anchor mapping for interrogation
    anchor_names = [selector.gene_names[a] for a in selector.anchors_idx]
    anchor_mapping = {
        'filtered_to_name': {k: anchor_names[k] for k in range(len(anchor_names))},
        'name_to_filtered': {name: k for k, name in enumerate(anchor_names)},
        'filtered_anchors': anchor_names,
        'original_anchors': fixed_anchors if fixed_anchors else anchor_names,
        'anchor_to_category': {}
    }
    
    # Add category information if available
    if gene_to_label:
        for name in anchor_names:
            if name in gene_to_label:
                anchor_mapping['anchor_to_category'][name] = gene_to_label[name]
            
    if verbose:
        print(f"[cluster_by_anchor] null cluster: {len(clusters['NULL'])} genes")
        for k, a in enumerate(selector.anchors_idx):
            anchor_name = selector.gene_names[a]
            # Add label if available
            display_name = anchor_name
            if gene_to_label and anchor_name in gene_to_label:
                display_name = f"{anchor_name} ({gene_to_label[anchor_name]})"
                
            print(f"[cluster_by_anchor] anchor {k} {display_name} → "
                  f"{len(clusters[anchor_name])} genes")
            
    return cluster_idx, clusters, D2, anchor_mapping

def create_cluster_summary(anchor_mapping: Dict, clusters: Dict[str, List[str]], 
                          original_markers: Optional[List[str]] = None) -> 'pd.DataFrame':
    """
    Create a summary table for cluster interrogation
    
    Args:
        anchor_mapping: Anchor mapping dictionary from cluster_by_anchor
        clusters: Clusters dictionary from cluster_by_anchor
        original_markers: Optional list of original markers for reference
        
    Returns:
        DataFrame with anchor information and cluster sizes
    """
    import pandas as pd
    
    summary_data = []
    
    for filtered_idx, anchor_name in anchor_mapping['filtered_to_name'].items():
        row = {
            'filtered_index': filtered_idx,
            'anchor_name': anchor_name,
            'cluster_size': len(clusters[anchor_name]),
            'category': anchor_mapping['anchor_to_category'].get(anchor_name, 'Unknown')
        }
        
        # Add original position if available
        if original_markers and anchor_name in original_markers:
            row['original_index'] = original_markers.index(anchor_name)
        else:
            row['original_index'] = None
            
        summary_data.append(row)
    
    # Add NULL cluster
    summary_data.append({
        'filtered_index': -1,
        'anchor_name': 'NULL',
        'cluster_size': len(clusters['NULL']),
        'category': 'Background',
        'original_index': None
    })
    
    return pd.DataFrame(summary_data)

def interrogate_cluster(cluster_idx: int, anchor_mapping: Dict, clusters: Dict[str, List[str]], 
                       analyzer: Optional['SpatialExpressionAnalyzer'] = None, top_genes: int = 10) -> Dict:
    """
    Get detailed information about a specific cluster
    
    Args:
        cluster_idx: The filtered cluster index (from cluster_by_anchor output)
        anchor_mapping: Anchor mapping from cluster_by_anchor
        clusters: Clusters dictionary from cluster_by_anchor
        analyzer: Optional SpatialExpressionAnalyzer for expression info
        top_genes: Number of top genes to return
        
    Returns:
        Dictionary with cluster information
    """
    if cluster_idx == -1:
        anchor_name = 'NULL'
        category = 'Background'
    else:
        if cluster_idx not in anchor_mapping['filtered_to_name']:
            raise ValueError(f"Cluster index {cluster_idx} not found. Available: {list(anchor_mapping['filtered_to_name'].keys())}")
        
        anchor_name = anchor_mapping['filtered_to_name'][cluster_idx]
        category = anchor_mapping['anchor_to_category'].get(anchor_name, 'Unknown')
    
    cluster_genes = clusters[anchor_name]
    
    result = {
        'filtered_index': cluster_idx,
        'anchor_name': anchor_name,
        'category': category,
        'cluster_size': len(cluster_genes),
        'genes': cluster_genes[:top_genes] if len(cluster_genes) > top_genes else cluster_genes
    }
    
    # Add expression information if analyzer provided
    if analyzer and analyzer.gene_info is not None:
        gene_info = analyzer.gene_info
        cluster_gene_info = gene_info[gene_info['name'].isin(cluster_genes)]
        if not cluster_gene_info.empty:
            top_expressed = cluster_gene_info.nlargest(top_genes, 'tot')
            result['top_expressed_genes'] = top_expressed[['name', 'tot']].to_dict('records')
    
    return result

def get_fixed_anchors(gene_info: 'pd.DataFrame', species: str) -> Dict[str, str]:
    """
    Get fixed anchor genes for a species with their spatial category labels.
    
    To customize the gene options, edit the control_genes.yaml file.
    
    Args:
        gene_info: DataFrame with 'name' and 'tot' columns
        species: Species name ('chick', 'human', 'mouse')
        
    Returns:
        Dictionary mapping category labels to selected gene names
        e.g., {'DV.stripes': 'BMP2', 'Fovea': 'FGF8', ...}
    """
    control_genes = CONTROL_GENES.get(species.lower(), {})
    anchors_with_labels = {}
    
    for category, genes in control_genes.items():
        available_genes = gene_info.loc[gene_info["name"].isin(genes), :]
        if not available_genes.empty:
            best_gene = available_genes.sort_values(by="tot", ascending=False).name.iloc[0]
            anchors_with_labels[category] = best_gene
            print(f"  {category}: {best_gene} (from options: {genes})")
        else:
            print(f"Warning: No genes found for {species} {category} category: {genes}")
            
    return anchors_with_labels

def reload_control_genes(config_file: str = "control_genes.yaml"):
    """
    Reload control genes from YAML file. 
    Call this after editing control_genes.yaml to pick up changes.
    """
    # If config_file is just a filename, look for it in the script directory
    if not os.path.isabs(config_file) and os.path.dirname(config_file) == "":
        script_dir = os.path.dirname(os.path.abspath(__file__))
        config_file = os.path.join(script_dir, config_file)
    
    config_file_abs = os.path.abspath(config_file)
    
    global CONTROL_GENES
    CONTROL_GENES = load_control_genes(config_file)
    print(f"Available species: {list(CONTROL_GENES.keys())}")
    return CONTROL_GENES

if __name__ == "__main__":
    print("Spatial Expression Analysis Module")
    print("==================================")
    print("Available classes:")
    print("- SpatialExpressionAnalyzer: Main analysis class")
    print("- SpatialAnalysisParams: Parameter configuration")
    print("- MarkerSelectorPy: Gene marker selection")
    print("Available functions:")
    print("- cluster_by_anchor: Gene clustering by anchors")
    print("- create_cluster_summary: Create summary table for clusters")
    print("- interrogate_cluster: Get detailed cluster information")
    print("- get_fixed_anchors: Get species-specific anchor genes")
    print("- reload_control_genes: Reload anchor gene configuration") 