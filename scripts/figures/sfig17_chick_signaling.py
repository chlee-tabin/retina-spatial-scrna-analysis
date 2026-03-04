# ---
# jupyter:
#   jupytext:
#     formats: py:percent
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Figure S17: Fgf, Bmp, and RA Signaling Pathway Maps in Chicken Retina

# %% [markdown]
# ## Utility: save_correlation_panels()

# %% tags=["cell-34"]
import sys
sys.path.append('..')
from spatial_expression_analysis import SpatialAnalysisParams, SpatialExpressionAnalyzer
import matplotlib.pyplot as plt
import numpy as np
import os
def save_correlation_panels(gene, analyzer, panel_top_n=12, compact_top_n=5, 
                         panel_ncols=4, output_dir="figures/correlation_panels"):
  """
  Save correlation panels for a gene with customizable parameters.
  Args:
      gene: Gene name to analyze
      analyzer: The SpatialExpressionAnalyzer object
      panel_top_n: Number of genes to show in panel view (default 12)
      compact_top_n: Number of genes to show in compact view (default 5)
      panel_ncols: Number of columns in panel view (default 4)
      output_dir: Directory to save figures
  """
  if gene not in analyzer.gene_names:
      print(f"Gene {gene} not found in dataset")
      return
  print(f"\n--- Saving Correlation Panels for {gene} ---")
  # Get correlations
  correlations = analyzer.get_gene_correlations(gene, top_n=max(panel_top_n, compact_top_n))
  # Gene to index mapping
  gene2idx = {g: i for i, g in enumerate(analyzer.gene_names)}
  # Calculate spatial max values
  spatial_max_values = {}
  for g, _ in correlations[:max(panel_top_n, compact_top_n)]:
      idx = gene2idx[g]
      spatial_max_values[g] = float(np.max(analyzer.images[idx]))
  # === PANEL VIEW (Grid) ===
  if panel_top_n > 0:
      # Calculate number of rows needed
      panel_nrows = (panel_top_n + panel_ncols - 1) // panel_ncols
      # Adjust figure size based on grid dimensions
      fig_width = 3 * panel_ncols
      fig_height = 3 * panel_nrows
      # Create figure with tighter spacing
      fig, axes = plt.subplots(panel_nrows, panel_ncols, figsize=(fig_width, fig_height))
      axes = axes.flatten() if panel_nrows * panel_ncols > 1 else [axes]
      # Plot each gene
      for i, (g, correlation) in enumerate(correlations[:panel_top_n]):
          if i >= len(axes):
              break
          ax = axes[i]
          idx = gene2idx[g]
          img = analyzer.images[idx].copy()
          # Apply mask if counts are available
          if analyzer.counts is not None:
              mask = (analyzer.counts < analyzer.params.mask_count_threshold)
              img = np.ma.masked_where(mask, img)
          # Plot the image
          im = ax.imshow(
              img,
              origin="lower",
              cmap="viridis",
              aspect="equal"
          )
          # Add reference lines
          if analyzer.dv_mid is not None and analyzer.nt_mid is not None:
              ax.axhline(y=analyzer.dv_mid, color="black", lw=0.5, alpha=0.7)
              ax.axvline(x=analyzer.nt_mid, color="black", lw=0.5, alpha=0.7)
          # Get spatial max value
          spatial_max = spatial_max_values[g]
          # Set title with correlation and spatial max
          if i == 0:  # Main gene
              title = f"{g}\n(main, max: {spatial_max:.2f})"
              ax.set_title(title, fontweight="bold", fontsize=10)
          else:
              title = f"{g}\n(r={correlation:.3f}, max:{spatial_max:.2f})"
              ax.set_title(title, fontsize=9)
          # Remove axis labels
          ax.set_xticks([])
          ax.set_yticks([])
      # Hide unused subplots
      for i in range(panel_top_n, len(axes)):
          axes[i].set_visible(False)
      # Add overall title
      fig.suptitle(f"Expression Patterns: {gene} and Top {panel_top_n-1} Correlated Genes",
                  fontsize=14, fontweight="bold", y=0.98)
      # Use tight_layout with reduced padding for tighter spacing
      plt.tight_layout(pad=0.5, h_pad=0.5, w_pad=0.5)
      plt.subplots_adjust(top=0.93)  # Adjust for suptitle
      # Save the panel figure
      panel_path = os.path.join(output_dir, f"{gene}_correlation_panel_top{panel_top_n}.png")
      plt.savefig(panel_path, dpi=300, bbox_inches='tight')
      plt.close()
      print(f"  Saved panel view (top {panel_top_n}) to {panel_path}")
  # === COMPACT VIEW (Horizontal) ===
  if compact_top_n > 0:
      # Get top correlations for compact view
      correlations_compact = correlations[:compact_top_n]
      # Adjust figure width based on number of genes
      fig_width = 3 * compact_top_n
      # Create horizontal subplot grid
      fig, axes = plt.subplots(1, compact_top_n, figsize=(fig_width, 4))
      if compact_top_n == 1:
          axes = [axes]
      # Plot each gene
      for i, (g, correlation) in enumerate(correlations_compact):
          ax = axes[i]
          idx = gene2idx[g]
          img = analyzer.images[idx].copy()
          # Apply mask if counts are available
          if analyzer.counts is not None:
              mask = (analyzer.counts < analyzer.params.mask_count_threshold)
              img = np.ma.masked_where(mask, img)
          # Plot the image
          im = ax.imshow(
              img,
              origin="lower",
              cmap="viridis",
              aspect="equal"
          )
          # Add reference lines
          if analyzer.dv_mid is not None and analyzer.nt_mid is not None:
              ax.axhline(y=analyzer.dv_mid, color="black", lw=0.5, alpha=0.7)
              ax.axvline(x=analyzer.nt_mid, color="black", lw=0.5, alpha=0.7)
          # Get spatial max value
          spatial_max = spatial_max_values[g]
          # Set title
          if i == 0:  # Main gene
              title = f"{g}\n(main gene)\nmax: {spatial_max:.2f}"
              ax.set_title(title, fontweight="bold", fontsize=10)
          else:
              title = f"{g}\nr = {correlation:.3f}\nmax: {spatial_max:.2f}"
              ax.set_title(title, fontsize=9)
          # Add axis labels to first subplot only
          if i == 0:
              ax.set_xlabel("Temporal ← NT → Nasal", fontsize=8)
              ax.set_ylabel("Ventral ← DV → Dorsal", fontsize=8)
          # Remove axis ticks
          ax.set_xticks([])
          ax.set_yticks([])
      # Add overall title
      plt.suptitle(f"Gene Expression Correlations: {gene} (Top {compact_top_n})",
                  fontsize=12, fontweight="bold", y=1.05)
      # Use tight_layout with reduced padding
      plt.tight_layout(pad=0.5, w_pad=0.3)
      plt.subplots_adjust(top=0.85)
      # Save the compact figure
      compact_path = os.path.join(output_dir, f"{gene}_correlation_compact_top{compact_top_n}.png")
      plt.savefig(compact_path, dpi=300, bbox_inches='tight')
      plt.close()
      print(f"  Saved compact view (top {compact_top_n}) to {compact_path}")

# %% [markdown]
# ## SF17A-D: Correlation panels for signaling pathway genes

# %% tags=["cell-35"]
# Create output directory if it doesn't exist
FIGURES_BASE = os.path.join(os.path.dirname(__file__), "..", "..", "figures")
output_dir = os.path.join(FIGURES_BASE, "Figure_SF17")
os.makedirs(output_dir, exist_ok=True)
# Example usage with different parameters
for gene in ["FGF8", ]:
  # # Default: 12 genes in panel (3x4), 5 genes in compact
  # save_correlation_panels(gene, chick_analyzer)
  # Custom: 20 genes in panel (4x5), 8 genes in compact
  save_correlation_panels(
      gene, chick_analyzer, panel_top_n=20, compact_top_n=8, panel_ncols=5,
      output_dir = output_dir
  )
  # # Only panel view with 16 genes (4x4)
  # save_correlation_panels(gene, chick_analyzer, panel_top_n=16, compact_top_n=0, panel_ncols=4)
  # # Only compact view with 10 genes
  # save_correlation_panels(gene, chick_analyzer, panel_top_n=0, compact_top_n=10)
print(f"\nAll correlation panels saved to {output_dir}/")

# %% tags=["cell-36"]
# Create output directory if it doesn't exist
output_dir = os.path.join(FIGURES_BASE, "Figure_SF17")
os.makedirs(output_dir, exist_ok=True)
# Example usage with different parameters
for gene in ["CYP26C1", ]:
  # # Default: 12 genes in panel (3x4), 5 genes in compact
  # save_correlation_panels(gene, chick_analyzer)
  # Custom: 20 genes in panel (4x5), 8 genes in compact
  save_correlation_panels(
      gene, chick_analyzer, panel_top_n=20, compact_top_n=8, panel_ncols=5,
      output_dir = output_dir
  )
  # # Only panel view with 16 genes (4x4)
  # save_correlation_panels(gene, chick_analyzer, panel_top_n=16, compact_top_n=0, panel_ncols=4)
  # # Only compact view with 10 genes
  # save_correlation_panels(gene, chick_analyzer, panel_top_n=0, compact_top_n=10)
print(f"\nAll correlation panels saved to {output_dir}/")

# %% tags=["cell-37"]
# Create output directory if it doesn't exist
output_dir = os.path.join(FIGURES_BASE, "Figure_SF17")
os.makedirs(output_dir, exist_ok=True)
# Example usage with different parameters
for gene in ["BMP2", ]:
  # # Default: 12 genes in panel (3x4), 5 genes in compact
  # save_correlation_panels(gene, chick_analyzer)
  # Custom: 20 genes in panel (4x5), 8 genes in compact
  save_correlation_panels(
      gene, chick_analyzer, panel_top_n=20, compact_top_n=8, panel_ncols=5,
      output_dir = output_dir
  )
  # # Only panel view with 16 genes (4x4)
  # save_correlation_panels(gene, chick_analyzer, panel_top_n=16, compact_top_n=0, panel_ncols=4)
  # # Only compact view with 10 genes
  # save_correlation_panels(gene, chick_analyzer, panel_top_n=0, compact_top_n=10)
print(f"\nAll correlation panels saved to {output_dir}/")

# %% tags=["cell-38"]
# Create output directory if it doesn't exist
output_dir = os.path.join(FIGURES_BASE, "Figure_SF17")
os.makedirs(output_dir, exist_ok=True)
# Example usage with different parameters
for gene in ["CYP1B1", ]:
  # # Default: 12 genes in panel (3x4), 5 genes in compact
  # save_correlation_panels(gene, chick_analyzer)
  # Custom: 20 genes in panel (4x5), 8 genes in compact
  save_correlation_panels(
      gene, chick_analyzer, panel_top_n=20, compact_top_n=8, panel_ncols=5,
      output_dir = output_dir
  )
  # # Only panel view with 16 genes (4x4)
  # save_correlation_panels(gene, chick_analyzer, panel_top_n=16, compact_top_n=0, panel_ncols=4)
  # # Only compact view with 10 genes
  # save_correlation_panels(gene, chick_analyzer, panel_top_n=0, compact_top_n=10)
print(f"\nAll correlation panels saved to {output_dir}/")

# %% tags=["cell-50"]
import matplotlib.pyplot as plt
for gene in ["BMP2", "CYP26C1"]:
  if gene in chick_analyzer.gene_names:
      print(f"\n--- Correlation Panel for {gene} ---")
      # Multi-panel grid view (main gene + top correlations)
      # Create figure instead of using show method
      correlations = chick_analyzer.get_gene_correlations(gene, top_n=12)
      # Create the panel manually to save it
      fig, axes = plt.subplots(3, 4, figsize=(12, 9))
      axes = axes.flatten()
      # Gene to index mapping
      gene2idx = {g: i for i, g in enumerate(chick_analyzer.gene_names)}
      # Calculate spatial max values
      spatial_max_values = {}
      for g, _ in correlations[:12]:
          idx = gene2idx[g]
          spatial_max_values[g] = float(np.max(chick_analyzer.images[idx]))
      # Plot each gene
      for i, (g, correlation) in enumerate(correlations[:12]):
          if i >= len(axes):
              break
          ax = axes[i]
          idx = gene2idx[g]
          img = chick_analyzer.images[idx].copy()
          # Apply mask if counts are available
          if chick_analyzer.counts is not None:
              mask = (chick_analyzer.counts < chick_analyzer.params.mask_count_threshold)
              img = np.ma.masked_where(mask, img)
          # Plot the image
          im = ax.imshow(
              img,
              origin="lower",
              cmap="viridis",
              aspect="equal"
          )
          # Add reference lines
          if chick_analyzer.dv_mid is not None and chick_analyzer.nt_mid is not None:
              ax.axhline(y=chick_analyzer.dv_mid, color="black", lw=0.5, alpha=0.7)
              ax.axvline(x=chick_analyzer.nt_mid, color="black", lw=0.5, alpha=0.7)
          # Get spatial max value
          spatial_max = spatial_max_values[g]
          # Set title with correlation and spatial max
          if i == 0:  # Main gene
              title = f"{g}\n(main gene, max: {spatial_max:.2f})"
              ax.set_title(title, fontweight="bold", fontsize=10)
          else:
              title = f"{g}\n(r = {correlation:.3f}, max: {spatial_max:.2f})"
              ax.set_title(title, fontsize=9)
          # Remove axis labels
          ax.set_xticks([])
          ax.set_yticks([])
      # Hide unused subplots
      for i in range(12, len(axes)):
          axes[i].set_visible(False)
      # Add overall title
      fig.suptitle(f"Expression Patterns: {gene} and Top Correlated Genes",
                  fontsize=14, fontweight="bold", y=1.0)
      # Add axis labels to first subplot
      axes[0].set_xlabel("Temporal ← NT → Nasal", fontsize=8)
      axes[0].set_ylabel("Ventral ← DV → Dorsal", fontsize=8)
      # Add parameters caption
      overall_max = max(spatial_max_values.values())
      caption_text = (f"Parameters: percentile={chick_analyzer.params.percentile_clip:.2f}, "
                    f"σ={chick_analyzer.params.smooth_sigma}, max_intensity={overall_max:.2f}")
      # Place caption on last visible subplot
      axes[11].text(0.98, 0.02, caption_text,
                   transform=axes[11].transAxes,
                   fontsize=8, ha='right', va='bottom',
                   bbox=dict(boxstyle='round,pad=0.3',
                            facecolor='white', alpha=0.8))
      plt.tight_layout()
      plt.subplots_adjust(top=0.93)
      # Save the figure
      plt.savefig(os.path.join(FIGURES_BASE, "Figure_SF17", f"SF17_{gene}_correlations.png"),
                  dpi=300, bbox_inches='tight')
      plt.close()
      print(f"Saved correlation panel for {gene}")
