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
# # Figure S19: Spatial Expression Clusters in Mouse and Human Retinas

# %% [markdown]
# ## Load human RPC data

# %% tags=["cell-74"]
import os
import sys
sys.path.append('..')

FIGURES_BASE = os.path.join(os.path.dirname(__file__), "..", "..", "figures")
from spatial_expression_analysis import (
    SpatialAnalysisParams, SpatialExpressionAnalyzer,
    reload_control_genes, get_fixed_anchors, MarkerSelectorPy,
    cluster_by_anchor, create_cluster_summary, CONTROL_GENES
)
import anndata as ad
human_adata = ad.read_h5ad("../data/20250604_human_RPC.h5ad")
print(f"Loaded human data: {human_adata.n_obs:,} cells × {human_adata.n_vars:,} genes")
# Derive parameters transparently based on the data characteristics
human_params = SpatialAnalysisParams.derive_parameters_from_data(
    n_cells=human_adata.n_obs, 
    species='human'
)
human_params

# %% [markdown]
# ## Create human spatial analyzer

# %% tags=["cell-78"]
# Initialize analyzer
human_analyzer = SpatialExpressionAnalyzer(human_params)
human_results = human_analyzer.run_full_analysis("../data/20250604_human_RPC.h5ad")

# %% [markdown]
# ## Human anchor selection

# %% [markdown]
# ## Selection of anchors
# With some fixed spatial patterns that we know of, we can maximize by choosing anchors by greedy approach. The anchors are selected based on maximum expression.

# %% tags=["cell-87"]
# To customize anchors, edit `control_genes.yaml` and call reload_control_genes()
reload_control_genes()

# %% tags=["cell-88"]
human_fixed_anchors_dict = get_fixed_anchors(human_analyzer.gene_info, 'human')
print("Fixed anchors from control_genes.yaml:")
for category, gene in human_fixed_anchors_dict.items():
    expression = human_analyzer.gene_info[human_analyzer.gene_info['name'] == gene]['tot'].iloc[0]
    print(f"  {category}: {gene} (expression: {expression:.1f})")
# Extract gene names for marker selection
human_fixed_anchors = list(human_fixed_anchors_dict.values())
print(f"\nUsing {len(human_fixed_anchors)} anchor genes for marker selection")
# Set up marker selector (replaces your MarkerSelectorPy setup)
valid_mask = human_analyzer.gene_info['tot'] >= 1000
human_selector = MarkerSelectorPy(
    gene_names=human_analyzer.gene_names,
    Q=human_analyzer.similarity_matrix,
    valid_mask=valid_mask,
    verbose=1,
    rownorm=False
)
# Select markers (replaces your select_markers call)
human_markers = human_selector.select_markers(K=20, fixed_anchors=human_fixed_anchors)
print("Selected markers:")
for i, marker in enumerate(human_markers, 1):
    expr = human_analyzer.gene_info[human_analyzer.gene_info['name'] == marker]['tot'].iloc[0]
    print(f"{i:2d}. {marker} ({expr:.1f})")

# %% tags=["cell-91"]
# First, create the gene-to-label mapping from control genes
human_gene_to_label = {}
human_control_genes = CONTROL_GENES.get('human', {})
for category, genes in human_control_genes.items():
    for gene in genes:
        if gene in human_analyzer.gene_names:
            human_gene_to_label[gene] = category

# %% tags=["cell-93"]
human_neighbors = human_selector.find_neighbors_to_anchors(
    m=20, # maximum candidate neighbors ("shortlist" size)
    max_rank_fraction=50/len(human_analyzer.gene_names), # acceptable level of rank quality. Only top 50 as rank
    gene_to_label=human_gene_to_label
)
# Keep only markers with neighbors
human_kept_markers = [nbr.name for nbr in human_neighbors if len(nbr.neighbors) > 0]
print(f"Kept {len(human_kept_markers)} markers with clear neighbors")
# Cluster genes (it generates a background expression pattern so incorporates the gene expression level)
human_Qavg = (human_analyzer.similarity_matrix * human_analyzer.gene_info['tot'].values.reshape(-1, 1)).sum(axis=0)
human_Qavg = human_Qavg / human_Qavg.sum()
human_cluster_idx, human_clusters, human_d2, human_anchor_mapping = cluster_by_anchor(
    human_selector,
    fixed_anchors=human_kept_markers,
    null_q=human_Qavg,
    gene_to_label=human_gene_to_label,
    verbose=1
)
# Add summary table
human_summary = create_cluster_summary(human_anchor_mapping, human_clusters, human_markers)
print("\n=== HUMAN CLUSTER SUMMARY ===")
print(human_summary)

# %% [markdown]
# ## SF19B: Human 20-anchor spatial patterns

# %% [markdown]
# ### Figure 5B

# %% tags=["cell-104"]
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
# The 20 anchors in order
key_genes = human_markers[:20]
# Get the anchors that actually formed clusters from the clustering results
cluster_forming_anchors = set(human_kept_markers)
# Create a dictionary mapping anchor names to cluster sizes from chick_summary
anchor_to_cluster_size = {}
if 'human_summary' in globals():
  for _, row in human_summary.iterrows():
      if row['anchor_name'] != 'NULL':
          anchor_to_cluster_size[row['anchor_name']] = row['cluster_size']
print(f"Cluster-forming anchors from data: {len(cluster_forming_anchors)} out of {len(key_genes)}")
# Create figure with subplots
fig, axes = plt.subplots(4, 5, figsize=(12.5, 10))
axes = axes.flatten()
# Gene to index mapping
gene2idx = {gene: i for i, gene in enumerate(human_analyzer.gene_names)}
# Calculate spatial max values for each gene
spatial_max_values = {}
for gene in key_genes:
  if gene in human_analyzer.gene_names:
      idx = gene2idx[gene]
      spatial_max_values[gene] = float(np.max(human_analyzer.images[idx]))
# Plot each gene
for i, gene in enumerate(key_genes):
  if gene in human_analyzer.gene_names:
      ax = axes[i]
      idx = gene2idx[gene]
      img = human_analyzer.images[idx].copy()
      # Apply mask if counts are available
      if human_analyzer.counts is not None:
          mask = (human_analyzer.counts < human_analyzer.params.mask_count_threshold)
          img = np.ma.masked_where(mask, img)
      # Plot the image
      im = ax.imshow(
          img,
          origin="lower",
          cmap="viridis",
          aspect="equal"
      )
      # Add reference lines if available
      if human_analyzer.dv_mid is not None and human_analyzer.nt_mid is not None:
          ax.axhline(y=human_analyzer.dv_mid, color="black", lw=0.5)
          ax.axvline(x=human_analyzer.nt_mid, color="black", lw=0.5)
      # Get spatial max value for this gene
      spatial_max = spatial_max_values[gene]
      # Check if this anchor formed a cluster and get cluster size
      forms_cluster = gene in cluster_forming_anchors
      # Create title with cluster size information
      if forms_cluster and gene in anchor_to_cluster_size:
          cluster_size = anchor_to_cluster_size[gene]
          # Format with comma for numbers >= 1000
          if cluster_size >= 1000:
              size_str = f"{cluster_size:,}"
          else:
              size_str = str(cluster_size)
          title_text = f"{gene} ({size_str})\nmax: {spatial_max:.2f}"
      else:
          # No cluster formed - just show gene name
          title_text = f"{gene}\nmax: {spatial_max:.2f}"
      ax.set_title(title_text, fontsize=11, fontweight='bold')
      # Remove axis ticks for cleaner look
      ax.set_xticks([])
      ax.set_yticks([])
      # Add red border for cluster-forming anchors
      if forms_cluster:
          # Create a rectangle patch for the border
          rect = patches.Rectangle((0, 0), 1, 1, linewidth=3,
                                  edgecolor='red', facecolor='none',
                                  transform=ax.transAxes)
          ax.add_patch(rect)
      # Add parameter info to bottom right of each subplot
      param_text = f"p={human_analyzer.params.percentile_clip:.2f}"
      ax.text(0.98, 0.02, param_text,
              transform=ax.transAxes,
              fontsize=8, ha='right', va='bottom',
              bbox=dict(boxstyle='round,pad=0.3',
                       facecolor='white', alpha=0.8))
      # Add spatial category label if available
      if gene in human_gene_to_label:
          category = human_gene_to_label[gene]
          ax.text(0.02, 0.98, category,
                 transform=ax.transAxes,
                 fontsize=8, ha='left', va='top',
                 bbox=dict(boxstyle='round,pad=0.3',
                          facecolor='yellow', alpha=0.5))
      # Optionally add axis labels to first subplot only
      if i == 0:
          ax.set_xlabel("Temporal ← NT → Nasal", fontsize=10)
          ax.set_ylabel("Ventral ← DV → Dorsal", fontsize=10)
# Add overall title with legend
# fig.suptitle("Spatial Expression Patterns of Top 20 Anchor Genes\n" +
#           "Red border = forms cluster, (n) = number of genes in cluster",
#           fontsize=14, fontweight='bold')
# Adjust layout and save
plt.tight_layout()
plt.subplots_adjust(top=0.90)  # Make room for suptitle with legend
os.makedirs(os.path.join(FIGURES_BASE, "Figure_SF19"), exist_ok=True)
plt.savefig(os.path.join(FIGURES_BASE, "Figure_SF19", "SF19_human_clusters.png"), dpi=300, bbox_inches='tight')
plt.show()
# Print detailed summary
print("\n" + "="*60)
print("CLUSTER FORMATION SUMMARY")
print("="*60)
print(f"Total anchors analyzed: {len(key_genes)}")
print(f"Anchors that formed clusters: {sum(1 for g in key_genes if g in cluster_forming_anchors)}")
print(f"Anchors that didn't form clusters: {sum(1 for g in key_genes if g not in cluster_forming_anchors)}")
print("\nCluster sizes for the 20 anchors:")
for gene in key_genes:
  if gene in anchor_to_cluster_size:
      size = anchor_to_cluster_size[gene]
      size_str = f"{size:,}" if size >= 1000 else str(size)
      category = human_gene_to_label.get(gene, "Unknown")
      print(f"  {gene:15} → {size_str:>6} genes ({category})")
  else:
      print(f"  {gene:15} → No cluster formed")
# Summary statistics
if anchor_to_cluster_size:
  sizes = [anchor_to_cluster_size[g] for g in key_genes if g in anchor_to_cluster_size]
  print(f"\nCluster size statistics for these 20 anchors:")
  print(f"  Mean: {np.mean(sizes):.1f} genes")
  print(f"  Median: {np.median(sizes):.1f} genes")
  print(f"  Range: {min(sizes)} - {max(sizes):,} genes")

# %% [markdown]
# ## Table S2: Human anchor gene pairs

# %% [markdown]
# ### Supplemental Table 2

# %% tags=["cell-106"]
import pandas as pd
# Create detailed list - ONLY genes actually in clusters
detailed_data = []
for _, row in human_summary.iterrows():
  if row['anchor_name'] != 'NULL':
      anchor = row['anchor_name']
      category = row['category']
      cluster_size = row['cluster_size']
      if anchor in human_analyzer.gene_names:
          # Get ALL cluster members
          cluster_genes = human_clusters.get(anchor, [])
          # Get enough correlations to cover all cluster members
          correlations = human_analyzer.get_gene_correlations(anchor,
                                                             top_n=len(human_analyzer.gene_names))
          corr_dict = dict(correlations)
          # Add each cluster member with its correlation
          for gene in cluster_genes:
              if gene in corr_dict:
                  # Find rank in correlation list
                  rank = next((i for i, (g, _) in enumerate(correlations, 1) if g == gene), None)
                  detailed_data.append({
                      'Anchor': anchor,
                      'Category': category,
                      'Cluster_Size': cluster_size,
                      'Correlation_Rank': rank,
                      'Gene': gene,
                      'Correlation': corr_dict[gene]
                  })
# Create DataFrame and sort
detailed_table = pd.DataFrame(detailed_data)
detailed_table = detailed_table.sort_values(['Anchor', 'Correlation_Rank'])
# Save to TSV
os.makedirs(os.path.join(FIGURES_BASE, "Tables"), exist_ok=True)
detailed_table.to_csv(os.path.join(FIGURES_BASE, "Tables", "TableS2_human_cluster_members.tsv"), sep='\t', index=False)
print(f"Saved {len(detailed_table)} cluster member entries")
print(f"This should match total of all cluster sizes: {human_summary[human_summary['anchor_name'] != 'NULL']['cluster_size'].sum()}")
# Show summary
print("\nAnchors included:")
print(detailed_table.groupby('Anchor')[['Category', 'Cluster_Size']].first())

# %% [markdown]
# ## Load mouse RPC data

# %% tags=["cell-120"]
import anndata as ad
# mouse_adata = ad.read_h5ad("../data/20240815_mouse_RPC.h5ad")
mouse_adata = ad.read_h5ad("../data/20250604_mouse_RPC.h5ad")
print(f"Loaded mouse data: {mouse_adata.n_obs:,} cells × {mouse_adata.n_vars:,} genes")
# Derive parameters transparently based on the data characteristics
mouse_params = SpatialAnalysisParams.derive_parameters_from_data(
    n_cells=mouse_adata.n_obs, 
    species='mouse'
)
mouse_params

# %% tags=["cell-121"]
mouse_params.explain_parameters()

# %% tags=["cell-122"]
# You could potentially change the parameters like this:
mouse_params.bin_size = 51
mouse_params.smooth_sigma=1.0
mouse_params.min_gene_count=30 # changing to match human/chick

# %% [markdown]
# ## Create mouse spatial analyzer

# %% tags=["cell-124"]
# Initialize analyzer
mouse_analyzer = SpatialExpressionAnalyzer(mouse_params)
# Run full analysis (this replaces all your manual preprocessing steps)
# Update path to point to data file location relative to notebooks/ directory
#mouse_results = mouse_analyzer.run_full_analysis("../data/20240815_mouse_RPC.h5ad")
mouse_results = mouse_analyzer.run_full_analysis("../data/20250604_mouse_RPC.h5ad")

# %% [markdown]
# ## Mouse anchor selection

# %% [markdown]
# ## Selection of anchors
# With some fixed spatial patterns that we know of, we can maximize by choosing anchors by greedy approach. The anchors are selected based on maximum expression.

# %% tags=["cell-131"]
# To customize anchors, edit `control_genes.yaml` and call reload_control_genes()
reload_control_genes()

# %% tags=["cell-132"]
mouse_fixed_anchors_dict = get_fixed_anchors(mouse_analyzer.gene_info, 'mouse')
print("Fixed anchors from control_genes.yaml:")
for category, gene in mouse_fixed_anchors_dict.items():
    expression = mouse_analyzer.gene_info[mouse_analyzer.gene_info['name'] == gene]['tot'].iloc[0]
    print(f"  {category}: {gene} (expression: {expression:.1f})")
# Extract gene names for marker selection
mouse_fixed_anchors = list(mouse_fixed_anchors_dict.values())
print(f"\nUsing {len(mouse_fixed_anchors)} anchor genes for marker selection")
# Set up marker selector (replaces your MarkerSelectorPy setup)
valid_mask = mouse_analyzer.gene_info['tot'] >= 1000
mouse_selector = MarkerSelectorPy(
    gene_names=mouse_analyzer.gene_names,
    Q=mouse_analyzer.similarity_matrix,
    valid_mask=valid_mask,
    verbose=1,
    rownorm=False
)
# Select markers (replaces your select_markers call)
mouse_markers = mouse_selector.select_markers(K=20, fixed_anchors=mouse_fixed_anchors)
print("Selected markers:")
for i, marker in enumerate(mouse_markers, 1):
    expr = mouse_analyzer.gene_info[mouse_analyzer.gene_info['name'] == marker]['tot'].iloc[0]
    print(f"{i:2d}. {marker} ({expr:.1f})")

# %% tags=["cell-135"]
# First, create the gene-to-label mapping from control genes
mouse_gene_to_label = {}
mouse_control_genes = CONTROL_GENES.get('mouse', {})
for category, genes in mouse_control_genes.items():
    for gene in genes:
        if gene in mouse_analyzer.gene_names:
            mouse_gene_to_label[gene] = category

# %% tags=["cell-137"]
mouse_neighbors = mouse_selector.find_neighbors_to_anchors(
    m=10, # maximum candidate neighbors ("shortlist" size)
    max_rank_fraction=50/len(mouse_analyzer.gene_names), # acceptable level of rank quality. Only top 50 as rank
    gene_to_label=mouse_gene_to_label
)
# Keep only markers with neighbors
mouse_kept_markers = [nbr.name for nbr in mouse_neighbors if len(nbr.neighbors) > 0]
print(f"Kept {len(mouse_kept_markers)} markers with clear neighbors")
# Cluster genes (it generates a background expression pattern so incorporates the gene expression level)
mouse_Qavg = (mouse_analyzer.similarity_matrix * mouse_analyzer.gene_info['tot'].values.reshape(-1, 1)).sum(axis=0)
mouse_Qavg = mouse_Qavg / mouse_Qavg.sum()
mouse_cluster_idx, mouse_clusters, mouse_d2, mouse_anchor_mapping = cluster_by_anchor(
    mouse_selector, 
    fixed_anchors=mouse_kept_markers, 
    null_q=mouse_Qavg, 
    gene_to_label=mouse_gene_to_label,
    verbose=1
)
# Add summary table
mouse_summary = create_cluster_summary(mouse_anchor_mapping, mouse_clusters, mouse_markers)
print("\n=== MOUSE CLUSTER SUMMARY ===")
print(mouse_summary)

# %% [markdown]
# ## SF19A: Mouse 20-anchor spatial patterns

# %% [markdown]
# ### Figure 5D

# %% tags=["cell-141"]
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
# The 20 anchors in order
key_genes = mouse_markers[:20]
# Get the anchors that actually formed clusters from the clustering results
cluster_forming_anchors = set(mouse_kept_markers)
# Create a dictionary mapping anchor names to cluster sizes from chick_summary
anchor_to_cluster_size = {}
if 'mouse_summary' in globals():
  for _, row in mouse_summary.iterrows():
      if row['anchor_name'] != 'NULL':
          anchor_to_cluster_size[row['anchor_name']] = row['cluster_size']
print(f"Cluster-forming anchors from data: {len(cluster_forming_anchors)} out of {len(key_genes)}")
# Create figure with subplots
fig, axes = plt.subplots(4, 5, figsize=(12.5, 10))
axes = axes.flatten()
# Gene to index mapping
gene2idx = {gene: i for i, gene in enumerate(mouse_analyzer.gene_names)}
# Calculate spatial max values for each gene
spatial_max_values = {}
for gene in key_genes:
  if gene in mouse_analyzer.gene_names:
      idx = gene2idx[gene]
      spatial_max_values[gene] = float(np.max(mouse_analyzer.images[idx]))
# Plot each gene
for i, gene in enumerate(key_genes):
  if gene in mouse_analyzer.gene_names:
      ax = axes[i]
      idx = gene2idx[gene]
      img = mouse_analyzer.images[idx].copy()
      # Apply mask if counts are available
      if mouse_analyzer.counts is not None:
          mask = (mouse_analyzer.counts < mouse_analyzer.params.mask_count_threshold)
          img = np.ma.masked_where(mask, img)
      # Plot the image
      im = ax.imshow(
          img,
          origin="lower",
          cmap="viridis",
          aspect="equal"
      )
      # Add reference lines if available
      if mouse_analyzer.dv_mid is not None and mouse_analyzer.nt_mid is not None:
          ax.axhline(y=mouse_analyzer.dv_mid, color="black", lw=0.5)
          ax.axvline(x=mouse_analyzer.nt_mid, color="black", lw=0.5)
      # Get spatial max value for this gene
      spatial_max = spatial_max_values[gene]
      # Check if this anchor formed a cluster and get cluster size
      forms_cluster = gene in cluster_forming_anchors
      # Create title with cluster size information
      if forms_cluster and gene in anchor_to_cluster_size:
          cluster_size = anchor_to_cluster_size[gene]
          # Format with comma for numbers >= 1000
          if cluster_size >= 1000:
              size_str = f"{cluster_size:,}"
          else:
              size_str = str(cluster_size)
          title_text = f"{gene} ({size_str})\nmax: {spatial_max:.2f}"
      else:
          # No cluster formed - just show gene name
          title_text = f"{gene}\nmax: {spatial_max:.2f}"
      ax.set_title(title_text, fontsize=11, fontweight='bold')
      # Remove axis ticks for cleaner look
      ax.set_xticks([])
      ax.set_yticks([])
      # Add red border for cluster-forming anchors
      if forms_cluster:
          # Create a rectangle patch for the border
          rect = patches.Rectangle((0, 0), 1, 1, linewidth=3,
                                  edgecolor='red', facecolor='none',
                                  transform=ax.transAxes)
          ax.add_patch(rect)
      # Add parameter info to bottom right of each subplot
      param_text = f"p={mouse_analyzer.params.percentile_clip:.2f}"
      ax.text(0.98, 0.02, param_text,
              transform=ax.transAxes,
              fontsize=8, ha='right', va='bottom',
              bbox=dict(boxstyle='round,pad=0.3',
                       facecolor='white', alpha=0.8))
      # Add spatial category label if available
      if gene in mouse_gene_to_label:
          category = mouse_gene_to_label[gene]
          ax.text(0.02, 0.98, category,
                 transform=ax.transAxes,
                 fontsize=8, ha='left', va='top',
                 bbox=dict(boxstyle='round,pad=0.3',
                          facecolor='yellow', alpha=0.5))
      # Optionally add axis labels to first subplot only
      if i == 0:
          ax.set_xlabel("Temporal ← NT → Nasal", fontsize=10)
          ax.set_ylabel("Ventral ← DV → Dorsal", fontsize=10)
# Add overall title with legend
# fig.suptitle("Spatial Expression Patterns of Top 20 Anchor Genes\n" +
#           "Red border = forms cluster, (n) = number of genes in cluster",
#           fontsize=14, fontweight='bold')
# Adjust layout and save
plt.tight_layout()
plt.subplots_adjust(top=0.90)  # Make room for suptitle with legend
plt.savefig(os.path.join(FIGURES_BASE, "Figure_SF19", "SF19_mouse_clusters.png"), dpi=300, bbox_inches='tight')
plt.show()
# Print detailed summary
print("\n" + "="*60)
print("CLUSTER FORMATION SUMMARY")
print("="*60)
print(f"Total anchors analyzed: {len(key_genes)}")
print(f"Anchors that formed clusters: {sum(1 for g in key_genes if g in cluster_forming_anchors)}")
print(f"Anchors that didn't form clusters: {sum(1 for g in key_genes if g not in cluster_forming_anchors)}")
print("\nCluster sizes for the 20 anchors:")
for gene in key_genes:
  if gene in anchor_to_cluster_size:
      size = anchor_to_cluster_size[gene]
      size_str = f"{size:,}" if size >= 1000 else str(size)
      category = mouse_gene_to_label.get(gene, "Unknown")
      print(f"  {gene:15} → {size_str:>6} genes ({category})")
  else:
      print(f"  {gene:15} → No cluster formed")
# Summary statistics
if anchor_to_cluster_size:
  sizes = [anchor_to_cluster_size[g] for g in key_genes if g in anchor_to_cluster_size]
  print(f"\nCluster size statistics for these 20 anchors:")
  print(f"  Mean: {np.mean(sizes):.1f} genes")
  print(f"  Median: {np.median(sizes):.1f} genes")
  print(f"  Range: {min(sizes)} - {max(sizes):,} genes")

# %% [markdown]
# ## Table S3: Mouse anchor gene pairs

# %% [markdown]
# ### Supplemental Table 3

# %% tags=["cell-144"]
import pandas as pd
# Create detailed list - ONLY genes actually in clusters
detailed_data = []
for _, row in mouse_summary.iterrows():
  if row['anchor_name'] != 'NULL':
      anchor = row['anchor_name']
      category = row['category']
      cluster_size = row['cluster_size']
      if anchor in mouse_analyzer.gene_names:
          # Get ALL cluster members
          cluster_genes = mouse_clusters.get(anchor, [])
          # Get enough correlations to cover all cluster members
          correlations = mouse_analyzer.get_gene_correlations(anchor,
                                                             top_n=len(mouse_analyzer.gene_names))
          corr_dict = dict(correlations)
          # Add each cluster member with its correlation
          for gene in cluster_genes:
              if gene in corr_dict:
                  # Find rank in correlation list
                  rank = next((i for i, (g, _) in enumerate(correlations, 1) if g == gene), None)
                  detailed_data.append({
                      'Anchor': anchor,
                      'Category': category,
                      'Cluster_Size': cluster_size,
                      'Correlation_Rank': rank,
                      'Gene': gene,
                      'Correlation': corr_dict[gene]
                  })
# Create DataFrame and sort
detailed_table = pd.DataFrame(detailed_data)
detailed_table = detailed_table.sort_values(['Anchor', 'Correlation_Rank'])
# Save to TSV
detailed_table.to_csv(os.path.join(FIGURES_BASE, "Tables", "TableS3_mouse_cluster_members.tsv"), sep='\t', index=False)
print(f"Saved {len(detailed_table)} cluster member entries")
print(f"This should match total of all cluster sizes: {mouse_summary[mouse_summary['anchor_name'] != 'NULL']['cluster_size'].sum()}")
# Show summary
print("\nAnchors included:")
print(detailed_table.groupby('Anchor')[['Category', 'Cluster_Size']].first())
