"""
retina_spatial_scrna: Spatial gene expression analysis for retinal scRNA-seq data.

Provides tools for spatial binning, DV/NT scoring, cross-species comparison,
and NMF-based pattern analysis across chick, human, and mouse retinal datasets.
"""

from retina_spatial_scrna.spatial_expression_analysis import (
    SpatialAnalysisParams,
    SpatialExpressionAnalyzer,
    MarkerSelectorPy,
)

__version__ = "0.1.0"
