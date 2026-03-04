#!/usr/bin/env python3
"""
Basic tests for the spatial expression analysis modules
"""

import unittest
import sys
import os
import numpy as np
import pandas as pd

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

class TestModuleImports(unittest.TestCase):
    """Test that all modules can be imported correctly"""
    
    def test_import_spatial_analysis(self):
        """Test importing spatial expression analysis module"""
        try:
            from spatial_expression_analysis import (
                SpatialExpressionAnalyzer,
                SpatialAnalysisParams,
                MarkerSelectorPy
            )
            self.assertTrue(True)
        except ImportError as e:
            self.fail(f"Failed to import spatial_expression_analysis: {e}")
    

class TestSpatialAnalysisParams(unittest.TestCase):
    """Test SpatialAnalysisParams functionality"""
    
    def setUp(self):
        from spatial_expression_analysis import SpatialAnalysisParams
        self.SpatialAnalysisParams = SpatialAnalysisParams
    
    def test_default_params(self):
        """Test default parameter creation"""
        params = self.SpatialAnalysisParams()
        self.assertIsInstance(params.bin_size, int)
        self.assertIsInstance(params.smooth_sigma, float)
        self.assertGreater(params.bin_size, 0)
        self.assertGreaterEqual(params.smooth_sigma, 0)
    
    def test_species_defaults(self):
        """Test species-specific parameter defaults"""
        species_list = ['chick', 'human', 'mouse']
        
        for species in species_list:
            params = self.SpatialAnalysisParams.get_species_defaults(species)
            self.assertIsInstance(params, self.SpatialAnalysisParams)
            self.assertGreater(params.bin_size, 0)
            self.assertGreaterEqual(params.smooth_sigma, 0)
    
    def test_invalid_species(self):
        """Test that invalid species raises error"""
        with self.assertRaises(ValueError):
            self.SpatialAnalysisParams.get_species_defaults('invalid_species')

class TestMarkerSelector(unittest.TestCase):
    """Test MarkerSelectorPy functionality with mock data"""
    
    def setUp(self):
        from spatial_expression_analysis import MarkerSelectorPy
        self.MarkerSelectorPy = MarkerSelectorPy
        
        # Create mock data
        np.random.seed(42)  # For reproducible tests
        self.n_genes = 100
        self.n_pixels = 50
        
        # Mock gene names
        self.gene_names = [f"GENE_{i:03d}" for i in range(self.n_genes)]
        
        # Mock similarity matrix (positive definite)
        random_matrix = np.random.randn(self.n_genes, self.n_genes)
        self.similarity_matrix = np.dot(random_matrix, random_matrix.T)
        
        # Normalize to make it look like cosine similarity
        norms = np.linalg.norm(self.similarity_matrix, axis=1, keepdims=True)
        self.similarity_matrix = self.similarity_matrix / (norms * norms.T)
        np.fill_diagonal(self.similarity_matrix, 1.0)  # Perfect self-similarity
    
    def test_marker_selector_creation(self):
        """Test creating MarkerSelectorPy instance"""
        selector = self.MarkerSelectorPy(
            gene_names=self.gene_names,
            Q=self.similarity_matrix,
            verbose=0
        )
        
        self.assertEqual(len(selector.gene_names), self.n_genes)
        self.assertEqual(selector.M, self.n_genes)
        self.assertIsNotNone(selector.name2idx)
    
    def test_marker_selection(self):
        """Test marker selection functionality"""
        selector = self.MarkerSelectorPy(
            gene_names=self.gene_names,
            Q=self.similarity_matrix,
            verbose=0
        )
        
        K = 5
        selected_markers = selector.select_markers(K)
        
        self.assertEqual(len(selected_markers), K)
        self.assertTrue(all(marker in self.gene_names for marker in selected_markers))
        self.assertEqual(len(set(selected_markers)), K)  # All unique

if __name__ == '__main__':
    # Create a test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    # Add test classes
    suite.addTests(loader.loadTestsFromTestCase(TestModuleImports))
    suite.addTests(loader.loadTestsFromTestCase(TestSpatialAnalysisParams))
    suite.addTests(loader.loadTestsFromTestCase(TestMarkerSelector))

    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)

    # Exit with non-zero code if tests failed
    if not result.wasSuccessful():
        sys.exit(1)