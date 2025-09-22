"""
Lung Cancer DGE Analysis Pipeline
==================================
A comprehensive pipeline for differential gene expression analysis
"""

from .data_loader import GEODataLoader
from .preprocessing import DataPreprocessor  
from .deg_analysis import DifferentialExpression
from .visualization import ExpressionVisualizer

__version__ = "1.0.0"
__author__ = "Ananaya"
__all__ = [
    "GEODataLoader",
    "DataPreprocessor",
    "DifferentialExpression", 
    "ExpressionVisualizer"
]

# Main pipeline class for easy access
class LungCancerPipeline:
    """Main pipeline wrapper for complete analysis"""
    
    def __init__(self):
        self.loader = GEODataLoader()
        self.preprocessor = DataPreprocessor()
        self.de_analyzer = DifferentialExpression()
        self.visualizer = ExpressionVisualizer()
    
    def analyze(self, dataset_id, output_dir='results'):
        """
        Run complete analysis on a dataset
        
        Parameters:
            dataset_id: GEO dataset ID
            output_dir: Output directory for results
            
        Returns:
            Dictionary with results
        """
        
        # Load data
        expr_df, metadata = self.loader.load_dataset(dataset_id)
        
        # Preprocess
        expr_filtered = self.preprocessor.filter_low_expression(expr_df)
        expr_normalized = self.preprocessor.normalize(expr_filtered)
        
        # Differential expression
        de_results = self.de_analyzer.analyze(expr_normalized, metadata)
        
        # Visualize
        self.visualizer.create_all_plots(expr_normalized, de_results, metadata)
        
        # Save results
        de_results.to_csv(f'{output_dir}/deg_results.csv', index=False)
        
        return {
            'expression': expr_normalized,
            'de_results': de_results,
            'metadata': metadata
        }
    
    def multi_dataset_analysis(self, dataset_ids):
        """Analyze multiple datasets"""
        
        results = {}
        for dataset_id in dataset_ids:
            results[dataset_id] = self.analyze(dataset_id)
        
        return results
