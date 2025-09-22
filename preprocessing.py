"""
Data preprocessing and normalization module
"""

import pandas as pd
import numpy as np
from scipy import stats

class DataPreprocessor:
    """Preprocess and normalize gene expression data"""
    
    def __init__(self):
        self.min_expression = 1.0
        self.min_sample_percent = 0.2
    
    def filter_low_expression(self, expr_df, min_expr=None, min_samples=None):
        """
        Filter genes with low expression
        
        Parameters:
            expr_df: Expression dataframe
            min_expr: Minimum expression threshold
            min_samples: Minimum fraction of samples
            
        Returns:
            Filtered expression dataframe
        """
        
        if min_expr is None:
            min_expr = self.min_expression
        if min_samples is None:
            min_samples = self.min_sample_percent
        
        # Calculate number of samples threshold
        n_samples_threshold = int(min_samples * expr_df.shape[1])
        
        # Filter genes
        genes_to_keep = (expr_df > min_expr).sum(axis=1) >= n_samples_threshold
        
        print(f"  Filtering: {expr_df.shape[0]} -> {genes_to_keep.sum()} genes")
        
        return expr_df.loc[genes_to_keep]
    
    def normalize(self, expr_df, method='log2'):
        """
        Normalize expression data
        
        Parameters:
            expr_df: Expression dataframe
            method: Normalization method ('log2', 'zscore', 'quantile')
            
        Returns:
            Normalized expression dataframe
        """
        
        if method == 'log2':
            # Log2 transformation with pseudocount
            normalized = np.log2(expr_df + 1)
            
        elif method == 'zscore':
            # Z-score normalization per gene
            normalized = expr_df.sub(expr_df.mean(axis=1), axis=0)
            normalized = normalized.div(expr_df.std(axis=1), axis=0)
            normalized = normalized.fillna(0)
            
        elif method == 'quantile':
            # Quantile normalization
            from scipy.stats import rankdata
            
            ranks = expr_df.apply(rankdata, axis=0)
            sorted_df = np.sort(expr_df.values, axis=0)
            mean_values = sorted_df.mean(axis=1)
            
            normalized_values = np.zeros_like(expr_df.values)
            for i in range(expr_df.shape[1]):
                normalized_values[ranks.iloc[:, i].astype(int) - 1, i] = mean_values
            
            normalized = pd.DataFrame(
                normalized_values,
                index=expr_df.index,
                columns=expr_df.columns
            )
        else:
            normalized = expr_df
        
        print(f"  Normalization method: {method}")
        
        return normalized
    
    def remove_batch_effects(self, expr_df, batch_info):
        """
        Remove batch effects using simple linear model
        
        Parameters:
            expr_df: Expression dataframe
            batch_info: Series with batch information
            
        Returns:
            Batch-corrected expression dataframe
        """
        
        corrected_df = expr_df.copy()
        
        for batch in batch_info.unique():
            batch_samples = batch_info[batch_info == batch].index
            batch_mean = expr_df[batch_samples].mean(axis=1)
            overall_mean = expr_df.mean(axis=1)
            
            for sample in batch_samples:
                if sample in corrected_df.columns:
                    corrected_df[sample] = corrected_df[sample] - batch_mean + overall_mean
        
        print("  Batch correction applied")
        
        return corrected_df
