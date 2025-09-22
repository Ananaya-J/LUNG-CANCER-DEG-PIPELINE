"""
Differential expression analysis module
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import false_discovery_control

class DifferentialExpression:
    """Perform differential expression analysis"""
    
    def __init__(self, fdr_threshold=0.05, log2fc_threshold=1.0):
        self.fdr_threshold = fdr_threshold
        self.log2fc_threshold = log2fc_threshold
    
    def analyze(self, expr_df, metadata):
        """
        Perform differential expression analysis
        
        Parameters:
            expr_df: Normalized expression dataframe
            metadata: Sample metadata with 'condition' column
            
        Returns:
            DataFrame with DE results
        """
        
        # Separate samples by condition
        normal_samples = metadata[metadata['condition'] == 'Normal']['sample_id'].tolist()
        tumor_samples = metadata[metadata['condition'] == 'Tumor']['sample_id'].tolist()
        
        # Filter to existing samples
        normal_samples = [s for s in normal_samples if s in expr_df.columns]
        tumor_samples = [s for s in tumor_samples if s in expr_df.columns]
        
        print(f"  Comparing {len(tumor_samples)} tumor vs {len(normal_samples)} normal samples")
        
        results = []
        
        # Calculate statistics for each gene
        for gene in expr_df.index:
            normal_vals = expr_df.loc[gene, normal_samples].values
            tumor_vals = expr_df.loc[gene, tumor_samples].values
            
            # Remove NaN values
            normal_vals = normal_vals[~np.isnan(normal_vals)]
            tumor_vals = tumor_vals[~np.isnan(tumor_vals)]
            
            if len(normal_vals) > 1 and len(tumor_vals) > 1:
                # Calculate means
                normal_mean = np.mean(normal_vals)
                tumor_mean = np.mean(tumor_vals)
                
                # Log2 fold change (data already log2 transformed)
                log2fc = tumor_mean - normal_mean
                
                # T-test (Welch's t-test for unequal variance)
                t_stat, p_value = stats.ttest_ind(
                    tumor_vals, 
                    normal_vals,
                    equal_var=False
                )
                
                # Effect size (Cohen's d)
                pooled_std = np.sqrt(
                    (np.var(normal_vals) + np.var(tumor_vals)) / 2
                )
                cohen_d = (tumor_mean - normal_mean) / pooled_std if pooled_std > 0 else 0
                
                results.append({
                    'gene': gene,
                    'baseMean': (normal_mean + tumor_mean) / 2,
                    'log2FoldChange': log2fc,
                    'pvalue': p_value,
                    'stat': t_stat,
                    'cohen_d': cohen_d,
                    'normal_mean': normal_mean,
                    'tumor_mean': tumor_mean
                })
        
        # Create results dataframe
        de_results = pd.DataFrame(results)
        
        # Multiple testing correction
        de_results['padj'] = false_discovery_control(de_results['pvalue'].values)
        
        # Determine significance
        de_results['significant'] = (
            (de_results['padj'] < self.fdr_threshold) & 
            (np.abs(de_results['log2FoldChange']) > self.log2fc_threshold)
        )
        
        # Add direction
        de_results['direction'] = 'NS'  # Not significant
        de_results.loc[
            (de_results['significant']) & (de_results['log2FoldChange'] > 0),
            'direction'
        ] = 'UP'
        de_results.loc[
            (de_results['significant']) & (de_results['log2FoldChange'] < 0),
            'direction'
        ] = 'DOWN'
        
        # Sort by adjusted p-value
        de_results = de_results.sort_values('padj')
        
        # Print summary
        n_sig = de_results['significant'].sum()
        n_up = (de_results['direction'] == 'UP').sum()
        n_down = (de_results['direction'] == 'DOWN').sum()
        
        print(f"  Significant genes: {n_sig} (UP: {n_up}, DOWN: {n_down})")
        
        return de_results
    
    def get_top_genes(self, de_results, n=20):
        """
        Get top differentially expressed genes
        
        Parameters:
            de_results: DE analysis results
            n: Number of top genes to return
            
        Returns:
            DataFrame with top genes
        """
        
        significant = de_results[de_results['significant']]
        
        # Get top upregulated and downregulated
        top_up = significant[significant['direction'] == 'UP'].head(n // 2)
        top_down = significant[significant['direction'] == 'DOWN'].head(n // 2)
        
        return pd.concat([top_up, top_down])
