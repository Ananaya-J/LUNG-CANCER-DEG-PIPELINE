"""
Visualization module for gene expression analysis
"""

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

class ExpressionVisualizer:
    """Create publication-quality visualizations"""
    
    def __init__(self, output_dir='results'):
        self.output_dir = output_dir
        
    def create_all_plots(self, expr_df, de_results, metadata):
        """
        Generate all standard plots
        
        Parameters:
            expr_df: Expression dataframe
            de_results: DE analysis results
            metadata: Sample metadata
        """
        
        # Create figure with subplots
        fig = plt.figure(figsize=(20, 15))
        
        # 1. Volcano plot
        ax1 = plt.subplot(2, 3, 1)
        self._volcano_plot(ax1, de_results)
        
        # 2. MA plot
        ax2 = plt.subplot(2, 3, 2)
        self._ma_plot(ax2, de_results)
        
        # 3. P-value histogram
        ax3 = plt.subplot(2, 3, 3)
        self._pvalue_histogram(ax3, de_results)
        
        # 4. PCA plot
        ax4 = plt.subplot(2, 3, 4)
        self._pca_plot(ax4, expr_df, metadata)
        
        # 5. Top genes heatmap
        ax5 = plt.subplot(2, 3, 5)
        self._gene_heatmap(ax5, expr_df, de_results)
        
        # 6. DEG summary
        ax6 = plt.subplot(2, 3, 6)
        self._deg_summary(ax6, de_results)
        
        plt.suptitle('Lung Cancer Differential Expression Analysis', fontsize=16, fontweight='bold')
        plt.tight_layout()
        
        # Save figure
        plt.savefig(f'{self.output_dir}/all_plots.png', dpi=300, bbox_inches='tight')
        plt.show()
        
    def _volcano_plot(self, ax, de_results):
        """Create volcano plot"""
        
        x = de_results['log2FoldChange']
        y = -np.log10(de_results['padj'] + 1e-300)
        
        # Plot by significance
        colors = {'UP': 'red', 'DOWN': 'blue', 'NS': 'gray'}
        
        for direction in ['NS', 'DOWN', 'UP']:
            mask = de_results['direction'] == direction
            ax.scatter(x[mask], y[mask], 
                      c=colors[direction], 
                      alpha=0.6 if direction == 'NS' else 0.8,
                      s=5, label=direction)
        
        # Add threshold lines
        ax.axhline(-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
        ax.axvline(-1, color='black', linestyle='--', alpha=0.5)
        ax.axvline(1, color='black', linestyle='--', alpha=0.5)
        
        ax.set_xlabel('Log2 Fold Change')
        ax.set_ylabel('-Log10(Adjusted P-value)')
        ax.set_title('Volcano Plot')
        ax.legend()
        
    def _ma_plot(self, ax, de_results):
        """Create MA plot"""
        
        A = de_results['baseMean']
        M = de_results['log2FoldChange']
        
        colors = {'UP': 'red', 'DOWN': 'blue', 'NS': 'gray'}
        
        for direction in ['NS', 'DOWN', 'UP']:
            mask = de_results['direction'] == direction
            ax.scatter(A[mask], M[mask], 
                      c=colors[direction], 
                      alpha=0.6 if direction == 'NS' else 0.8,
                      s=5)
        
        ax.axhline(0, color='black', linestyle='-', alpha=0.5)
        ax.axhline(-1, color='black', linestyle='--', alpha=0.5)
        ax.axhline(1, color='black', linestyle='--', alpha=0.5)
        
        ax.set_xlabel('Average Expression')
        ax.set_ylabel('Log2 Fold Change')
        ax.set_title('MA Plot')
        ax.set_xscale('log')
        
    def _pvalue_histogram(self, ax, de_results):
        """Create p-value distribution histogram"""
        
        ax.hist(de_results['pvalue'], bins=50, alpha=0.7, color='skyblue', edgecolor='black')
        ax.set_xlabel('P-value')
        ax.set_ylabel('Frequency')
        ax.set_title('P-value Distribution')
        ax.axvline(0.05, color='red', linestyle='--', alpha=0.5, label='p=0.05')
        ax.legend()
        
    def _pca_plot(self, ax, expr_df, metadata):
        """Create PCA plot"""
        
        # Prepare data
        samples = metadata['sample_id'].tolist()
        samples_in_expr = [s for s in samples if s in expr_df.columns]
        X = expr_df[samples_in_expr].T.values
        
        # Standardize
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        
        # PCA
        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(X_scaled)
        
        # Plot
        conditions = []
        for s in samples_in_expr:
            cond = metadata[metadata['sample_id'] == s]['condition'].values[0]
            conditions.append(cond)
        
        for condition in set(conditions):
            indices = [i for i, c in enumerate(conditions) if c == condition]
            color = 'blue' if condition == 'Normal' else 'red'
            ax.scatter(pca_result[indices, 0], 
                      pca_result[indices, 1],
                      label=condition, alpha=0.7, s=50, color=color)
        
        ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%})')
        ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%})')
        ax.set_title('PCA Plot')
        ax.legend()
        
    def _gene_heatmap(self, ax, expr_df, de_results):
        """Create heatmap of top DEGs"""
        
        # Get top genes
        sig_genes = de_results[de_results['significant']]
        
        if len(sig_genes) > 0:
            # Select top genes
            top_up = sig_genes[sig_genes['direction'] == 'UP'].head(10)
            top_down = sig_genes[sig_genes['direction'] == 'DOWN'].head(10)
            top_genes = pd.concat([top_up, top_down])['gene'].tolist()
            
            # Get expression data
            heatmap_data = expr_df.loc[top_genes]
            
            # Z-score normalize
            heatmap_zscore = (heatmap_data.T - heatmap_data.T.mean()) / heatmap_data.T.std()
            heatmap_zscore = heatmap_zscore.T
            
            # Plot
            im = ax.imshow(heatmap_zscore, cmap='RdBu_r', aspect='auto', vmin=-2, vmax=2)
            ax.set_xlabel('Samples')
            ax.set_ylabel('Genes')
            ax.set_title(f'Top {len(top_genes)} DEGs Heatmap')
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        else:
            ax.text(0.5, 0.5, 'No significant DEGs', 
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_title('Heatmap')
            
    def _deg_summary(self, ax, de_results):
        """Create DEG summary bar plot"""
        
        up_count = (de_results['direction'] == 'UP').sum()
        down_count = (de_results['direction'] == 'DOWN').sum()
        ns_count = (de_results['direction'] == 'NS').sum()
        
        categories = ['Upregulated', 'Downregulated', 'Not Significant']
        counts = [up_count, down_count, ns_count]
        colors = ['red', 'blue', 'gray']
        
        bars = ax.bar(categories, counts, color=colors, alpha=0.7)
        
        # Add value labels on bars
        for bar, count in zip(bars, counts):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{count:,}', ha='center', va='bottom')
        
        ax.set_ylabel('Number of Genes')
        ax.set_title('DEG Summary')
        ax.set_ylim(0, max(counts) * 1.1)
