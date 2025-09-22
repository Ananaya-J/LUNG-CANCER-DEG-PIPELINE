"""
Data loading module for GEO datasets
"""

import pandas as pd
import numpy as np
import urllib.request
import gzip
import io

class GEODataLoader:
    """Load and parse GEO datasets"""
    
    def __init__(self):
        self.supported_datasets = {
            'GSE31210': 'Lung adenocarcinoma',
            'GSE19188': 'Non-small cell lung cancer',
            'GSE32863': 'Lung adenocarcinoma',
            'GSE43458': 'Lung adenocarcinoma'
        }
    
    def load_dataset(self, dataset_id):
        """
        Load a GEO dataset
        
        Parameters:
            dataset_id: GEO accession (e.g., 'GSE31210')
            
        Returns:
            expr_df: Expression dataframe (genes x samples)
            metadata: Sample metadata
        """
        
        # For demonstration, create example data
        # In production, this would download real GEO data
        
        np.random.seed(42)
        
        # Generate example data based on dataset
        if dataset_id == 'GSE31210':
            n_genes, n_tumor, n_normal = 20000, 226, 20
        elif dataset_id == 'GSE19188':
            n_genes, n_tumor, n_normal = 20000, 91, 65
        else:
            n_genes, n_tumor, n_normal = 15000, 60, 30
        
        # Create gene IDs
        gene_ids = [f"GENE_{i:05d}" for i in range(n_genes)]
        
        # Known cancer genes (will be made differential)
        cancer_genes_idx = {
            0: 'EGFR', 1: 'KRAS', 2: 'ALK', 3: 'MET',
            4: 'ERBB2', 5: 'PIK3CA', 6: 'BRAF',
            7: 'TP53', 8: 'STK11', 9: 'KEAP1',
            10: 'CDKN2A', 11: 'RB1', 12: 'PTEN'
        }
        
        # Generate expression matrix
        n_samples = n_tumor + n_normal
        expression_matrix = np.random.normal(8, 2, (n_genes, n_samples))
        
        # Make cancer genes differential
        for idx, gene in cancer_genes_idx.items():
            if idx < 7:  # Oncogenes - upregulated
                expression_matrix[idx, :n_normal] -= 2
                expression_matrix[idx, n_normal:] += 2
            else:  # Tumor suppressors - downregulated
                expression_matrix[idx, :n_normal] += 2
                expression_matrix[idx, n_normal:] -= 2
        
        # Add some random DEGs
        random_degs = np.random.choice(range(13, n_genes), 500, replace=False)
        for idx in random_degs:
            if np.random.random() > 0.5:
                expression_matrix[idx, :n_normal] -= 1.5
                expression_matrix[idx, n_normal:] += 1.5
            else:
                expression_matrix[idx, :n_normal] += 1.5
                expression_matrix[idx, n_normal:] -= 1.5
        
        # Convert to actual expression values
        expression_matrix = 2 ** expression_matrix
        
        # Create dataframe
        sample_ids = [f"Sample_{i:03d}" for i in range(n_samples)]
        expr_df = pd.DataFrame(expression_matrix, index=gene_ids, columns=sample_ids)
        
        # Create metadata
        metadata = pd.DataFrame({
            'sample_id': sample_ids,
            'condition': ['Normal'] * n_normal + ['Tumor'] * n_tumor
        })
        
        return expr_df, metadata
