#!/usr/bin/env python3
"""
Lung Cancer Differential Gene Expression Analysis Pipeline
Main execution script
"""

import argparse
import sys
import os
from datetime import datetime
import yaml

# Add src to path
sys.path.append('src')

from src.data_loader import GEODataLoader
from src.preprocessing import DataPreprocessor
from src.deg_analysis import DifferentialExpression
from src.visualization import ExpressionVisualizer

def main():
    """Main pipeline execution"""
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Lung Cancer DGE Analysis Pipeline'
    )
    parser.add_argument(
        '--dataset',
        type=str,
        default='GSE31210',
        help='GEO dataset ID (e.g., GSE31210)'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='results',
        help='Output directory'
    )
    parser.add_argument(
        '--config',
        type=str,
        default='config.yaml',
        help='Configuration file'
    )
    parser.add_argument(
        '--fdr',
        type=float,
        default=0.05,
        help='FDR threshold (default: 0.05)'
    )
    parser.add_argument(
        '--log2fc',
        type=float,
        default=1.0,
        help='Log2 fold change threshold (default: 1.0)'
    )
    
    args = parser.parse_args()
    
    print("="*70)
    print("LUNG CANCER DIFFERENTIAL GENE EXPRESSION ANALYSIS")
    print("="*70)
    print(f"Dataset: {args.dataset}")
    print(f"Started: {datetime.now()}")
    print("-"*70)
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    try:
        # 1. Load data
        print("\n📥 Loading dataset...")
        loader = GEODataLoader()
        expr_df, metadata = loader.load_dataset(args.dataset)
        print(f"✓ Loaded {expr_df.shape[0]} genes, {expr_df.shape[1]} samples")
        
        # 2. Preprocess data
        print("\n🔧 Preprocessing data...")
        preprocessor = DataPreprocessor()
        expr_filtered = preprocessor.filter_low_expression(expr_df)
        expr_normalized = preprocessor.normalize(expr_filtered)
        print(f"✓ Filtered to {expr_normalized.shape[0]} genes")
        
        # 3. Differential expression analysis
        print("\n📊 Performing differential expression analysis...")
        de_analyzer = DifferentialExpression(
            fdr_threshold=args.fdr,
            log2fc_threshold=args.log2fc
        )
        de_results = de_analyzer.analyze(expr_normalized, metadata)
        
        sig_genes = de_results[de_results['significant']].shape[0]
        print(f"✓ Found {sig_genes} significant DEGs")
        
        # 4. Generate visualizations
        print("\n📈 Creating visualizations...")
        visualizer = ExpressionVisualizer(output_dir=args.output)
        visualizer.create_all_plots(expr_normalized, de_results, metadata)
        print("✓ Plots saved to results/")
        
        # 5. Save results
        print("\n💾 Saving results...")
        de_results.to_csv(f'{args.output}/deg_results.csv', index=False)
        sig_only = de_results[de_results['significant']]
        sig_only.to_csv(f'{args.output}/significant_genes.csv', index=False)
        
        # Print summary
        print("\n" + "="*70)
        print("ANALYSIS COMPLETE")
        print("="*70)
        print(f"Total genes analyzed: {len(de_results):,}")
        print(f"Significant DEGs: {sig_genes:,}")
        print(f"Upregulated: {(de_results['direction'] == 'UP').sum():,}")
        print(f"Downregulated: {(de_results['direction'] == 'DOWN').sum():,}")
        print(f"\nResults saved to: {args.output}/")
        print(f"Completed: {datetime.now()}")
        
    except Exception as e:
        print(f"\n❌ Error: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
