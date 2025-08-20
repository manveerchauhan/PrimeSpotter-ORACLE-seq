#!/usr/bin/env python3
"""
Oracle Analysis Script

This script analyzes the output from the Oracle PrimeSpotter integration pipeline,
generating comparative plots and statistics for internal priming analysis.

Author: Oracle Analysis Pipeline
"""

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
import argparse
import logging
from scipy import stats
import json


class OracleAnalyzer:
    """Main class for Oracle analysis and visualization."""
    
    def __init__(self, results_dir):
        """Initialize analyzer with results directory."""
        self.results_dir = Path(results_dir)
        self.output_dir = self.results_dir / 'visualizations'
        self.setup_plotting()
        self.load_data()
    
    def setup_plotting(self):
        """Setup matplotlib and seaborn styling."""
        plt.style.use('default')
        sns.set_palette("husl")
        plt.rcParams['figure.figsize'] = (12, 8)
        plt.rcParams['font.size'] = 12
        plt.rcParams['axes.titlesize'] = 14
        plt.rcParams['axes.labelsize'] = 12
        plt.rcParams['xtick.labelsize'] = 10
        plt.rcParams['ytick.labelsize'] = 10
        plt.rcParams['legend.fontsize'] = 11
    
    def load_data(self):
        """Load summary statistics from pipeline output."""
        summary_file = self.results_dir / 'summary_statistics.csv'
        
        if not summary_file.exists():
            raise FileNotFoundError(f"Summary statistics not found: {summary_file}")
        
        self.data = pd.read_csv(summary_file)
        logging.info(f"Loaded data for {len(self.data)} samples")
        
        # Ensure required columns exist
        required_cols = ['dataset', 'depth', 'ip_rate', 'platform', 'sample_type']
        missing_cols = [col for col in required_cols if col not in self.data.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")
    
    def plot_platform_comparison(self):
        """Create comparative plots between ONT and PacBio platforms."""
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('Platform Comparison: ONT vs PacBio Internal Priming', fontsize=16)
        
        # 1. Box plot by platform
        ax1 = axes[0, 0]
        sns.boxplot(data=self.data, x='platform', y='ip_rate', ax=ax1)
        ax1.set_title('IP Rate Distribution by Platform')
        ax1.set_ylabel('Internal Priming Rate')
        
        # Add statistical test
        ont_rates = self.data[self.data['platform'] == 'ONT']['ip_rate']
        pb_rates = self.data[self.data['platform'] == 'PacBio']['ip_rate']
        if len(ont_rates) > 0 and len(pb_rates) > 0:
            stat, p_value = stats.mannwhitneyu(ont_rates, pb_rates, alternative='two-sided')
            ax1.text(0.5, 0.95, f'Mann-Whitney U p-value: {p_value:.4f}', 
                    transform=ax1.transAxes, ha='center', va='top',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        # 2. Violin plot by sample type and platform
        ax2 = axes[0, 1]
        sns.violinplot(data=self.data, x='sample_type', y='ip_rate', hue='platform', ax=ax2)
        ax2.set_title('IP Rate by Sample Type and Platform')
        ax2.set_ylabel('Internal Priming Rate')
        ax2.tick_params(axis='x', rotation=45)
        
        # 3. Scatter plot: depth vs IP rate
        ax3 = axes[1, 0]
        for platform in self.data['platform'].unique():
            platform_data = self.data[self.data['platform'] == platform]
            ax3.scatter(platform_data['depth'], platform_data['ip_rate'], 
                       label=platform, alpha=0.7, s=60)
        ax3.set_xlabel('Sequencing Depth')
        ax3.set_ylabel('Internal Priming Rate')
        ax3.set_title('IP Rate vs Sequencing Depth')
        ax3.set_xscale('log')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        # 4. Heatmap of IP rates
        ax4 = axes[1, 1]
        pivot_data = self.data.pivot_table(
            values='ip_rate', 
            index='sample_type', 
            columns='platform', 
            aggfunc='mean'
        )
        sns.heatmap(pivot_data, annot=True, fmt='.4f', cmap='YlOrRd', ax=ax4)
        ax4.set_title('Average IP Rate Heatmap')
        
        plt.tight_layout()
        
        # Save plot
        output_file = self.output_dir / 'ip_rates_comparison.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logging.info(f"Platform comparison plot saved: {output_file}")
    
    def plot_rarefaction_curves(self):
        """Create rarefaction curves showing IP rate vs sequencing depth."""
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('Rarefaction Analysis: IP Rate vs Sequencing Depth', fontsize=16)
        
        # 1. Overall rarefaction curves by platform
        ax1 = axes[0, 0]
        for platform in self.data['platform'].unique():
            platform_data = self.data[self.data['platform'] == platform]
            depth_grouped = platform_data.groupby('depth')['ip_rate'].agg(['mean', 'std']).reset_index()
            
            ax1.plot(depth_grouped['depth'], depth_grouped['mean'], 
                    marker='o', label=f'{platform} (mean)', linewidth=2)
            ax1.fill_between(depth_grouped['depth'], 
                           depth_grouped['mean'] - depth_grouped['std'],
                           depth_grouped['mean'] + depth_grouped['std'],
                           alpha=0.3)
        
        ax1.set_xlabel('Sequencing Depth')
        ax1.set_ylabel('Internal Priming Rate')
        ax1.set_title('Overall Rarefaction Curves')
        ax1.set_xscale('log')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # 2. Rarefaction by sample type
        ax2 = axes[0, 1]
        sample_types = self.data['sample_type'].unique()
        colors = plt.cm.Set1(np.linspace(0, 1, len(sample_types)))
        
        for i, sample_type in enumerate(sample_types):
            sample_data = self.data[self.data['sample_type'] == sample_type]
            depth_grouped = sample_data.groupby('depth')['ip_rate'].mean().reset_index()
            
            ax2.plot(depth_grouped['depth'], depth_grouped['ip_rate'], 
                    marker='s', label=sample_type, color=colors[i], linewidth=2)
        
        ax2.set_xlabel('Sequencing Depth')
        ax2.set_ylabel('Internal Priming Rate')
        ax2.set_title('Rarefaction by Sample Type')
        ax2.set_xscale('log')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # 3. Individual dataset trajectories
        ax3 = axes[1, 0]
        for dataset in self.data['dataset'].unique():
            dataset_data = self.data[self.data['dataset'] == dataset].sort_values('depth')
            if len(dataset_data) > 1:  # Only plot if multiple depths available
                ax3.plot(dataset_data['depth'], dataset_data['ip_rate'], 
                        marker='o', label=dataset, alpha=0.7, linewidth=1)
        
        ax3.set_xlabel('Sequencing Depth')
        ax3.set_ylabel('Internal Priming Rate')
        ax3.set_title('Individual Dataset Trajectories')
        ax3.set_xscale('log')
        ax3.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        ax3.grid(True, alpha=0.3)
        
        # 4. Coefficient of variation analysis
        ax4 = axes[1, 1]
        cv_data = []
        for dataset in self.data['dataset'].unique():
            dataset_data = self.data[self.data['dataset'] == dataset]
            if len(dataset_data) > 1:
                cv = dataset_data['ip_rate'].std() / dataset_data['ip_rate'].mean()
                platform = dataset_data['platform'].iloc[0]
                sample_type = dataset_data['sample_type'].iloc[0]
                cv_data.append({
                    'dataset': dataset,
                    'cv': cv,
                    'platform': platform,
                    'sample_type': sample_type
                })
        
        if cv_data:
            cv_df = pd.DataFrame(cv_data)
            sns.boxplot(data=cv_df, x='platform', y='cv', ax=ax4)
            ax4.set_title('Coefficient of Variation by Platform')
            ax4.set_ylabel('CV of IP Rate')
        
        plt.tight_layout()
        
        # Save plot
        output_file = self.output_dir / 'rarefaction_curves.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logging.info(f"Rarefaction curves saved: {output_file}")
    
    def plot_depth_saturation(self):
        """Analyze saturation effects at different sequencing depths."""
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        fig.suptitle('Sequencing Depth Saturation Analysis', fontsize=16)
        
        # 1. IP rate vs total reads
        ax1 = axes[0]
        for platform in self.data['platform'].unique():
            platform_data = self.data[self.data['platform'] == platform]
            ax1.scatter(platform_data['total_reads'], platform_data['ip_rate'], 
                       label=platform, alpha=0.7, s=60)
        
        ax1.set_xlabel('Total Reads')
        ax1.set_ylabel('Internal Priming Rate')
        ax1.set_title('IP Rate vs Total Reads')
        ax1.set_xscale('log')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # 2. Absolute IP read count vs depth
        ax2 = axes[1]
        for platform in self.data['platform'].unique():
            platform_data = self.data[self.data['platform'] == platform]
            ax2.scatter(platform_data['depth'], platform_data['ip_reads'], 
                       label=platform, alpha=0.7, s=60)
        
        ax2.set_xlabel('Sequencing Depth')
        ax2.set_ylabel('IP Read Count')
        ax2.set_title('IP Read Count vs Depth')
        ax2.set_xscale('log')
        ax2.set_yscale('log')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        # Save plot
        output_file = self.output_dir / 'depth_saturation.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logging.info(f"Depth saturation plot saved: {output_file}")
    
    def generate_statistical_report(self):
        """Generate comprehensive statistical report."""
        report_file = self.results_dir / 'statistical_report.txt'
        
        with open(report_file, 'w') as f:
            f.write("ORACLE PRIMESPOTTER STATISTICAL REPORT\n")
            f.write("=" * 50 + "\n\n")
            
            # Overall summary
            f.write("OVERALL SUMMARY\n")
            f.write("-" * 20 + "\n")
            f.write(f"Total samples analyzed: {len(self.data)}\n")
            f.write(f"Platforms: {', '.join(self.data['platform'].unique())}\n")
            f.write(f"Sample types: {', '.join(self.data['sample_type'].unique())}\n")
            f.write(f"Depth range: {self.data['depth'].min():,} - {self.data['depth'].max():,}\n")
            f.write(f"Overall IP rate: {self.data['ip_rate'].mean():.4f} ± {self.data['ip_rate'].std():.4f}\n\n")
            
            # Platform comparison
            f.write("PLATFORM COMPARISON\n")
            f.write("-" * 20 + "\n")
            platform_stats = self.data.groupby('platform')['ip_rate'].agg(['count', 'mean', 'std', 'min', 'max'])
            f.write(platform_stats.to_string())
            f.write("\n\n")
            
            # Statistical tests
            ont_rates = self.data[self.data['platform'] == 'ONT']['ip_rate']
            pb_rates = self.data[self.data['platform'] == 'PacBio']['ip_rate']
            
            if len(ont_rates) > 0 and len(pb_rates) > 0:
                stat, p_value = stats.mannwhitneyu(ont_rates, pb_rates, alternative='two-sided')
                f.write("STATISTICAL TESTS\n")
                f.write("-" * 17 + "\n")
                f.write(f"Mann-Whitney U test (ONT vs PacBio):\n")
                f.write(f"  Statistic: {stat:.2f}\n")
                f.write(f"  P-value: {p_value:.6f}\n")
                f.write(f"  Significant (α=0.05): {'Yes' if p_value < 0.05 else 'No'}\n\n")
            
            # Sample type analysis
            f.write("SAMPLE TYPE ANALYSIS\n")
            f.write("-" * 20 + "\n")
            sample_stats = self.data.groupby('sample_type')['ip_rate'].agg(['count', 'mean', 'std'])
            f.write(sample_stats.to_string())
            f.write("\n\n")
            
            # Depth analysis
            f.write("DEPTH ANALYSIS\n")
            f.write("-" * 15 + "\n")
            depth_stats = self.data.groupby('depth')['ip_rate'].agg(['count', 'mean', 'std'])
            f.write(depth_stats.to_string())
            f.write("\n\n")
            
            # Correlation analysis
            f.write("CORRELATION ANALYSIS\n")
            f.write("-" * 20 + "\n")
            corr_depth = stats.spearmanr(self.data['depth'], self.data['ip_rate'])
            corr_total = stats.spearmanr(self.data['total_reads'], self.data['ip_rate'])
            f.write(f"Spearman correlation (depth vs IP rate): r={corr_depth.correlation:.4f}, p={corr_depth.pvalue:.6f}\n")
            f.write(f"Spearman correlation (total reads vs IP rate): r={corr_total.correlation:.4f}, p={corr_total.pvalue:.6f}\n")
        
        logging.info(f"Statistical report saved: {report_file}")
    
    def run_analysis(self):
        """Run complete analysis pipeline."""
        logging.info("Starting Oracle analysis")
        
        # Generate all plots
        self.plot_platform_comparison()
        self.plot_rarefaction_curves()
        self.plot_depth_saturation()
        
        # Generate statistical report
        self.generate_statistical_report()
        
        logging.info("Analysis completed successfully")
        
        # Print summary to console
        self.print_summary()
    
    def print_summary(self):
        """Print analysis summary to console."""
        print("\n" + "=" * 60)
        print("ORACLE ANALYSIS SUMMARY")
        print("=" * 60)
        print(f"Samples analyzed: {len(self.data)}")
        print(f"Platforms: {', '.join(self.data['platform'].unique())}")
        print(f"Overall IP rate: {self.data['ip_rate'].mean():.4f} ± {self.data['ip_rate'].std():.4f}")
        
        print("\nPlatform comparison:")
        platform_summary = self.data.groupby('platform')['ip_rate'].agg(['mean', 'std'])
        for platform, stats in platform_summary.iterrows():
            print(f"  {platform}: {stats['mean']:.4f} ± {stats['std']:.4f}")
        
        print(f"\nVisualizations saved to: {self.output_dir}")
        print(f"Statistical report: {self.results_dir / 'statistical_report.txt'}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Oracle Analysis Script"
    )
    parser.add_argument(
        'results_dir',
        help='Directory containing Oracle pipeline results'
    )
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose logging'
    )
    
    args = parser.parse_args()
    
    # Setup logging
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    # Validate results directory
    if not os.path.exists(args.results_dir):
        print(f"Error: Results directory not found: {args.results_dir}")
        sys.exit(1)
    
    # Run analysis
    try:
        analyzer = OracleAnalyzer(args.results_dir)
        analyzer.run_analysis()
        
        print("\nAnalysis completed successfully!")
        
    except Exception as e:
        logging.error(f"Analysis failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()