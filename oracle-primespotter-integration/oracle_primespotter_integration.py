#!/usr/bin/env python3
"""
Oracle PrimeSpotter Integration Pipeline

This script processes rarefied BAM files through PrimeSpotter to detect internal priming
events and generates comprehensive statistics for comparison across different sequencing
platforms and depths.

Author: Oracle Analysis Pipeline
"""

import os
import sys
import json
import subprocess
import logging
import pandas as pd
import pysam
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
from collections import defaultdict
import argparse


class OraclePrimeSpotterPipeline:
    """Main pipeline class for Oracle PrimeSpotter integration."""
    
    def __init__(self, config_file):
        """Initialize pipeline with configuration."""
        self.config = self.load_config(config_file)
        self.setup_logging()
        self.setup_directories()
        
    def load_config(self, config_file):
        """Load and validate configuration file."""
        try:
            with open(config_file, 'r') as f:
                config = json.load(f)
            
            # Validate required fields
            required_fields = ['primespotter_dir', 'genome_ref', 'gtf_file', 
                             'bam_input_dir', 'output_dir']
            for field in required_fields:
                if field not in config or config[field].startswith('[PLACEHOLDER'):
                    raise ValueError(f"Please configure {field} in {config_file}")
            
            return config
        except Exception as e:
            logging.error(f"Error loading config: {e}")
            sys.exit(1)
    
    def setup_logging(self):
        """Setup logging configuration."""
        log_file = os.path.join(self.config['output_dir'], 'oracle_pipeline.log')
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler()
            ]
        )
    
    def setup_directories(self):
        """Create output directory structure."""
        base_dir = Path(self.config['output_dir'])
        subdirs = ['read_ids', 'gene_level_counts', 'visualizations', 'primespotter_output']
        
        for subdir in subdirs:
            (base_dir / subdir).mkdir(parents=True, exist_ok=True)
    
    def find_bam_files(self):
        """Find all rarefied BAM files matching expected patterns."""
        bam_dir = Path(self.config['bam_input_dir'])
        bam_files = []
        
        for dataset in self.config['datasets']:
            for depth in self.config['rarefaction_depths']:
                pattern = f"{dataset}_{depth}_rarefied.bam"
                bam_path = bam_dir / pattern
                if bam_path.exists():
                    bam_files.append({
                        'path': str(bam_path),
                        'dataset': dataset,
                        'depth': depth,
                        'basename': pattern.replace('.bam', '')
                    })
                else:
                    logging.warning(f"BAM file not found: {bam_path}")
        
        logging.info(f"Found {len(bam_files)} BAM files to process")
        return bam_files
    
    def run_primespotter(self, bam_info):
        """Run PrimeSpotter on a single BAM file."""
        try:
            bam_path = bam_info['path']
            output_prefix = os.path.join(
                self.config['output_dir'], 'primespotter_output', 
                bam_info['basename']
            )
            
            # Construct PrimeSpotter command
            primespotter_script = os.path.join(
                self.config['primespotter_dir'], 'PrimeSpotter.py'
            )
            
            cmd = [
                'python3', primespotter_script,
                '-i', bam_path,
                '-g', self.config['genome_ref'],
                '-a', self.config['gtf_file'],
                '-o', output_prefix,
                '--window_size', str(self.config['primespotter_params']['window_size']),
                '--min_oligo_length', str(self.config['primespotter_params']['min_oligo_length']),
                '--min_quality', str(self.config['primespotter_params']['min_quality'])
            ]
            
            logging.info(f"Running PrimeSpotter on {bam_info['basename']}")
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            return {
                'bam_info': bam_info,
                'output_prefix': output_prefix,
                'success': True,
                'stdout': result.stdout,
                'stderr': result.stderr
            }
            
        except subprocess.CalledProcessError as e:
            logging.error(f"PrimeSpotter failed for {bam_info['basename']}: {e}")
            return {
                'bam_info': bam_info,
                'success': False,
                'error': str(e)
            }
        except Exception as e:
            logging.error(f"Unexpected error for {bam_info['basename']}: {e}")
            return {
                'bam_info': bam_info,
                'success': False,
                'error': str(e)
            }
    
    def parse_primespotter_output(self, result):
        """Parse PrimeSpotter SAM output to extract IP statistics and read IDs."""
        if not result['success']:
            return None
        
        bam_info = result['bam_info']
        sam_file = f"{result['output_prefix']}_annotated.sam"
        
        if not os.path.exists(sam_file):
            logging.error(f"PrimeSpotter output not found: {sam_file}")
            return None
        
        try:
            ip_reads = []
            normal_reads = []
            total_reads = 0
            ip_count = 0
            
            with pysam.AlignmentFile(sam_file, 'r') as sam:
                for read in sam:
                    if read.is_unmapped:
                        continue
                    
                    total_reads += 1
                    
                    # Check for IP tag
                    ip_tag = read.get_tag('IP') if read.has_tag('IP') else 'F'
                    
                    if ip_tag == 'T':
                        ip_reads.append(read.query_name)
                        ip_count += 1
                    else:
                        normal_reads.append(read.query_name)
            
            # Calculate IP rate
            ip_rate = ip_count / total_reads if total_reads > 0 else 0
            
            # Save read IDs
            self.save_read_ids(bam_info, ip_reads, normal_reads)
            
            return {
                'dataset': bam_info['dataset'],
                'depth': bam_info['depth'],
                'total_reads': total_reads,
                'ip_reads': ip_count,
                'normal_reads': len(normal_reads),
                'ip_rate': ip_rate,
                'ip_read_ids': ip_reads,
                'normal_read_ids': normal_reads
            }
            
        except Exception as e:
            logging.error(f"Error parsing {sam_file}: {e}")
            return None
    
    def save_read_ids(self, bam_info, ip_reads, normal_reads):
        """Save read IDs to text files."""
        base_name = f"{bam_info['dataset']}_{bam_info['depth']}"
        read_ids_dir = Path(self.config['output_dir']) / 'read_ids'
        
        # Save IP reads
        ip_file = read_ids_dir / f"{base_name}_ip_reads.txt"
        with open(ip_file, 'w') as f:
            for read_id in ip_reads:
                f.write(f"{read_id}\n")
        
        # Save normal reads
        normal_file = read_ids_dir / f"{base_name}_normal_reads.txt"
        with open(normal_file, 'w') as f:
            for read_id in normal_reads:
                f.write(f"{read_id}\n")
        
        logging.info(f"Saved {len(ip_reads)} IP reads and {len(normal_reads)} normal reads for {base_name}")
    
    def process_all_bams(self):
        """Process all BAM files through PrimeSpotter pipeline."""
        bam_files = self.find_bam_files()
        
        if not bam_files:
            logging.error("No BAM files found to process")
            return []
        
        results = []
        
        # Use parallel processing if configured
        if self.config.get('parallel_processes', 1) > 1:
            with ProcessPoolExecutor(max_workers=self.config['parallel_processes']) as executor:
                # Submit all jobs
                future_to_bam = {
                    executor.submit(self.run_primespotter, bam_info): bam_info 
                    for bam_info in bam_files
                }
                
                # Process completed jobs with progress bar
                for future in tqdm(as_completed(future_to_bam), 
                                 total=len(bam_files), 
                                 desc="Processing BAM files"):
                    result = future.result()
                    parsed_result = self.parse_primespotter_output(result)
                    if parsed_result:
                        results.append(parsed_result)
        else:
            # Sequential processing
            for bam_info in tqdm(bam_files, desc="Processing BAM files"):
                result = self.run_primespotter(bam_info)
                parsed_result = self.parse_primespotter_output(result)
                if parsed_result:
                    results.append(parsed_result)
        
        return results
    
    def generate_summary_statistics(self, results):
        """Generate summary statistics CSV file."""
        if not results:
            logging.warning("No results to summarize")
            return
        
        # Convert to DataFrame
        df = pd.DataFrame(results)
        
        # Add platform column
        df['platform'] = df['dataset'].apply(lambda x: 'ONT' if x.startswith('ONT') else 'PacBio')
        df['sample_type'] = df['dataset'].apply(lambda x: x.split('_')[1] if '_' in x else 'unknown')
        
        # Save summary
        summary_file = os.path.join(self.config['output_dir'], 'summary_statistics.csv')
        df.to_csv(summary_file, index=False)
        
        logging.info(f"Summary statistics saved to {summary_file}")
        
        # Print summary to console
        print("\n=== ORACLE PRIMESPOTTER SUMMARY ===")
        print(f"Total datasets processed: {len(df)}")
        print(f"Average IP rate: {df['ip_rate'].mean():.4f}")
        print("\nBy Platform:")
        platform_summary = df.groupby('platform')['ip_rate'].agg(['mean', 'std', 'count'])
        print(platform_summary)
        print("\nBy Sample Type:")
        sample_summary = df.groupby('sample_type')['ip_rate'].agg(['mean', 'std', 'count'])
        print(sample_summary)
    
    def run_pipeline(self):
        """Execute the complete pipeline."""
        logging.info("Starting Oracle PrimeSpotter integration pipeline")
        
        # Process all BAM files
        results = self.process_all_bams()
        
        # Generate summary statistics
        self.generate_summary_statistics(results)
        
        logging.info("Pipeline completed successfully")
        return results


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Oracle PrimeSpotter Integration Pipeline"
    )
    parser.add_argument(
        'config', 
        help='JSON configuration file'
    )
    parser.add_argument(
        '--verbose', '-v', 
        action='store_true',
        help='Enable verbose logging'
    )
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Validate config file exists
    if not os.path.exists(args.config):
        print(f"Error: Configuration file not found: {args.config}")
        sys.exit(1)
    
    # Run pipeline
    try:
        pipeline = OraclePrimeSpotterPipeline(args.config)
        results = pipeline.run_pipeline()
        
        print(f"\nPipeline completed successfully!")
        print(f"Results saved to: {pipeline.config['output_dir']}")
        
    except Exception as e:
        logging.error(f"Pipeline failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()