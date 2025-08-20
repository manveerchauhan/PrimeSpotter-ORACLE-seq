# Oracle PrimeSpotter Integration Pipeline

A comprehensive pipeline for analyzing internal priming events across different long-read sequencing platforms and sample types using PrimeSpotter.

## Overview

This pipeline processes rarefied BAM files through PrimeSpotter to detect internal priming (IP) events, extracts statistics, and generates comparative analyses between ONT (Oxford Nanopore Technologies) and PacBio platforms across various sample types and sequencing depths.

## Features

- **Batch Processing**: Process multiple rarefied BAM files in parallel
- **Platform Comparison**: Compare IP rates between ONT and PacBio platforms
- **Sample Type Analysis**: Analyze different sample types (scRNA, snRNA, bulk, directRNA)
- **Rarefaction Analysis**: Study IP rate changes across sequencing depths
- **Statistical Reports**: Generate comprehensive statistical comparisons
- **Visualization**: Create publication-ready plots and figures
- **Modular Design**: Easy to extend and customize

## Requirements

### Dependencies

```bash
# Python packages
pip install pandas numpy matplotlib seaborn scipy tqdm pysam

# External tools
# - PrimeSpotter (https://github.com/yourusername/PrimeSpotter)
# - samtools (for BAM file handling)
```

### Input Data Structure

Your rarefied BAM files should follow this naming convention:
```
{DATASET}_{DEPTH}_rarefied.bam
```

Where:
- `DATASET`: ONT_scRNA, PB_scRNA, ONT_snRNA, PB_snRNA, ONT_bulk, PB_bulk, ONT_directRNA
- `DEPTH`: 1000, 5000, 10000, 50000, 100000 (or your chosen depths)

Example: `ONT_scRNA_1000_rarefied.bam`

## Quick Start

### 1. Setup Configuration

Copy and edit the configuration template:

```bash
cp oracle_config_template.json oracle_config.json
```

Edit `oracle_config.json` with your paths:

```json
{
    "primespotter_dir": "/path/to/PrimeSpotter/",
    "oracle_seq_dir": "/path/to/ORACLE-seq/",
    "genome_ref": "/path/to/reference.fa",
    "gtf_file": "/path/to/annotations.gtf",
    "bam_input_dir": "/path/to/rarefied_bams/",
    "output_dir": "/path/to/oracle_results/",
    "datasets": [
        "ONT_scRNA", "PB_scRNA", "ONT_snRNA", 
        "PB_snRNA", "ONT_bulk", "PB_bulk", "ONT_directRNA"
    ],
    "rarefaction_depths": [1000, 5000, 10000, 50000, 100000],
    "parallel_processes": 4,
    "primespotter_params": {
        "window_size": 10,
        "min_oligo_length": 18,
        "min_quality": 20
    }
}
```

### 2. Run the Pipeline

```bash
# Run the main integration pipeline
python oracle_primespotter_integration.py oracle_config.json

# Optional: enable verbose logging
python oracle_primespotter_integration.py oracle_config.json --verbose
```

### 3. Generate Analysis and Plots

```bash
# Run analysis on the results
python oracle_analysis.py /path/to/oracle_results/

# Optional: enable verbose logging
python oracle_analysis.py /path/to/oracle_results/ --verbose
```

## Output Structure

After running the pipeline, your output directory will contain:

```
oracle_results/
├── summary_statistics.csv           # IP rates for all datasets/depths
├── oracle_pipeline.log             # Detailed pipeline log
├── statistical_report.txt          # Comprehensive statistical analysis
├── read_ids/                       # Read ID files for each dataset/depth
│   ├── ONT_scRNA_1000_ip_reads.txt
│   ├── ONT_scRNA_1000_normal_reads.txt
│   ├── PB_scRNA_1000_ip_reads.txt
│   └── ... (for each dataset/depth combination)
├── gene_level_counts/              # Gene-level IP statistics (future feature)
├── visualizations/                 # Generated plots and figures
│   ├── ip_rates_comparison.png     # Platform comparison plots
│   ├── rarefaction_curves.png      # Rarefaction analysis
│   └── depth_saturation.png       # Saturation analysis
└── primespotter_output/           # Raw PrimeSpotter outputs
    ├── ONT_scRNA_1000_rarefied_annotated.sam
    └── ... (for each processed BAM)
```

## Key Files Description

### `oracle_primespotter_integration.py`
Main pipeline script that:
- Processes rarefied BAM files through PrimeSpotter
- Extracts IP statistics and read IDs
- Generates summary statistics
- Supports parallel processing

### `oracle_analysis.py`
Analysis script that:
- Creates comparative visualizations
- Performs statistical tests
- Generates comprehensive reports
- Produces publication-ready figures

### `oracle_config_template.json`
Configuration template with:
- All required path specifications
- PrimeSpotter parameter settings
- Processing options

## Usage Examples

### Basic Usage

```bash
# 1. Configure your paths
cp oracle_config_template.json my_config.json
# Edit my_config.json with your specific paths

# 2. Run the pipeline
python oracle_primespotter_integration.py my_config.json

# 3. Generate analysis
python oracle_analysis.py /path/to/your/oracle_results/
```

### Advanced Usage

```bash
# Run with specific configuration and verbose logging
python oracle_primespotter_integration.py my_config.json --verbose

# Analyze results with detailed logging
python oracle_analysis.py /path/to/results/ --verbose
```

### Processing Subset of Data

To process only specific datasets, modify the `datasets` list in your config file:

```json
{
    "datasets": ["ONT_scRNA", "PB_scRNA"],
    "rarefaction_depths": [1000, 10000, 100000]
}
```

## Configuration Options

### Core Paths
- `primespotter_dir`: Path to PrimeSpotter installation
- `genome_ref`: Reference genome FASTA file
- `gtf_file`: Gene annotation GTF file
- `bam_input_dir`: Directory containing rarefied BAM files
- `output_dir`: Output directory for results

### Processing Options
- `parallel_processes`: Number of parallel processes (default: 4)
- `datasets`: List of dataset types to process
- `rarefaction_depths`: List of depths to analyze

### PrimeSpotter Parameters
- `window_size`: Window size for priming detection (default: 10)
- `min_oligo_length`: Minimum oligo-dT length (default: 18)
- `min_quality`: Minimum base quality (default: 20)

## Expected Results

### Summary Statistics
The pipeline generates comprehensive statistics including:
- IP rates by platform and sample type
- Total read counts and IP read counts
- Statistical comparisons between platforms

### Visualizations
Generated plots include:
- Platform comparison box plots and violin plots
- Rarefaction curves showing IP rate vs depth
- Saturation analysis plots
- Heatmaps of IP rates across conditions

### Statistical Analysis
The analysis includes:
- Mann-Whitney U tests for platform comparisons
- Spearman correlations for depth relationships
- Coefficient of variation analysis
- Comprehensive statistical reporting

## Troubleshooting

### Common Issues

1. **PrimeSpotter not found**
   ```
   Error: PrimeSpotter.py not found in specified directory
   ```
   - Verify `primespotter_dir` path in config
   - Ensure PrimeSpotter is properly installed

2. **BAM files not found**
   ```
   Warning: BAM file not found: /path/to/file.bam
   ```
   - Check BAM file naming convention
   - Verify `bam_input_dir` path
   - Ensure BAM files are indexed

3. **Missing dependencies**
   ```
   ModuleNotFoundError: No module named 'pysam'
   ```
   - Install required Python packages:
     ```bash
     pip install pandas numpy matplotlib seaborn scipy tqdm pysam
     ```

4. **Memory issues with large files**
   - Reduce `parallel_processes` in config
   - Process datasets individually
   - Ensure sufficient system memory

### Log Files

Check the pipeline log for detailed information:
```bash
tail -f oracle_results/oracle_pipeline.log
```

## Contributing

To extend the pipeline:

1. **Add new sample types**: Modify the `datasets` list in config
2. **Add new analysis**: Extend the `OracleAnalyzer` class
3. **Customize plots**: Modify plotting functions in `oracle_analysis.py`
4. **Add new statistics**: Extend the statistical reporting functions

## Citation

If you use this pipeline in your research, please cite:
- PrimeSpotter: [Citation for PrimeSpotter]
- ORACLE-seq: [Citation for ORACLE-seq methodology]

## Support

For questions or issues:
1. Check the troubleshooting section above
2. Review log files for detailed error messages
3. Ensure all dependencies are properly installed
4. Verify input file formats and naming conventions

## License

This pipeline is provided under the same license as the associated research project.