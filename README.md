# genTile Shiny App

A web interface for designing CRISPR guides with interactive visualization and multi-modal selection.

## Overview

genTile Shiny provides a user-friendly web interface to design CRISPR guides for different applications through your browser. Upload gene lists or coordinates, select guide modes, and visualize results interactively.

## Features

- **Input Methods**: Gene symbols, Ensembl IDs, or genomic coordinates
- **Cell-Type Specific TSS**: Uses ENCODE CAGE data for 13+ cell lines including K562, HeLa, HepG2
- **Guide Selection Modes**:
  - **Tiling**: Evenly spaced guides across regions
  - **CRISPRi**: Knockdown guides near TSS (-50 to +300bp)
  - **CRISPRa**: Activation guides upstream of TSS (-400 to -50bp)
- **Filtering**: Remove guides with restriction sites or genomic variants
- **Visualization**: Interactive genome browser plots with guide positioning
- **Export**: Guide files (.txt) and genome browser format (.bed)

## Usage

### Input

**Gene Mode:**
```
TP53
BRCA1
ENSG00000141510
```

**Coordinates Mode:**
```
enhancer1,chr8:127736230,+
promoter2,chr17:7687546,-
```

### Parameters

- **Cell Line**: Choose from K562, HeLa, HepG2, etc. for CAGE-based TSS
- **Sequence Range**: Upstream (default: 1500bp) and downstream (default: 500bp) from TSS
- **Guide Selection**: Choose one or more modes (Tiling, CRISPRi, CRISPRa)
- **Filtering**: Optional restriction enzyme and K562 variant filtering

### Results

- **Interactive table** with guide details and scoring
- **Genome visualization** showing guide positions and target regions
- **Download options** for selected guides

## File Structure

```
genTile_shiny/
├── app.R                    # Main Shiny application
├── functions/
│   ├── pipeline_runner.R   # Backend pipeline execution  
│   └── plotting_functions.R # Visualization functions
└── run_app.R               # Command-line launcher
```

## Running the App

```r
# From R console
shiny::runApp("app.R")

# From command line  
Rscript run_app.R
```

## Requirements

- R with shiny, plotly, DT packages
- genTile pipeline installed and configured
- Reference genome (hg38) and FlashFry database

## Support

Report issues at: [GitHub Issues](https://github.com/nathanaelandrews/genTile_shiny/issues)
