# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Educational tutorial for glycoproteomic (proteomics with post-translational modifications) data analysis using R and Python. The project teaches bioinformatics workflows for processing, normalizing, analyzing, and visualizing mass spectrometry proteomics data.

## Project Structure

- **Chapter 1**: R basics and package introduction
- **Chapter 2**: Data normalization (sample loading, TMM) with UMAP visualization
- **Chapter 3**: Differential expression analysis using limma
- **Chapter 4**: Gene ontology and pathway enrichment with clusterProfiler
- **Chapter 5**: Structural analysis using localCIDER and StructureMap

Each chapter contains:
- `.Rmd` files (R Markdown source - primary code)
- `.md` files (rendered output for GitHub)
- `training_data/` folder with example datasets

## Running Code

**R Markdown files**: Open in RStudio and knit, or use:
```r
rmarkdown::render("Chapter_X/file.Rmd")
```

**Python scripts**: Called from R via reticulate:
```r
library(reticulate)
use_condaenv(condaenv = '/opt/anaconda3/envs/structure_analysis', required = TRUE)
source_python("script.py")
```

## Key R Packages

- **Data manipulation**: tidyverse (dplyr, tidyr, purrr, stringr)
- **Data I/O**: readxl, writexl, readr
- **Visualization**: ggplot2, ggpubr, ComplexHeatmap
- **Statistics**: rstatix, limma, edgeR
- **Bioinformatics**: clusterProfiler, org.Hs.eg.db, AnnotationDbi
- **R-Python bridge**: reticulate

## Key Python Packages

- Core: pandas, numpy, matplotlib, seaborn
- Dimensionality reduction: umap
- Structural analysis (Chapter 5): localcider, structuremap

## Code Style

- Uses tidyverse conventions with pipe operators (`|>`)
- 2-space indentation (configured in .Rproj)
- UTF-8 encoding
- Data follows tidy data principles (each variable a column, each observation a row)

## Data Flow Pattern

1. Raw mass spectrometry results (CSV/TSV from search software)
2. PSM filtering and quality control (XCorr, PPM thresholds)
3. Protein-level aggregation
4. Normalization (Sample Loading → TMM)
5. Statistical analysis (limma for differential expression)
6. Visualization and enrichment analysis

## Important Notes

- Large data files (.csv, .xlsx, .tsv) are gitignored - training data must be downloaded separately or generated
- Python integration requires conda environment setup for Chapter 5
- Column names in proteomics data often contain special characters - use `name_repair = 'universal'` when reading
