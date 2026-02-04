# Glycoproteomic Data Analysis using R and Python

A comprehensive tutorial for bioinformaticians and proteomics researchers to learn glycoproteomic (proteomics with post-translational modifications) data analysis workflows using R and Python.

## Target Audience

This tutorial is designed for:
- Graduate students and researchers in proteomics/glycoproteomics
- Bioinformaticians learning mass spectrometry data analysis
- Scientists with basic R/Python knowledge wanting to analyze TMT-based quantitative proteomics data

**Prerequisites:**
- Basic familiarity with R and/or Python
- Understanding of proteomics concepts (proteins, peptides, mass spectrometry)
- Familiarity with statistical concepts (p-values, fold changes, normalization)

## Tutorial Structure

| Chapter | Topic | Learning Outcomes |
|---------|-------|-------------------|
| **Chapter 1** | [R Basics](Chapter_1_R_Basics/R_Basics.md) | Install R/RStudio, understand tidyverse, perform basic statistical tests (K-S test, Wilcoxon), create publication-quality plots |
| **Chapter 2** | [Data Normalization](Chapter_2_Data_Normalization/Data_Normalization.md) | Filter PSMs, aggregate to protein level, apply sample loading and TMM normalization, visualize with UMAP |
| **Chapter 3** | [Differential Analysis](Chapter_3_Differential_Analysis/Differential_Expression_Analysis.md) | Set up limma design matrices, perform differential expression analysis, interpret results |
| **Chapter 4** | [Enrichment Analysis](Chapter_4_Enrichment_Analysis/Enrichment_Analysis.md) | Run GSEA with GO/KEGG, analyze protein complexes (CORUM) and domains (Pfam) |
| **Chapter 5** | [Structure Analysis](Chapter_5_Structure_Related_Analysis/Structure_Related_Analysis.md) | Analyze peptide physicochemical properties, integrate AlphaFold structural data |

## Installation & Setup

### R Environment

1. **Install R** (version 4.3+): https://cran.r-project.org/
2. **Install RStudio**: https://posit.co/download/rstudio-desktop/
3. **Install required packages:**

```r
# CRAN packages
install.packages(c(
  "tidyverse", "readxl", "writexl", "readr",
  "showtext", "rstatix", "ggpubr", "reticulate"
))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "limma", "edgeR", "clusterProfiler",
  "org.Hs.eg.db", "AnnotationDbi", "ComplexHeatmap"
))
```

### Python Environment (for Chapters 2 & 5)

```bash
# Create conda environment for UMAP (Chapter 2)
conda create -n UMAP_env python=3.9
conda activate UMAP_env
pip install umap-learn pandas numpy matplotlib seaborn

# Create conda environment for structure analysis (Chapter 5)
conda create -n structure_analysis python=3.9
conda activate structure_analysis
pip install pandas numpy matplotlib seaborn localcider structuremap
```

### Platform-Specific Notes

- **Windows**: Ensure conda is added to PATH during installation
- **macOS**: If using Apple Silicon, some packages may require Rosetta 2
- **Linux**: May require additional system libraries for some R packages

## Quick Start

1. Clone this repository:
```bash
git clone https://github.com/lfu46/Glycoproteomic-Data-Analysis-using-R-and-Python.git
cd Glycoproteomic-Data-Analysis-using-R-and-Python
```

2. Open the `.Rproj` file in RStudio

3. Start with Chapter 1 to verify your R setup works:
```r
rmarkdown::render("Chapter_1_R_Basics/R_Basics.Rmd")
```

4. Each chapter builds on previous ones, so work through them sequentially.

## Data Access

See [DATA_ACCESS.md](DATA_ACCESS.md) for detailed information about:
- Sample datasets included for quick testing
- Where to download full datasets
- Expected data formats and column descriptions

## Key R Packages Used

| Category | Packages |
|----------|----------|
| Data manipulation | tidyverse (dplyr, tidyr, purrr, stringr) |
| Data I/O | readxl, writexl, readr |
| Visualization | ggplot2, ggpubr, ComplexHeatmap |
| Statistics | rstatix, limma, edgeR |
| Bioinformatics | clusterProfiler, org.Hs.eg.db, AnnotationDbi |
| R-Python bridge | reticulate |

## Code Style

This tutorial follows tidyverse conventions:
- Pipe operator (`|>`) for chaining operations
- Tidy data principles (each variable a column, each observation a row)
- 2-space indentation
- UTF-8 encoding

## Typical Data Flow

```
Raw MS Results (CSV/TSV from search software)
    ↓
PSM Filtering (XCorr, PPM thresholds)
    ↓
Protein-level Aggregation
    ↓
Normalization (Sample Loading → TMM)
    ↓
Statistical Analysis (limma)
    ↓
Visualization & Enrichment Analysis
```

## Citation

If you use this tutorial in your research, please cite:

```
Fu, L. (2025). Glycoproteomic Data Analysis using R and Python: A Practical Tutorial.
GitHub: https://github.com/lfu46/Glycoproteomic-Data-Analysis-using-R-and-Python
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Support

- **Issues**: Report bugs or request features via [GitHub Issues](https://github.com/lfu46/Glycoproteomic-Data-Analysis-using-R-and-Python/issues)
- **Questions**: For questions about the tutorial content, open a GitHub Discussion

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
