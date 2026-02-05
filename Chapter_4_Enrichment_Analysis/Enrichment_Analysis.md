Chapter 4: Enrichment Analysis
================
Longping Fu
03/05/2025

``` r
library(tidyverse)
```

    ## Warning: package 'purrr' was built under R version 4.4.1

    ## Warning: package 'lubridate' was built under R version 4.4.1

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.4     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(clusterProfiler)
```

    ## Warning: package 'clusterProfiler' was built under R version 4.4.2

    ## 
    ## clusterProfiler v4.14.4 Learn more at https://yulab-smu.top/contribution-knowledge-mining/
    ## 
    ## Please cite:
    ## 
    ## Guangchuang Yu, Li-Gen Wang, Yanyan Han and Qing-Yu He.
    ## clusterProfiler: an R package for comparing biological themes among
    ## gene clusters. OMICS: A Journal of Integrative Biology. 2012,
    ## 16(5):284-287
    ## 
    ## Attaching package: 'clusterProfiler'
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     simplify
    ## 
    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

``` r
library(org.Hs.eg.db)
```

    ## Loading required package: AnnotationDbi

    ## Warning: package 'AnnotationDbi' was built under R version 4.4.1

    ## Loading required package: stats4
    ## Loading required package: BiocGenerics

    ## Warning: package 'BiocGenerics' was built under R version 4.4.1

    ## 
    ## Attaching package: 'BiocGenerics'
    ## 
    ## The following objects are masked from 'package:lubridate':
    ## 
    ##     intersect, setdiff, union
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union
    ## 
    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, saveRDS, setdiff,
    ##     table, tapply, union, unique, unsplit, which.max, which.min
    ## 
    ## Loading required package: Biobase

    ## Warning: package 'Biobase' was built under R version 4.4.1

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.
    ## 
    ## Loading required package: IRanges

    ## Warning: package 'IRanges' was built under R version 4.4.2

    ## Loading required package: S4Vectors

    ## Warning: package 'S4Vectors' was built under R version 4.4.1

    ## 
    ## Attaching package: 'S4Vectors'
    ## 
    ## The following object is masked from 'package:clusterProfiler':
    ## 
    ##     rename
    ## 
    ## The following objects are masked from 'package:lubridate':
    ## 
    ##     second, second<-
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename
    ## 
    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand
    ## 
    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname
    ## 
    ## 
    ## Attaching package: 'IRanges'
    ## 
    ## The following object is masked from 'package:clusterProfiler':
    ## 
    ##     slice
    ## 
    ## The following object is masked from 'package:lubridate':
    ## 
    ##     %within%
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     reduce
    ## 
    ## 
    ## Attaching package: 'AnnotationDbi'
    ## 
    ## The following object is masked from 'package:clusterProfiler':
    ## 
    ##     select
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

``` r
# import differential analysis result
Iron_protein_toptable_tb <- read_csv(
  'training_data/Iron_protein_toptable_tb.csv'
)
```

    ## Rows: 7107 Columns: 7
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (1): UniProt_Accession
    ## dbl (6): logFC, AveExpr, t, P.Value, adj.P.Val, B
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# generate ranked gene list
# tutorial: https://yulab-smu.top/biomedical-knowledge-mining-book/faq.html
# How to prepare your own geneList
genelist <- Iron_protein_toptable_tb$logFC
names(genelist) <- Iron_protein_toptable_tb$UniProt_Accession
genelist <- sort(genelist, decreasing = TRUE)

# GSEA for GO
Iron_GO_GSEA <- gseGO(
  geneList = genelist,
  # 'ALL' for all 'BP', 'MF', and 'CC'
  ont = 'ALL',
  OrgDb = org.Hs.eg.db,
  # check the manual for org.Hs.eg.db for all available identifiers
  # https://bioconductor.org/packages/release/data/annotation/manuals/org.Hs.eg.db/man/org.Hs.eg.db.pdf
  keyType = 'UNIPROT',
  pvalueCutoff = 1
)
```

    ## using 'fgsea' for GSEA analysis, please cite Korotkevich et al (2019).
    ## 
    ## preparing geneSet collections...
    ## GSEA analysis...

    ## Warning in fgseaMultilevel(pathways = pathways, stats = stats, minSize =
    ## minSize, : For some of the pathways the P-values were likely overestimated. For
    ## such pathways log2err is set to NA.

    ## Warning in fgseaMultilevel(pathways = pathways, stats = stats, minSize =
    ## minSize, : For some pathways, in reality P-values are less than 1e-10. You can
    ## set the `eps` argument to zero for better estimation.

    ## leading edge analysis...
    ## done...

``` r
# check GSEA GO result
Iron_GO_GSEA@result |> 
  as_tibble() |> 
  arrange(p.adjust)
```

    ## # A tibble: 5,436 × 12
    ##    ONTOLOGY ID         Description setSize enrichmentScore   NES pvalue p.adjust
    ##    <chr>    <chr>      <chr>         <int>           <dbl> <dbl>  <dbl>    <dbl>
    ##  1 CC       GO:0022626 cytosolic …     104           0.652  2.79  1e-10  2.09e-8
    ##  2 BP       GO:0002181 cytoplasmi…     145           0.616  2.78  1e-10  2.09e-8
    ##  3 CC       GO:1904813 ficolin-1-…      98           0.605  2.57  1e-10  2.09e-8
    ##  4 CC       GO:0101002 ficolin-1-…     123           0.570  2.50  1e-10  2.09e-8
    ##  5 CC       GO:0034774 secretory …     185           0.534  2.50  1e-10  2.09e-8
    ##  6 CC       GO:0060205 cytoplasmi…     187           0.529  2.48  1e-10  2.09e-8
    ##  7 CC       GO:0031983 vesicle lu…     188           0.525  2.46  1e-10  2.09e-8
    ##  8 MF       GO:0003735 structural…     148           0.527  2.38  1e-10  2.09e-8
    ##  9 MF       GO:0045296 cadherin b…     273           0.479  2.35  1e-10  2.09e-8
    ## 10 CC       GO:0044391 ribosomal …     169           0.499  2.31  1e-10  2.09e-8
    ## # ℹ 5,426 more rows
    ## # ℹ 4 more variables: qvalue <dbl>, rank <dbl>, leading_edge <chr>,
    ## #   core_enrichment <chr>

``` r
# GSEA for KEGG
Iron_KEGG_GSEA <- gseKEGG(
  geneList = genelist,
  # 'hsa' for Homo sapiens (human)
  organism = 'hsa',
  # 'uniprot' for UniProt_Accession
  keyType = 'uniprot',
  pvalueCutoff = 1
)
```

    ## Reading KEGG annotation online: "https://rest.kegg.jp/link/hsa/pathway"...
    ## Reading KEGG annotation online: "https://rest.kegg.jp/list/pathway/hsa"...
    ## Reading KEGG annotation online: "https://rest.kegg.jp/conv/uniprot/hsa"...
    ## using 'fgsea' for GSEA analysis, please cite Korotkevich et al (2019).
    ## 
    ## preparing geneSet collections...
    ## GSEA analysis...

    ## Warning in fgseaMultilevel(pathways = pathways, stats = stats, minSize =
    ## minSize, : For some of the pathways the P-values were likely overestimated. For
    ## such pathways log2err is set to NA.
    ## Warning in fgseaMultilevel(pathways = pathways, stats = stats, minSize =
    ## minSize, : For some pathways, in reality P-values are less than 1e-10. You can
    ## set the `eps` argument to zero for better estimation.

    ## leading edge analysis...
    ## done...

``` r
# check GSEA KEGG result
Iron_KEGG_GSEA@result |> 
  as_tibble() |> 
  arrange(p.adjust)
```

    ## # A tibble: 308 × 11
    ##    ID       Description  setSize enrichmentScore   NES   pvalue p.adjust  qvalue
    ##    <chr>    <chr>          <int>           <dbl> <dbl>    <dbl>    <dbl>   <dbl>
    ##  1 hsa03010 Ribosome         124           0.582  2.61 1   e-10  7.70e-9 5.71e-9
    ##  2 hsa05012 Parkinson d…     188           0.472  2.22 1   e-10  7.70e-9 5.71e-9
    ##  3 hsa05016 Huntington …     216           0.454  2.19 1   e-10  7.70e-9 5.71e-9
    ##  4 hsa05022 Pathways of…     291           0.421  2.10 1   e-10  7.70e-9 5.71e-9
    ##  5 hsa05014 Amyotrophic…     253           0.431  2.12 3.15e-10  1.94e-8 1.44e-8
    ##  6 hsa05171 Coronavirus…     134           0.504  2.28 1.12e- 9  5.77e-8 4.28e-8
    ##  7 hsa03050 Proteasome        41           0.708  2.57 1.96e- 9  8.61e-8 6.39e-8
    ##  8 hsa05020 Prion disea…     185           0.443  2.09 7.60e- 9  2.93e-7 2.17e-7
    ##  9 hsa01200 Carbon meta…      90           0.528  2.25 4.99e- 8  1.71e-6 1.27e-6
    ## 10 hsa05010 Alzheimer d…     243           0.386  1.88 3.17e- 7  9.77e-6 7.25e-6
    ## # ℹ 298 more rows
    ## # ℹ 3 more variables: rank <dbl>, leading_edge <chr>, core_enrichment <chr>

``` r
# download CORUM database
# https://mips.helmholtz-muenchen.de/corum/download
# UniProt-CORUM Mapping, Corum 5.1 release (2025-01-21)
corum_human <- read_delim(
  'training_data/corum_uniprotCorumMapping.txt',
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  mutate(
    corum_id = paste('corum_id', corum_id, sep = '_')
  ) |> 
  dplyr::select(corum_id, UniProtKB_accession_number)
```

    ## Rows: 25877 Columns: 2
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): UniProtKB_accession_number
    ## dbl (1): corum_id
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# one column for TERM, one column for GENE
# for proteomic result, most of the time GENE refers to UniProt_Accession
colnames(corum_human) <- c('TERM', 'GENE')

# download Pfam database
# https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.0/proteomes/
# 9606.tsv.gz, 2024-05-28 13:25 2.8M
pfam_human <- read_tsv(
  'training_data/9606.tsv',
  skip = 3,
  col_names = FALSE
) |> 
  dplyr::select(X6, X1)
```

    ## Rows: 128058 Columns: 14
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (5): X1, X6, X7, X8, X14
    ## dbl (9): X2, X3, X4, X5, X9, X10, X11, X12, X13
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# one column for TERM, one column for GENE
colnames(pfam_human) <- c('TERM', 'GENE')

# GSEA for protein complex
Iron_CORUM_GSEA <- GSEA(
  geneList = genelist,
  pvalueCutoff = 1,
  # parameter TERM2GENE, one column for TERM, one column for GENE
  TERM2GENE = corum_human
)
```

    ## using 'fgsea' for GSEA analysis, please cite Korotkevich et al (2019).
    ## 
    ## preparing geneSet collections...
    ## GSEA analysis...

    ## Warning in fgseaMultilevel(pathways = pathways, stats = stats, minSize =
    ## minSize, : For some of the pathways the P-values were likely overestimated. For
    ## such pathways log2err is set to NA.

    ## Warning in fgseaMultilevel(pathways = pathways, stats = stats, minSize =
    ## minSize, : For some pathways, in reality P-values are less than 1e-10. You can
    ## set the `eps` argument to zero for better estimation.

    ## leading edge analysis...
    ## done...

``` r
# check GSEA CORUM result
Iron_CORUM_GSEA@result |> 
  as_tibble() |> 
  arrange(p.adjust)
```

    ## # A tibble: 176 × 11
    ##    ID        Description setSize enrichmentScore   NES   pvalue p.adjust  qvalue
    ##    <chr>     <chr>         <int>           <dbl> <dbl>    <dbl>    <dbl>   <dbl>
    ##  1 corum_id… corum_id_3…      76           0.689  2.81 1   e-10  8.80e-9 7.37e-9
    ##  2 corum_id… corum_id_3…     100           0.649  2.79 1   e-10  8.80e-9 7.37e-9
    ##  3 corum_id… corum_id_1…      36           0.748  2.62 1.85e-10  1.09e-8 9.09e-9
    ##  4 corum_id… corum_id_3…      45           0.696  2.56 1.08e- 9  4.76e-8 3.98e-8
    ##  5 corum_id… corum_id_1…      22           0.765  2.44 2.16e- 7  7.62e-6 6.38e-6
    ##  6 corum_id… corum_id_3…      31           0.689  2.36 5.40e- 7  1.36e-5 1.14e-5
    ##  7 corum_id… corum_id_3…      31           0.689  2.36 5.40e- 7  1.36e-5 1.14e-5
    ##  8 corum_id… corum_id_1…      16           0.792  2.32 4.69e- 6  1.03e-4 8.64e-5
    ##  9 corum_id… corum_id_32      20           0.743  2.30 5.67e- 6  1.11e-4 9.29e-5
    ## 10 corum_id… corum_id_5…      14           0.816  2.32 6.70e- 6  1.18e-4 9.87e-5
    ## # ℹ 166 more rows
    ## # ℹ 3 more variables: rank <dbl>, leading_edge <chr>, core_enrichment <chr>

``` r
# GSEA for protein domain
Iron_Pfam_GSEA <- GSEA(
  geneList = genelist,
  pvalueCutoff = 1,
  # parameter TERM2GENE, one column for TERM, one column for GENE
  TERM2GENE = pfam_human
)
```

    ## using 'fgsea' for GSEA analysis, please cite Korotkevich et al (2019).
    ## 
    ## preparing geneSet collections...
    ## GSEA analysis...
    ## leading edge analysis...
    ## done...

``` r
# check GSEA Pfam result
Iron_Pfam_GSEA@result |> 
  as_tibble() |> 
  arrange(p.adjust)
```

    ## # A tibble: 174 × 11
    ##    ID      Description setSize enrichmentScore   NES   pvalue   p.adjust  qvalue
    ##    <chr>   <chr>         <int>           <dbl> <dbl>    <dbl>      <dbl>   <dbl>
    ##  1 PF00071 PF00071          76           0.605  2.45 6.11e-10    1.06e-7 8.48e-8
    ##  2 PF00096 PF00096          90          -0.462 -2.04 2.33e- 6    2.03e-4 1.62e-4
    ##  3 PF00025 PF00025          15           0.768  2.20 4.97e- 5    2.88e-3 2.30e-3
    ##  4 PF00104 PF00104          16          -0.731 -2.17 1.04e- 4    3.02e-3 2.41e-3
    ##  5 PF00105 PF00105          16          -0.731 -2.17 1.04e- 4    3.02e-3 2.41e-3
    ##  6 PF00227 PF00227          16           0.740  2.16 9.04e- 5    3.02e-3 2.41e-3
    ##  7 PF00153 PF00153          34          -0.543 -1.96 2.53e- 4    6.29e-3 5.02e-3
    ##  8 PF00118 PF00118          11           0.765  2.02 3.53e- 4    7.68e-3 6.14e-3
    ##  9 PF00651 PF00651          25          -0.594 -1.96 5.63e- 4    1.09e-2 8.69e-3
    ## 10 PF00856 PF00856          23          -0.607 -1.98 6.43e- 4    1.11e-2 8.90e-3
    ## # ℹ 164 more rows
    ## # ℹ 3 more variables: rank <dbl>, leading_edge <chr>, core_enrichment <chr>

``` r
sessionInfo()
```

    ## R version 4.4.0 (2024-04-24)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Sonoma 14.6.1
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/New_York
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] org.Hs.eg.db_3.20.0    AnnotationDbi_1.68.0   IRanges_2.40.1        
    ##  [4] S4Vectors_0.44.0       Biobase_2.66.0         BiocGenerics_0.52.0   
    ##  [7] clusterProfiler_4.14.4 lubridate_1.9.4        forcats_1.0.0         
    ## [10] stringr_1.5.1          dplyr_1.1.4            purrr_1.0.4           
    ## [13] readr_2.1.5            tidyr_1.3.1            tibble_3.2.1          
    ## [16] ggplot2_3.5.1          tidyverse_2.0.0       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] DBI_1.2.3               gson_0.1.0              rlang_1.1.5            
    ##  [4] magrittr_2.0.3          DOSE_4.0.0              compiler_4.4.0         
    ##  [7] RSQLite_2.3.9           png_0.1-8               vctrs_0.6.5            
    ## [10] reshape2_1.4.4          pkgconfig_2.0.3         crayon_1.5.3           
    ## [13] fastmap_1.2.0           XVector_0.46.0          utf8_1.2.4             
    ## [16] rmarkdown_2.29          tzdb_0.4.0              enrichplot_1.26.6      
    ## [19] UCSC.utils_1.2.0        bit_4.5.0.1             xfun_0.50              
    ## [22] zlibbioc_1.52.0         cachem_1.1.0            aplot_0.2.4            
    ## [25] GenomeInfoDb_1.42.3     jsonlite_1.8.9          blob_1.2.4             
    ## [28] BiocParallel_1.40.0     parallel_4.4.0          R6_2.6.1               
    ## [31] stringi_1.8.4           RColorBrewer_1.1-3      GOSemSim_2.32.0        
    ## [34] Rcpp_1.0.14             knitr_1.49              ggtangle_0.0.6         
    ## [37] R.utils_2.12.3          Matrix_1.7-2            splines_4.4.0          
    ## [40] igraph_2.1.4            timechange_0.3.0        tidyselect_1.2.1       
    ## [43] qvalue_2.38.0           rstudioapi_0.17.1       yaml_2.3.10            
    ## [46] codetools_0.2-20        lattice_0.22-6          plyr_1.8.9             
    ## [49] treeio_1.30.0           withr_3.0.2             KEGGREST_1.46.0        
    ## [52] evaluate_1.0.3          gridGraphics_0.5-1      Biostrings_2.74.1      
    ## [55] pillar_1.10.1           ggtree_3.14.0           ggfun_0.1.8            
    ## [58] generics_0.1.3          vroom_1.6.5             hms_1.1.3              
    ## [61] munsell_0.5.1           scales_1.3.0            tidytree_0.4.6         
    ## [64] glue_1.8.0              lazyeval_0.2.2          tools_4.4.0            
    ## [67] data.table_1.16.4       fgsea_1.32.2            fs_1.6.5               
    ## [70] fastmatch_1.1-6         cowplot_1.1.3           grid_4.4.0             
    ## [73] ape_5.8-1               colorspace_2.1-1        nlme_3.1-167           
    ## [76] GenomeInfoDbData_1.2.13 patchwork_1.3.0         cli_3.6.4              
    ## [79] gtable_0.3.6            R.methodsS3_1.8.2       yulab.utils_0.2.0      
    ## [82] digest_0.6.37           ggrepel_0.9.6           ggplotify_0.1.2        
    ## [85] farver_2.1.2            memoise_2.0.1           htmltools_0.5.8.1      
    ## [88] R.oo_1.27.0             lifecycle_1.0.4         httr_1.4.7             
    ## [91] GO.db_3.20.0            bit64_4.6.0-1
