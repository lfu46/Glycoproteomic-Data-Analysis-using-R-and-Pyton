Chapter 5: Structure Related Analysis
================
Longping Fu
03/20/2025

To perform structure-related analysis of peptide sequences and
post-translational modifications, you first need to install two Python
packages: [localCIDER](https://pappulab.github.io/localCIDER/) and
[StructureMap](https://github.com/MannLabs/structuremap). Please refer
to their documentation for installation instructions and usage details.

Additionally, you can use
[Anaconda](https://docs.anaconda.com/anaconda/install/) to set up a
virtual environment for managing dependencies efficiently.

``` r
# import packages
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
library(reticulate)
```

    ## Warning: package 'reticulate' was built under R version 4.4.1

``` r
# use specific virtual environment you just created
use_condaenv(
  condaenv = '/opt/anaconda3/envs/structure_analysis',
  required = TRUE
)
```

``` r
# execute the python script for structure analysis
source_python("structure_analysis.py")
```

    ## Warning in py_to_r.pandas.core.frame.DataFrame(<environment>): index contains
    ## duplicated values: row names not set

``` r
## check the result
# peptide physicochemical property analysis
HEK_Nterm_Kd_half_life_sequence |> 
  as_tibble() |> 
  select(Nterm_13mer, hydropathy, isoelectric_point) |> 
  head(10)
```

    ## # A tibble: 10 × 3
    ##    Nterm_13mer   hydropathy isoelectric_point
    ##    <chr>              <dbl>             <dbl>
    ##  1 VVISAYRKALDDM       4.86              7   
    ##  2 SSVATTSGKAPPN       3.96             12.2 
    ##  3 WQQNNVQRLKQML       3.22             13.9 
    ##  4 EFASGFASEQCPE       4.05              1.75
    ##  5 LAGYDPTPTMRDV       4.03              3.94
    ##  6 AIEPPPLDAVIEA       5.05              1.75
    ##  7 LGSTPHNLTDANI       4.19              5.25
    ##  8 SSSTSFMSSSSSS       4.19              7   
    ##  9 QYEEAVRDYEKVY       2.99              4.05
    ## 10 AAGHYASDEVREK       3.35              5.25

``` r
# protein modification structure analysis
common_Nterm_alphafold_N_terminus |>
  as_tibble() |>
  select(protein_id, AA, position, structure_group, nAA_24_180_pae, nAA_12_70_pae, high_acc_5, low_acc_5, IDR) |> 
  head(10)
```

    ## # A tibble: 10 × 9
    ##    protein_id AA    position structure_group nAA_24_180_pae nAA_12_70_pae
    ##    <chr>      <chr>    <dbl> <chr>                    <dbl>         <dbl>
    ##  1 E9PAV3     M            1 unstructured                 5             0
    ##  2 E9PAV3     P            2 unstructured                 6             0
    ##  3 E9PAV3     G            3 unstructured                 7             0
    ##  4 E9PAV3     E            4 unstructured                 8             0
    ##  5 E9PAV3     A            5 unstructured                 9             0
    ##  6 E9PAV3     T            6 unstructured                10             0
    ##  7 E9PAV3     E            7 unstructured                10             0
    ##  8 E9PAV3     T            8 unstructured                10             0
    ##  9 E9PAV3     V            9 unstructured                10             0
    ## 10 E9PAV3     P           10 unstructured                10             0
    ## # ℹ 3 more variables: high_acc_5 <dbl>, low_acc_5 <dbl>, IDR <dbl>

``` r
enrichment_N_terminus |>
  as_tibble() |>
  select(ptm, roi, oddsr, p, p_adj_bf, p_adj_bh) |> 
  head(10)
```

    ## # A tibble: 7 × 6
    ##   ptm          roi        oddsr       p p_adj_bf p_adj_bh
    ##   <chr>        <chr>      <dbl>   <dbl>    <dbl>    <dbl>
    ## 1 common_Nterm BEND       0.796 0.00749   0.0524   0.0262
    ## 2 common_Nterm HELX       0.998 0.969     1        0.969 
    ## 3 common_Nterm STRN       1.16  0.00596   0.0417   0.0262
    ## 4 common_Nterm TURN       0.901 0.145     1        0.253 
    ## 5 common_Nterm IDR        1.10  0.0148    0.103    0.0345
    ## 6 common_Nterm high_acc_5 1.04  0.420     1        0.490 
    ## 7 common_Nterm low_acc_5  0.963 0.420     1        0.490

``` r
common_Nterm_proximity |>
  as_tibble() |>
  head(10)
```

    ## # A tibble: 10 × 7
    ##    protein_id ptm   n_ptms pvalue_1d pvalue_3d pvalue_1d_adj_bh pvalue_3d_adj_bh
    ##    <chr>      <chr>  <dbl>     <dbl>     <dbl>            <dbl>            <dbl>
    ##  1 E9PAV3     comm…      9    0         0                0                0     
    ##  2 O00151     comm…      5    0.441     0.758            0.693            0.924 
    ##  3 O00193     comm…      7    0.0117    0.0093           0.0819           0.0747
    ##  4 O00273     comm…      5    0.588     0.411            0.799            0.712 
    ##  5 O00299     comm…      2    0.731     0.985            0.869            0.998 
    ##  6 O14545     comm…      2    0.0072    0.0078           0.0652           0.0689
    ##  7 O14556     comm…      3    0.108     0.0786           0.318            0.261 
    ##  8 O14561     comm…      4    0.0265    0.175            0.129            0.432 
    ##  9 O14737     comm…      4    0.386     0.327            0.643            0.651 
    ## 10 O14745     comm…      2    0.129     0.166            0.353            0.416

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
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] reticulate_1.40.0 lubridate_1.9.4   forcats_1.0.0     stringr_1.5.1    
    ##  [5] dplyr_1.1.4       purrr_1.0.4       readr_2.1.5       tidyr_1.3.1      
    ##  [9] tibble_3.2.1      ggplot2_3.5.1     tidyverse_2.0.0  
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Matrix_1.7-2      jsonlite_1.8.9    gtable_0.3.6      compiler_4.4.0   
    ##  [5] Rcpp_1.0.14       tidyselect_1.2.1  png_0.1-8         scales_1.3.0     
    ##  [9] yaml_2.3.10       fastmap_1.2.0     lattice_0.22-6    R6_2.6.1         
    ## [13] generics_0.1.3    knitr_1.49        munsell_0.5.1     pillar_1.10.1    
    ## [17] tzdb_0.4.0        rlang_1.1.5       utf8_1.2.4        stringi_1.8.4    
    ## [21] xfun_0.50         timechange_0.3.0  cli_3.6.4         withr_3.0.2      
    ## [25] magrittr_2.0.3    digest_0.6.37     grid_4.4.0        rstudioapi_0.17.1
    ## [29] hms_1.1.3         lifecycle_1.0.4   vctrs_0.6.5       evaluate_1.0.3   
    ## [33] glue_1.8.0        colorspace_2.1-1  rmarkdown_2.29    tools_4.4.0      
    ## [37] pkgconfig_2.0.3   htmltools_0.5.8.1
