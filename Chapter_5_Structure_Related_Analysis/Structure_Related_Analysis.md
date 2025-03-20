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
  head(10)
```

    ## # A tibble: 10 × 18
    ##    Index       UniProt_Accession Protein.Start Gene   Entry.Name    Kd half_life
    ##    <chr>       <chr>                     <dbl> <chr>  <chr>      <dbl>     <dbl>
    ##  1 P49368_121  P49368                      121 CCT3   TCPG_HUMAN  2.36     0.298
    ##  2 Q14157_778  Q14157                      778 UBAP2L UBP2L_HUM…  2.15     0.328
    ##  3 Q15154_1067 Q15154                     1067 PCM1   PCM1_HUMAN  2.09     0.339
    ##  4 Q15393_740  Q15393                      740 SF3B3  SF3B3_HUM…  1.96     0.360
    ##  5 O75436_252  O75436                      252 VPS26A VP26A_HUM…  1.94     0.365
    ##  6 Q08211_820  Q08211                      820 DHX9   DHX9_HUMAN  1.94     0.366
    ##  7 Q9UN37_307  Q9UN37                      307 VPS4A  VPS4A_HUM…  1.86     0.381
    ##  8 P50402_185  P50402                      185 EMD    EMD_HUMAN   1.75     0.406
    ##  9 Q99615_343  Q99615                      343 DNAJC7 DNJC7_HUM…  1.74     0.408
    ## 10 Q13813_429  Q13813                      429 SPTAN1 SPTN1_HUM…  1.71     0.415
    ## # ℹ 11 more variables: RSS <dbl>, Percentile <dbl>, category <chr>,
    ## #   Sequence <chr>, Full_Protein_Length <dbl>, Nterm_sequence <chr>,
    ## #   Nterm_13mer <chr>, Nterm_terminus <chr>,
    ## #   Nterm_13mers_sequence_parameter <list>, hydropathy <dbl>,
    ## #   isoelectric_point <dbl>

``` r
# protein modification structure analysis
common_Nterm_alphafold_N_terminus |>
  as_tibble() |>
  head(10)
```

    ## # A tibble: 10 × 33
    ##    protein_id protein_number AA    position quality x_coord_c x_coord_ca
    ##    <chr>               <dbl> <chr>    <dbl>   <dbl>     <dbl>      <dbl>
    ##  1 E9PAV3                  1 M            1    43.7      32.1       31.5
    ##  2 E9PAV3                  1 P            2    36.9      32.6       31.8
    ##  3 E9PAV3                  1 G            3    37.2      31.5       32.7
    ##  4 E9PAV3                  1 E            4    38.0      32.2       31.8
    ##  5 E9PAV3                  1 A            5    42.9      31.0       31.4
    ##  6 E9PAV3                  1 T            6    41.5      28.7       30.2
    ##  7 E9PAV3                  1 E            7    43.5      27.9       28.8
    ##  8 E9PAV3                  1 T            8    44.0      27.5       27.5
    ##  9 E9PAV3                  1 V            9    39.2      26.7       26.3
    ## 10 E9PAV3                  1 P           10    41.3      27.2       28.1
    ## # ℹ 26 more variables: x_coord_cb <dbl>, x_coord_n <dbl>, y_coord_c <dbl>,
    ## #   y_coord_ca <dbl>, y_coord_cb <dbl>, y_coord_n <dbl>, z_coord_c <dbl>,
    ## #   z_coord_ca <dbl>, z_coord_cb <dbl>, z_coord_n <dbl>,
    ## #   secondary_structure <chr>, structure_group <chr>, BEND <dbl>, HELX <dbl>,
    ## #   STRN <dbl>, TURN <dbl>, unstructured <dbl>, nAA_24_180_pae <dbl>,
    ## #   nAA_12_70_pae <dbl>, high_acc_5 <dbl>, low_acc_5 <dbl>,
    ## #   nAA_24_180_pae_smooth10 <dbl>, IDR <dbl>, common_Nterm <dbl>, …

``` r
enrichment_N_terminus |>
  as_tibble() |>
  head(10)
```

    ## # A tibble: 7 × 13
    ##   quality_cutoff ptm       roi   n_aa_ptm n_aa_roi n_ptm_in_roi n_ptm_not_in_roi
    ##            <dbl> <chr>     <chr>    <dbl>    <dbl>        <dbl>            <dbl>
    ## 1              0 common_N… BEND      2869    23359          140             2729
    ## 2              0 common_N… HELX      2869   144289         1070             1799
    ## 3              0 common_N… STRN      2869    51865          436             2433
    ## 4              0 common_N… TURN      2869    32157          217             2652
    ## 5              0 common_N… IDR       2869   151391         1188             1681
    ## 6              0 common_N… high…     2869   298954         2238              631
    ## 7              0 common_N… low_…     2869    87472          631             2238
    ## # ℹ 6 more variables: n_naked_in_roi <dbl>, n_naked_not_in_roi <dbl>,
    ## #   oddsr <dbl>, p <dbl>, p_adj_bf <dbl>, p_adj_bh <dbl>

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
