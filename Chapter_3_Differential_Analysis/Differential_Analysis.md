Chapter 3: Differential Analysis
================
Longping Fu
03/05/2025

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
library(limma)
```

    ## Warning: package 'limma' was built under R version 4.4.2

``` r
# experimental model
Experiment_Model <- model.matrix(~ 0 + factor(rep(c("case1", "case2"), each = 3), levels = c("case1", "case2")))
colnames(Experiment_Model) <- c("case1", "case2")
Contrast_Matrix <- makeContrasts(case1_case2 = case1 - case2, levels = Experiment_Model)
```

``` r
# import normalized data
protein_quantification_raw_sl_tmm <- read_csv(
  'training_data/protein_quantification_raw_sl_tmm.csv'
)
```

    ## Rows: 7107 Columns: 30
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr  (3): UniProt_Accession, Gene.Symbol, Annotation
    ## dbl (27): Iron_1, Iron_2, Iron_3, Copper_4, Copper_5, Copper_6, Ctrl_7, Ctrl...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# calculate log2 transformed Iron experiment TMT channel intensity
Iron_protein_raw_sl_tmm_log2 <- protein_quantification_raw_sl_tmm |> 
  mutate(
    log2_Iron_1_sl_tmm = log2(Iron_1_sl_tmm),
    log2_Iron_2_sl_tmm = log2(Iron_2_sl_tmm),
    log2_Iron_3_sl_tmm = log2(Iron_3_sl_tmm),
    log2_Ctrl_7_sl_tmm = log2(Ctrl_7_sl_tmm),
    log2_Ctrl_8_sl_tmm = log2(Ctrl_8_sl_tmm),
    log2_Ctrl_9_sl_tmm = log2(Ctrl_9_sl_tmm)
  )

# differential analysis
Iron_protein_raw_sl_tmm_log2_data_matrix <- data.matrix(
  Iron_protein_raw_sl_tmm_log2 |> select(starts_with('log2'))
)
rownames(Iron_protein_raw_sl_tmm_log2_data_matrix) <- Iron_protein_raw_sl_tmm_log2$UniProt_Accession

Iron_protein_lmfit <- lmFit(Iron_protein_raw_sl_tmm_log2_data_matrix, Experiment_Model)
Iron_protein_lmFit_contrast <- contrasts.fit(Iron_protein_lmfit, Contrast_Matrix)
Iron_protein_lmFit_contrast_eBayes <- eBayes(Iron_protein_lmFit_contrast)
Iron_protein_toptable <- topTable(
  Iron_protein_lmFit_contrast_eBayes,
  number = Inf,
  adjust.method = 'BH'
)

rownames_Iron_protein_toptable <- rownames(Iron_protein_toptable)
Iron_protein_toptable_tb <- tibble(Iron_protein_toptable)
Iron_protein_toptable_tb$UniProt_Accession <- rownames_Iron_protein_toptable

write_csv(
  Iron_protein_toptable_tb,
  file = 'training_data/Iron_protein_toptable_tb.csv'
)
```

``` r
Iron_protein_toptable_tb |> 
  mutate(
    log10_adjpvalue = -log10(adj.P.Val)
  ) |> 
  ggplot() +
  geom_point(
    aes(
      x = logFC,
      y = log10_adjpvalue
    )
  ) +
  labs(x = 'log2 fold change', y = '-log10 adj.P.value')
```

![](Differential_Analysis_files/figure-gfm/volcano%20plot-1.png)<!-- -->

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
    ##  [1] limma_3.62.2    lubridate_1.9.4 forcats_1.0.0   stringr_1.5.1  
    ##  [5] dplyr_1.1.4     purrr_1.0.4     readr_2.1.5     tidyr_1.3.1    
    ##  [9] tibble_3.2.1    ggplot2_3.5.1   tidyverse_2.0.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] bit_4.5.0.1       gtable_0.3.6      crayon_1.5.3      compiler_4.4.0   
    ##  [5] tidyselect_1.2.1  parallel_4.4.0    scales_1.3.0      statmod_1.5.0    
    ##  [9] yaml_2.3.10       fastmap_1.2.0     R6_2.6.1          labeling_0.4.3   
    ## [13] generics_0.1.3    knitr_1.49        munsell_0.5.1     pillar_1.10.1    
    ## [17] tzdb_0.4.0        rlang_1.1.5       stringi_1.8.4     xfun_0.50        
    ## [21] bit64_4.6.0-1     timechange_0.3.0  cli_3.6.4         withr_3.0.2      
    ## [25] magrittr_2.0.3    digest_0.6.37     grid_4.4.0        vroom_1.6.5      
    ## [29] rstudioapi_0.17.1 hms_1.1.3         lifecycle_1.0.4   vctrs_0.6.5      
    ## [33] evaluate_1.0.3    glue_1.8.0        farver_2.1.2      colorspace_2.1-1 
    ## [37] rmarkdown_2.29    tools_4.4.0       pkgconfig_2.0.3   htmltools_0.5.8.1
