Differential Analysis
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
