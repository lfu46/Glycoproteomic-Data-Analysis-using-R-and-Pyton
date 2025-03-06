Data Normalization
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
# load the peptides result from searching software
peptide_result_raw <- read_csv(
  # the csv file is too large to be uploaded to github
  # change the path to your result files
  'training_data/ronghuwulab_1741203378.csv',
  col_names = TRUE,
  name_repair = 'universal'
)
```

    ## New names:
    ## Rows: 119682 Columns: 105
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "," chr
    ## (14): Obs.m.z.Link, SrchName, ..Ions, ..Ions.Link, Reference, Reference.... dbl
    ## (91): ScanF, Obs.m.z, z, SrchID, RunID, ScansID, PPM, Dalton, XCorr, ..9...
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `Obs m/z` -> `Obs.m.z`
    ## • `Obs m/z Link` -> `Obs.m.z.Link`
    ## • `&#916;Corr` -> `..916.Corr`
    ## • `# Ions` -> `..Ions`
    ## • `# Ions Link` -> `..Ions.Link`
    ## • `Reference Link` -> `Reference.Link`
    ## • `Redun Link` -> `Redun.Link`
    ## • `Peptide Link` -> `Peptide.Link`
    ## • `Trimmed Peptide` -> `Trimmed.Peptide`
    ## • `Pept. Length` -> `Pept..Length`
    ## • `Start Position` -> `Start.Position`
    ## • `End Position` -> `End.Position`
    ## • `Gene Symbol` -> `Gene.Symbol`
    ## • `Protein MWT(kDa)` -> `Protein.MWT.kDa.`
    ## • `127nto126` -> `..127nto126`
    ## • `127nto126 Link` -> `..127nto126.Link`
    ## • `131cto126` -> `..131cto126`
    ## • `127cto126` -> `..127cto126`
    ## • `128nto126` -> `..128nto126`
    ## • `128cto126` -> `..128cto126`
    ## • `129nto126` -> `..129nto126`
    ## • `129cto126` -> `..129cto126`
    ## • `130nto126` -> `..130nto126`
    ## • `130cto126` -> `..130cto126`
    ## • `131to126` -> `..131to126`
    ## • `126 Adjusted Intensity` -> `..126.Adjusted.Intensity`
    ## • `127n Adjusted Intensity` -> `..127n.Adjusted.Intensity`
    ## • `127c Adjusted Intensity` -> `..127c.Adjusted.Intensity`
    ## • `128n Adjusted Intensity` -> `..128n.Adjusted.Intensity`
    ## • `128c Adjusted Intensity` -> `..128c.Adjusted.Intensity`
    ## • `129n Adjusted Intensity` -> `..129n.Adjusted.Intensity`
    ## • `129c Adjusted Intensity` -> `..129c.Adjusted.Intensity`
    ## • `130n Adjusted Intensity` -> `..130n.Adjusted.Intensity`
    ## • `130c Adjusted Intensity` -> `..130c.Adjusted.Intensity`
    ## • `131 Adjusted Intensity` -> `..131.Adjusted.Intensity`
    ## • `126 Sn` -> `..126.Sn`
    ## • `127n Sn` -> `..127n.Sn`
    ## • `127c Sn` -> `..127c.Sn`
    ## • `128n Sn` -> `..128n.Sn`
    ## • `128c Sn` -> `..128c.Sn`
    ## • `129n Sn` -> `..129n.Sn`
    ## • `129c Sn` -> `..129c.Sn`
    ## • `130n Sn` -> `..130n.Sn`
    ## • `130c Sn` -> `..130c.Sn`
    ## • `131 Sn` -> `..131.Sn`
    ## • `126 Raw Intensity` -> `..126.Raw.Intensity`
    ## • `127n Raw Intensity` -> `..127n.Raw.Intensity`
    ## • `127c Raw Intensity` -> `..127c.Raw.Intensity`
    ## • `128n Raw Intensity` -> `..128n.Raw.Intensity`
    ## • `128c Raw Intensity` -> `..128c.Raw.Intensity`
    ## • `129n Raw Intensity` -> `..129n.Raw.Intensity`
    ## • `129c Raw Intensity` -> `..129c.Raw.Intensity`
    ## • `130n Raw Intensity` -> `..130n.Raw.Intensity`
    ## • `130c Raw Intensity` -> `..130c.Raw.Intensity`
    ## • `131 Raw Intensity` -> `..131.Raw.Intensity`
    ## • `126 Noise` -> `..126.Noise`
    ## • `127n Noise` -> `..127n.Noise`
    ## • `127c Noise` -> `..127c.Noise`
    ## • `128n Noise` -> `..128n.Noise`
    ## • `128c Noise` -> `..128c.Noise`
    ## • `129n Noise` -> `..129n.Noise`
    ## • `129c Noise` -> `..129c.Noise`
    ## • `130n Noise` -> `..130n.Noise`
    ## • `130c Noise` -> `..130c.Noise`
    ## • `131 Noise` -> `..131.Noise`
    ## • `126 Mz Error` -> `..126.Mz.Error`
    ## • `127n Mz Error` -> `..127n.Mz.Error`
    ## • `127c Mz Error` -> `..127c.Mz.Error`
    ## • `128n Mz Error` -> `..128n.Mz.Error`
    ## • `128c Mz Error` -> `..128c.Mz.Error`
    ## • `129n Mz Error` -> `..129n.Mz.Error`
    ## • `129c Mz Error` -> `..129c.Mz.Error`
    ## • `130n Mz Error` -> `..130n.Mz.Error`
    ## • `130c Mz Error` -> `..130c.Mz.Error`
    ## • `131 Mz Error` -> `..131.Mz.Error`
    ## • `Sum Sn` -> `Sum.Sn`
    ## • `Ion Injection Time` -> `Ion.Injection.Time`
    ## • `Ms2 Ion Inj Time` -> `Ms2.Ion.Inj.Time`
    ## • `Elapsed Scan Time` -> `Elapsed.Scan.Time`
    ## • `Precursor Intensity` -> `Precursor.Intensity`
    ## • `Precursor Max Intensity` -> `Precursor.Max.Intensity`
    ## • `Precursor Max Sn` -> `Precursor.Max.Sn`
    ## • `Isolation Mz` -> `Isolation.Mz`
    ## • `Isolation Specificity` -> `Isolation.Specificity`
    ## • `Time To Max` -> `Time.To.Max`
    ## • `Difference Max Intensity` -> `Difference.Max.Intensity`
    ## • `Peak Width` -> `Peak.Width`
    ## • `Num Precursor Peaks` -> `Num.Precursor.Peaks`

``` r
peptide_result_raw_no_decoy_contaminants <- peptide_result_raw |> 
  # only keep the psm from human proteome
  filter(str_detect(Reference, 'HUMAN')) |> 
  # remove contaminant
  filter(!str_detect(Reference, 'contaminant')) |> 
  # remove decoy
  filter(!str_detect(Reference, '##')) |> 
  # remove low quality PSMs
  # for Sequest, XCorr and PPM are used
  # for other searching software, please use corresponding searching score
  # for example, Hyperscore in MSFragger
  filter(XCorr > 1.2) |> 
  filter(PPM > -10, PPM < 10) |> 
  filter(Sum.Sn > 45) |> 
  separate(Reference, into = c('sp', 'UniProt_Accession', 'protein_name'), sep = '\\|') |> 
  select(
    # keep the columns which are useful in the following data analysis
    Peptide,
    Trimmed.Peptide, 
    Pept.Length = Pept..Length, 
    Start.Position, 
    End.Position, 
    Trypticity, 
    MissedCleav,
    Obs.m.z, 
    XCorr, 
    PPM, 
    UniProt_Accession, 
    Gene.Symbol, 
    Annotation, 
    # rename the TMT channel base on your experimental design
    Iron_1 = ..126.Sn, Iron_2 = ..127n.Sn, Iron_3 = ..127c.Sn,
    Copper_4 = ..128n.Sn, Copper_5 = ..128c.Sn, Copper_6 = ..129n.Sn,
    Ctrl_7 = ..129c.Sn, Ctrl_8 = ..130n.Sn, Ctrl_9 = ..130c.Sn, 
    Sum.Sn
  )
```

``` r
# quantification on protein level
protein_quantification <- peptide_result_raw_no_decoy_contaminants |> 
  group_by(
    UniProt_Accession, Gene.Symbol, Annotation
  ) |> 
  summarise(
    Iron_1 = sum(Iron_1),
    Iron_2 = sum(Iron_2),
    Iron_3 = sum(Iron_3),
    Copper_4 = sum(Copper_4),
    Copper_5 = sum(Copper_5),
    Copper_6 = sum(Copper_6),
    Ctrl_7 = sum(Ctrl_7),
    Ctrl_8 = sum(Ctrl_8),
    Ctrl_9 = sum(Ctrl_9)
  ) |> 
  ungroup()
```

    ## `summarise()` has grouped output by 'UniProt_Accession', 'Gene.Symbol'. You can
    ## override using the `.groups` argument.

``` r
# check the sum of each TMT channel
colSums(
  protein_quantification |> select(Iron_1:Ctrl_9)
)
```

    ##   Iron_1   Iron_2   Iron_3 Copper_4 Copper_5 Copper_6   Ctrl_7   Ctrl_8 
    ## 49494037 34450910 42939901 40187880 32316432 37936857 49340660 27321462 
    ##   Ctrl_9 
    ## 24539529

``` r
# sample loading normalization
target_mean_protein <- mean(colSums(protein_quantification |> select(Iron_1:Ctrl_9)))
norm_facs_protein <- target_mean_protein/colSums(protein_quantification |> select(Iron_1:Ctrl_9))
protein_sl <- tibble(
  sweep(
    protein_quantification |> select(Iron_1:Ctrl_9), 
    2, 
    norm_facs_protein, 
    FUN = '*'
  )
)
colnames(protein_sl) <- c(
  'Iron_1_sl', 'Iron_2_sl', 'Iron_3_sl', 
  'Copper_4_sl', 'Copper_5_sl', 'Copper_6_sl',
  'Ctrl_7_sl', 'Ctrl_8_sl', 'Ctrl_9_sl'
)
protein_quantification_raw_sl <- bind_cols(protein_quantification, protein_sl)
```

``` r
# check the sum of each TMT channel
colSums(
  protein_quantification_raw_sl |> select(Iron_1_sl:Ctrl_9_sl)
)
```

    ##   Iron_1_sl   Iron_2_sl   Iron_3_sl Copper_4_sl Copper_5_sl Copper_6_sl 
    ##    37614186    37614186    37614186    37614186    37614186    37614186 
    ##   Ctrl_7_sl   Ctrl_8_sl   Ctrl_9_sl 
    ##    37614186    37614186    37614186

``` r
# TMM normalization
library(edgeR)
```

    ## Warning: package 'edgeR' was built under R version 4.4.2

    ## Loading required package: limma

    ## Warning: package 'limma' was built under R version 4.4.2

``` r
norm_facs_protein_sl_tmm <- calcNormFactors(protein_quantification_raw_sl |> select(Iron_1_sl:Ctrl_9_sl))
protein_tmm <- tibble(
  sweep(
    protein_quantification_raw_sl |> select(Iron_1_sl:Ctrl_9_sl), 
    2, 
    norm_facs_protein_sl_tmm, 
    FUN = "/"
  )
)
colnames(protein_tmm) <- c(
  'Iron_1_sl_tmm', 'Iron_2_sl_tmm', 'Iron_3_sl_tmm', 
  'Copper_4_sl_tmm', 'Copper_5_sl_tmm', 'Copper_6_sl_tmm',
  'Ctrl_7_sl_tmm', 'Ctrl_8_sl_tmm', 'Ctrl_9_sl_tmm'
)
protein_quantification_raw_sl_tmm <- bind_cols(protein_quantification_raw_sl, protein_tmm)
```

``` r
# check the distribution of the intensity of each channel
protein_quantification_raw_sl_tmm |> 
  select(Iron_1_sl_tmm:Ctrl_9_sl_tmm) |> 
  pivot_longer(cols = Iron_1_sl_tmm:Ctrl_9_sl_tmm, names_to = 'Exp', values_to = 'Intensity') |> 
  mutate(
    log2_intensity = log2(Intensity)
  ) |> 
  ggplot() +
  geom_boxplot(
    aes(
      x = Exp,
      y = log2_intensity
    )
  ) +
  labs(x = '', y = 'log2 intensity') +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
  )
```

![](Data_Normalization_files/figure-gfm/check%20intensity%20distribution-1.png)<!-- -->
