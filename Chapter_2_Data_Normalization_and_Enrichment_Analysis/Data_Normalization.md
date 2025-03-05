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
peptide_result <- read_csv(
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
