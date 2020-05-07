
<!-- README.md is generated from README.Rmd. Please edit that file -->

# stimcellselector

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/rohitfarmer/stimcellselector.svg?branch=master)](https://travis-ci.com/rohitfarmer/stimcellselector)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/rohitfarmer/stimcellselector?branch=master&svg=true)](https://ci.appveyor.com/project/rohitfarmer/stimcellselector)
[![Codecov test
coverage](https://codecov.io/gh/rohitfarmer/stimcellselector/branch/master/graph/badge.svg)](https://codecov.io/gh/rohitfarmer/stimcellselector?branch=master)
<!-- badges: end -->

The goal of this package is to select cells from stimulated samples that
have responded to the stimulant in CyTOF/Flow cytometry stimulation
assays. Starting from the identified cell populations either through
automated clustering such as FlowSOM or traditional cell gating, the
primary function `stim_cell_selector()` follows a heuristic approach to
group cells into responding and non-responding.

For a combination of cell population and stimulation type (e.g., CD127+
T-helper cells and interferon-alpha), `stim_cell_selector()` starts by
performing k-means clustering on the combined set of cells from
stimulated and unstimulated samples. K-means clustering is performed on
expression data of all the state markers combined. Upon clustering using
a contingency table as drawn below, a Fisher’s exact test determines the
effect size and the statistical significance of partitioning. Cells form
the combinations that pass the Fisher’s exact test are considered as
responding. An optional UMAP plot can also be generated to verify the
cell partitioning in responding and non-responding groups
visually.

``` r
matrix(nrow = 2, ncol = 2, dimnames = list(c("Cluster1", "Cluster2"), c("Stim", "Unstim")))
#>          Stim Unstim
#> Cluster1   NA     NA
#> Cluster2   NA     NA
```

An optional UMAP plot can also be generated to verify the cell
partitioning in responding and non-responding groups visually.

## Installation

You can install the released version of stimcellselector from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("stimcellselector")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("rohitfarmer/stimcellselector")
```

## Example

An example using the sample CyTOF data set `chi11_1k` included in the
package. `suppressMessages()` is used to reduce the verbosity of the
`stim_cell_selector` function.

``` r
library(stimcellselector)

selected_data <- suppressMessages(stim_cell_selector(chi11_1k$expr_data, chi11_1k$state_markers,
                  chi11_1k$cluster_col, chi11_1k$stim_label,
                  chi11_1k$unstim_label, seed_val = 123, umap = TRUE, umap_cells = 50))

class(selected_data)
#> [1] "list"

attributes(selected_data)
#> $names
#> [1] "selected_expr_data"    "summary"               "stacked_bar_plot_data"
#> [4] "umap_plot_data"
```

### Output

`stim_cell_selector()` returns a list with the selected expression data,
a summary, data to plot stacked bar plots to visualize the k-means and
Fisher’s exact test results and data to plot the optional UMAPs.

``` r
head(selected_data$selected_expr_data)
#> # A tibble: 6 x 36
#>   cluster_id sample_id condition patient_id stim_type  CD45   CD7   CD19
#>   <fct>      <fct>     <fct>     <fct>      <fct>     <dbl> <dbl>  <dbl>
#> 1 42         CHI-011_… CHI       CHI-011    A          2.77  3.71 1.20  
#> 2 24         CHI-011_… CHI       CHI-011    A          2.59  4.42 0.435 
#> 3 31         CHI-011_… CHI       CHI-011    A          3.17  3.89 0     
#> 4 31         CHI-011_… CHI       CHI-011    A          3.58  2.78 0     
#> 5 32         CHI-011_… CHI       CHI-011    A          2.93  3.59 0     
#> 6 11         CHI-011_… CHI       CHI-011    A          3.29  3.76 0.0348
#> # … with 28 more variables: pPLCg2 <dbl>, CD4 <dbl>, IgD <dbl>,
#> #   CD20 <dbl>, CD25 <dbl>, pSTAT5 <dbl>, CD123 <dbl>, AKT <dbl>,
#> #   pSTAT1 <dbl>, CD27 <dbl>, pP38 <dbl>, CD24 <dbl>, pSTAT3 <dbl>,
#> #   CD11c <dbl>, CD14 <dbl>, CD56 <dbl>, IkBa <dbl>, pCREB <dbl>,
#> #   CD16 <dbl>, CD38 <dbl>, CD8 <dbl>, CD45RA <dbl>, CD3 <dbl>,
#> #   pERK1_2 <dbl>, HLA_DR <dbl>, pS6 <dbl>, CD127 <dbl>, merging1 <fct>
```

``` r
head(selected_data$summary)
#> # A tibble: 6 x 14
#>   cluster stim_type f_p_value fold_change pPLCg2 pSTAT5   AKT pSTAT1
#>   <fct>   <fct>         <dbl>       <dbl>  <dbl>  <dbl> <dbl>  <dbl>
#> 1 CD3 CD… A          2.16e-20       0.64   0.372  1.43  0      1.79 
#> 2 CD3 CD… A          1.36e-22       0.67   0.487  1.83  0      2.16 
#> 3 CD11c … A          2.22e-28       0.72   0.408  1.46  0      2.66 
#> 4 CD19 C… A          2.50e-26       0.72   0.239  1.83  0      2.12 
#> 5 CD3 CD… A          9.20e-14       0.570  0.597  1.16  0.118  1.56 
#> 6 CD3 CD… A          3.08e- 2      -0.99   0.113  0.210 0      0.425
#> # … with 6 more variables: pP38 <dbl>, pSTAT3 <dbl>, IkBa <dbl>,
#> #   pCREB <dbl>, pERK1_2 <dbl>, pS6 <dbl>
```

``` r
head(selected_data$stacked_bar_plot_data)
#> # A tibble: 6 x 5
#>   cluster              stim_type stim_status k_cluster count
#>   <fct>                <fct>     <fct>       <fct>     <dbl>
#> 1 CD3 CD8 CD127 CD45RA A         unstim      cluster1   1474
#> 2 CD3 CD8 CD127 CD45RA A         unstim      cluster2    178
#> 3 CD3 CD8 CD127 CD45RA A         stim        cluster1    423
#> 4 CD3 CD8 CD127 CD45RA A         stim        cluster2   1182
#> 5 CD3 CD8 CD27 CD127   A         unstim      cluster1   3370
#> 6 CD3 CD8 CD27 CD127   A         unstim      cluster2    316
```

``` r
head(selected_data$umap_plot_data)
#> # A tibble: 6 x 4
#>   cluster              stim_type UMAP1              UMAP2            
#>   <fct>                <fct>     <fct>              <fct>            
#> 1 CD3 CD8 CD127 CD45RA A         1.18542357206345   0.660537173897028
#> 2 CD3 CD8 CD127 CD45RA A         -0.844023268222809 -1.40735436424613
#> 3 CD3 CD8 CD127 CD45RA A         -0.61260942697525  -0.80203414902091
#> 4 CD3 CD8 CD127 CD45RA A         0.720577914714813  1.86169802203774 
#> 5 CD3 CD8 CD127 CD45RA A         1.13555618047714   -1.41719735607505
#> 6 CD3 CD8 CD127 CD45RA A         -1.50886527776718  -1.08557809337974
```

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

## Contact

Rohit Farmer: <rohit.farmer@gmail.com>
