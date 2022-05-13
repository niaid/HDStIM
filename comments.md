# For Testing
## Main HDStIM Function
```
mapped_data <- HDStIM(chi11$expr_data, chi11$state_markers,
                   chi11$cluster_col, chi11$stim_label,
                   chi11$unstim_label, seed_val = 123, umap = FALSE, umap_cells = NULL,
                   verbose = FALSE)
```
## Marker Ranking by Boruta
```
attribute_stats <- marker_ranking_boruta(mapped_data, path = NULL, n_cells = NULL,
                                         max_runs = 1000, seed_val = 123,
                                         verbose = FALSE)
```

# For Building
## README.Rmd
Changes made to README.Rmd are overwritten on README.md during build. Therefore, do not edit README.md directly. After making changes to README.Rmd run knitr.

## Updating Documentation
If the option to update Roxygen2 documentation is enabled in Rstudio settings then it should get triggered by default on package build/install/restart/check. However, if not then use the following command.

`devtools::document()`

## Pkgdown Website
If changes are made to README.Rmd or vignettes, build pkgdown website before pushing the code to GitHub.

`pkgdown::build_site_github_pages()`
