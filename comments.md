# For Building and Testing
## Main HDStIM Function
mapped_data <- HDStIM(chi11$expr_data, chi11$state_markers,
                   chi11$cluster_col, chi11$stim_label,
                   chi11$unstim_label, seed_val = 123, umap = FALSE, umap_cells = NULL,
                   verbose = FALSE)

## Pkgdown Website
pkgdown::build_site_github_pages()
