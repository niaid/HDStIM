## code to prepare `DATASET` dataset goes here

chi11_1k <- readRDS(file.path("data-raw", "chi11_1k.rds"))
usethis::use_data(chi11_1k, compress = "xz")
