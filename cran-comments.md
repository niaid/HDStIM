# Resubmission on
## Response to the email from Gregor Seyer

> Gregor Seyer <gregor.seyer@wu.ac.at>	21 June 2022 at 04:38
> To: Rohit Farmer <rohit.farmer@gmail.com>, CRAN <cran-submissions@r-project.org>

Dear Gregor,

Thank you for your suggestions. I am resubmitting the package with the corrections. 

> If there are references describing the methods in your package, please add these in the description field of your DESCRIPTION file in the form
> authors (year) <doi:...>
> authors (year) <arXiv:...>
> authors (year, ISBN:...)
> or if those are not available: <https:...>
> with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for auto-linking.
> (If you want to add a title as well please put it in quotes: "Title")

We are still working on the manuscript. I will update the description as soon as I have the relevant information.

> Please always write package names, software names and API (application programming interface) names in single quotes in title and description. e.g: --> 'HDStIM'

Done.

> \dontrun{} should only be used if the example really cannot be executed (e.g. because of missing additional software, missing API keys, ...) by the user. That's why wrapping examples in \dontrun{} adds the comment ("# Not run:") as a warning for the user.
> Does not seem necessary.
> Please unwrap the examples if they are executable in < 5 sec, or replace \dontrun{} with \donttest{}.

I unwraped the examples from most of the functions except `marker_ranking_boruta()` and `plot_umap()` where I replaced \dontrun{} with \donttest{}.

> Please ensure that you do not use more than 2 cores in your examples, vignettes, etc.

There is only one function, "marker_ranking_boruta()" that calls the Boruta (https://cran.r-project.org/package=Boruta) package, which might be using multiple cores. However, Boruta doesn't provide a way to control the number of cores used. Please let me know if I have misunderstood this comment.

Regards,
Rohit

# Checks
## devtools::check(args = c('--as-cran'))
── R CMD check results ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── HDStIM 0.1.0 ────
Duration: 3m 49.7s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

## Win builder log
* using log directory 'd:/RCompile/CRANguest/R-release/HDStIM.Rcheck'
* using R version 4.2.0 (2022-04-22 ucrt)
* using platform: x86_64-w64-mingw32 (64-bit)
* using session charset: UTF-8
* checking for file 'HDStIM/DESCRIPTION' ... OK
* checking extension type ... Package
* this is package 'HDStIM' version '0.1.0'
* package encoding: UTF-8
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Rohit Farmer <rohit.farmer@gmail.com>'

New submission

Possibly misspelled words in DESCRIPTION:
  cytometry (21:93, 21:292)
  intracellular (21:322, 21:608)
  unstimulated (21:477)
  vitro (21:246)
* checking package namespace information ... OK
* checking package dependencies ... OK
* checking if this is a source package ... OK
* checking if there is a namespace ... OK
* checking for hidden files and directories ... OK
* checking for portable file names ... OK
* checking whether package 'HDStIM' can be installed ... OK
* checking installed package size ... OK
* checking package directory ... OK
* checking for future file timestamps ... OK
* checking 'build' directory ... OK
* checking DESCRIPTION meta-information ... OK
* checking top-level files ... OK
* checking for left-over files ... OK
* checking index information ... OK
* checking package subdirectories ... OK
* checking R files for non-ASCII characters ... OK
* checking R files for syntax errors ... OK
* checking whether the package can be loaded ... OK
* checking whether the package can be loaded with stated dependencies ... OK
* checking whether the package can be unloaded cleanly ... OK
* checking whether the namespace can be loaded with stated dependencies ... OK
* checking whether the namespace can be unloaded cleanly ... OK
* checking loading without being on the library search path ... OK
* checking use of S3 registration ... OK
* checking dependencies in R code ... OK
* checking S3 generic/method consistency ... OK
* checking replacement functions ... OK
* checking foreign function calls ... OK
* checking R code for possible problems ... [18s] OK
* checking Rd files ... [1s] OK
* checking Rd metadata ... OK
* checking Rd line widths ... OK
* checking Rd cross-references ... OK
* checking for missing documentation entries ... OK
* checking for code/documentation mismatches ... OK
* checking Rd \usage sections ... OK
* checking Rd contents ... OK
* checking for unstated dependencies in examples ... OK
* checking contents of 'data' directory ... OK
* checking data for non-ASCII characters ... OK
* checking LazyData ... OK
* checking data for ASCII and uncompressed saves ... OK
* checking installed files from 'inst/doc' ... OK
* checking files in 'vignettes' ... OK
* checking examples ... [10s] OK
* checking for unstated dependencies in 'tests' ... OK
* checking tests ... [44s] OK
  Running 'testthat.R' [44s]
* checking for unstated dependencies in vignettes ... OK
* checking package vignettes in 'inst/doc' ... OK
* checking re-building of vignette outputs ... [72s] OK
* checking PDF version of manual ... OK
* checking for detritus in the temp directory ... OK
* DONE
Status: 1 NOTE

