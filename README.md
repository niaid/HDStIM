
<!-- README.md is generated from README.Rmd. Please edit that file -->

# stimcellselector <a href='https://niaid.github.io/stimcellselector/'><img src='man/figures/sticker.png' align="right" height="200" /></a>

<!-- badges: start -->

[![Actions
Status](https://github.com/niaid/stimcellselector/workflows/R-CMD-check/badge.svg)](https://github.com/niaid/stimcellselector/actions?query=workflow%3AR-CMD-check)
[![Codecov test
coverage](https://codecov.io/gh/rohitfarmer/stimcellselector/branch/master/graph/badge.svg?token=IXR3EVFLXA)](https://codecov.io/gh/rohitfarmer/stimcellselector?branch=master)
[![](https://img.shields.io/github/languages/code-size/rohitfarmer/stimcellselector.svg)](https://github.com/niaid/stimcellselector)
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
responding.

## Installation

You can install the released version of stimcellselector from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("stimcellselector")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("niaid/stimcellselector")
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
```

For more information on the interpredation of the output and how to
generate relevant figures please visit the package website at
<https://niaid.github.io/stimcellselector/>.

## Contact

Rohit Farmer: <rohit.farmer@nih.gov>

## Disclaimer

THE SOFTWARE IS PROVIDED “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT
ARE DISCLAIMED. IN NO EVENT SHALL THE NATIONAL CANCER INSTITUTE (THE
PROVIDER), THE NATIONAL INSTITUTES OF HEALTH, THE U.S. GOVERNMENT OR THE
INDIVIDUAL DEVELOPERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
