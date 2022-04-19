
<!-- README.md is generated from README.Rmd. Please edit that file -->

# HDStIM <img src='man/figures/sticker.png' align="right" height = 200 />

<!-- badges: start -->

[![Actions
Status](https://github.com/niaid/HDStIM/workflows/R-CMD-check/badge.svg)](https://github.com/niaid/HDStIM/actions?query=workflow%3AR-CMD-check)
<!-- badges: end -->

The goal of this package is to identify response to a stimulant in
CyTOF/Flow cytometry stimulation assays by labeling cells as responded
or not based on an unsupervised high dimensional approach. Starting from
the annotated cell populations either through automated clustering such
as FlowSOM or traditional cell gating, the primary function `HDStIM()`
follows a heuristic approach to label cells as responding or
non-responding.

For a combination of cell population and stimulation type (e.g., CD127+
T-helper cells and interferon-alpha), `HDStIM()` starts by performing
k-means clustering on the combined set of cells from stimulated and
unstimulated samples. K-means clustering is performed on expression data
of all the state markers combined. Upon clustering using a contingency
table, a Fisher’s exact test determines the effect size and the
statistical significance of partitioning. Cells form the combinations
that pass the Fisher’s exact test are labelled as responding.

## Installation

You can install the released version of stimcellselector from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("HDStIM")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("niaid/HDStIM")
```

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
