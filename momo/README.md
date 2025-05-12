
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

## Overview

*momo* is a new R package for tagging-based *mo*vement *mo*deling,
specifically designed to estimate fine-scale animal movement patterns
using various tagging data, such as mark-recapture, mark-resight, and
data-logging/archival tags. The movement model is based on the
advection-diffusion equation and incorporates habitat preference
functions, enabling the reconstruction of individual movement paths and
inference of population-level movement dynamics using a small set of
interpretable parameters. Movement rates are estimated as a function of
environmental conditions or geographic locations.

The current version of the package allows to estimate habitat preference
functions and movement patterns. Further package development will allow
to estimate spatiotemporally varying recapture probability, estimate
natural and fishing mortality rates, and relative biomass indices and
distribution.

To get started with *momo*, install and load the package, simulate a
tagging dataset, and fit the movement model:

``` r
## devtools::install_github("tokami/momo/momo")
library(momo)

## Simulate a small tagging data set
sim <- sim.momo()

## Fit momo to the simulated data
fit <- fit.momo(sim)

## Plot simulated data and model predictions
plotmomo.compare(sim = sim, fit = fit)
```

## More information

More detailed examples and documentation for *momo* can be found at
<https://tokami.github.io/momo/>. The *pkgdown* page includes links to
articles, vignettes, functions descriptions, information to version
updates, and much more. In case, your question is not answered by the
package documentation and on the *pkgdown* pages, please write an email
to the maintainer: [Tobias
Mildenberger](mailto:t.k.mildenberger@gmail.com). In case you find bugs,
please post an issue on [here](https://github.com/tokami/momo/issues).

### Citation

Please use the R command `citation("momo")` to receive information on
how to cite this package.

### Funding

The development of *momo* was cofunded by the European Union.

------------------------------------------------------------------------

<img src="man/figures/EN_Co-fundedbytheEU_RGB_POS.png" width="250" />
