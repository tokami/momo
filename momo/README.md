
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

momo is a new R package for tagging-based **mo**vement **mo**deling,
specifically designed to estimate fine-scale animal movement patterns
using various tagging data, such as mark-recapture, mark-resight, and
data-logging/archival tags. The movement model is based on the
advection-diffusion equation and incorporates habitat preference
functions, enabling the reconstruction of individual movement paths and
inference of population-level movement dynamics using a small set of
interpretable parameters. Movement rates are estimated as a function of
environmental conditions or geographic locations.

## Installation

<div class=".momo-release">

``` r
# install.packages("remotes")
remotes::install_github("tokami/momo")
```

</div>

<div class=".momo-devel">

``` r
# install.packages("remotes")
remotes::install_github("tokami/momo", ref = "dev")
```

</div>

## Usage

Get started with:

## Citation

Please use the R command `citation("momo")` to receive information on
how to cite this package.

## Documentation

The [tutorial](link) vignette demonstrates the use of the main functions
of momo for estimating habitat suitability and movement patterns for
simulated example data.

## Questions / Issues

In case you have questions or find bugs, please write an email to
[Tobias Mildenberger](mailto:t.k.mildenberger@gmail.com) or post on
[momo/issues](https://github.com/tokami/momo/issues).

## References

1.  Mildenberger, T. K., Maunder, M., Nielsen, A. (in process). momo: an
    R package for the estimation of fine-scale animal movement based on
    tagging data [link](link)
2.  Mildenberger, T. K., et al. (in process). Habitat suitability and
    movement of Skipjack tuna in the Eastern Pacific Ocean [link](link)
