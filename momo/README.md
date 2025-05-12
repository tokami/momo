
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

*momo* is a new R package for tagging-based *mo*vement *mo*deling,
specifically designed to estimate fine-scale animal movement patterns
using various tagging data, such as mark-recapture, mark-resight, and
data-logging/archival tags. The movement model is based on the
advection-diffusion equation and incorporates habitat preference
functions, enabling the reconstruction of individual movement paths and
inference of population-level movement dynamics using a small set of
interpretable parameters. Movement rates are estimated as a function of
environmental conditions or geographic locations.

## Installation

The current version of *momo* (v0.0.1) requires R $>= 4.0.0$ and can be
downloaded from github as follows.

<div class=".momo-release">

``` r
# install.packages("remotes")
remotes::install_github("tokami/momo")
```

</div>

The development version can be installed with:

<div class=".momo-devel">

``` r
remotes::install_github("tokami/momo", ref = "dev")
```

</div>

## Usage

To get started with *momo*, load the package, simulate a tagging
dataset, and fit the movement model:

    #> Tags recaptured outside of spatial domain: 5. Removing them.
    #> Building the model, that can take a few minutes.
    #> Model built (0.83min). Minimizing neg. loglik.
    #>   0:    -3295.9999:  0.00000  0.00000 -4.60517 -4.60517
    #>   1:    -3343.7354: 0.252603 -0.221448 -4.44081 -4.45347
    #>   2:    -5136.5907: -0.0151293 -0.152162 -4.07660 -4.14436
    #>   3:    -5390.7222: -0.00219382 -0.173481 -3.68787 -3.75329
    #>   4:    -5395.0961: -0.0433644 -0.146010 -3.67570 -3.77448
    #>   5:    -5409.3956: -0.0209431 -0.171092 -3.64322 -3.80381
    #>   6:    -5420.5467: -0.0383744 -0.159877 -3.55973 -3.87300
    #>   7:    -5421.5147: -0.0178165 -0.175019 -3.55520 -3.87860
    #>   8:    -5423.6450: -0.0185156 -0.154777 -3.54043 -3.88729
    #>   9:    -5424.7615: -0.0206635 -0.167597 -3.52638 -3.90567
    #>  10:    -5425.0283: -0.0275750 -0.162817 -3.52337 -3.90817
    #>  11:    -5425.3096: -0.0207830 -0.162037 -3.51784 -3.91112
    #>  12:    -5425.6381: -0.0314227 -0.169760 -3.50555 -3.91558
    #>  13:    -5425.7780: -0.0142813 -0.148726 -3.48326 -3.92751
    #>  14:    -5426.0136: -0.0152051 -0.157547 -3.48227 -3.93040
    #>  15:    -5426.0942: -0.0242188 -0.157902 -3.48153 -3.93287
    #>  16:    -5426.3003: -0.0243042 -0.166718 -3.48044 -3.93588
    #>  17:    -5426.3443: -0.0262978 -0.167359 -3.47155 -3.93800
    #>  18:    -5426.3599: -0.0242056 -0.165899 -3.47177 -3.94126
    #>  19:    -5426.3641: -0.0247623 -0.165485 -3.47169 -3.94139
    #>  20:    -5426.3660: -0.0245207 -0.165744 -3.47076 -3.94273
    #>  21:    -5426.3675: -0.0244556 -0.165201 -3.46959 -3.94283
    #>  22:    -5426.3678: -0.0247224 -0.165322 -3.46953 -3.94312
    #>  23:    -5426.3681: -0.0246854 -0.165413 -3.46918 -3.94333
    #>  24:    -5426.3683: -0.0246647 -0.165477 -3.46894 -3.94395
    #>  25:    -5426.3683: -0.0247056 -0.165505 -3.46873 -3.94392
    #>  26:    -5426.3683: -0.0247165 -0.165518 -3.46881 -3.94382
    #>  27:    -5426.3683: -0.0247115 -0.165510 -3.46880 -3.94384
    #> Minimization done (0.022min). Model converged. Estimating uncertainty.

<img src="man/figures/unnamed-chunk-4-1.png" alt="Simulated data and model predictions"  />

This example illustrates the basic workflow: prepare data (here
simulate), fit the model, and access results. For real applications,
momo supports multiple types of tagging data, customizable environmental
covariates, and both estimation and prediction features.

Detailed examples and guidance are provided in the package vignettes
(see <https://tokami.github.io/momo/>).

## Getting help

You can find more information about *momo* on its *pkgdown* page at
<https://tokami.github.io/momo/>. The page includes links to vignettes,
functions descriptions, version updates, and many more. In case, your
question is not answered on the *pkgdown* webpage, please write an email
to the maintainer: [Tobias
Mildenberger](mailto:t.k.mildenberger@gmail.com). In case you find bugs,
please post an issue on [here](https://github.com/tokami/momo/issues).

## Citation

Please use the R command `citation("momo")` to receive information on
how to cite this package.

## References

1.  Mildenberger, T. K., Maunder, M., Nielsen, A. (in prep). momo: an R
    package for the estimation of fine-scale animal movement based on
    tagging data
2.  Mildenberger, T. K., et al. (in prep). Habitat suitability and
    movement of Skipjack tuna in the Eastern Pacific Ocean
