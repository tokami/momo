##' def.par
##'
##' @description Get a list with parameters and default initial values.
##'
##' @param dat data frame with input data as produced by the function
##'     [check.momo.data].
##' @param conf configuration list as produced by the function [def.conf].
##'
##' @return A named list.
##'
##' @examples
##' data(skjepo)
##'
##' par <- with(skjepo, def.par(dat, conf))
##'
##' @export
def.par <- function(dat, conf){

    par <- list()

    ## Taxis
    ## par$alpha <- matrix(rep(0, length(dat$knots.tax)),
    ##                     nrow(dat$knots.tax),
    ##                     ncol(dat$knots.tax))

    ## TRY:
    par$alpha <- array(rep(0, length(dat$knots.tax)),
                       dim = c(nrow(dat$knots.tax),
                               ncol(dat$knots.tax),
                               max(conf$ienvS$tax)))

    ## Diffusion
    ## par$beta <- matrix(rep(log(0.01), length(dat$knots.dif)),
    ##                     nrow(dat$knots.dif),
    ##                    ncol(dat$knots.dif))

    ## TRY:
    par$beta <- array(rep(log(0.01), length(dat$knots.dif)),
                        dim = c(nrow(dat$knots.dif),
                                ncol(dat$knots.dif),
                                max(conf$ienvS$dif)))

    ## Advection
    ## par$gamma <- matrix(rep(0, length(dat$env)),
    ##                     2,  ## 1: x direction; 2: y direction
    ##                     length(dat$env))

    ## TRY:
    par$gamma <- array(rep(0, length(dat$env)),
                        dim = c(2,  ## 1: x direction; 2: y direction
                                length(dat$env),
                                max(conf$ienvS$adv)))

    ## Observation uncertainty
    par$logSdObsATS <- log(0.01)
    par$logSdObsSTS <- log(0.01)

    ## Survival
    par$logLambda <- matrix(0, length(dat$effort), 1)
    par$lambdaEC <- matrix(0, length(dat$effort), 1)
    par$logM <- log(0.1)

    return(par)
}



##' def.map
##'
##' @description Get a list with default mappings for parameters.
##'
##' @param dat data frame with input data as produced by the function
##'     [check.momo.data].
##' @param conf configuration list as produced by the function [def.conf].
##' @param par parameter list with initial values as produced by the function
##'     [def.par].
##'
##' @examples
##' data(skjepo)
##'
##' map <- with(skjepo, def.map(dat, conf, par))
##'
##' @export
def.map <- function(dat, conf, par){

    map <- list()

    ## Taxis
    map$alpha <- get.map(as.numeric(apply(par$alpha, 2,
                                          function(x) c(NA, rep(1, length(x)-1)))))

    ## Diffusion
    tmp <- apply(par$beta, 2,
                 function(x) c(NA, rep(1, length(x)-1)))
    tmp <- matrix(tmp, nrow(par$beta), ncol(par$beta))
    tmp[1,1] <- 1 ## one intercept
    map$beta <- get.map(as.numeric(tmp))

    ## Advection
    if(conf$use.advection){
        ## Link x and y direction by default
        map$gamma <- factor(sapply(1:ncol(par$gamma),
                                   function(x) rep(x,nrow(par$gamma))))
    }else{
        map$gamma <- factor(rep(NA, length(par$gamma)))
    }

    ## Observation error
    if(!conf$use.atags){
        map$logSdObsATS <- factor(rep(NA, length(par$logSdObsATS)))
    }

    if(!conf$use.stags){
        map$logSdObsSTS <- factor(rep(NA, length(par$logSdObsSTS)))
    }

    ## Survival
    if(is.null(dat$effort) || !conf$use.effort){
        map$logM <- factor(rep(NA, length(par$logM)))
        map$logLambda <- factor(rep(NA, length(par$logLambda)))
    }
    map$lambdaEC <- factor(rep(NA, length(par$lambdaEC)))


return(map)
}




##' Lower bounds
##'
##' @description Get lower bounds for parameters.
##'
##' @param par parameter list with initial values as produced by the function
##'     [def.par].
##'
##' @return List with lower bounds for parameters.
get.lower.bounds <- function(par){
    list()
}



##' Upper bounds
##'
##' @description Get upper bounds for parameters.
##'
##' @param par parameter list with initial values as produced by the function
##'     [def.par].
##'
##' @return List with upper bounds for parameters.
get.upper.bounds <- function(par){
    list()
}
