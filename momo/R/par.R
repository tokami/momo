
##' def.par
##'
##' @param dat Data list with required information
##' @param conf Configuration list
##'
##' @export
def.par <- function(dat, conf){

    par <- list()

    ## Taxis
    par$alpha <- matrix(rep(0, length(conf$knots.tax)),
                        nrow(conf$knots.tax),
                        ncol(conf$knots.tax))

    ## Diffusion
    par$beta <- matrix(rep(log(0.01), length(conf$knots.dif)),
                        nrow(conf$knots.dif),
                       ncol(conf$knots.dif))

    ## Advection
    par$gamma <- matrix(rep(0, length(dat$env)),
                        2,  ## 1: x direction; 2: y direction
                        length(dat$env))

    ## Observation uncertainty
    par$logSdObsATS <- log(0.01)

    ## Survival
    par$logLambda <- matrix(0, length(dat$effort), 1)
    par$lambdaEC <- matrix(0, length(dat$effort), 1)
    par$logM <- log(0.1)

    return(par)
}



##' def.map
##'
##' @param dat Data list with required information
##' @param conf Configuration list
##' @param par parameter list
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

    ## Survival
    if(is.null(dat$effort) || !conf$use.effort){
        map$logM <- factor(rep(NA, length(par$logM)))
        map$logLambda <- factor(rep(NA, length(par$logLambda)))
    }
    map$lambdaEC <- factor(rep(NA, length(par$lambdaEC)))


return(map)
}
