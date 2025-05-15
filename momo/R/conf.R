##' Generate default configuration list
##'
##' @description `def.conf` generates a list with default configurations
##'     required by the function [fit.momo].
##'
##' @param dat data frame with input data as produced by the function
##'     [check.momo.data].
##'
##' @return A list with configurations.
##'
##' @examples
##' data(skjepo)
##' conf <- def.conf(skjepo$dat)
##'
##' @export
def.conf <- function(dat){

    ## TODO: for user choices, manipulate,
    ## flags to turn stuff on, off, couple parameters, etc.
    ## everything else (data prep for RTMB -> to setup.momo.data)

    conf <- list()

    ## Flags
    conf$use.ctags <- ifelse(!is.null(dat$ctags), TRUE, FALSE)
    conf$use.atags <- ifelse(!is.null(dat$atags), TRUE, FALSE)
    conf$use.stags <- ifelse(!is.null(dat$stags), TRUE, FALSE)
    conf$use.expm <- FALSE
    conf$use.effort <- ifelse(!is.null(dat$effort), TRUE, FALSE)
    conf$use.catch <- FALSE
    conf$use.taxis <- TRUE
    conf$use.advection <- FALSE
    conf$use.length <- FALSE
    conf$est.var.ctags <- FALSE
    conf$est.var.atags <- FALSE
    conf$est.n <- FALSE
    conf$pred.move <- TRUE
    conf$use.boundaries <- FALSE
    conf$use.rel.events <- FALSE

    ## Env
    if(!is.null(dat$env)){
        nenv <- length(dat$env)
        ## OLDER:
        ## ## TODO: how to do this automatic? env could be annually, by quarter, season, require a time stamp? attributes(dat$env)
        ## ## ienv <- rep(sapply(dat$env, function(x) dim(x)[3]),
        ## ##             each = length(dat$time.cont))
        ## ienv <- sapply(dat$env, function(x){
        ##     if(dim(x)[3] == diff(dat$trange)){
        ##         res <- rep(1:dim(x)[3], each = (length(dat$time.cont)-1)/dim(x)[3])
        ##         if(length(res) != length(dat$time.cont)){
        ##             res <- c(res, rep(dim(x)[3], length(dat$time.cont) - length(res)))
        ##         }
        ##     }else{
        ##         res <- rep(1, length(dat$time.cont))
        ##     }
        ##     return(res)
        ## })
        ## ienv.tax <- matrix(ienv, nenv, length(dat$time.cont))
        tmp <- sapply(dat$env,
                      function(x)
                          if(!is.null(attributes(x)$dimnames[[3]])){
                              as.numeric(attributes(x)$dimnames[[3]])
                          }else{
                              seq(dat$trange[1], dat$trange[2],
                                  length.out = dim(x)[3]+1)
                          })
        ienv.tax <- t(apply(tmp, 2,
                          function(x){
                              tmp <- as.integer(cut(dat$time.cont,
                                                    unique(c(x,dat$trange[2])),
                                                    right = FALSE,
                                                    include.lowest = TRUE))
                              tmp[is.na(tmp)] <- min(tmp, na.rm = TRUE)
                              return(tmp)
                          }))
        if(dat$const.dif){
            ienv.dif <- matrix(1, nenv, length(dat$time.cont))
        }else{
            ienv.dif <- matrix(ienv, nenv, length(dat$time.cont))
        }

        ienv.adv.x <- matrix(0, nenv, length(dat$time.cont))
        ienv.adv.y <- matrix(0, nenv, length(dat$time.cont))
        conf$ienv <- list(tax   = ienv.tax,
                          dif   = ienv.dif,
                          adv.x = ienv.adv.x,
                          adv.y = ienv.adv.y)
    }


    ## Seasonally varying prefs
    if(!is.null(dat$env)){

        nenv <- length(dat$env)
        ienvS <- matrix(1, nenv, length(dat$time.cont))

        conf$ienvS <- list(tax   = ienvS,
                           dif   = ienvS,
                           adv = ienvS)
    }

    ## REMOVE:
    ## if(!is.null(dat$ctags)){
    ##     rec <- !is.na(dat$ctags$x1) & !is.na(dat$ctags$y1) & !is.na(dat$ctags$t1)
    ## }else{
    ##     rec <- NULL
    ## }
    ## conf$rec <- rec


    ## Boundaries
    ## tmp <- array(1, c(dat$nx+2, dat$ny+2,1))
    ## tmp[c(1,dat$nx+2),,1] <- 0
    ## tmp[,c(1,dat$ny+2),1] <- 0
    tmp <- array(1, c(dat$nx, dat$ny,1))
    tmp[c(1,dat$nx),,1] <- 0
    tmp[,c(1,dat$ny),1] <- 0
    conf$boundaries <- list(tmp)

    conf$ibound <- matrix(1, 1, length(dat$time.cont))


    return(conf)
}
