##' def.conf
##'
##' @param dat Dat frame with required information created by function \code{\link{check.momo.data}}
##' @param const.dif Constant diffusion? logical
##'
##' @export
def.conf <- function(dat, const.dif = TRUE){

    ## TODO: for user choices, manipulate,
    ## flags to turn stuff on, off, couple parameters, etc.
    ## everything else (data prep for RTMB -> to setup.momo.data)

    conf <- list()

    ## Flags
    conf$use.ctags <- ifelse(!is.null(dat$ctags), TRUE, FALSE)
    conf$use.atags <- ifelse(!is.null(dat$atags), TRUE, FALSE)
    conf$use.kf <- TRUE
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

    ## Env
    if(!is.null(dat$env)){
        nenv <- length(dat$env)
        ## TODO: how to do this automatic? env could be annually, by quarter, season, require a time stamp? attributes(dat$env)
        ## ienv <- rep(sapply(dat$env, function(x) dim(x)[3]),
        ##             each = length(dat$time.cont))
        ienv <- sapply(dat$env, function(x){
            if(dim(x)[3] == diff(dat$trange)){
                res <- rep(1:dim(x)[3], each = (length(dat$time.cont)-1)/dim(x)[3])
                if(length(res) != length(dat$time.cont)) res <- c(res, rep(dim(x)[3], length(dat$time.cont) - length(res)))
            }else{
                res <- rep(1, length(dat$time.cont))
            }
            return(res)
        })
        ienv.tax <- matrix(ienv, nenv, length(dat$time.cont))

        if(const.dif){
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


    ## Tags
    if(!is.null(dat$ctags)){
        itrec <- as.integer(cut(dat$ctags$t1, dat$time.cont,
                               include.lowest = TRUE))
        icrec <- dat$celltable[cbind(as.integer(cut(dat$ctags$x1, dat$xgr)),
                                     as.integer(cut(dat$ctags$y1, dat$ygr)))]

    }else{
        itrec <- NULL
        icrec <- NULL
    }
    conf$itrec <- itrec
    conf$icrec <- icrec

    ## Release events
    if(!is.null(dat$ctags)){
        itrel <- as.integer(cut(dat$ctags$t0, dat$time.cont,
                               include.lowest = TRUE))
        icrel <- dat$celltable[cbind(as.integer(cut(dat$ctags$x0, dat$xgr)),
                                     as.integer(cut(dat$ctags$y0, dat$ygr)))]
        id <- paste0(itrel,":",icrel)
        uni <- unique(id)
        irel.event <- match(id, uni)
        rel.events <- data.frame(itrel = as.numeric(sapply(strsplit(uni, ":"),"[[",1)),
                                 icrel = as.numeric(sapply(strsplit(uni, ":"),"[[",2)),
                                 itrec = as.numeric(by(itrec, irel.event, max)))
        rownames(rel.events) <- NULL
    }else{
        rel.events <- NULL
        irel.event <- NULL
    }
    conf$rel.events <- rel.events
    conf$irel.event <- irel.event


    if(!is.null(dat$ctags)){
        excl.ctags <- rep(0, nrow(dat$ctags))
        excl.ctags[(itrec - itrel) <= 0] <- 1
    }else{
        excl.ctags <- NULL
    }
    conf$excl.ctags <- excl.ctags

    if(!is.null(dat$ctags)){
        rec <- !is.na(dat$ctags$x1) & !is.na(dat$ctags$y1) & !is.na(dat$ctags$t1)
    }else{
        rec <- NULL
    }
    conf$rec <- rec


    ## Knots
    env.obs <- get.env(dat, conf)
    if(is.null(env.obs)){
        env.obs <- dat$env
    }
    conf$knots.tax <- sapply(env.obs,
                             function(x)
                                 quantile(as.numeric(x),
                                          c(0.05, 0.5, 0.95), na.rm = TRUE))

    if(any(apply(conf$knots.tax,2,duplicated))) warning("Some knots are the same! This will likely give an error!")

    if(const.dif){
        conf$knots.dif <- matrix(0,
                                 1, ## constant diffusion by default
                                 length(dat$env))
    }else{
        conf$knots.dif <- sapply(dat$env,
                                 function(x)
                                     quantile(as.numeric(x),
                                              c(0.05, 0.5, 0.95), na.rm = TRUE))
    }


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
