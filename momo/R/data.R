##' Simulated tags
##'
##' @name simTags
##' @docType data
##' @author Tobias Mildenberger \email{tobm@dtu.dk}
##' @references \url{bla bla}
##' @keywords data
NULL


##' prep.ctags
##'
##' @param x Data frame with required information
##' @param names Names of columns that contain in following order: release time (t0), recapture time (t1), release location x (x0), recapture location x (x1), release location y (y0), recapture location y (y1)
##' @param origin Optional, if time column not provided as decimal date, origin can be used to convert date columns to dates
##' @param speed.limit Optional speed in km/d
##' @param verbose Print messages?
##'
##' @export
prep.ctags <- function(x,
                       names,
                       origin = NULL,
                       speed.limit = NULL,
                       verbose = TRUE){

    ## Select columns
    idx <- names %in% colnames(x)
    if(sum(idx) != length(names)){
        stop(paste0(paste(names[!idx],collapse = ","),
                    " not found in x!"))
    }
    res <- x[,match(names,colnames(x))]
    colnames(res) <- c("t0","t1","x0","x1","y0","y1")

    ## Convert dates
    if(!is.null(origin)){
        res$t0 <- lubridate::decimal_date(as.Date(res$t0, origin = origin))
        res$t1 <- lubridate::decimal_date(as.Date(res$t1, origin = origin))
    }

    ## Remove tags without release location
    idx <- which(is.na(res$x0) | is.na(res$y0))
    if(length(idx) > 0){
        if(verbose) writeLines(paste0("Removing ", length(idx),
                          "tags without release location."))
        res <- res[-idx,]
    }

    ## Remove tags that are recaptured but miss date or location info
    recaptured <- !is.na(res$t1) | !is.na(res$x1) | !is.na(res$y1)
    idx <- which(recaptured & (is.na(res$t1) |
                               is.na(res$x1) |
                               is.na(res$y1)))
    if(length(idx) > 0){
        if(verbose) writeLines(paste0("Removing ", length(idx),
                          " recaptured tags with incomplete recapture information (date, location)."))
        res <- res[-idx,]
    }

    ## Omit tags with t0 >= t1
    idx <- which(!is.na(res$t1) & res$t0 >= res$t1)
    if(length(idx) > 0){
        if(verbose) writeLines(paste0("Removing ", length(idx)," recaptured tags for which release time is after recapture time (t0 > t1)."))
        res <- res[-idx,]
    }


    ## Speed filter
    if(!is.null(speed.limit)){
        tdiff <- as.numeric(res$t1 - res$t0) * 365
        dist <- rep(NA, nrow(res))
        idx <- which(!is.na(res$x1) & !is.na(res$y1))
        for(i in 1:length(idx)){
            dist[idx[i]] <- geosphere::distm(cbind(res$x0[idx[i]], res$y0[idx[i]]),
                                                 cbind(res$x1[idx[i]], res$y1[idx[i]]),
                                                 fun = geosphere::distGeo) / 1000 ## km
        }
        speed <- dist / tdiff

        idx <- which(speed > speed.limit)
        if(length(idx) > 0){
            writeLines(paste0("Removing ", length(idx)," recaptured tags that are faster than the provided speed limit(",speed.limit," km/d)."))
            res <- res[-idx,]
        }
    }

    res <- add.class(res, "momo.ctags")

    return(res)
}



##' prep.atags
##'
##' @param x Data frame with required information
##' @param names Names of columns that contain in following order: release time (t0), recapture time (t1), release location x (x0), recapture location x (x1), release location y (y0), recapture location y (y1)
##' @param origin Optional, if time column not provided as decimal date, origin can be used to convert date columns to dates
##' @param speed.limit Optional speed in km/d
##' @param verbose Print messages?
##'
##' @export
prep.atags <- function(x,
                       names,
                       origin = NULL,
                       speed.limit = NULL,
                       verbose = TRUE){


    x0 <- x
    x <- lapply(seq_along(x), function(i) {
        x[[i]]$id <- i
        return(x[[i]])
    })
    x <- do.call(rbind, x)


    ## Select columns
    idx <- names %in% colnames(x)
    if(sum(idx) != length(names)){
        stop(paste0(paste(names[!idx],collapse = ","),
                    " not found in x!"))
    }
    dat <- x[,match(names,colnames(x))]
    colnames(dat) <- c("t","x","y")
    id <- x$id

    ## Convert dates
    if(!is.null(origin)){
        dat$t <- lubridate::decimal_date(as.Date(dat$t, origin = origin))
    }

    res <- split(dat, id)

    res <- add.class(res, "momo.atags")

    return(res)
}


##' prep.env
##'
##' @param x Data frame with required information
##' @param names Names of columns that contain in following order: release time (t0), recapture time (t1), release location x (x0), recapture location x (x1), release location y (y0), recapture location y (y1)
##' @param origin Optional, if time column not provided as decimal date, origin can be used to convert date columns to dates
##' @param verbose Print messages?
##'
##' @export
prep.env <- function(x,
                     names,
                     origin = NULL,
                     verbose = TRUE){

    res <- x

    if(inherits(x, "list")){
        res <- simplify2array(x)
    }

    if(length(dim(x)) == 2){
        res <- array(x, c(dim(x),1))
    }

    ## if(length(dim(x)) == 3){
    ##     res <- lapply(seq(dim(x)[3]), function(i) x[,,i])
    ## }

    res <- add.class(res, "momo.env")

    return(res)
}


##' prep.effort
##'
##' @param x Data frame with required information
##' @param names Names of columns that contain in following order: release time (t0), recapture time (t1), release location x (x0), recapture location x (x1), release location y (y0), recapture location y (y1)
##' @param origin Optional, if time column not provided as decimal date, origin can be used to convert date columns to dates
##' @param verbose Print messages?
##'
##' @export
prep.effort <- function(x,
                     names,
                     origin = NULL,
                     verbose = TRUE){

    res <- x

    if(inherits(x, "list")){
        res <- simplify2array(x)
    }

    if(length(dim(x)) == 2){
        res <- array(x, c(dim(x),1))
    }

    ## if(length(dim(x)) == 3){
    ##     res <- lapply(seq(dim(x)[3]), function(i) x[,,i])
    ## }

    res <- add.class(res, "momo.effort")

    return(res)
}



##' setup.momo.data
##'
##' @export
setup.momo.data <- function(grid,
                            env,
                            ctags = NULL,
                            atags = NULL,
                            effort = NULL,
                            knots.tax = NULL,
                            knots.dif = NULL,
                            const.dif = TRUE,
                            trange = NULL,
                            dt = NULL,
                            time.cont = NULL,
                            verbose = TRUE){

    res <- list()

    env <- check.that.list(env)

    ## Time
    if(is.null(trange)){
        if(!is.null(ctags) || !is.null(atags)){
            dim.tags <- get.dim(ctags, atags)
            trange <- dim.tags$trange
        }else{
            trange <- c(0, max(1,max(sapply(env, function(x) dim(x)[3]))-1))
        }
    }
    res$trange[1] <- floor(trange[1])
    res$trange[2] <- ceiling(trange[2])
    res$dt <- ifelse(is.null(dt), 0.1, dt)
    if(is.null(time.cont)){
        res$time.cont <- seq(res$trange[1],
                             res$trange[2],
                             res$dt)
    }else{
        res$time.cont <- time.cont
    }


    nt <- length(res$time.cont)

    ## Space
    res$xygrid <- grid$xygrid
    res$igrid <- grid$igrid
    res$celltable <- grid$celltable
    res$xrange <- attr(grid, "xrange")
    res$yrange <- attr(grid, "yrange")
    res$dxdy <- attr(grid, "dxdy")
    res$xgr <- attr(grid, "xgr")
    res$ygr <- attr(grid, "ygr")
    res$xcen <- attr(grid, "xcen")
    res$ycen <- attr(grid, "ycen")
    res$nx <- attr(grid, "nx")
    res$ny <- attr(grid, "ny")
    res$nextTo <- get.neighbours(res$celltable)
    res$next.dist <- c(res$dxdy[1],res$dxdy[1],
                       res$dxdy[2],res$dxdy[2])
    ## TODO: use sum(res$dxdy^2)/2 for all diagonal neighbours

    ## Env
    res$env <- env
    if(!is.null(env)){
        res$xranges <- matrix(res$xrange, nrow = length(env), ncol = 2, byrow = TRUE)
        res$yranges <- matrix(res$yrange, nrow = length(env), ncol = 2, byrow = TRUE)
    }

    ## ctags
    res$ctags <- ctags

    if(!is.null(ctags)){

        ## Checks
        ## TODO: outsource in check function?

        ## Remove tags released outside of spatial domain
        idx <- which(res$ctags$x0 < res$xrange[1] |
                     res$ctags$x0 > res$xrange[2] |
                     res$ctags$y0 < res$yrange[1] |
                     res$ctags$y0 > res$yrange[2])
        if(length(idx) > 0){
            if(verbose) writeLines(paste0("Tags released outside of spatial domain: ", length(idx), ". Removing them."))
            res$ctags <- res$ctags[-idx,]
        }

        ## Remove tags recpatured outside of spatial domain
        idx <- which(!is.na(res$ctags$x1) & (res$ctags$x1 < res$xrange[1] |
                                             res$ctags$x1 > res$xrange[2]) |
                     !is.na(res$ctags$y1) & (res$ctags$y1 < res$yrange[1] |
                                             res$ctags$y1 > res$yrange[2]))
        if(length(idx) > 0){
            if(verbose) writeLines(paste0("Tags recaptured outside of spatial domain: ", length(idx),". Removing them."))
            res$ctags <- res$ctags[-idx,]
        }

        ## Remove tags that are in grid = NA
        idx <- which(is.na(res$celltable[cbind(cut(res$ctags$x0, res$xgr),
                                               cut(res$ctags$y0, res$ygr))]))
        if(length(idx) > 0){
            if(verbose) writeLines(paste0("Tags released in grid cells that are NA: ", length(idx),". Removing them."))
            res$ctags <- res$ctags[-idx,]
        }
        idx <- which(!is.na(res$ctags$x1) & !is.na(res$ctags$y1) &
                     is.na(res$celltable[cbind(cut(res$ctags$x1, res$xgr),
                                               cut(res$ctags$y1, res$ygr))]))
        if(length(idx) > 0){
            if(verbose) writeLines(paste0("Tags recaptured in grid cells that are NA: ", length(idx),". Removing them."))
            res$ctags <- res$ctags[-idx,]
        }

        ## Remove tags outside of time period
        idx <- which(res$ctags$t0 < res$trange[1] |
                     res$ctags$t0 > res$trange[2] |
                     (!is.na(res$ctags$t1) &
                      (res$ctags$t1 < res$trange[1] |
                       res$ctags$t1 > res$trange[2])))
        if(length(idx) > 0){
            if(verbose) writeLines(paste0("Tags released or recaptured before or after time period: ", length(idx),". Removing them."))
            res$ctags <- res$ctags[-idx,]
        }

        ## CHECK: same time step recovery allowed for expm?
        res$ctags <- res$ctags[is.na(res$ctags$t1) |
                               res$ctags$t0 != res$ctags$t1,]
        ## TODO: cut to time step than check if equal?

        ## Remove all non recovered if effort is zero
        if(is.null(effort)){
            idx <- which(is.na(res$ctags$t1) |
                         is.na(res$ctags$x1) |
                         is.na(res$ctags$y1))
            if(length(idx) > 0){
                res$ctags <- res$ctags[-idx,]
            }
        }

        ## Discretised
        res$ctags$itrel <- as.integer(cut(res$ctags$t0, res$time.cont,
                                          include.lowest = TRUE))
        res$ctags$itrec <- as.integer(cut(res$ctags$t1, res$time.cont,
                                          include.lowest = TRUE))
        res$ctags$icrel <- res$celltable[cbind(as.integer(cut(res$ctags$x0, res$xgr)),
                                               as.integer(cut(res$ctags$y0, res$ygr)))]
        res$ctags$icrec <- res$celltable[cbind(as.integer(cut(res$ctags$x1, res$xgr)),
                                               as.integer(cut(res$ctags$y1, res$ygr)))]

        ## Remove tags recovered in same time step as released
        ind <- which(res$ctags$itrec - res$ctags$itrel <= 0)
        if(length(ind) > 0){
            res$ctags <- res$ctags[-ind,]
        }
        ## FUTURE: not needed for KF, could remove quite a few tags!
        ## make flag to swith this off, possible to have this later and only use when running expm?

        ## Release events
        tmp <- get.release.events(res)
        res$rel.events <- tmp$rel.events
        res$ctags$rel.event <- tmp$idx

    }


    ## atags
    res$atags <- atags

    if(!is.null(atags)){

        na <- length(atags)
        res$atags <- vector("list", na)
        for(a in 1:na){
            tag <- as.data.frame(atags[[a]])
            tag$it <- as.integer(cut(tag$t, res$time.cont,
                                     include.lowest = TRUE))
            tag$ic <- res$celltable[cbind(as.integer(cut(tag$x, res$xgr)),
                                          as.integer(cut(tag$y, res$ygr)))]
            res$atags[[a]] <- tag
        }

        ## Checks
        ## Deactivate any archival tag that leaves the grid (NA in ic)
        for(a in 1:na){
            tag <- as.data.frame(res$atags[[a]])
            if(any(is.na(tag$ic))){
                tag$use <- 0
            }else{
                tag$use <- 1
            }
            res$atags[[a]] <- tag
        }

    }

    res$len.mins <- matrix(0, 2, 2)

    ## Effort
    res$effort <- effort

    if(!is.null(res$effort)){

        if(!inherits(res$effort, "list")) res$effort <- list(res$effort)

        res$xranges.eff <- matrix(res$xrange, nrow = length(res$effort),
                                  ncol = 2, byrow = TRUE)
        res$yranges.eff <- matrix(res$yrange, nrow = length(res$effort),
                                  ncol = 2, byrow = TRUE)
        res$ne <- 1
        res$time.cont.eff <- rep(0, 10)
        res$nm <- 1
    }else{
        res$xranges.eff <- matrix(res$xrange, nrow = 1, ncol = 2, byrow = TRUE)
        res$yranges.eff <- matrix(res$yrange, nrow = 1, ncol = 2, byrow = TRUE)

        res$ne <- 1
        res$time.cont.eff <- rep(0, 10)
        res$nm <- 1
    }

    ieff <- rep(sapply(res$effort, function(x) dim(x)[3]),
                    each = length(res$time.cont))
    res$ieff <- matrix(ieff, length(res$effort), length(res$time.cont))

    res$log2steps <- 0


    res$const.dif <- const.dif

    ## Knots
    res$knots.tax <- knots.tax
    res$knots.dif <- knots.dif
    ## env.obs <- get.env(dat, conf)
    ## if(is.null(env.obs)){
    ##     env.obs <- dat$env
    ## }
    env.obs <- env
    if(is.null(res$knots.tax)){
        res$knots.tax <- sapply(env.obs,
                                 function(x)
                                     quantile(as.numeric(x),
                                              c(0.05, 0.5, 0.95), na.rm = TRUE))
    }

    if(any(apply(res$knots.tax,2,duplicated))) warning("Some knots are the same! This will likely give an error!")

    if(is.null(res$knots.dif)){
        if(const.dif){
            res$knots.dif <- matrix(0,
                                     1, ## constant diffusion by default
                                     length(env))
        }else{
            res$knots.dif <- sapply(env,
                                     function(x)
                                         quantile(as.numeric(x),
                                                  c(0.05, 0.5, 0.95), na.rm = TRUE))
        }
    }


    ## Prediction
    res$xygrid.pred <- res$xygrid
    res$igrid.pred <- res$igrid
    res$time.cont.pred <- res$time.cont
    res$env.pred <- sapply(env, function(x) seq(min(x, na.rm = TRUE),
                                                max(x, na.rm = TRUE),
                                                length.out = 100))

    res$ddt <- 0.01
    res$eps <- 0.000001

    ## Return
    res <- add.class(res, "momo.data")
    return(res)
}
