##' Full habitat function
##'
##' @description Creates a light habitat class that is being used within
##'     \emph{momo}.
##'
##' @param FIELDS A list with 3-D array fields that are used for interpolation,
##'     where the first 2 dimensions span the x and y direction of the spatial
##'     field and the third dimension indicates the time dimension.
##' @param XR Limits of the 2-D fields in x direction.
##' @param YR Limits of the 2-D fields in y direction.
##' @param ienv Indicator matrix mapping each model time step to the time steps
##'     of each field. The first dimension corresponds to the number of fields
##'     and the second to the number of model time steps.
##' @param time.cont Vector with continuous model time steps.
##' @param S List of smooth functions relating fields to habitat preference
##'     functions.
##' @param dS List of derivatives of smooth functions relating fields to habitat
##'     preference functions.
##' @param ienvS Indicator matrix mapping each model time step to the time step
##'     of each spline. The first dimension corresponds to the number of splines
##'     and the second to the number of model time steps.
##'
##' @return List with functions.
##'
##' @importFrom abind abind
##'
habi.full <- function(FIELDS, XR, YR, ienv, time.cont, S, dS, ienvS){

    nenv <- length(FIELDS)

    dxfield <- function(field, xr){
        nr <- nrow(field)
        nc <- ncol(field)
        ret <- left <- right <- matrix(NA_real_, nr, nc)
        onedx <- (xr[2] - xr[1]) / nr ## before / (nr - 1) WHY:?
        ret[2:(nr-1),] <- (field[3:nr,] - field[1:(nr-2),]) / (2 * onedx)
        ## Account for edges and NA in maps (land, islands)
        left[1:(nr-1),] <- (field[2:nr,] - field[1:(nr-1),]) / onedx
        right[2:nr,] <- (field[2:nr,] - field[1:(nr-1),]) / onedx
        comb <- rowMeans(abind::abind(left, right, along = 3), dims = 2, na.rm = TRUE)
        ret[is.na(ret)] <- comb[is.na(ret)]
        return(ret)
    }

    dyfield <- function(field, yr){
        t(dxfield(t(field), yr))
    }

    LIV <- lapply(1:nenv,
                  function(i)
                      lapply(1:dim(FIELDS[[i]])[3],
                             function(j)
                                 RTMB::interpol2Dfun(FIELDS[[i]][,,j],
                                                     xlim = XR[i,],
                                                     ylim = YR[i,],
                                                     R = 1)))
    LIVdx <- lapply(1:nenv,
                    function(i)
                        lapply(1:dim(FIELDS[[i]])[3],
                               function(j)
                                   RTMB::interpol2Dfun(dxfield(FIELDS[[i]][,,j],
                                                               XR[i,]),
                                                       xlim = XR[i,],
                                                       ylim = YR[i,],
                                                       R = 1)))
    LIVdy <- lapply(1:nenv,
                    function(i)
                        lapply(1:dim(FIELDS[[i]])[3],
                               function(j)
                                   RTMB::interpol2Dfun(dyfield(FIELDS[[i]][,,j],
                                                               YR[i,]),
                                                       xlim = XR[i,],
                                                       ylim = YR[i,],
                                                       R = 1)))

    val <- function(xy, t){
        "c" <- ADoverload("c")
        "[<-" <- ADoverload("[<-")
        if(!inherits(xy, "matrix") && !inherits(xy, "data.frame")) xy <- RTMB::matrix(xy, 1, 2)
        h <- rep(0, nrow(xy))
        for(i in 1:nenv){
            it <- ienv[i,as.integer(cut(t, time.cont, include.lowest = TRUE))]
            is <- ienvS[i,as.integer(cut(t, time.cont, include.lowest = TRUE))]
            if(!is.empty(S[[i]][[is]]) && it != 0){
                for(c in 1:nrow(xy)){
                    h[c] <- h[c] + S[[i]][[is]](LIV[[i]][[it]](
                        smooth.identity(xy[c,1], XR[i,1], XR[i,2]),
                        smooth.identity(xy[c,2], YR[i,1], YR[i,2])))
                }
            }
        }
        return(h)
    }


    grad <- function(xy, t){
        "c" <- ADoverload("c")
        "[<-" <- ADoverload("[<-")
        if(!inherits(xy, "matrix") && !inherits(xy, "data.frame")) xy <- RTMB::matrix(xy, 1, 2)
        dh <- RTMB::matrix(0, nrow(xy), 2)
        dxytmp <- c(0,0)
        for(i in 1:nenv){
            it <- ienv[i,as.integer(cut(t, time.cont, include.lowest = TRUE))]
            is <- ienvS[i,as.integer(cut(t, time.cont, include.lowest = TRUE))]
            if(!is.empty(dS[[i]][[is]]) && it != 0){
                for(c in 1:nrow(xy)){
                    dxytmp[1] <- LIVdx[[i]][[it]](
                        smooth.identity(xy[c,1], XR[i,1], XR[i,2]),
                        smooth.identity(xy[c,2], YR[i,1], YR[i,2]))
                    dxytmp[2] <- LIVdy[[i]][[it]](
                        smooth.identity(xy[c,1], XR[i,1], XR[i,2]),
                        smooth.identity(xy[c,2], YR[i,1], YR[i,2]))
                    dh[c,] <- dh[c,] + dS[[i]][[is]](LIV[[i]][[it]](
                        smooth.identity(xy[c,1], XR[i,1], XR[i,2]),
                        smooth.identity(xy[c,2], YR[i,1], YR[i,2]))) * dxytmp
                }
            }
        }
        return(dh)
    }

    valF <- function(xy, t){
        h <- 0
        for(i in 1:nenv){
            it <- ienv[i,as.integer(cut(t, time.cont, include.lowest = TRUE))]
            h <- h + LIV[[i]][[it]](
                smooth.identity(xy[1], XR[i,1], XR[i,2]),
                smooth.identity(xy[2], YR[i,1], YR[i,2]))
        }
        return(h)
    }

    env2val <- function(env, combine = FALSE){
        "c" <- ADoverload("c")
        "[<-" <- ADoverload("[<-")
        h <- RTMB::matrix(0, nrow(env), ncol(env))
        ## HERE: env for each season?
        is <- 1
        for(i in 1:ncol(env)){
            if(!is.empty(S[[i]][[is]])){
                h[,i] <- h[,i] + S[[i]][[is]](env[,i])
            }
        }
        if(combine){
            h <- RTMB::apply(h, 1, sum)
        }
        return(h)
    }

    res <- list(
        LIV = LIV,
        LIVdx = LIVdx,
        LIVdy = LIVdy,
        S = S,
        dS = dS,
        val = val,
        grad = grad,
        valF = valF,
        env2val = env2val
    )
    return(res)
}


##' Light habitat function
##'
##' @description Creates a light habitat class that is being used within
##'     \emph{momo}.
##'
##' @param FIELDS A list with 3-D array fields that are used for interpolation,
##'     where the first 2 dimensions span the x and y direction of the spatial
##'     field and the third dimension indicates the time dimension.
##' @param XR Limits of the 2-D fields in x direction.
##' @param YR Limits of the 2-D fields in y direction.
##' @param ienv Indicator matrix mapping each model time step to the time steps
##'     of each field. The first dimension corresponds to the number of fields
##'     and the second to the number of model time steps.
##' @param time.cont Vector with continuous model time steps.
##'
##' @return List with functions.
##'
habi.light <- function(FIELDS, XR, YR, ienv, time.cont){

    nenv <- length(FIELDS)

    LIV <- lapply(1:nenv,
                  function(i)
                      lapply(1:dim(FIELDS[[i]])[3],
                             function(j)
                                 RTMB::interpol2Dfun(FIELDS[[i]][,,1],
                                                     xlim = XR[i,],
                                                     ylim = YR[i,],
                                                     R = 1)))

    val <- function(xy, t){
        "c" <- ADoverload("c")
        "[<-" <- ADoverload("[<-")
        if(!inherits(xy, "matrix") && !inherits(xy, "data.frame")) xy <- RTMB::matrix(xy, 1, 2)
        h <- rep(0, nrow(xy))
        for(i in 1:nenv){
            it <- ienv[i,as.integer(cut(t, time.cont, include.lowest = TRUE))]
            if(it != 0){
                for(c in 1:nrow(xy)){
                    h[c] <- h[c] + LIV[[i]][[it]](
                        smooth.identity(xy[c,1], XR[i,1], XR[i,2]),
                        smooth.identity(xy[c,2], YR[i,1], YR[i,2]))
                }
            }
        }
        return(h)
    }

    res <- list(LIV = LIV,
                val = val)
    return(res)
}
