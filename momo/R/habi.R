


get.liv <- function(FIELDS, XR, YR) {

    nenv <- length(FIELDS)

    dxfield <- function(field, xr) {
        nr <- nrow(field)
        nc <- ncol(field)
        ret <- matrix(NA_real_, nr, nc)

        onedx <- (xr[2] - xr[1]) / (nr - 1) ## as xr cell centers

        for (i in 1:nr) {
            for (j in 1:nc) {
                val <- field[i, j]
                if (is.na(val)) next

                left <- if (i > 1) field[i - 1, j] else NA_real_
                right <- if (i < nr) field[i + 1, j] else NA_real_

                if (!is.na(left) && !is.na(right)) {
                    ret[i, j] <- (right - left) / (2 * onedx)
                } else if (!is.na(right)) {
                    ret[i, j] <- (right - val) / onedx
                } else if (!is.na(left)) {
                    ret[i, j] <- (val - left) / onedx
                }
            }
        }

       ret[is.na(ret)] <- 0

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

    LIV0 <- lapply(1:nenv,
                  function(i)
                      lapply(1:dim(FIELDS[[i]])[3],
                             function(j){
                                 fieldi <- FIELDS[[i]][,,j]
                                 fieldi[is.na(fieldi)] <- 0
                                 RTMB::interpol2Dfun(fieldi,
                                                     xlim = XR[i,],
                                                     ylim = YR[i,],
                                                     R = 1)}))

    res <- list(LIV = LIV,
                LIV0 = LIV0,
                LIVdx = LIVdx,
                LIVdy = LIVdy)

    return(res)
}



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
habi.full <- function(liv, XR, YR, ienv, time.cont, S, dS, ienvS){
    "c" <- ADoverload("c")
    "[<-" <- ADoverload("[<-")

    LIV <- liv$LIV
    LIV0 <- liv$LIV0
    LIVdx <- liv$LIVdx
    LIVdy <- liv$LIVdy

    nenv <- length(LIV)

    val <- function(xy, t){
        t.index <- t2index(t, time.cont)
        h <- rep(0, nrow(xy))
        for(i in 1:nenv){
            it <- ienv[i,t.index]
            is <- ienvS[i,t.index]
            if(!is.empty(S[[i]]) && !is.empty(S[[i]][[is]]) && it != 0){
                h <- h + S[[i]][[is]](LIV0[[i]][[it]](xy[,1], xy[,2]))  ## HERE: LIV
                    ## smooth.identity(xy[,1], XR[i,1], XR[i,2]),
                    ## smooth.identity(xy[,2], YR[i,1], YR[i,2])))
            }
        }
        return(h)
    }

    grad <- function(xy, t){
        t.index <- t2index(t, time.cont)
        dh <- dxytmp <- RTMB::matrix(0, nrow(xy), 2)
        for(i in 1:nenv){
            it <- ienv[i,t.index]
            is <- ienvS[i,t.index]
            if(!is.empty(dS[[i]]) && !is.empty(dS[[i]][[is]]) && it != 0){
                dxytmp[,1] <- LIVdx[[i]][[it]](xy[,1], xy[,2])
                ## smooth.identity(xy[,1], XR[i,1], XR[i,2]),
                ## smooth.identity(xy[,2], YR[i,1], YR[i,2]))
                dxytmp[,2] <- LIVdy[[i]][[it]](xy[,1], xy[,2])
                ## smooth.identity(xy[,1], XR[i,1], XR[i,2]),
                ## smooth.identity(xy[,2], YR[i,1], YR[i,2]))
                dh <- dh + dS[[i]][[is]](LIV0[[i]][[it]](xy[,1], xy[,2])) * dxytmp  ## HERE: LIV
                ## smooth.identity(xy[,1], XR[i,1], XR[i,2]),
                ## smooth.identity(xy[,2], YR[i,1], YR[i,2]))) * dxytmp
            }
        }
        return(dh)
    }

    valF <- function(xy, t){
        t.index <- t2index(t, time.cont)
        h <- 0
        for(i in 1:nenv){
            it <- ienv[i,t.index]
            h <- h + LIV[[i]][[it]](xy[,1], xy[,2])
                ## smooth.identity(xy[,1], XR[i,1], XR[i,2]),
                ## smooth.identity(xy[,2], YR[i,1], YR[i,2]))
        }
        return(h)
    }

    env2val <- function(env, combine = FALSE){
        h <- matrix(0, nrow(env), ncol(env))
        for(i in 1:ncol(env)){
            is <- ienvS[i,1] ## HERE: find better solution for seasonal splines
            if(!is.empty(S[[i]]) && !is.empty(S[[i]][[is]])){
                h[,i] <- h[,i] + S[[i]][[is]](env[,i])
            }
        }
        if(combine){
            h <- apply(h, 1, sum)
        }
        return(h)
    }

    res <- list(
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
habi.light <- function(liv, XR, YR, ienv, time.cont){
    "c" <- ADoverload("c")
    "[<-" <- ADoverload("[<-")

    LIV <- liv$LIV

    nenv <- length(LIV)

    val <- function(xy, t){
        t.index <- t2index(t, time.cont)
        h <- rep(0, nrow(xy))
        for(i in 1:nenv){
            it <- ienv[i,t.index]
            if(it != 0){
                h <- h + LIV[[i]][[it]](xy[,1], xy[,2])
                    ## smooth.identity(xy[,1], XR[i,1], XR[i,2]),
                    ## smooth.identity(xy[,2], YR[i,1], YR[i,2]))
            }
        }
        return(h)
    }

    res <- list(val = val)
    return(res)
}
