
habi.full <- function(FIELDS, XR, YR, ienv, time.cont, S, dS){

    nenv <- length(FIELDS)

    dxfield <- function(field, xr){
        ret <- matrix(0, nrow(field), ncol(field))
        onedx <- (xr[2] - xr[1]) / nrow(field) ## WHY: (nrow(field) - 1)
        twodx <- 2 * onedx
        for(i in 2:(nrow(field)-1)){
            ret[i,] <- (field[i+1,] - field[i-1,]) / twodx
        }
        ## TODO: make this an option!
        ## ## right edge
        ## i <- nrow(field)
        ## ret[i,] <- (field[i,] - field[i-1,]) / onedx
        ## ## left edge
        ## i <- 1
        ## ret[i,] <- (field[i+1,] - field[i,]) / onedx
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
            if(!is.empty(S[[i]]) && it != 0){
                for(c in 1:nrow(xy)){
                    h[c] <- h[c] + S[[i]](LIV[[i]][[it]](
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
            if(!is.empty(dS[[i]]) && it != 0){
                for(c in 1:nrow(xy)){
                    dxytmp[1] <- LIVdx[[i]][[it]](
                        smooth.identity(xy[c,1], XR[i,1], XR[i,2]),
                        smooth.identity(xy[c,2], YR[i,1], YR[i,2]))
                    dxytmp[2] <- LIVdy[[i]][[it]](
                        smooth.identity(xy[c,1], XR[i,1], XR[i,2]),
                        smooth.identity(xy[c,2], YR[i,1], YR[i,2]))
                    dh[c,] <- dh[c,] + dS[[i]](LIV[[i]][[it]](
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

    env2val <- function(env){
        "c" <- ADoverload("c")
        "[<-" <- ADoverload("[<-")
        h <- RTMB::matrix(0, nrow(env), ncol(env))
        for(i in 1:ncol(env)){
            if(!is.empty(S[[i]])){
                h <- h + S[[i]](env[,i])
            }
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
