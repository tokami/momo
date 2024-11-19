##' Negative log likelihood function
##'
##' @param par a list of initial parameters
##' @param dat the data set
##'
##' @details exported, but mainly intended for internal use
##'
##' @return a scalar function value
##'
##' @importFrom Matrix expm
##'
##' @export
nll <- function(par, dat){

    ## R's JIT does not preserve a few specific definitions made by RTMB
    "c" <- RTMB::ADoverload("c")
    "[<-" <- RTMB::ADoverload("[<-")
    "diag<-" <- RTMB::ADoverload("diag<-")

    nll <- loglik.ctags <- loglik.atags <- 0

    ## Flags
    flag.const.dif <- ifelse(nrow(dat$knots.dif) == 1, TRUE, FALSE)
    flag.effort <- dat$use.effort

    ## Conversions
    sdObsATS <- exp(par$logSdObsATS)
    lambda <- exp(par$logLambda)
    nat.mort <- exp(par$logM)

    ## Dimensions
    nr <- nrow(dat$rel.events)
    nenv <- length(dat$env)
    nt <- length(dat$time.cont)
    ntp <- length(dat$time.cont.pred)
    nc <- nrow(dat$igrid)
    ncp <- nrow(dat$igrid.pred)
    ne <- length(dat$effort)
    ncem <- nc + ne + 1


    ## Setup environmental functions --------------------------
    env.func.tax <- env.func.dif <-
        env.func.adv.x <- env.func.adv.y <-
            env.dfunc.tax <- env.dfunc.dif <-
                env.dfunc.adv <- vector("list", nenv)
    for(i in 1:nenv){
        env.func.tax[[i]] <- poly.fun(dat$knots.tax[,i], par$alpha[,i])
        env.dfunc.tax[[i]] <- poly.fun(dat$knots.tax[,i], par$alpha[,i], deriv = TRUE)
        env.func.dif[[i]] <- poly.fun(dat$knots.dif[,i], par$beta[,i])
        env.func.adv.x[[i]] <- poly.fun(NULL, par$gamma[1,i], adv = TRUE)
        env.func.adv.y[[i]] <- poly.fun(NULL, par$gamma[2,i], adv = TRUE)
    }

    ## Setup habitat objects ---------------------------------
    habitat.tax <- habi.full(dat$env, dat$xranges, dat$yranges,
                             dat$ienv$tax, dat$time.cont,
                             env.func.tax, env.dfunc.tax)
    habitat.dif <- habi.full(dat$env, dat$xranges, dat$yranges,
                             dat$ienv$dif, dat$time.cont,
                             env.func.dif, env.dfunc.dif)
    habitat.adv.x <- habi.full(dat$env, dat$xranges, dat$yranges,
                               dat$ienv$adv.x, dat$time.cont,
                               env.func.adv.x, env.dfunc.adv)
    habitat.adv.y <- habi.full(dat$env, dat$xranges, dat$yranges,
                               dat$ienv$adv.y, dat$time.cont,
                               env.func.adv.y, env.dfunc.adv)

    if(flag.effort){
        effort <- habi.light(dat$effort, dat$xranges.eff, dat$yranges.eff,
                             dat$ieff, dat$time.cont)
    }


    if(dat$use.kf){

        ## Kalman filter  ---------------------------------

        ## Conventional tags
        if(dat$use.ctags){

            nctags <- nrow(dat$ctags)
            res.cx <- res.cy <- res.ct <- rep(NA, nctags)

            mxy <- vxy <- c(0,0)

            for(i in 1:nctags){

                mxy[1] <- dat$ctags$x0[i]
                mxy[2] <- dat$ctags$y0[i]
                vxy[1] <- 0
                vxy[2] <- 0

                ts <- seq(dat$ctags$t0[i], dat$ctags$t1[i], by = dat$ddt)
                nts <- length(ts)

                if(nts > 0){

                    lpC <- lpC.cum <- lpS <- rep(0, nts+1)

                    for(t in 1:nts){
                        if(dat$use.taxis){
                            mxy <- mxy + habitat.tax$grad(mxy, ts[t]) * dat$ddt
                        }
                        if(dat$use.advection){
                            mxy[1] <- mxy[1] + habitat.adv.x$val(mxy, ts[t]) * dat$ddt
                            mxy[2] <- mxy[2] + habitat.adv.y$val(mxy, ts[t]) * dat$ddt
                        }
                        vxy <- vxy + 2 * exp(habitat.dif$val(mxy, ts[t])) * dat$ddt

                        if(flag.effort){
                            e <- 1 ## fleet LATER:
                            fish.mort <- (lambda[1,e] + par$lambdaEC[1,e] * ts[t]) *
                                effort$val(mxy, ts[t])
                            tot.mort <- fish.mort + nat.mort
                            lpS[t+1] <- lpS[t] - tot.mort * dat$ddt
                            lpC[t+1] <- lpS[t] + log(fish.mort) - log(tot.mort) +
                                log1p(-exp(-tot.mort * dat$ddt))
                            lpC.cum[t+1] <- RTMBconvenience::logspace_add(lpC.cum[t],
                                                                          lpC[t+1])
                        }
                    }
                    if(flag.effort){
                        lpC.cum <- RTMBconvenience::logspace_sub(lpC.cum, 0)
                    }

                    if(dat$excl.ctags[i] == 0){

                        if(!is.na(dat$ctags$x1[i])){

                            loglik.ctags <- loglik.ctags +
                                dnorm(dat$ctags$x1[i], mxy[1], sqrt(vxy[1]), TRUE)
                            res.cx[i] <- (dat$ctags$x1[i] - mxy[1])/sqrt(vxy[1])
                            loglik.ctags <- loglik.ctags +
                                dnorm(dat$ctags$y1[i], mxy[2], sqrt(vxy[2]), TRUE)
                            res.cy[i] <- (dat$ctags$y1[i] - mxy[2])/sqrt(vxy[2])

                            if(flag.effort){
                                loglik.ctags <- loglik.ctags + lpC[t+1]
                                if(!RTMB:::ad_context()){
                                    res.ct[i] <- qnorm(runif(1, exp(lpC.cum[t]),
                                                             exp(lpC.cum[t+1])))
                                }
                            }
                        }else if(flag.effort){
                            loglik.ctags <- loglik.ctags + log1p(-exp(lpC.cum[t+1]))
                            if(!RTMB:::ad_context()){
                                res.ct[i] <- qnorm(runif(1, exp(lpC.cum[t+1]), 1))
                            }
                        }
                    }
                }
            }

            nll <- nll - loglik.ctags

            REPORT(res.cx)
            REPORT(res.cy)
            REPORT(res.ct)
        }

        ## Archival tags
        if(dat$use.atags){

            natags <- length(dat$atags)

            res.axy <- vector("list", natags)
            res.at <- rep(NA, natags)

            for(at in 1:natags){

                thisAT <- dat$atags[[at]]

                res.axy[[at]] <- matrix(NA, nrow(thisAT), 2)

                ## Use tag?
                if(nrow(thisAT) > 3 && all(thisAT$use == 1)){

                    obsvar <- exp(2 * par$logSdObsATS)
                    lastt <- as.numeric(thisAT[1,1])
                    lastxy <- as.numeric(thisAT[1,2:3])
                    P <- D <- A <- c(0,0)

                    lpS <- lpC <- lpC.cum <- rep(0, nrow(thisAT))

                    for(r in 2:nrow(thisAT)){
                        thist <- thisAT[r,1]
                        dt <- (thist - lastt)
                        D <- exp(habitat.dif$val(lastxy, thisAT[r,1]))
                        if(dat$use.taxis){
                            T <- habitat.tax$grad(lastxy, thisAT[r,1])
                        }
                        if(dat$use.advection){
                            A[1] <- habitat.adv.x$val(lastxy, thisAT[r,1])
                            A[2] <- habitat.adv.y$val(lastxy, thisAT[r,1])
                        }
                        predxy <- lastxy + (A + T) * dt
                        PP <- P + 2 * D * dt
                        if(r == nrow(thisAT)){
                            F <- PP
                        }else{
                            F <- PP + obsvar
                        }
                        thisxy <- as.numeric(thisAT[r,2:3])
                        w <- thisxy - predxy
                        lastxy <- predxy + PP / F * w
                        P <- PP - PP / F * PP
                        lastt <- thist
                        loglik.atags <- loglik.atags + dnorm(w[1], 0, sqrt(F[1]), TRUE)
                        loglik.atags <- loglik.atags + dnorm(w[2], 0, sqrt(F[2]), TRUE)
                        res.axy[[at]][r,1] <- w[1] / sqrt(F[1])
                        res.axy[[at]][r,2] <- w[2] / sqrt(F[2])

                        if(flag.effort){
                            e <- 1 ## fleet LATER:
                            fish.mort <- (lambda[1,e] + par$lambdaEC[1,e] * thist) * effort$val(predxy, thist)
                            tot.mort <- fish.mort + nat.mort
                            lpS[r] <- lpS[r-1] - tot.mort * dt
                            lpC[r] <- lpS[r-1] + log(fish.mort) - log(tot.mort) +
                                log1p(-exp(-tot.mort * dt))
                            lpC.cum[r] <- RTMBconvenience::logspace_add(lpC.cum[r-1],
                                                                        lpC[r])
                        } ## end effort
                    }
                    if(flag.effort){
                        loglik.atags <- loglik.atags + lpC[r]
                        if(!RTMB:::ad_context()){
                            res.at[at] <- qnorm(runif(1,
                                                      exp(lpC.cum[r-1]),
                                                      exp(lpC.cum[r-1]) +
                                                      fish.mort / tot.mort *
                                                      (1 - exp(-tot.mort * dt))))
                        }
                    }
                }
            }

            nll <- nll - loglik.atags

            REPORT(res.axy)
            REPORT(res.at)
        }

    }else{

        ## Matrix exponential ---------------------------------

        ## Conventional tags
        if(dat$use.ctags){

            nctags <- nrow(dat$ctags)
            res.cx <- res.cy <- res.ct <- rep(NA, nctags)

            ## Loop over release events
            for(r in 1:nr){

                itrel <- dat$rel.events$itrel[r]
                icrel <- dat$rel.events$icrel[r]
                itrec <- dat$rel.events$itrec[r]
                itmax <- itrec - itrel

                if(itmax > 0){

                    dt.ctags <- (dat$time.cont[itrec] - dat$time.cont[itrel]) / itmax
                    itmaxp1 <- itmax + 1

                    ## Distribution probability
                    dist.prob <- RTMB::matrix(0, itmaxp1, ncem)
                    dist.prob[1, icrel] <- 1


                    ## Loop over time
                    for(t in 2:itmaxp1){

                        itabs <- itrel + t - 1

                        ## Set to zero
                        Zstar <- Dstar <- Astar <-
                            Cstar <- Ostar <- RTMB::matrix(0, ncem, ncem)

                        ## Taxis rate -----------------------
                        tmp <- habitat.tax$grad(dat$xygrid, dat$time.cont[itabs])
                        j = 2
                        ind <- which(!is.na(dat$nextTo[,j]))
                        Zstar[cbind(ind, dat$nextTo[ind,j])] <-
                            0.5 * tmp[ind,1] / dat$next.dist[j-1]
                        j = 3
                        ind <- which(!is.na(dat$nextTo[,j]))
                        Zstar[cbind(ind, dat$nextTo[ind,j])] <-
                            -0.5 * tmp[ind,1] / dat$next.dist[j-1]
                        j = 4
                        ind <- which(!is.na(dat$nextTo[,j]))
                        Zstar[cbind(ind, dat$nextTo[ind,j])] <-
                            0.5 * tmp[ind,2] / dat$next.dist[j-1]
                        j = 5
                        ind <- which(!is.na(dat$nextTo[,j]))
                        Zstar[cbind(ind, dat$nextTo[ind,j])] <-
                            - 0.5 * tmp[ind,2] / dat$next.dist[j-1]

                        ## Diffusion rate --------------------
                        hD <- habitat.dif$val(dat$xygrid, dat$time.cont[itabs])
                        for(j in 2:ncol(dat$nextTo)){
                            ind <- which(!is.na(dat$nextTo[,j]))
                            Dstar[cbind(ind, dat$nextTo[ind,j])] <-
                                1.5 * exp(hD[ind]) / dat$next.dist[j-1]^2
                        }

                        ## Advection rate -------------------------
                        hAx <- habitat.adv.x$val(dat$xygrid, dat$time.cont[itabs])
                        hAy <- habitat.adv.y$val(dat$xygrid, dat$time.cont[itabs])
                        ## right neighbour
                        j = 2
                        ind <- which(!is.na(dat$nextTo[,j]))
                        Astar[cbind(ind, dat$nextTo[ind,j])] <-
                            0.5 * hAx[ind] / dat$next.dist[j-1]
                        ## left neighbour
                        j = 3
                        ind <- which(!is.na(dat$nextTo[,j]))
                        Astar[cbind(ind, dat$nextTo[ind,j])] <-
                            - 0.5 * hAx[ind] / dat$next.dist[j-1]
                        ## top neighbour
                        j = 4
                        ind <- which(!is.na(dat$nextTo[,j]))
                        Astar[cbind(ind, dat$nextTo[ind,j])] <-
                            0.5 * hAy[ind] / dat$next.dist[j-1]
                        ## down neighbour
                        j = 5
                        ind <- which(!is.na(dat$nextTo[,j]))
                        Astar[cbind(ind, dat$nextTo[ind,j])] <-
                            - 0.5 * hAy[ind] / dat$next.dist[j-1]

                        ## Mass balance ---------------------------
                        Zstar[cbind(1:ncem,1:ncem)] <- - RTMB::rowSums(Zstar)
                        Dstar[cbind(1:ncem,1:ncem)] <- - RTMB::rowSums(Dstar)

                        ## Movement rates ------------------------
                        Mstar <- Zstar + Dstar + Astar

                        ## Add natural and fishing mortality
                        if(flag.effort){

                            Cstar <- Ostar <- RTMB::matrix(0, ncem, ncem)

                            ## Instantaneous fishing mortality rate
                            for(e in 1:ne){
                                fish.mort <- (lambda[1,e] + par$lambdaEC[1,e] * dat$time.cont[itabs]) *
                                    effort$val(dat$xygrid, dat$time.cont[itabs]) ## how does this work with multiple fleets?
                                Cstar[1:nc,nc+e] <- fish.mort
                            }

                            ## Instantaneous natural mortality rate
                            Ostar[1:nc,ncem] <- rep(nat.mort, nc)

                            ## Mass balance
                            Cstar[cbind(1:ncem,1:ncem)] <- - RTMB::rowSums(Cstar)
                            Ostar[cbind(1:ncem,1:ncem)] <- - RTMB::rowSums(Ostar)

                            ## Movement rates
                            Mstar <- Mstar + Cstar + Ostar
                        }

                        ## Account for time step
                        Mstar <- Mstar * dt.ctags

                        ## Matrix exponential
                        M <- matexpo(Mstar, dat$log2steps)

                        ## Distribution probability
                        dist.prob[t,] <- RTMB:::matrix(RTMB:::matrix(dist.prob[t-1,], 1, ncem) %*% M, 1, ncem)

                    } ## end of time loop

                    ## Likelihood contribution
                    thisCT <- dat$ctags[dat$irel.event == r,]
                    nh <- nrow(thisCT)
                    if(flag.effort){
                        lpC <- rep(0, nh)
                    }

                    for(h in 1:nh){
                        ind.tag <- which(dat$irel.event == r)[h]

                        if(dat$excl.ctags[ind.tag] == 0){
                            e <- 1 ## fleet LATER:
                            ## eind <- thisCT[h,"fleet"]

                            if(flag.effort){
                                lpC.cum <- log(dist.prob[itmaxp1, nc + 1])
                                if(ne > 1){
                                    for(e in 2:ne){
                                        lpC.cum <- RTMBconvenience::logspace_add(lpC.cum, dist.prob[itmaxp1, nc + e])
                                    }
                                }
                            }

                            if(dat$rec[ind.tag]){

                                ## Recovered tags
                                recapTime2 <- dat$itrec[ind.tag]
                                recapTime <- as.integer(recapTime2 - itrel)
                                recapLoc <- as.integer(dat$icrec[ind.tag])
                                ## print(recapTime+1)
                                ## print(dim(dist.prob))
                                loglik.ctags <- loglik.ctags +
                                    log(dist.prob[recapTime+1, recapLoc])

                                if(flag.effort){
                                    fish.mort <- (lambda[1,e] + par$lambdaEC[1,e] * dat$time.cont[recapTime2]) *
                                        effort$val(dat$xygrid[recapLoc,], dat$time.cont[recapTime2]) ## how does this work with multiple fleets?
                                    tot.mort <- fish.mort + nat.mort
                                    lpC[h] <- log(fish.mort) - log(tot.mort) + log1p(-exp(-tot.mort * dt.ctags))
                                    loglik.ctags <- loglik.ctags + lpC[h]
                                }

                                ## Residuals ----------------------------
                                if(!RTMB:::ad_context()){

                                    n <- max(dat$igrid$idx)
                                    xres <- rep(0, n)
                                    for(i in 1:n){
                                        xres[i] <- sum(dist.prob[recapTime+1,1:nc][dat$igrid$idx == i])
                                    }
                                    xres.cum <- c(0, cumsum(xres))
                                    xres.cum <- xres.cum / xres.cum[length(xres.cum)]
                                    res.cx[ind.tag] <- qnorm(runif(1,
                                                                   xres.cum[dat$igrid$idx[recapLoc]],
                                                                   xres.cum[dat$igrid$idx[recapLoc]+1]))
                                    yres.cum <- c(0, cumsum(dist.prob[recapTime+1, dat$igrid$idx == dat$igrid$idx[recapLoc]]))
                                    yres.cum <- yres.cum / yres.cum[length(yres.cum)]
                                    res.cy[ind.tag] <- qnorm(runif(1,
                                                                   yres.cum[dat$igrid$idy[recapLoc]],
                                                                   yres.cum[dat$igrid$idy[recapLoc]+1]))

                                    if(flag.effort){
                                        res.ct[ind.tag] <- qnorm(runif(1,
                                                                       exp(lpC.cum),
                                                                       exp(lpC.cum) + exp(lpC[h])))
                                        ## TODO: that should be the probability of t-1 and t
                                    }
                                }

                            }else if(flag.effort){

                                ## Not recaptured tags: prob of being in any cell + prob of M = (1 - prob being caught)
                                lpNoC <- RTMBconvenience::logspace_sub(0, lpC.cum)
                                ## prob still in one of the cells or dead due to nat mort
                                loglik.ctags <- loglik.ctags + lpNoC

                                ## Only time residual for non-recovered tags
                                if(!RTMB:::ad_context()){
                                    res.ct[ind.tag] <- qnorm(runif(1, exp(lpC.cum), 1))
                                }

                            }
                        }
                    } ## end tag loop

                } ## end if itmax == 0

            } ## end release location loop



            nll <- nll - loglik.ctags

            REPORT(res.cx)
            REPORT(res.cy)
            REPORT(res.ct)

        }

        if(dat$use.atags){

            natags <- length(dat$atags)

            res.axy <- vector("list", natags)
            res.at <- rep(NA, natags)

            ## Loop over atags
            for(a in 1:natags){

                tag <- dat$atags[[a]]
                itrel <- tag$it[1]
                icrel <- tag$ic[1]
                itmax <- nrow(tag)
                itrec <- tag$it[itmax]

                res.axy[[a]] <- matrix(NA, itmax, 2)

                ## Use tag?
                if(itmax > 3 && all(tag$use == 1)){

                    ## Distribution probability
                    dist.prob <- RTMB::matrix(0, itmax, ncem)
                    dist.prob[1, icrel] <- 1

                    ## Loop over time
                    for(t in 2:itmax){

                        itabs <- tag$it[t]
                        dt.atags <- tag$t[t] - tag$t[t-1]

                        ## Set to zero
                        Zstar <- Dstar <- Astar <-
                            Cstar <- Ostar <- RTMB::matrix(0, ncem, ncem)

                        ## Taxis rate ----------------------
                        tmp <- habitat.tax$grad(dat$xygrid, dat$time.cont[itabs])
                        j = 2
                        ind <- which(!is.na(dat$nextTo[,j]))
                        Zstar[cbind(ind, dat$nextTo[ind,j])] <-
                            0.5 * tmp[ind,1] / dat$next.dist[j-1]
                        j = 3
                        ind <- which(!is.na(dat$nextTo[,j]))
                        Zstar[cbind(ind, dat$nextTo[ind,j])] <-
                            -0.5 * tmp[ind,1] / dat$next.dist[j-1]
                        j = 4
                        ind <- which(!is.na(dat$nextTo[,j]))
                        Zstar[cbind(ind, dat$nextTo[ind,j])] <-
                            0.5 * tmp[ind,2] / dat$next.dist[j-1]
                        j = 5
                        ind <- which(!is.na(dat$nextTo[,j]))
                        Zstar[cbind(ind, dat$nextTo[ind,j])] <-
                            - 0.5 * tmp[ind,2] / dat$next.dist[j-1]

                        ## Diffusion rate ---------------------
                        hD <- habitat.dif$val(dat$xygrid, dat$time.cont[itabs])
                        for(j in 2:ncol(dat$nextTo)){
                            ind <- which(!is.na(dat$nextTo[,j]))
                            Dstar[cbind(ind, dat$nextTo[ind,j])] <-
                                1.5 * exp(hD[ind]) / dat$next.dist[j-1]^2
                        }

                        ## Advection rate ---------------------
                        hAx <- habitat.adv.x$val(dat$xygrid, dat$time.cont[itabs])
                        hAy <- habitat.adv.y$val(dat$xygrid, dat$time.cont[itabs])
                        ## right neighbour
                        j = 2
                        ind <- which(!is.na(dat$nextTo[,j]))
                        Astar[cbind(ind, dat$nextTo[ind,j])] <-
                            0.5 * hAx[ind] / dat$next.dist[j-1]
                        ## left neighbour
                        j = 3
                        ind <- which(!is.na(dat$nextTo[,j]))
                        Astar[cbind(ind, dat$nextTo[ind,j])] <-
                            -0.5 * hAx[ind] / dat$next.dist[j-1]
                        ## top neighbour
                        j = 4
                        ind <- which(!is.na(dat$nextTo[,j]))
                        Astar[cbind(ind, dat$nextTo[ind,j])] <-
                            0.5 * hAy[ind] / dat$next.dist[j-1]
                        ## down neighbour
                        j = 5
                        ind <- which(!is.na(dat$nextTo[,j]))
                        Astar[cbind(ind, dat$nextTo[ind,j])] <-
                            -0.5 * hAy[ind] / dat$next.dist[j-1]

                        ## Mass balance ----------------------
                        Zstar[cbind(1:ncem,1:ncem)] <- - RTMB::rowSums(Zstar)
                        Dstar[cbind(1:ncem,1:ncem)] <- - RTMB::rowSums(Dstar)

                        ## Movement rates ------------------
                        Mstar <- Zstar + Dstar + Astar

                        ## Add natural and fishing mortality
                        if(flag.effort){

                            Cstar <- Ostar <- RTMB::matrix(0, ncem, ncem)

                            ## Instantaneous fishing mortality rate
                            for(e in 1:ne){
                                fish.mort <- (lambda[1,e] + par$lambdaEC[1,e] * dat$time.cont[itabs]) *
                                    effort$val(dat$xygrid, dat$time.cont[itabs]) ## how does this work with multiple fleets?
                                Cstar[1:nc,nc+e] <- fish.mort
                            }

                            ## Instantaneous natural mortality rate
                            Ostar[1:nc,ncem] <- rep(nat.mort, nc)

                            ## Mass balance
                            Cstar[cbind(1:ncem,1:ncem)] <- - RTMB::rowSums(Cstar)
                            Ostar[cbind(1:ncem,1:ncem)] <- - RTMB::rowSums(Ostar)

                            ## Movement rates
                            Mstar <- Mstar + Cstar + Ostar
                        }

                        ## Account for time step
                        Mstar <- Mstar * dt.atags

                        ## Matrix exponential
                        M <- matexpo(Mstar, dat$log2steps)

                        ## Distribution probability
                        dist.prob[t,] <- RTMB:::matrix(RTMB:::matrix(dist.prob[t-1,], 1, ncem) %*% M, 1, ncem)

                    } ## end of time loop


                    ## Likelihood contribution -----------------------------
                    for(t in 2:(itmax-1)){
                        ## px <- dnorm(dat$xygrid[,1], tag$x[t], sdObsATS)
                        ## px <- px/sum(px)
                        ## py <- dnorm(dat$xygrid[,2], tag$y[t], sdObsATS)
                        ## py <- py/sum(py)
                        ## tmp <- log(sum(dist.prob[t,1:nc] * px * py))
                        ## tmp <- log(sum(dist.prob[t,1:nc] *
                        ##                dnorm(dat$xygrid[,1], tag$x[t], sdObsATS) *
                        ##                dnorm(dat$xygrid[,2], tag$y[t], sdObsATS)))
                        ## tmp <- log(
                        ##     sum(dist.prob[t,1:nc] *
                        ##         (pnorm(dat$xgr[dat$igrid$idx+1], tag$x[t], sdObsATS) -
                        ##          pnorm(dat$xgr[dat$igrid$idx], tag$x[t], sdObsATS)) *
                        ##         (pnorm(dat$ygr[dat$igrid$idy+1], tag$y[t], sdObsATS) -
                        ##          pnorm(dat$ygr[dat$igrid$idy], tag$y[t], sdObsATS))))
                        tmp <- log(dist.prob[t, tag$ic[t]])
                        loglik.atags <- loglik.atags + tmp
                    }
                    loglik.atags <- loglik.atags +
                        log(dist.prob[itmax,tag$ic[itmax]])


                    if(!RTMB:::ad_context()){

                        n <- max(dat$igrid$idx)

                        for(t in 2:itmax){
                            recapLoc <- tag$ic[t]
                            xres <- rep(0, n)
                            for(i in 1:n){
                                xres[i] <- sum(dist.prob[t,1:nc][dat$igrid$idx == i])
                            }
                            xres.cum <- c(0, cumsum(xres))
                            xres.cum <- xres.cum / xres.cum[length(xres.cum)]
                            res.axy[[a]][t,1] <- qnorm(runif(1,
                                                             xres.cum[dat$igrid$idx[recapLoc]],
                                                             xres.cum[dat$igrid$idx[recapLoc] + 1]))
                            yres.cum <- c(0,
                                          cumsum(dist.prob[t,
                                                           dat$igrid$idx ==
                                                           dat$igrid$idx[recapLoc]]))
                            yres.cum <- yres.cum / yres.cum[length(yres.cum)]
                            res.axy[[a]][t,2] <- qnorm(runif(1,
                                                             yres.cum[dat$igrid$idy[recapLoc]],
                                                             yres.cum[dat$igrid$idy[recapLoc] + 1]))
                        }
                    }

                    if(flag.effort){
                        e <- 1 ## fleet LATER:
                        fish.mort <- (lambda[1,e] + par$lambdaEC[1,e] * dat$time.cont[itabs]) *
                            effort$val(dat$xygrid[recapLoc,], dat$time.cont[itabs]) ## how does this work with multiple fleets?
                        tot.mort <- fish.mort + nat.mort
                        lpC <- log(fish.mort) - log(tot.mort) + log1p(-exp(-tot.mort * dt.atags))
                        loglik.atags <- loglik.atags + lpC

                        if(!RTMB:::ad_context()){
                            lpC.cum <- log(dist.prob[itmax, nc + 1])
                            if(ne > 1){
                                for(e in 2:ne){
                                    lpC.cum <- RTMBconvenience::logspace_add(lpC.cum, dist.prob[itmax, nc + e])
                                }
                            }
                            res.at[a] <- runif(1,
                                               exp(lpC.cum),
                                               exp(lpC.cum) + exp(lpC))
                        }
                    }

                } ## use this a-tag
            } ## end of atag loop

            nll <- nll - loglik.atags

            REPORT(res.axy)
            REPORT(res.at)

        } ## use.atags

    } ## expm

    REPORT(loglik.ctags)
    REPORT(loglik.atags)

    ## Prediction ---------------------------------
    prefT.pred <- habitat.tax$env2val(dat$env.pred)
    prefD.pred <- habitat.dif$env2val(dat$env.pred)

    REPORT(prefT.pred)
    ADREPORT(prefT.pred)
    REPORT(prefD.pred)
    ADREPORT(prefD.pred)


    hT.pred <- hTdx.pred <- hTdy.pred <- hD.pred <-
        hAx.pred <- hAy.pred <- matrix(0, ncp, ntp)
    for(t in 1:ntp){
        hT.pred[,t] <- habitat.tax$val(dat$xygrid.pred, dat$time.cont.pred[t])
        tmp <- habitat.tax$grad(dat$xygrid.pred, dat$time.cont.pred[t])
        hTdx.pred[,t] <- tmp[,1]
        hTdy.pred[,t] <- tmp[,2]
        hD.pred[,t] <- habitat.dif$val(dat$xygrid.pred, dat$time.cont.pred[t])
        hAx.pred[,t] <- habitat.adv.x$val(dat$xygrid.pred, dat$time.cont.pred[t])
        hAy.pred[,t] <- habitat.adv.y$val(dat$xygrid.pred, dat$time.cont.pred[t])
    }

    REPORT(hT.pred)
    REPORT(hTdx.pred)
    REPORT(hTdy.pred)
    REPORT(hD.pred)
    REPORT(hAx.pred)
    REPORT(hAy.pred)

    return(nll)
}
