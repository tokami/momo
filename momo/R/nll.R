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

    if(dat$use.boundaries){
    bound <- habi.light(dat$boundaries,
                             dat$xranges + c(-1,1) * dat$dxdy[1],
                             dat$yranges + c(-1,1) * dat$dxdy[2],
                        dat$ibound, dat$time.cont)
    }


    if(dat$use.kf){

        ## Kalman filter  ---------------------------------

        ## Conventional tags
        if(dat$use.ctags){

            nctags <- nrow(dat$ctags)
            resid.ctags <- matrix(NA, nctags, 4)

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
                        move <- c(0,0)
                        if(dat$use.taxis){
                            move <- move + habitat.tax$grad(mxy, ts[t]) * dat$ddt
                        }
                        if(dat$use.advection){
                            move <- move + c(habitat.adv.x$val(mxy, ts[t]),
                                             habitat.adv.y$val(mxy, ts[t])) * dat$ddt
                        }
                        if(dat$use.boundaries){
                            move <- move * bound$val(mxy, ts[t])
                        }
                        dif <- exp(habitat.dif$val(mxy, ts[t]))
                        if(dat$use.boundaries){
                            dif <- dif * bound$val(mxy, ts[t])
                        }
                        mxy <- mxy + move
                        vxy <- vxy + 2 * dif * dat$ddt

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

                            ## Likelihood
                            loglik.ctags <- loglik.ctags +
                                dnorm(dat$ctags$x1[i], mxy[1], sqrt(vxy[1]), TRUE)
                            loglik.ctags <- loglik.ctags +
                                dnorm(dat$ctags$y1[i], mxy[2], sqrt(vxy[2]), TRUE)
                            if(flag.effort){
                                loglik.ctags <- loglik.ctags + lpC[t+1]
                            }

                            ## Residuals
                            if(!RTMB:::ad_context()){

                                resid.ctags[i,1] <- (dat$ctags$x1[i] - mxy[1]) /
                                    sqrt(vxy[1])
                                resid.ctags[i,2] <- (dat$ctags$y1[i] - mxy[2]) /
                                    sqrt(vxy[2])
                                resid.ctags[i,3] <- sqrt((dat$ctags$x1[i] - mxy[1])^2 +
                                                         (dat$ctags$y1[i] - mxy[2])^2) /
                                    sqrt((vxy[1] + vxy[2])/2)

                                if(flag.effort){
                                    resid.ctags[i,4] <- qnorm(runif(1, exp(lpC.cum[t]),
                                                                    exp(lpC.cum[t+1])))
                                }

                            }



                        }else if(flag.effort){

                            ## Likelihood
                            loglik.ctags <- loglik.ctags + log1p(-exp(lpC.cum[t+1]))

                            ## Residuals
                            if(!RTMB:::ad_context()){
                                resid.ctags[i,4] <- qnorm(runif(1,
                                                                exp(lpC.cum[t+1]), 1))
                            }
                        }
                    }
                }
            }

            nll <- nll - loglik.ctags

            REPORT(resid.ctags)
        }

        ## Archival tags
        if(dat$use.atags){

            natags <- length(dat$atags)

            resid.atags.fine <- vector("list", natags)
            resid.atags <- matrix(NA, natags, 4)

            for(a in 1:natags){

                tag <- dat$atags[[a]]

                resid.atags.fine[[a]] <- matrix(NA, nrow(tag), 3)

                ## Use tag?
                if(nrow(tag) > 3 && all(tag$use == 1)){

                    obsvar <- exp(2 * par$logSdObsATS)
                    lastt <- as.numeric(tag[1,1])
                    lastxy <- as.numeric(tag[1,2:3])
                    P <- D <- A <- c(0,0)

                    lpS <- lpC <- lpC.cum <- rep(0, nrow(tag))

                    for(r in 2:nrow(tag)){
                        move <- c(0,0)
                        thist <- tag[r,1]
                        dt <- (thist - lastt)
                        if(dat$use.taxis){
                            move <- move + habitat.tax$grad(lastxy, tag[r,1]) * dt
                        }
                        if(dat$use.advection){
                            move <- move + c(habitat.adv.x$val(lastxy, tag[r,1]),
                                             habitat.adv.y$val(lastxy, tag[r,1])) * dt
                        }
                        if(dat$use.boundaries){
                            move <- move * bound$val(lastxy, tag[r,1])
                        }
                        dif <- exp(habitat.dif$val(lastxy, tag[r,1]))
                        if(dat$use.boundaries){
                            dif <- dif * bound$val(lastxy, tag[r,1])
                        }
                        predxy <- lastxy + move
                        PP <- P + 2 * dif * dt
                        if(r == nrow(tag)){
                            F <- PP
                        }else{
                            F <- PP + obsvar
                        }
                        thisxy <- as.numeric(tag[r,2:3])
                        w <- thisxy - predxy
                        lastxy <- predxy + PP / F * w
                        P <- PP - PP / F * PP
                        lastt <- thist
                        loglik.atags <- loglik.atags + dnorm(w[1], 0, sqrt(F[1]), TRUE)
                        loglik.atags <- loglik.atags + dnorm(w[2], 0, sqrt(F[2]), TRUE)


                        if(!RTMB:::ad_context()){
                            resid.atags.fine[[a]][r,1] <- w[1] / sqrt(F[1])
                            resid.atags.fine[[a]][r,2] <- w[2] / sqrt(F[2])
                            resid.atags.fine[[a]][r,3] <- sqrt((thisxy[1] -
                                                                predxy[1])^2 +
                                                               (thisxy[2] -
                                                                predxy[2])^2) /
                                sqrt((F[1] + F[2])/2)
                        }

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

                    if(!RTMB:::ad_context()){
                        resid.atags[r,1] <- resid.atags.fine[[a]][nrow(tag),1]
                        resid.atags[r,2] <- resid.atags.fine[[a]][nrow(tag),2]
                        resid.atags[r,3] <- resid.atags.fine[[a]][nrow(tag),3]
                    }

                    if(flag.effort){
                        ## Likelihood
                        loglik.atags <- loglik.atags + lpC[r]

                        ## Residuals
                        if(!RTMB:::ad_context()){
                            resid.atags[a,4] <- qnorm(runif(1,
                                                             exp(lpC.cum[r-1]),
                                                             exp(lpC.cum[r-1]) +
                                                             fish.mort / tot.mort *
                                                             (1 - exp(-tot.mort * dt))))
                        }
                    }
                }
            }

            nll <- nll - loglik.atags

            REPORT(resid.atags)
            REPORT(resid.atags.fine)
        }

    }else{

        ## Matrix exponential ---------------------------------

        ## Conventional tags
        if(dat$use.ctags){

            nctags <- nrow(dat$ctags)
            resid.ctags <- matrix(NA, nctags, 4)

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

                        ## Taxis -----------------------
                        if(dat$use.taxis){
                            move <- habitat.tax$grad(dat$xygrid,
                                                     dat$time.cont[itabs]) * dt.ctags
                            if(dat$use.boundaries){
                                move <- move * bound$val(dat$xygrid,
                                                         dat$time.cont[itabs])
                            }
                            xyind <- c(2,2,1,1)
                            posneg <- c(0.5,-0.5,-0.5,0.5)
                            for(j in 2:ncol(dat$nextTo)){
                                ind <- which(!is.na(dat$nextTo[,j]))
                                Zstar[cbind(ind, dat$nextTo[ind,j])] <-
                                    posneg[j-1] * move[ind,xyind[j-1]] /
                                    dat$next.dist[j-1]
                            }
                        }

                        ## Diffusion rate --------------------
                        hD <- exp(habitat.dif$val(dat$xygrid,
                                                  dat$time.cont[itabs])) * dt.ctags
                        if(dat$use.boundaries){
                            hD <- hD * bound$val(dat$xygrid,
                                                 dat$time.cont[itabs])
                        }
                        for(j in 2:ncol(dat$nextTo)){
                            ind <- which(!is.na(dat$nextTo[,j]))
                            Dstar[cbind(ind, dat$nextTo[ind,j])] <-
                                2 * hD[ind] / dat$next.dist[j-1]^2
                        }

                        ## Advection -----------------------
                        if(dat$use.advection){
                            move <- c(habitat.adv.x$val(dat$xygrid,
                                                        dat$time.cont[itabs]),
                                      habitat.adv.y$val(dat$xygrid,
                                                        dat$time.cont[itabs])) * dt.ctags
                            if(dat$use.boundaries){
                                move <- move * bound$val(dat$xygrid,
                                                         dat$time.cont[itabs])
                            }
                            xyind <- c(2,2,1,1)
                            posneg <- c(0.5,-0.5,-0.5,0.5)
                            for(j in 2:ncol(dat$nextTo)){
                                ind <- which(!is.na(dat$nextTo[,j]))
                                Astar[cbind(ind, dat$nextTo[ind,j])] <-
                                    posneg[j-1] * move[ind,xyind[j-1]] /
                                    dat$next.dist[j-1]
                            }
                        }

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
                                    effort$val(dat$xygrid, dat$time.cont[itabs]) * dt.ctags
                                ## how does this work with multiple fleets?
                                Cstar[1:nc,nc+e] <- fish.mort
                            }

                            ## Instantaneous natural mortality rate
                            Ostar[1:nc,ncem] <- rep(nat.mort, nc) * dt.ctags

                            ## Mass balance
                            Cstar[cbind(1:ncem,1:ncem)] <- - RTMB::rowSums(Cstar)
                            Ostar[cbind(1:ncem,1:ncem)] <- - RTMB::rowSums(Ostar)

                            ## Movement rates
                            Mstar <- Mstar + Cstar + Ostar
                        }

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
                                    resid.ctags[ind.tag,1] <- qnorm(runif(1,
                                                                   xres.cum[dat$igrid$idx[recapLoc]],
                                                                   xres.cum[dat$igrid$idx[recapLoc]+1]))
                                    yres.cum <- c(0, cumsum(dist.prob[recapTime+1, dat$igrid$idx == dat$igrid$idx[recapLoc]]))
                                    yres.cum <- yres.cum / yres.cum[length(yres.cum)]
                                    resid.ctags[ind.tag,2] <- qnorm(runif(1,
                                                                   yres.cum[dat$igrid$idy[recapLoc]],
                                                                   yres.cum[dat$igrid$idy[recapLoc]+1]))

                                    if(flag.effort){
                                        resid.ctags[ind.tag,4] <- qnorm(runif(1,
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
                                    resid.ctags[ind.tag,4] <- qnorm(runif(1, exp(lpC.cum), 1))
                                }

                            }
                        }
                    } ## end tag loop

                } ## end if itmax == 0

            } ## end release location loop



            nll <- nll - loglik.ctags

            REPORT(resid.ctags)
        }

        if(dat$use.atags){

            natags <- length(dat$atags)

            resid.atags.fine <- vector("list", natags)
            resid.atags <- matrix(NA, natags, 4)

            ## Loop over atags
            for(a in 1:natags){

                tag <- dat$atags[[a]]
                itrel <- tag$it[1]
                icrel <- tag$ic[1]
                itmax <- nrow(tag)
                itrec <- tag$it[itmax]

                resid.atags.fine[[a]] <- matrix(NA, itmax, 4)

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

                        ## Taxis -----------------------
                        if(dat$use.taxis){
                            move <- habitat.tax$grad(dat$xygrid,
                                                 dat$time.cont[itabs]) * dt.atags
                            if(dat$use.boundaries){
                                move <- move * bound$val(dat$xygrid,
                                                         dat$time.cont[itabs])
                            }
                            xyind <- c(2,2,1,1)
                            posneg <- c(0.5,-0.5,-0.5,0.5)
                            for(j in 2:ncol(dat$nextTo)){
                                ind <- which(!is.na(dat$nextTo[,j]))
                                Zstar[cbind(ind, dat$nextTo[ind,j])] <-
                                    posneg[j-1] * move[ind,xyind[j-1]] /
                                    dat$next.dist[j-1]
                            }
                        }

                        ## Diffusion rate --------------------
                        hD <- exp(habitat.dif$val(dat$xygrid,
                                                  dat$time.cont[itabs])) * dt.atags
                        if(dat$use.boundaries){
                            hD <- hD * bound$val(dat$xygrid,
                                                 dat$time.cont[itabs])
                        }
                        for(j in 2:ncol(dat$nextTo)){
                            ind <- which(!is.na(dat$nextTo[,j]))
                            Dstar[cbind(ind, dat$nextTo[ind,j])] <-
                                2 * hD[ind] / dat$next.dist[j-1]^2
                        }

                        ## Advection -----------------------
                        if(dat$use.advection){
                            move <- c(habitat.adv.x$val(dat$xygrid,
                                                        dat$time.cont[itabs]),
                                      habitat.adv.y$val(dat$xygrid,
                                                        dat$time.cont[itabs])) * dt.atags
                            if(dat$use.boundaries){
                                move <- move * bound$val(dat$xygrid,
                                                         dat$time.cont[itabs])
                            }
                            xyind <- c(2,2,1,1)
                            posneg <- c(0.5,-0.5,-0.5,0.5)
                            for(j in 2:ncol(dat$nextTo)){
                                ind <- which(!is.na(dat$nextTo[,j]))
                                Astar[cbind(ind, dat$nextTo[ind,j])] <-
                                    posneg[j-1] * move[ind,xyind[j-1]] /
                                    dat$next.dist[j-1]
                            }
                        }

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
                                    effort$val(dat$xygrid, dat$time.cont[itabs]) * dt.atags
                                ## how does this work with multiple fleets?
                                Cstar[1:nc,nc+e] <- fish.mort
                            }

                            ## Instantaneous natural mortality rate
                            Ostar[1:nc,ncem] <- rep(nat.mort, nc) * dt.atags

                            ## Mass balance
                            Cstar[cbind(1:ncem,1:ncem)] <- - RTMB::rowSums(Cstar)
                            Ostar[cbind(1:ncem,1:ncem)] <- - RTMB::rowSums(Ostar)

                            ## Movement rates
                            Mstar <- Mstar + Cstar + Ostar
                        }

                        ## Matrix exponential
                        M <- matexpo(Mstar, dat$log2steps)

                        ## Distribution probability
                        dist.prob[t,] <- RTMB:::matrix(RTMB:::matrix(dist.prob[t-1,], 1, ncem) %*% M, 1, ncem)

                    } ## end of time loop


                    ## Likelihood contribution -----------------------------
                    for(t in 2:(itmax-1)){
                        oprob <- dnorm(dat$xygrid[,1], tag$x[t],
                                       sdObsATS / dat$dxdy[1]) *
                            dnorm(dat$xygrid[,2], tag$y[t],
                                  sdObsATS/ dat$dxdy[2])
                        ## oprob <- oprob / (sum(oprob) * prod(dat$dxdy))
                        oprob <- oprob / sum(oprob)
                        loglik.atags <- loglik.atags +
                            log(sum(dist.prob[t,1:nc] * oprob))
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
                            resid.atags.fine[[a]][t,1] <- qnorm(runif(1,
                                                             xres.cum[dat$igrid$idx[recapLoc]],
                                                             xres.cum[dat$igrid$idx[recapLoc] + 1]))
                            yres.cum <- c(0,
                                          cumsum(dist.prob[t,
                                                           dat$igrid$idx ==
                                                           dat$igrid$idx[recapLoc]]))
                            yres.cum <- yres.cum / yres.cum[length(yres.cum)]
                            resid.atags.fine[[a]][t,2] <- qnorm(runif(1,
                                                             yres.cum[dat$igrid$idy[recapLoc]],
                                                             yres.cum[dat$igrid$idy[recapLoc] + 1]))
                        }
                        resid.atags[a,1] <- resid.atags.fine[[a]][nrow(tag),1]
                        resid.atags[a,2] <- resid.atags.fine[[a]][nrow(tag),2]
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
                            resid.atags[a,4] <- runif(1,
                                                      exp(lpC.cum),
                                                      exp(lpC.cum) + exp(lpC))
                        }
                    }

                } ## use this a-tag
            } ## end of atag loop

            nll <- nll - loglik.atags

            REPORT(resid.atags)
            REPORT(resid.atags.fine)

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
