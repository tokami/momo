##' Fit the movement model
##'
##' @param dat list of dat objects
##' @param conf list of initial model parameters
##' @param par list of initial model parameters
##'
##' @details ...
##'
##' @return list of result of minimization conducted via nlminb
##'
##' @export
##'
##' @examples
##' 3 + 3
fit.momo <- function(dat,
                     conf = NULL,
                     par = NULL,
                     map = NULL,
                     newtonsteps = 3,
                     rm.unidentified = FALSE,
                     run = TRUE,
                     lower = get.lower.bounds(par),
                     upper = get.upper.bounds(par),
                     sim.condRE = TRUE,
                     ignore.parm.uncertainty = FALSE,
                     rel.tol = 1e-10,
                     verbose = TRUE,
                     do.sdreport = TRUE,
                     do.report = TRUE,
                     use.expm = NULL,
                     use.rel.events = FALSE,
                     ...){

    ## Flags
    sim.flag <- ifelse(inherits(dat, "momo.sim"), TRUE, FALSE)

    ## RTMB does not allow you to pass the dat to MakeADFun
    cmb <- function(f, d) function(p) f(p, d)

    if(sim.flag){
        sim <- dat
        dat <- sim$dat
        conf <- sim$conf
        par <- sim$par
        map <- sim$map
    }

    if(is.null(conf)) conf <- def.conf(dat)
    if(!is.null(use.expm)) conf$use.expm <- use.expm
    if(!is.null(use.rel.events)) conf$use.rel.events <- use.rel.events
    if(is.null(par)) par <- def.par(dat, conf)
    if(is.null(map)) map <- def.map(dat, conf, par)

    ## Do not use release events for KF
    if(!conf$use.rel.events && !conf$use.expm){
        dat$rel.events <- dat$ctags
        dat$ctags$rel.event <- 1:nrow(dat$ctags)
    }

    ## Combine conf and dat
    tmb.all <- c(dat, conf)

    if(verbose) writeLines("Building the model, that can take a few minutes.")


    t1 <- Sys.time()
    obj <- RTMB::MakeADFun(func = cmb(nll, tmb.all),
                           parameters = par,
                           map = map,
                           silent = TRUE)
    t2 <- Sys.time()

    lower2 <- rep(-Inf, length(obj$par))
    upper2 <- rep(Inf, length(obj$par))
    for(nn in names(lower)) lower2[names(obj$par) == nn] <- lower[[nn]]
    for(nn in names(upper)) upper2[names(obj$par) == nn] <- upper[[nn]]

    if(!run) return(list(sdrep = NA,
                         pl = parameters,
                         plsd = NA,
                         dat = dat,
                         conf = conf,
                         opt = NA,
                         obj = obj))

    if(verbose) writeLines(paste0("Model built (",
                                  signif(as.numeric(difftime(t2, t1,
                                                             units = "mins")),2),
                                  "min). Minimizing neg. loglik."))

    opt <- nlminb(obj$par, obj$fn, obj$gr,
                  control = list(trace = as.integer(verbose),
                                 eval.max = 2000,
                                 iter.max = 1000,
                                 rel.tol = rel.tol),
                  lower = lower2,
                  upper = upper2)
    t3 <- Sys.time()

    if(verbose) writeLines(paste0("Minimization done (",
                                  signif(as.numeric(difftime(t3, t2,
                                                             units = "mins")),2),
                                  "min). Model ", "not "[opt$convergence],
                                  "converged. Estimating uncertainty."))

    res <- list(dat = dat,
                conf = conf,
                map = map,
                opt = opt,
                obj = obj,
                low = lower,
                hig = upper)

    times <- c(makeadfun = signif(as.numeric(difftime(t2, t1, units = "mins")),2), nlminb = signif(as.numeric(difftime(t3, t2, units = "mins")),2))


    if(do.report){
        rep <- obj$report()
        res$rep <- rep
    }

    if(do.sdreport){
        sdrep <- sdreport(obj = obj,
                          par.fixed = opt$par,
                          ignore.parm.uncertainty = ignore.parm.uncertainty)
        t4 <- Sys.time()

        pl <- as.list(sdrep, "Est")
        plsd <- as.list(sdrep, "Std")

        ## TODO: make optional
        sdrep$cov <- NULL ## save memory

        res <- c(list(sdrep = sdrep,
                      pl = pl,
                      plsd = plsd),
                 res)

        times <- c(times, sdreport = signif(as.numeric(difftime(t4, t3, units = "mins")),2))
    }

    res$times <- times

    attr(res, "RemoteSha") <- substr(packageDescription("momo")$RemoteSha, 1, 12)
    attr(res, "Version") <- packageDescription("momo")$Version
    res <- add.class(res, "momo.fit")

    return(res)
}


##' Lower bounds
##' @param parameters initial values for the model in a format similar to what
##'     is returned from the defpar function
##' @return a named list
get.lower.bounds <- function(parameters){
    list()
}

##' Upper bounds
##' @param parameters initial values for the model in a format similar to what
##'     is returned from the defpar function
##' @return a named list
get.upper.bounds <- function(parameters){
    list()
}
