##' Fit the movement model
##'
##' @description The movement model, momo, allows estimation of habitat
##'     preferences and fine-scale movement rates using various types of tagging
##'     data. It supports flexible spatiotemporal modeling with options for
##'     Kalman filter or matrix exponential approaches, and includes tools for
##'     data preparation, model fitting, and visualization.
##'
##' @param dat data frame with input data as produced by the function
##'     [check.momo.data].
##' @param conf configuration list as produced by the function [def.conf]. If
##'     `NULL` (default), [def.conf] is used to generate the default
##'     configuration list based on the input data.
##' @param par parameter list with initial values as produced by the function
##'     [def.par]. If `NULL` (default), [def.par] is used to generate the
##'     default parameter list based on the input data and configuration list.
##' @param map list with parameter mapping as produced by the function
##'     [def.map]. If `NULL` (default), [def.map] is used to generate the
##'     default mapping list based on the input data and configuration and
##'     parameter lists.
##' @param run logical; If `FALSE`, the AD object is returned with the
##'     optimization. Default: `FALSE`.
##' @param lower by default, [get.lower.bounds] is used.
##' @param upper by default, [get.upper.bounds] is used.
##' @param rel.tol option passed to [stats::nlminb] sets the convergence
##'     criteria. Default: `1e-10`.
##' @param do.sdreport logical; If `FALSE`, [RTMB::sdreport] is not run and no
##'     parameter uncertainties are returned. Default: `TRUE`.
##' @param do.report logical; If `FALSE`, `obj$report()` is not run and RTMB
##'     variables might not be reported. Default: `TRUE`.
##' @param use.expm logical; allows to overwrite setting in configuration list
##'     (`conf$use.expm`). If `NULL`, setting from configuration list is used.
##'     Default: `NULL`.
##' @param use.rel.events logical; allows to overwrite setting in configuration
##'     list (`conf$use.rel.events`). If `NULL`, setting from configuration list
##'     is used. Default: `NULL`.
##' @param verbose if `TRUE`, print information to console. Default: `TRUE`.
##' @param ... extra arguments to [RTMB::MakeADFun].
##'
##' @return A list of class `momo.fit`.
##'
##' @examples
##' data(skjepo)
##'
##' fit <- fit.momo(skjepo)
##'
##' @export
fit.momo <- function(dat,
                     conf = NULL,
                     par = NULL,
                     map = NULL,
                     run = TRUE,
                     lower = get.lower.bounds(par),
                     upper = get.upper.bounds(par),
                     rel.tol = 1e-10,
                     do.sdreport = TRUE,
                     do.report = TRUE,
                     use.expm = NULL,
                     use.rel.events = NULL,
                     verbose = TRUE,
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
    if(!conf$use.rel.events && !conf$use.expm && !is.null(dat$ctags)){
        dat$rel.events <- dat$ctags
        dat$ctags$rel.event <- 1:nrow(dat$ctags)
    }

    ## Combine conf and dat
    tmb.all <- c(dat, conf)
    tmb.all$use.ukf <- FALSE ## LATER: make option later

    if(verbose) writeLines("Building the model, that can take a few minutes.")

    t1 <- Sys.time()
    obj <- RTMB::MakeADFun(func = cmb(nll, tmb.all),
                           parameters = par,
                           map = map,
                           silent = TRUE,
                           ...)
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

    res$times <- c(makeadfun = signif(as.numeric(difftime(t2, t1, units = "mins")),2), nlminb = signif(as.numeric(difftime(t3, t2, units = "mins")),2))

    attr(res, "RemoteSha") <- substr(packageDescription("momo")$RemoteSha, 1, 12)
    attr(res, "Version") <- packageDescription("momo")$Version
    res <- add.class(res, "momo.fit")

    if(do.sdreport){
        res <- add.sdreport(res)
    }

    if(do.report){
        res <- add.report(res)
    }

    return(res)
}



##' Add TMB's sdreport
##'
##' @description Run and add TMB's sdreporting function and add it to the fitted
##'     list.
##'
##' @param fit A fitted momo list of class `momo.fit` as returned by the
##'     function [fit.momo].
##'
##' @return A list of class `momo.fit` as returned by the function [fit.momo]
##'     with sdreport.
##'
##' @export
add.sdreport <- function(fit){

    check.class(fit, "momo.fit")

    res <- fit

    t3 <- Sys.time()
    sdrep <- RTMB::sdreport(obj = fit$obj)
    t4 <- Sys.time()

    pl <- as.list(sdrep, "Est")
    plsd <- as.list(sdrep, "Std")

    sdrep$cov <- NULL ## save memory

    res <- c(list(sdrep = sdrep,
                  pl = pl,
                  plsd = plsd),
             res)

    res$times <- c(res$times,
                   sdreport = signif(as.numeric(difftime(t4, t3, units = "mins")),2))


    res <- add.class(res, "momo.fit")
    return(res)
}

##' Add TMB's report
##'
##' @description Run and add TMB's reporting function and add it to the fitted
##'     list.
##'
##' @param fit A fitted momo list of class `momo.fit` as returned by the
##'     function [fit.momo].
##'
##' @return A list of class `momo.fit` as returned by the function [fit.momo]
##'     with reported quantites.
##'
##' @export
add.report <- function(fit){

    check.class(fit, "momo.fit")

    rep <- fit$obj$report()
    fit$rep <- rep

    fit <- add.class(fit, "momo.fit")
    return(fit)
}
