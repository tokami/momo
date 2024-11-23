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
                     conf,
                     par,
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
                     ...){

    ## RTMB does not allow you to pass the dat to MakeADFun
    cmb <- function(f, d) function(p) f(p, d)

    if(is.null(map)) map <- def.map(dat, conf, par)

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
                                  "mins). Minimizing neg. loglik."))

    t1 <- Sys.time()
    opt <- nlminb(obj$par, obj$fn, obj$gr,
                  control = list(trace = as.integer(verbose),
                                 eval.max = 2000,
                                 iter.max = 1000,
                                 rel.tol = rel.tol),
                  lower = lower2,
                  upper = upper2)
    t2 <- Sys.time()

    if(verbose) writeLines(paste0("Minimization done (",
                                  signif(as.numeric(difftime(t2, t1,
                                                             units = "mins")),2),
                                  "mins). Model ", "not "[opt$convergence],
                                  "converged. Estimating uncertainty."))

    res <- list(dat = dat,
                conf = conf,
                opt = opt,
                obj = obj,
                low = lower,
                hig = upper)

    if(do.report){
        rep <- obj$report()
        res$rep <- rep
    }


    if(do.sdreport){
        sdrep <- sdreport(obj = obj,
                          par.fixed = opt$par,
                          ignore.parm.uncertainty = ignore.parm.uncertainty)

        pl <- as.list(sdrep, "Est")
        plsd <- as.list(sdrep, "Std")

        sdrep$cov <- NULL ## save memory

        res <- c(list(sdrep = sdrep,
                      pl = pl,
                      plsd = plsd),
                 res)
    }

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
