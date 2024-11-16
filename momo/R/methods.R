##' Do sd report
##'
##' # @method sdreport momo.fit
##' ## TODO: why does the method sdreport.momo.fit not work?
##'
##' @param fit fit
##'
##' @details ...
##'
##' @return list of result of minimization conducted via nlminb
##'
##' @export
##'
##' @examples
##' 3 + 3
do.sdreport <- function(fit,
                        ignore.parm.uncertainty = FALSE){

    res <- fit

    sdrep <- RTMB::sdreport(obj = fit$obj,
                            par.fixed = fit$opt$par,
                            ignore.parm.uncertainty = ignore.parm.uncertainty)

    pl <- as.list(sdrep, "Est")
    plsd <- as.list(sdrep, "Std")

    sdrep$cov <- NULL ## save memory

    res <- c(list(sdrep = sdrep,
                  pl = pl,
                  plsd = plsd),
             res)

    attr(res, "RemoteSha") <- substr(packageDescription("momo")$RemoteSha, 1, 12)
    attr(res, "Version") <- packageDescription("momo")$Version
    res <- add.class(res, "momo.fit")

    return(res)
}
