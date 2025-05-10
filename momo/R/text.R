##' Description of model
##'
##' @description Writes a string to install the version of the package which was
##'     used to run the model.
##'
##' @param fit A fitted momo list of class `momo.fit` as returned by the
##'     function [fit.momo].
##' @param ... Additional parameters to be passed to [writeLines].
##'
##' @export
model.version.info <-function(fit, ...){
    check.class(fit, "momo.fit")

    ret <- c(
        '# The fit was run with a specific version of movemod package.',
        '# If in the mean time version on your system has been updated',
        '# you can revert back to the version used by inserting this:',
        '',
        paste0('devtools::install_github("tokami/momo/momo@',
               attr(fit, "RemoteSha"),'")'),
        '',
        '# right before the momo package is loaded'
    )

    writeLines(ret, ...)
}
