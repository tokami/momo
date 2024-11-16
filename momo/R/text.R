##' Description of model
##' @param fit returned object from fit.movemod
##' @param ... Additional parameters to be passed to ...
##' @details Writes a string to install the version of the package which was
##'     used to run the model.
##' @export
model.version.info <-function(fit, ...){
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
