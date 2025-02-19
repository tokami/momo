##' Plot momo grid
##'
##' @param grid grid
##'
##' @export
plotmomo.grid <- function(grid,
                          main = "Grid",
                          labels = TRUE,
                          plot.land = FALSE,
                          keep.par = FALSE,
                          ...){

    xlims <- attributes(grid)$xrange
    ylims <- attributes(grid)$yrange

    if(!keep.par){
        opar <- par(no.readonly = TRUE)
        on.exit(par(opar))
        par(mfrow = c(1,1))
    }

    plot(xlims, ylims,
         xlim = xlims,
         ylim = ylims,
         main = main,
         ty = "n",
         xaxs = "i", yaxs = "i",
         xlab = "x", ylab = "y")
    c0 <- t(grid$celltable)
    c0[c0 > 0] <- 1
    image(attributes(grid)$xgr,
          attributes(grid)$ygr,
          t(c0),
          col = adjustcolor("dodgerblue2",0.2),
          xlim = xlims,
          ylim = ylims,
          add = TRUE)
    if(plot.land){
        try(maps::map("world",
                  xlim = xlims,
                  ylim = ylims,
                  fill = TRUE, plot = TRUE, add = TRUE,
                  col = grey(0.95), border = grey(0.5)), silent = TRUE)
    }
    labs <- as.numeric(grid$celltable)
    labs <- labs[!is.na(labs)]
    if(labels) text(grid$xygrid[,1], grid$xygrid[,2], labs)
    abline(v = attributes(grid)$xgr)
    abline(h = attributes(grid)$ygr)
    box(lwd = 1.5)

    return(invisible(NULL))
}


##' Plot momo env
##' @export
plotmomo.env <- function(x,
                         select = NULL,
                         main = "Environmental fields",
                         labels = TRUE,
                         keep.par = FALSE,
                         plot.land = FALSE,
                         xlab = "x",
                         ylab = "y",
                         ...){

    if(inherits(x, "momo.sim")){
        env <- x$env
    }else{
        env <- x
    }

    if(inherits(env, "list")){
        stop("Not implemented yet for lists, plot one list element after the other!")
    }

    if(!is.null(select)){
        env <- env[,,select]
    }

    nt <- dim(env)[3]

    if(any(names(attributes(env)) == "dimnames")){
        xlims <- range(as.numeric(attributes(env)$dimnames[[1]]))
        ylims <- range(as.numeric(attributes(env)$dimnames[[2]]))
    }else{
        xlims <- c(1,nrow(env))
        ylims <- c(1,nrow(env))
    }

    if(!keep.par){
        opar <- par(no.readonly = TRUE)
        on.exit(par(opar))
        par(mfrow = n2mfrow(nt, asp = 2),
            mar = c(1.5,1.5,1.5,1.5),
            oma = c(3,3,ifelse(main == "", 0, 1.5),0))
    }
    for(i in 1:nt){
        x <- as.numeric(rownames(env[,,i]))
        if(length(x) == 0) x <- 1:nrow(env[,,i])
        y <- as.numeric(colnames(env[,,i]))
        if(length(y) == 0) y <- 1:ncol(env[,,i])
        plot(1,1, type = "n",
             xlim = xlims, ylim = ylims,
             xlab = "",
             ylab = "")
        image(x, y, env[,,i], col = terrain.colors(100), add = TRUE)
        if(plot.land){
            try(maps::map("world",
                      xlim = xlims,
                      ylim = ylims,
                      fill = TRUE, plot = TRUE, add = TRUE,
                      col = adjustcolor(grey(0.7),0.5),
                      border = grey(0.5)), silent = TRUE)
        }
        contour(x, y, env[,,i], add = TRUE)
        if(nt > 1) legend("topleft", legend = paste0("Field ", i),
                          bg = "white", pch = NA)
        box(lwd = 1.5)
    }
    if(!keep.par){
        mtext(main, 3, 0, outer = TRUE)
        mtext(xlab, 1, 1, outer = TRUE)
        mtext(ylab, 2, 1.5, outer = TRUE)
    }


    return(invisible(NULL))
}

##' Plot momo effort
##' @export
plotmomo.effort <- function(effort,
                            main = "Effort fields",
                            labels = TRUE,
                            keep.par = FALSE,
                            plot.land = FALSE,
                            xlab = "x",
                            ylab = "y",
                            ...){

    nt <- dim(effort)[3]

    if(any(names(attributes(effort)) == "dimnames")){
        xlims <- range(as.numeric(attributes(effort)$dimnames[[1]]))
        ylims <- range(as.numeric(attributes(effort)$dimnames[[2]]))
    }else{
        xlims <- c(1,nrow(effort))
        ylims <- c(1,nrow(effort))
    }

    if(!keep.par){
        opar <- par(no.readonly = TRUE)
        on.exit(par(opar))
        par(mfrow = n2mfrow(nt),
            mar = c(1.5,1.5,1.5,1.5),
            oma = c(3,3,1.5,0))
    }
    for(i in 1:nt){
        x <- as.numeric(rownames(effort[,,i]))
        if(length(x) == 0) x <- 1:nrow(effort[,,i])
        y <- as.numeric(colnames(effort[,,i]))
        if(length(y) == 0) y <- 1:ncol(effort[,,i])
        plot(1,1, type = "n",
             xlim = xlims, ylim = ylims,
             xlab = "",
             ylab = "")
        image(x, y, effort[,,i], col = terrain.colors(100), add = TRUE)
        contour(x, y, effort[,,i], add = TRUE)
        legend("topleft", legend = paste0("Field ", i),
               bg = "white", pch = NA)
        if(plot.land){
            try(maps::map("world",
                      xlim = xlims,
                      ylim = ylims,
                      fill = TRUE, plot = TRUE, add = TRUE,
                      col = grey(0.95), border = grey(0.5)), silent = TRUE)
        }
        box(lwd = 1.5)
    }
    mtext(xlab, 1, 1, outer = TRUE)
    mtext(ylab, 2, 1, outer = TRUE)
    mtext(main, 3, 0, outer = TRUE)


    return(invisible(NULL))
}

##' Plot momo ctags
##' @export
plotmomo.ctags <- function(x,
                           main = "Conventional tags",
                           plot.land = FALSE,
                           keep.par = FALSE,
                           xlim = NULL,
                           ylim = NULL,
                           add = FALSE,
                           ...){

    if(inherits(x, "momo.sim") || inherits(x, "momo.data")){
        ctags <- x$ctags
    }else{
        ctags <- x
    }


    if(is.null(xlim)){
        xlims <- range(ctags$x0, ctags$x1, na.rm = TRUE)
    }else{
        xlims <- xlim
    }
    if(is.null(ylim)){
        ylims <- range(ctags$y0, ctags$y1, na.rm = TRUE)
    }else{
        ylims <- ylim
    }

    if(!keep.par & !add){
        opar <- par(no.readonly = TRUE)
        on.exit(par(opar))
        par(mfrow = c(1,1))
    }

    if(!add){
        plot(0, 0, ty = "n", main = main,
             xlim = xlims,
             ylim = ylims,
             xlab = "x", ylab = "y")
    }
    arrows(ctags$x0, ctags$y0, ctags$x1, ctags$y1,
           col = adjustcolor("grey60",0.4),
           length = 0.1)
    points(ctags$x0, ctags$y0, pch = 16, col = "grey30", cex = 0.8)
    if(plot.land && !add){
        try(maps::map("world",
                  xlim = xlims,
                  ylim = ylims,
                  fill = TRUE, plot = TRUE, add = TRUE,
                  col = grey(0.95), border = grey(0.5)), silent = TRUE)
    }
    if(!add){
        box(lwd = 1.5)
    }

    return(invisible(NULL))
}


##' Plot momo atags
##' @export
plotmomo.atags <- function(x,
                           main = "Archival tags",
                           keep.par = FALSE,
                           plot.land = FALSE,
                           xlim = NULL,
                           ylim = NULL,
                           xlab = "x",
                           ylab = "y",
                           leg.pos = "topright",
                           ...){

    if(inherits(x, "momo.sim") || inherits(x, "momo.data")){
        atags <- x$atags
    }else{
        atags <- x
    }

    if(is.null(xlim)){
        xlims <- range(sapply(atags, function(x) range(x[,2], na.rm = TRUE)))
    }else{
        xlims <- xlim
    }
    if(is.null(ylim)){
        ylims <- range(sapply(atags, function(x) range(x[,3], na.rm = TRUE)))
    }else{
        ylims <- ylim
    }

    if(!keep.par){
        opar <- par(no.readonly = TRUE)
        on.exit(par(opar))
        par(mfrow = c(1,1))
    }

    cols <- get.momo.cols(2)

    plot(0,0, ty = "n", main = main,
         xlim = xlims,
         ylim = ylims,
         xlab = xlab, ylab = ylab)
    if(plot.land){
        try(maps::map("world",
                  xlim = xlims,
                  ylim = ylims,
                  fill = TRUE, plot = TRUE, add = TRUE,
                  col = grey(0.95), border = grey(0.5)), silent = TRUE)
    }
    for(i in 1:length(atags)){
        points(atags[[i]][1,2], atags[[i]][1,3], col = cols[1], pch = 1)
        lines(atags[[i]][,2], atags[[i]][,3], col = adjustcolor("grey60",0.3))
        points(atags[[i]][nrow(atags[[i]]),2], atags[[i]][nrow(atags[[i]]),3],
               col = cols[2], pch = 0)
    }
    legend(leg.pos, legend = c("Release", "Recapture"),
           pch = c(1,0),
           col = cols,
           bg = "white")
    box(lwd = 1.5)

    return(invisible(NULL))
}


##' Plot preference
##' @export
plotmomo.pref <- function(x,
                          type = "taxis",
                          col = "black",
                          lwd = 1,
                          ci = 0.95,
                          par = NULL,
                          funcs = NULL,
                          env = NULL,
                          keep.par = FALSE,
                          add = FALSE,
                          ...){

    if(!keep.par){
        opar <- par(no.readonly = TRUE)
        on.exit(par(opar))
        par(mfrow = c(1,1))
    }


    if(inherits(x, "momo.fit")){

        sdr <- x$sdr
        env.pred <- x$dat$env.pred

        ## TODO: make matrices if more than one env field!
        ## improve this code, convert sdr$values into matrices!
        ## TODO: make polygon and CI plotting conditional if fit was run with sdreport=TRUE!

        i = 1

        if(type == "taxis"){

            if(!is.null(sdr)){
                ind <- which(names(sdr$value) == "prefT.pred")
                par.est <- x$pl$alpha[,i,]
            }else{
                ind <- which(names(x$rep) == "prefT.pred")
                par.est <- c(0,x$opt$par[names(x$opt$par) == "alpha"])
            }
            knots <- x$dat$knots.tax[,i]


        }else if(type == "diffusion"){

            if(!is.null(sdr)){
                ind <- which(names(sdr$value) == "prefD.pred")
                par.est <- x$pl$beta[,i]
            }else{
                ind <- which(names(x$rep) == "prefD.pred")
                par.est <- x$opt$par[names(x$opt$par) == "beta"]
            }
            knots <- x$dat$knots.dif[,i]

        }else stop("only taxis and diffusion implemented yet.")

        if(!is.null(sdr)){
            pref <- sdr$value[ind]
            prefsd <- sdr$sd[ind]
            preflow <- pref - qnorm(ci + (1 - ci)/2) * prefsd
            prefup <- pref + qnorm(ci + (1 - ci)/2) * prefsd
        }else{
            pref <- x$rep[[ind]]
            prefsd <- preflow <- prefup <- NULL
        }

        xlims <- apply(env.pred, 2, range)

        ylims <- range(pref, preflow, prefup, na.rm = TRUE) ## if more env fields this should be matrices
        alpha <- 0.3

        if(!add){
            plot(NA, ty = 'n',
                 xlim = xlims,
                 ylim = ylims,
                 ylab = "Preference",
                 xlab = "Covariate")
        }
        polygon(c(env.pred[,i], rev(env.pred[,i])),
                c(preflow, rev(prefup)),
                border = NA,
                col = rgb(t(col2rgb(col))/255, alpha=alpha))
        ## rug(x$dat$env$env.obs[,inp$env$var[i]]) ## TODO:

        points(knots, par.est, pch = 16, cex = 1.2, col = col)
        lines(env.pred[,i], pref, col = col, lwd = lwd)

        if(!add) box(lwd = 1.5)

    }else if(inherits(x, "momo.sim")){

        i = 1

        grid <- x$grid
        env <- x$env
        par <- x$par.sim
        dat <- x$dat
        funcs <- NULL

        if(is.null(par)) stop("No parameters provided! Use par = list() to specify parameters for taxis.")

        par <- get.sim.par(par)
        env <- check.that.list(env)
        dat <- setup.momo.data(grid,
                               env,
                               trange = c(0,
                                          max(sapply(env,
                                                     function(x) dim(x)[3]))),
                               knots.tax = dat$knots.tax,
                               knots.dif = dat$knots.dif)
        conf <- def.conf(dat)
        funcs <- get.sim.funcs(funcs, dat, conf, env, par)

        env.pred <- dat$env.pred

        xlims <- apply(dat$env.pred, 2, range)

        if(type == "taxis"){
            knots <- dat$knots.tax[,i]
            par <- par$alpha[,i,]
        }else if(type == "diffusion"){
            knots <- dat$knots.dif[,i]
            par <- par$beta[,i]
        }

        get.true.pref <- momo:::poly.fun(as.numeric(knots),
                                         as.numeric(par))

        pref <- get.true.pref(dat$env.pred[,i])

        ylims <- range(pref, par)

        alpha <- 0.3


        if(!add){
            plot(NA, ty = 'n',
                 xlim = xlims,
                 ylim = ylims,
                 ylab = "Preference",
                 xlab = "Covariate")
        }
        lines(env.pred[,i], pref,
              col = col,
              lwd = lwd)
        points(knots, par,
               col = col,
               pch = 15, cex = 1.2)
        if(!add) box(lwd = 1.5)
    }
}


##' Plot spatial preference
##' @export
plotmomo.pref.spatial <- function(x,
                                  type = "taxis",
                                  col = hcl.colors(14, "YlOrRd", rev = TRUE),
                                  ci = 0.95,
                                  par = NULL,
                                  funcs = NULL,
                                  env = NULL,
                                  keep.par = FALSE,
                                  knots.tax = NULL,
                                  add = FALSE,
                                  ...){

    if(!keep.par){
        opar <- par(no.readonly = TRUE)
        on.exit(par(opar))
        par(mfrow = c(1,1))
    }

    if(inherits(x, "momo.fit")){

        sdr <- x$sdr
        env.pred <- x$dat$env.pred

        ## TODO: make matrices if more than one env field!
        ## improve this code, convert sdr$values into matrices!
        ## TODO: make polygon and CI plotting conditional if fit was run with sdreport=TRUE!

        i = 1

        if(type == "taxis"){

            if(!is.null(sdr)){
                ind <- which(names(sdr$value) == "prefT.pred")
                par.est <- x$pl$alpha[,i,]
            }else{
                ind <- which(names(x$rep) == "prefT.pred")
                par.est <- c(0,x$opt$par[names(x$opt$par) == "alpha"])
            }
            knots <- x$dat$knots.tax[,i]
            if(!is.null(par)) par.true <- par$alpha[,i,]


        }else if(type == "diffusion"){

            if(!is.null(sdr)){
                ind <- which(names(sdr$value) == "prefD.pred")
                par.est <- x$pl$beta[,i]
            }else{
                ind <- which(names(x$rep) == "prefD.pred")
                par.est <- x$opt$par[names(x$opt$par) == "beta"]
            }
            knots <- x$dat$knots.dif[,i]
            if(!is.null(par)) par.true <- par$beta[,i]

        }else stop("only taxis and diffusion implemented yet.")

        if(!is.null(sdr)){
            pref <- sdr$value[ind]
            prefsd <- sdr$sd[ind]
            preflow <- pref - qnorm(ci + (1 - ci)/2) * prefsd
            prefup <- pref + qnorm(ci + (1 - ci)/2) * prefsd
        }else{
            pref <- x$rep[[ind]]
            prefsd <- preflow <- prefup <- NULL
        }

        xlims <- apply(env.pred, 2, range)

        ylims <- range(pref, preflow, prefup) ## if more env fields this should be matrices
        if(!is.null(par)) ylims <- range(ylims, par.true)

        cols <- get.momo.cols(2)
        alpha <- 0.3

        get.true.pref <- momo:::poly.fun(as.numeric(knots),
                                         as.numeric(par.est))

        i = 1
        pref.pred <- get.true.pref(as.numeric(x$dat$env[[i]][,,1]))

        mat <- x$dat$env[[i]][,,1]
        mat[] <- pref.pred

        if(!add){
            plot(NA,
                 xlim = x$dat$xrange,
                 ylim = x$dat$yrange,
                 xlab = "x",
                 ylab = "y")
        }

        image(as.numeric(rownames(mat)),
              as.numeric(colnames(mat)),
              mat,
            xlim = x$dat$xrange,
            ylim = x$dat$yrange,
            add = TRUE)

        box(lwd=1.5)

    }else{

        grid <- x$grid
        env <- x$env
        par <- x$par
        dat <- x$dat
        knots.tax <- x$knots
        funcs <- NULL

        if(is.null(par)) stop("No parameters provided! Use par = list() to specify parameters for taxis.")

        if(!add){
            plot(NA,
                 xlim = attr(grid,"xrange"),
                 ylim = attr(grid,"yrange"),
                 xlab = "x",
                 ylab = "y")
        }

        par <- get.sim.par(par)
        env <- check.that.list(env)
        dat <- setup.momo.data(grid, env, trange = c(0,
                                                     max(sapply(env,
                                                                function(x) dim(x)[3]))))
        conf <- def.conf(dat)
        funcs <- get.sim.funcs(funcs, dat, conf, env, par)

        ## uv.true <- t(apply(x$dat$xygrid, 1, function(xy)
        ##     funcs$tax(xy,NA)))

        get.true.pref <- momo:::poly.fun(as.numeric(x$knots),
                                         as.numeric(x$par$alpha))

        i = 1
        pref.pred <- get.true.pref(as.numeric(env[[i]]))

        image(
            matrix(pref.pred, attr(grid, "nx"),
                   attr(grid,"ny")),
            xlim = attr(grid,"xrange"),
            ylim = attr(grid,"yrange"),
            add = TRUE)

        if(!add) box(lwd = 1.5)

    }
}



##' Plot taxis
##' @export
plotmomo.taxis <- function(x,
                           select = NULL,
                           average = TRUE,
                           cor = 0.1,
                           col = "black",
                           alpha = 0.5,
                           lwd = 1,
                           main = "Taxis",
                           xlab = "x",
                           ylab = "y",
                           xaxt = "s",
                           yaxt = "s",
                           plot.land = FALSE,
                           keep.par = FALSE,
                           add = FALSE,
                           ...){

    if(is.null(select)) select <- 1:length(x$dat$time.cont.pred)

    if(!keep.par){
        opar <- par(no.readonly = TRUE)
        on.exit(par(opar))
        if(average || length(select) == 1){
            par(mfrow = c(1,1))
        }else{
            par(mfrow = n2mfrow(length(select), asp = 2))
        }
    }

    if(inherits(x, "momo.fit")){

        if(average){
            if(length(select) > 1){
                tax.x <- apply(x$rep$hTdx.pred[,select], 1, mean, na.rm = TRUE)
                tax.y <- apply(x$rep$hTdy.pred[,select], 1, mean, na.rm = TRUE)
            }else{
                tax.x <- x$rep$hTdx.pred[,select]
                tax.y <- x$rep$hTdy.pred[,select]
            }
        }else{
            tax.x <- x$rep$hTdx.pred[,select]
            tax.y <- x$rep$hTdy.pred[,select]
        }

        if(!inherits(tax.x, "matrix")){
            tax.x <- as.matrix(tax.x)
            tax.y <- as.matrix(tax.y)
        }

        for(i in 1:ncol(tax.x)){

            if(!add){
                plot(NA,
                     xlim = x$dat$xrange,
                     ylim = x$dat$yrange,
                     xlab = xlab,
                     ylab = ylab,
                     xaxt = xaxt,
                     yaxt = yaxt,
                     main = "")
            }
            if(plot.land){
                try(maps::map("world",
                     xlim = x$dat$xrange,
                     ylim = x$dat$yrange,
                          fill = TRUE, plot = TRUE, add = TRUE,
                     col = adjustcolor(grey(0.95),0.5),
                     border = grey(0.7)), silent = TRUE)
            }

            arrows(x$dat$xygrid.pred[,1],
                   x$dat$xygrid.pred[,2],
                   x$dat$xygrid.pred[,1] + tax.x[,i] * cor,
                   x$dat$xygrid.pred[,2] + tax.y[,i] * cor,
                   col = col,
                   lwd = lwd,
                   length = .1)

            if(!add) box(lwd = 1.5)

        }

    }else if(inherits(x, "momo.sim")){

        grid <- x$grid
        env <- x$env
        par <- x$par.sim
        dat <- x$dat
        funcs <- NULL

        if(is.null(par)) stop("No parameters provided! Use par = list() to specify parameters for taxis.")

        if(!add){
            plot(NA,
                 xlim = attr(grid,"xrange"),
                 ylim = attr(grid,"yrange"),
                 xlab = "x",
                 ylab = "y",
                 main = main)
        }

        par <- get.sim.par(par)
        env <- check.that.list(env)
        dat <- setup.momo.data(grid, env,
                               trange = c(0,
                                          max(sapply(env,
                                                     function(x) dim(x)[3]))),
                               knots.tax = dat$knots.tax,
                               knots.dif = dat$knots.dif)

        conf <- def.conf(dat)
        funcs <- get.sim.funcs(funcs, dat, conf, env, par)
        ## uv.true <- t(apply(x$dat$xygrid, 1, function(xy)
        ##     funcs$tax(xy,NA)))
        hTdx.true <- sapply(dat$time.cont.pred,
                            function(t) apply(dat$xygrid.pred, 1, function(x) funcs$tax(x,t)[1]))
        hTdy.true <- sapply(dat$time.cont.pred,
                            function(t) apply(dat$xygrid.pred, 1, function(x) funcs$tax(x,t)[2]))
        uv.true <- matrix(NA, nrow(dat$xygrid), 2)
        uv.true[,1] <- rowMeans(hTdx.true)
        uv.true[,2] <- rowMeans(hTdy.true)
        arrows(grid$xygrid[,1],
               grid$xygrid[,2],
               grid$xygrid[,1]+uv.true[,1] * cor,
               grid$xygrid[,2]+uv.true[,2] * cor,
               ## col = get.momo.cols(1, alpha, type = "true"),
               col = col,
               lwd = lwd,
               length = .1)
        ## uv.true <- t(apply(x$xygrid, 1, function(xy)
        ##     taxis.fun(xy,1)))
        ## arrows(x$xygrid[,1],
        ##        x$xygrid[,2],
        ##        x$xygrid[,1]+uv.true[,1],
        ##        x$xygrid[,2]+uv.true[,2],
        ##        length=.1)

        if(!add) box(lwd = 1.5)
    }
}


##' Plot diffusion
##' @export
plotmomo.dif <- function(x,
                         cor = NULL,
                         main = "Diffusion",
                         col = "black",
                         alpha = 0.5,
                         lwd = 1,
                         keep.par = FALSE,
                         add = FALSE,
                         ...){

    if(!keep.par){
        opar <- par(no.readonly = TRUE)
        on.exit(par(opar))
        par(mfrow = c(1,1))
    }

    if(inherits(x, "momo.fit")){

        if(!add){
            plot(NA,
                 xlim = x$dat$xrange,
                 ylim = x$dat$yrange,
                 xlab = "x",
                 ylab = "y",
                 main = "")
        }

        dif.est <- exp(apply(x$rep$hD.pred, 1, mean))
        dif.est * x$dat$dxdy[1]

        if(is.null(cor)){
            cor <- min(diff(x$dat$xrange)/2,
                       diff(x$dat$yrange)/2) / max(sqrt(dif.est)) / 20
        }

        points(x$dat$xygrid[,1],
               x$dat$xygrid[,2],
               col = col,
               lwd = lwd,
               cex = sqrt(dif.est) * cor)

    }else if(inherits(x, "momo.sim")){

        grid <- x$grid
        env <- x$env
        par <- x$par.sim
        dat <- x$dat
        funcs <- NULL

        if(is.null(par)) stop("No parameters provided! Use par = list() to specify parameters for taxis.")

        if(!add){
            plot(NA,
                 xlim = attr(grid,"xrange"),
                 ylim = attr(grid,"yrange"),
                 xlab = "x",
                 ylab = "y",
                 main = main)
        }

        par <- get.sim.par(par)
        env <- check.that.list(env)
        dat <- setup.momo.data(grid, env,
                               trange = c(0,
                                          max(sapply(env,
                                                     function(x) dim(x)[3]))),
                               knots.tax = dat$knots.tax,
                               knots.dif = dat$knots.dif)

        conf <- def.conf(dat)
        funcs <- get.sim.funcs(funcs, dat, conf, env, par)

        D.true <- sapply(dat$time.cont.pred,
                         function(t) apply(dat$xygrid.pred, 1,
                                           function(x) exp(funcs$dif(x,t)[1])))

     if(is.null(cor)){
         cor <- min(diff(dat$xrange)/2,
                    diff(dat$yrange)/2) / max(sqrt(rowMeans(D.true))) / 20
        }

        points(dat$xygrid[,1],
               dat$xygrid[,2],
               ## col = get.momo.cols(1, alpha, type = "true"),
               col = col,
               lwd = lwd,
               cex = sqrt(rowMeans(D.true)) * cor)

    }
    if(!add) box(lwd = 1.5)
}


##' Plot residuals
##' @export
plotmomo.resid <- function(fit,
                           add.dist = FALSE,
                           keep.par = FALSE,
                           ...){

    if(!keep.par){
        opar <- par(no.readonly = TRUE)
        on.exit(par(opar))
        par(mfrow = c(1,1))
    }

    if(!inherits(fit, "momo.fit")) stop("fit must be a fitted momo object ('momo.fit')!")

    if(!any(names(fit) == "diag")){
        fit <- get.diag(fit)
    }

    if(fit$conf$use.ctags && fit$conf$use.atags){
        if(fit$conf$use.effort){
            nmati <- 18
        }else{
            nmati <- 12
        }
    }else{
        if(fit$conf$use.effort){
            nmati <- 9
        }else{
            nmati <- 6
        }
    }
    if(add.dist){
        nmati <- nmati + 3
    }
    layout(matrix(1:nmati, 3))
    par(mar = c(4,4,3,2), oma = c(1,1,1,1))

    plot.resid.bias <- function(resid, x,
                                pval = NULL,
                                ylab = "Residuals"){
        plot(x, resid,
             ty = "n",
             xlab = "Time of recapture",
             ylab = ylab)
        abline(h = 0, lty = 3)
        points(x, resid, pch = 1)
        box(lwd = 1.5)
        if(!is.null(pval) && !is.na(pval)){
            title(paste0("Bias p-val: ", signif(pval,5)),
                  col.main = ifelse(pval >= 0.05,
                                    get.momo.cols(1, type = "notsig"),
                                    get.momo.cols(1, type = "sig")))
        }
    }

    plot.resid.acf <- function(resid, pval = NULL, lag.max = 4,
                               ylab = "ACF"){
        ## inds <- which(acf.signf(resid, lag.max = lag.max))
        ## if(length(inds) > 0){
        ##     txt <- paste0("lag.signf: ", paste0(inds, collapse=","))
        ## }else{
        ##     txt <- ""
        ## }
        acf(resid, main='', lag.max = lag.max, ylab = ylab)
        ## legend('topright', legend=NA, title = txt, col=2,
        ##        bty='n', pt.cex=0, text.col=2)
        box(lwd = 1.5)
        if(!is.null(pval) && !is.na(pval)){
            title(paste0("LBox p-val: ", signif(pval,5)),
                  col.main = ifelse(pval >= 0.05,
                                    get.momo.cols(1, type = "notsig"),
                                    get.momo.cols(1, type = "sig")))
        }
    }

    plot.resid.normal <- function(resid, pval = NULL){
        qqnorm(resid,
               main = "")  ## ylim = range(resi[is.finite(resi)]),
        qqline(resid)
        box(lwd = 1.5)
        if(!is.null(pval) && !is.na(pval)){
            title(paste0("Shapiro p-val: ", signif(pval,5)),
                  col.main = ifelse(pval >= 0.05,
                                    get.momo.cols(1, type = "notsig"),
                                    get.momo.cols(1, type = "sig")))
        }
    }

    plot.resid.spatial <- function(resid, x, y,
                                   pval = NULL,
                                   xlab = "x",
                                   ylab = "y",
                                   cor = 1){
        plot(x, y,
             ty = "n",
             xlab = xlab,
             ylab = ylab)
        points(x, y, cex = abs(resid) * cor,
               col = ifelse(resid > 0,
                            get.momo.cols(1, type = "pos"),
                            get.momo.cols(1, type = "neg")))
        box(lwd = 1.5)
        ## if(!is.null(pval)){
        ##     coli <- ifelse(pval >= 0.05,
        ##                        get.momo.cols(1, type = "notsig"),
        ##                        get.momo.cols(1, type = "sig"))
        ##     title(paste0("Moran's I p-val: ", signif(pval,5)),
        ##           col.main = ifelse(is.na(coli), "grey50", coli))
        ## }
    }


    if(fit$conf$use.ctags){

        ## x
        tmp <- data.frame(fit$dat$ctags,
                          resid = fit$rep$resid.ctags[,1])
        tmp <- tmp[!is.na(tmp$resid),]
        tmp <- tmp[order(tmp$t1),]
        plot.resid.bias(tmp$resid, tmp$t1,
                        fit$diag$ctags.x$ttest$p.value,
                        ylab = expression("ctags x residuals"))
        ## hist(resi, main = "", xlab = "Residuals")
        plot.resid.spatial(tmp$resid, tmp$x1, tmp$y1,
                           fit$diag$ctags.x$moransI$p.value)
        ## plot.resid.acf(tmp$resid,
        ##                fit$diag$ctags.x$box$p.value)
        plot.resid.normal(tmp$resid,
                          fit$diag$ctags.x$shapiro$p.value)

        ## y
        tmp <- data.frame(fit$dat$ctags,
                          resid = fit$rep$resid.ctags[,2])
        tmp <- tmp[!is.na(tmp$resid),]
        tmp <- tmp[order(tmp$t1),]
        plot.resid.bias(tmp$resid, tmp$t1,
                        fit$diag$ctags.y$ttest$p.value,
                        ylab = expression("ctags y residuals"))
        ## hist(resi, main = "", xlab = "Residuals")
        plot.resid.spatial(tmp$resid, tmp$x1, tmp$y1,
                           fit$diag$ctags.y$moransI$p.value)
        ## plot.resid.acf(tmp$resid,
        ##                fit$diag$ctags.y$box$p.value)
        plot.resid.normal(tmp$resid,
                          fit$diag$ctags.y$shapiro$p.value)

        if(add.dist){
            ## dist
            tmp <- data.frame(fit$dat$ctags,
                              resid = fit$rep$resid.ctags[,3])
            tmp <- tmp[!is.na(tmp$resid),]
            tmp <- tmp[order(tmp$t1),]
            plot.resid.bias(tmp$resid, tmp$t1,
                            fit$diag$ctags.x$ttest$p.value,
                            ylab = expression("ctags dist residuals"))
            ## hist(resi, main = "", xlab = "Residuals")
            plot.resid.spatial(tmp$resid, tmp$x1, tmp$y1,
                               fit$diag$ctags.x$moransI$p.value)
            ## plot.resid.acf(tmp$resid,
            ##                fit$diag$ctags.x$box$p.value)
            plot.resid.normal(tmp$resid,
                              fit$diag$ctags.x$shapiro$p.value)

        }

        ## t
        if(fit$conf$use.effort){
            tmp <- data.frame(fit$dat$ctags,
                              resid = fit$rep$resid.ctags[,4])
            tmp <- tmp[!is.na(tmp$resid),]
            tmp <- tmp[order(tmp$t0),]
            plot.resid.bias(tmp$resid, tmp$t0,
                            fit$diag$ctags.t$ttest$p.value,
                            ylab = expression("ctags t residuals"))
            plot.resid.spatial(tmp$resid, tmp$x1, tmp$y1,
                               fit$diag$ctags.t$moransI$p.value)
            ## plot.resid.acf(tmp$resid,
            ##                fit$diag$ctags.t$box$p.value)
            plot.resid.normal(tmp$resid,
                              fit$diag$ctags.t$shapiro$p.value)
        }
    }
    if(fit$conf$use.atags){
        ## x
        resi <- lapply(fit$rep$resid.atags.fine, function(x) x[,1])
        xlims <- c(0, max(sapply(resi,length)))
        ylims <- range(unlist(resi), na.rm = TRUE)
        plot(NA, ty = "n",
             xlim = xlims, ylim = ylims,
             xlab = "Index", ylab = "Residuals")
        for(i in 1:length(resi)){
            points(seq_along(resi[[i]]), resi[[i]],
                   ty = "b", col = adjustcolor("grey60",0.4),
                   pch = NA)
            points(seq_along(resi[[i]]), resi[[i]],
                   col = adjustcolor("grey30",0.4))
        }
        abline(h=0)
        mtext("Archival tags: x", 3, 0.5, font = 2)
        ## Histogram
        hist(unlist(resi), main = "",
             xlab = "Residuals")
        ## QQplot
        qqnorm(unlist(resi), main = "")
        qqline(unlist(resi))
        ## y
        resi <- lapply(fit$rep$resid.atags.fine, function(x) x[,2])
        xlims <- c(0, max(sapply(resi,length)))
        ylims <- range(unlist(resi), na.rm = TRUE)
        plot(NA, ty = "n",
             xlim = xlims, ylim = ylims,
             xlab = "Index", ylab = "Residuals")
        for(i in 1:length(resi)){
            points(seq_along(resi[[i]]), resi[[i]],
                   ty = "b", col = adjustcolor("grey60",0.4),
                   pch = NA)
            points(seq_along(resi[[i]]), resi[[i]],
                   col = adjustcolor("grey30",0.4))
        }
        abline(h=0)
        mtext("Archival tags: y", 3, 0.5, font = 2)
        ## Histogram
        hist(unlist(resi), main = "",
             xlab = "Residuals")
        ## QQplot
        qqnorm(unlist(resi), main = "")
        qqline(unlist(resi))
        ## t
        if(fit$conf$use.effort){
            resi <- fit$rep$resid.atags[,4]
            plot(resi,
                 xlab = "Index", ylab = "Residuals")
            abline(h=0)
            mtext("Archival tags: t", 3, 0.5, font = 2)
            ## Histogram
            hist(resi, 30, main = "",
                 xlab = "Residuals")
            ## QQplot
            qqnorm(resi, main = "")
            qqline(resi)
        }
    }
}


##' Plot taxis2
##' @export
plotmomo.taxis2 <- function(x,
                            cor = 0.1,
                            par = NULL,
                            funcs = NULL,
                            env = NULL,
                            vals = FALSE,
                            digits = 2,
                            plot.numbers = TRUE,
                            plot.cors = TRUE,
                            keep.par = FALSE,
                            ...){

    if(inherits(x, "momo.fit")){
        dat <- x$dat
        fit <- x
    }else{
        dat <- x
        fit <- NULL
    }

    if(!keep.par){
        opar <- par(no.readonly = TRUE)
        on.exit(par(opar))
        if(is.null(par)) par(mfrow = c(1,2)) else par(mfrow = c(3,2))
    }

    if(is.null(par)){
        mati <- matrix(1:2,1,2)
    }else{
        if(plot.numbers && plot.cors){
            mati <- matrix(1:6,3,2, byrow = TRUE)
        }else if(plot.numbers && !plot.cors){
            mati <- matrix(1:4,2,2)
        }else{
            mati <- matrix(1:2,1,2)
        }
    }
    layout(mati)
    par(mar = c(2,2,1,1), oma = c(3,3,2,2))

    if(!is.null(par)){
        par <- get.sim.par(par)
        env <- check.that.list(env)
        dat$env <- env
        dat$env.pred <- NULL
        conf <- def.conf(dat)
        funcs <- get.sim.funcs(funcs, dat, x$conf, env, par)
        hTdx.true <- sapply(dat$time.cont.pred,
                            function(t) apply(dat$xygrid.pred, 1, function(x) funcs$tax(x, t)[1]))
        hTdy.true <- sapply(dat$time.cont.pred,
                            function(t) apply(dat$xygrid.pred, 1, function(x) funcs$tax(x, t)[2]))
        ## uv.true <- t(apply(dat$xygrid, 1, function(xy) funcs$tax(xy, 1)))
        umat.true <- vmat.true <-
            matrix(NA, nrow(dat$celltable), ncol(dat$celltable))
        ## umat.true[cbind(dat$igrid[,1],dat$igrid[,2])] <- uv.true[,1]
        ## vmat.true[cbind(dat$igrid[,1],dat$igrid[,2])] <- uv.true[,2]
        umat.true[cbind(dat$igrid[,1],dat$igrid[,2])] <- rowMeans(hTdx.true)
        vmat.true[cbind(dat$igrid[,1],dat$igrid[,2])] <- rowMeans(hTdy.true)
    }

    if(!is.null(fit)){
        umat.est <- vmat.est <-
            matrix(NA, nrow(dat$celltable), ncol(dat$celltable))
        umat.est[cbind(grid$igrid[,1],grid$igrid[,2])] <- rowMeans(fit$rep$hTdx.pred)
        vmat.est[cbind(grid$igrid[,1],grid$igrid[,2])] <- rowMeans(fit$rep$hTdy.pred)
    }

    if(!is.null(par) && !is.null(fit)){
        umat.re <- (umat.est - umat.true) / umat.true
        vmat.re <- (vmat.est - vmat.true) / vmat.true
    }

    plot.single <- function(mati, dat, vals = TRUE,
                            main = "", cor = 1){
        plot(1, 1, ty = "n",
             xlim = xlims,
             ylim = ylims,
             main = main,
             xaxs = "i", yaxs = "i",
             xlab = "", ylab = "")
        abline(v = dat$xgr, col = "grey70")
        abline(h = dat$ygr, col = "grey70")
        if(vals){
            text(dat$xcen[dat$igrid[,1]],
                 dat$ycen[dat$igrid[,2]],
                 labels = as.numeric(mati))
        }else{
            points(dat$xcen[dat$igrid[,1]],
                   dat$ycen[dat$igrid[,2]],
                   cex = abs(as.numeric(mati)) * cor,
                   lwd = 1.5,
                   col = ifelse(as.numeric(mati) > 0, get.momo.cols(type = "pos"),
                                get.momo.cols(type = "neg")))
        }
        box(lwd = 1.5)
    }

    xlims <- dat$xrange
    ylims <- dat$yrange

    if(!is.null(par) && plot.numbers){
        plot.single(signif(umat.true, digits), dat)
        mtext("u", 3, 1, font = 2)
        plot.single(signif(vmat.true, digits), dat)
        mtext("v", 3, 1, font = 2)
        mtext("True", 4, 1, font = 2)
    }

    if(!is.null(fit) && plot.numbers){
        plot.single(signif(umat.est, digits), dat)
        if(is.null(par)) mtext("u", 3, 1, font = 2)
        plot.single(signif(vmat.est, digits), dat)
        if(is.null(par)) mtext("v", 3, 1, font = 2)
        mtext("Est", 4, 1, font = 2)
    }

    if(!is.null(par) && !is.null(fit) && plot.cors){
        plot.single(signif(umat.re, digits), dat, vals = vals, cor = cor)
        if(!plot.numbers) mtext("u", 3, 1, font = 2)
        plot.single(signif(vmat.re, digits), dat, vals = vals, cor = cor)
        if(!plot.numbers) mtext("v", 3, 1, font = 2)
        if(plot.numbers) mtext("Relative error", 4, 1, font = 2)
    }

    mtext("x", 1, 1, outer = TRUE)
    mtext("y", 2, 1, outer = TRUE)
}


##' Plot residuals2
##' @export
plotmomo.resid2 <- function(fit,
                            keep.par = FALSE,
                            ...){

    if(!keep.par){
        opar <- par(no.readonly = TRUE)
        on.exit(par(opar))
        par(mfrow = c(1,1))
    }

    if(!inherits(fit, "momo.fit")) stop("fit must be a fitted momo object ('momo.fit')!")

    if(!any(names(fit) == "diag")){
        fit <- get.diag(fit)
    }

    if(fit$conf$use.ctags && fit$conf$use.atags){
        if(fit$conf$use.effort){
            layout(matrix(1:30,5,6))
        }else{
            layout(matrix(1:20,5,4))
        }
    }else{
        if(fit$conf$use.effort){
            layout(matrix(1:15,5,3))
        }else{
            layout(matrix(1:10,5,2))
        }
    }
    par(mar = c(4,4,1,1), oma = c(1,1,1,1))

    plot.resid.single <- function(resid, x,
                                  xlab = "",
                                  ylab = "Residuals"){
        plot(x, resid,
             ty = "n",
             xlab = xlab,
             ylab = ylab)
        abline(h = 0, lty = 3)
        points(x, resid, pch = 1)
        mod <- loess(resid ~ x)
        x.new <- seq(min(x), max(x), 0.01)
        lines(x.new, predict(mod, newdata = data.frame(x = x.new)),
              lwd = 3, col = get.momo.cols())
        box(lwd = 1.5)
    }

    if(fit$conf$use.ctags){

        ## x
        tmp <- data.frame(fit$dat$ctags,
                          resid = fit$rep$resid.ctags[,1])
        tmp <- tmp[!is.na(tmp$resid),]
        tmp <- tmp[order(tmp$t1),]
        tmp$tdiff <- tmp$t1 - tmp$t0
        tmp$dist <- sqrt((tmp$x1 - tmp$x0)^2 + (tmp$y1 - tmp$y0)^2)
        plot.resid.single(tmp$resid, tmp$x0,
                          xlab = "Release time t0",
                          ylab = expression("ctags x residuals"))
        ## plot.resid.single(tmp$resid, tmp$x0,
        ##                   xlab = "Recapture time t1",
        ##                   ylab = expression("ctags x residuals"))
        plot.resid.single(tmp$resid, tmp$tdiff,
                          xlab = "Time at sea",
                          ylab = expression("ctags x residuals"))
        plot.resid.single(tmp$resid, tmp$x0,
                          xlab = "Release location x0",
                          ylab = expression("ctags x residuals"))
        plot.resid.single(tmp$resid, tmp$y0,
                          xlab = "Release location y0",
                          ylab = expression("ctags x residuals"))
        ## plot.resid.single(tmp$resid, tmp$x1,
        ##                   xlab = "Recapture location x1",
        ##                   ylab = expression("ctags x residuals"))
        ## plot.resid.single(tmp$resid, tmp$y1,
        ##                   xlab = "Recapture location y1",
        ##                   ylab = expression("ctags x residuals"))
        plot.resid.single(tmp$resid, tmp$dist,
                          xlab = "Distance",
                          ylab = expression("ctags x residuals"))

        ## y
        tmp <- data.frame(fit$dat$ctags,
                          resid = fit$rep$resid.ctags[,2])
        tmp <- tmp[!is.na(tmp$resid),]
        tmp <- tmp[order(tmp$t1),]
        tmp$tdiff <- tmp$t1 - tmp$t0
        tmp$dist <- sqrt((tmp$x1 - tmp$x0)^2 + (tmp$y1 - tmp$y0)^2)
        plot.resid.single(tmp$resid, tmp$x0,
                          xlab = "Release time t0",
                          ylab = expression("ctags y residuals"))
        ## plot.resid.single(tmp$resid, tmp$x0,
        ##                   xlab = "Recapture time t1",
        ##                   ylab = expression("ctags y residuals"))
        plot.resid.single(tmp$resid, tmp$tdiff,
                          xlab = "Time at sea",
                          ylab = expression("ctags y residuals"))
        plot.resid.single(tmp$resid, tmp$x0,
                          xlab = "Release location x0",
                          ylab = expression("ctags y residuals"))
        plot.resid.single(tmp$resid, tmp$y0,
                          xlab = "Release location y0",
                          ylab = expression("ctags y residuals"))
        ## plot.resid.single(tmp$resid, tmp$x1,
        ##                   xlab = "Recapture location x1",
        ##                   ylab = expression("ctags y residuals"))
        ## plot.resid.single(tmp$resid, tmp$y1,
        ##                   xlab = "Recapture location y1",
        ##                   ylab = expression("ctags y residuals"))
        plot.resid.single(tmp$resid, tmp$dist,
                          xlab = "Distance",
                          ylab = expression("ctags y residuals"))

        ## t
        if(fit$conf$use.effort){
            tmp <- data.frame(fit$dat$ctags,
                              resid = fit$rep$resid.ctags[,4])
            tmp <- tmp[!is.na(tmp$resid),]
            tmp <- tmp[order(tmp$t0),]

        }
    }

    if(fit$conf$use.atags){
        ## x
        resi <- lapply(fit$rep$res.axy, function(x) x[,1])
        xlims <- c(0, max(sapply(resi,length)))
        ylims <- range(unlist(resi), na.rm = TRUE)
        plot(NA, ty = "n",
             xlim = xlims, ylim = ylims,
             xlab = "Index", ylab = "Residuals")
        for(i in 1:length(resi)){
            points(seq_along(resi[[i]]), resi[[i]],
                   ty = "b", col = adjustcolor("grey60",0.4),
                   pch = NA)
            points(seq_along(resi[[i]]), resi[[i]],
                   col = adjustcolor("grey30",0.4))
        }
        abline(h=0)
        mtext("Archival tags: x", 3, 0.5, font = 2)
        ## Histogram
        hist(unlist(resi), main = "",
             xlab = "Residuals")
        ## QQplot
        qqnorm(unlist(resi), main = "")
        qqline(unlist(resi))
        ## y
        resi <- lapply(fit$rep$res.axy, function(x) x[,2])
        xlims <- c(0, max(sapply(resi,length)))
        ylims <- range(unlist(resi), na.rm = TRUE)
        plot(NA, ty = "n",
             xlim = xlims, ylim = ylims,
             xlab = "Index", ylab = "Residuals")
        for(i in 1:length(resi)){
            points(seq_along(resi[[i]]), resi[[i]],
                   ty = "b", col = adjustcolor("grey60",0.4),
                   pch = NA)
            points(seq_along(resi[[i]]), resi[[i]],
                   col = adjustcolor("grey30",0.4))
        }
        abline(h=0)
        mtext("Archival tags: y", 3, 0.5, font = 2)
        ## Histogram
        hist(unlist(resi), main = "",
             xlab = "Residuals")
        ## QQplot
        qqnorm(unlist(resi), main = "")
        qqline(unlist(resi))
        ## t
        if(fit$conf$use.effort){
            resi <- fit$rep$res.at
            plot(resi,
                 xlab = "Index", ylab = "Residuals")
            abline(h=0)
            mtext("Archival tags: t", 3, 0.5, font = 2)
            ## Histogram
            hist(resi, 30, main = "",
                 xlab = "Residuals")
            ## QQplot
            qqnorm(resi, main = "")
            qqline(resi)
        }
    }
}


##' Plot simulated data
##' @export
plotmomo.sim <- function(x,
                         keep.par = FALSE,
                         ...){

    if(!keep.par){
        opar <- par(no.readonly = TRUE)
        on.exit(par(opar))
        par(mfrow = c(2,3), mar = c(4,4,1,1), oma = c(1,1,1,1))
    }


    i = 1
    plotmomo.env(x, keep.par = TRUE, main = "")
    add.lab(LETTERS[i])
    i = i + 1
    plotmomo.pref(x, keep.par = TRUE, main = "")
    add.lab(LETTERS[i])
    i = i + 1
    ## plotmomo.pref.spatial(x, keep.par = TRUE, main = "")
    plotmomo.taxis(x, keep.par = TRUE, main = "")
    add.lab(LETTERS[i])
    i = i + 1
    plotmomo.dif(x, keep.par = TRUE, main = "")
    add.lab(LETTERS[i])
    i = i + 1
    if(!is.null(x$ctags)){
        plotmomo.ctags(x$ctags, keep.par = TRUE, main = "")
        add.lab(LETTERS[i])
        i = i + 1
    }
    if(!is.null(x$atags)){
        plotmomo.atags(x$atags, keep.par = TRUE, main = "")
        add.lab(LETTERS[i])
        i = i + 1
    }

}

##' Plot data
##' @export
plotmomo.data <- function(x,
                         keep.par = FALSE,
                         ...){

    if(!keep.par){
        opar <- par(no.readonly = TRUE)
        on.exit(par(opar))
        par(mfrow = c(2,3), mar = c(4,4,1,1), oma = c(1,1,1,1))
    }

    i = 1
    plotmomo.env(x$env, keep.par = TRUE, main = "")
    add.lab(LETTERS[i])
    i = i + 1
    plotmomo.pref(x, keep.par = TRUE, main = "")
    add.lab(LETTERS[i])
    i = i + 1
    ## plotmomo.pref.spatial(x, keep.par = TRUE, main = "")
    plotmomo.taxis(x, keep.par = TRUE, main = "")
    add.lab(LETTERS[i])
    i = i + 1
    plotmomo.dif(x, keep.par = TRUE, main = "")
    add.lab(LETTERS[i])
    i = i + 1
    if(!is.null(x$ctags)){
        plotmomo.ctags(x$ctags, keep.par = TRUE, main = "")
        add.lab(LETTERS[i])
        i = i + 1
    }
    if(!is.null(x$atags)){
        plotmomo.atags(x$atags, keep.par = TRUE, main = "")
        add.lab(LETTERS[i])
        i = i + 1
    }

}


##' Plot to compare single quantity
##'
##' @param fit A result report as generated by running \code{\link{fit.momo}} or
##'     a list containing objects fitted with spict.
##' @param ... Optional additional momo fits.
##' @param quantity Quantities that can be compared. The following options are
#'     currently implemented: \code{"taxis"} for the taxis, \code{"dif"} for the
#'     diffusion, \code{"par"} for the parameters.
##' @param asp Positive number deining the target aspect ratio (columns / rows)
##'     of the plot arrangement.
plotmomo.compare.one <- function(fit, ...,
                                 quantity = c("pref","taxis",
                                              "dif","par"),
                                 keep.par = FALSE,
                                 asp = 2,
                                 col = momo:::get.momo.cols(10),
                                 lty = 1:10,
                                 plot.legend = 1
                                 ){

    if("momo.fit" %in% class(fit) || "momo.sim" %in% class(fit)){
        fitlist <- list(fit, ...)
    }else if(inherits(fit, "list")){
        fitlist <- c(fit, ...)
    }else stop("Please provide fitted momo objects either individually or as list.")

    sim.ind <- lapply(fitlist, function(x) inherits(x, "momo.sim"))

    quantity <- match.arg(quantity)
    n <- length(fitlist)

    if(quantity == "pref"){
        plotmomo.pref(fitlist[[1]], col = col[1],
                       main = "",
                       keep.par = TRUE)
        if(n > 1){
            for(i in 2:n){
                plotmomo.pref(fitlist[[i]], add = TRUE,
                               col = col[i], lty = lty[i],
                               keep.par = TRUE)
            }
        }
    }

    if(quantity == "taxis"){
        plotmomo.taxis(fitlist[[1]], col = col[1],
                       main = "",
                       keep.par = TRUE)
        if(n > 1){
            for(i in 2:n){
                plotmomo.taxis(fitlist[[i]], add = TRUE,
                               col = col[i], lty = lty[i],
                               keep.par = TRUE)
            }
        }
    }

    if(quantity == "dif"){
        plotmomo.dif(fitlist[[1]],
                     col = col[1], lty = lty[1],
                     main = "",
                     keep.par = TRUE)
        if(n > 1){
            for(i in 2:n){
                plotmomo.dif(fitlist[[i]], add = TRUE,
                             col = col[i], lty = lty[i],
                             keep.par = TRUE)
            }
        }
    }

    if(quantity == "par"){

        idx <- which(sapply(fitlist, function(x) inherits(x, "momo.fit")))
        if(length(idx) > 0){
            pars <- unique(unlist(lapply(fitlist[idx], get.par.names)))
        }

        tmp <- lapply(fitlist, function(x){
            if(inherits(x, "momo.sim")){
                nam <- names(x$par.sim)
                map <- names(x$map)[match(nam,names(x$map))]
                map <- map[!is.na(map)]
                mapped <- unlist(x$map[map])
                mapped <- is.na(mapped)
                pars <- unlist(x$par.sim)
                pars <- pars[!names(pars) %in% names(mapped)[mapped]]
                ind <- unlist(sapply(c("beta","logSdObsATS"),
                                     function(x) grep(x, names(pars))))
                if(length(ind) > 0){
                    pars[ind] <- exp(pars[ind])
                }
                lo <- pars
                hi <- pars
            }else if(inherits(x, "momo.fit")){
                ## TODO: make a get.par function!
                nam <- unique(names(x$opt$par))
                map <- names(x$map)[match(nam,names(x$map))]
                map <- map[!is.na(map)]
                mapped <- unlist(x$map[map])
                mapped <- is.na(mapped)
                pars <- unlist(x$pl[nam])
                pars <- pars[!names(pars) %in% names(mapped)[mapped]]
                sds <- unlist(x$plsd[nam])
                sds <- sds[!names(sds) %in% names(mapped)[mapped]]
                lo <- pars - 1.96 * sds
                hi <- pars + 1.96 * sds
                ind <- unlist(sapply(c("beta","logSdObsATS"),
                                     function(x) grep(x, names(pars))))
                if(length(ind) > 0){
                    lo[ind] <- exp(pars[ind] - 1.96 * sds[ind])
                    hi[ind] <- exp(pars[ind] + 1.96 * sds[ind])
                    pars[ind] <- exp(pars[ind])
                }
            }
            return(c(pars,lo,hi))
        })

        ylim <- range(unlist(tmp), na.rm = TRUE)
        xlim <- c(1, max(unique(unlist(lapply(tmp, function(x)
            length(unique(names(x)))))))) + 0.5 * c(-1,1)

        i = 1
        if(inherits(fitlist[[i]], "momo.sim")){

            nam <- names(fitlist[[i]]$par.sim)

            map <- names(fitlist[[i]]$map)[match(nam,names(fitlist[[i]]$map))]
            map <- map[!is.na(map)]
            mapped <- unlist(fitlist[[i]]$map[map])
            mapped <- is.na(mapped)

            pars <- unlist(fitlist[[i]]$par.sim[nam])
            pars <- pars[!names(pars) %in% names(mapped)[mapped]]

            ind <- which(names(pars) %in% c("beta","logSdObsATS"))
            if(length(ind) > 0){
                pars[ind] <- exp(pars[ind])
            }

        }else if(inherits(fitlist[[i]], "momo.fit")){

            ## TODO: make a get.par function!
            nam <- unique(names(fitlist[[i]]$opt$par))

            map <- names(fitlist[[i]]$map)[match(nam,names(fitlist[[i]]$map))]
            map <- map[!is.na(map)]
            mapped <- unlist(fitlist[[i]]$map[map])
            mapped <- is.na(mapped)

            pars <- unlist(fitlist[[i]]$pl[nam])
            pars <- pars[!names(pars) %in% names(mapped)[mapped]]
            sds <- unlist(fitlist[[i]]$plsd[nam])
            sds <- sds[!names(sds) %in% names(mapped)[mapped]]
            lo <- pars - 1.96 * sds
            hi <- pars + 1.96 * sds
            ind <- which(names(pars) %in% c("beta","logSdObsATS"))
            if(length(ind) > 0){
                lo[ind] <- exp(pars[ind] - 1.96 * sds[ind])
                hi[ind] <- exp(pars[ind] + 1.96 * sds[ind])
                pars[ind] <- exp(pars[ind])
            }
        }

        plot(seq(pars), pars,
             ty = "n",
             xlim = xlim,
             xaxt = "n",
             ylim = ylim,
             xlab = "Parameter",
             ylab = "Value")
        axis(1, at = seq(pars), labels = names(pars))

        addi <- seq(-0.1, 0.1, length.out = n)

        if(inherits(fitlist[[i]], "momo.fit")){
            arrows(seq(pars) + addi[i], lo,
                   seq(pars) + addi[i], hi,
                   length = 0.1,
                   angle = 90,
                   code = 3,
                   col = col[i])
        }
        points(seq(pars) + addi[i], pars, col = col[i])

        if(n > 1){
            for(i in 2:n){
                if(inherits(fitlist[[i]], "momo.sim")){
                    pars <- unlist(fitlist[[i]]$par.sim)
                    ind <- which(names(pars) %in% c("beta","logSdObsATS"))
                    if(length(ind) > 0){
                        pars[ind] <- exp(pars[ind])
                    }
                }else if(inherits(fitlist[[i]], "momo.fit")){

                    ## TODO: make a get.par function!
                    nam <- unique(names(fitlist[[i]]$opt$par))

                    map <- names(fitlist[[i]]$map)[match(nam,names(fitlist[[i]]$map))]
                    map <- map[!is.na(map)]
                    mapped <- unlist(fitlist[[i]]$map[map])
                    mapped <- is.na(mapped)

                    pars <- unlist(fitlist[[i]]$pl[nam])
                    pars <- pars[!names(pars) %in% names(mapped)[mapped]]
                    sds <- unlist(fitlist[[i]]$plsd[nam])
                    sds <- sds[!names(sds) %in% names(mapped)[mapped]]
                    lo <- pars - 1.96 * sds
                    hi <- pars + 1.96 * sds
                    ind <- which(names(pars) %in% c("beta","logSdObsATS"))
                    if(length(ind) > 0){
                        lo[ind] <- exp(pars[ind] - 1.96 * sds[ind])
                        hi[ind] <- exp(pars[ind] + 1.96 * sds[ind])
                        pars[ind] <- exp(pars[ind])
                    }
                    arrows(seq(pars) + addi[i], lo,
                           seq(pars) + addi[i], hi,
                           length = 0.1,
                           angle = 90,
                           code = 3,
                           col = col[i])
                }
                points(seq(pars) + addi[i], pars, col = col[i])
            }
        }
        box(lwd = 1.5)
    }
}

##' Plot simulated data
##'
##' @param fit A result report as generated by running \code{\link{fit.momo}} or
##'     a list containing objects fitted with spict.
##' @param ... Optional additional momo fits.
##' @param quantity Quantities that can be compared. The following options are
#'     currently implemented: \code{"taxis"} for the taxis, \code{"dif"} for the
#'     diffusion, \code{"par"} for the parameters.
##' @param asp Positive number deining the target aspect ratio (columns / rows)
##'     of the plot arrangement.
##'
##' @export
plotmomo.compare <- function(fit, ...,
                             quantity = c("pref","taxis",
                                          "dif","par"),
                             keep.par = FALSE,
                             col = momo:::get.momo.cols(10),
                             asp = 2,
                             plot.legend = 1
                             ){

    if("momo.fit" %in% class(fit) || "momo.sim" %in% class(fit)){
        fitlist <- list(fit, ...)
    }else if(inherits(fit, "list")){
        fitlist <- c(fit, ...)
    }else stop("Please provide fitted momo objects either individually or as list.")

    sim.ind <- lapply(fitlist, function(x) inherits(x, "momo.sim"))

    quantity <- match.arg(quantity, several.ok = TRUE)
    nq <- length(quantity)

    if(!keep.par){
        opar <- par(no.readonly = TRUE)
        on.exit(par(opar))
        mfrow <- n2mfrow(nq, asp = asp)
        par(mar = c(4.5,4,1,1)+0.1, oma = c(1,1,1,1))
        if(as.integer(plot.legend) == 1){
            layout(rbind(matrix(1:max(nq,prod(mfrow)),
                                nrow = mfrow[1],
                                ncol = mfrow[2],
                                byrow = TRUE),
                         rep(nq+1, mfrow[2])),
                   heights = c(rep(1, mfrow[1]), 0.15))
        }else{
            layout(matrix(1:nq,
                          nrow = mfrow[1],
                          ncol = mfrow[2],
                          byrow = TRUE))
        }
    }

    for(i in 1:nq){
        plotmomo.compare.one(fitlist,
                             quantity = quantity[i],
                             col = col,
                             plot.legend = as.integer(plot.legend) == 2 && i == nq,
                             keep.par = TRUE)
        add.lab(LETTERS[i])
    }

    if(as.integer(plot.legend) == 1){
        nfit <- sum(!unlist(sim.ind))
        if(is.null(names(fitlist))){
            leg.text <- ifelse(unlist(sim.ind), "Sim",
                               paste0("Fit ",
                                      cumsum(!unlist(sim.ind))))
        }else{
            leg.text <- names(fitlist)
        }

        par(mar = c(1,5,0,0))
        plot.new()
        legend("center", legend = leg.text,
               lwd = 2,
               horiz = TRUE,
               bty = "s",
               box.lwd = 1.5,
               col = col[1:length(fitlist)],
               bg = "white")
    }
}
