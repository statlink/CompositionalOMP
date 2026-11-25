boot.gomp <- function(y, x, type = "areg", a = 1, sig = 0.05, B = 500, ncores = 1) {

  runtime <- proc.time()
  n <- dim(y)[1]

  bootfun <- function(i) {
    ina <- rangen::Sample.int(n, n, replace = TRUE)
    yb <- y[ina, , drop = FALSE]   ;   xb <- x[ina, , drop = FALSE]
    if ( type == "diri" ) {
      mod <- CompositionalOMP::diri.gomp(yb, xb, sig = sig)
      pvalue <- mod$pvalue.res[-1, 1]  ;  ic <- mod$ebic.res[-1, 1]
    } else if ( type == "areg" ) {
      mod <- CompositionalOMP::areg.gomp(yb, xb, a = a, sig = sig)
      pvalue <- mod$pvalue.res[-1, 1]  ;  ic <- mod$bic.res[-1, 1]
    }
    list(pvalue = pvalue, ic = ic)
  }

  if ( ncores <= 1 ) {
    sel <- lapply(1:B, bootfun)
  } else {
    ncores <- min(ncores, B)
    environment(bootfun) <- globalenv()   ## workers look up everything in their own global env
    cl <- parallel::makeCluster(ncores)
    on.exit( parallel::stopCluster(cl), add = TRUE )
    parallel::clusterExport(cl, c("y", "x", "n", "type", "a", "sig"), envir = environment())
    sel <- parallel::parLapply(cl, 1:B, bootfun)
  }

  proptab <- function(s) {
    s <- unlist(s)
    if ( length(s) == 0 ) {
      res <- matrix(nrow = 0, ncol = 2)
    } else {
      tab <- Rfast::Table(s)
      res <- cbind( as.numeric( names(tab) ), as.numeric(tab) / B )
    }
    colnames(res) <- c("Variable", "Selection proportion")
    rownames(res) <- NULL
    res
  }

  res.pvalue <- proptab( lapply(sel, `[[`, "pvalue") )
  res.ic <- proptab( lapply(sel, `[[`, "ic") )

  runtime <- proc.time() - runtime
  list(runtime = runtime, pvalue = res.pvalue, ic = res.ic)
}
