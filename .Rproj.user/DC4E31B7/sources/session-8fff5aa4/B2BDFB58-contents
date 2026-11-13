boot.gomp <- function(y, x, type = "diri", a = 1, xstand = TRUE, sig = 0.05,
                      B = 500, ncores = 1) {

  runtime <- proc.time()
  n <- dim(y)[1]

  bootfun <- function(i) {
    ina <- rangen::Sample.int(n, n, replace = TRUE)
    yb <- y[ina, ]   ;   xb <- x[ina, ]
    if ( type == "diri" ) {
      CompositionalOMP::diri.gomp(yb, xb, xstand = xstand, sig = sig)$pvalue.res[-1, 1]
    } else if ( type == "areg" ) {
      CompositionalOMP::areg.gomp(yb, xb, a = a, xstand = xstand, sig = sig)$res[-1, 1]
    }
  }

  if ( ncores <= 1 ) {
    sel <- lapply(1:B, bootfun)
  } else {
    ncores <- min(ncores, B)
    environment(bootfun) <- globalenv()   ## workers look up everything in their own global env
    cl <- parallel::makeCluster(ncores)
    on.exit( parallel::stopCluster(cl), add = TRUE )
    parallel::clusterExport(cl, c("y", "x", "n", "type", "a", "xstand", "sig", "R"),
                            envir = environment())
    sel <- parallel::parLapply(cl, 1:B, bootfun)
  }

  sel <- unlist(sel)
  tab <- Rfast::Table(sel)
  res <- cbind( as.numeric( names(tab) ), tab / B )
  colnames(res) <- c("Variable", "Selection proportion")
  rownames(res) <- NULL
  runtime <- proc.time() - runtime
  list(runtime = runtime, res = res)
}
