areg.gomp.path <- function(y, x, a = seq(-1, 1, by = 0.1), xstand = TRUE, sig = 0.05, R = 499,
                           model = FALSE, ncores = 1) {

  runtime <- proc.time()

  if ( xstand )  x <- Rfast::standardise(x)
  if ( min(y) == 0 )  a <- a[a > 0]
  lena <- length(a)

  fun <- function(i) {
    CompositionalOMP::areg.gomp(y, x, a = a[i], xstand = FALSE, sig = sig, R = R, model = model)
  }

  if ( ncores <= 1 ) {
    path <- lapply(1:lena, fun)
  } else {
    ncores <- min(ncores, lena)
    environment(fun) <- globalenv()   ## workers look up everything in their own global env
    cl <- parallel::makeCluster(ncores)
    on.exit( parallel::stopCluster(cl), add = TRUE )
    parallel::clusterExport(cl, c("y", "x", "a", "sig", "R", "model"), envir = environment())
    path <- parallel::parLapply(cl, 1:lena, fun)
  }

  names(path) <- paste("alfa=", a, sep = "")
  runtime <- proc.time() - runtime

  list(runtime = runtime, path = path)
}
