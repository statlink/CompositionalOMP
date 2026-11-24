cv.diri.gomp <- function(y, x, sig = 0.05, nfolds = 10, folds = NULL, ncores = 1) {

  x <- as.matrix(x)
  y <- as.matrix(y)
  force(sig)
  n <- dim(y)[1]   ;   d <- dim(y)[2]
  if ( is.null(folds) )  folds <- Compositional::makefolds(1:n, nfolds = nfolds, stratified = FALSE)
  nfolds <- length(folds)
  apa <- proc.time()

  kldfun <- function(ytest, est) {
    est <- pmax(est, .Machine$double.eps)
    kl <- ytest * log(ytest / est)
    sum(kl, na.rm = TRUE) / dim(ytest)[1]   ## na.rm removes the 0 * log(0) NaNs
  }

  predfun <- function(res, be, xtest, ytrain, ntest) {
    if ( length(res) == 0 ) {
      est <- Compositional::frechet(ytrain, 1)
      est <- matrix(est, nrow = ntest, ncol = d, byrow = TRUE)
    } else {
      est <- cbind(1, exp( cbind(1, xtest[, res, drop = FALSE]) %*% be ) )
      est <- est / Rfast::rowsums(est)
    }
    est
  }

  runfold <- function(i) {
    xtest <- x[ folds[[ i ]], , drop = FALSE]
    ytest <- y[ folds[[ i ]], , drop = FALSE]
    xtrain <- x[ -folds[[ i ]], , drop = FALSE]
    ytrain <- y[ -folds[[ i ]], , drop = FALSE]
    mod <- CompositionalOMP::diri.gomp(ytrain, xtrain, sig = sig, model = TRUE)
    ntest <- dim(ytest)[1]
    ## p-value
    est <- predfun(mod$pvalue.res[-1, 1], mod$mod$be, xtest, ytrain, ntest)
    kla <- kldfun(ytest, est)
    ## eBIC
    est <- predfun(mod$ebic.res[-1, 1], mod$mod.ebic$be, xtest, ytrain, ntest)
    klb <- kldfun(ytest, est)
    c(kla, klb)
  }

  if ( ncores <= 1 ) {
    kula <- lapply(1:nfolds, runfold)
  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    on.exit( parallel::stopCluster(cl), add = TRUE )
    parallel::clusterExport(cl, varlist = c("x", "y", "sig", "folds", "d", "kldfun", "predfun"),
                            envir = environment())
    kula <- parallel::parLapply(cl, 1:nfolds, runfold)
  }
  kula <- do.call(rbind, kula)
  colnames(kula) <- c("p-value", "eBIC")

  kl <- Rfast::colmeans(kula)
  names(kl) <- c("p-value", "eBIC")
  apa <- proc.time() - apa

  list(runtime = apa, perf = kl)
}
