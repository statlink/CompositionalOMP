cv.diri.gomp <- function(y, x, sig = 0.05, nfolds = 10, folds = NULL, seed = NULL) {

  n <- dim(y)[1]   ;   p <- dim(x)[2]
  ina <- 1:n
  if ( is.null(folds) )  folds <- Compositional::makefolds(ina, nfolds = nfolds,
                                                           stratified = FALSE, seed = seed)
  nfolds <- length(folds)
  apa <- proc.time()
  kula <- matrix(nrow = nfolds, ncol = 2)
  colnames(kula) <- c("p-value", "eBIC")

  for ( i in 1:nfolds ) {
    xtest <- x[ folds[[ i ]], ]
    ytest <- y[ folds[[ i ]], ]
    xtrain <- x[ -folds[[ i ]], ]
    ytrain <- y[-folds[[ i ]], ]
    mod <- CompositionalOMP::diri.gomp(ytrain, xtrain, sig = sig, model = TRUE)
    ntest <- dim(ytest)[1]
    res <- mod$pvalue.res[-1, 1]
    if ( length(res) == 0 ) {
      est <- Compositional::frechet(ytrain, 1)
      est <- matrix( rep(est, ntest), nrow = ntest, byrow = TRUE)
    } else {
      est <- cbind(1, exp( cbind(1, xtest[, res]) %*% mod$mod$be ) )
      est <- est / Rfast::rowsums(est)
    }
    kl <- ytest * log(ytest / est)
    kl[ is.infinite(kl) ] <- NA
    kula[i, 1] <- sum(kl, na.rm = TRUE) / dim(est)[1]
    ### eBIC now
    res <- mod$ebic.res[-1, 1]
    if ( length(res) == 0 ) {
      est <- Compositional::frechet(ytrain, 1)
      est <- matrix( rep(est, ntest), nrow = ntest, byrow = TRUE)
    } else  est <- Compositional::diri.reg.irls(ytrain, xtrain[, res], xnew = xtest[, res])$est
    kl <- ytest * log(ytest / est)
    kl[ is.infinite(kl) ] <- NA
    kula[i, 2] <- sum(kl, na.rm = TRUE) / dim(est)[1]
  }

  kl <- Rfast::colmeans(kula)
  names(kl) <- c("p-value", "eBIC")
  apa <- proc.time() - apa

  list(runtime = apa, perf = kl)
}
