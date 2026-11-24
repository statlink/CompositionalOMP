cv.areg.gomp <- function( y, x, a = seq(0.1, 1, by = 0.1), sig = 0.05,
                          nfolds = 10, folds = NULL, ncores = 1) {
  if ( min(y) == 0 )  a <- a[a > 0]
  lena <- length(a)
  n <- dim(y)[1]
  d <- dim(y)[2]
  if ( is.null(folds) )  folds <- Compositional::makefolds(1:n, nfolds = nfolds, stratified = FALSE)
  nfolds <- length(folds)
  apa <- proc.time()

  kldfun <- function(ytest, est) {
    kl <- ytest * log(ytest / est)
    kl[ is.infinite(kl) ] <- NA
    sum(kl, na.rm = TRUE) / dim(ytest)[1]
  }

  runfold <- function(i) {
    xtest <- x[ folds[[ i ]], , drop = FALSE]
    ytest <- y[ folds[[ i ]], , drop = FALSE]
    xtrain <- x[ -folds[[ i ]], , drop = FALSE]
    ytrain <- y[ -folds[[ i ]], , drop = FALSE]
    path <- CompositionalOMP::areg.gomp.path(ytrain, xtrain, a = a, sig = sig, model = TRUE)$path
    ntest <- dim(ytest)[1]
    kla <- klb <- numeric(lena)
    for ( j in 1:lena ) {
      mod <- path[[ j ]]
      ## p-value
      res <- mod$pvalue.res[-1, 1]
      if ( length(res) == 0 ) {
        est <- matrix( Compositional::frechet(ytrain, a[j]), nrow = ntest, ncol = d, byrow = TRUE )
      } else {
        est <- cbind(1, exp( cbind(1, xtest[, res, drop = FALSE]) %*% mod$be ) )
        est <- est / Rfast::rowsums(est)
      }
      kla[j] <- kldfun(ytest, est)
      ## BIC
      res <- mod$bic.res[-1, 1]
      if ( length(res) == 0 ) {
        est <- matrix( Compositional::frechet(ytrain, a[j]), nrow = ntest, ncol = d, byrow = TRUE )
      } else {
        est <- cbind(1, exp( cbind(1, xtest[, res, drop = FALSE]) %*% mod$be.bic ) )
        est <- est / Rfast::rowsums(est)
      }
      klb[j] <- kldfun(ytest, est)
    }
    c(kla, klb)
  }

  if ( ncores <= 1 ) {
    kk <- lapply(1:nfolds, runfold)
  } else {
    cl <- parallel::makePSOCKcluster(ncores)
    on.exit( parallel::stopCluster(cl), add = TRUE )
    parallel::clusterExport(cl, varlist = c("x", "y", "a", "sig", "folds", "lena", "d", "kldfun"),
                            envir = environment())
    kk <- parallel::parLapply(cl, 1:nfolds, runfold)
  }
  kk <- do.call(rbind, kk)
  kula <- kk[, 1:lena, drop = FALSE]
  kulb <- kk[, -(1:lena), drop = FALSE]

  kl <- rbind( Rfast::colmeans(kula), Rfast::colmeans(kulb) )
  colnames(kl) <- paste("alfa=", a, sep = "")
  rownames(kl) <- c("p-value", "BIC")
  opt <- rbind( c( min(kl[1, ]), a[ which.min(kl[1, ]) ] ),
                c( min(kl[2, ]), a[ which.min(kl[2, ]) ] ) )
  colnames(opt) <- c("KLD", "alpha")
  rownames(opt) <- c("p-value", "BIC")
  apa <- proc.time() - apa
  list(runtime = apa, perf = kl, opt = opt)
}
