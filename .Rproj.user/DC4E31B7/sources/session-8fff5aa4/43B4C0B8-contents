cv.alfagomp <- function( y, x, a = seq(0.1, 1, by = 0.1), sig = 0.05, R = 499,
                         nfolds = 10, folds = NULL, nc = 1, seed = NULL ) {
  if ( min(y) == 0 )  a <- a[a>0]
  lena <- length(a)
  n <- dim(y)[1]
  ina <- 1:n
  if ( is.null(folds) )  folds <- Compositional::makefolds(ina, nfolds = nfolds,
                                                           stratified = FALSE, seed = seed)
  nfolds <- length(folds)

  apa <- proc.time()
  kula <- matrix(nrow = nfolds, ncol = lena)
  for ( i in 1:nfolds ) {
    xtest <- x[ folds[[ i ]], ]
    ytest <- y[ folds[[ i ]], ]
    xtrain <- x[ -folds[[ i ]], ]
    ytrain <- y[-folds[[ i ]], ]
    path <- CompositionalOMP::areg.gomp.path(ytrain, xtrain, a = a, sig = sig, R = R,
                                             model = TRUE, ncores = ncores)$path
    for ( j in 1:lena ) {


      kl <- ytest * log(ytest / yest)
      kl[ is.infinite(kl) ] <- NA
      kula[i, j] <- sum(kl, na.rm = TRUE) / dim(yest)[1]
    }

    kl <- Rfast::colmeans(kula)
    opt <- c( min(kl), a[ which.min(kl) ] )
    names(opt) <- c( "KLD", "alpha")
    apa <- proc.time() - apa


    kula <- kula[, 1:la]
    kl <- Rfast::colmeans(kula)
    opt <- c( min(kl), a[ which.min(kl) ] )
    names(opt) <- c( "KLD", "alpha")
    apa <- proc.time() - apa

  }

  list(runtime = apa, perf = kl, opt = opt)
}
