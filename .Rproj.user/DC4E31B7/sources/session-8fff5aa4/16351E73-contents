cv.areg.gomp <- function( y, x, a = seq(0.1, 1, by = 0.1), sig = 0.05,
                         nfolds = 10, folds = NULL) {
  if ( min(y) == 0 )  a <- a[a>0]
  lena <- length(a)
  n <- dim(y)[1]
  ina <- 1:n
  if ( is.null(folds) )  folds <- Compositional::makefolds(ina, nfolds = nfolds, stratified = FALSE)
  nfolds <- length(folds)
  apa <- proc.time()
  kula <- kulb <- matrix(nrow = nfolds, ncol = lena)
  d <- dim(y)[2]

  for ( i in 1:nfolds ) {
    xtest <- x[ folds[[ i ]], ]
    ytest <- y[ folds[[ i ]], ]
    xtrain <- x[ -folds[[ i ]], ]
    ytrain <- y[-folds[[ i ]], ]
    path <- CompositionalOMP::areg.gomp.path(ytrain, xtrain, a = a, sig = sig, model = TRUE)$path
    ntest <- dim(ytest)[1]
    for ( j in 1:lena ) {
      mod <- path[[ j ]]
      ### p-value
      res <- mod$pvalue.res[-1, 1]
      if ( length(res) == 0 ) {
        est <- Compositional::frechet(ytrain, a[j])
        est <- matrix(est, nrow = ntest, ncol = d, byrow = TRUE)
      } else {
        est <- cbind(1, exp( cbind(1, xtest[, res]) %*% mod$be ) )
        est <- est / Rfast::rowsums(est)
      }
      kl <- ytest * log(ytest / est)
      kl[ is.infinite(kl) ] <- NA
      kula[i, j] <- sum(kl, na.rm = TRUE) / ntest
      ## BIC
      res <- mod$bic.res[-1, 1]
      if ( length(res) == 0 ) {
        est <- Compositional::frechet(ytrain, a[j])
        est <- matrix(est, nrow = ntest, ncol = d, byrow = TRUE)
      } else {
        est <- cbind(1, exp( cbind(1, xtest[, res]) %*% mod$be.bic ) )
        est <- est / Rfast::rowsums(est)
      }
      kl <- ytest * log(ytest / est)
      kl[ is.infinite(kl) ] <- NA
      kulb[i, j] <- sum(kl, na.rm = TRUE) / ntest
    }
  }

  kl <- rbind( Rfast::colmeans(kula), Rfast::colmeans(kulb) )
  colnames(kl) <- paste("alfa=", a, sep = "")
  rownames(kl) <- c("p-value", "BIC")
  opt <- rbind( c( min(kl[1, ]), a[ which.min(kl[1, ]) ] ), c( min(kl[2, ]), a[ which.min(kl[2, ]) ] ) )
  colnames(opt) <- c( "KLD", "alpha")
  rownames(opt) <- c("p-value", "BIC")
  apa <- proc.time() - apa

  list(runtime = apa, perf = kl, opt = opt)
}
