areg.gomp <- function(y, x, a = 1, sig = 0.05, model = FALSE) {

  runtime <- proc.time()

  mx <- Rfast::colmeans(x)
  sx <- Rfast::colVars(x, std = TRUE)
  x <- Rfast::standardise(x)
  d <- dim(y)[2] - 1
  dm <- dim(x)   ;   n <- dm[1]
  p <- dm[2]     ;   ind <- 1:p

  ya <- Compositional::alfa(y, a)$aff
  sse <- Rfast::eachrow( ya, Rfast::colmeans(ya), oper = "-" )
  sse <- sum(sse^2)
  ela <- Rfast::eachcol.apply(x, ya[, 1])^2
  for ( j in 2:d )  ela <- ela + Rfast::eachcol.apply(x, ya[, j])^2
  sel <- which.max(ela)
  sela <- NULL
  #rho <- c(0, CompositionalSR::alfa.betest(y, x[, sel], a, index = "all", R = R, ya = ya))
  rho <- c(0, CompositionalOMP::areg.ftest(y, x[, sel], a, index = 1, ya = ya)[2] )
  i <- 2
  while ( rho[i] < sig  &  length(sela) < p - 1  &  i < n ) {
    i <- i + 1
    sela <- c(sela, sel)
    ind[sela] <- 0
    mod <- CompositionalSR::areg(y, x[, sela, drop = FALSE], a = a, xnew = x[, sela, drop = FALSE], yb = ya)
    sse[i - 1] <- mod$dev
    res <- ya - Compositional::alfa(mod$est, a)$aff
    r <- numeric(p)
    ii <- ind[ind > 0]
    r[ind] <- Rfast::eachcol.apply(x, res[, 1], indices = ii)^2
    for ( j in 2:d )  r[ind] <- r[ind] + Rfast::eachcol.apply(x, res[, j], indices = ii)^2
    sel <- which.max(r)
    #rho[i] <- CompositionalSR::alfa.betest(y, x[, c(sel, sela)], a, index = 1, R = R, ya = ya)
    rho[i] <- CompositionalOMP::areg.ftest(y, x[, c(sel, sela)], a, index = 1, ya = ya)[2]
  }

  len <- length(sela)
  res <- cbind( c(0, sela), sse[1:(len + 1)], c(NA, rho[seq_len(len) + 1]) )
  colnames(res) <- c("Selected", "SSE", "p-value")

  mod <- NULL
  if ( model &  dim(res)[1] > 1 )  {
    #be[-1, ] <- be[-1, , drop = FALSE] / sx[sela]
    #be[1, ] <- be[1, ] - Rfast::colsums( mx[sela] * be[-1, , drop = FALSE] )
    xsela <- x[, sela, drop = FALSE]
    xsela <- t( t(xsela) * sx[sela]  + mx[sela] )
    mod <- CompositionalSR::areg(y, xsela, a, yb = ya)
  }

  runtime <- proc.time() - runtime

  list(runtime = runtime, res = res, mod = mod)
}



