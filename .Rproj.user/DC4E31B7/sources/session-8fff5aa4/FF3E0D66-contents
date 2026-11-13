diri.gomp <- function(y, x, xstand = TRUE, sig = 0.05, model = FALSE) {

  runtime <- proc.time()

  if ( xstand )  x <- Rfast::standardise(x)
  d <- dim(y)[2] - 1
  dm <- dim(x)   ;   n <- dm[1]
  p <- dm[2]     ;   ind <- 1:p
  ly <- log(y)

  mod <- Compositional::diri.nr(y)
  lik <- mod$loglik
  ela <- Rfast::eachcol.apply(x, y[, 2])^2
  for ( j in 3:d )  ela <- ela + Rfast::eachcol.apply(x, y[, j])^2
  sel <- which.max(ela)
  sela <- NULL
  mod <- Compositional::diri.reg.irls(y, x[, sel], xnew = x[, sel])
  lik[2] <- mod$loglik
  stat <- 2 * ( lik[2] - lik[1] )
  pvalue <- pchisq(stat, d, lower.tail = FALSE)
  rho <- c(0, pvalue)

  i <- 2
  while ( rho[i] < sig  &  length(sela) < p - 1  &  i < n ) {
    i <- i + 1
    sela <- c(sela, sel)
    ind[sela] <- 0

    mu <- mod$est
    psi1_phi <- trigamma(mod$phi)
    alpha <- mod$phi * mu
    h <- ly - digamma(alpha)
    hbar <- Rfast::rowsums(mu * h)
    res <- mod$phi * mu[, -1] * (h[, -1] - hbar)

    r <- numeric(p)
    ii <- ind[ind > 0]
    r[ind] <- Rfast::eachcol.apply(x, res[, 1], indices = ii)^2
    for ( j in 2:d )  r[ind] <- r[ind] + Rfast::eachcol.apply(x, res[, j], indices = ii)^2
    sel <- which.max(r)
    mod <- Compositional::diri.reg.irls(y, x[, c(sela, sel)], xnew = x[, c(sela, sel)])
    lik[i] <- mod$loglik
    stat <- 2 * ( lik[i] - lik[i - 1] )
    rho[i] <- pchisq(stat, d, lower.tail = FALSE)
  }

  len <- length(sela)
  pvalue.res <- cbind( c(0, sela), lik[1:(len + 1)], c(NA, rho[seq_len(len) + 1]) )
  colnames(pvalue.res) <- c("Selected", "log-lik", "p-value")

  res <- pvalue.res   ;   ebic.res <- NULL
  if ( dim(pvalue.res)[1] > 1 ) {
    sel <- 1:c(len + 1)
    if ( n > p )  gam <- 0 else gam <- 0.5
    for ( i in 1:len ) {
      ebic <- 2 * log(n) + 2 * gam * (p - i) / (i + 1)
      if ( 2 * (res[i + 1, 2] - res[i, 2] ) < ebic ) {
        sel[i + 1] <- 0
      }
    }
    ebic.res <- res[sel, , drop = FALSE]
  }

  mod <- NULL
  if ( model  &  dim(res)[1] > 1 )  mod <- Compositional::diri.reg.irls(y, x[, sela, drop = FALSE])

  runtime <- proc.time() - runtime

  list(runtime = runtime, pvalue.res = pvalue.res, ebic.res = ebic.res, mod = mod)
}



