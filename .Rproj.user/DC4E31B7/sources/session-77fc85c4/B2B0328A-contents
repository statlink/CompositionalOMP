areg.gomp <- function(y, x, a = 1, sig = 0.05, model = FALSE) {

  runtime <- proc.time()

  mx <- Rfast::colmeans(x)
  sx <- Rfast::colVars(x, std = TRUE)
  x <- Rfast::standardise(x)
  D <- dim(y)[2]   ;   d <- D - 1
  dm <- dim(x)   ;   n <- dm[1]
  p <- dm[2]     ;   ind <- 1:p
  nd <- n * d
  H <- Compositional::helm(D)   ## the same d x D Helmert sub-matrix used by alfa()

  ## ---- Gallant's (1987) F-test, two-stage estimator ----
  ## Sigma from the SSE fit of the full model; then both models are re-fitted by minimising
  ## S(B) = sum_i e_i' Sigma^{-1} e_i, with nls.lm on the whitened residuals vec(E U'),
  ## where Sigma^{-1} = U'U.

  ## fitted values in alpha-space and nu_i = softmax(a * eta_i)
  afit <- function(be, X) {
    mu <- cbind( 1, exp(X %*% matrix(be, ncol = d)) )
    mu <- mu / Rfast::rowsums(mu)
    if ( abs(a) <= 1e-6 ) {
      lmu <- log(mu)
      z <- lmu - Rfast::rowmeans(lmu)
      v <- matrix(1 / D, n, D)
    } else {
      v <- mu^a
      v <- v / Rfast::rowsums(v)
      z <- (D * v - 1) / a
    }
    list(m = z %*% t(H), v = v)
  }
  resfun <- function(be, X, yw, tU)  as.vector( yw - afit(be, X)$m %*% tU )
  ## block (j', k) of the Jacobian = -diag( (w_k U')_{.j'} ) X,  w_ijk = D nu_ik (h_jk - hbar_ij)
  jacfun <- function(be, X, yw, tU) {
    v <- afit(be, X)$v
    vh <- v %*% t(H)
    k <- dim(X)[2]
    J <- matrix(0, nd, k * d)
    for ( l in 2:D ) {
      awl <- ( D * v[, l] * ( matrix(H[, l], n, d, byrow = TRUE) - vh ) ) %*% tU
      for ( j in 1:d )  J[(j - 1) * n + 1:n, (l - 2) * k + 1:k] <- - awl[, j] * X
    }
    J
  }
  wfit <- function(be, X, yw, tU) {
    mod <- minpack.lm::nls.lm( par = as.vector(be), fn = resfun, jac = jacfun,
                               X = X, yw = yw, tU = tU,
                               control = minpack.lm::nls.lm.control(maxiter = 1000) )
    list( be = matrix(mod$par, ncol = d), value = sum(mod$fvec^2) )
  }
  ## xs: predictors of the full model, the tested predictor is its LAST column
  ## m1, e1: SSE fit of the full model and its alpha-residuals
  ftest <- function(xs, m1, e1) {
    tU <- t( chol( solve( crossprod(e1) / n ) ) )
    yw <- ya %*% tU
    X1 <- cbind(1, xs)
    k1 <- dim(X1)[2]
    f1 <- wfit(m1$be, X1, yw, tU)                   ## full model, warm start: SSE estimates
    if ( k1 == 2 ) {                                ## restricted model is intercept only
      e0 <- Rfast::eachrow( ya, Rfast::colmeans(ya), oper = "-" )
      S0 <- sum( (e0 %*% tU)^2 )
    } else {                                        ## restricted model, same Sigma
      b0 <- f1$be[-k1, , drop = FALSE]
      S0 <- wfit(b0, X1[, -k1, drop = FALSE], yw, tU)$value
    }
    p1 <- k1 * d   ;   r <- d
    stat <- ( (S0 - f1$value) / r ) / ( f1$value / (nd - p1) )
    pf(stat, r, nd - p1, lower.tail = FALSE)
  }

  powcl <- function(z) {
    z <- z^a
    z / Rfast::rowsums(z)
  }

  ## derivative-based residuals of the alpha-regression (n x d matrix):
  dres <- function(e, v) {
    g <- e %*% H                                   ## n x D
    rt <- D * v * ( g - Rfast::rowsums(v * g) )    ## D * nu_ik * (g_ik - gbar_i)
    rt[, -1, drop = FALSE]                         ## the first component is the reference
  }

  crit <- function(rt, ii)   Rfast::colsums( crossprod(rt, x[, ii, drop = FALSE])^2 )

  fit.test <- function(xs) {
    m1 <- tryCatch( CompositionalSR::areg(y, xs, a = a, xnew = xs, yb = ya),
                    error = function(e) NULL )
    if ( is.null(m1) )  return( list(pval = NA, mod = NULL, e = NULL, v = NULL) )
    e1 <- ya - Compositional::alfa(m1$est, a)$aff
    pval <- tryCatch( ftest(xs, m1, e1), error = function(e) NA )
    list(pval = pval, mod = m1, e = e1, v = powcl(m1$est))
  }

  ya <- Compositional::alfa(y, a)$aff
  e0 <- Rfast::eachrow( ya, Rfast::colmeans(ya), oper = "-" )  ## intercept-only residuals
  p0 <- d
  sse <- sum(e0^2)
  phi <- sse / (nd - d)
  v0 <- matrix( Rfast::colmeans( powcl(y) ), n, D, byrow = TRUE )
  ela <- crit( dres(e0, v0), ind )
  sel <- which.max(ela)
  sela <- NULL

  ft <- fit.test(x[, sel, drop = FALSE])
  test <- ft$pval
  if ( !is.finite(test) )  test <- 1
  rho <- c(0, test)

  i <- 2
  while ( rho[i] < sig  &  length(sela) < p - 1  &  i < n ) {
    i <- i + 1
    sela <- c(sela, sel)
    ind[sela] <- 0
    xs <- x[, sela, drop = FALSE]

    ## the alternative model of the previous test is the current model, no refitting
    mod <- ft$mod
    if ( is.null(mod) || !is.finite(mod$dev) ) {
      sela <- sela[-length(sela)]
      rho[i] <- 1
      break
    }
    sse[i - 1] <- mod$dev
    res <- ft$e                 ## alpha-residuals of the current model
    rt <- dres(res, ft$v)       ## derivative-based residuals, used for the selection
    p0 <- length(mod$be)
    ## number of parameters and dispersion of the current model, computed along the path
    phi[i - 1] <- sse[i - 1] / (nd - p0)

    r <- numeric(p)
    ii <- ind[ind > 0]
    r[ii] <- crit(rt, ii)
    cand <- ii[ order(r[ii], decreasing = TRUE) ]
    cand <- cand[ r[cand] > 0 ]
    qrs <- qr(xs)
    sel <- NA
    for ( k in cand ) {
      e <- qr.resid( qrs, x[, k] )
      if ( sum(e^2) / (n - 1) > 1e-6 ) {
        sel <- k   ;   break
      }
      ind[k] <- 0   ## collinear with the selected set now, so also later (the span only grows)
    }

    if ( is.na(sel) ) {
      rho[i] <- 1
    } else {
      ## columns ordered as c(sela, sel), so that if sel is accepted this fit is the next model
      ft <- fit.test(x[, c(sela, sel), drop = FALSE])
      test <- ft$pval
      if ( is.finite(test) ) {
        rho[i] <- test
      } else {
        rho[i] <- 1
      }
    }
  }

  len <- length(sela)
  sse <- sse[1:(len + 1)]
  phi <- phi[1:(len + 1)]

  dbic <- NULL
  len.bic <- 0
  if ( len > 0 ) {
    k <- 1:len
    dbic <- d * log(n) - ( sse[k] - sse[k + 1] ) / phi[k + 1]
    stp <- which( !(dbic < 0) )
    len.bic <- if ( length(stp) == 0 )  len  else  stp[1] - 1
  }
  sela.bic <- sela[ seq_len(len.bic) ]

  pvalue.res <- cbind( c(0, sela), sse, c(NA, rho[seq_len(len) + 1]) )
  colnames(pvalue.res) <- c("Selected", "SSE", "p-value")
  kb <- 1:(len.bic + 1)
  bic.res <- cbind( c(0, sela)[kb], sse[kb], phi[kb], c(NA, dbic)[kb] )
  colnames(bic.res) <- c("Selected", "SSE", "phi", "BIC change")

  be <- NULL
  be.bic <- NULL
  if ( model  &  len > 0 )  {
    be <- CompositionalSR::areg(y, x[, sela, drop = FALSE], a, yb = ya)$be
    be[-1, ] <- be[-1, , drop = FALSE] / sx[sela]
    be[1, ] <- be[1, ] - Rfast::colsums( mx[sela] * be[-1, , drop = FALSE] )
    if ( len.bic == len ) {
      be.bic <- be
    } else if ( len.bic > 0 ) {
      be.bic <- CompositionalSR::areg(y, x[, sela.bic, drop = FALSE], a, yb = ya)$be
      be.bic[-1, ] <- be.bic[-1, , drop = FALSE] / sx[sela.bic]
      be.bic[1, ] <- be.bic[1, ] - Rfast::colsums( mx[sela.bic] * be.bic[-1, , drop = FALSE] )
    }
  }

  runtime <- proc.time() - runtime

  list(runtime = runtime, pvalue.res = pvalue.res, be = be, bic.res = bic.res, be.bic = be.bic)
}