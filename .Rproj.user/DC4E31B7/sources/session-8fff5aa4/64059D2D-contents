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

  ftest <- function(e0, e1, p0, p1) {
    stat <- n * sum( diag( solve( crossprod(e1), crossprod(e0) ) ) )
    r <- p1 - p0
    stat <- (nd - p1) / nd * (stat - nd) / r
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

  fit.test <- function(xs, e0, p0) {
    m1 <- tryCatch( CompositionalSR::areg(y, xs, a = a, xnew = xs, yb = ya),
                    error = function(e) NULL )
    if ( is.null(m1) )  return( list(pval = NA, mod = NULL, e = NULL, v = NULL) )
    e1 <- ya - Compositional::alfa(m1$est, a)$aff
    pval <- tryCatch( ftest(e0, e1, p0, length(m1$be)), error = function(e) NA )
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

  ft <- fit.test(x[, sel, drop = FALSE], e0, p0)
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
    res <- ft$e                 ## alpha-residuals, used by the F-test
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
      ft <- fit.test(x[, c(sela, sel), drop = FALSE], res, p0)
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
