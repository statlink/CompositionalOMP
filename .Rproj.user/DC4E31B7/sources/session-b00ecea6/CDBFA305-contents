areg.ftest <- function(y, x, a, index, ya = NULL) {
  ## Gallant's (1987) F-test for multivariate non-linear regression, alpha-regression.
  ## Stage 1: SSE fit of the full model -> Sigma.
  ## Stage 2: both models re-fitted by minimising S(B) = sum_i e_i' Sigma^{-1} e_i
  ##          with minpack.lm::nls.lm on the whitened residuals e_i' U', where Sigma^{-1} = U'U.
  x <- as.matrix(x)
  if ( is.null(ya) )  ya <- Compositional::alfa(y, a)$aff
  n <- dim(ya)[1]   ;   d <- dim(ya)[2]   ;   D <- d + 1
  h <- Compositional::helm(D)

  ## stage 1
  m1 <- CompositionalSR::areg(y, x, a, yb = ya, xnew = x)
  e1 <- ya - Compositional::alfa(m1$est, a)$aff
  sigma <- crossprod(e1) / n
  tU <- t( chol( solve(sigma) ) )   ## Sigma^{-1} = U'U, whitened residuals: E %*% t(U)
  yw <- ya %*% tU

  ## fitted values in alpha-space and nu_i = softmax(a * eta_i)
  fitted <- function(be, X) {
    mu <- cbind( 1, exp(X %*% matrix(be, ncol = d)) )
    mu <- mu / Rfast::rowsums(mu)
    if ( a == 0 ) {
      lmu <- log(mu)
      z <- lmu - Rfast::rowmeans(lmu)
      v <- matrix(1 / D, n, D)
    } else {
      v <- mu^a
      v <- v / Rfast::rowsums(v)
      z <- (D * v - 1) / a
    }
    list(m = z %*% t(h), v = v)
  }
  ## whitened residuals, vec(E U')
  resfun <- function(be, X)  as.vector( yw - fitted(be, X)$m %*% tU )
  ## Jacobian: block (j', k) = -diag( (w_k U')_{.j'} ) X,
  ## with w_k the n x d matrix w_ijk = D nu_ik (h_jk - hbar_ij), k = 2, ..., D
  jacfun <- function(be, X) {
    v <- fitted(be, X)$v
    vh <- v %*% t(h)
    k <- dim(X)[2]
    J <- matrix(0, n * d, k * d)
    for ( l in 2:D ) {
      wl <- D * v[, l] * ( matrix(h[, l], n, d, byrow = TRUE) - vh )
      awl <- wl %*% tU
      for ( j in 1:d )  J[(j - 1) * n + 1:n, (l - 2) * k + 1:k] <- - awl[, j] * X
    }
    J
  }
  wfit <- function(be, X) {
    mod <- minpack.lm::nls.lm( par = as.vector(be), fn = resfun, jac = jacfun, X = X,
                               control = minpack.lm::nls.lm.control(maxiter = 1000) )
    if ( !mod$info %in% 1:3 )  warning("nls.lm: ", mod$message)
    list( be = matrix(mod$par, ncol = d), value = sum(mod$fvec^2) )
  }

  ## stage 2: full model, warm start at the SSE estimates
  X1 <- cbind(1, x)
  b1 <- as.matrix(m1$be)
  f1 <- wfit(b1, X1)
  p1 <- length(b1)

  ## stage 2: restricted model, same Sigma
  x0 <- x[, -index, drop = FALSE]
  if ( dim(x0)[2] >= 1 ) {
    X0 <- cbind(1, x0)
    b0 <- f1$be[-(index + 1), , drop = FALSE]   ## warm start
    S0 <- wfit(b0, X0)$value
    p0 <- length(b0)
  } else {   ## intercept only: fitted value is the sample mean
    e0 <- Rfast::eachrow(ya, Rfast::colmeans(ya), oper = "-")
    S0 <- sum( (e0 %*% tU)^2 )
    p0 <- d
  }

  r <- p1 - p0
  S1 <- f1$value
  stat <- ( (S0 - S1) / r ) / ( S1 / (n * d - p1) )
  pvalue <- pf(stat, r, n * d - p1, lower.tail = FALSE)
  res <- c(stat, pvalue)
  names(res) <- c("statistic", "p-value")
  res
}
