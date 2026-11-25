areg.ftest <- function(y, x, a, index) {
  x <- as.matrix(x)
  ya <- Compositional::alfa(y, a)$aff
  n <- dim(ya)[1]   ;   d <- dim(ya)[2]
  m1 <- CompositionalSR::areg(y, x, a, yb = ya, xnew = x)
  e1 <- ya - Compositional::alfa(m1$est, a)$aff
  sigma <- crossprod(e1) / n
  p1 <- length(m1$be)

  if ( dim(x[, -index, drop = FALSE])[2] >= 1 ) {
    m0 <- CompositionalSR::areg(y, x[, -index], a, yb = ya, xnew = x[, -index])
    p0 <- length(m0$be)
    e0 <- ya - Compositional::alfa(m0$est, a)$aff
  } else {
    e0 <- t( t(ya) - Rfast::colmeans(ya, a) )
    p0 <- d
  }

  r <- p1 - p0
  stat <- sum( Rfast::mahala(e0, numeric(d), sigma) )
  stat <- (n * d - p1) / (n * d) * (stat - n * d) / r
  pvalue <- pf(stat, r, n * d - p1, lower.tail = FALSE)
  res <- c(stat, pvalue)
  names(res) <- c("statistic", "p-value")
  res
}
