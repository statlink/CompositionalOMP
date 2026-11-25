diri.reg.irls <- function(y, x, xnew = NULL, tol = 1e-6, maxit = 100) {

  runtime <- proc.time()
  y <- as.matrix(y)
  D <- dim(y)[2]
  x <- model.matrix( ~., data = as.data.frame(x) )

  mod <- .diri_reg_irls_cpp( x, log(y), matrix(0, 0, 0), 1.0, tol, maxit )
  beta <- mod$be
  colnames(beta) <- paste0("Y", 2:D)
  rownames(beta) <- colnames(x)

  est <- NULL
  if ( !is.null(xnew) ) {
    xnew <- model.matrix(~., data.frame(xnew) )
    mu <- cbind( 1, exp(xnew %*% beta) )
    est <- mu / Rfast::rowsums(mu)
  }

  runtime <- proc.time() - runtime

  list( runtime = runtime, iters = mod$iters, loglik = mod$loglik, phi = mod$phi, be = beta, est = est)
}
