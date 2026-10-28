areg.gomp.path <- function(y, x, a = seq(-1, 1, by = 0.1), xstand = TRUE, sig = 0.05,
                           model = FALSE) {

  runtime <- proc.time()

  if ( xstand )  x <- Rfast::standardise(x)
  if ( min(y) == 0 )  a <- a[a > 0]
  lena <- length(a)

  fun <- function(i) {
    CompositionalOMP::areg.gomp(y, x, a = a[i], xstand = FALSE, sig = sig, model = model)
  }

  path <- lapply(1:lena, fun)
  names(path) <- paste("alfa=", a, sep = "")
  runtime <- proc.time() - runtime

  list(runtime = runtime, path = path)
}
