// Dirichlet regression (softmax link for the mean, constant precision phi)
// fitted by Fisher scoring / IRLS for beta and Newton-Raphson for phi.
// C++ port of the R function diri.reg.irls().

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// digamma and trigamma for x > 0: recurrence up to x >= 10, then asymptotic series
// (absolute error ~1e-15). Much faster than R::digamma / R::trigamma.
static inline double fast_digamma(double x) {
  double acc = 0.0;
  while ( x < 10.0 ) { acc -= 1.0 / x;  x += 1.0; }
  const double z = 1.0 / (x * x);
  return acc + std::log(x) - 0.5 / x -
    z * ( 1.0/12 - z * ( 1.0/120 - z * ( 1.0/252 - z * ( 1.0/240 - z * ( 1.0/132 -
    z * ( 691.0/32760 - z / 12.0 ) ) ) ) ) );
}

static inline double fast_trigamma(double x) {
  double acc = 0.0;
  while ( x < 10.0 ) { acc += 1.0 / (x * x);  x += 1.0; }
  const double ix = 1.0 / x, z = ix * ix;
  return acc + ix + 0.5 * z +
    ix * z * ( 1.0/6 - z * ( 1.0/30 - z * ( 1.0/42 - z * ( 1.0/30 - z * ( 5.0/66 -
    z * ( 691.0/2730 - z * 7.0/6 ) ) ) ) ) );
}

static inline void apply_digamma(const arma::mat& a, arma::mat& out) {
  const double* pa = a.memptr();
  double* po = out.memptr();
  const arma::uword N = a.n_elem;
  for ( arma::uword i = 0; i < N; ++i )  po[i] = fast_digamma( pa[i] );
}

static inline void apply_trigamma(const arma::mat& a, arma::mat& out) {
  const double* pa = a.memptr();
  double* po = out.memptr();
  const arma::uword N = a.n_elem;
  for ( arma::uword i = 0; i < N; ++i )  po[i] = fast_trigamma( pa[i] );
}

static inline double sum_lgamma(const arma::mat& a) {
  const double* pa = a.memptr();
  const arma::uword N = a.n_elem;
  double s = 0.0;
  for ( arma::uword i = 0; i < N; ++i )  s += R::lgammafn( pa[i] );
  return s;
}

// x     : n x K design matrix (intercept included)
// ly    : n x D matrix of log(y)
// beta0 : K x (D-1) starting values; if it has 0 rows, ALR least squares are used
// [[Rcpp::export(name = ".diri_reg_irls_cpp")]]
List diri_reg_irls_cpp(const arma::mat& x, const arma::mat& ly, arma::mat beta0,
                       const double phi0 = 1.0, const double tol = 1e-6,
                       const int maxit = 100) {

  const arma::uword n = ly.n_rows, D = ly.n_cols, p = D - 1, K = x.n_cols;
  if ( x.n_rows != n )  stop("x and y must have the same number of rows.");
  if ( D < 2 )  stop("y must have at least two components.");

  arma::mat beta;
  if ( beta0.n_rows == 0 ) {
    // additive log-ratio least squares (first component is the reference)
    arma::mat z = ly.cols(1, p);
    z.each_col() -= ly.col(0);
    bool ok = arma::solve( beta, x.t() * x, x.t() * z, arma::solve_opts::likely_sympd + arma::solve_opts::no_approx );
    if ( !ok )  stop("The design matrix is singular.");
  } else {
    if ( beta0.n_rows != K || beta0.n_cols != p )  stop("Starting values have wrong dimensions.");
    beta = beta0;
  }

  // products x_ik * x_il (k <= l), fixed across iterations: n x K(K+1)/2
  const arma::uword Kq = K * (K + 1) / 2, pq = p * (p + 1) / 2;
  arma::mat XX(n, Kq);
  {
    arma::uword c = 0;
    for ( arma::uword k = 0; k < K; ++k )
      for ( arma::uword l = k; l < K; ++l, ++c )  XX.col(c) = x.col(k) % x.col(l);
  }

  double phi = phi0, loglik = 0.0, loglik_old = -arma::datum::inf;
  int iter = 0;

  arma::mat eta(n, p), mu(n, D), alpha(n, D), dg(n, D), tg(n, D),
            h(n, D), a(n, D), R(n, p), m(n, p), a1(n, p), W(n, pq), M(Kq, pq),
            I_mat(K * p, K * p);
  arma::vec den(n), hbar(n), A(n), S_vec(K * p), step(K * p);

  for ( iter = 1; iter <= maxit; ++iter ) {

    // ---- mean, precision and log-likelihood ----
    eta = x * beta;
    eta = arma::exp(eta);
    den = 1.0 + arma::sum(eta, 1);
    mu.col(0) = 1.0 / den;
    mu.cols(1, p) = eta.each_col() / den;
    alpha = phi * mu;
    loglik = n * R::lgammafn(phi) - sum_lgamma(alpha) + arma::accu( (alpha - 1.0) % ly );

    if ( !std::isfinite(loglik) )  stop("Non-finite log-likelihood.");
    if ( std::abs(loglik - loglik_old) < tol )  break;

    // ---- score for beta ----
    apply_digamma(alpha, dg);
    h = ly - dg;
    hbar = arma::sum(mu % h, 1);
    R = h.cols(1, p);
    R.each_col() -= hbar;
    R %= phi * mu.cols(1, p);
    S_vec = arma::vectorise( x.t() * R );

    // ---- expected information for beta ----
    apply_trigamma(alpha, tg);
    a = tg % arma::square(mu);
    A = arma::sum(a, 1);
    m = mu.cols(1, p);
    a1 = a.cols(1, p);
    {
      arma::uword c = 0;
      for ( arma::uword j = 0; j < p; ++j ) {
        for ( arma::uword l = j; l < p; ++l, ++c ) {
          W.col(c) = m.col(j) % m.col(l) % A - a1.col(j) % m.col(l) - m.col(j) % a1.col(l);
          if ( j == l )  W.col(c) += a1.col(j);
        }
      }
    }
    M = XX.t() * W;     // M(kl, jl') = sum_i x_ik x_il W_i(j, l')
    M *= phi * phi;
    {
      arma::uword c = 0;
      for ( arma::uword j = 0; j < p; ++j ) {
        for ( arma::uword l = j; l < p; ++l, ++c ) {
          arma::uword r = 0;
          for ( arma::uword k1 = 0; k1 < K; ++k1 ) {
            for ( arma::uword k2 = k1; k2 < K; ++k2, ++r ) {
              const double v = M(r, c);
              I_mat(j * K + k1, l * K + k2) = v;  I_mat(j * K + k2, l * K + k1) = v;
              I_mat(l * K + k1, j * K + k2) = v;  I_mat(l * K + k2, j * K + k1) = v;
            }
          }
        }
      }
    }

    bool ok = arma::solve( step, I_mat, S_vec, arma::solve_opts::likely_sympd + arma::solve_opts::no_approx );
    if ( !ok )  stop("The Fisher information matrix is singular.");
    beta += arma::reshape(step, K, p);

    // ---- update phi (direct Newton, no link) ----
    const double S_phi = n * fast_digamma(phi) - arma::accu(mu % dg) + arma::accu(mu % ly);
    const double H_phi = n * fast_trigamma(phi) - arma::accu( arma::square(mu) % tg );
    double phi_new = phi - S_phi / H_phi;
    if ( phi_new <= 0 )  phi_new = 1e-4;
    phi = phi_new;

    loglik_old = loglik;
  }
  if ( iter > maxit )  iter = maxit;

  return List::create( _["iters"] = iter, _["loglik"] = loglik,
                       _["phi"] = phi, _["be"] = beta );
}
