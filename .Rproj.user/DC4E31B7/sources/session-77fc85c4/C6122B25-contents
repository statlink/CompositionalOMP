\name{diri.reg.irls}
\alias{diri.reg.irls}
\title{
Dirichlet regression via Fisher scoring (IRLS), implemented in C++
}

\description{
Dirichlet regression with a softmax (additive log-ratio) link for the mean and a constant precision parameter, fitted by Fisher scoring (IRLS) for the regression coefficients and Newton-Raphson for the precision.
}

\usage{
diri.reg.irls(y, x, xnew = NULL, tol = 1e-6, maxit = 100)
}

\arguments{
\item{y}{
A matrix with the compositional data (response). Zero values are not allowed.
}
\item{x}{
A matrix or a data.frame with the predictor variables.
}
\item{xnew}{
If you have new data use it, otherwise leave it NULL.
}
\item{tol}{
The tolerance value to terminate the algorithm; the change in the log-likelihood between successive iterations.
}
\item{maxit}{
The maximum number of iterations.
}
}

\details{
The mean of the Dirichlet distribution is linked to the predictors through \eqn{\mu_{ij} = e^{\bm{x}_i^\top\bm{\beta}_j} / (1 + \sum_{k=2}^D e^{\bm{x}_i^\top\bm{\beta}_k})}, with the first component as the reference (\eqn{\bm{\beta}_1 = \bm{0}}), and the precision \eqn{\phi} is common to all observations. At every iteration, the coefficients are updated with a Fisher scoring (IRLS) step using the expected information matrix, followed by a Newton-Raphson step for \eqn{\phi}. Starting values for the coefficients are the least squares estimates of the additive log-ratio transformed data, and \eqn{\phi} starts at 1. The whole algorithm is implemented in C++.
}

\value{
A list including:
\item{runtime}{
The time required by the regression.
}
\item{iters}{
The number of iterations required.
}
\item{loglik}{
The value of the log-likelihood.
}
\item{phi}{
The estimated precision parameter.
}
\item{be}{
The estimated regression coefficients, a matrix with \eqn{D-1} columns, one for each of the components 2 to \eqn{D}.
}
\item{est}{
The fitted or the predicted values (if xnew is not NULL).
}
}

\references{
Maier M. J. (2014). DirichletReg: Dirichlet Regression for Compositional Data in R. Research Report Series/Department of Statistics and Mathematics, 125. WU Vienna University of Economics and Business, Vienna.

Tsagris M. and Stewart C. (2018). A Dirichlet regression model for compositional data with zeros. Lobachevskii Journal of Mathematics, 39(3): 398--412.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

\seealso{
\code{\link{diri.gomp}, \link{cv.diri.gomp}}
}

\examples{
x <- as.vector(iris[, 4])
y <- as.matrix(iris[, 1:3])
y <- y / rowSums(y)
mod <- diri.reg.irls(y, x)
}
