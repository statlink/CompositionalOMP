\name{dirigomp.ebic.update}
\alias{dirigomp.ebic.update}
\title{
Update the results of the \eqn{gamma}-Orthogonal Matching Pursuit (gamma-OMP) via eBIC
}

\description{
Update the results of the \eqn{\gamma}-OMP via eBIC.
}

\usage{
dirigomp.ebic.update(res, n, p)
}

\arguments{
\item{res}{
The matrix with the results of the \eqn{\gamma}-OMP for the Dirichlet regression.
}
\item{n}{
The sample size.
}
\item{p}{
The number of predictor variables.
}
}

\details{
Suppose you ran \eqn{\gamma}-OMP with Dirichlet regression and now want to
run it again using the eBIC as an extra selection cterion. Instead of running
it again, you may simply update the results using eBIC.
}

\value{
The updated "res"" matrix with the results of the \eqn{\gamma}-OMP with the
Dirichlet regression.
}

\references{
Tsagris M., Papadovasilakis Z., Lakiotaki K. and Tsamardinos I. (2022).
The \eqn{\gamma}-OMP algorithm for feature selection with application to gene
expression data. IEEE/ACM Transactions on Computational Biology and
Bioinformatics, 19(2): 1214--1224.

Chen J. and Chen Z. (2012). Extended BIC for small-n-large-P sparse GLM.
Statistica Sinica, 22(2): 555--574.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{diri.gomp}
}
}

\examples{
y <- as.matrix(iris[, 1:4])
y <- y / rowSums(y)
x <- matrix(rnorm(150 * 30), ncol = 30)
a <- diri.gomp(y, x)
b <- diri.gomp(y, x, ebic = TRUE)
a2 <- dirigomp.ebic.update(a$res, dim(y)[1], dim(x)[2])
}

