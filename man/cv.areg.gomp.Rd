\name{cv.areg.gomp}
\alias{cv.areg.gomp}
\title{
Cross-validation for \eqn{gamma}-OMP algorithm for the \eqn{alpha}-regression
}

\description{
Cross-validation for \eqn{gamma}-OMP algorithm for the \eqn{alpha}-regression.
}

\usage{
cv.areg.gomp(y, x, a = seq(0.1, 1, by = 0.1), sig = 0.05, R = 499,
nfolds = 10, folds = NULL, seed = NULL, ncores = 1)
}

\arguments{
\item{y}{
The compositional responses, zero values are allowed.
}
\item{x}{
A matrix with continuous data, where the rows denote the observations and the
columns are the variables.
}
\item{a}{
A vector with \eqn{\alpha} values for the \eqn{\alpha}-transformation.
}
\item{sig}{
The significance level to terminate the algorithm, default value is 0.05.
}
\item{R}{
How many permutations to perform to test the significance of the candidate selected variable.
}
\item{nfolds}{
The number of folds to create.
}
\item{folds}{
If you aleady have the list with the folds pass it here.
}
\item{seed}{
If you want the folds to be the same type a seed number here.
}
\item{ncores}{
If you want parallel computations specify the number of cores to use.
}
}

\value{
A list including
\item{runtime}{
The runtime of the algorithm.
}
\item{perf}{
The estimated performance, predictive Kullback-Leibler divergence,
for each value of \eqn{\alpha}.
}
\item{opt}{
A vector with the minimum Kullback-Leibler divergence and the corresponding
value of \eqn{\alpha}.
}
}

\references{
Pati Y. C., Rezaiifar R. and Krishnaprasad P. S. (1993). Orthogonal matching pursuit: Recursive function approximation with applications to wavelet decomposition. In Signals, Systems and Computers. 1993 Conference Record of The Twenty-Seventh Asilomar Conference on. IEEE.

Davis G. (1994). Adaptive Nonlinear Approximations. PhD thesis.
http://www.geoffdavis.net/papers/dissertation.pdf

Mallat S. G. and Zhang Z. (1993). Matching pursuits with time-frequency dictionaries. IEEE Transactions on signal processing, 41(12), 3397--3415.
https://www.di.ens.fr/~mallat/papiers/MallatPursuit93.pdf

Tsagris M., Papadovasilakis Z., Lakiotaki K. and Tsamardinos I. (2018). Efficient feature selection on gene expression data: Which algorithm to use? BioRxiv.

Tsagris M., Papadovasilakis Z., Lakiotaki K. and Tsamardinos I. (2022).
The \eqn{\gamma}-OMP algorithm for feature selection with application to gene expression data.
IEEE/ACM Transactions on Computational Biology and Bioinformatics, 19(2): 1214--1224.

Tsagris M. and Pantazis Y. (2026). The \eqn{\alpha}--regression for compositional data:
a unified framework for standard, temporal and spatial regression models including compositional
predictors. https://arxiv.org/pdf/2510.12663

Tsagris M. (2015). Regression analysis with compositional data containing zero values.
Chilean Journal of Statistics, 6(2): 47-57. https://arxiv.org/pdf/1508.01913v1.pdf

Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for
compositional data. In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
https://arxiv.org/pdf/1106.1451.pdf
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{areg.gomp}, \link{areg.gomp.path} }
}

\examples{
y <- as.matrix(iris[, 1:4])
y <- y / rowSums(y)
x <- matrix(rnorm(150 * 10), ncol = 10)
}

