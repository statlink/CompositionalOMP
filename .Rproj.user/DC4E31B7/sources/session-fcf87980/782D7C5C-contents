\name{areg.ftest}
\alias{areg.ftest}
\title{
F-test for hypothesis testing in the \eqn{alpha}-regression
}

\description{
F-test for hypothesis testing in the \eqn{alpha}-regression.
}

\usage{
areg.ftest(y, x, a, index)
}

\arguments{
\item{y}{
The compositional responses, zero values are allowed.
}
\item{x}{
A matrix with continuous data, where the rows denote the observations and the columns are the variables.
}
\item{a}{
The value of \eqn{\alpha} for the \eqn{\alpha}-transformation.
}
\item{index}{
A number specifying which variable to test.
}
}

\value{
A vector with the F-test statistic and its associated p-value.
}

\references{
Gallant A. R. (1987). Nonlinear Statistical Models. Wiley, New York.

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
\code{ \link{areg.gomp.path}, \link{diri.gomp} }
}

\examples{
y <- as.matrix(iris[, 1:4])
y <- y / rowSums(y)
x <- rnorm(150)
areg.ftest(y, x, a = 1, index = 1)
}

