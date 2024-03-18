\name{alpha-RF}
\alias{alfa.rf}
\title{
\eqn{\alpha}-RF
}

\description{
\eqn{\alpha}-RF.
}

\usage{
alfa.rf(xnew, y, x, a = seq(-1, 1, by = 0.1), size = c(1, 2, 3),
depth = c(0, 1), splits = 2:5, R = 500)
}

\arguments{
\item{xnew}{
A matrix with the new compositional data whose group is to be predicted.
Zeros are allowed, but you must be careful to choose strictly positive vcalues of \eqn{\alpha}.
}
\item{y}{
The response variable, it can either be a factor (for classification) or a
numeric vector (for regression). Depending on the nature of the response
variable, the function will proceed with the necessary task.
}
\item{x}{
A matrix with the compositional data.
}
\item{a}{
A vector with a grid of values of the power transformation, it has to be
between -1 and 1. If zero values are present it has to be greater than 0.
If a=0, the isometric log-ratio transformation is applied.
}
\item{size}{
The minimal node size to split at.
}
\item{depth}{
The maximal tree depth. A value of NULL or 0 corresponds to unlimited depth,
1 to tree stumps (1 split per tree).
}
\item{splits}{
The number of random splits to consider for each candidate splitting variable.
}
\item{R}{
The number of trees.
}
}

\details{
For each value of \eqn{\alpha}, the compositional data are transformed and then
the random forest (RF) is applied for one or more combinations of the hyper-parameters.
}

\value{
A list including:
\item{mod}{
A list with the results of the RF model for each value of \eqn{\alpha} that
includes the RF output (a ranger class object) as provided by the function "ranger"
of the package \bold{ranger}, the configurations used and the predicted values of "xnew".
}
}

\references{
Wright M. N. and Ziegler A. (2017). ranger: A fast implementation of random
forests for high dimensional data in C++ and R.
Journal of Statisrical Software 77:1-17. doi:10.18637/jss.v077.i01.

Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation
for compositional data. In Proceedings of the 4th Compositional Data Analysis
Workshop, Girona, Spain. https://arxiv.org/pdf/1106.1451.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{alfarf.tune}
}
}

\examples{
x <- as.matrix(iris[, 1:4])
x <- x/ rowSums(x)
y <- iris[, 5]
mod <- alfa.rf(x, y, x, a = c(0, 0.5, 1), size = 3, depth = 1, splits = 2:3, R = 500)
mod
}

