\name{Tuning the parameters of the alpha-RF}
\alias{alfarf.tune}
\title{
Tuning the parameters of the \eqn{\alpha}-RF
}

\description{
Tuning the parameters of the \eqn{\alpha}-RF.
}

\usage{
alfarf.tune(y, x, a = seq(-1, 1, by = 0.1), size = c(1, 2, 3),
depth = c(0, 1), splits = 2:5, R = 500, ncores = 1, folds = NULL,
nfolds = 10, stratified = TRUE, seed = NULL, graph = FALSE)
}

\arguments{
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
\item{ncores}{
The number of cores to use. If more than 1, parallel computing will take place.
It is advisable to use it if you have many observations and or many variables,
otherwise it will slow down the process.
}
\item{folds}{
If you have the list with the folds supply it here. You can also leave it
NULL and it will create folds.
}
\item{nfolds}{
The number of folds in the cross validation.
}
\item{stratified}{
Do you want the folds to be created in a stratified way? TRUE or FALSE.
}
\item{seed}{
You can specify your own seed number here or leave it NULL.
}
\item{graph}{
If graph is TRUE (default value) a plot will appear.
}
}

\details{
K-fold cross-validation of the \eqn{\alpha}-RF with compositional predictor
variables is performed to select the optimal value of \eqn{\alpha} and the optimal
configutrations of the random forest (RF).
}

\value{
If graph is true, a graph with the estimated performance for each value of \eqn{\alpha}.
A list including:
\item{per}{
A vector with the estimated performance for each value of \eqn{\alpha}.
}
\item{performance}{
A vector with the optimal performance and the optimal combinations of the
hyper-parameters of the RF.
}
\item{best_a}{
The value of \eqn{\alpha} corresponding to the optimal performance.
}
\item{runtime}{
The time required by the cross-validation procedure.
}
}

\references{
Wright M. N. and Ziegler A. (2017). ranger: A fast implementation of random
forests for high dimensional data in C++ and R.
Journal of Statisrical Software 77:1-17. doi:10.18637/jss.v077.i01.

Friedman Jerome, Trevor Hastie and Robert Tibshirani (2009).
The elements of statistical learning, 2nd edition. Springer, Berlin.

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
\code{\link{alfa.rf}
}
}

\examples{
x <- as.matrix(iris[, 1:4])
x <- x/ rowSums(x)
y <- iris[, 5]
mod <- alfa.rf(x, y, x, a = c(0, 0.5, 1), size = 3, depth = 1, splits = 2:3, R = 500)
mod
}

