\name{Tuning the parameters of the alfa-PPR}
\alias{alfappr.tune}
\title{
Tuning the parameters of the\eqn{\alpha}-PPR
}

\description{
Tuning the parameters of the\eqn{\alpha}-PPR.
}

\usage{
alfappr.tune(y, x, a = seq(-1, 1, by = 0.1), nterms = 1:10, ncores = 1,
folds = NULL, nfolds = 10, seed = NULL, graph = FALSE)
}

\arguments{
\item{y}{
The response variable, a numerical vector.
}
\item{x}{
A matrix with the compositional data.
}
\item{a}{
A vector with a grid of values of the power transformation, it has to be
between -1 and 1. If zero values are present it has to be greater than 0.
If a=0, the isometric log-ratio transformation is applied.
}
\item{nterms}{
The number of terms to include in the model.
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
\item{seed}{
You can specify your own seed number here or leave it NULL.
}
\item{graph}{
If graph is TRUE (default value) a plot will appear.
}
}

\details{
K-fold cross-validation of the \eqn{\alpha}-PPR with compositional predictor
variables is performed to select the optimal value of \eqn{\alpha} and the numer of
terms in the PPR.
}

\value{
If graph is true, a graph with the estimated performance for each value of \eqn{\alpha}.
A list including:
\item{per}{
A vector with the estimated performance for each value of \eqn{\alpha}.
}
\item{performance}{
A vector with the optimal performance and the optimal number of terms.
}
\item{best_a}{
The value of \eqn{\alpha} corresponding to the optimal performance.
}
\item{runtime}{
The time required by the cross-validation procedure.
}
}

\references{
Friedman J. H. and Stuetzle W. (1981). Projection pursuit regression.
Journal of the American Statistical Association, 76, 817-823. doi: 10.2307/2287576.

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
\code{\link{alfa.ppr}
}
}

\examples{
x <- as.matrix(iris[, 1:3])
x <- x/ rowSums(x)
y <- iris[, 4]
mod <- alfappr.tune(y, x, a = c(0, 0.5, 1), nterms = c(2, 3))
mod
}

