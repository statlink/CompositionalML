\name{alfa-PPR with compositional predictor variables}
\alias{alfa.ppr}
\title{
\eqn{\alpha}-PPR with compositional predictor variables
}

\description{
\eqn{\alpha}-PPR with compositional predictor variables.
}

\usage{
alfa.ppr(xnew, y, x, a = seq(-1, 1, by = 0.1), nterms = 1:10)
}

\arguments{
\item{xnew}{
A matrix with the new compositional data whose group is to be predicted.
Zeros are allowed, but you must be careful to choose strictly positive vcalues of \eqn{\alpha}.
}
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
}

\details{
This is the standard projection pursuit regression (PPR) applied to the
\eqn{\alpha}-transformed compositional predictors.
See the built-in function "ppr" for more details.
}

\value{
A list including:
\item{mod}{
A list with the results of the PPR model for each value of \eqn{\alpha} that
includes the PPR output as provided by the function "ppr", for each value of "nterms".
}
\item{est}{
A list with the predicted response values of "xnew" for each value of \eqn{\alpha}
and number of "nterms".
}
}

\references{
Friedman J. H. and Stuetzle W. (1981). Projection pursuit regression.
Journal of the American Statistical Association, 76, 817-823. doi: 10.2307/2287576.

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
\code{\link{alfappr.tune}
}
}

\examples{
x <- as.matrix(iris[, 1:3])
x <- x/ rowSums(x)
y <- iris[, 4]
mod <- alfa.ppr(x, y, x, a = c(0, 0.5, 1), nterms = c(2, 3))
mod
}

