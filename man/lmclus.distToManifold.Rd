\name{lmclus.distToManifold}
\alias{lmclus.distToManifold}
\alias{distToManifoldPure}

\title{
    Calculates distance from point to manifold
}
\description{
    \code{lmclus.distToManifold} calculates distance from point to manifold defined by basis and origin.
}
\usage{
lmclus.distToManifold(point, basis, origin)
distToManifoldPure(point, basis_t)
}

\arguments{
  \item{point}{point in N-dimensional space}
  \item{basis}{a basis matrix}
  \item{basis_t}{a transpose of basis matrix}
  \item{origin}{an origin of manifold or translation vector for manifold from Euclidian origin point}
}

\details{
    \code{distToManifoldPure} accepts point vector with value translated to origin and transposed basis matrix.
    \code{lmclus.distToManifold} before calculating distance would translate point to origin if specified.
}
\value{
    Returns scalar distance value.
}

\author{Art Diky}

\examples{
    basis <- matrix(c(1,0,0,0,1,0,0,0,1), 3, 3)
    point <- c(2,2,2)
    origin <- c(0,0,0)
    lmclus.distToManifold(point, basis, origin)
}

\keyword{manifold}
\keyword{distance}
\keyword{linear algebra}
