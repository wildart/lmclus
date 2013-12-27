\name{lmclus.kittler}
\alias{lmclus.kittler}

\title{
Kittler threasholding
}
\description{
Performs Kittler threasholding on collection
}
\usage{
    lmclus.kittler(v,bins)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v}{vector of scalar values}
  \item{bins}{number of bins in histogram}
}
\value{
  \item{threshold}{Histogram separation threshold}
  \item{separation}{Separation value}
}
\author{Art Diky}

\examples{
    dists <- sample(1:10,100,replace=T)
    res <- lmclus.kittler(dist)
    show(res)
}

\keyword{ threshold }
