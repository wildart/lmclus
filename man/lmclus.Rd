\name{lmclus}
\alias{lmclus}
\alias{lmclusPure}
\alias{lmclus.default}

\title{Linear Manifold Clustering}
\description{
  \code{lmclus} clusters the n-dimensional dataset using linear manifold clustering technique.
}

\usage{
lmclusPure( X, maxDim, numOfClus, noiseSize, bestBound,
            errorBound, maxBinPortion, hisSampling, hisConstSize,
            sampleHeuristic, sampleFactor, randomSeed, showLog )

lmclus(X, \dots)
\method{lmclus}{default}(X, params, \dots)
}

\arguments{
\item{X}{N-dimensional input dataset matrix}
\item{maxDim}{The maximum dimensionality of linear manifold clusters the algorithm will search for, starting with dimensionality 1.}
\item{numOfClus}{The maximum number of clusters we expect to see in the data. }
\item{noiseSize}{Determine the minimum number of points that may be considered as a cluster.}
\item{bestBound}{A histogram 'goodness of separation' threshold, which the higher it is set the more coarse the clustering will be, and the lower it is set the finer the clustering will be.}
\item{errorBound}{A bound on the probability (epsilon) of not sampling points coming from the same cluster using the sampling heuristic.}
\item{maxBinPortion}{Maximum number of points (portion) in a histogram's bin. Does not matter if hisConstSize>0.}
\item{hisSampling}{Sample points to create the distance histograms to improve the efficiency of the algorithm.}
\item{hisConstSize}{Create histograms with a constant number of bins. if value equals zero, histograms of variable sizes will be generated, according to some heuristic. If value greater than zero than that value will be the number of bins.}
\item{sampleHeuristic}{Sampling heuristic mode (1,2,3).}
\item{sampleFactor}{Value will be multiplied by the number of points in the dataset for heuristic 3.}
\item{randomSeed}{Specifies a random generator seed.}
\item{showLog}{Show execution log.}

\item{params}{Parameters of linear manifold clustering, see \link{lmclus.params-class}}

\item{\ldots}{not used}
}

\value{
  \code{lmclusPure} returns a list with three components:
  \item{clusters}{a vector of indexes}
  \item{cluster_dimensions}{a vector of cluster dimensions}
  \item{thresholds}{a vector of a separation thresholds}
 
  \code{lmclus} returns data similar to \code{lmclusPure}, but only accepts params object.
}

\author{Rave Harpaz, Art Diky}

\keyword{manifold}
\keyword{clustering}

\references{"Linear manifold clustering in high dimensional spaces by stochastic search", Pattern Recognition (2007), vol. 40(10), pp 2672-2684.}

\examples{
    ## Create dataset
    ds <- replicate(10, rnorm(100));
    
    # setup parameters
    params <- new("lmclus.params", maxDim = 9, showLog = TRUE)
    print(params)
    
    # cluster dataset
    res <- lmclus( ds, params )
    print(res)
    
    # get cluster 1 data points
    c1_idx <- res$clusters[[1]]
    ds_c1 <- ds[c1_idx, ]
}