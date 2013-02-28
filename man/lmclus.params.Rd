\name{lmclus.params-class}
\docType{class}
\alias{lmclus.params-class}
\title{Class for LMCLUS parameters}
\description{
This class holds parameters for Linear Manifold Clustering method (\link{lmclus}).
}
\section{Creating Objects}{
\code{new(lmclus.params,}\cr
\code{maxDim = ...., # Max dimension}\cr
\code{numberOfClusters = ...., # Number of clusters}\cr
\code{noiseSize = ...., # Noise size}\cr
\code{bestBound = ...., # Best bound}\cr
\code{errorBound = ...., # Error bound}\cr
\code{hisSampling = ...., # Sample points for distance histogram}\cr
\code{hisConstSize = ...., # Histogram bins}\cr
\code{maxBinPortion = ...., # Maximum number of points in a histogram bin}\cr
\code{sampleHeuristic = ...., # Sampling heuristic}\cr
\code{sampleFactor = ...., # Sampling factor}\cr
\code{randomSeed = ...., # RNG seed}\cr
\code{showLog = ...., # Show execution log}\cr
\code{)}}

\section{Slots}{
\describe{
\item{\code{maxDim}:}{The maximum dimensionality of linear manifold clusters the algorithm will search for, starting with dimensionality 1.}
\item{\code{numberOfClusters}:}{The maximum number of clusters we expect to see in the data.}
\item{\code{noiseSize}:}{Determine the minimum number of points that may be considered as a cluster.}
\item{\code{bestBound}:}{A histogram 'goodness of separation' threshold, which the higher it is set the more coarse the clustering will be, and the lower it is set the finer the clustering will be.}
\item{\code{errorBound}:}{A bound on the probability (epsilon) of not sampling points coming from the same cluster using the sampling heuristic.}
\item{\code{hisSampling}:}{Sample points to create the distance histograms to improve the efficiency of the algorithm.}
\item{\code{hisConstSize}:}{Create histograms with a constant number of bins. if value equals zero, histograms of variable sizes will be generated, according to some heuristic. If value greater than zero than that value will be the number of bins.}
\item{\code{maxBinPortion}:}{Maximum number of points (portion) in a histogram's bin. Does not matter if hisConstSize>0.}
\item{\code{sampleHeuristic}:}{Sampling heuristic mode (1,2,3).}
\item{\code{sampleFactor}:}{Value will be multiplied by the number of points in the dataset for heuristic 3.}
\item{\code{randomSeed}:}{Specifies a random generator seed.}
\item{\code{showLog}:}{Show execution log.}
}
}

\section{Prototype}{
\describe{
\item{integer maxDim}{ = 1 }
\item{integer numberOfClusters}{ = 1 }
\item{integer noiseSize}{ = 20 }
\item{double bestBound}{ = 1.0 }
\item{double errorBound}{ = 0.0001 }
\item{logical hisSampling}{ = FALSE }
\item{integer hisConstSize}{ = 0 }
\item{double maxBinPortion}{ = 0.1 }
\item{integer sampleHeuristic}{ = 3 }
\item{double sampleFactor}{ = 0.003 }
\item{integer randomSeed}{ = 0 }
\item{logical showLog}{ = FALSE }
}
}
\keyword{ parameters }
