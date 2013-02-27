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
\item{\code{maxDim}:}{integer}
\item{\code{numberOfClusters}:}{integer}
\item{\code{noiseSize}:}{integer}
\item{\code{bestBound}:}{double}
\item{\code{errorBound}:}{double}
\item{\code{hisSampling}:}{logical}
\item{\code{hisConstSize}:}{integer}
\item{\code{maxBinPortion}:}{double}
\item{\code{sampleHeuristic}:}{integer}
\item{\code{sampleFactor}:}{double}
\item{\code{randomSeed}:}{integer}
\item{\code{showLog}:}{logical}
}
}

\section{Prototype}{
\describe{
\item{integer maxDim}{ = 1 }
\item{integer numberOfClusters}{ = 1 }
\item{integer noiseSize}{ = 20 }
\item{integer bestBound}{ = 1.0 }
\item{integer errorBound}{ = 0.0001 }
\item{integer hisSampling}{ = FALSE }
\item{integer hisConstSize}{ = 0 }
\item{integer maxBinPortion}{ = 0.1 }
\item{integer sampleHeuristic}{ = 3 }
\item{integer sampleFactor}{ = 0.003 }
\item{integer randomSeed}{ = 0 }
\item{integer showLog}{ = FALSE }
}
}
\keyword{ parameters }
