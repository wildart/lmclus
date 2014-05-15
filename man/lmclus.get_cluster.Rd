\name{lmclus.get_cluster}
\alias{lmclus.get_cluster}
\title{Get cluster from LMCLUS result}
\description{Extracts particular cluster from result of \code{lmclus} call.}
\usage{lmclus.get_cluster(clusters, i)}

\arguments{
  \item{clusters}{Clusters as a result of \code{lmclus} function call}
  \item{i}{Cluster index}
}
\value{
    Returns list of following values:
    \item{id}{Cluster id}
    \item{dim}{Cluster manifold dimensionality}
    \item{origin}{Cluster manifold origin point}
    \item{basis}{Cluster manifold basis vectors (as rows in matrix)}, 
    \item{labels}{Dataset point indexes that are in this cluster}
    \item{threshold}{Cluster separation threshold}
}
\author{Art Diky}

\examples{
cluster10 <- lmclus.get_cluster(clusters, 10)
cluster10$basis
}
