## LMCLUS is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## See <http://www.gnu.org/licenses/>.

setClass(
    Class="lmclus.params",
    representation=representation(
        maxDim = "numeric", 
        numberOfClusters = "numeric", 
        noiseSize  = "numeric", 
        bestBound  = "numeric", 
        errorBound = "numeric", 
        
        hisSampling = "logical",
        hisConstSize = "numeric", 
        maxBinPortion = "numeric",
        
        sampleHeuristic = "numeric", 
        sampleFactor  = "numeric", 
        randomSeed = "numeric",
        
        showLog = "logical"
    ),
    prototype=prototype(
        
        maxDim = 1, 
        numberOfClusters = 1, 
        noiseSize  = 20, 
        bestBound  = 1.0, 
        errorBound = 0.0001, 

        hisSampling = FALSE,
        hisConstSize = 0,
        maxBinPortion = 0.1,
        
        sampleHeuristic = 3, 
        sampleFactor  = 0.003, 
        randomSeed = 0,
        
        showLog = FALSE
    )
)

setMethod ("show", "lmclus.params",
function (object)
{
    cat("Linear Manifold Clustering Parameters \n")
    
    cat("Max dimension:", object@maxDim, "\n")
    cat("Number of clusters:", object@numberOfClusters, "\n")
    cat("Noise size:" , object@noiseSize, "\n")
    cat("Best bound:" , object@bestBound, "\n")
    cat("Error bound:" , object@errorBound , "\n")
    
    cat("Sample points for distance histogram:" , object@hisSampling , "\n")
    cat("Histogram bins:" , object@hisConstSize , "\n")
    cat("Maximum number of points in a histogram's bin:" , object@maxBinPortion , "\n")
    
    cat("Sampling heuristic:" , object@sampleHeuristic , "\n")
    cat("Sampling factor:" , object@sampleFactor , "\n")
    cat("Random seed:" , object@randomSeed , "\n")
    
    cat("Show log:" , object@showLog , "\n")
    cat("\n")
})

lmclusPure <- function(X, maxDim, numOfClus, noiseSize, bestBound, errorBound, maxBinPortion,
                       hisSampling, hisConstSize, sampleHeuristic, sampleFactor, randomSeed, showLog) 
{
    .Call("lmclus", X, maxDim, numOfClus, noiseSize, bestBound, errorBound, maxBinPortion, 
          hisSampling, hisConstSize, sampleHeuristic, sampleFactor, randomSeed, showLog,
          package = "lmclus")
}

lmclus <- function(X, ...) UseMethod("lmclus")

lmclus.default <- function(X, params, ...) 
{
    stopifnot(isClass(params), class(params)[1] == "lmclus.params")
    
    X <- as.matrix(X)
    
    maxDim <- as.integer(params@maxDim)
    numOfClus <- as.integer(params@numberOfClusters)
    noiseSize <- as.integer(params@noiseSize)
    bestBound <- as.double(params@bestBound)
    errorBound <- as.double(params@errorBound)
    maxBinPortion <- as.double(params@maxBinPortion)
    hisSampling <- as.integer(params@hisSampling)
    hisConstSize <- as.integer(params@hisConstSize)
    sampleHeuristic <- as.integer(params@sampleHeuristic)
    sampleFactor <- as.double(params@sampleFactor)
    randomSeed <- as.integer(params@randomSeed)
    showLog <- as.integer(params@showLog)
    
    res <- lmclusPure(X, maxDim, numOfClus, noiseSize, bestBound, errorBound, maxBinPortion,
                      hisSampling, hisConstSize, sampleHeuristic, sampleFactor, randomSeed, showLog)
    
    res$call <- match.call()
    
    class(res) <- "lmclus"
    res
}

print.lmclus <- function(x, ...) {
    cat("\nCall:\n")
    print(x$call)
    cat("\nClusters:\n")
    print(x$clusters)
    cat("\nCluster dimensions:\n")
    print(x$cluster_dimensions)
}