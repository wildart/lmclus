/*
 * =====================================================================================
 *
 *       Filename:  lmclus_r.cpp
 *
 *    Description:  LMCLUS R-interface
 *
 *        Version:  1.0
 *        Created:  02/20/2013 03:21:42 PM
 *       Revision:  none
 *       Compiler:  gcc 4.7.2
 *
 *         Author:  Art Wild (wildart@gmail.com) 
 *        Company:  Patter Recognition Lab, GC CUNY
 *
 * =====================================================================================
 */
#include "lmclus.hpp"

#include <R.h>
#include <Rdefines.h>

extern "C" {
SEXP lmclus(SEXP Xs, SEXP maxDim, SEXP numOfClus, SEXP noiseSize, SEXP bestBound, SEXP errorBound, 
	    SEXP maxBinPortion, SEXP hisSampling, SEXP hisConstSize, SEXP sampleHeuristic, 
	    SEXP sampleFactor, SEXP randomSeed, SEXP showLog) 
{
    Rprintf("Linear manifold clustering...\n"); 
    int n, k, nprotect = 0, show_log;
    size_t i, j;
    
    cpplog::StringLogger log;
    try{   
    
    // Set parameters
    clustering::lmclus::Parameters params;
    params.MAX_DIM = INTEGER(maxDim)[0];
    params.NUM_OF_CLUS = INTEGER(numOfClus)[0];
    params.NOISE_SIZE = static_cast<unsigned int>(INTEGER(noiseSize)[0]);
    params.BEST_BOUND = REAL(bestBound)[0];
    params.BEST_BOUND = REAL(errorBound)[0];
    params.MAX_BIN_PORTION = REAL(maxBinPortion)[0];
    params.HIS_SAMPLING = static_cast<bool>(INTEGER(hisSampling)[0]);
    params.CONST_SIZE_HIS = INTEGER(hisConstSize)[0];
    params.SAMPLING_HEURISTIC = INTEGER(sampleHeuristic)[0];
    params.SAMPLING_FACTOR = REAL(sampleFactor)[0];
    params.RANDOM_SEED = static_cast<unsigned int>(INTEGER(randomSeed)[0]);
    
    show_log = INTEGER(showLog)[0];
    
    // get dataset
    SEXP Rdim = getAttrib(Xs, R_DimSymbol);
    n = INTEGER(Rdim)[0];
    k = INTEGER(Rdim)[1];
    Xs = AS_NUMERIC(Xs);
    arma::mat data(REAL(Xs), n, k, false);
    Rprintf("Data dims: (%d, %d)\n", n, k);
    
    std::vector<arma::uvec> labels;
    std::vector<double> thresholds; 
    std::vector<arma::mat> basises; 
    std::vector<int> clusterDims;
    
    clustering::lmclus::LMCLUS lmclus(&log); 
    lmclus.cluster(data, params, labels, thresholds, basises, clusterDims);
    if (show_log)
        Rprintf("%s", log.getString().c_str());
    
    Rprintf("Clusters found: %d\n", labels.size());
    
    SEXP Return_lst, Rnames, Rthresholds, RclusterDims, Rlabels;
    
     // Thresholds
    PROTECT(Rthresholds = allocVector(REALSXP,thresholds.size())); nprotect++;
    for (i = 0; i < thresholds.size(); ++i)
        REAL(Rthresholds)[i] = thresholds[i];
    
    // Dimensions
    PROTECT(RclusterDims = allocVector(INTSXP,thresholds.size())); nprotect++;
    for (i = 0; i < clusterDims.size(); ++i)
      INTEGER(RclusterDims)[i] = clusterDims[i];
    
    // Labels
    PROTECT(Rlabels = allocVector(VECSXP,labels.size())); nprotect++;
    for (i = 0; i < labels.size(); ++i){
        SEXP lbls;
        PROTECT(lbls = allocVector(INTSXP,labels[i].n_elem)); nprotect++;
        for (j = 0; j < labels[i].n_elem; ++j)
            INTEGER(lbls)[j] = labels[i][j];
        SET_VECTOR_ELT(Rlabels, i, lbls);
    }
    
    // Result list
    PROTECT(Return_lst = allocVector(VECSXP,3)); nprotect++;
    
    /* set names */
    PROTECT(Rnames = NEW_CHARACTER(3)); nprotect++;
    SET_STRING_ELT(Rnames, 0, mkChar("thresholds"));
    SET_STRING_ELT(Rnames, 1, mkChar("cluster_dimensions"));
    SET_STRING_ELT(Rnames, 2, mkChar("clusters"));
    //SET_STRING_ELT(Rnames, 3, mkChar("basis"));
    SET_NAMES(Return_lst, Rnames);

    /* set values */
    SET_VECTOR_ELT(Return_lst, 0, Rthresholds);
    SET_VECTOR_ELT(Return_lst, 1, RclusterDims);
    SET_VECTOR_ELT(Return_lst, 2, Rlabels);
    
    UNPROTECT(nprotect);
    return Return_lst;

    } catch(...) { 
        Rprintf("Log:\n%s", log.getString().c_str());
        ::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return R_NilValue;
}
}