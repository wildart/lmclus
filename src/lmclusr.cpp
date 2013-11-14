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

void output_callback(const char *msg){
    Rprintf("%s\n", msg);
}

extern "C" {
SEXP lmclus(SEXP Xs, SEXP maxDim, SEXP numOfClus, SEXP noiseSize, SEXP bestBound, SEXP errorBound, 
	    SEXP maxBinPortion, SEXP hisSampling, SEXP hisConstSize, SEXP sampleHeuristic, 
	    SEXP sampleFactor, SEXP randomSeed, SEXP showLog) 
{
    Rprintf("Linear manifold clustering...\n"); 
    int n, m, nprotect = 0, show_log;
    size_t i, j, k;
    
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
    m = INTEGER(Rdim)[1];
    Xs = AS_NUMERIC(Xs);
    arma::mat data(REAL(Xs), n, m, false);
    Rprintf("Data dims: (%d, %d)\n", n, m);

    if(params.MAX_DIM >= m){ 
        Rprintf("Linear manifold dimension must be less than the dimension of the data !!!\n");
        return R_NilValue;
    }
    
    std::vector<arma::uvec> labels;
    std::vector<double> thresholds; 
    std::vector<arma::mat> bases; 
    std::vector<int> clusterDims;
    
    clustering::lmclus::LMCLUS lmclus(&log); 
    lmclus.cluster(data, params, labels, thresholds, bases, clusterDims, output_callback);
    if (show_log)
        Rprintf("%s", log.getString().c_str());
    
    Rprintf("Clusters found: %d\n", labels.size());
    
    SEXP Return_lst, Rnames, Rthresholds, RclusterDims, Rlabels, Rbases;
    
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
    
    // bases
    PROTECT(Rbases = allocVector(VECSXP,bases.size())); nprotect++;
    for (i = 0; i < bases.size(); ++i){
        int r = bases[i].n_rows, c = bases[i].n_cols;
        SEXP bss, dimnames;
        PROTECT(bss = allocMatrix(REALSXP, r, c)); nprotect++;
        //PROTECT(dimnames = allocVector(VECSXP, 2)); nprotect++;
        //SET_VECTOR_ELT(dimnames, 0, r);
        //SET_VECTOR_ELT(dimnames, 1, c);
        //setAttrib(bss, R_DimNamesSymbol, dimnames);
        for (j = 0; j < r; ++j)
            for (k = 0; k < c; ++k)
                REAL(bss)[j+r*k] = bases[i].at(j, k);
        SET_VECTOR_ELT(Rbases, i, bss);
    }
    
    // Result list
    PROTECT(Return_lst = allocVector(VECSXP,4)); nprotect++;
    
    /* set names */
    PROTECT(Rnames = NEW_CHARACTER(4)); nprotect++;
    SET_STRING_ELT(Rnames, 0, mkChar("thresholds"));
    SET_STRING_ELT(Rnames, 1, mkChar("cluster_dimensions"));
    SET_STRING_ELT(Rnames, 2, mkChar("clusters"));
    SET_STRING_ELT(Rnames, 3, mkChar("bases"));
    SET_NAMES(Return_lst, Rnames);

    /* set values */
    SET_VECTOR_ELT(Return_lst, 0, Rthresholds);
    SET_VECTOR_ELT(Return_lst, 1, RclusterDims);
    SET_VECTOR_ELT(Return_lst, 2, Rlabels);
    SET_VECTOR_ELT(Return_lst, 3, Rbases);
    
    UNPROTECT(nprotect);
    return Return_lst;

    } catch(...) { 
        Rprintf("Log:\n%s", log.getString().c_str());
        ::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return R_NilValue;
}
}
