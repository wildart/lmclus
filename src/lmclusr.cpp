/* Copyright 2005-2014 Pattern Recognition Laboratory, The City University of New York
 * =====================================================================================
 *
 *       Filename:  lmclusr.cpp
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
#include <armadillo>
#include "lmclus.hpp"

#include <R.h>
#include <Rdefines.h>


void output_callback(const char *msg) {
    Rprintf("%s\n", msg);
}

extern "C" {

SEXP distToManifold(SEXP x, SEXP b_t) {
    int l, n, m;
    try {
    // Get point
    l = LENGTH(x);
    x = AS_NUMERIC(x);
    arma::rowvec point(REAL(x), l, false);
    // Get basis
    SEXP Rdim = getAttrib(b_t, R_DimSymbol);
    n = INTEGER(Rdim)[0];
    m = INTEGER(Rdim)[1];
    b_t = AS_NUMERIC(b_t);
    arma::mat B_T(REAL(b_t), n, m, false);
    // Calculate distance
    double dist = clustering::lmclus::LMCLUS::distanceToManifold(point, B_T);
    return ScalarReal(dist);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    return R_NilValue;
}

SEXP kittler(SEXP nH, SEXP minX, SEXP maxX) {
    unsigned int i, n, nprotect = 0;
    try {
    n = LENGTH(nH);
    nH = AS_NUMERIC(nH);
    arma::vec hist(REAL(nH), n, false);
    double RHmin = REAL(minX)[0];
    double RHmax = REAL(maxX)[0];

    Kittler K;
    bool res = K.FindThreshold(hist, RHmin, RHmax);
    if (!res)
        Rprintf("no minimum, unimode histogram\n");

    SEXP Return_lst, Rnames, Rthreshold, Rdiscriminability, Rdepth,
         Rcriterion, Rglobalmin, Rminidx;
    PROTECT(Rthreshold = ScalarReal(K.GetThreshold())); nprotect++;
    PROTECT(Rdiscriminability = ScalarReal(K.GetDiscrim())); nprotect++;
    PROTECT(Rdepth = ScalarReal(K.GetDepth())); nprotect++;
    PROTECT(Rglobalmin = ScalarReal(K.GetGlobalMin())); nprotect++;
    PROTECT(Rminidx = ScalarInteger(K.GetMinIndex())); nprotect++;
    // Criterion function
    auto cf = K.GetCriterionFunc();
    PROTECT(Rcriterion = allocVector(REALSXP, cf.size())); nprotect++;
    for (i = 0; i < cf.size(); ++i)
        REAL(Rcriterion)[i] = cf[i];

    // Result
    PROTECT(Return_lst = allocVector(VECSXP, 6)); nprotect++;

    /* set names */
    PROTECT(Rnames = NEW_CHARACTER(6)); nprotect++;
    SET_STRING_ELT(Rnames, 0, mkChar("threshold"));
    SET_STRING_ELT(Rnames, 1, mkChar("discriminability"));
    SET_STRING_ELT(Rnames, 2, mkChar("depth"));
    SET_STRING_ELT(Rnames, 3, mkChar("criterion"));
    SET_STRING_ELT(Rnames, 4, mkChar("globalmin"));
    SET_STRING_ELT(Rnames, 5, mkChar("minindex"));
    SET_NAMES(Return_lst, Rnames);

    /* set values */
    SET_VECTOR_ELT(Return_lst, 0, Rthreshold);
    SET_VECTOR_ELT(Return_lst, 1, Rdiscriminability);
    SET_VECTOR_ELT(Return_lst, 2, Rdepth);
    SET_VECTOR_ELT(Return_lst, 3, Rcriterion);
    SET_VECTOR_ELT(Return_lst, 4, Rglobalmin);
    SET_VECTOR_ELT(Return_lst, 5, Rminidx);

    UNPROTECT(nprotect);
    return Return_lst;
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    return R_NilValue;
}

SEXP lmclus(SEXP Xs, SEXP maxDim, SEXP numOfClus, SEXP noiseSize, SEXP bestBound,
    SEXP errorBound, SEXP maxBinPortion, SEXP hisSampling, SEXP hisConstSize,
    SEXP sampleHeuristic, SEXP sampleFactor, SEXP randomSeed, SEXP showLog,
    SEXP hisThr, SEXP algnBasis, SEXP zdSearch) {

    Rprintf("Linear manifold clustering...\n");
    int n, m, show_log, nprotect = 0;
    size_t i, j, k;

    try {
    // Set parameters
    clustering::lmclus::Parameters params;
    params.MAX_DIM = INTEGER(maxDim)[0];
    params.NUM_OF_CLUS = INTEGER(numOfClus)[0];
    params.NOISE_SIZE = static_cast<unsigned int>(INTEGER(noiseSize)[0]);
    params.BEST_BOUND = REAL(bestBound)[0];
    params.ERROR_BOUND = REAL(errorBound)[0];
    params.MAX_BIN_PORTION = REAL(maxBinPortion)[0];
    params.HIS_SAMPLING = static_cast<bool>(INTEGER(hisSampling)[0]);
    params.CONST_SIZE_HIS = INTEGER(hisConstSize)[0];
    params.SAMPLING_HEURISTIC = INTEGER(sampleHeuristic)[0];
    params.SAMPLING_FACTOR = REAL(sampleFactor)[0];
    params.RANDOM_SEED = static_cast<unsigned int>(INTEGER(randomSeed)[0]);
    params.HIS_THR = INTEGER(hisThr)[0];
    params.ALIGN_BASIS = INTEGER(algnBasis)[0];
    params.ZEROD_SEARCH = INTEGER(zdSearch)[0];

    show_log = INTEGER(showLog)[0];
    cpplog::BaseLogger *log;
    if (!show_log)
        log = new cpplog::NullLogger();
    else
        log = new cpplog::FileLogger("/tmp/lmclus.log");

    LOG_INFO(log) << params;

    // get dataset
    SEXP Rdim = getAttrib(Xs, R_DimSymbol);
    n = INTEGER(Rdim)[0];
    m = INTEGER(Rdim)[1];
    Xs = AS_NUMERIC(Xs);
    arma::mat data(REAL(Xs), n, m, false);
    Rprintf("Data dims: (%d, %d)\n", n, m);

    if (params.MAX_DIM >= m) {
        Rprintf("Linear manifold dimension must be less than "
                "the dimension of the data !!!\n");
        return R_NilValue;
    }

    std::vector<arma::uvec> labels;
    std::vector<int> clusterDims;
    std::vector<clustering::lmclus::Separation> separations;

    clustering::lmclus::LMCLUS lmclus(log);
    lmclus.cluster(data, params, labels, clusterDims, separations, output_callback);

    Rprintf("Clusters found: %d\n", labels.size());

    SEXP Return_lst, Rnames, Rthresholds, RclusterDims, Rlabels,
        Rbases, Rorigins, Rhistograms, Rhmins, Rdistances;

     // Thresholds
    PROTECT(Rthresholds = allocVector(REALSXP, separations.size())); nprotect++;
    for (i = 0; i < separations.size(); ++i)
        REAL(Rthresholds)[i] = separations[i].get_threshold();

    // Dimensions
    PROTECT(RclusterDims = allocVector(INTSXP, clusterDims.size())); nprotect++;
    for (i = 0; i < clusterDims.size(); ++i)
      INTEGER(RclusterDims)[i] = clusterDims[i];

    // Labels
    PROTECT(Rlabels = allocVector(VECSXP, labels.size())); nprotect++;
    for (i = 0; i < labels.size(); ++i) {
        labels[i] += 1; // adjust indexes
        SEXP lbls;
        PROTECT(lbls = allocVector(INTSXP, labels[i].n_elem)); nprotect++;
        for (j = 0; j < labels[i].n_elem; ++j)
            INTEGER(lbls)[j] = labels[i][j];
        SET_VECTOR_ELT(Rlabels, i, lbls);
    }

    // Bases
    PROTECT(Rbases = allocVector(VECSXP, separations.size())); nprotect++;
    for (i = 0; i < separations.size(); ++i) {
        arma::mat basis = separations[i].get_projection();
        unsigned int r = basis.n_rows, c = basis.n_cols;

        SEXP bss;
        PROTECT(bss = allocMatrix(REALSXP, r, c)); nprotect++;
        for (j = 0; j < r; ++j)
            for (k = 0; k < c; ++k)
                REAL(bss)[j+r*k] = basis.at(j, k);
        SET_VECTOR_ELT(Rbases, i, bss);
    }

    // Origins
    PROTECT(Rorigins = allocVector(VECSXP, separations.size())); nprotect++;
    for (i = 0; i < separations.size(); ++i) {
        arma::rowvec origin = separations[i].get_origin();
        SEXP orgn;
        PROTECT(orgn = allocVector(REALSXP, origin.n_elem)); nprotect++;
        for (j = 0; j < origin.n_elem; ++j)
            REAL(orgn)[j] = origin[j];
        SET_VECTOR_ELT(Rorigins, i, orgn);
    }

    // Distance histogram
    PROTECT(Rhistograms = allocVector(VECSXP, separations.size())); nprotect++;
    for (i = 0; i < separations.size(); ++i) {
        arma::vec histogram = separations[i].get_histogram();
        SEXP hist;
        PROTECT(hist = allocVector(REALSXP, histogram.n_elem)); nprotect++;
        for (j = 0; j < histogram.n_elem; ++j)
            REAL(hist)[j] = histogram[j];
        SET_VECTOR_ELT(Rhistograms, i, hist);
    }

    // Histogram distances
    PROTECT(Rdistances = allocVector(VECSXP, separations.size())); nprotect++;
    for (i = 0; i < separations.size(); ++i) {
        arma::vec distances = separations[i].get_distances();
        SEXP dist;
        PROTECT(dist = allocVector(REALSXP, distances.n_elem)); nprotect++;
        for (j = 0; j < distances.n_elem; ++j)
            REAL(dist)[j] = distances[j];
        SET_VECTOR_ELT(Rdistances, i, dist);
    }

    // Histogram minimum
    PROTECT(Rhmins = allocVector(INTSXP, separations.size())); nprotect++;
    for (i = 0; i < separations.size(); ++i)
      INTEGER(Rhmins)[i] = separations[i].get_global_min();

    // Result list
    PROTECT(Return_lst = allocVector(VECSXP, 8)); nprotect++;

    /* set names */
    PROTECT(Rnames = NEW_CHARACTER(8)); nprotect++;
    SET_STRING_ELT(Rnames, 0, mkChar("thresholds"));
    SET_STRING_ELT(Rnames, 1, mkChar("cluster_dimensions"));
    SET_STRING_ELT(Rnames, 2, mkChar("clusters"));
    SET_STRING_ELT(Rnames, 3, mkChar("bases"));
    SET_STRING_ELT(Rnames, 4, mkChar("origins"));
    SET_STRING_ELT(Rnames, 5, mkChar("histograms"));
    SET_STRING_ELT(Rnames, 6, mkChar("global_mins"));
    SET_STRING_ELT(Rnames, 7, mkChar("distances"));
    SET_NAMES(Return_lst, Rnames);

    /* set values */
    SET_VECTOR_ELT(Return_lst, 0, Rthresholds);
    SET_VECTOR_ELT(Return_lst, 1, RclusterDims);
    SET_VECTOR_ELT(Return_lst, 2, Rlabels);
    SET_VECTOR_ELT(Return_lst, 3, Rbases);
    SET_VECTOR_ELT(Return_lst, 4, Rorigins);
    SET_VECTOR_ELT(Return_lst, 5, Rhistograms);
    SET_VECTOR_ELT(Return_lst, 6, Rhmins);
    SET_VECTOR_ELT(Return_lst, 7, Rdistances);

    UNPROTECT(nprotect);

    // Clear log
    delete log;

    return Return_lst;
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    return R_NilValue;
}
}
