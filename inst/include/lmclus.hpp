 /* ********************************************************************************
  * * TITLE:        lmclus.hpp
  * *
  * * PURPOSE:      A linear manifold clustering algorithm based on the paper
  * *               "Linear manifold clustering in high dimensional spaces by stochastic search", 
  * *                Pattern Recognition (2007), vol. 40(10), pp 2672-2684.
  * *
  * * AUTHOR:       Rave Harpaz 
  * *               Pattern Recognition Laboratory
  * *               Department of Computer Science
  * *               The Graduate Center
  * *               The City University of New York
  * *               365 Fifth Avenue, New York, New York 10016
  * *               email: rbharpaz@sci.brooklyn.cuny.edu
  * * DATE:         01/01/2005
  * *
  * * VERSION:      1.00
  * *
  * * LANGUAGE:     C++
  * *
  * * SYSTEM:       Ubuntu 64 bit linux workstation using kernal 2.6.15-26-amd64
  * *
  * * COMPILER:     gcc/g++ compiler version 4.0.3
  * *
  * * REFERENCES:   Rave Harpaz and Robert Haralick, 
  * *               "Linear manifold clustering in high dimensional spaces by stochastic search", 
  * *               Pattern Recognition (2007), vol. 40(10), pp 2672-2684.
  * *
  * * REVISIONS:    Art Diky
  * *               Compiled as R-module with armadillo and boost
  * *               email: adiky@gc.cuny.edu
  * *               02/18/2013
  * *
  * * Copyright 2005-2013 Pattern Recognition Laboratory, The City University of New York
  * **********************************************************************************/

#include <iostream>
#include <armadillo>

#include "separation.hpp"
#include "Kittler.h"

#define CPPLOG_FILTER_LEVEL 2
#include "cpplog.hpp"

struct Parameters
{
    int MAX_DIM;
    int NUM_OF_CLUS;
    int LABELED_DATA;
    int CONST_SIZE_HIS;
    unsigned int NOISE_SIZE;
    double BEST_BOUND;
    double ERROR_BOUND; 
    double MAX_BIN_PORTION;
    
    unsigned long int RANDOM_SEED;
    int SAMPLING_HEURISTIC;
    double SAMPLING_FACTOR;
    bool HIS_SAMPLING;
    bool SAVE_RESULT;
    
    friend ostream & operator<<(ostream &o, Parameters &p)
    {
        o<<"MAX_DIM="<<p.MAX_DIM<<endl;
        o<<"NUM_OF_CLUS="<<p.NUM_OF_CLUS<<endl;
        o<<"BEST_BOUND="<<p.BEST_BOUND<<endl;
        o<<"ERROR_BOUND="<<p.ERROR_BOUND<<endl;
        o<<"LABELED_DATA="<<p.LABELED_DATA<<endl;
        o<<"CONST_SIZE_HIS="<<p.CONST_SIZE_HIS<<endl;
        o<<"MAX_BIN_PORTION="<<p.MAX_BIN_PORTION<<endl;
        o<<"NOISE_SIZE="<<p.NOISE_SIZE<<endl;
        o<<"RANDOM_SEED="<<p.RANDOM_SEED<<endl;
        o<<"SAMPLING_HEURISTIC="<<p.SAMPLING_HEURISTIC<<endl;
        o<<"SAMPLING_FACTOR="<<p.SAMPLING_FACTOR<<endl;
        o<<"HIS_SAMPLING="<<p.HIS_SAMPLING<<endl;
        o<<"SAVE_RESULT="<<p.SAVE_RESULT<<endl;

        return o;
    }

};


class LMCLUS
{
private:
    LMCLUS(const LMCLUS& rhs) = delete;
    void operator=(const LMCLUS& rhs) = delete;
    
    // sampling functions
    int sampleQuantity(int lmDim, int fullSpcDim, const int dataSize, const Parameters &para);
    arma::uvec samplePoints(const arma::mat &data, const int lmDim);
    arma::uvec sample(const int n, const int k);

    // basis generation functions
    std::pair<arma::rowvec, arma::mat> formBasis(const arma::mat &points);
    arma::mat gramSchmidtOrthogonalization(const arma::mat &M);
        
    // spearation detection functions
    Separation findBestSeparation(const arma::mat &data, const int SubSpaceDim, const Parameters &para);
    //std::pair<arma::uvec, arma::uvec> findBestPoints(const arma::mat &data, const Separation &best_sep);
    
    arma::mat findBestPoints(const arma::mat &data, const Separation &sep, arma::mat &nonClusterPoints);
    arma::vec determineDistances(const arma::mat &data, const arma::mat &P, const arma::rowvec &origin, const Parameters &para);
    
    unsigned int randromNumber();
    
    cpplog::OstreamLogger *log;
    bool logCreated;
    std::mt19937 engine;
    std::uniform_int_distribution<unsigned int> dist;
    
public:
    LMCLUS() : logCreated(true), engine(std::random_device{}()), dist(std::uniform_int_distribution<unsigned int>())
    {
        // Setup logger
        log = new cpplog::StdErrLogger();
    }    
    
    LMCLUS(cpplog::OstreamLogger *mlog) : logCreated(false), engine(std::random_device{}()), dist(std::uniform_int_distribution<unsigned int>())
    {
        log = mlog;
    }    
    
    ~LMCLUS()
    {
        if (logCreated)
            delete log;
    }
    
    void cluster(const arma::mat &data, const Parameters &para, 
                 std::vector<arma::uvec> &labels, std::vector<double> &thresholds, std::vector<arma::mat> &basises, std::vector<int> &clusterDims);
};
