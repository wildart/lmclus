/*
 * =====================================================================================
 *
 *       Filename:  lmclus_demo.cpp
 *
 *    Description:  Linear Manifold Clustering
 *
 *        Version:  1.0
 *        Created:  02/19/2013 08:16:54 PM
 *       Revision:  none
 *       Compiler:  gcc 4.7.2
 *
 *         Author:  Art Wild (wildart@gmail.com)
 *        Company:  Patter Recognition Lab, GC CUNY
 *
 * =====================================================================================
 */
#include "lmclus.hpp"

#define TEST(log, cond, msg) if(!(cond)){ LOG_ERROR(log) << msg; return EXIT_FAILURE;}

int main ( int argc, char *argv[] )
{
    cpplog::OstreamLogger *log;
    log = new cpplog::StdErrLogger();
    
    // Load dataset
    arma::mat data{ -1,1,0,1, 
                    0,0,1,0,
                    0,0,1,2,
                    0,0,1,0};
    data.reshape(4,4);
    LOG_INFO(log) << "data:\n" << data;
    
    arma::colvec exp_origin{-1,0,0,0};
    arma::mat exp_basis1{2,1,2, 
                         0,1,0,
                         0,1,2,
                         0,1,0};
    exp_basis1.reshape(3,4);
    
    // Create basis
    clustering::lmclus::LMCLUS lmclus(log);
    arma::rowvec origin = data.row(0);    
    auto basis = lmclus.formBasis(data, origin);
    LOG_INFO(log) << "origin:\n" << origin;
    LOG_INFO(log) << "basis:\n" << basis;    
    TEST(log, sum(exp_origin == origin) == 4, "Invalid origin")
    TEST(log, sum(sum(exp_basis1 == basis)) == 12, "Invalid basis1")

    // Orthoganalization
    arma::mat exp_basis2{1,0,0, 
                         0,1/sqrt(3),-1/sqrt(6),
                         0,1/sqrt(3),2/sqrt(6),
                         0,1/sqrt(3),-1/sqrt(6)};
    exp_basis2.reshape(3,4);
    auto bo = lmclus.gramSchmidtOrthogonalization(basis);
    LOG_INFO(log) << "basis:\n" << exp_basis2;
    LOG_INFO(log) << "basis:\n" << bo;
    TEST(log, sum(sum(round(exp_basis2 - bo))) == 0, "Invalid basis2")



    arma::vec dists{0.02310271, 0, 0.03380884, 0.04338802, 0.04507846, 0.03718973, 0.04843846, 0.07005137, 0.06697438, 0.06818648, 0.06165015, 0.07331724, 0.0572927, 0.06115168, 0.06136745, 0.06916228, 0.07505889, 0.04025576, 0.0348092, 0.06977559, 0.07038524, 0.06660204, 0.06956347, 0.06107903, 0.06637905, 0.05116174, 0.067156, 0.07132017, 0.074072, 0.07345591, 0.06706556, 0.03416455, 0.07011909, 0.04506517, 0.0574775, 0.06259683, 0.03340868, 0.04202428, 0.05425579, 0.05954518, 0.06063532, 0.04872525, 0.072886, 0.04538827, 0.06488559, 0.07271835, 0.06768831, 0.07327199, 0.06472511, 0.07290822, 0.06789454, 0.0654445, 0.06324504, 0.06488559, 0.06893446, 0.06855627, 0.06228776, 0.06645411, 0.06816398, 0.06424475, 0.06488559, 0.06645411, 0.04788272, 0.0503709, 0.08044789, 0.08274521, 0.08274075, 0.000229518, 0.01966824, 0.02400097, 0.03040416, 0.01700431, 0.02230689, 0.009251605, 0.02558243, 0.00778831, 0.01411486, 0.02450959, 0.02142606, 0.02875149, 0.0009508195, 0.02315529, 0.003659238, 0.009273629, 0.01588809, 0.02306605, 0.009173657, 0.02709979, 0.02265084, 0.01876146, 0.01276366, 0.01776059, 0.005630435, 0.003261429, 0.01021569, 0.008393755, 0.02532516, 0.03487231, 0.02801513, 0.007998879, 0.02221921, 0.00717457, 0.01726988, 0.01852432, 0.01742978, 0.0213977, 0.01059003, 0.01074981, 0.001320267, 0.001012724, 0.0009508195, 0.01603255, 0.01672346, 0.0004692429, 0.01419651, 0.001320267, 0.01335964, 0.01544578, 0.002465553, 0.01034083, 0.001146503, 0, 0.01852432, 0.02111023, 0.002860147, 0.01778593, 0.02954474, 0.0286408, 0.06894352, 0.06043474, 0.06406747, 0.06518222, 0.06632843, 0.06237962};
    auto hist = clustering::lmclus::histBootstrapping(dists, 30);
    LOG_INFO(log) << "hist:\n" << hist.t();
    
    // Clear log
    delete log;
    
    return EXIT_SUCCESS;
}
