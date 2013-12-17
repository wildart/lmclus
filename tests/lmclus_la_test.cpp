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
    
    // Clear log
    delete log;
    
    return EXIT_SUCCESS;
}
