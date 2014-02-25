 /* Copyright 2005-2014 Pattern Recognition Laboratory, The City University of New York
  * *******************************************************************************
  * * TITLE:        lmclus.cpp
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
  * *               Rewrote with armadillo
  * *               email: adiky@sci.brooklyn.cuny.edu
  * *               02/18/2013
  * *  
  * **********************************************************************************/

#include <vector>
#include <random>
#include <functional>
#include <limits>
#include <chrono>

#include "Kittler.hpp"
#include "lmclus.hpp"

#define EPS 1e-8

/* Random number generator
 */
unsigned int clustering::lmclus::LMCLUS::randromNumber() {
    return dist(engine);
}

// sampling functions ----------------------------------------------------------------------------------------------

/* SampleQuantity
 * --------------
 * determine the number of times to sample the data in order to guaranty
 * that the points sampled are from the same cluster with probability
 * of error that does not exceed an error bound. 3 different types of heuristics may be used
 * depending on LMCLUS's input parameters.
 */
int clustering::lmclus::LMCLUS::sampleQuantity(int LMDim, int fullSpcDim, const int DataSize, const Parameters &para)
{
    double k=static_cast<double>(para.NUM_OF_CLUS);

    if(k==1)                               // case where there is only one cluster
        return 1;

    double p=1/k;                          // p=probability that 1 point comes from a certain cluster

    double P=pow(p,LMDim);                   // P=probability that "k+1" points are from the same cluster

    double N=(log10(para.ERROR_BOUND))/(log10(1-P));

    int NumOfSamples = 0;

    LOG_DEBUG(log) << "number of samples by first heuristic="<<N<<", by second heuristic="<<DataSize*para.SAMPLING_FACTOR;

    switch(para.SAMPLING_HEURISTIC){
        case 1: {
            NumOfSamples=static_cast<int>(N);
            break;
        }
        case 2:{
            NumOfSamples=DataSize*para.SAMPLING_FACTOR;
            break;
        }
        case 3: {
            if( N < (DataSize*para.SAMPLING_FACTOR))
                NumOfSamples=static_cast<int>(N);
            else
                NumOfSamples=DataSize*para.SAMPLING_FACTOR; 
        }
    }
            
    LOG_DEBUG(log) << "number of samples="<<NumOfSamples;

    return NumOfSamples;
}

/* samplePoints
 * ------------
 * Sample randomly LMDim+1 points from the dataset, making sure
 * that the same point is not sampled twice. the function will return
 * a index vector, specifying the index of the sampled points.
 */
arma::uvec clustering::lmclus::LMCLUS::samplePoints(const arma::mat &data, const int LMDim)
{
    size_t NumOfPoints = LMDim+1, empty = data.n_rows+1;
    arma::uvec point_index(NumOfPoints);  
    bool match, zero;
    
    std::uniform_int_distribution<unsigned int>::param_type newParams{0, data.n_rows-1};
    dist.param(newParams);
    
    point_index.fill(empty);
    size_t count = 0;
    while ( count < NumOfPoints ) {
        size_t index = randromNumber();
        LOG_TRACE(log) << "Point index: " << index;
        size_t found = arma::sum(point_index == index);  
        if (found == 0)
        {
            match = false;
            zero = true;
            for(size_t i=0; i<NumOfPoints && point_index(i) != empty; i++)
            {
                match = true;
                zero = true;
                size_t ii = point_index(i);
                LOG_TRACE(log) << "Check points: " << i << " <-> " << ii;
                for(size_t j=0; j < data.n_cols; j++)
                {
                   match = match && ((data(index,j) - data(ii,j)) < EPS); 
                   if (!match) break;
                   zero  = zero && data(index,j) == 0.0;
                }
                if (!match || zero) break;
            }
            // new point does not match or zero
            if (!match || zero)
            {
                LOG_TRACE(log) << "Added index " << index << ", # " << count;
                point_index(count) = index;
                count++;
            }
        }
    }
    return point_index;
}

/* sample
 * ------------
 * Sample uniformly k integers from the integer range 0:n-1, making sure that
 * the same integer is not sampled twice. the function returns an intger vector
 * containing the sampled integers.
 */
arma::uvec clustering::lmclus::LMCLUS::sample(const int n, const int k)
{
    arma::uvec index = arma::zeros<arma::uvec>(n);
    arma::uvec SampleIndex(k);

    // setup rng limits
    std::uniform_int_distribution<unsigned int>::param_type newParams(0, n);
    dist.param(newParams);

    for(int i=0; i<k; i++) {
        auto rn = randromNumber();

        while(index(rn) != 0 )                       // make sure it was not already chosen
            rn= randromNumber();
        index(rn)=1;
        SampleIndex(i) = rn;
    }

    return SampleIndex;
}

// Basis generation functions ----------------------------------------------------------------------------------------------

/* FormBasis
 * ---------
 * the idea is to pick a point (origin) from the sampled points and generate
 * the basis vectors by subtracting all other points from the origin,
 * creating a basis matrix with one less vector than the number of sampled points.
 * the basis matrix generated in this function is the transpose of the basis matrix. 
 * the function also updates the origin of the basis vectors.
 */
arma::mat clustering::lmclus::LMCLUS::formBasis(const arma::mat &points, arma::rowvec& origin)
{
    origin = points.row(0);  // let the origin be the fisrt point sampled
    arma::mat basis(points.n_rows-1, points.n_cols);   // create the B (Basis) transpose matrix

    // find the b_i-th basis vector by subtracting each other point (vector) from 'x_1'
    arma::rowvec b_i;
    for (size_t i=1; i < points.n_rows ; i++) {
        b_i = points.row(i) - origin;
        basis.row(i-1) = b_i;   // set the i-th row in the transpose basis matrix to its b_i-th basis vector
    }
    
    return basis;
}

/* GramSchmidtOrthogonalization
 * -----------------------------
 * the idea of the Gram-Schmidt orthogonalization process is to subtract from
 * every new (unsettled basis vector) its components in the directions that
 * are already settled, and then make it a unit vector (orthonormal).
 *
 * more formally:
 * --------------
 * b_i  : is the already settled orthogonal i-th basis vector
 * b_j  : is the already settled orthogonal j-th basis vector
 * b_i' : is the transpose of b_i
 * m_i  : is an unsettled i-th basis vector (coming from the original set of unsettled basis vectors)
 * sum(...) : the sum is taken over all j's from 0 to i-1
 * || b_i || : is the norm (length of the vector)
 * to calculate b_i use : b_i = ( m_i - sum( (b_j'm_i)b_j ) ) / || b_i ||
 */
arma::mat clustering::lmclus::LMCLUS::gramSchmidtOrthogonalization(const arma::mat &M)
{
    arma::mat B( M.n_rows, M.n_cols );       // create an uninitialized orthogonal basis matrix

    for(unsigned int i=0; i< M.n_rows; i++) {// for each vector in the original basis convert to an orthogonal vector
        auto m_i = M.row(i);                 // the i-th original basis vector
        
        arma::rowvec x_i( M.n_cols );
        x_i.zeros();

        double x_j;
        for (unsigned int j=0; j < i; j++) { // calculate the sum of i-th vector's components in the directions of the
            // already settled orthogonal basis vectors
            auto b_j = B.row(j);
            x_j = arma::dot( b_j, m_i );
            x_i += (  b_j * x_j );
        }
        arma::rowvec b_i = m_i - x_i;        // subtract the sum calculated above (stored in x_i) from the i-th unsettled vector
        double n= arma::norm(b_i, 2);        // make b_i a unit vector
        if(n!=0.0)
            b_i = b_i * (1/n);        
        B.row(i) = b_i;
    }
    //LOG_TRACE(log) << "B: \n" << B;
    return B;
}    

// Separation detection functions ----------------------------------------------------------------------------------------------

/* DetermineDistances
 * ------------------
 * determine the distance of each point in the data set from to a linear manifold,
 *  using: d_n= || (I-P)(z_n-origin) ||
 * where P is the projection operator matrix, and z_n is the point.
 * (the distance will not be determined for points that were already sampled to create the linear manifold).
 * depending on LMCLUS's input parameters only a sample of distances will be computed to enhance efficiency.
 */
arma::vec clustering::lmclus::LMCLUS::determineDistances(
    const arma::mat &data, const arma::mat &P, 
    const arma::rowvec &origin, const Parameters &para)
{
    arma::mat data1 = data;

    if(para.HIS_SAMPLING) 
    {
        double Z_01=2.576;      // Z random variable, confidence interval 0.99
        double delta_p=0.2;
        double delta_mu=0.1;
        double P=1/static_cast<double>(para.NUM_OF_CLUS);
        double Q=1-P;
        double n1=(Q/P)*((Z_01*Z_01)/(delta_p*delta_p));
        double p=( P<=Q ? P : Q );
        double n2=(Z_01*Z_01)/(delta_mu*delta_mu*p);
        double n3= ( n1 >= n2 ? n1 : n2 );
        unsigned int n4= static_cast<int> (n3);
        int n= ( data.n_rows <= n4 ? data.n_rows-1 : n4 );

        auto sampleIndex = sample(data.n_rows, n);
        data1 = data.rows(sampleIndex);
    }
    
    // vector to hold distances of points from basis
    arma::vec Distances = arma::zeros<arma::vec>(data1.n_rows);
    unsigned int i;
    #pragma omp parallel for private(i) shared(data1, origin, P, Distances)
    for (i=0; i < Distances.n_rows; i++) {
        Distances(i) = distanceToManifold(data1.row(i) - origin, P);
    }

    return Distances;
}

/* distanceToManifold
 * ------------------
 * calculates distance from point to manifold defined by basis and origin
 * using: d_n= || (I-P)x || (1)
 * where P is the projection operator matrix, and z_n is the point.
 * But calculations are optimized in a following way
 * || x ||^2 = || Px + (I-P)x ||^2
 * || x ||^2 = || Px ||^2 + || (I-P)x ||^2
 * because P = BB' and P^2 = BB'BB' = BIB' = BB' = P then
 * || x ||^2 = || BB'x ||^2 + || (I-P)x ||^2
 * but || BB'x ||^2 = (BB'x)'BB'x = x'BB'BB'x = x'BIB'x = x'BB'x = = (B'x)'B'x = || B'x ||^2
 * so || x ||^2 = || B'x ||^2 + || (I-P)x ||^2
 * || x ||^2 = || B'x ||^2 + d_n from (1)
 * thus d_n = || x ||^2 -|| B'x ||^2
 * where B is a basis matrix and B' - its transpose.
 * Note: depending on LMCLUS's input parameters only a sample of distances will be computed to enhance efficiency.
 */
double clustering::lmclus::LMCLUS::distanceToManifold(
    const arma::rowvec &point, const arma::mat &B_T)
{
    double d_n = 0.0, c, b;
    arma::vec d_v = B_T * point.t();
    c=arma::norm(point, 2);
    b=arma::norm(d_v, 2);
    d_n=sqrt((c*c)-(b*b));
    if(d_n<0 && d_n>1000000000)
        d_n = 0.0;
    return d_n;
}

/* findBestSeparation
 * --------------------
 * LMCLUS's main function:
 * 1- sample trial linear manifolds by sampling points from the data
 * 2- create distance histograms of the data points to each trial linear manifold
 * 3- of all the linear manifolds sampled select the one whose associated distance histogram shows the best separation between to modes.
 */
clustering::lmclus::Separation clustering::lmclus::LMCLUS::findBestSeparation (const arma::mat &data, const int LMDim, const Parameters &params)
{
    int DataSize = data.n_rows;
    int FullSpcDim = data.n_cols;

    LOG_INFO(log)<<"data size="<<DataSize<<"   linear manifold dim="<<LMDim<<"   space dim="<<FullSpcDim<<"   searching for separation ...";
    LOG_INFO(log)<<"------------------------------------------------------------";
    Separation best_sep;                                // contains info about best separation

    // determine number of samples of "LMDim+1" points
    int i;
    int Q = LMCLUS::sampleQuantity( LMDim, FullSpcDim, DataSize, params );
    LOG_TRACE(log) << "Collect " << Q << " sample(s)";
    
    #pragma omp parallel for private(i) shared(best_sep, data, Q)
    for (i=0; i < Q; i++) {                         // sample Q times SubSpaceDim+1 points
        LOG_TRACE(log) << "Iteration: " << i;
        arma::uvec points_idx = samplePoints(data, LMDim);   //  sample LMDim+1 points
        if (points_idx.n_elem < static_cast<unsigned>(LMDim+1)) continue;        

        // form basis (transpose) of the linear manifold spanned by the sampled points (vectors)
        arma::mat sample = data.rows(points_idx);
        arma::rowvec origin = sample.row(0); 
        auto basis = formBasis(sample, origin);
        LOG_TRACE(log) << "Origin: \n" << origin << "Collected basis :\n" << basis <<std::endl;
        // orthogonalize Basis ( with orthonormal basis-vectors)
        auto B_T = gramSchmidtOrthogonalization(basis);           
        LOG_TRACE(log) << "Orthogonal basis:\n" << basis;
        
        // determine distances of points to the linear manifold
        arma::vec distances = determineDistances(data, B_T, origin, params);
        double RHmin=distances.min();
        double RHmax=distances.max(); 
        LOG_TRACE(log) << "Distances from " << RHmin << " to " << RHmax;

        // generate a distances histogram
        arma::uvec hist = arma::hist(distances, params.CONST_SIZE_HIS>0 ? params.CONST_SIZE_HIS : static_cast<int>(distances.n_elem * params.MAX_BIN_PORTION));
        
        // Cannot calculate distances
        auto filledBins = find(hist > 0).eval();
        LOG_TRACE(log) << "Non-empty bins: \n" << filledBins.n_elem;
        if ( filledBins.n_elem < 2)
            continue;
        
        // Threshold histogram and determine goodness of separation
        arma::vec histNorm = arma::conv_to< arma::vec >::from(hist) / static_cast<double>(distances.n_elem);
        LOG_TRACE(log) << "Create histogram: \n" << histNorm.t();
        
        Kittler K;
        K.FindThreshold(histNorm, RHmin, RHmax);
        LOG_TRACE(log) << "Found threshold: " << K.GetDiscrim()*K.GetDepth();

        Separation sep(K.GetDiscrim(), K.GetDepth(), K.GetThreshold(), origin, B_T, hist);  // store separation info
        
        // keep track of best histogram/separation
        #pragma omp critical(find_separation)
        if ( (K.GetDiscrim()*K.GetDepth() ) > best_sep.get_criteria()  ){
            LOG_TRACE(log) << "Found orthogonal basis:\n" << sep.get_projection();
            best_sep=sep;
        }
    }

    if (best_sep.get_criteria()==0)
        LOG_DEBUG(log) << "no good histograms to separate data !!!";
    else {
        LOG_DEBUG(log) <<"sep width="<<best_sep.get_width()<<"  sep depth="<<best_sep.get_depth() << "  sep criteria="<<best_sep.get_criteria();
    }
    return best_sep;
}


void clustering::lmclus::LMCLUS::find_manifold(const arma::mat &data, const Parameters &para, 
                 arma::uvec &points_index, std::vector<unsigned int> &nonClusterPoints,
                 std::vector<Separation> &separations, bool &Noise, int &SepDim)
{
    for(int lm_dim=1; lm_dim < para.MAX_DIM+1 && !Noise; lm_dim++) {
        while(true) {
            // find the best fit of a set of points to a linear manifold of dimensionality lm_dim
            Separation best_sep=findBestSeparation(data.rows(points_index), lm_dim, para);
            //LOG_TRACE(log) << "Proj: \n"<< best_sep.get_projection();
            LOG_DEBUG(log) << "BEST_BOUND: "<< best_sep.get_criteria() << "(" << para.BEST_BOUND << ")";
            
            if (best_sep.get_criteria() < para.BEST_BOUND ) break;
            SepDim=lm_dim;
            
            // find which points are close to manifold 
            std::vector<unsigned int> best_points;
            double threshold = best_sep.get_threshold();
            arma::mat Projection = best_sep.get_projection();
            arma::rowvec origin = best_sep.get_origin();
            
            LOG_TRACE(log) << "Threshold: " << threshold;
            LOG_TRACE(log) << "Origin: \n" << origin;

            unsigned int i, idx;
            //#pragma omp parallel for private(i, idx) shared(points_index, data)
            for(i=0; i < points_index.n_rows; i++) {
                idx = points_index(i);
                double d_n = distanceToManifold(data.row(idx) - origin, Projection);
                if(d_n < threshold)    // point i has distances less than the threshold value
                    best_points.push_back(idx);          // add to best points
                else
                    nonClusterPoints.push_back(idx);
            }
            
            LOG_DEBUG(log) << "Separated points: "<< best_points.size();
            points_index = arma::conv_to<arma::uvec>::from(best_points);                
                            
            if(points_index.n_rows< para.NOISE_SIZE) {       // small amount of points is considered noise
                Noise=true;
                LOG_INFO(log)<<"noise less than "<<para.NOISE_SIZE<<" points";
                break;
            }
            
            separations.push_back(best_sep);
            LOG_DEBUG(log) << "Best basis:"<< std::endl << best_sep.get_projection();
        }

        LOG_INFO(log) << "# of points = " << points_index.n_rows;

        if(lm_dim<para.MAX_DIM && !Noise)
            LOG_INFO(log)<<"no separation, increasing dim ...";
        else
            LOG_INFO(log)<<"no separation ...";

    }
}

void clustering::lmclus::LMCLUS::cluster(const arma::mat &data, const Parameters &para, 
    std::vector<arma::uvec> &labels, std::vector<double> &thresholds, 
    std::vector<arma::mat> &bases, std::vector<int> &clusterDims, 
    std::vector<arma::vec> &origins, callback_t progress)
{
    // Initialize random generator
    if(para.RANDOM_SEED>0)
        engine.seed(para.RANDOM_SEED);
    else 
    {
        // obtain a seed from the system clock:
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        engine.seed(seed);
    }
    
    std::vector<unsigned int> noise;
    bool groupNoise = true;

    std::vector<Separation> separations;
    arma::uvec points_index = arma::linspace<arma::uvec>(0, data.n_rows, data.n_rows);

    int ClusterNum=0;
    do {
        bool Noise=false;
        std::vector<unsigned int> nonClusterPoints;
        int SepDim=0;                             // dimension in which separation was found
        
        //find_manifold(data, para, points_index, nonClusterPoints, separations, Noise, SepDim);

        for(int lm_dim=1; lm_dim < para.MAX_DIM+1 && !Noise; lm_dim++) {
            while(true) {
                // find the best fit of a set of points to a linear manifold of dimensionality lm_dim
                Separation best_sep=findBestSeparation(data.rows(points_index), lm_dim, para);
                //LOG_TRACE(log) << "Proj: \n"<< best_sep.get_projection();
                LOG_DEBUG(log) << "BEST_BOUND: "<< best_sep.get_criteria() << "(" << para.BEST_BOUND << ")";
                if (best_sep.get_criteria() < para.BEST_BOUND ) break;
                SepDim=lm_dim;
                
                // find which points are close to manifold 
                std::vector<unsigned int> best_points;
                double threshold = best_sep.get_threshold();
                arma::mat Projection = best_sep.get_projection();
                arma::rowvec origin = best_sep.get_origin();
                
                LOG_TRACE(log) << "Threshold: " << threshold;
                LOG_TRACE(log) << "Origin: \n" << origin;

                unsigned int i, idx;
                /*
                for(i=0; i < points_index.n_rows; i++) {
                    idx = points_index(i);
                    arma::rowvec Z_n  = data.row(idx) - origin; // point with respect to basis
                    // calculate distance of point from basis:  d_n= || (I-P)(z_n-origin) ||
                    arma::vec d_v = Projection * Z_n.t();
                    double d_n=0;
                    double c = norm(Z_n, 2);
                    double b = norm(d_v, 2);
                    d_n=sqrt(fabs((c*c)-(b*b)));
                    if(d_n < threshold)    // point i has distances less than the threshold value
                        best_points.push_back(idx);          // add to best points
                    else
                        nonClusterPoints.push_back(idx);
                } */
                for(i=0; i < points_index.n_rows; i++) {
                    idx = points_index(i);
                    double d_n = distanceToManifold(data.row(idx) - origin, Projection);
                    if(d_n < threshold)    // point i has distances less than the threshold value
                        best_points.push_back(idx);          // add to best points
                    else
                        nonClusterPoints.push_back(idx);
                }
                LOG_DEBUG(log) << "Separated points: "<< best_points.size();
                points_index = arma::conv_to<arma::uvec>::from(best_points);                
                                
                if(points_index.n_rows< para.NOISE_SIZE) {       // small amount of points is considered noise
                    Noise=true;
                    LOG_INFO(log)<<"noise less than "<<para.NOISE_SIZE<<" points";
                    break;
                }
                
                separations.push_back(best_sep);
                LOG_DEBUG(log) << "Best basis:"<< std::endl << best_sep.get_projection();
            }

            LOG_INFO(log) << "# of points = " << points_index.n_rows;

            if(lm_dim<para.MAX_DIM && !Noise)
                LOG_INFO(log)<<"no separation, increasing dim ...";
            else
                LOG_INFO(log)<<"no separation ...";

        }

        // second phase (steps 9-12)
        ClusterNum++;
        std::string msg("found cluster # "+to_string(ClusterNum)+
            ", size="+to_string(points_index.n_rows)+", dim="+to_string(SepDim));
        LOG_INFO(log)<<msg;
        if(progress != nullptr) 
            progress(msg.c_str());
        
        if (Noise && groupNoise) {
            for(size_t i = 0; i < points_index.n_rows; ++i)
                noise.push_back(points_index(i));
            ClusterNum--;
        } else {
            clusterDims.push_back(SepDim);
            //labels.push_back(points_index);     
            if (SepDim > 0)
                labels.push_back(points_index);     
            else
                // if dimension is 0 then dataset is unseparable so list it as cluster
                labels.push_back(points_index);
            
            // Save cluster basis
            if (separations.size() > 0)
            {
                Separation bs = separations[separations.size()-1];
                bases.push_back(bs.get_projection());
                thresholds.push_back(bs.get_threshold());
                origins.push_back(bs.get_origin());
            
                LOG_TRACE(log) << "Basis: \n"<< bs.get_projection();
                LOG_TRACE(log) << "Origin: " << bs.get_origin();
                LOG_TRACE(log) << "Threshold: " << bs.get_threshold();
            }
            else
            {
                thresholds.push_back(0.0);
                bases.push_back(arma::zeros<arma::mat>(1,para.MAX_DIM+1));
                origins.push_back(arma::zeros<arma::vec>(para.MAX_DIM+1));
            }    
        }
        
        // separate cluster points from the rest of the data
        points_index = arma::conv_to<arma::uvec>::from(nonClusterPoints);
        separations.clear();
        
        // Stop clustering if we reached limit
        if (ClusterNum == para.NUM_OF_CLUS) break;

    } while(points_index.n_rows > para.NOISE_SIZE);
    
    // take care of remaining points
    if(points_index.n_rows>0)
    {
        for(unsigned int c=0; c<points_index.n_rows; ++c)
            noise.push_back(points_index(c));
    }
    
    // Add noise as cluster
    if (noise.size() > 0)
    {
        auto noise_index = arma::conv_to<arma::uvec>::from(noise);
        clusterDims.push_back(0);
        thresholds.push_back(0.0);        
        bases.push_back(arma::zeros<arma::mat>(1,para.MAX_DIM+1));
        origins.push_back(arma::zeros<arma::vec>(para.MAX_DIM+1));
        labels.push_back(noise_index); 
        LOG_INFO(log)<<"found cluster #(noise) size="<<noise.size()<<" !!!";
    }
    
    return;
}

