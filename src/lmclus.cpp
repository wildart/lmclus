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
#include <iomanip>
#include <cassert>

#include "lmclus.hpp"

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
int clustering::lmclus::LMCLUS::sampleQuantity(int LMDim, int fullSpcDim, const int DataSize, const Parameters &params)
{
    double k=static_cast<double>(params.NUM_OF_CLUS);

    if(k==1)                               // case where there is only one cluster
        return 1;

    double p=1/k;                          // p=probability that 1 point comes from a certain cluster

    double P= LMDim == 0 ? 2 : pow(p,LMDim);                   // P=probability that "k+1" points are from the same cluster

    double N=(log10(params.ERROR_BOUND))/(log10(1-P));

    int NumOfSamples = 0;

    LOG_DEBUG(log) << "number of samples by first heuristic="<<N<<", by second heuristic="<<DataSize*params.SAMPLING_FACTOR;

    switch(params.SAMPLING_HEURISTIC){
        case 1: {
            NumOfSamples=static_cast<int>(N);
            break;
        }
        case 2:{
            NumOfSamples=DataSize*params.SAMPLING_FACTOR;
            break;
        }
        case 3: {
            if( N < (DataSize*params.SAMPLING_FACTOR))
                NumOfSamples=static_cast<int>(N);
            else
                NumOfSamples=DataSize*params.SAMPLING_FACTOR;
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
    size_t pindex, pid, NumOfPoints = LMDim+1;
    vector<unsigned int> selected_points;
    vector<unsigned int> point_index(data.n_rows);
    bool zero, unique;

    LOG_TRACE(log) << "Total points: " << data.n_rows;
    std::uniform_int_distribution<unsigned int>::param_type newParams{0, data.n_rows-1};
    dist.param(newParams);

    // make index
    for(size_t i=0; i<data.n_rows-1;++i)
        point_index[i] = i;

    while (point_index.size() > 0) {
        // select random point
        pindex = randromNumber();
        pid = pindex % point_index.size();
        LOG_TRACE(log) << "Point index: " << pid <<
            "(" << pindex << ", " << point_index[pid] << ") - Left: " << point_index.size();

        // check if it zero and matches
        zero = true;
        for(size_t j=0; j < data.n_cols; j++)
            zero  = zero && data(point_index[pid],j) == 0.0;

        unique = true;
        for(size_t i=0; i<selected_points.size(); i++) {
            LOG_TRACE(log) << "Check points: " << selected_points[i] << " <-> " << point_index[pid];
            bool match = true;
            for(size_t j=0; j < data.n_cols; j++)
                match = match && ((data(selected_points[i],j) - data(point_index[pid],j)) < arma::datum::eps);
            unique = unique && !match;
        }

        // select non-zero unique point
        if (!zero && unique) {
            selected_points.push_back(point_index[pid]);
            LOG_TRACE(log) << "Selected: " << selected_points.size();
            if (selected_points.size() == NumOfPoints)
                break;
        }

        // remove selected point from index
        point_index.erase(point_index.begin()+pid);
    }

    return arma::conv_to<arma::uvec>::from(selected_points);
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
    const arma::mat &data, const arma::mat &basis,
    const arma::rowvec &origin, const Parameters &params)
{
    arma::mat data1 = data;

    if(params.HIS_SAMPLING)
    {
        double Z_01=2.576;      // Z random variable, confidence interval 0.99
        double delta_p=0.2;
        double delta_mu=0.1;
        double P=1/static_cast<double>(params.NUM_OF_CLUS);
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
    #pragma omp parallel for private(i) shared(data1, origin, basis, Distances)
    for (i=0; i < Distances.n_rows; i++) {
        Distances(i) = distanceToManifold(data1.row(i) - origin, basis);
    }

    return Distances;
}

/* distanceToManifold
 * ------------------
 * calculates distance from point to manifold defined by basis and origin
 * using: d_n= || (I-P)x || (1)
 * where P is the projection operator matrix, and x is the point.
 * But calculations are optimized in a following way
 * || x ||^2 = || Px + (I-P)x ||^2
 * || x ||^2 = || Px ||^2 + || (I-P)x ||^2
 * because P = BB' and P^2 = BB'BB' = BIB' = BB' = P then
 * || x ||^2 = || BB'x ||^2 + || (I-P)x ||^2
 * but || BB'x ||^2 = (BB'x)'BB'x = x'BB'BB'x = x'BIB'x = x'BB'x = = (B'x)'B'x = || B'x ||^2
 * so || x ||^2 = || B'x ||^2 + || (I-P)x ||^2
 * || x ||^2 = || B'x ||^2 + d_n^2 from (1)
 * thus d_n^2 = || x ||^2 -|| B'x ||^2
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
    d_n=(c*c)-(b*b);
    d_n = (d_n<0) ? 0.0 : sqrt(d_n); //&& d_n>1000000000
    return d_n;
}

double clustering::lmclus::LMCLUS::projectTo1D(
    const arma::rowvec &point, const arma::mat &B_T)
{
    arma::vec d_v = B_T * point.t();
    return d_v(0);
}


/* findBestSeparation
 * --------------------
 * 1- sample trial linear manifolds by sampling points from the data
 * 2- create distance histograms of the data points to each trial linear manifold
 * 3- of all the linear manifolds sampled select the one whose associated distance histogram shows the best separation between to modes.
 */
clustering::lmclus::Separation clustering::lmclus::LMCLUS::findBestSeparation(
    const arma::mat &data, const int LMDim, const Parameters &params) {

    int DataSize = data.n_rows;
    int FullSpcDim = data.n_cols;
    size_t histBins = params.CONST_SIZE_HIS>0 ? params.CONST_SIZE_HIS :
                static_cast<size_t>(DataSize * params.MAX_BIN_PORTION);

    LOG_INFO(log)<<"data size="<<DataSize<<"   linear manifold dim="<<LMDim
        <<"   space dim="<<FullSpcDim<<"   hist bins="<<histBins
        <<"   searching for separation ...";
    LOG_INFO(log)<<"------------------------------------------------------------";
    Separation best_sep;    // contains info about best separation

    // determine number of samples of "LMDim+1" points
    int i;
    int Q = LMCLUS::sampleQuantity( LMDim, FullSpcDim, DataSize, params );
    LOG_TRACE(log) << "Collect " << Q << " sample(s)";

    //#pragma omp parallel for private(i) shared(best_sep, data, Q)
    for (i=0; i < Q; i++) {     // sample Q times SubSpaceDim+1 points
        LOG_TRACE(log) << "Iteration: " << i;
        arma::uvec points_idx = samplePoints(data, LMDim);  //  sample LMDim+1 points
        if (points_idx.n_elem < static_cast<unsigned>(LMDim+1))
            continue;

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
        arma::uvec hist = arma::hist(distances, histBins);

        // Cannot calculate distances
        auto filledBins = find(hist > 0).eval();
        LOG_TRACE(log) << "Non-empty bins: \n" << filledBins.n_elem;
        if ( filledBins.n_elem < 2)
            continue;

        arma::vec histNorm;
        if (histBins >= params.HIS_THR || params.HIS_THR == 0) {
        // Threshold histogram and determine goodness of separation
            histNorm = arma::conv_to< arma::vec >::from(hist) / static_cast<double>(distances.n_elem);
        }
        // Calculate density function
        else {
            LOG_DEBUG(log) << "Start histogram bootstrapping...";
            histNorm = histBootstrapping(distances, params.HIS_THR);
            if ( histNorm.n_elem < 2)
                continue;
        }
        LOG_TRACE(log) << "Create histogram: \n" << histNorm.t();

        Kittler K;
        K.FindThreshold(histNorm, RHmin, RHmax);
        LOG_TRACE(log) << "Found threshold: " << K.GetDiscrim()*K.GetDepth();

        // store separation info
        Separation sep(K.GetDiscrim(), K.GetDepth(), K.GetThreshold(),
                origin, B_T, histNorm, K.GetMinIndex(), distances);

        // keep track of best histogram/separation
        //#pragma omp critical(find_separation)
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

/* findBestZeroManifoldSeparation
 * --------------------
 * 1- sample 0-dimensioanl manifolds origins (points)
 * 2- calculate distance to origin from each point of dataset
 * 3- calculate the densify function (actualy mass function) of distances from origin
      and determine the best separation if it is multimodal
 */
clustering::lmclus::Separation
    clustering::lmclus::LMCLUS::findBestZeroManifoldSeparation(
    const arma::mat &data, const Parameters &params, const Separation &sep) {

    int DataSize = data.n_rows;
    int FullSpcDim = data.n_cols;
    size_t histBins = params.CONST_SIZE_HIS>0 ? params.CONST_SIZE_HIS :
                static_cast<size_t>(DataSize * params.MAX_BIN_PORTION);
    Separation best_sep(FullSpcDim);

    LOG_INFO(log)<<"data size="<<DataSize<<"   linear manifold dim=0"
        <<"   space dim="<<FullSpcDim<<"   hist bins="<<histBins
        <<"   searching for separation ...";
    LOG_INFO(log)<<"------------------------------------------------------------";

    // determine number of samples of "LMDim+1" points
    int i, Q = LMCLUS::sampleQuantity( 1, FullSpcDim, DataSize, params );
    LOG_TRACE(log) << "Collect " << Q << " sample(s)";

    for (i = 0; i < Q; ++i) {
        unsigned int point_idx = rand() % data.n_rows;

        // Calculate distances to origin
        arma::vec distances = arma::zeros<arma::vec>(data.n_rows);
        arma::mat B_T = sep.get_projection();
        arma::rowvec origin = sep.get_origin();
        double zd_origin = dot(B_T.row(0), (data.row(point_idx) - origin));
        for (size_t j=0; j < data.n_rows; j++) {
            distances[j] = abs(zd_origin - dot(B_T.row(0), (data.row(j) - origin)));
        }
        double RHmin=distances.min();
        double RHmax=distances.max();

        arma::vec histNorm;
        // Build histogram for 0-dim manifold search
        if (histBins >= params.HIS_THR || params.HIS_THR == 0) {
            // generate a distances histogram
            arma::uvec hist = arma::hist(distances, histBins);

            // Cannot calculate distances
            auto filledBins = find(hist > 0).eval();
            LOG_TRACE(log) << "Non-empty bins: \n" << filledBins.n_elem;
            if ( filledBins.n_elem < 2)
                continue;

            // Threshold histogram and determine goodness of separation
            histNorm = arma::conv_to< arma::vec >::from(hist) / static_cast<double>(distances.n_elem);
        }
        // Calculate density function
        else {
            LOG_DEBUG(log) << "Start histogram bootstrapping...";
            histNorm = histBootstrapping(distances, params.HIS_THR);
            if ( histNorm.n_elem < 2)
                continue;

            origin = data.row(point_idx);
        }
        LOG_TRACE(log) << "Create histogram: \n" << histNorm.t();

        // Find threshold
        Kittler K;
        K.FindThreshold(histNorm, RHmin, RHmax);
        LOG_TRACE(log) << "Found threshold: " << K.GetDepth();

        // store separation info
        Separation sep(K.GetDiscrim(), K.GetDepth(), K.GetThreshold(),
                origin, B_T, histNorm, K.GetMinIndex(), distances);

        //#pragma omp critical(find_separation)
        if ((K.GetDiscrim()*K.GetDepth()) > best_sep.get_criteria()) {
            LOG_TRACE(log) << "Depth: " << K.GetDepth();
            LOG_TRACE(log) << "0D Origin: " << sep.get_origin();
            LOG_TRACE(log) << "Found orthogonal basis:\n" << sep.get_projection();
            best_sep=sep;
        }
    }

    if (best_sep.get_criteria() == 0)
        LOG_DEBUG(log) << "no good histograms to separate data !!!";
    else{
        LOG_DEBUG(log) <<"sep width="<<best_sep.get_width()<<"  sep depth="
            <<best_sep.get_depth() << "  sep criteria="<<best_sep.get_criteria();
        LOG_DEBUG(log)<< "his_size=" << best_sep.get_histogram().n_elem
            << ", bins=" << histBins;
    }

    return best_sep;
}


/* findManifold
 * --------------------
 * 1- calculate best separation of a manifold
 * 2- filter point that belong to found manifold
 */
bool clustering::lmclus::LMCLUS::findManifold(const arma::mat &data,
        const Parameters &params,
        arma::uvec &points_index, std::vector<unsigned int> &nonClusterPoints,
        Separation &best_sep, bool &Noise, int lm_dim)
{
    bool separated = false;
    while (true) {
        LOG_DEBUG(log) << "lm_dim: "<< lm_dim;
        // find the best fit of a set of points to a linear manifold of dimensionality lm_dim
        Separation sep;
        if (lm_dim == 0) {
            sep = findBestZeroManifoldSeparation(data.rows(points_index), params, best_sep);
        } else {
            sep = findBestSeparation(data.rows(points_index), lm_dim, params);
        }
        LOG_INFO(log) << "BEST_BOUND: "<< sep.get_criteria() << " >= " << params.BEST_BOUND;
        if (sep.get_criteria() < params.BEST_BOUND ) break;

        // find which points are close to manifold
        separated = true;
        best_sep = sep;
        std::vector<unsigned int> best_points;
        double threshold = best_sep.get_threshold();
        arma::mat Projection = best_sep.get_projection();
        arma::rowvec origin = best_sep.get_origin();
        arma::vec distances = best_sep.get_distances();

        LOG_DEBUG(log) << "Threshold: " << threshold;
        LOG_DEBUG(log) << "Origin: \n" << origin;

        LOG_TRACE(log) << "Indexes: ";
        unsigned int i, idx;
        std::vector<unsigned int> discard;
        for(i=0; i < points_index.n_rows; i++) {
            idx = points_index(i);
            //TODO: This is a temporary. See determineDistances()
            // if (lm_dim == 0)
            //     d_n = abs(zd_origin - dot(Projection.row(0), (data.row(idx) - origin)));
            // else
            //     d_n = distanceToManifold(data.row(idx) - origin, Projection);

            if (distances(i) < threshold) {    // point i has distances less than the threshold value
                best_points.push_back(idx);  // add to best points
                LOG_TRACE(log) << distances(i) << "(" << idx << ")";
            }
            else
                discard.push_back(idx);
        }

        LOG_DEBUG(log) << "Separated points: "<< best_points.size();
        if (best_points.size() < params.NOISE_SIZE) {       // small amount of points is considered noise
            Noise=true;
            LOG_INFO(log)<<"cluster size is less then noise threshold: "<<params.NOISE_SIZE<<" points";
            break;
        }
        nonClusterPoints.insert (nonClusterPoints.end(), discard.begin(), discard.end());
        points_index = arma::conv_to<arma::uvec>::from(best_points);

        LOG_DEBUG(log) << "Best basis:"<< std::endl << best_sep.get_projection();
    }

    LOG_INFO(log) << "# of points = " << points_index.n_rows;
    if(lm_dim > 0 && lm_dim < params.MAX_DIM && !Noise)
        LOG_INFO(log)<<"no separation, increasing dim ...";
    else
        LOG_INFO(log)<<"no separation ...";

    return separated;
}


/* cluster
 * --------------------
 * LMCLUS's main function:
 * 1- find manifold for specific dimension
 * 2- if could not be found increase dimension
 * 3- separate found manifold from data and begin search again
 */
void clustering::lmclus::LMCLUS::cluster(const arma::mat &data,
    const Parameters &params, std::vector<arma::uvec> &labels,
    std::vector<int> &clusterDims,
    std::vector<clustering::lmclus::Separation> &separations,
    callback_t progress) {

    // Initialize random generator
    if (params.RANDOM_SEED>0) {
        engine.seed(params.RANDOM_SEED);
    } else {
        // obtain a seed from the system clock:
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        engine.seed(seed);
    }

    // Build index
    arma::uvec points_index = arma::zeros<arma::uvec>(data.n_rows);
    for (size_t i = 0; i < points_index.n_elem; ++i)
        points_index[i] = i;

    int ClusterNum = 0;
    do {
        Separation best_sep(data.n_cols);
        bool Noise = false;
        std::vector<unsigned int> nonClusterPoints;
        int SepDim = 0;  // dimension in which separation was found
        // Bottom-up search for proper dimension
        for (int i = 1; i < params.MAX_DIM+1; i++) {
            if (findManifold(data, params, points_index, nonClusterPoints, best_sep, Noise, i))
                SepDim = i;
            if (Noise) break;
        }

        // Top-down search for embedded manifolds
        if (SepDim>0 && params.ZEROD_SEARCH) {
            LOG_DEBUG(log) << "Top-Down search d=" << SepDim-1;
            if (findManifold(data, params, points_index, nonClusterPoints, best_sep, Noise, SepDim-1)) {
                SepDim--;
                LOG_DEBUG(log) << "criteria: " << best_sep.get_criteria();
                LOG_DEBUG(log) << "origin: " << best_sep.get_origin();
                LOG_DEBUG(log) << "density: " << best_sep.get_histogram().t();
                LOG_DEBUG(log) << "distances: " << best_sep.get_distances().t();
            }
        }

        // second phase (steps 9-12)
        ClusterNum++;
        std::string msg("found cluster # "+to_string(ClusterNum)+
            ", size="+to_string(points_index.n_rows)+", dim="+to_string(SepDim));
        LOG_INFO(log)<<msg;
        if (progress != nullptr)
            progress(msg.c_str());

        clusterDims.push_back(SepDim);
        labels.push_back(points_index);

        // Realign basis of best separation
        if (params.ALIGN_BASIS){
            auto manifold  = alignBasis(data, points_index, SepDim == 0 ? 1 : SepDim);
            best_sep.set_projection(get<0>(manifold));
            best_sep.set_origin(get<1>(manifold));
            //LOG_DEBUG(log)<< "Basis after alignment:" << best_sep.get_projection();
            //LOG_DEBUG(log)<< "Origin after alignment:" << best_sep.get_origin();
        }
        separations.push_back(best_sep);

        // separate cluster points from the rest of the data
        points_index = arma::conv_to<arma::uvec>::from(nonClusterPoints);

        // Stop clustering if we reached limit
        if (ClusterNum == params.NUM_OF_CLUS) break;

    } while (points_index.n_rows > params.NOISE_SIZE);

    // Add noise as cluster
    if (points_index.n_rows > 0)
    {
        // make basis for noise
        arma::mat basis;
        arma::rowvec origin;
        // if (params.ALIGN_BASIS){
        //     auto manifold  = alignBasis(data, points_index, 1);
        //     basis = get<0>(manifold);
        //     origin = get<1>(manifold);
        // } else {
            basis = arma::zeros<arma::mat>(1, params.MAX_DIM+1);
            origin = arma::zeros<arma::rowvec>(1);
        //}

        Separation sep(0., 0., 0.,
            origin, basis,
            arma::zeros<arma::vec>(1), 0,
            arma::zeros<arma::vec>(1)
            );
        clusterDims.push_back(0);
        labels.push_back(points_index);
        separations.push_back(sep);
        LOG_INFO(log)<<"found cluster #(noise) size="<< points_index.n_rows <<" !!!";
    }

    return;
}

/* Performing PCA on cluster
    data - MxN matrix
    where  M - number of points in cluster, N - space dimensionality
*/
tuple<arma::mat, arma::rowvec> clustering::lmclus::LMCLUS::alignBasis(
    const arma::mat &data, const arma::uvec &labels, int d)
{
    arma::mat cluster = data.rows(labels);
    arma::rowvec origin = arma::mean( cluster, 0 );

    // substract mean from every point
    for (size_t i = 0; i < cluster.n_rows; ++i)
        cluster.row(i) -= origin;

    // arma::mat cov = (cluster.t()*cluster)/(cluster.n_rows-1);
    // arma::vec eigval;
    // arma::mat eigvec;
    // arma::eig_sym(eigval, eigvec, cov);
    // arma::uvec indices = sort_index(eigval, "descend");
    // LOG_DEBUG(log) << "EVAL: " << eigval;
    // LOG_DEBUG(log) << "EVEC:\n" << eigvec;
    // LOG_DEBUG(log) << "IDX:\n" << indices.rows( 0, d-1 );
    // arma::mat basis = eigvec.cols(indices.rows( 0, d-1 )).t();

    // Perform svd
    cluster /= sqrt(cluster.n_rows-1);
    arma::mat U, V;
    arma::vec s;
    arma::svd(U, s, V, cluster);
    arma::uvec indices = sort_index(s, "descend");
    LOG_DEBUG(log) << "s: " << s.t();
    LOG_DEBUG(log) << "U:\n" << U.t();
    LOG_DEBUG(log) << "V:\n" << V.t();
    arma::mat basis = V.cols( indices.rows( 0, d-1 ) ).t();

    return make_tuple(basis, origin);
}

arma::vec clustering::lmclus::LMCLUS::histBootstrapping(
    const arma::vec &distances, size_t bins)
{
    vector<size_t> counts;
    arma::vec dsorted = sort(distances);

    // count points
    size_t i = 0;
    double p = dsorted[i];
    counts.push_back(1);
    i++;
    while (i < dsorted.n_elem){
        if ((dsorted[i] - p) < arma::datum::eps){
            counts[counts.size()-1]++;
            dsorted.shed_row(i);
        } else {
            p = dsorted[i];
            i++;
            counts.push_back(1);
        }
    }
    assert(counts.size() == dsorted.n_elem);

    int S = dsorted.n_elem;
    int L = static_cast<size_t>(sqrt(S)/2.);
    if (L < 1)
        return arma::mat();

    // generate a mass function
    arma::vec emf_x = arma::vec(S-2*L);
    arma::vec emf_y = arma::vec(S-2*L);
    for (int i = L; i < S-L; i++){
        emf_x[i-L] = dsorted[i];
        size_t c = 0; // Cout all points in interval
        for (int j = -L; j <= L; j++)
            c += counts[i+j];
        emf_y[i-L] = c/(dsorted[i+L]-dsorted[i-L]);
    }

    // generate bins boundaries
    double min_d = emf_x[0], max_d = emf_x[emf_x.n_elem-1];
    arma::vec bbins = arma::vec(bins+1);
    for (size_t i = 0; i <= bins; i++){
        bbins[i] = min_d + i*(max_d-min_d)/bins;
    }

    // interpolate and integrate linear piecewise PDF
    double ilppdf, lx, ly, gx, gy;
    arma::vec lppdf = arma::vec(bins+1);
    lppdf[0] = emf_y[0];
    ilppdf = emf_y[0];
    for (size_t i=1; i<bins; i++){
        size_t tail = arma::find(emf_x < bbins[i]).eval().n_elem;
        ly = emf_y[tail-1];
        lx = emf_x[tail-1];
        gy = emf_y[tail];
        gx = emf_x[tail];
        lppdf[i] = (bbins[i] - lx)*(gy-ly)/(gx-lx)+ly;
        ilppdf += 2*lppdf[i];
        //std::cout << i << ": " << bbins[i] << " = " << lppdf[i] << "(" << ilppdf << ", " << tail << ")"<< std::endl;
    }
    lppdf[bins] = emf_y[emf_y.n_elem-1];
    ilppdf += lppdf[bins];
    ilppdf *=(max_d-min_d)/(2*bins);
    lppdf /= ilppdf;
    lppdf /= arma::sum(lppdf);

    return lppdf;
}
