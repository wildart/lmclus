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
#include <cstdlib>
#include <string>

#define CFG "config"
#define IN "input"
#define OUT "output"
#define AUX "aux"
#define DIST "dist"
#define SEP "separator"
#define LOGF "log"
#define SKIP "skip"
#include "cmdline.hpp"

#define SI_SUPPORT_IOSTREAMS
#define SECTION "lmclus"
#include "SimpleIni.h"

#include "lmclus.hpp"

#define TEST(log, cond, msg) if(!(cond)){ LOG_ERROR(log) << msg; return EXIT_FAILURE;}

void specifyCommandLineParameters(cmdline::parser &cmdOpt)
{
    cmdOpt.add<std::string>(CFG, 'c', "configuration file", true, "");
    cmdOpt.add<std::string>(IN, 'i', "dataset input file", true, "");
    cmdOpt.add<std::string>(OUT, 'o', "cluster output file", false, "");
    cmdOpt.add<std::string>(LOGF, 'l', "log file", false, ""); 
    cmdOpt.add<char>(SEP, 's', "separator in output file", false, ','); 
    cmdOpt.add(SKIP, 'S', "skip faulty records");   
}
    
int main ( int argc, char *argv[] )
{
    cmdline::parser cmdOpt;
    specifyCommandLineParameters(cmdOpt);
    if (!cmdOpt.parse(argc, argv)) {
        std::cout << cmdOpt.error()<< std::endl << cmdOpt.usage() << std::endl;
        return EXIT_FAILURE;
    }
 
    // reading command line
    std::string parametersFile(cmdOpt.get<std::string>(CFG));
    std::string inputDataFile(cmdOpt.get<std::string>(IN));
    
    // Setup logger
    cpplog::OstreamLogger *log;
    if (cmdOpt.exist(LOGF))
    {
      log = new cpplog::FileLogger(cmdOpt.get<std::string>(LOGF)+".log", true);  
    }
    else 
    {
      log = new cpplog::StdErrLogger();
    }
    
    // load lmclus parameters from parameters file 
    CSimpleIniA ini(false, false, false);
    SI_Error rc = ini.LoadFile(parametersFile.c_str());
    if (rc < 0){
        LOG_ERROR(log) << "Problem reading parameters file: " << parametersFile;
        return EXIT_FAILURE;
    }
    
    // Set parameters from ini-file
    Parameters params;
    params.MAX_DIM = static_cast<int>(ini.GetLongValue(SECTION, "MAX_DIM", 2));
    params.NUM_OF_CLUS = static_cast<int>(ini.GetLongValue(SECTION, "NUM_OF_CLUS", 2));
    params.LABELED_DATA = static_cast<int>(ini.GetLongValue(SECTION, "LABELED_DATA", 1));
    params.CONST_SIZE_HIS = static_cast<int>(ini.GetLongValue(SECTION, "CONST_SIZE_HIS", 0));
    params.NOISE_SIZE = static_cast<unsigned int>(ini.GetLongValue(SECTION, "NOISE_SIZE", 2));
    params.BEST_BOUND = ini.GetDoubleValue(SECTION, "BEST_BOUND", 1.0);
    params.ERROR_BOUND = ini.GetDoubleValue(SECTION, "ERROR_BOUND", 0.0001);
    params.MAX_BIN_PORTION = ini.GetDoubleValue(SECTION, "MAX_BIN_PORTION", 0.1);
    
    params.RANDOM_SEED = static_cast<unsigned int>(ini.GetLongValue(SECTION, "RANDOM_SEED", 0));
    params.SAMPLING_HEURISTIC = static_cast<int>(ini.GetLongValue(SECTION, "SAMPLING_HEURISTIC", 3));
    params.SAMPLING_FACTOR = ini.GetDoubleValue(SECTION, "SAMPLING_FACTOR", 0.003);
    params.HIS_SAMPLING = ini.GetBoolValue(SECTION, "HIS_SAMPLING", false);
    params.SAVE_RESULT = ini.GetBoolValue(SECTION, "SAVE_RESULT", false);
    
    // Load dataset
    arma::mat data;
    data.load(inputDataFile, arma::csv_ascii);
    //data = arma::randu<arma::mat>(100,11);
    
    LOG_INFO(log) << "Parameters file: " << parametersFile;
    LOG_INFO(log) << "Parameters: " << params;
    LOG_INFO(log) << "Input: " << inputDataFile;
    LOG_INFO(log) << "Input rows: " << data.n_rows;
    LOG_INFO(log) << "Input cols: " << data.n_cols;
    
    TEST(log, data.n_cols == 11, "Invalid column number")
    TEST(log, data.n_rows == 3000, "Invalid row number")
    
    auto original_labels = data.col(data.n_cols-1);
    data.shed_col(data.n_cols-1);
    TEST(log, data.n_cols == 10, "Invalid column number")
    
    TEST(log, params.MAX_DIM < static_cast<int>(data.n_cols), "Linear manifold dimension must be less than the dimension of the data !!!")

    // Process data
    LMCLUS lmclus(log);
    std::vector<arma::uvec> labels;
    std::vector<double> thresholds; 
    std::vector<arma::mat> basises; 
    std::vector<int> clusterDims;
    lmclus.cluster(data, params, labels, thresholds, basises, clusterDims);
    
    // Process results
    size_t clusterNum = labels.size();
    TEST(log, clusterDims.size() == clusterNum, "Invalid number of clusters")
    for(size_t i=0; i < clusterNum; ++i)
    {
        LOG_INFO(log) << "Cluster " << i << " dimension: " << clusterDims[i] << ", size: " << labels[i].n_elem;
    }
    
    LOG_INFO(log) << "labels found: " << labels.size();
    LOG_INFO(log) << "thresholds found: " << thresholds.size();
    LOG_INFO(log) << "basises found: " << basises.size();
    LOG_INFO(log) << "clusterDims found: " << clusterDims.size();
    
    // Clear log
    delete log;
    
    return EXIT_SUCCESS;
}
