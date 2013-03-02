# build library
add_library(lmclus Kittler.cpp lmclus.cpp)

# build R library
add_library(lmclusr SHARED lmclusr.cpp)

# link with library
target_link_libraries(lmclusr lmclus ${BLAS_LIBS} ${R_LIBRARIES})