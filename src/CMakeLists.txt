# build library
add_library(lmclus Kittler.cpp lmclus.cpp)

# build R library
find_package(R REQUIRED)
include_directories(${R_INCLUDE_DIR})
add_library(lmclusr SHARED lmclusr.cpp)
target_link_libraries(lmclusr lmclus ${BLAS_LIBS} ${R_LIBRARIES})