project(lmclus)
cmake_minimum_required (VERSION 2.8.0)

message(STATUS "Configuring LMCLUS")

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules/")

# gcc flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -Wall -Wno-deprecated -O2")

include_directories(inst/include/armadillo)
include_directories(inst/include/cpplog)
include_directories(inst/include/)

if(NOT WIN32)
    find_package(MKL)
    if (MKL_FOUND)
        set(BLAS_LIBS ${MKL_LIBRARIES})
    else(MKL_FOUND)
        find_package(BLAS REQUIRED)
        message(STATUS "Found BLAS libraries: ${BLAS_LIBRARIES}")
        set(BLAS_LIBS ${BLAS_LIBRARIES})
    endif(MKL_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
endif(NOT WIN32)

if(WIN32)
	link_directories(${CMAKE_SOURCE_DIR}/contrib/libs)
	set(BLAS_LIBS blas lapack)
endif(WIN32)

#find_package(Armadillo REQUIRED)
#include_directories(${ARMADILLO_INCLUDE_DIRS})
find_package(R REQUIRED)
include_directories(${R_INCLUDE_DIR})

option(LMCLUS_BUILD_TESTING "Build the test suite." ON)

if (LMCLUS_BUILD_TESTING)
  enable_testing()
endif (LMCLUS_BUILD_TESTING)

# include sources
add_subdirectory(src)
if (LMCLUS_BUILD_TESTING)
  add_subdirectory (tests)
endif (LMCLUS_BUILD_TESTING)
