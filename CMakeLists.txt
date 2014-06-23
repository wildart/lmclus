project(lmclus)
cmake_minimum_required (VERSION 2.8.0)

message(STATUS "Configuring LMCLUS")

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules/")

# gcc flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall ")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wno-deprecated -O3")


include_directories(include/armadillo)
include_directories(include/cpplog)
include_directories(include/)

if(NOT WIN32)
    find_package(MKL)
    if (MKL_FOUND)
        set(BLAS_LIBS ${MKL_LIBRARIES})
    else(MKL_FOUND)
        find_package(LAPACK REQUIRED)
        message(STATUS "Found LAPACK libraries: ${LAPACK_LIBRARIES}")
    endif(MKL_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
endif(NOT WIN32)

if(WIN32)
	link_directories(${CMAKE_SOURCE_DIR}/contrib/libs)
	set(BLAS_LIBS blas lapack)
endif(WIN32)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

option(LMCLUS_BUILD_TESTING "Build the test suite." ON)

if (LMCLUS_BUILD_TESTING)
  enable_testing()
endif (LMCLUS_BUILD_TESTING)

# include sources
add_subdirectory(src)
if (LMCLUS_BUILD_TESTING)
  add_subdirectory (tests)
endif (LMCLUS_BUILD_TESTING)
