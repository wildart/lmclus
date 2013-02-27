# add some includes
include_directories(../inst/include/third-party)
# build executable
add_executable(lmclus_test lmclus_test.cpp)
# link demo with library
target_link_libraries(lmclus_test lmclus ${BLAS_LIBS})
# run test
add_test(lmclus_test lmclus_test -c lmclus.ini -i testData)
