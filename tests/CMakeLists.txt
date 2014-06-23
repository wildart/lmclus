# add some includes
include_directories(../include/third-party)
# build executable
add_executable(lmclus_test lmclus_test.cpp)
add_executable(lmclus_la_test lmclus_la_test.cpp)
# link test with library
target_link_libraries(lmclus_test lmclus ${LAPACK_LIBRARIES})
target_link_libraries(lmclus_la_test lmclus ${LAPACK_LIBRARIES})
# run test
add_test(lmclus_la_test lmclus_la_test)
add_test(lmclus_test lmclus_test -c ${CMAKE_SOURCE_DIR}/data/lmclus.ini -i ${CMAKE_SOURCE_DIR}/data/testData)
