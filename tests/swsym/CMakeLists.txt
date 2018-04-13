cmake_minimum_required(VERSION 2.8.2)
include(CTest)


project(TestingSym)
# Trivial example using gtest and gmock
add_executable(swsym_test swsymTest.cpp)

target_link_libraries(swsym_test
        gtest gmock_main
        ${ARMADILLO_LIBRARIES}
        ${MKL_LIBRARIES} pthread dl)

add_test(NAME swsym_test COMMAND swsym_test)
