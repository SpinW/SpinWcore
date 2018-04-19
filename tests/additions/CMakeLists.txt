cmake_minimum_required(VERSION 2.8.2)
include(CTest)


project(TestingAdditions)
# Trivial example using gtest and gmock
add_executable(additions_test additionsTest.cpp)

target_link_libraries(additions_test
        SpinW
        gtest gmock_main
        ${ARMADILLO_LIBRARIES}
        ${MKL_LIBRARIES} pthread dl)

add_test(NAME swsym_test COMMAND swsym_test)
