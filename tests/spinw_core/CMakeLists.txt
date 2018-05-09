cmake_minimum_required(VERSION 2.8.2)
include(CTest)


project(TestingSpinw_core)
# Trivial example using gtest and gmock
add_executable(sw_core
        swcore_atom.cpp
        swcore_intmatrix.cpp)

target_link_libraries(sw_core
        gtest gmock_main
        SpinW
        ${ARMADILLO_LIBRARIES}
        ${MKL_LIBRARIES} pthread dl)

add_test(NAME swcore COMMAND sw_core)
