cmake_minimum_required(VERSION 3.9)
project(SpinW)

set(CMAKE_CXX_STANDARD 11)

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

#find_package(LAPACK REQUIRED)
#include_directories(${LAPACK_INCLUDE_DIRS})
#set(PROJECT_LIBS "${PROJECT_LIBS} ${LAPACK_LIBRARIES}")


find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
find_package(MKL REQUIRED)
include_directories(${MKL_INCLUDE_DIRS})

add_library(SpinW SHARED spinw.cpp spinw.h Hspinw.cpp Hspinw.h sw_helper.h sw_additions.cpp sw_additions.h spinw_external.cpp)
target_link_libraries(SpinW ${ARMADILLO_LIBRARIES} ${MKL_LIBRARIES} pthread dl)

