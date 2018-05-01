cmake_minimum_required(VERSION 3.6)

project(SpinW)
set(DO_TESTS "TRUE")
set(MAKE_IDIVIDUAL_TESTS "FALSE")

set(CMAKE_CXX_STANDARD 11)
#Make highly optimised
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# I've had problems with gcc < 5.3 and OpenMP does not work correctly.
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # require at least gcc 4.8
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8)
        message(FATAL_ERROR "GCC version must be at least 4.8, but might not work with MATLAB.")
    endif()
        if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.3)
        message(WARNING "GCC version must be at least 5.3 to be fully compatible!")
    endif()
else()
    message(WARNING "You are using an unsupported compiler! Compilation has only been tested with GCC.")
endif()

# Armadillo is the base package so it is required
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

# If we want threading then we must link with OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# We use Intel MKL as a numeric library but fall back to LAPACK and BLAS if this fails.
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
find_package(MKL)
#if (MKL_FOUND)
include_directories(${MKL_INCLUDE_DIR})
#    set(LAPACK_INCLUDE_DIRS '')
#else()
#    # This should find the BLAS libraries as well.
#    find_package(LAPACK REQUIRED)
#    include_directories(${LAPACK_INCLUDE_DIRS})
#endif()

# Make the shared library
add_library(
        SpinW SHARED
        src/spinw.cpp include/spinw.h
        src/atom.cpp include/atom.h
        src/intmatrix.cpp include/intmatrix.h
        src/Hspinw.cpp include/Hspinw.h # The main files and C interface
        include/sw_structs.h # Helpful structures
        src/sw_additions.cpp include/sw_additions.h
        src/spinw_external.cpp # Helpful functions
        include/symmetryDat.h # Symmetry files
        swsym/src/swsym.cpp swsym/include/swsym.h
        swsym/src/position.cpp swsym/include/position.h
        swsym/src/bond.cpp swsym/include/bond.h
        swsym/src/symOperator.cpp swsym/include/symOperator.h
        swsym/src/oporder.cpp swsym/include/oporder.h
        src/templateFuncs.tpp include/templateFuncs.h # Template functions
        )

# Link with the required libraries.
target_link_libraries(SpinW
        ${ARMADILLO_LIBRARIES}
        ${MKL_LIBRARIES} pthread dl)

# We can add tests here...
if(DO_TESTS)
    add_subdirectory(tests)
endif()