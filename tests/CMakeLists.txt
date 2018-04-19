cmake_minimum_required(VERSION 2.8.2)
include(CTest)

if (CMAKE_VERSION VERSION_LESS 3.2)
    set(UPDATE_DISCONNECTED_IF_AVAILABLE "")
else()
    set(UPDATE_DISCONNECTED_IF_AVAILABLE "UPDATE_DISCONNECTED 1")
endif()

include(DownloadProject.cmake)
download_project(
        PROJ                googletest
        GIT_REPOSITORY      https://github.com/google/googletest.git
        GIT_TAG             master
        ${UPDATE_DISCONNECTED_IF_AVAILABLE}
        )

# Prevent GoogleTest from overriding our compiler/linker options
# when building with Visual Studio
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})

# When using CMake 2.8.11 or later, header path dependencies
# are automatically added to the gtest and gmock targets.
# For earlier CMake versions, we have to explicitly add the
# required directories to the header search path ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
    include_directories("${gtest_SOURCE_DIR}/include"
            "${gmock_SOURCE_DIR}/include")
endif()

# Use gtest and gmock
if (MAKE_IDIVIDUAL_TESTS)
    add_subdirectory(templates)
    add_subdirectory(swsym)
    add_subdirectory(additions)
else()
    project(Testing)
    add_executable(all_test
            templates/templateTest.cpp
            swsym/swsymTest.cpp
            additions/additionsTest.cpp
            )
    target_link_libraries(all_test
            gtest gmock_main
            SpinW
            ${ARMADILLO_LIBRARIES}
            ${MKL_LIBRARIES} pthread dl)
    add_test(NAME all_test COMMAND all_test)
endif()
