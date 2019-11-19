cmake_minimum_required(VERSION 2.8)

project(KnightKing)

set(KTK_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
option(WITH_TESTS "Build unit test programs" OFF)

include(knightking.cmake)

find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(MPI REQUIRED)
if(MPI_FOUND)
    set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
    set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
    set(MPI_EXECUTABLE_SUFFIX .openmpi)
    include_directories(${MPI_INCLUDE_PATH})
endif()

if(WITH_TESTS)
    enable_testing()
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/ext/gtest/googletest/include")
    add_subdirectory(ext/gtest)
    link_directories(${CMAKE_BINARY_DIR}/lib)
    set(GTEST_LIBRARIES "gtest" "gtest_main")
endif()

#use c++11
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${KTK_RUNTIME_OUTPUT_DIRECTORY})

include_directories(include)
include_directories(ext)

add_subdirectory(src)
