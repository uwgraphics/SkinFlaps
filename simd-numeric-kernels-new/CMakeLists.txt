cmake_minimum_required(VERSION 3.5.1)
set (CMAKE_CXX_STANDARD 11)

find_package(Eigen3)

project(simd_numeric-kernels)

include_directories(Common)
include_directories(TestDeps)
include_directories(${EIGEN3_INCLUDE_DIRS})

set (CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-qopenmp -std=gnu++11")


add_subdirectory(Tests)
