cmake_minimum_required(VERSION 3.17)

project(StraightECMC)

set(CMAKE_CXX_STANDARD 11)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pedantic-errors -Wall -Wextra -Wconversion -Wsign-conversion -Wshadow")
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release) # Use CMAKE_CXX_FLAGS_RELEASE for compilation (usually -O3 -DNDEBUG).
endif()

find_package(Boost REQUIRED)
include_directories(SYSTEM ${Boost_INCLUDE_DIR})

# Set options for the HighFive package (see HighFive/CMakeLists.txt)
set(USE_BOOST ON)
set(USE_EIGEN OFF)
set(USE_XTENSOR OFF)
set(USE_OPENCV OFF)
set(HIGHFIVE_UNIT_TESTS OFF)
set(HIGHFIVE_EXAMPLES OFF)
set(HIGHFIVE_PARALLEL_HDF5 OFF)
set(HIGHFIVE_BUILD_DOCS OFF)
add_subdirectory(HighFive)

add_executable(StraightEcmc main.cc
        cell.cc
        StraightEcmc.h
        straightEcmc.cc
        run.cc
        init.cc
        inout.cc
        Param.h
        refresh.cc)

target_link_libraries(StraightEcmc HighFive)
