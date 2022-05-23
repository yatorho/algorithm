cmake_minimum_required (VERSION 3.8)

project(ANT_COLONY_ALGORITHM)

set(CMAKE_CXX_STANDARD 14)
# src\ant_colony_algorithm.cpp
# aux_source_directory(./src SRC)
FILE(GLOB SOURCES ./src/*.cpp)


add_executable(out_put ${SOURCES})