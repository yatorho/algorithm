cmake_minimum_required(VERSION 3.10)

project(PARTICLE)

aux_source_directory("pco/src" SRC)

include_directories("pco/include")


add_executable(PCO_OUT ${SRC} test.cc)
# add_executable(PCO_OUT main.cpp)
