# CMakeList.txt: genetic_algorithm 的 CMake 项目，在此处包括源代码并定义
# 项目特定的逻辑。
#
cmake_minimum_required (VERSION 3.8)

project ("genetic_algorithm")

set(INCLUDEFILES ./)

aux_source_directory(./*.cpp SRC)

include_directories(./)

set(MAINFILE "test.cc")

add_executable (genetic_algorithm  ${MAINFILE} genetic.cpp message.cpp)
# add_executable(genetic_algorithm ${MAINFILE} )

