cmake_minimum_required (VERSION 2.6)
project (Diplo)
set (CMAKE_CXX_STANDARD 11)

add_library(sequence sequence.cpp)

add_library(matchfinder matchfinder.cpp)
target_link_libraries(matchfinder sequence)

add_executable(main main.cpp)
target_link_libraries(main sequence)
target_link_libraries(main matchfinder)