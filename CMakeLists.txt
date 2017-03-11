cmake_minimum_required (VERSION 2.6)
project (Diplo)
set (CMAKE_CXX_STANDARD 11)

add_library(sequence sequence.cpp)

add_library(matchfinder matchfinder.cpp)
target_link_libraries(matchfinder sequence)

add_library(matchfinder2 matchfinder2.cpp)
target_link_libraries(matchfinder2 sequence)
target_link_libraries(matchfinder2 matchfinder)

add_library(tests tests.cpp)
target_link_libraries(tests matchfinder)

add_library(utils utils.cpp)
target_link_libraries(utils matchfinder)

add_executable(main main.cpp)
target_link_libraries(main sequence)
target_link_libraries(main matchfinder)
target_link_libraries(main matchfinder2)
target_link_libraries(main tests)
target_link_libraries(main utils)