cmake_minimum_required (VERSION 2.6)
project (Diplo)
set (CMAKE_CXX_STANDARD 11)

message("CMAKE_CXX_FLAGS_DEBUG is ${CMAKE_CXX_FLAGS_DEBUG}")
message("CMAKE_CXX_FLAGS_RELEASE is ${CMAKE_CXX_FLAGS_RELEASE}")

add_library(sequence sequence.cpp)

add_library(matchfinder matchfinder.cpp)
target_link_libraries(matchfinder sequence)

add_library(matchfinder2 matchfinder2.cpp)
target_link_libraries(matchfinder2 sequence)
target_link_libraries(matchfinder2 matchfinder)

add_library(sequencegraph sequencegraph.cpp)
target_link_libraries(sequencegraph sequence)
target_link_libraries(sequencegraph matchfinder)

add_library(utils utils.cpp)
target_link_libraries(utils matchfinder)

add_library(reconstruction reconstruction.cpp)
target_link_libraries(reconstruction sequencegraph utils)


add_library(tests tests.cpp)
target_link_libraries(tests matchfinder)

add_library(filestorage filestorage.cpp)
target_link_libraries(filestorage matchfinder matchfinder2)

add_executable(main main.cpp)
target_link_libraries(main sequence)
target_link_libraries(main matchfinder)
target_link_libraries(main matchfinder2)
target_link_libraries(main tests)
target_link_libraries(main utils)
target_link_libraries(main sequencegraph)
target_link_libraries(main reconstruction)
target_link_libraries(main filestorage)