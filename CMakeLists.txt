cmake_minimum_required(VERSION 3.19)
project(elgamal)

set(CMAKE_CXX_STANDARD 14)

add_executable(elgamal main.cpp)

target_link_libraries(elgamal gmp)