cmake_minimum_required(VERSION 3.28)
project(coarsing)

set(CMAKE_CXX_STANDARD 17)

add_executable(coarsing main.cpp
        src/coarsing.cpp)
