cmake_minimum_required(VERSION 3.22)
project(DLA_cpp)

set(CMAKE_CXX_STANDARD 20)

add_executable(DLA_cpp main.cpp convex_hull.cpp point.h convex_hull.h)
