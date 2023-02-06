cmake_minimum_required(VERSION 3.22)
project(DLA_cpp)

set(CMAKE_CXX_STANDARD 20)

find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})

add_executable(DLA_cpp main.cpp convex_hull.cpp point.h convex_hull.h)
target_link_libraries(DLA_cpp PUBLIC ${OpenCV_LIBS})
