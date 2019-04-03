cmake_minimum_required(VERSION 3.10)
project(tof_efficiency)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra")

find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})
include_directories(${CMAKE_SOURCE_DIR})

add_executable(tof_efficiency efficiency.cpp Access.cpp)
target_link_libraries(tof_efficiency ${ROOT_LIBRARIES})