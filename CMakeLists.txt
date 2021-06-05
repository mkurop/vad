cmake_minimum_required(VERSION 3.10)

project(vad VERSION 1.0.0)

find_package(Armadillo REQUIRED)

find_package(Boost REQUIRED)

if(${CMAKE_CLANG})
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
  set(CMAKE_CXX_COMPILER clang++)
else()
  set(CMAKE_CXX_COMPILER /usr/bin/g++)
endif()

add_subdirectory(src)
add_subdirectory(example)


