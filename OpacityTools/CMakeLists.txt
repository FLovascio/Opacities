#support for apple silicon is mandated
cmake_minimum_required(VERSION 3.19) 

# set the project name and version
project(OpacityTools VERSION 0.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
add_subdirectory(src)