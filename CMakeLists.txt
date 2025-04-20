cmake_minimum_required(VERSION 3.15...3.26)

project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(PYBIND11_NEWPYTHON ON)

add_subdirectory(externals/pybind11)
add_subdirectory(src/pytherm)
add_subdirectory(src/binds)