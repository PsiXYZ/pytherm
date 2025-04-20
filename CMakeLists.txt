# cmake_minimum_required(VERSION 3.30)
cmake_minimum_required(VERSION 3.15...3.26)

project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(PYBIND11_NEWPYTHON ON)

message("zzzzzzzzzzzz")
message("${CMAKE_BUILD_TYPE}")
# set(PYTHON_IS_DEBUG ON)

# find_package(Python 3.13 COMPONENTS Interpreter Development REQUIRED)
# find_package(pybind11 CONFIG REQUIRED)

# if(CMAKE_BUILD_TYPE STREQUAL "Debug")
# string(TOUPPER "${CMAKE_BUILD_TYPE}" uppercase_CMAKE_BUILD_TYPE)
# if("${uppercase_CMAKE_BUILD_TYPE}" MATCHES DEBUG|RELWITHDEBINFO)
#     message("debug mode")
#     set(PYTHON_IS_DEBUG ON)
# endif()

# find_package(Python COMPONENTS Interpreter Development)
# find_package(pybind11 CONFIG)
# set_target_properties(Python::Module PROPERTIES  MAP_IMPORTED_CONFIG_DEBUG ";RELEASE")


add_subdirectory(externals/pybind11)
add_subdirectory(src/pytherm)
add_subdirectory(src/binds)