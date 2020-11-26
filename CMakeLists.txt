cmake_minimum_required(VERSION 3.17)
project(rebin)

find_package(ROOT 6.22 REQUIRED COMPONENTS)
include(${ROOT_USE_FILE}) # Enables CLion to find root headers

set(CMAKE_CXX_STANDARD 17)

add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME} ${ROOT_LIBRARIES})