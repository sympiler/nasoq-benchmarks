cmake_minimum_required(VERSION 3.9)

set(CMAKE_CXX_STANDARD 11)
enable_language(C)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

find_package(MKL REQUIRED)

add_subdirectory(nasoq)
add_subdirectory(osqp)

include_directories(${CMAKE_CURRENT_LIST_DIR}/nasoq/common/ ${CMAKE_CURRENT_LIST_DIR}/smp-format/)
set(OSQP_INC osqp/include/)

add_subdirectory(drivers)

