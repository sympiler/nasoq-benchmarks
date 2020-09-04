cmake_minimum_required(VERSION 3.9)

set(CMAKE_CXX_STANDARD 11)
enable_language(C)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

find_package(MKL REQUIRED)

find_package(OpenMP)

set(MOSEK_INC "$ENV{MOSEKROOT}/h")
set(MOSEK_LIB "$ENV{MOSEKROOT}/bin")

set(GUROBI_INC "$ENV{GUROBIROOT}/include")
set(GUROBI_LIB "$ENV{GUROBIROOT}/lib")


add_subdirectory(nasoq)
# add_subdirectory(osqp)


include_directories(${CMAKE_CURRENT_LIST_DIR}/nasoq/common/ ${CMAKE_CURRENT_LIST_DIR}/smp-format/)
set(OSQP_INC osqp/include/)

add_subdirectory(drivers)
