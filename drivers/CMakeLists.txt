

project(osqp-bench)
set(SOURCE_FILES osqp.cpp)
add_executable(osqp-bench ${SOURCE_FILES})
target_include_directories(osqp-bench PRIVATE ${OSQP_INC} ${MKL_INCLUDE_DIR})
target_link_libraries(osqp-bench osqp ${MKL_LIBRARIES})
