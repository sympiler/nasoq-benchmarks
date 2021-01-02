project(osqp-bench)
set(SOURCE_FILES osqp.cpp)
add_executable(osqp-bench ${SOURCE_FILES})
target_include_directories(osqp-bench PRIVATE ${OSQP_INC} ${MKL_INCLUDE_DIR})
target_link_libraries(osqp-bench osqp ${MKL_LIBRARIES})


if(EXISTS ${GUROBI_INC})
project(gurobi-bench)
set(GR_SOURCE_FILES gurobi.cpp Gurobi/Gurobi.cpp)
ADD_DEFINITIONS(-DGRB_LICENSE_FILE="$ENV{HOME}/gurobi.lic")
add_executable(gurobi-bench ${GR_SOURCE_FILES})
target_include_directories(gurobi-bench PRIVATE ${GUROBI_INC}  "Gurobi/"
        ${EIGEN_INC}
        ##FIXME MKL include is added for a mess in header files but is not required here!
        "../nasoq/QP/" ${MKL_INCLUDE_DIR} "../nasoq/matrixVector/")
target_link_directories(gurobi-bench PRIVATE ${GUROBI_LIB} )

target_link_libraries(gurobi-bench
        gurobi_g++5.2
        gurobi81
        GurobiJni81
        OpenMP::OpenMP_CXX
        )
endif()


if(EXISTS ${MOSEK_INC})
project(mosek-bench)
set(MSK_SOURCE_FILES mosek.cpp Mosek/QuadProgMosek.cpp)
ADD_DEFINITIONS(-DkMosekLicensePath="$ENV{HOME}/mosek/mosek.lic")
add_executable(mosek-bench ${MSK_SOURCE_FILES})
target_include_directories(mosek-bench PRIVATE ${MOSEK_INC} "Mosek/"
         ${EIGEN_INC}
##FIXME MKL include is added for a mess in header files but is not required here!
        "../nasoq/QP/" ${MKL_INCLUDE_DIR} "../nasoq/matrixVector/")
target_link_directories(mosek-bench PRIVATE ${MOSEK_LIB}  )
target_link_libraries(mosek-bench
        mosek64
        )
endif()