CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(affine-g2o)

set( CMAKE_BUILD_TYPE "Release" )
set( CMAKE_CXX_FLAGS "-std=c++11 -O3" )

list( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules )

# Eigen
include_directories( "/usr/include/eigen3" "./include" )
#include_directories

find_package( G2O REQUIRED )
include_directories( ${G2O_INCLUDE_DIRS} )

find_package( Cholmod REQUIRED )
include_directories( ${CHOLMOD_INCLUDE_DIR} )

add_executable( affine_g2o ../src/affine_g2o.cpp )
target_link_libraries( affine_g2o
    g2o_core g2o_stuff 
    ${CHOLMOD_LIBRARIES}
)

add_executable( sim_generate ../src/simulation_data_creation.cpp )
target_link_libraries( sim_generate
    g2o_core g2o_stuff 
    ${CHOLMOD_LIBRARIES}
)
