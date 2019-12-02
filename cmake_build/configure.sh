#!/bin/sh

# this is an example file...please DO NOT MODIFY if you don't want to do it for everyone
# to use it, copy it to another file and run it

# additional compiler flags could be added customizing the corresponding var, for example
# -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -msse3 ". Note that the "defaults are already correctly coded"
# so we should add here only machine specific stuff

# an effort is made to autodetect all of the libraries needed HOWEVER:
# METIS_APPLICATION needs the var PARMETIS_ROOT_DIR to be specified by the user (not needed if the app is set to OFF)
# TRILINOS_APPLICATION needs the var PARMETIS_ROOT_DIR to be specified by the user (not needed if the app is set to OFF)

#the user should also note that the symbol "\" marks that the command continues on the next line. IT SHOULD ONLY BE FOLLOWED
#BY the "ENTER" and NOT by any space!!

clear

# you may want to decomment this the first time you compile
rm CMakeCache.txt
rm *.cmake
rm -rf CMakeFiles\

cmake ..                                                                            \
-DCMAKE_C_COMPILER=/usr/bin/gcc                                                     \
-DCMAKE_CXX_COMPILER=/usr/bin/g++                                                   \
-DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -msse3 -std=c++11 "                           \
-DCMAKE_C_FLAGS="${CMAKE_C_FLAGS} -msse3 "                                          \
-DBOOST_ROOT="/home/bsaridar/Software/boost_1_69_0"                     		    \
-DPYTHON_EXECUTABLE="/usr/bin/python3"                                              \
-DCMAKE_BUILD_TYPE=FullDebug                                                          \
-DMESHING_APPLICATION=ON                                                            \
-DEXTERNAL_SOLVERS_APPLICATION=ON                                                   \
-DSTRUCTURAL_MECHANICS_APPLICATION=ON                                               \
-DCONVECTION_DIFFUSION_APPLICATION=OFF                                              \
-DSOLID_MECHANICS_APPLICATION=OFF                                                   \
-DCONSTITUTIVE_MODELS_APPLICATION=OFF                                               \
-DFLUID_DYNAMICS_APPLICATION=ON                                                     \
-DMESH_MOVING_APPLICATION=ON                                                        \
-DFSI_APPLICATION=OFF                                                               \
-DKRATOS_BUILD_TESTING=ON                                                           \
-DDEM_APPLICATION=OFF                                                               \
-DSWIMMING_DEM_APPLICATION=OFF                                                      \
-DMIXED_ELEMENT_APPLICATION=OFF                                                     \
-DSHAPE_OPTIMIZATION_APPLICATION=OFF                                                \
-DTOPOLOGY_OPTIMIZATION_APPLICATION=OFF                                             \
-DMETIS_APPLICATION=ON                                                             \
-DTRILINOS_INCLUDE_DIR="/usr/include/trilinos"                                      \
-DTRILINOS_LIBRARY_DIR="/usr/lib/x86_64-linux-gnu"                                  \
-DTRILINOS_LIBRARY_PREFIX="trilinos_"                                               \
-DUSE_METIS_5=ON                                                                    \
-DPARMETIS_ROOT_DIR="/home/youruser/compiled_libraries/ParMetis-3.1.1"              \
-DTRILINOS_APPLICATION=OFF                                                          \
-DTRILINOS_ROOT="/home/youruser/compiled_libraries/trilinos-10.2.0"                 \
-DEMPIRE_APPLICATION=ON 							    \
-DMAPPING_APPLICATION=ON 							    \
-DCOMPRESSIBLE_POTENTIAL_FLOW_APPLICATION=ON 					    \
-DINSTALL_EMBEDDED_PYTHON=ON                                         \
-DCO_SIMULATION_APPLICATION=ON


# decomment this to have it verbose
# make VERBOSE=1 -j4
make -j4
make install
