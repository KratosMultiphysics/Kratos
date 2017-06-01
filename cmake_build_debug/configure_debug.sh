#!/bin/sh

# this is an example file...please DO NOT MODIFY if you don't want to do it for everyone
#to use it, copy it to another file and run it

# additional compiler flags could be added customizing the corresponding var, for example
# -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -msse3 ". Note that the "defaults are already correctly coded"
#so we should ass here only machine specific stuff

#an effort is made to autodetect all of the libraries needed HOWEVER:
#METIS_APPLICATION needs the var PARMETIS_ROOT_DIR to be specified by the user (not needed if the app is set to OFF)
#TRILINOS_APPLICATION needs the var PARMETIS_ROOT_DIR to be specified by the user (not needed if the app is set to OFF)
#MKL_SOLVERS_APPLICATION needs the var MKLSOLVER_INCLUDE_DIR and MKLSOLVER_LIB_DIR to be specified by the user (not needed if the app is set to OFF)
#note that the MKLSOLVER_LIB_DIR should include /lib/em64t. This is needed as intel is changing location of mkl at every update of the compiler!!

#the user should also note that the symbol "\" marks that the command continues on the next line. IT SHOULD ONLY BE FOLLOWED
#BY the "ENTER" and NOT by any space!!

clear

#you may want to decomment this the first time you compile
rm CMakeCache.txt
rm *.cmake
rm -rf CMakeFiles\

cmake ..                            						                               	                \
-DCMAKE_C_COMPILER=/usr/bin/gcc                         			                   		                \
-DCMAKE_INSTALL_RPATH="/home/rzorrilla/Kratos/libs"                              			                        \
-DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE                                        			                        \
-DCMAKE_CXX_COMPILER=/usr/bin/g++						                                                \
-DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -msse3 -std=c++11" 					                                \
-DCMAKE_C_FLAGS="${CMAKE_C_FLAGS} -msse3 " 					                            	                \
-DBOOST_ROOT="/home/rzorrilla/Compiled_libraries/boost_1_59_0"				                                        \
-DPYTHON_LIBRARY="/usr/lib/python3.5/config-3.5m-x86_64-linux-gnu/libpython3.5m.so" 		                                \
-DPYTHON_INCLUDE_DIR="/usr/include/python3.5"                                       		                                \
-DCMAKE_BUILD_TYPE=Debug  							                                   	        \
-DKRATOS_DEBUG=ON														\
-DDO_NOT_DEFINE_NDEBUG=ON													\
-DINCOMPRESSIBLE_FLUID_APPLICATION=OFF 						                              	                \
-DMESHING_APPLICATION=ON 							                                                \
-DEXTERNAL_SOLVERS_APPLICATION=ON						                                                \
-DPFEM_APPLICATION=OFF 								                                                \
-DSTRUCTURAL_APPLICATION=OFF 							                                                \
-DCONVECTION_DIFFUSION_APPLICATION=ON 						                              	                \
-DFLUID_DYNAMICS_APPLICATION=ON 						                                                \
-DALE_APPLICATION=ON 								                                                \
-DFSI_APPLICATION=ON 								                                                \
-DMAPPING_APPLICATION=ON                                                                                                        \
-DOPENCL_APPLICATION=OFF     						                                                        \
-DMIXED_ELEMENT_APPLICATION=OFF							                                                \
-DMKL_SOLVERS_APPLICATION=OFF							                                  	        \
-DMKLSOLVER_INCLUDE_DIR="/opt/intel/Compiler/11.1/072/mkl/include"	            			                        \
-DMKLSOLVER_LIB_DIR="/opt/intel/Compiler/11.1/072/mkl/lib/em64t"	              			                        \
-DMETIS_APPLICATION=ON								                                                \
-DMETIS_INCLUDE_DIR="/usr/include/"                                                                                             \
-DUSE_METIS_5=ON                                                                                                                \
-DPARMETIS_ROOT_DIR="/usr/lib/"                                                   			                        \
-DTRILINOS_APPLICATION=ON							                                                \
-DTRILINOS_LIBRARY_DIR="/usr/lib/x86_64-linux-gnu/"                                                                             \
-DTRILINOS_INCLUDE_DIR="/usr/include/trilinos/"                                                                                 \
-DTRILINOS_LIBRARY_PREFIX="trilinos_"                                                                                           \
-DINCLUDE_MMG=ON                                                                                                                \
-DMMG_INCLUDE_DIR="/home/rzorrilla/Compiled_libraries/mmg/include/"                                                             \
-DMMG2D_INCLUDE_DIR="/home/rzorrilla/Compiled_libraries/mmg/include/mmg/mmg2d/"                                                 \
-DMMG3D_INCLUDE_DIR="/home/rzorrilla/Compiled_libraries/mmg/include/mmg/mmg3d/"                                                 \
-DMMGS_INCLUDE_DIR="/home/rzorrilla/Compiled_libraries/mmg/include/mmg/mmgs/"                                                   \
-DMMG_LIBRARY="/home/rzorrilla/Compiled_libraries/mmg/lib/libmmg.a"                                                             \
-DMMG2D_LIBRARY="/home/rzorrilla/Compiled_libraries/mmg/lib/libmmg2d.a"                                                         \
-DMMG3D_LIBRARY="/home/rzorrilla/Compiled_libraries/mmg/lib/libmmg3d.a"                                                         \
-DMMGS_LIBRARY="/home/rzorrilla/Compiled_libraries/mmg/lib/libmmgs.a"                                                           \
-DDEM_APPLICATION=OFF								                                                \
-DSWIMMING_DEM_APPLICATION=OFF                                             					                \
-DSOLID_MECHANICS_APPLICATION=ON   				                                    		                \
-DSTRUCTURAL_MECHANICS_APPLICATION=ON 				                                		                \
-DCOMPRESSIBLE_POTENTIAL_FLOW_APPLICATION=ON 				                                		        \
-DCONTACT_STRUCTURAL_MECHANICS_APPLICATION=OFF		                                		                        \
-DPFEM_SOLID_MECHANICS_APPLICATION=OFF 				                                		                \
-DFSI_TRILINOS_APPLICATION=OFF                                                                                                  \
-DINSTALL_EMBEDDED_PYTHON=ON                                               					                \
-DKRATOS_BUILD_TESTING=ON                                                                                                       \

#decomment this to have it verbose
# make VERBOSE=1 -j4
make -j4 install
