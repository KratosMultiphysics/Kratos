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


clear all
clear

#you may want to decomment this the first time you compile
rm CMakeCache.txt

cmake ..  									\
-DCMAKE_C_COMPILER=icc								\
-DCMAKE_CXX_COMPILER=icpc						        \
-DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -msse2 " 					\
-DCMAKE_C_FLAGS="${CMAKE_C_FLAGS} -msse2 " 					\
-DCMAKE_BUILD_TYPE=Release							\
-DINCOMPRESSIBLE_FLUID_APPLICATION=ON  				                \
-DMESHING_APPLICATION=ON 							\
-DEXTERNAL_SOLVERS_APPLICATION=ON						\
-DPFEM_APPLICATION=OFF 								\
-DDEM_APPLICATION=OFF 								\
-DSTRUCTURAL_APPLICATION=ON 							\
-DFEM_DEM_APPLICATION=ON 							\
-DCONSTITUTIVE_LAWS_APPLICATION=OFF						\
-DTHERMO_MECHANICAL_APPLICATION=OFF 						\
-DCONVECTION_DIFFUSION_APPLICATION=OFF 						\
-DFLUID_DYNAMICS_APPLICATION=OFF 						\
-DALE_APPLICATION=OFF 								\
-DFSI_APPLICATION=OFF 								\
-DOPENCL_APPLICATION=OFF							\
-DMPI_SEARCH_APPLICATION=OFF 							\
-DTURBULENT_FLOW_APPLICATION=OFF						\
-DMIXED_ELEMENT_APPLICATION=ON							\
-DMKL_SOLVERS_APPLICATION=ON							\
-DMKLSOLVER_INCLUDE_DIR="/home/nelson/Compilers/mkl/include"                    \
-DMKLSOLVER_LIB_DIR="/home/nelson/Compilers/mkl/lib/em64t"                      \
-DMETIS_APPLICATION=ON								\
-DPARMETIS_ROOT_DIR="/home/nelson/compiled_libraries/parmetis-4.0.3" 		\
-DTRILINOS_APPLICATION=ON							\
-DTRILINOS_ROOT="/usr/local/trilinos-11.2.3"                                    \
#-DTRILINOS_LIB_DIRS="/usr/local/trilinos-11.2.3/lib"                            \
#-DTRILINOS_INC_DIRS="/usr/local/trilinos-11.2.3/include"                        \


#-DUSE_INTEL_GREATER_THAN_12=TRUE                                                
#decomment this to have it verbose
# make VERBOSE=1 -j4
make -j4
make install
                      
