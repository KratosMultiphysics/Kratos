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
rm -rf CMakeFiles/

export BLA_VENDOR="ATLAS"

# Release with debug information
#-DCMAKE_BUILD_TYPE=RelWithDebInfo 
#-DCMAKE_BUILD_TYPE=Release

cmake ..  									\
-DCMAKE_C_COMPILER=gcc								\
-DCMAKE_CXX_COMPILER=g++							\
-DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -msse3 " 					\
-DCMAKE_C_FLAGS="${CMAKE_C_FLAGS} -msse3 " 					\
-DCMAKE_BUILD_TYPE=RelWithDebInfo     					        \
-DINCOMPRESSIBLE_FLUID_APPLICATION=OFF  					\
-DMESHING_APPLICATION=OFF 							\
-DEXTERNAL_SOLVERS_APPLICATION=ON						\
-DPFEM_APPLICATION=OFF 								\
-DSTRUCTURAL_APPLICATION=ON 							\
-DSOLID_MECHANICS_APPLICATION=ON   					        \
-DPFEM_SOLID_MECHANICS_APPLICATION=ON  					        \
-DPFEM_SOLID_APPLICATION=OFF 							\
-DCONVECTION_DIFFUSION_APPLICATION=OFF 						\
-DFLUID_DYNAMICS_APPLICATION=OFF 						\
-DALE_APPLICATION=OFF								\
-DFSI_APPLICATION=OFF								\
-DOPENCL_APPLICATION=OFF							\
-DMIXED_ELEMENT_APPLICATION=OFF							\
-DMKL_SOLVERS_APPLICATION=OFF							\
-DMKLSOLVER_INCLUDE_DIR="/opt/intel/mkl/include"              		        \
-DMKLSOLVER_LIB_DIR="/opt/intel/mkl/lib/intel64"		                \
-DMETIS_APPLICATION=OFF								\
-DTRILINOS_APPLICATION=OFF							\
-DDEM_APPLICATION=OFF								\
-DINCLUDE_PASTIX=ON                                                             \
-DKRATOS_SUPERLU_FORTRAN_MANGLING="-DAdd_"                                      \
-DPASTIX_INSTALL_DIR="/home/jmaria/kratos/external_libraries/pastix_release_4030/install" \
-DSCOTCH_INSTALL_DIR="/home/jmaria/kratos/external_libraries/scotch_6.0.0/lib"            \
-DEXCLUDE_ITSOL=ON                                                                        \
-DINSTALL_EMBEDDED_PYTHON=ON                                                              \

#-DPARMETIS_ROOT_DIR="/home/rrossi/compiled_libraries/ParMetis-3.1.1" 		\
#-DTRILINOS_ROOT="/home/rrossi/compiled_libraries/trilinos-10.2.0"		\
#-DPYTHON_INCLUDE_DIRS="/usr/include/python2.7"                                  \
#-DPYTHON_LIBRARIES="/usr/lib/python2.7/config/"                                 \
#-DPYTHON_LIBRARY:FILEPATH="/usr/lib/python2.7/config/libpython2.7.so"           \



#decomment this to have it verbose
# make VERBOSE=1 -j4
make -j4
make install

