### Using MPI in PoromechanicsApplication ###

Note: For the moment, MPI only works in Linux and requires compiling METIS_APPLICATION and TRILINOS_APPLICATION. Non-local Damage and Fracture Propagation features do not work in MPI.

## Instructions to compile PoromechanicsApplication for MPI (tested in Ubuntu 16.04) ##

1. Make sure that the following lines are properly set in the configure.sh file:

-DMETIS_APPLICATION=ON								                                        \
-DMETIS_INCLUDE_DIR="/usr/include/"                                                         \
-DUSE_METIS_5=ON                                                                            \
-DPARMETIS_ROOT_DIR="/usr/lib/"                                                   			\
-DTRILINOS_APPLICATION=ON							                                        \
-DTRILINOS_LIBRARY_DIR="/usr/lib/x86_64-linux-gnu/"                                         \
-DTRILINOS_INCLUDE_DIR="/usr/include/trilinos/"                                             \
-DTRILINOS_LIBRARY_PREFIX="trilinos_"                                                       \
-DEXTERNAL_SOLVERS_APPLICATION=ON						                                    \
-DSOLID_MECHANICS_APPLICATION=ON   					                                    \
-DFLUID_DYNAMICS_APPLICATION=ON   					                                    \
-DPOROMECHANICS_APPLICATION=ON \
-DUSE_PORO_MPI=ON \

2. Uncomment (remove #~ ) the following line in PoromechanicsApplication/custom_problemtype/Poromechanics_Application.gid/poromechanics_main.py

#~ import KratosMultiphysics.TrilinosApplication as TrilinosApplication
