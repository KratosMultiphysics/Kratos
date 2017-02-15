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
-DPOROMECHANICS_APPLICATION=ON \

2. Uncomment (remove #~ ) the following line in PoromechanicsApplication/custom_problemtype/Poromechanics_Application.gid/poromechanics_main.py

#~ import KratosMultiphysics.TrilinosApplication as TrilinosApplication

3. Uncomment (remove //~ ) the following lines in PoromechanicsApplication/custom_python/add_custom_strategies_to_python.cpp

//~ #include "mpi.h"
//~ #include "Epetra_FECrsMatrix.h"
//~ #include "Epetra_FEVector.h"
//~ #include "trilinos_space.h"

//~ #include "custom_strategies/schemes/trilinos_newmark_quasistatic_U_Pw_scheme.hpp"
//~ #include "custom_strategies/schemes/trilinos_newmark_quasistatic_damped_U_Pw_scheme.hpp"
//~ #include "custom_strategies/schemes/trilinos_newmark_dynamic_U_Pw_scheme.hpp"

//~ typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
//~ typedef Scheme< TrilinosSparseSpaceType, LocalSpaceType > TrilinosBaseSchemeType;
//~ typedef TrilinosNewmarkQuasistaticUPwScheme<TrilinosSparseSpaceType, LocalSpaceType> TrilinosNewmarkQuasistaticUPwSchemeType;
//~ typedef TrilinosNewmarkQuasistaticDampedUPwScheme<TrilinosSparseSpaceType, LocalSpaceType> TrilinosNewmarkQuasistaticDampedUPwSchemeType;
//~ typedef TrilinosNewmarkDynamicUPwScheme<TrilinosSparseSpaceType, LocalSpaceType> TrilinosNewmarkDynamicUPwSchemeType;

//~ class_< TrilinosNewmarkQuasistaticUPwSchemeType, bases<TrilinosBaseSchemeType>, boost::noncopyable >( "TrilinosNewmarkQuasistaticUPwScheme", 
    //~ init< double, double, double >() );

//~ class_< TrilinosNewmarkQuasistaticDampedUPwSchemeType,bases< TrilinosBaseSchemeType >, boost::noncopyable >("TrilinosNewmarkQuasistaticDampedUPwScheme",
    //~ init<  double, double, double, double, double >());

//~ class_< TrilinosNewmarkDynamicUPwSchemeType,bases< TrilinosBaseSchemeType >, boost::noncopyable >("TrilinosNewmarkDynamicUPwScheme",
    //~ init<  double, double, double, double, double >());

4. Uncomment the following lines in PoromechanicsApplication/CMakeLists.txt

#~ include_directories( ${CMAKE_SOURCE_DIR}/applications/trilinos_application )
#~ include_directories(${TRILINOS_INCLUDE_DIR})

#~ target_link_libraries(KratosPoromechanicsApplication KratosCore KratosTrilinosApplication KratosSolidMechanicsApplication)

and comment the line just above the previous one:

target_link_libraries(KratosPoromechanicsApplication KratosCore KratosSolidMechanicsApplication)
