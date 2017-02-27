### Using MPI in DamApplication ###

Note: For the moment, MPI only works in Linux and requires compiling METIS_APPLICATION and TRILINOS_APPLICATION. Non-local Damage does not work in MPI.

## Instructions to compile DamApplication for MPI (tested in Ubuntu 16.04) ##

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
-DDAM_APPLICATION=ON \

2. Uncomment (remove #~ ) the following line in GiDInterface/Kratos.gid/apps/Dam/python/dam_main.py

#~ import KratosMultiphysics.TrilinosApplication as TrilinosApplication

3. Uncomment (remove //~ ) the following lines in DamApplication/custom_python/add_custom_strategies_to_python.cpp

//~ #include "mpi.h"
//~ #include "Epetra_FECrsMatrix.h"
//~ #include "Epetra_FEVector.h"
//~ #include "trilinos_space.h"

//~ #include "custom_strategies/schemes/trilinos_incrementalupdate_static_damped_scheme.hpp"
//~ #include "custom_strategies/schemes/trilinos_dam_UP_scheme.hpp"

//~ typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
//~ typedef Scheme< TrilinosSparseSpaceType, LocalSpaceType > TrilinosBaseSchemeType;
//~ typedef TrilinosIncrementalUpdateStaticDampedScheme<TrilinosSparseSpaceType, LocalSpaceType> TrilinosIncrementalUpdateStaticDampedSchemeType;
//~ typedef TrilinosDamUPScheme<TrilinosSparseSpaceType, LocalSpaceType> TrilinosDamUPSchemeType;

//~ class_< TrilinosIncrementalUpdateStaticDampedSchemeType, bases<TrilinosBaseSchemeType>, boost::noncopyable >( "TrilinosIncrementalUpdateStaticDampedScheme", 
    //~ init< double >() );

//~ class_< TrilinosDamUPSchemeType, bases< TrilinosBaseSchemeType >,  boost::noncopyable >("TrilinosDamUPScheme",
    //~ init< double, double, double, double >());

4. Uncomment the following lines in DamApplication/CMakeLists.txt

#~ include_directories( ${CMAKE_SOURCE_DIR}/applications/trilinos_application )
#~ include_directories(${TRILINOS_INCLUDE_DIR})

#~ target_link_libraries(KratosDamApplication KratosCore KratosTrilinosApplication KratosConvectionDiffusionApplication KratosSolidMechanicsApplication KratosPoromechanicsApplication)

and comment the line just above the previous one:

target_link_libraries(KratosDamApplication KratosCore KratosConvectionDiffusionApplication KratosSolidMechanicsApplication KratosPoromechanicsApplication)
