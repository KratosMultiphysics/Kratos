//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes
#if defined(KRATOS_PYTHON)

// External includes
#include <boost/python.hpp>

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"
#include "custom_python/add_cross_sections_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"

#include "solid_mechanics_application.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;



BOOST_PYTHON_MODULE(KratosSolidMechanicsApplication)
{

    class_<KratosSolidMechanicsApplication,
           KratosSolidMechanicsApplication::Pointer,
           bases<KratosApplication>, boost::noncopyable >("KratosSolidMechanicsApplication")
           ;

    AddCustomUtilitiesToPython();
    AddCustomStrategiesToPython();
    AddCustomConstitutiveLawsToPython(); 
    AddCrossSectionsToPython();
    AddCustomProcessesToPython();

    //registering variables in python ( if must to be seen from python )

    // Generalized eigenvalue problem
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( BUILD_LEVEL )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( EIGENVALUE_VECTOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( EIGENVECTOR_MATRIX )

    //For process information
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( LOAD_RESTART )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( COMPUTE_LUMPED_MASS_MATRIX )
   
    //For explicit schemes
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MIDDLE_VELOCITY )
            
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( POINT_LOAD )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( POINT_MOMENT )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( LINE_LOAD )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( SURFACE_LOAD )

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( FOLLOWER_POINT_MOMENT )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( FOLLOWER_POINT_LOAD )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( FOLLOWER_LINE_LOAD )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( FOLLOWER_SURFACE_LOAD )
      
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( LINE_LOADS_VECTOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SURFACE_LOADS_VECTOR )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( WRITE_ID )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( RAYLEIGH_ALPHA )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( RAYLEIGH_BETA )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( EXTERNAL_FORCES_VECTOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( INTERNAL_FORCES_VECTOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONTACT_FORCES_VECTOR )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( STABILIZATION_FACTOR )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( PRESSURE_REACTION )
      
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( VON_MISES_STRESS )

    //For beams
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( CROSS_AREA )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( MEAN_RADIUS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SECTION_SIDES )
      
    //For shells cross section
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHELL_CROSS_SECTION )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHELL_CROSS_SECTION_OUTPUT_PLY_ID )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHELL_CROSS_SECTION_OUTPUT_PLY_LOCATION )
    
    //For shell generalized variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHELL_STRAIN )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHELL_STRAIN_GLOBAL )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHELL_CURVATURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHELL_CURVATURE_GLOBAL )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHELL_FORCE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHELL_FORCE_GLOBAL )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHELL_MOMENT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHELL_MOMENT_GLOBAL )


}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
