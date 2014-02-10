//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
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

    // cross section
    AddCrossSectionsToPython();
    
    
    //registering variables in python ( if must to be seen from python )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONSTITUTIVE_LAW_NAME );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( WRITE_ID );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( PREVIOUS_DELTA_TIME );

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( INTERNAL_FORCE );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( EXTERNAL_FORCE );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONTACT_FORCE );

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( EXTERNAL_FORCES_VECTOR );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( INTERNAL_FORCES_VECTOR );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONTACT_FORCES_VECTOR );

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( CROSS_AREA );

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( IMPOSED_DISPLACEMENT );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( REACTION_PRESSURE );

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( VON_MISES_STRESS );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( PLASTIC_STRAIN );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( DELTA_PLASTIC_STRAIN );

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( DETERMINANT_F );
    
    // cross section
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHELL_CROSS_SECTION )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHELL_CROSS_SECTION_OUTPUT_PLY_ID )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHELL_CROSS_SECTION_OUTPUT_PLY_LOCATION )
    
    // shell generalized variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHELL_STRAIN )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHELL_FORCE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHELL_STRAIN_GLOBAL )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHELL_FORCE_GLOBAL )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHELL_CURVATURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHELL_MOMENT )
    
    // material orientation
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( MATERIAL_ORIENTATION_DX )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( MATERIAL_ORIENTATION_DY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( MATERIAL_ORIENTATION_DZ )
    
    // orthotropic/anisotropic constants
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( YOUNG_MODULUS_X )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( YOUNG_MODULUS_Y )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( YOUNG_MODULUS_Z )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHEAR_MODULUS_XY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHEAR_MODULUS_YZ )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHEAR_MODULUS_XZ )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( POISSON_RATIO_XY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( POISSON_RATIO_YZ )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( POISSON_RATIO_XZ )
    
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( PLASTIC_DISSIPATION );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( DELTA_PLASTIC_DISSIPATION );

}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
