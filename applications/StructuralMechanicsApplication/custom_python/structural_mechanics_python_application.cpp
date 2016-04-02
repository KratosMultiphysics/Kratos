// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "structural_mechanics_application_variables.h"
#include "structural_mechanics_application.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_cross_sections_to_python.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;



BOOST_PYTHON_MODULE(KratosStructuralMechanicsApplication)
{

    class_<KratosStructuralMechanicsApplication,
           KratosStructuralMechanicsApplication::Pointer,
           bases<KratosApplication>, boost::noncopyable >("KratosStructuralMechanicsApplication")
           ;

    AddCustomStrategiesToPython();
    AddCustomUtilitiesToPython();
    AddCrossSectionsToPython();

    //registering variables in python
    // cross section
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHELL_CROSS_SECTION )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHELL_CROSS_SECTION_OUTPUT_PLY_ID )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHELL_CROSS_SECTION_OUTPUT_PLY_LOCATION )

    // shell generalized variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHELL_STRAIN )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHELL_STRAIN_GLOBAL )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHELL_CURVATURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHELL_CURVATURE_GLOBAL )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHELL_FORCE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHELL_FORCE_GLOBAL )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHELL_MOMENT )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHELL_MOMENT_GLOBAL )

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( POINT_MOMENT )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( LOCAL_POINT_MOMENT )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( CROSS_AREA )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( MEAN_RADIUS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SECTION_SIDES )

    // material orientation
//     KRATOS_REGISTER_IN_PYTHON_VARIABLE( MATERIAL_ORIENTATION_DX )
//     KRATOS_REGISTER_IN_PYTHON_VARIABLE( MATERIAL_ORIENTATION_DY )
//     KRATOS_REGISTER_IN_PYTHON_VARIABLE( MATERIAL_ORIENTATION_DZ )

    // orthotropic/anisotropic constants
//     KRATOS_REGISTER_IN_PYTHON_VARIABLE( YOUNG_MODULUS_X )
//     KRATOS_REGISTER_IN_PYTHON_VARIABLE( YOUNG_MODULUS_Y )
//     KRATOS_REGISTER_IN_PYTHON_VARIABLE( YOUNG_MODULUS_Z )
//     KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHEAR_MODULUS_XY )
//     KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHEAR_MODULUS_YZ )
//     KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHEAR_MODULUS_XZ )
//     KRATOS_REGISTER_IN_PYTHON_VARIABLE( POISSON_RATIO_XY )
//     KRATOS_REGISTER_IN_PYTHON_VARIABLE( POISSON_RATIO_YZ )
//     KRATOS_REGISTER_IN_PYTHON_VARIABLE( POISSON_RATIO_XZ )

    /* Adding the SPRISM EAS variables */
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(ALPHA_EAS);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(EAS_IMP);

    /* Adding the SPRISM additional variables */
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(ANG_ROT);

    /* Adding the SPRISM number of transversal integration points */
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(NINT_TRANS);

    /* Adding the SPRISM variable to deactivate the quadratic interpolation */
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(QUAD_ON);

    /* Hencky strain */
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(HENCKY_STRAIN_VECTOR);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(HENCKY_STRAIN_TENSOR);
}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
