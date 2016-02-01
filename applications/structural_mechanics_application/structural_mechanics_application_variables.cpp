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
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
typedef array_1d<double,3> Vector3;

//geometrical
KRATOS_CREATE_VARIABLE( double, AREA )
KRATOS_CREATE_VARIABLE( double, IX )
KRATOS_CREATE_VARIABLE( double, IY )
KRATOS_CREATE_VARIABLE( double, IZ )
KRATOS_CREATE_VARIABLE( double, CROSS_AREA )
KRATOS_CREATE_VARIABLE( double, MEAN_RADIUS )
KRATOS_CREATE_VARIABLE( int,    SECTION_SIDES )
KRATOS_CREATE_VARIABLE( Matrix , GEOMETRIC_STIFFNESS )

//othotropic/anisotropic constants
// KRATOS_CREATE_VARIABLE( double, YOUNG_MODULUS_X )
// KRATOS_CREATE_VARIABLE( double, YOUNG_MODULUS_Y )
// KRATOS_CREATE_VARIABLE( double, YOUNG_MODULUS_Z )
// KRATOS_CREATE_VARIABLE( double, SHEAR_MODULUS_XY )
// KRATOS_CREATE_VARIABLE( double, SHEAR_MODULUS_YZ )
// KRATOS_CREATE_VARIABLE( double, SHEAR_MODULUS_XZ )
// KRATOS_CREATE_VARIABLE( double, POISSON_RATIO_XY )
// KRATOS_CREATE_VARIABLE( double, POISSON_RATIO_YZ )
// KRATOS_CREATE_VARIABLE( double, POISSON_RATIO_XZ )

//material orientation
// KRATOS_CREATE_VARIABLE( Vector3, MATERIAL_ORIENTATION_DX )
// KRATOS_CREATE_VARIABLE( Vector3, MATERIAL_ORIENTATION_DY )
// KRATOS_CREATE_VARIABLE( Vector3, MATERIAL_ORIENTATION_DZ )

//shell generalized variables
KRATOS_CREATE_VARIABLE( Matrix, SHELL_STRAIN )
KRATOS_CREATE_VARIABLE( Matrix, SHELL_STRAIN_GLOBAL )
KRATOS_CREATE_VARIABLE( Matrix, SHELL_CURVATURE )
KRATOS_CREATE_VARIABLE( Matrix, SHELL_CURVATURE_GLOBAL )
KRATOS_CREATE_VARIABLE( Matrix, SHELL_FORCE )
KRATOS_CREATE_VARIABLE( Matrix, SHELL_FORCE_GLOBAL )
KRATOS_CREATE_VARIABLE( Matrix, SHELL_MOMENT )
KRATOS_CREATE_VARIABLE( Matrix, SHELL_MOMENT_GLOBAL )

//cross section
KRATOS_CREATE_VARIABLE( ShellCrossSection::Pointer, SHELL_CROSS_SECTION )
KRATOS_CREATE_VARIABLE( int,          SHELL_CROSS_SECTION_OUTPUT_PLY_ID )
KRATOS_CREATE_VARIABLE( double,	SHELL_CROSS_SECTION_OUTPUT_PLY_LOCATION )

}
