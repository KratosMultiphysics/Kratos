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

// Geometrical
KRATOS_CREATE_VARIABLE( double, AREA )
KRATOS_CREATE_VARIABLE( double, IX )
KRATOS_CREATE_VARIABLE( double, IY )
KRATOS_CREATE_VARIABLE( double, IZ )
KRATOS_CREATE_VARIABLE( double, CROSS_AREA )
KRATOS_CREATE_VARIABLE( double, MEAN_RADIUS )
KRATOS_CREATE_VARIABLE( int,    SECTION_SIDES )
KRATOS_CREATE_VARIABLE( Matrix , GEOMETRIC_STIFFNESS )

// Orthotropic/anisotropic constants
// KRATOS_CREATE_VARIABLE( double, YOUNG_MODULUS_X )
// KRATOS_CREATE_VARIABLE( double, YOUNG_MODULUS_Y )
// KRATOS_CREATE_VARIABLE( double, YOUNG_MODULUS_Z )
// KRATOS_CREATE_VARIABLE( double, SHEAR_MODULUS_XY )
// KRATOS_CREATE_VARIABLE( double, SHEAR_MODULUS_YZ )
// KRATOS_CREATE_VARIABLE( double, SHEAR_MODULUS_XZ )
// KRATOS_CREATE_VARIABLE( double, POISSON_RATIO_XY )
// KRATOS_CREATE_VARIABLE( double, POISSON_RATIO_YZ )
// KRATOS_CREATE_VARIABLE( double, POISSON_RATIO_XZ )

// Material orientation
// KRATOS_CREATE_VARIABLE( Vector3, MATERIAL_ORIENTATION_DX )
// KRATOS_CREATE_VARIABLE( Vector3, MATERIAL_ORIENTATION_DY )
// KRATOS_CREATE_VARIABLE( Vector3, MATERIAL_ORIENTATION_DZ )

// Shell generalized variables
KRATOS_CREATE_VARIABLE( Matrix, SHELL_STRAIN )
KRATOS_CREATE_VARIABLE( Matrix, SHELL_STRAIN_GLOBAL )
KRATOS_CREATE_VARIABLE( Matrix, SHELL_CURVATURE )
KRATOS_CREATE_VARIABLE( Matrix, SHELL_CURVATURE_GLOBAL )
KRATOS_CREATE_VARIABLE( Matrix, SHELL_FORCE )
KRATOS_CREATE_VARIABLE( Matrix, SHELL_FORCE_GLOBAL )
KRATOS_CREATE_VARIABLE( Matrix, SHELL_MOMENT )
KRATOS_CREATE_VARIABLE( Matrix, SHELL_MOMENT_GLOBAL )

// Cross section
KRATOS_CREATE_VARIABLE( ShellCrossSection::Pointer, SHELL_CROSS_SECTION )
KRATOS_CREATE_VARIABLE( int,          SHELL_CROSS_SECTION_OUTPUT_PLY_ID )
KRATOS_CREATE_VARIABLE( double,	SHELL_CROSS_SECTION_OUTPUT_PLY_LOCATION )

// Conditions
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( POINT_MOMENT )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( LOCAL_POINT_MOMENT )

// Adding the SPRISM EAS variables 
KRATOS_CREATE_VARIABLE(double, ALPHA_EAS);
KRATOS_CREATE_VARIABLE(bool, EAS_IMP);

// Adding the SPRISM additional variables 
KRATOS_CREATE_VARIABLE(double, ANG_ROT);

// Adding the Sprism number of transversal integration points 
KRATOS_CREATE_VARIABLE(int, NINT_TRANS);

// Additional strain measures 
KRATOS_CREATE_VARIABLE(Vector, HENCKY_STRAIN_VECTOR);
KRATOS_CREATE_VARIABLE(Matrix, HENCKY_STRAIN_TENSOR);

}
