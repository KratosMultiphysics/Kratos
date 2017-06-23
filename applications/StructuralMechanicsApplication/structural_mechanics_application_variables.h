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

#if !defined(KRATOS_STRUCTURAL_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED)
#define KRATOS_STRUCTURAL_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED

// System includes
#include<map>

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"
#include "solid_mechanics_application.h"
#include "structural_mechanics_application.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/shell_cross_section.hpp"
#include "custom_utilities/multipoint_constraint_data.hpp"

namespace Kratos
{
typedef array_1d<double,3> Vector3;
typedef std::string string;
typedef MpcData::Pointer MpcDataPointerType;


// Generalized eigenvalue problem
KRATOS_DEFINE_VARIABLE(int, BUILD_LEVEL)
KRATOS_DEFINE_VARIABLE(Vector, EIGENVALUE_VECTOR)
KRATOS_DEFINE_VARIABLE(Matrix, EIGENVECTOR_MATRIX)

// Geometrical
KRATOS_DEFINE_VARIABLE(double, AREA)
KRATOS_DEFINE_VARIABLE(double, IT)
KRATOS_DEFINE_VARIABLE(double, IY)
KRATOS_DEFINE_VARIABLE(double, IZ)
KRATOS_DEFINE_VARIABLE(double, CROSS_AREA)
KRATOS_DEFINE_VARIABLE(double, MEAN_RADIUS)
KRATOS_DEFINE_VARIABLE(int, SECTION_SIDES)
KRATOS_DEFINE_VARIABLE(Matrix, GEOMETRIC_STIFFNESS)


// Othotropic/anisotropic constants
// KRATOS_DEFINE_VARIABLE( double, YOUNG_MODULUS_X )
// KRATOS_DEFINE_VARIABLE( double, YOUNG_MODULUS_Y )
// KRATOS_DEFINE_VARIABLE( double, YOUNG_MODULUS_Z )
// KRATOS_DEFINE_VARIABLE( double, SHEAR_MODULUS_XY )
// KRATOS_DEFINE_VARIABLE( double, SHEAR_MODULUS_YZ )
// KRATOS_DEFINE_VARIABLE( double, SHEAR_MODULUS_XZ )
// KRATOS_DEFINE_VARIABLE( double, POISSON_RATIO_XY )
// KRATOS_DEFINE_VARIABLE( double, POISSON_RATIO_YZ )
// KRATOS_DEFINE_VARIABLE( double, POISSON_RATIO_XZ )

// Material orientation
// KRATOS_DEFINE_VARIABLE( Vector3, MATERIAL_ORIENTATION_DX )
// KRATOS_DEFINE_VARIABLE( Vector3, MATERIAL_ORIENTATION_DY )
// KRATOS_DEFINE_VARIABLE( Vector3, MATERIAL_ORIENTATION_DZ )

// Truss generalized variables
KRATOS_DEFINE_VARIABLE(double, TRUSS_PRESTRESS_PK2)
KRATOS_DEFINE_VARIABLE(bool, TRUSS_IS_CABLE)
// Beam generalized variables
KRATOS_DEFINE_VARIABLE(double, AREA_EFFECTIVE_Y)
KRATOS_DEFINE_VARIABLE(double, AREA_EFFECTIVE_Z)
KRATOS_DEFINE_VARIABLE(double, INERTIA_ROT_Y)
KRATOS_DEFINE_VARIABLE(double, INERTIA_ROT_Z)
KRATOS_DEFINE_VARIABLE(Vector, LOCAL_AXES_VECTOR)
KRATOS_DEFINE_VARIABLE(bool, LUMPED_MASS_MATRIX)

//ROCO NET
KRATOS_DEFINE_VARIABLE(double, KB_RING)
KRATOS_DEFINE_VARIABLE(double, DIAMETER_RING)
KRATOS_DEFINE_VARIABLE(double, WIRE_THICKNESS_RING)
KRATOS_DEFINE_VARIABLE(int, WIRE_WINDINGS_RING)
KRATOS_DEFINE_VARIABLE(double, NODAL_MASS_RING)
KRATOS_DEFINE_VARIABLE(double, KT_RING)



// Shell generalized variables
KRATOS_DEFINE_VARIABLE(Matrix, SHELL_STRAIN)
KRATOS_DEFINE_VARIABLE(Matrix, SHELL_STRAIN_GLOBAL)
KRATOS_DEFINE_VARIABLE(Matrix, SHELL_CURVATURE)
KRATOS_DEFINE_VARIABLE(Matrix, SHELL_CURVATURE_GLOBAL)
KRATOS_DEFINE_VARIABLE(Matrix, SHELL_FORCE)
KRATOS_DEFINE_VARIABLE(Matrix, SHELL_FORCE_GLOBAL)
KRATOS_DEFINE_VARIABLE(Matrix, SHELL_MOMENT)
KRATOS_DEFINE_VARIABLE(Matrix, SHELL_MOMENT_GLOBAL)

// Membrane1 vairiables
KRATOS_DEFINE_VARIABLE(double, PRESTRESS_11)
KRATOS_DEFINE_VARIABLE(double, PRESTRESS_22)
KRATOS_DEFINE_VARIABLE(double, PRESTRESS_12)


// Cross section
KRATOS_DEFINE_VARIABLE(ShellCrossSection::Pointer, SHELL_CROSS_SECTION)
KRATOS_DEFINE_VARIABLE(int, SHELL_CROSS_SECTION_OUTPUT_PLY_ID)
KRATOS_DEFINE_VARIABLE(double, SHELL_CROSS_SECTION_OUTPUT_PLY_LOCATION)

// Nodal stiffness for the nodal concentrated element
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( NODAL_STIFFNESS )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( NODAL_DAMPING_RATIO)

// CONDITIONS
/* Beam conditions */
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(POINT_MOMENT)
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(LOCAL_POINT_MOMENT)
/* Torque conditions */
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(POINT_TORQUE)

// Adding the SPRISM EAS variables
KRATOS_DEFINE_VARIABLE(double, ALPHA_EAS);
KRATOS_DEFINE_VARIABLE(bool, EAS_IMP);
KRATOS_DEFINE_VARIABLE(bool, SPRISM_TL_UL);

// Adding the SPRISM additional variables
KRATOS_DEFINE_VARIABLE(double, ANG_ROT); // TODO: Transform into a vector

// Adding the Sprism number of transversal integration points
KRATOS_DEFINE_VARIABLE(int, NINT_TRANS);

// Adding the SPRISM variable to deactivate the quadratic interpolation
KRATOS_DEFINE_VARIABLE(bool, QUAD_ON);

// Additional strain measures
KRATOS_DEFINE_VARIABLE(Vector, HENCKY_STRAIN_VECTOR);
KRATOS_DEFINE_VARIABLE(Matrix, HENCKY_STRAIN_TENSOR);

// For MPC implementations
KRATOS_DEFINE_VARIABLE(MpcDataPointerType, MPC_POINTER);
}

#endif /* KRATOS_STRUCTURAL_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED */
