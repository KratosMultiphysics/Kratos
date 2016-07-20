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

#if !defined(KRATOS_STRUCTURAL_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_STRUCTURAL_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED

// System includes

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

namespace Kratos
{
typedef array_1d<double,3> Vector3;

struct contact_container 
{
    Condition::Pointer                 condition;
    std::vector<double>  local_coordinates_slave;
  
    ~contact_container(){}
    
    void clear()
    {
        condition = nullptr;
        local_coordinates_slave.clear();
    }
    
    void print()
    {
       KRATOS_WATCH(condition);
       std::cout << std::endl;
       for (unsigned int i = 0; i < local_coordinates_slave.size(); i++)
       {
            std::cout << "local_coordinates_slave_" << i << " " << local_coordinates_slave[i] << std::endl;
       }
       std::cout << std::endl;
       
    }
    
    void save( Serializer& rSerializer ) const
    {
        rSerializer.save("condition",                               condition);
        rSerializer.save("local_coordinates_slave",   local_coordinates_slave);
    }

    void load( Serializer& rSerializer )
    {
        rSerializer.load("condition",                              condition);
        rSerializer.load("local_coordinates_slave",  local_coordinates_slave);
    }
};

// Geometrical
KRATOS_DEFINE_VARIABLE( double, AREA )
KRATOS_DEFINE_VARIABLE( double, IX )
KRATOS_DEFINE_VARIABLE( double, IY )
KRATOS_DEFINE_VARIABLE( double, IZ )
KRATOS_DEFINE_VARIABLE( double, CROSS_AREA )
KRATOS_DEFINE_VARIABLE( double, MEAN_RADIUS )
KRATOS_DEFINE_VARIABLE( int,    SECTION_SIDES )
KRATOS_DEFINE_VARIABLE( Matrix , GEOMETRIC_STIFFNESS )

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

// Shell generalized variables
KRATOS_DEFINE_VARIABLE( Matrix, SHELL_STRAIN )
KRATOS_DEFINE_VARIABLE( Matrix, SHELL_STRAIN_GLOBAL )
KRATOS_DEFINE_VARIABLE( Matrix, SHELL_CURVATURE )
KRATOS_DEFINE_VARIABLE( Matrix, SHELL_CURVATURE_GLOBAL )
KRATOS_DEFINE_VARIABLE( Matrix, SHELL_FORCE )
KRATOS_DEFINE_VARIABLE( Matrix, SHELL_FORCE_GLOBAL )
KRATOS_DEFINE_VARIABLE( Matrix, SHELL_MOMENT )
KRATOS_DEFINE_VARIABLE( Matrix, SHELL_MOMENT_GLOBAL )

// Cross section
KRATOS_DEFINE_VARIABLE( ShellCrossSection::Pointer, SHELL_CROSS_SECTION )
KRATOS_DEFINE_VARIABLE( int,          SHELL_CROSS_SECTION_OUTPUT_PLY_ID )
KRATOS_DEFINE_VARIABLE( double,	SHELL_CROSS_SECTION_OUTPUT_PLY_LOCATION )

// Nodal stiffness for the nodal concentrated element
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( NODAL_STIFFNESS )

// CONDITIONS
/* Beam conditions */
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( POINT_MOMENT )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( LOCAL_POINT_MOMENT )
/* Torque conditions */
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( POINT_TORQUE )
/* Mortar method */ 
KRATOS_DEFINE_VARIABLE( std::vector<contact_container>*, CONTACT_CONTAINERS ) // A vector of which contains the structure which defines the contact conditions
KRATOS_DEFINE_VARIABLE( Element::Pointer, ELEMENT_POINTER )                   // A pointer to the element belonging to this condition
KRATOS_DEFINE_VARIABLE( int , INTEGRATION_ORDER_CONTACT )                     // The integration order computed in the contact
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( LAGRANGE_MULTIPLIER )              // The lagrange multiplier used to enforce the mortar constraint
KRATOS_DEFINE_VARIABLE( Matrix, MORTAR_CONTACT_OPERATOR )                     // Mortar Contact Operator
KRATOS_DEFINE_VARIABLE( double, ACTIVE_CHECK_FACTOR )                         // The factor employed to consider an active/inactive node
KRATOS_DEFINE_VARIABLE( double, WEIGHTED_GAP )                                // The integrated gap employed in mortar formulation

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
}

#endif	/* KRATOS_STRUCTURAL_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED */
