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
    Condition*            condition;
    double             contact_area;
    std::vector<double> contact_gap;
    std::vector<bool>  active_nodes;
  
    void print()
    {
       KRATOS_WATCH(condition);
       KRATOS_WATCH(contact_area);
       for (unsigned int i = 0; i < contact_gap.size(); i++)
       {
            std::cout << "contact_gap_" << i << " " << contact_gap[i] << std::endl;
       }
       
       for (unsigned int i = 0; i < active_nodes.size(); i++)
       {
            std::cout << "active_nodes_" << i << " " << active_nodes[i] << std::endl;
       }
    }
    
    void save( Serializer& rSerializer ) const
    {
        rSerializer.save("condition",       condition);
        rSerializer.save("contact_area", contact_area);
        rSerializer.save("contact_gap",   contact_gap);
//         rSerializer.save("active_nodes", active_nodes); // NOTE: Error, ask
    }

    void load( Serializer& rSerializer )
    {
        rSerializer.load("condition",       condition);
        rSerializer.load("contact_area", contact_area);
        rSerializer.load("contact_gap",   contact_gap);
//         rSerializer.load("active_nodes", active_nodes); // NOTE: Error, ask 
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
/* Mortar method */ // TODO: Remove unnecessary variables
KRATOS_DEFINE_VARIABLE( std::vector<contact_container>*, CONTACT_CONTAINERS ) // A vector of which contains the structure which defines the contact conditions
KRATOS_DEFINE_VARIABLE( Element::Pointer, ELEMENT_POINTER )                   // A pointer to the element belonging to this condition
KRATOS_DEFINE_VARIABLE( int , INTEGRATION_ORDER_CONTACT )                     // The integration order computed in the contact
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( CONTACT_MESH_TYING_FORCE )         // The "force" resulting from contact
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( LAGRANGE_MULTIPLIER )              // The lagrange multiplier used to enforce the mortar constraint
KRATOS_DEFINE_VARIABLE( Matrix, MORTAR_CONTACT_OPERATOR )                     // Mortar Contact Operator
KRATOS_DEFINE_VARIABLE( Matrix, CONTACT_STIFFNESS_MATRIX )                    // LHS contribution
KRATOS_DEFINE_VARIABLE( Matrix, GAP_DERIVATIVES_MATRIX )                      // LHS contribution
KRATOS_DEFINE_VARIABLE( Vector, CONTACT_FORCES_VECTOR )                       // RHS contribution
KRATOS_DEFINE_VARIABLE( Vector, NORMAL_GAPS_VECTOR )                          // RHS contribution

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
