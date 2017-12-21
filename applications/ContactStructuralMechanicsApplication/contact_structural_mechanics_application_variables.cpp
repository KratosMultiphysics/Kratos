// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix
//

// System includes

// External includes

// Project includes
#include "contact_structural_mechanics_application_variables.h"

namespace Kratos
{
typedef array_1d<double,3> Vector3;

// VARIABLES
/* Mortar method general variables */
KRATOS_CREATE_VARIABLE( Element::Pointer , ELEMENT_POINTER )                         // A pointer  the element belonging to this condition
KRATOS_CREATE_VARIABLE( int , INTEGRATION_ORDER_CONTACT )                            // The integration order computed in the contact
KRATOS_CREATE_VARIABLE( double, ACTIVE_CHECK_FACTOR )                                // The factor employed to search an active/inactive node

/* Weighted values */
KRATOS_CREATE_VARIABLE( double, WEIGHTED_GAP )                                       // The integrated gap employed in mortar formulation
KRATOS_CREATE_VARIABLE( double, WEIGHTED_SLIP )                                      // The integrated slip employed in mortar formulation
KRATOS_CREATE_VARIABLE( double, WEIGHTED_FRICTION )                                  // The integrated friction coefficient
KRATOS_CREATE_VARIABLE( double, WEIGHTED_SCALAR_RESIDUAL )                           // The integrated scalar residual  
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( WEIGHTED_VECTOR_RESIDUAL )                // The integrated vector residual    

/* For ALM mortar condition */
KRATOS_CREATE_VARIABLE( double, AUGMENTED_NORMAL_CONTACT_PRESSURE )                  // The resultant augmented pressure in the normal direction
KRATOS_CREATE_VARIABLE( double, AUGMENTED_TANGENT_CONTACT_PRESSURE )                 // The resultant augmented pressure in the tangent direction
KRATOS_CREATE_VARIABLE( double, TANGENT_FACTOR )                                     // The proportion between the tangent and normal penalty
KRATOS_CREATE_VARIABLE( int, CONSIDER_NORMAL_VARIATION )                             // A value used to check if consider normal variation or not
KRATOS_CREATE_VARIABLE( bool, CONSIDER_PAIR_VARIATION )                              // A value used to check if consider variation or not in the active inactive pairs
KRATOS_CREATE_VARIABLE( bool, ADAPT_PENALTY )                                        // To set if the penalty is recalculated or not
KRATOS_CREATE_VARIABLE( double, MAX_GAP_FACTOR )                                     // The factor between the nodal H and the max gap considered to recalculate the penalty

/* For mesh tying mortar condition */
KRATOS_CREATE_VARIABLE( std::string, TYING_VARIABLE )                                // The variable name for the mesh tying
}
