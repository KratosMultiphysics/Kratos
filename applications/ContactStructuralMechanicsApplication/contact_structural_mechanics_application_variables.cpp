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
KRATOS_CREATE_VARIABLE( boost::shared_ptr<ConditionSet>, CONTACT_SETS )              // An unordered map of which contains the structure which defines the contact conditions
KRATOS_CREATE_VARIABLE( Element::Pointer , ELEMENT_POINTER )                         // A pointer  the element belonging to this condition
KRATOS_CREATE_VARIABLE( int , INTEGRATION_ORDER_CONTACT )                            // The integration order computed in the contact
KRATOS_CREATE_VARIABLE( Matrix, MORTAR_CONTACT_OPERATOR )                            // Mortar Contact Operator
KRATOS_CREATE_VARIABLE( double, ACTIVE_CHECK_FACTOR )                                // The factor employed to search an active/inactive node

/* Weighted values */
KRATOS_CREATE_VARIABLE( double, WEIGHTED_GAP )                                       // The integrated gap employed in mortar formulation
KRATOS_CREATE_VARIABLE( double, WEIGHTED_SLIP )                                      // The integrated slip employed in mortar formulation
KRATOS_CREATE_VARIABLE( double, WEIGHTED_FRICTION )                                  // The integrated friction coefficient
KRATOS_CREATE_VARIABLE( double, WEIGHTED_SCALAR_RESIDUAL )                           // The integrated scalar residual  
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( WEIGHTED_VECTOR_RESIDUAL )                // The integrated vector residual    

/* Matrix to store the derivatives of the normal */
KRATOS_CREATE_VARIABLE( Matrix, DELTA_NORMAL )                                       // Directional derivative of the normal

/* For ALM mortar condition */
KRATOS_CREATE_VARIABLE( double, AUGMENTED_NORMAL_CONTACT_PRESSURE )                  // The resultant augmented pressure in the normal direction
KRATOS_CREATE_VARIABLE( double, AUGMENTED_TANGENT_CONTACT_PRESSURE )                 // The resultant augmented pressure in the tangent direction
KRATOS_CREATE_VARIABLE( double, PENALTY_PARAMETER )                                  // The penalty factor for the ALM
KRATOS_CREATE_VARIABLE( double, SCALE_FACTOR )                                       // The scale factor for the ALM
KRATOS_CREATE_VARIABLE( double, TANGENT_FACTOR )                                     // The proportion between the tangent and normal penalty
KRATOS_CREATE_VARIABLE( bool, CONSIDER_NORMAL_VARIATION )                            // A value used to check if consider normal variation or not
KRATOS_CREATE_VARIABLE( bool, CONSIDER_PAIR_VARIATION )                              // A value used to check if consider variation or not in the active inactive pairs

/* For mesh tying mortar condition */
KRATOS_CREATE_VARIABLE( std::string, TYING_VARIABLE )                                // The variable name for the mesh tying

/* Processes utilities */
// KRATOS_CREATE_VARIABLE( boost::shared_ptr<ProcessFactoryUtility>, PROCESSES_LIST )   // A pointer to the processes list
}
