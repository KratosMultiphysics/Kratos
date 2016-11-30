// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix
//
#include "contact_structural_mechanics_application_variables.h"

namespace Kratos
{
typedef array_1d<double,3> Vector3;

// CONDITIONS
/* Mortar method */
KRATOS_CREATE_VARIABLE( std::vector<contact_container>*, CONTACT_CONTAINERS )     // A vector of which contains the structure which defines the contact conditions
KRATOS_CREATE_VARIABLE( Element::Pointer , ELEMENT_POINTER )                      // A pointer to the element belonging to this condition
KRATOS_CREATE_VARIABLE( int , INTEGRATION_ORDER_CONTACT )                         // The integration order computed in the contact
KRATOS_CREATE_VARIABLE( Matrix, MORTAR_CONTACT_OPERATOR )                         // Mortar Contact Operator
KRATOS_CREATE_VARIABLE( double, ACTIVE_CHECK_FACTOR )                             // The factor employed to search an active/inactive node
KRATOS_CREATE_VARIABLE( double, NORMAL_AUGMENTATION_FACTOR )                      // The constant that is considered for the check of active or inactive (when 0 it doesn't accept traction)
KRATOS_CREATE_VARIABLE( double, TANGENT_AUGMENTATION_FACTOR )                     // The constant that is considered for the check if the node is slip/stick
KRATOS_CREATE_VARIABLE( double, WEIGHTED_GAP )                                    // The integrated gap employed in mortar formulation
KRATOS_CREATE_VARIABLE( double, WEIGHTED_SLIP )                                   // The integrated slip employed in mortar formulation
KRATOS_CREATE_VARIABLE( double, WEIGHTED_FRICTION )                               // The integrated friction employed in mortar formulation
KRATOS_CREATE_VARIABLE( Matrix, DELTA_NORMAL )                                    // Directional derivative of the normal
KRATOS_CREATE_VARIABLE( bool, AUXILIAR_ACTIVE )                                   // Auxiliar boolean to check if the node is active or not
KRATOS_CREATE_VARIABLE( bool, AUXILIAR_SLIP )                                     // Auxiliar boolean to check if the node is stick or not
KRATOS_CREATE_VARIABLE( double, GAP_GP )                                          // A double storing the gap of the GP
KRATOS_CREATE_VARIABLE( double, SLIP_GP )                                         // A double storing the slip of the GP
KRATOS_CREATE_VARIABLE( double, DOUBLE_LM_FACTOR )                                // The double LM parameter
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( DOUBLE_LM )                            // The double LM
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( NORMAL_CONTACT_STRESS_GP )             // For getting the normal contact stress in the GP
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( TANGENTIAL_CONTACT_STRESS_GP )         // For getting the tangential contact stress in the GP
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( NORMAL_GP )                            // For getting the normal in the GP
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( TANGENT_GP )                           // For getting the tangent in the GP
}
