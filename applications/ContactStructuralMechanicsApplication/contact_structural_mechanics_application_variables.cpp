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
#include "contact_structural_mechanics_application_variables.h"

namespace Kratos
{
typedef array_1d<double,3> Vector3;

// TODO: Clean all this!!!

// CONDITIONS
/* Mortar method general variables */
KRATOS_CREATE_VARIABLE( std::vector<contact_container>*, CONTACT_CONTAINERS )     // A vector of which contains the structure which defines the contact conditions // TODO: Remove this, deprecated
KRATOS_CREATE_VARIABLE( ConditionMap*, CONTACT_SETS )                              // An unordered map of which contains the structure which defines the contact conditions
KRATOS_CREATE_VARIABLE( Element::Pointer , ELEMENT_POINTER )                      // A pointer to the element belonging to this condition
KRATOS_CREATE_VARIABLE( int , INTEGRATION_ORDER_CONTACT )                         // The integration order computed in the contact
KRATOS_CREATE_VARIABLE( Matrix, MORTAR_CONTACT_OPERATOR )                         // Mortar Contact Operator
KRATOS_CREATE_VARIABLE( double, ACTIVE_CHECK_FACTOR )                             // The factor employed to search an active/inactive node

/* The complementary values */
// NOTE: This will be eventually not necessary
KRATOS_CREATE_VARIABLE( double, NORMAL_AUGMENTATION_FACTOR )                      // The constant that is considered for the check of active or inactive (when 0 it doesn't accept traction)
KRATOS_CREATE_VARIABLE( double, TANGENT_AUGMENTATION_FACTOR )                     // The constant that is considered for the check if the node is slip/stick

/* Weighted values */
KRATOS_CREATE_VARIABLE( double, WEIGHTED_GAP )                                    // The integrated gap employed in mortar formulation
KRATOS_CREATE_VARIABLE( double, WEIGHTED_SLIP )                                   // The integrated slip employed in mortar formulation

/* Matrix to store the derivatives of the normal */
KRATOS_CREATE_VARIABLE( Matrix, DELTA_NORMAL )                                    // Directional derivative of the normal

/* Auxiliar booleans to store the change in active/inactive slip/stick */
KRATOS_CREATE_VARIABLE( bool, AUXILIAR_ACTIVE )                                   // Auxiliar boolean to check if the node is active or not
KRATOS_CREATE_VARIABLE( bool, AUXILIAR_SLIP )                                     // Auxiliar boolean to check if the node is stick or not

/* For ALM mortar condition */
KRATOS_CREATE_VARIABLE( double, PENALTY_FACTOR )                                  // The penalty factor for the ALM
KRATOS_CREATE_VARIABLE( double, SCALE_FACTOR )                                    // The scale factor for the ALM

/* For mesh tying mortar condition */
KRATOS_CREATE_VARIABLE( std::string, TYING_VARIABLE )                             // The variable name for the mesh tying 
}
