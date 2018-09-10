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
typedef Geometry<Node<3>> GeometryType;

// VARIABLES
/* Mortar method general variables */
KRATOS_CREATE_VARIABLE( int , INTEGRATION_ORDER_CONTACT )                            // The integration order computed in the contact
KRATOS_CREATE_VARIABLE( double, DISTANCE_THRESHOLD )                                 // The distance threshold considered
KRATOS_CREATE_VARIABLE( double, ACTIVE_CHECK_FACTOR )                                // The factor employed to search an active/inactive node
KRATOS_CREATE_VARIABLE( GeometryType::Pointer, PAIRED_GEOMETRY )                     // The paired geometry with the current condition
KRATOS_CREATE_VARIABLE( Vector3, PAIRED_NORMAL )                                     // The normal of the paired geometry
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( AUXILIAR_COORDINATES )                    // Auxiliar coordinates used to map
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( DELTA_COORDINATES )                       // Delta coordinates used to map
KRATOS_CREATE_VARIABLE( double, NORMAL_GAP )                                         // The normal gap employed in contact formulation

/* Weighted values */
KRATOS_CREATE_VARIABLE( double, WEIGHTED_GAP )                                       // The integrated gap employed in mortar formulation
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( WEIGHTED_SLIP )                           // The integrated slip employed in mortar formulation
KRATOS_CREATE_VARIABLE( double, WEIGHTED_SCALAR_RESIDUAL )                           // The integrated scalar residual
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( WEIGHTED_VECTOR_RESIDUAL )                // The integrated vector residual

/* For ALM mortar condition */
KRATOS_CREATE_VARIABLE( bool, ACTIVE_SET_CONVERGED )                                 // To know if the active set has converged
KRATOS_CREATE_VARIABLE( double, DYNAMIC_FACTOR )                                     // The factor considered for dynamic problems (in order to take intro account the gap evolution)
KRATOS_CREATE_VARIABLE( double, LAGRANGE_MULTIPLIER_CONTACT_PRESSURE )               // The lagrange multiplier for normal contact pressure
KRATOS_CREATE_VARIABLE( double, AUGMENTED_NORMAL_CONTACT_PRESSURE )                  // The resultant augmented pressure in the normal direction
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( AUGMENTED_TANGENT_CONTACT_PRESSURE )      // The resultant augmented pressure in the tangent direction
KRATOS_CREATE_VARIABLE( int, CONSIDER_NORMAL_VARIATION )                             // A value used to check if consider normal variation or not
KRATOS_CREATE_VARIABLE( bool, ADAPT_PENALTY )                                        // To set if the penalty is recalculated or not
KRATOS_CREATE_VARIABLE( double, MAX_GAP_FACTOR )                                     // The factor between the nodal H and the max gap considered to recalculate the penalty

/* For mesh tying mortar condition */
KRATOS_CREATE_VARIABLE( std::string, TYING_VARIABLE )                                // The variable name for the mesh tying
}
