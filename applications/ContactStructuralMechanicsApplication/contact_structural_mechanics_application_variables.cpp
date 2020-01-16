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
typedef Geometry<Node<3>> GeometryType;

// VARIABLES
// MPC Contact related variables
KRATOS_CREATE_VARIABLE( MasterSlaveConstraint::Pointer, CONSTRAINT_POINTER )      // Pointer to the constraint of the condition
KRATOS_CREATE_VARIABLE( double, REACTION_CHECK_STIFFNESS_FACTOR )                 // The reaction factor to be considered on the tension check

/* Mortar method general variables */
KRATOS_CREATE_VARIABLE( int , INNER_LOOP_ITERATION )                              // The number of loops in the simplified semi-smooth inner iteration
KRATOS_CREATE_VARIABLE( int , INTEGRATION_ORDER_CONTACT )                         // The integration order computed in the contact
KRATOS_CREATE_VARIABLE( double, DISTANCE_THRESHOLD )                              // The distance threshold considered
KRATOS_CREATE_VARIABLE( double, ZERO_TOLERANCE_FACTOR )                           // The epsilon factor considered
KRATOS_CREATE_VARIABLE( double, ACTIVE_CHECK_FACTOR )                             // The factor employed to search an active/inactive node
KRATOS_CREATE_VARIABLE( double, SLIP_THRESHOLD )                                  // The threshold employed to search an slip/stick node
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( AUXILIAR_COORDINATES )                 // Auxiliar coordinates used to map
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( DELTA_COORDINATES )                    // Delta coordinates used to map
KRATOS_CREATE_VARIABLE( double, NORMAL_GAP )                                      // The normal gap employed in contact formulation
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( TANGENT_SLIP )                         // The tangent slip gap employed in contact formulation

/* Weighted values */
KRATOS_CREATE_VARIABLE( double, WEIGHTED_GAP )                                    // The integrated gap employed in mortar formulation
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( WEIGHTED_SLIP )                        // The integrated slip employed in mortar formulation
KRATOS_CREATE_VARIABLE( double, WEIGHTED_SCALAR_RESIDUAL )                        // The integrated scalar residual
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( WEIGHTED_VECTOR_RESIDUAL )             // The integrated vector residual

/* For ALM mortar condition */
KRATOS_CREATE_VARIABLE( bool, ACTIVE_SET_COMPUTED )                               // To know if the active set has been computed
KRATOS_CREATE_VARIABLE( bool, ACTIVE_SET_CONVERGED )                              // To know if the active set has converged
KRATOS_CREATE_VARIABLE( bool, SLIP_SET_CONVERGED )                                // To know if the slip set has converged
KRATOS_CREATE_VARIABLE( double, OPERATOR_THRESHOLD )                              // Consider objetive/non-objetive formulation threshold
KRATOS_CREATE_VARIABLE( double, SLIP_AUGMENTATION_COEFFICIENT )                   // Coefficient to improve the slip computation convergence (augmented part related)
KRATOS_CREATE_VARIABLE( double, DYNAMIC_FACTOR )                                  // The factor considered for dynamic problems (in order to take intro account the gap evolution)
KRATOS_CREATE_VARIABLE( double, LAGRANGE_MULTIPLIER_CONTACT_PRESSURE )            // The lagrange multiplier for normal contact pressure
KRATOS_CREATE_VARIABLE( double, AUGMENTED_NORMAL_CONTACT_PRESSURE )               // The resultant augmented pressure in the normal direction
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( AUGMENTED_TANGENT_CONTACT_PRESSURE )   // The resultant augmented pressure in the tangent direction
KRATOS_CREATE_VARIABLE( int, CONSIDER_NORMAL_VARIATION )                          // A value used to check if consider normal variation or not
KRATOS_CREATE_VARIABLE( bool, ADAPT_PENALTY )                                     // To set if the penalty is recalculated or not
KRATOS_CREATE_VARIABLE( double, MAX_GAP_FACTOR )                                  // The factor between the nodal H and the max gap considered to recalculate the penalty

/* For mesh tying mortar condition */
KRATOS_CREATE_VARIABLE( std::string, TYING_VARIABLE )                             // The variable name for the mesh tying

/* Explicit simulation */
KRATOS_CREATE_VARIABLE( double, MAX_GAP_THRESHOLD )                               // The gap considered as threshold to rescale penalty

/* Frictional laws */
KRATOS_CREATE_VARIABLE( FrictionalLaw::Pointer, FRICTIONAL_LAW )                  // The frictional law considered
KRATOS_CREATE_VARIABLE( double, TRESCA_FRICTION_THRESHOLD )                       // The threshold value for Tresca frictional contact
}
