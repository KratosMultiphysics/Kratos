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

#if !defined(KRATOS_CONTACT_STRUCTURAL_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_CONTACT_STRUCTURAL_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"
#include "includes/mapping_variables.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    typedef array_1d<double,3> Vector3;
    typedef Geometry<Node<3>> GeometryType;

///@}
///@name  Enum's
///@{

    enum NormalDerivativesComputation {NO_DERIVATIVES_COMPUTATION = 0, ELEMENTAL_DERIVATIVES = 1, NODAL_ELEMENTAL_DERIVATIVES = 2, NO_DERIVATIVES_COMPUTATION_WITH_NORMAL_UPDATE = 3};

///@}
///@name  Functions
///@{

// VARIABLES
/* Mortar method */
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, int , INTEGRATION_ORDER_CONTACT )                           // The integration order computed in the contact
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, double, DISTANCE_THRESHOLD )                                // The distance threshold considered
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, double, ACTIVE_CHECK_FACTOR )                               // The factor employed to search an active/inactive node
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, GeometryType::Pointer, PAIRED_GEOMETRY )                    // The paired geometry with the current condition
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, Vector3, PAIRED_NORMAL )                                    // The normal of the paired geometry
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(CONTACT_STRUCTURAL_MECHANICS_APPLICATION, AUXILIAR_COORDINATES )                    // Auxiliar coordinates used to map
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(CONTACT_STRUCTURAL_MECHANICS_APPLICATION, DELTA_COORDINATES )                       // Delta coordinates used to map
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, double, NORMAL_GAP )                                        // The normal gap employed in contact formulation

/* Weighted values */
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, double, WEIGHTED_GAP )                                      // The integrated gap employed in mortar formulation
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(CONTACT_STRUCTURAL_MECHANICS_APPLICATION, WEIGHTED_SLIP )                           // The integrated slip employed in mortar formulation
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, double, WEIGHTED_SCALAR_RESIDUAL )                          // The integrated scalar residual
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(CONTACT_STRUCTURAL_MECHANICS_APPLICATION, WEIGHTED_VECTOR_RESIDUAL )                // The integrated vector residual

/* For ALM mortar condition */
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, bool, ACTIVE_SET_CONVERGED )                                // To know if the active set has converged
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, double, DYNAMIC_FACTOR )                                    // The factor considered for dynamic problems (in order to take intro account the gap evolution)
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, double, LAGRANGE_MULTIPLIER_CONTACT_PRESSURE )              // The lagrange multiplier for normal contact pressure
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, double, AUGMENTED_NORMAL_CONTACT_PRESSURE )                 // The resultant augmented pressure in the normal direction
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(CONTACT_STRUCTURAL_MECHANICS_APPLICATION, AUGMENTED_TANGENT_CONTACT_PRESSURE )      // The resultant augmented pressure in the tangent direction
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, int, CONSIDER_NORMAL_VARIATION )                            // A value used to check if consider normal variation or not
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, bool, ADAPT_PENALTY )                                       // To set if the penalty is recalculated or not
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, double, MAX_GAP_FACTOR )                                    // The factor between the nodal H and the max gap considered to recalculate the penalty

/* For mesh tying mortar condition */
KRATOS_DEFINE_APPLICATION_VARIABLE( CONTACT_STRUCTURAL_MECHANICS_APPLICATION, std::string, TYING_VARIABLE )                               // The variable name for the mesh tying
}

#endif	/* KRATOS_CONTACT_STRUCTURAL_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED */
