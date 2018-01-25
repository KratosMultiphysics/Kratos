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
#include "includes/model_part.h"
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
    
    enum NormalDerivativesComputation {NO_DERIVATIVES_COMPUTATION = 0, ELEMENTAL_DERIVATIVES = 1, NODAL_ELEMENTAL_DERIVATIVES = 2};
    
///@}
///@name  Functions
///@{

// VARIABLES
/* Mortar method */ 
KRATOS_DEFINE_VARIABLE( int , INTEGRATION_ORDER_CONTACT )                            // The integration order computed in the contact
KRATOS_DEFINE_VARIABLE( double, ACTIVE_CHECK_FACTOR )                                // The factor employed to search an active/inactive node
KRATOS_DEFINE_VARIABLE( GeometryType::Pointer, PAIRED_GEOMETRY )                     // The paired geometry with the current condition
KRATOS_DEFINE_VARIABLE( Vector3, PAIRED_NORMAL )                                     // The normal of the paired geometry
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( AUXILIAR_COORDINATES )                    // Auxiliar coordinates used to map
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( DELTA_COORDINATES )                       // Delta coordinates used to map
KRATOS_DEFINE_VARIABLE( double, NORMAL_GAP )                                         // The normal gap employed in contact formulation

/* Weighted values */
KRATOS_DEFINE_VARIABLE( double, WEIGHTED_GAP )                                       // The integrated gap employed in mortar formulation
KRATOS_DEFINE_VARIABLE( double, WEIGHTED_SLIP )                                      // The integrated slip employed in mortar formulation
KRATOS_DEFINE_VARIABLE( double, WEIGHTED_FRICTION )                                  // The integrated friction coefficient
KRATOS_DEFINE_VARIABLE( double, WEIGHTED_SCALAR_RESIDUAL )                           // The integrated scalar residual  
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( WEIGHTED_VECTOR_RESIDUAL )                // The integrated vector residual         

/* For ALM mortar condition */
KRATOS_DEFINE_VARIABLE( double, AUGMENTED_NORMAL_CONTACT_PRESSURE )                 // The resultant augmented pressure in the normal direction
KRATOS_DEFINE_VARIABLE( double, AUGMENTED_TANGENT_CONTACT_PRESSURE )                // The resultant augmented pressure in the tangent direction
KRATOS_DEFINE_VARIABLE( double, TANGENT_FACTOR )                                    // The proportion between the tangent and normal penalty
KRATOS_DEFINE_VARIABLE( int, CONSIDER_NORMAL_VARIATION )                            // A value used to check if consider normal variation or not
KRATOS_DEFINE_VARIABLE( bool, ADAPT_PENALTY )                                       // To set if the penalty is recalculated or not
KRATOS_DEFINE_VARIABLE( double, MAX_GAP_FACTOR )                                    // The factor between the nodal H and the max gap considered to recalculate the penalty

/* For mesh tying mortar condition */
KRATOS_DEFINE_VARIABLE( std::string, TYING_VARIABLE )                               // The variable name for the mesh tying  
}       

#endif	/* KRATOS_CONTACT_STRUCTURAL_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED */
