// KRATOS  ___|  |       |       |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//           | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License: BSD License
//   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Vicente Mataix Ferrandiz
//

// System includes
#include <cmath>

// External includes

// Project includes
/* Mortar includes */
#include "custom_conditions/ALM_frictional_mortar_contact_axisym_condition.h"

/* Utilities */

namespace Kratos 
{
/************************************* OPERATIONS **********************************/
/***********************************************************************************/

template< std::size_t TNumNodes, bool TNormalVariation >
Condition::Pointer AugmentedLagrangianMethodFrictionalMortarContactAxisymCondition<TNumNodes,TNormalVariation>::Create( 
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesPointerType pProperties ) const
{
    return Kratos::make_shared<  AugmentedLagrangianMethodFrictionalMortarContactAxisymCondition<TNumNodes, TNormalVariation > >( NewId, this->GetGeometry().Create( rThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TNumNodes, bool TNormalVariation >
Condition::Pointer AugmentedLagrangianMethodFrictionalMortarContactAxisymCondition<TNumNodes,TNormalVariation>::Create(
    IndexType NewId,
    GeometryPointerType pGeom,
    PropertiesPointerType pProperties) const
{
    return Kratos::make_shared< AugmentedLagrangianMethodFrictionalMortarContactAxisymCondition<TNumNodes, TNormalVariation> >( NewId, pGeom, pProperties );
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

template< std::size_t TNumNodes, bool TNormalVariation >
AugmentedLagrangianMethodFrictionalMortarContactAxisymCondition<TNumNodes, TNormalVariation>::~AugmentedLagrangianMethodFrictionalMortarContactAxisymCondition( )
= default;

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TNumNodes, bool TNormalVariation >
double AugmentedLagrangianMethodFrictionalMortarContactAxisymCondition<TNumNodes,TNormalVariation>::GetAxisymmetricCoefficient(const GeneralVariables& rVariables) const
{
    const double radius = CalculateRadius(rVariables);
    const double thickness = this->GetProperties()[THICKNESS];
    return (2.0 * Globals::Pi * radius/thickness);
}

/*************************COMPUTE AXYSIMMETRIC RADIUS*******************************/
/***********************************************************************************/

template< std::size_t TNumNodes, bool TNormalVariation >
double AugmentedLagrangianMethodFrictionalMortarContactAxisymCondition<TNumNodes,TNormalVariation>::CalculateRadius(const GeneralVariables& rVariables) const
{
    KRATOS_TRY;

    double current_radius = 0.0;
//     double reference_radius = 0.0;

    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node)
    {
        // Displacement from the reference to the current configuration
//         const array_1d<double, 3 > DeltaDisplacement = this->GetGeometry()[i_node].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[i_node].FastGetSolutionStepValue(DISPLACEMENT,1);  
	    const array_1d<double, 3 > current_position = this->GetGeometry()[i_node].Coordinates();
// 	    const array_1d<double, 3 > ReferencePosition = current_position - DeltaDisplacement;
	    
	    current_radius   += current_position[0] * rVariables.NSlave[i_node];
// 	    reference_radius += ReferencePosition[0] * rVariables.NSlave[i_node];
    }
    
    return current_radius;
//     return reference_radius;
        
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template class AugmentedLagrangianMethodFrictionalMortarContactAxisymCondition<2, false>;
template class AugmentedLagrangianMethodFrictionalMortarContactAxisymCondition<2, true>;

} // Namespace Kratos
