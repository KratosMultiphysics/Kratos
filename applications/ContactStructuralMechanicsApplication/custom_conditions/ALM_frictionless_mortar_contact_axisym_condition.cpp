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
#include "custom_conditions/ALM_frictionless_mortar_contact_axisym_condition.h"

/* Utilities */

namespace Kratos 
{
/************************************* OPERATIONS **********************************/
/***********************************************************************************/

template< unsigned int TNumNodes, bool TNormalVariation >
Condition::Pointer AugmentedLagrangianMethodFrictionlessMortarContactAxisymCondition<TNumNodes,TNormalVariation>::Create( 
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesPointerType pProperties ) const
{
    return boost::make_shared< AugmentedLagrangianMethodFrictionlessMortarContactAxisymCondition<TNumNodes, TNormalVariation> >( NewId, this->GetGeometry().Create( rThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TNumNodes, bool TNormalVariation >
Condition::Pointer AugmentedLagrangianMethodFrictionlessMortarContactAxisymCondition<TNumNodes,TNormalVariation>::Create(
    IndexType NewId,
    GeometryPointerType pGeom,
    PropertiesPointerType pProperties) const
{
    return boost::make_shared< AugmentedLagrangianMethodFrictionlessMortarContactAxisymCondition<TNumNodes, TNormalVariation> >( NewId, pGeom, pProperties );
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

template< unsigned int TNumNodes, bool TNormalVariation >
AugmentedLagrangianMethodFrictionlessMortarContactAxisymCondition<TNumNodes, TNormalVariation>::~AugmentedLagrangianMethodFrictionlessMortarContactAxisymCondition( )
= default;

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TNumNodes, bool TNormalVariation >
double AugmentedLagrangianMethodFrictionlessMortarContactAxisymCondition<TNumNodes,TNormalVariation>::GetIntegrationWeight(
    GeneralVariables& rVariables,
    const GeometryType::IntegrationPointsArrayType& ThisIntegrationMethod,
    const unsigned int& PointNumber
    )
{
    const double Radius = CalculateRadius(rVariables);
    const double Thickness = (this->GetValue(ELEMENT_POINTER))->GetProperties()[THICKNESS];
    const double AxiSymCoefficient = 2.0 * M_PI * Radius/Thickness;
    return MortarBaseType::GetIntegrationWeight(rVariables, ThisIntegrationMethod, PointNumber) * AxiSymCoefficient;
}

/*************************COMPUTE AXYSIMMETRIC RADIUS*******************************/
/***********************************************************************************/

template< unsigned int TNumNodes, bool TNormalVariation >
double AugmentedLagrangianMethodFrictionlessMortarContactAxisymCondition<TNumNodes,TNormalVariation>::CalculateRadius(GeneralVariables& rVariables)
{
    KRATOS_TRY;

    double CurrentRadius = 0.0;
//     double ReferenceRadius = 0.0;

    for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
    {
        // Displacement from the reference to the current configuration
//         const array_1d<double, 3 > DeltaDisplacement = this->GetGeometry()[iNode].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[iNode].FastGetSolutionStepValue(DISPLACEMENT,1);  
	    const array_1d<double, 3 > CurrentPosition = this->GetGeometry()[iNode].Coordinates();
// 	    const array_1d<double, 3 > ReferencePosition = CurrentPosition - DeltaDisplacement;
	    
	    CurrentRadius   += CurrentPosition[0] * rVariables.NSlave[iNode];
// 	    ReferenceRadius += ReferencePosition[0] * rVariables.NSlave[iNode];
    }
    
    return CurrentRadius;
//     return ReferenceRadius;
        
    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template class AugmentedLagrangianMethodFrictionlessMortarContactAxisymCondition<2, false>;
template class AugmentedLagrangianMethodFrictionlessMortarContactAxisymCondition<2, true>;

} // Namespace Kratos
