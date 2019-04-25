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
#include "custom_conditions/penalty_frictional_mortar_contact_axisym_condition.h"

/* Utilities */

namespace Kratos
{
/************************************* OPERATIONS **********************************/
/***********************************************************************************/

template< std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster >
Condition::Pointer PenaltyMethodFrictionalMortarContactAxisymCondition<TNumNodes,TNormalVariation, TNumNodesMaster>::Create(
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesPointerType pProperties ) const
{
    return Kratos::make_shared< PenaltyMethodFrictionalMortarContactAxisymCondition<TNumNodes, TNormalVariation, TNumNodesMaster > >( NewId, this->GetGeometry().Create( rThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster >
Condition::Pointer PenaltyMethodFrictionalMortarContactAxisymCondition<TNumNodes,TNormalVariation, TNumNodesMaster>::Create(
    IndexType NewId,
    GeometryPointerType pGeom,
    PropertiesPointerType pProperties) const
{
    return Kratos::make_shared< PenaltyMethodFrictionalMortarContactAxisymCondition<TNumNodes, TNormalVariation, TNumNodesMaster> >( NewId, pGeom, pProperties );
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

template< std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster >
PenaltyMethodFrictionalMortarContactAxisymCondition<TNumNodes, TNormalVariation, TNumNodesMaster>::~PenaltyMethodFrictionalMortarContactAxisymCondition( )
= default;

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster >
bool PenaltyMethodFrictionalMortarContactAxisymCondition<TNumNodes,TNormalVariation, TNumNodesMaster>::IsAxisymmetric() const
{
    return true;
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster >
double PenaltyMethodFrictionalMortarContactAxisymCondition<TNumNodes,TNormalVariation, TNumNodesMaster>::GetAxisymmetricCoefficient(const GeneralVariables& rVariables) const
{
    const double radius = CalculateRadius(rVariables);
    const double thickness = this->GetProperties()[THICKNESS];
    return (2.0 * Globals::Pi * radius/thickness);
}

/*************************COMPUTE AXYSIMMETRIC RADIUS*******************************/
/***********************************************************************************/

template< std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster >
double PenaltyMethodFrictionalMortarContactAxisymCondition<TNumNodes,TNormalVariation, TNumNodesMaster>::CalculateRadius(const GeneralVariables& rVariables) const
{
    KRATOS_TRY;

    double current_radius = 0.0;

    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        // Displacement from the reference to the current configuration
        const array_1d<double, 3 >& r_current_position = this->GetGeometry()[i_node].Coordinates();

        current_radius += r_current_position[0] * rVariables.NSlave[i_node];
    }

    return current_radius;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template class PenaltyMethodFrictionalMortarContactAxisymCondition<2, false>;
template class PenaltyMethodFrictionalMortarContactAxisymCondition<2, true>;

} // Namespace Kratos
