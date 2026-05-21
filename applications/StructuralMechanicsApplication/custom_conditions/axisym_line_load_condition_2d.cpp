// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "custom_conditions/axisym_line_load_condition_2d.h"
#include "custom_utilities/structural_mechanics_math_utilities.hpp"

namespace Kratos
{
/******************************* CONSTRUCTOR ***************************************/
/***********************************************************************************/

AxisymLineLoadCondition2D::AxisymLineLoadCondition2D(
    IndexType NewId,
    GeometryType::Pointer pGeometry
    ) : LineLoadCondition<2>( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

AxisymLineLoadCondition2D::AxisymLineLoadCondition2D(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties
    ) : LineLoadCondition<2>( NewId, pGeometry, pProperties )
{
}

/********************************* CREATE ******************************************/
/***********************************************************************************/

Condition::Pointer AxisymLineLoadCondition2D::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<AxisymLineLoadCondition2D>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer AxisymLineLoadCondition2D::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<AxisymLineLoadCondition2D>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer AxisymLineLoadCondition2D::Clone (
    IndexType NewId,
    NodesArrayType const& ThisNodes
    ) const
{
    KRATOS_TRY

    Condition::Pointer p_new_cond = Kratos::make_intrusive<AxisymLineLoadCondition2D>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    p_new_cond->SetData(this->GetData());
    p_new_cond->Set(Flags(*this));
    return p_new_cond;

    KRATOS_CATCH("");
}

/******************************* DESTRUCTOR ****************************************/
/***********************************************************************************/

AxisymLineLoadCondition2D::~AxisymLineLoadCondition2D()
{
}

/***********************************************************************************/
/********************************* PROTECTED ***************************************/
/***********************************************************************************/

double AxisymLineLoadCondition2D::GetIntegrationWeight(
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const SizeType PointNumber,
    const double detJ
    ) const
{
    // We calculate the axisymmetric coefficient
    Vector N;
    N = GetGeometry().ShapeFunctionsValues( N, IntegrationPoints[PointNumber].Coordinates() );
    const double Radius = StructuralMechanicsMathUtilities::CalculateRadius(N, GetGeometry());
    const double Thickness = (GetProperties().Has( THICKNESS ) == true) ? this->GetProperties()[THICKNESS] : 1.0;
    const double AxiSymCoefficient = 2.0 * Globals::Pi * Radius/Thickness;

    return AxiSymCoefficient * IntegrationPoints[PointNumber].Weight() * detJ;
}


/***********************************************************************************/
/***********************************************************************************/

void AxisymLineLoadCondition2D::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LineLoadCondition<2> );
}

/***********************************************************************************/
/***********************************************************************************/

void AxisymLineLoadCondition2D::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LineLoadCondition<2> );
}

} // Namespace Kratos


