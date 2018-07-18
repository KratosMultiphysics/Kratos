// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "custom_elements/axisym_total_lagrangian.h"
#include "custom_utilities/structural_mechanics_math_utilities.hpp"

namespace Kratos
{
//******************************* CONSTRUCTOR ****************************************
//************************************************************************************

AxisymTotalLagrangian::AxisymTotalLagrangian(
    IndexType NewId,
    GeometryType::Pointer pGeometry
    )
        : TotalLagrangian( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

AxisymTotalLagrangian::AxisymTotalLagrangian(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties
    )
        : TotalLagrangian( NewId, pGeometry, pProperties )
{
}

//********************************* CREATE *******************************************
//************************************************************************************

Element::Pointer AxisymTotalLagrangian::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_shared<AxisymTotalLagrangian>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

//************************************************************************************
//************************************************************************************

Element::Pointer AxisymTotalLagrangian::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_shared<AxisymTotalLagrangian>( NewId, pGeom, pProperties );
}

//******************************* DESTRUCTOR *****************************************
//************************************************************************************

AxisymTotalLagrangian::~AxisymTotalLagrangian()
{
}

//************************************************************************************
//********************************* PROTECTED ****************************************
//************************************************************************************

double AxisymTotalLagrangian::GetIntegrationWeight(
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const IndexType PointNumber,
    const double detJ
    )
{
    // We calculate the axisymmetric coefficient
    Vector N;
    N = GetGeometry().ShapeFunctionsValues( N, IntegrationPoints[PointNumber].Coordinates() );
    const double radius = StructuralMechanicsMathUtilities::CalculateRadius(N, GetGeometry());
    const double thickness = (GetProperties().Has( THICKNESS ) == true) ? this->GetProperties()[THICKNESS] : 1.0;
    const double axi_sym_coefficient = 2.0 * Globals::Pi * radius/thickness;

    return axi_sym_coefficient * IntegrationPoints[PointNumber].Weight() * detJ;
}

//************************************************************************************
//************************************************************************************

void AxisymTotalLagrangian::save( Serializer& rSerializer ) const
{
    rSerializer.save( "Name", "AxisymTotalLagrangian" );
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, TotalLagrangian );
}

//************************************************************************************
//************************************************************************************

void AxisymTotalLagrangian::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, TotalLagrangian );
}

} // Namespace Kratos


