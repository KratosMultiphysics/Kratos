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
#include "custom_elements/axisym_updated_lagrangian.h"
#include "custom_utilities/structural_mechanics_math_utilities.hpp"

namespace Kratos
{
//******************************* CONSTRUCTOR ****************************************
//************************************************************************************

AxisymUpdatedLagrangian::AxisymUpdatedLagrangian(
    IndexType NewId,
    GeometryType::Pointer pGeometry
    )
        : UpdatedLagrangian( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

AxisymUpdatedLagrangian::AxisymUpdatedLagrangian(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties
    )
        : UpdatedLagrangian( NewId, pGeometry, pProperties )
{
}

//********************************* CREATE *******************************************
//************************************************************************************

Element::Pointer AxisymUpdatedLagrangian::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_shared<AxisymUpdatedLagrangian>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

//************************************************************************************
//************************************************************************************

Element::Pointer AxisymUpdatedLagrangian::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_shared<AxisymUpdatedLagrangian>( NewId, pGeom, pProperties );
}

//******************************* DESTRUCTOR *****************************************
//************************************************************************************

AxisymUpdatedLagrangian::~AxisymUpdatedLagrangian()
{
}

//************************************************************************************
//********************************* PROTECTED ****************************************
//************************************************************************************

double AxisymUpdatedLagrangian::GetIntegrationWeight(
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

void AxisymUpdatedLagrangian::save( Serializer& rSerializer ) const
{
    rSerializer.save( "Name", "AxisymUpdatedLagrangian" );
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, UpdatedLagrangian );
}

//************************************************************************************
//************************************************************************************

void AxisymUpdatedLagrangian::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, UpdatedLagrangian );
}

} // Namespace Kratos


