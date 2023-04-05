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
#include "custom_elements/axisym_small_displacement.h"
#include "custom_utilities/structural_mechanics_math_utilities.hpp"

namespace Kratos
{

//******************************* CONSTRUCTOR ****************************************
//************************************************************************************

AxisymSmallDisplacement::AxisymSmallDisplacement(
    IndexType NewId,
    GeometryType::Pointer pGeometry
    )
        : SmallDisplacement( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

AxisymSmallDisplacement::AxisymSmallDisplacement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties
    )
        : SmallDisplacement( NewId, pGeometry, pProperties )
{
}

//********************************* CREATE *******************************************
//************************************************************************************

Element::Pointer AxisymSmallDisplacement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<AxisymSmallDisplacement>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

//************************************************************************************
//************************************************************************************

Element::Pointer AxisymSmallDisplacement::Create(
    IndexType NewId,
        GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<AxisymSmallDisplacement>( NewId, pGeom, pProperties );
}

//******************************* DESTRUCTOR *****************************************
//************************************************************************************

AxisymSmallDisplacement::~AxisymSmallDisplacement()
{
}

//************************************************************************************
//********************************* PROTECTED ****************************************
//************************************************************************************

void AxisymSmallDisplacement::CalculateB(
    Matrix& rB,
    const Matrix& DN_DX,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const IndexType PointNumber
    ) const
{
    KRATOS_TRY;

    const unsigned int NumberOfNodes = GetGeometry().PointsNumber();

    Vector N;
    N = GetGeometry().ShapeFunctionsValues( N, IntegrationPoints[PointNumber].Coordinates() );
    const double Radius = StructuralMechanicsMathUtilities::CalculateRadius(N, GetGeometry());

    rB.clear();

    for ( unsigned int i = 0; i < NumberOfNodes; i++ )
    {
        const unsigned int index = 2 * i;

        rB(0, index + 0) = DN_DX(i, 0);
        rB(1, index + 1) = DN_DX(i, 1);
        rB(2, index + 0) = N[i]/Radius;
        rB(3, index + 0) = DN_DX(i, 1);
        rB(3, index + 1) = DN_DX(i, 0);
    }

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void AxisymSmallDisplacement::ComputeEquivalentF(
    Matrix& rF,
    const Vector& rStrainVector
    ) const
{
    rF(0,0) = 1.0 + rStrainVector[0];
    rF(0,1) = 0.5 * rStrainVector[3];
    rF(1,0) = 0.5 * rStrainVector[3];
    rF(1,1) = 1.0 + rStrainVector[1];
}

//************************************************************************************
//************************************************************************************

double AxisymSmallDisplacement::GetIntegrationWeight(
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const IndexType PointNumber,
    const double detJ
    ) const
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

void AxisymSmallDisplacement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, SmallDisplacement );
}

//************************************************************************************
//************************************************************************************

void AxisymSmallDisplacement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, SmallDisplacement );
}

} // Namespace Kratos


