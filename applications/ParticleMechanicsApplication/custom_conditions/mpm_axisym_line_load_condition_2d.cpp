//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


// System includes

// External includes

// Project includes
#include "custom_conditions/mpm_axisym_line_load_condition_2d.h"
#include "custom_utilities/particle_mechanics_math_utilities.hpp"

namespace Kratos
{
/******************************* CONSTRUCTOR ***************************************/
/***********************************************************************************/

MPMAxisymLineLoadCondition2D::MPMAxisymLineLoadCondition2D(
    IndexType NewId,
    GeometryType::Pointer pGeometry
    )
        : MPMLineLoadCondition2D( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

MPMAxisymLineLoadCondition2D::MPMAxisymLineLoadCondition2D(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties
    )
        : MPMLineLoadCondition2D( NewId, pGeometry, pProperties )
{
}

/********************************* CREATE ******************************************/
/***********************************************************************************/

Condition::Pointer MPMAxisymLineLoadCondition2D::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_shared<MPMAxisymLineLoadCondition2D>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer MPMAxisymLineLoadCondition2D::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_shared<MPMAxisymLineLoadCondition2D>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/******************************* DESTRUCTOR ****************************************/
/***********************************************************************************/

MPMAxisymLineLoadCondition2D::~MPMAxisymLineLoadCondition2D()
{
}

/***********************************************************************************/
/********************************* PROTECTED ***************************************/
/***********************************************************************************/

double MPMAxisymLineLoadCondition2D::GetIntegrationWeight(
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const SizeType PointNumber,
    const double detJ
    )
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

void MPMAxisymLineLoadCondition2D::save( Serializer& rSerializer ) const
{
    rSerializer.save( "Name", "MPMAxisymLineLoadCondition2D" );
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MPMLineLoadCondition2D );
}

/***********************************************************************************/
/***********************************************************************************/

void MPMAxisymLineLoadCondition2D::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MPMLineLoadCondition2D );
}

} // Namespace Kratos


