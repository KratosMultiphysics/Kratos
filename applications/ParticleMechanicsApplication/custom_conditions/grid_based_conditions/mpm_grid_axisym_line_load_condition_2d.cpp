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
#include "custom_conditions/grid_based_conditions/mpm_grid_axisym_line_load_condition_2d.h"
#include "custom_utilities/particle_mechanics_math_utilities.h"

namespace Kratos
{
/******************************* CONSTRUCTOR ***************************************/
/***********************************************************************************/

MPMGridAxisymLineLoadCondition2D::MPMGridAxisymLineLoadCondition2D(
    IndexType NewId,
    GeometryType::Pointer pGeometry
    )
        : MPMGridLineLoadCondition2D( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

MPMGridAxisymLineLoadCondition2D::MPMGridAxisymLineLoadCondition2D(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties
    )
        : MPMGridLineLoadCondition2D( NewId, pGeometry, pProperties )
{
}

/********************************* CREATE ******************************************/
/***********************************************************************************/

Condition::Pointer MPMGridAxisymLineLoadCondition2D::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_shared<MPMGridAxisymLineLoadCondition2D>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer MPMGridAxisymLineLoadCondition2D::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_shared<MPMGridAxisymLineLoadCondition2D>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/******************************* DESTRUCTOR ****************************************/
/***********************************************************************************/

MPMGridAxisymLineLoadCondition2D::~MPMGridAxisymLineLoadCondition2D()
{
}

/***********************************************************************************/
/********************************* PROTECTED ***************************************/
/***********************************************************************************/

double MPMGridAxisymLineLoadCondition2D::GetIntegrationWeight(
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const unsigned int PointNumber,
    const double detJ
    )
{
    // We calculate the axisymmetric coefficient
    Vector N;
    N = GetGeometry().ShapeFunctionsValues( N, IntegrationPoints[PointNumber].Coordinates() );
    const double Radius = ParticleMechanicsMathUtilities<double>::CalculateRadius(N, GetGeometry());
    const double Thickness = (GetProperties().Has( THICKNESS ) == true) ? this->GetProperties()[THICKNESS] : 1.0;
    const double AxiSymCoefficient = 2.0 * Globals::Pi * Radius/Thickness;

    return AxiSymCoefficient * IntegrationPoints[PointNumber].Weight() * detJ;
}


/***********************************************************************************/
/***********************************************************************************/

void MPMGridAxisymLineLoadCondition2D::save( Serializer& rSerializer ) const
{
    rSerializer.save( "Name", "MPMGridAxisymLineLoadCondition2D" );
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MPMGridLineLoadCondition2D );
}

/***********************************************************************************/
/***********************************************************************************/

void MPMGridAxisymLineLoadCondition2D::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MPMGridLineLoadCondition2D );
}

} // Namespace Kratos


