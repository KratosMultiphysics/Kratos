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
#include "custom_conditions/grid_based_conditions/mpm_grid_axisym_point_load_condition.h"
#include "custom_utilities/particle_mechanics_math_utilities.h"

namespace Kratos
{
//******************************* CONSTRUCTOR ****************************************
//************************************************************************************

MPMGridAxisymPointLoadCondition::MPMGridAxisymPointLoadCondition(
    IndexType NewId,
    GeometryType::Pointer pGeometry
    )
        : MPMGridPointLoadCondition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

MPMGridAxisymPointLoadCondition::MPMGridAxisymPointLoadCondition(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties
    )
        : MPMGridPointLoadCondition( NewId, pGeometry, pProperties )
{
}

//********************************* CREATE *******************************************
//************************************************************************************

Condition::Pointer MPMGridAxisymPointLoadCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_shared<MPMGridAxisymPointLoadCondition>( NewId, pGeom, pProperties );
}

//************************************************************************************
//************************************************************************************

Condition::Pointer MPMGridAxisymPointLoadCondition::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_shared<MPMGridAxisymPointLoadCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

//******************************* DESTRUCTOR *****************************************
//************************************************************************************

MPMGridAxisymPointLoadCondition::~MPMGridAxisymPointLoadCondition()
{
}

//************************************************************************************
//********************************* PROTECTED ****************************************
//************************************************************************************

double MPMGridAxisymPointLoadCondition::GetPointLoadIntegrationWeight()
{
    // We calculate the axisymmetric coefficient
    const double Radius = ParticleMechanicsMathUtilities<double>::CalculateRadiusPoint(GetGeometry());
    const double Thickness = (GetProperties().Has( THICKNESS ) == true) ? this->GetProperties()[THICKNESS] : 1.0;
    const double AxiSymCoefficient = 2.0 * Globals::Pi * Radius/Thickness;

    return AxiSymCoefficient;
}


//************************************************************************************
//************************************************************************************

void MPMGridAxisymPointLoadCondition::save( Serializer& rSerializer ) const
{
    rSerializer.save( "Name", "MPMGridAxisymPointLoadCondition" );
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MPMGridPointLoadCondition );
}

//************************************************************************************
//************************************************************************************

void MPMGridAxisymPointLoadCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MPMGridPointLoadCondition );
}

} // Namespace Kratos


