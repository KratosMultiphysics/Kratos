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
#include "custom_conditions/particle_based_conditions/mpm_particle_axisym_penalty_dirichlet_condition_2d.h"

namespace Kratos
{
/******************************* CONSTRUCTOR ***************************************/
/***********************************************************************************/

MPMParticleAxisymPenaltyDirichletCondition2D::MPMParticleAxisymPenaltyDirichletCondition2D(
    IndexType NewId,
    GeometryType::Pointer pGeometry
    )
        : MPMParticlePenaltyDirichletCondition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

MPMParticleAxisymPenaltyDirichletCondition2D::MPMParticleAxisymPenaltyDirichletCondition2D(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties
    )
        : MPMParticlePenaltyDirichletCondition( NewId, pGeometry, pProperties )
{
}

/********************************* CREATE ******************************************/
/***********************************************************************************/

Condition::Pointer MPMParticleAxisymPenaltyDirichletCondition2D::Create(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<MPMParticleAxisymPenaltyDirichletCondition2D>(NewId, pGeometry, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer MPMParticleAxisymPenaltyDirichletCondition2D::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<MPMParticleAxisymPenaltyDirichletCondition2D>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/******************************* DESTRUCTOR ****************************************/
/***********************************************************************************/

MPMParticleAxisymPenaltyDirichletCondition2D::~MPMParticleAxisymPenaltyDirichletCondition2D()
{
}

/***********************************************************************************/
/********************************* PROTECTED ***************************************/
/***********************************************************************************/

double MPMParticleAxisymPenaltyDirichletCondition2D::GetIntegrationWeight()
{
    const array_1d<double,3>& xg_c = this->GetValue(MPC_COORD);
    const double integration_weight = MPMParticlePenaltyDirichletCondition::GetIntegrationWeight();
    const double thickness = (GetProperties().Has( THICKNESS ) == true) ? this->GetProperties()[THICKNESS] : 1.0;
    const double axis_symmetric_weight = integration_weight * 2.0 * Globals::Pi * xg_c[0] / thickness;

    return axis_symmetric_weight;
}


/***********************************************************************************/
/***********************************************************************************/

void MPMParticleAxisymPenaltyDirichletCondition2D::save( Serializer& rSerializer ) const
{
    rSerializer.save( "Name", "MPMParticleAxisymPenaltyDirichletCondition2D" );
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MPMParticlePenaltyDirichletCondition );
}

/***********************************************************************************/
/***********************************************************************************/

void MPMParticleAxisymPenaltyDirichletCondition2D::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MPMParticlePenaltyDirichletCondition );
}

} // Namespace Kratos


