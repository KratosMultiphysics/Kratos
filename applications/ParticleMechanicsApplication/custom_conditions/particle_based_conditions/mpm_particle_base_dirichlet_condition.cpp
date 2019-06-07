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
#include "custom_conditions/particle_based_conditions/mpm_particle_base_dirichlet_condition.h"

namespace Kratos
{
//************************************************************************************
//************************************************************************************

void MPMParticleBaseDirichletCondition::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    /* NOTE:
    In the InitializeSolutionStep of each time step the nodal initial conditions are evaluated.
    This function is called by the base scheme class.*/
    // Here MPC_DISPLACEMENT is updated in terms of velocity and acceleration is added
    array_1d<double,3>& MPC_Displacement = this->GetValue(MPC_DISPLACEMENT);
    array_1d<double,3>& MPC_Velocity = this->GetValue(MPC_VELOCITY);
    const array_1d<double,3>& MPC_Acceleration = this->GetValue(MPC_ACCELERATION);
    const double& delta_time = rCurrentProcessInfo[DELTA_TIME];

    MPC_Displacement += (MPC_Velocity * delta_time) + (0.5 * MPC_Acceleration * delta_time * delta_time);
    MPC_Velocity += (MPC_Acceleration * delta_time);
}

//************************************************************************************
//************************************************************************************

void MPMParticleBaseDirichletCondition::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // Update the MPC Position
    const array_1d<double,3> & xg_c = this->GetValue(MPC_COORD);
    array_1d<double,3> & displacement = this->GetValue(MPC_DISPLACEMENT);
    const array_1d<double,3> & new_xg_c = xg_c + displacement ;
    this->SetValue(MPC_COORD,new_xg_c);

    // Set displacement to zero (NOTE: to use incremental displacement, use MPC_VELOCITY)
    displacement.clear();

    KRATOS_CATCH( "" )
}

} // Namespace Kratos


