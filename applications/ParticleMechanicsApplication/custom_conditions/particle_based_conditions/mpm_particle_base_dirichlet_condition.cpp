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
    // Here MPC_IMPOSED_DISPLACEMENT is updated in terms of velocity and acceleration is added
    array_1d<double,3>& MPC_Imposed_Displacement = this->GetValue(MPC_IMPOSED_DISPLACEMENT);
    const array_1d<double,3>& MPC_Imposed_Velocity = this->GetValue(MPC_IMPOSED_VELOCITY);
    const array_1d<double,3>& MPC_Imposed_Acceleration = this->GetValue(MPC_IMPOSED_ACCELERATION);
    const double& delta_time = rCurrentProcessInfo[DELTA_TIME];

    // Convert imposition of velocity and acceleration to displacement
    // NOTE: This only consider translational velocity and acceleration: no angular
    MPC_Imposed_Displacement += (MPC_Imposed_Velocity * delta_time) + (0.5 * MPC_Imposed_Acceleration * delta_time * delta_time);

    // Prepare variables
    GeneralVariables Variables;
    const array_1d<double, 3 > & xg_c = this->GetValue(MPC_COORD);
    Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg_c);

    // Get NODAL_AREA from MPC_Area
    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    const double & r_mpc_area = this->GetIntegrationWeight();
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        if (r_geometry[i].SolutionStepsDataHas(NODAL_AREA))
        {
            r_geometry[i].SetLock();
            r_geometry[i].FastGetSolutionStepValue(NODAL_AREA, 0) += Variables.N[i] * r_mpc_area;
            r_geometry[i].UnSetLock();
        }
        else break;
    }

}

//************************************************************************************
//************************************************************************************

void MPMParticleBaseDirichletCondition::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const array_1d<double,3> & xg_c = this->GetValue(MPC_COORD);
    array_1d<double,3> & delta_xg_c = this->GetValue(MPC_IMPOSED_DISPLACEMENT);

    // Update the MPC Position
    const array_1d<double,3> new_xg_c = xg_c + delta_xg_c;
    this->SetValue(MPC_COORD,new_xg_c);

    // Update total MPC Displacement
    array_1d<double,3> & MPC_Displacement = this->GetValue(MPC_DISPLACEMENT);
    MPC_Displacement += delta_xg_c;
    this->SetValue(MPC_DISPLACEMENT,MPC_Displacement);

    // Reset value of incremental displacement
    delta_xg_c.clear();

    KRATOS_CATCH( "" )
}


Vector& MPMParticleBaseDirichletCondition::MPMShapeFunctionPointValues( Vector& rResult, const array_1d<double,3>& rPoint )
{
    KRATOS_TRY

    rResult = MPMParticleBaseCondition::MPMShapeFunctionPointValues(rResult, rPoint);

    // Additional check to modify zero shape function values
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();

    double denominator = 1.0;
    const double small_cut_instability_tolerance = 0.01;
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        if (rResult[i] < small_cut_instability_tolerance){
            rResult[i] = small_cut_instability_tolerance;
            denominator += rResult[i];
        }
    }

    rResult = rResult/denominator;

    return rResult;

    KRATOS_CATCH( "" )
}

} // Namespace Kratos


