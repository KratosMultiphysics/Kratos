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
#include "includes/checks.h"

namespace Kratos
{

void MPMParticleBaseDirichletCondition::InitializeSolutionStep( const ProcessInfo& rCurrentProcessInfo )
{
    /* NOTE:
    In the InitializeSolutionStep of each time step the nodal initial conditions are evaluated.
    This function is called by the base scheme class.*/
    // Here MPC_IMPOSED_DISPLACEMENT is updated in terms of velocity and acceleration is added
    const double& delta_time = rCurrentProcessInfo[DELTA_TIME];

    // Convert imposition of velocity and acceleration to displacement
    // NOTE: This only consider translational velocity and acceleration: no angular
    m_imposed_displacement += (m_imposed_velocity * delta_time) + (0.5 * m_imposed_acceleration * delta_time * delta_time);

    // Prepare variables
    GeneralVariables Variables;
    MPMShapeFunctionPointValues(Variables.N);

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

void MPMParticleBaseDirichletCondition::FinalizeSolutionStep( const ProcessInfo& rCurrentProcessInfo )
{
    // Update the MPC Position
    m_xg += m_imposed_displacement;

    // Update total MPC Displacement
    m_displacement += m_imposed_displacement;

    m_imposed_displacement = ZeroVector(3);
}

void MPMParticleBaseDirichletCondition::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3 > >& rVariable,
    std::vector<array_1d<double, 3 > >& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);

    if (rVariable == MPC_IMPOSED_DISPLACEMENT) {
        rValues[0] = m_imposed_displacement;
    }
    else if (rVariable == MPC_IMPOSED_VELOCITY) {
        rValues[0] = m_imposed_velocity;
    }
    else if (rVariable == MPC_IMPOSED_ACCELERATION) {
        rValues[0] = m_imposed_acceleration;
    }
    else if (rVariable == MPC_CONTACT_FORCE) {
        rValues[0] = m_contact_force;
    }
    else {
        MPMParticleBaseCondition::CalculateOnIntegrationPoints(
            rVariable, rValues, rCurrentProcessInfo);
    }
}

void MPMParticleBaseDirichletCondition::SetValuesOnIntegrationPoints(
    const Variable<array_1d<double, 3 > >& rVariable,
    const std::vector<array_1d<double, 3 > >& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR_IF(rValues.size() > 1)
        << "Only 1 value per integration point allowed! Passed values vector size: "
        << rValues.size() << std::endl;

    if (rVariable == MPC_IMPOSED_DISPLACEMENT) {
        m_imposed_displacement = rValues[0];
    }
    else if (rVariable == MPC_IMPOSED_VELOCITY) {
        m_imposed_velocity = rValues[0];
    }
    else if (rVariable == MPC_IMPOSED_ACCELERATION) {
        m_imposed_acceleration = rValues[0];
    }
    else if (rVariable == MPC_CONTACT_FORCE) {
        m_contact_force = rValues[0];
    }
    else {
        MPMParticleBaseCondition::SetValuesOnIntegrationPoints(
            rVariable, rValues, rCurrentProcessInfo);
    }
}

void MPMParticleBaseDirichletCondition::MPMShapeFunctionPointValues( Vector& rResult ) const
{
    KRATOS_TRY

    MPMParticleBaseCondition::MPMShapeFunctionPointValues(rResult);

    // Additional check to modify zero shape function values
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();

    double denominator = 1.0;
    const double small_cut_instability_tolerance = 0.01;
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        if (rResult[i] < small_cut_instability_tolerance){
            denominator += (small_cut_instability_tolerance - rResult[i]);
            rResult[i] = small_cut_instability_tolerance;
        }
    }

    rResult = rResult / denominator;

    KRATOS_CATCH( "" )
}

int MPMParticleBaseDirichletCondition::Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    MPMParticleBaseCondition::Check(rCurrentProcessInfo);

    // Verify that the dofs exist
    for (const auto& r_node : this->GetGeometry().Points()){
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NORMAL,r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NODAL_AREA,r_node)
    }

    return 0;
}

} // Namespace Kratos


