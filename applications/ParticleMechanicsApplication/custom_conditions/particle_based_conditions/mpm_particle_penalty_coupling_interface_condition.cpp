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
#include "includes/define.h"
#include "custom_conditions/particle_based_conditions/mpm_particle_penalty_coupling_interface_condition.h"
#include "includes/kratos_flags.h"

namespace Kratos
{
//******************************* CONSTRUCTOR ****************************************
//************************************************************************************

MPMParticlePenaltyCouplingInterfaceCondition::MPMParticlePenaltyCouplingInterfaceCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : MPMParticlePenaltyDirichletCondition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

MPMParticlePenaltyCouplingInterfaceCondition::MPMParticlePenaltyCouplingInterfaceCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : MPMParticlePenaltyDirichletCondition( NewId, pGeometry, pProperties )
{
}

//********************************* CREATE *******************************************
//************************************************************************************

Condition::Pointer MPMParticlePenaltyCouplingInterfaceCondition::Create(IndexType NewId,GeometryType::Pointer pGeom,PropertiesType::Pointer pProperties) const
{
    return Kratos::make_shared<MPMParticlePenaltyCouplingInterfaceCondition>(NewId, pGeom, pProperties);
}

//************************************************************************************
//************************************************************************************

Condition::Pointer MPMParticlePenaltyCouplingInterfaceCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_shared<MPMParticlePenaltyCouplingInterfaceCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

//******************************* DESTRUCTOR *****************************************
//************************************************************************************

MPMParticlePenaltyCouplingInterfaceCondition::~MPMParticlePenaltyCouplingInterfaceCondition()
{
}

//************************************************************************************
//************************************************************************************

void MPMParticlePenaltyCouplingInterfaceCondition::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    MPMParticlePenaltyDirichletCondition::FinalizeSolutionStep(rCurrentProcessInfo);

    // Estimating the contact forces at the boundary
    if (Is(INTERFACE))
    {
        this->CalculateContactForce( rCurrentProcessInfo );
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void MPMParticlePenaltyCouplingInterfaceCondition::CalculateContactForce( ProcessInfo& rCurrentProcessInfo )
{
    GeometryType& rGeom = GetGeometry();
    const unsigned int number_of_nodes = rGeom.PointsNumber();
    const unsigned int dimension = rGeom.WorkingSpaceDimension();
    const unsigned int block_size = this->GetBlockSize();

    // Prepare variables
    GeneralVariables Variables;
    const array_1d<double, 3 > & xg_c = this->GetValue(MPC_COORD);
    Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg_c);

    // Compute nodal force: this has considered CONTACT condition
    Vector nodal_force;
    this->CalculateRightHandSide(nodal_force, rCurrentProcessInfo);

    // Interpolate the force to MPC_Force assuming linear shape function
    array_1d<double, 3 > MPC_Force = ZeroVector(3);
    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        for ( unsigned int j = 0; j < dimension; j++)
        {
            MPC_Force[j] += Variables.N[i] *  nodal_force[block_size * i + j];
        }
    }
    MPC_Force *= -1.0;

    // Set Contact Force
    this->SetValue(MPC_CONTACT_FORCE, MPC_Force);

}


} // Namespace Kratos


