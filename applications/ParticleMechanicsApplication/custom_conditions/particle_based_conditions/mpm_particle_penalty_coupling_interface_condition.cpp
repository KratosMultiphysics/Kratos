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
#include "utilities/math_utils.h"
#include "custom_utilities/particle_mechanics_math_utilities.h"

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
    return Kratos::make_intrusive<MPMParticlePenaltyCouplingInterfaceCondition>(NewId, pGeom, pProperties);
}

//************************************************************************************
//************************************************************************************

Condition::Pointer MPMParticlePenaltyCouplingInterfaceCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<MPMParticlePenaltyCouplingInterfaceCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

//******************************* DESTRUCTOR *****************************************
//************************************************************************************

MPMParticlePenaltyCouplingInterfaceCondition::~MPMParticlePenaltyCouplingInterfaceCondition()
{
}

//************************************************************************************
//************************************************************************************

void MPMParticlePenaltyCouplingInterfaceCondition::InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
    if (Is(INTERFACE))
    {
        GeometryType& rGeom = GetGeometry();
        const unsigned int number_of_nodes = rGeom.PointsNumber();

        // At the beginning of NonLinearIteration, REACTION has to be reset to zero
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            rGeom[i].SetLock();
            rGeom[i].FastGetSolutionStepValue(REACTION).clear();
            rGeom[i].UnSetLock();
        }

        mReactionIsAdded = false;
    }
}

//************************************************************************************
//************************************************************************************

void MPMParticlePenaltyCouplingInterfaceCondition::FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Estimating the contact forces at the boundary
    if (Is(INTERFACE))
    {
        this->CalculateInterfaceContactForce( rCurrentProcessInfo );
    }

    KRATOS_CATCH( "" )
}


void MPMParticlePenaltyCouplingInterfaceCondition::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    MPMParticlePenaltyDirichletCondition::FinalizeSolutionStep(rCurrentProcessInfo);

    // Estimating the contact forces at the boundary
    if (Is(INTERFACE))
    {
        this->CalculateInterfaceContactForce( rCurrentProcessInfo );
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void MPMParticlePenaltyCouplingInterfaceCondition::CalculateAll(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo,
    bool CalculateStiffnessMatrixFlag,
    bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY

    MPMParticlePenaltyDirichletCondition::CalculateAll(rLeftHandSideMatrix, rRightHandSideVector,
                                                        rCurrentProcessInfo, CalculateStiffnessMatrixFlag,
                                                        CalculateResidualVectorFlag);

    // Append penalty force at the nodes
    if (Is(INTERFACE) && !mReactionIsAdded)
    {
        this->CalculateNodalContactForce( rRightHandSideVector, rCurrentProcessInfo, CalculateResidualVectorFlag );
        mReactionIsAdded = true;
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void MPMParticlePenaltyCouplingInterfaceCondition::CalculateNodalContactForce( const VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, const bool CalculateResidualVectorFlag )
{
    if ( CalculateResidualVectorFlag == true )
    {
        GeometryType& rGeom = GetGeometry();
        const unsigned int number_of_nodes = rGeom.size();
        const unsigned int dimension = rGeom.WorkingSpaceDimension();
        const unsigned int block_size = this->GetBlockSize();

        // Calculate nodal forces
        Vector nodal_force = ZeroVector(3);
        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            for (unsigned int j = 0; j < dimension; j++)
            {
                nodal_force[j] = rRightHandSideVector[block_size * i + j];
            }

            // Check whether there nodes are active and associated to material point elements
            const double& nodal_mass = rGeom[i].FastGetSolutionStepValue(NODAL_MASS, 0);
            if (nodal_mass > std::numeric_limits<double>::epsilon())
            {
                rGeom[i].SetLock();
                rGeom[i].FastGetSolutionStepValue(REACTION) += nodal_force;
                rGeom[i].UnSetLock();
            }

        }
    }
}

//************************************************************************************
//************************************************************************************

void MPMParticlePenaltyCouplingInterfaceCondition::CalculateInterfaceContactForce( ProcessInfo& rCurrentProcessInfo )
{
    GeometryType& rGeom = GetGeometry();
    const unsigned int number_of_nodes = rGeom.PointsNumber();
    // const unsigned int dimension = rGeom.WorkingSpaceDimension();

    // Prepare variables
    GeneralVariables Variables;
    const array_1d<double, 3 > & xg_c = this->GetValue(MPC_COORD);
    const double & MPC_Area = this->GetValue(MPC_AREA);
    Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg_c);

    // Interpolate the force to MPC_Force assuming linear shape function
    array_1d<double, 3 > MPC_Force = ZeroVector(3);
    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        // Check whether there is material point inside the node
        const double& nodal_mass = rGeom[i].FastGetSolutionStepValue(NODAL_MASS, 0);
        const double nodal_area  = rGeom[i].FastGetSolutionStepValue(NODAL_AREA, 0);
        const Vector nodal_force = rGeom[i].FastGetSolutionStepValue(REACTION);

        if (nodal_mass > std::numeric_limits<double>::epsilon())
        {
            MPC_Force += Variables.N[i] * nodal_force * MPC_Area / nodal_area;
        }
    }

    // Apply in the normal contact direction and allow releasing motion
    if (Is(CONTACT))
    {
        // Apply only in the normal direction
        array_1d<double, 3 > & unit_normal_vector = this->GetValue(MPC_NORMAL);
        ParticleMechanicsMathUtilities<double>::Normalize(unit_normal_vector);
        const double normal_force = MathUtils<double>::Dot(MPC_Force,unit_normal_vector);

        // This check is done to avoid sticking forces
        if (normal_force > 0.0)
            MPC_Force = -1.0 * normal_force * unit_normal_vector;
        else
            MPC_Force = ZeroVector(3);
    }
    // Apply a sticking contact
    else{
        MPC_Force *= -1.0;
    }

    // Set Contact Force
    this->SetValue(MPC_CONTACT_FORCE, MPC_Force);

}


} // Namespace Kratos


