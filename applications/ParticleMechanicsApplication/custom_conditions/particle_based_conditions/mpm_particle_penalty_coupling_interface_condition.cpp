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

// void MPMParticlePenaltyCouplingInterfaceCondition::CalculateNodalContactForceAlternative( const MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo, const bool CalculateStiffnessMatrixFlag )
// {
//     if ( CalculateStiffnessMatrixFlag == true )
//     {
//         GeometryType& rGeom = GetGeometry();
//         const unsigned int number_of_nodes = rGeom.size();
//         const unsigned int dimension = rGeom.WorkingSpaceDimension();
//         const unsigned int block_size = this->GetBlockSize();
//         const unsigned int matrix_size = number_of_nodes * block_size;

//         // Prepare variables
//         GeneralVariables Variables;
//         Variables.CurrentDisp = this->CalculateCurrentDisp(Variables.CurrentDisp, rCurrentProcessInfo);

//         // Obtain nodal_displacement
//         Vector nodal_displacement = ZeroVector(matrix_size);
//         for (unsigned int i = 0; i < number_of_nodes; i++)
//         {
//             for ( unsigned int j = 0; j < dimension; j++)
//             {
//                 nodal_displacement[block_size * i + j] = Variables.CurrentDisp(i,j);
//             }
//         }

//         // Calculate local reaction
//         Vector local_reaction = ZeroVector(matrix_size);
//         local_reaction = prod(rLeftHandSideMatrix, nodal_displacement);

//         // if((this->GetId() == 8554 ||this->GetId() == 8560 )){
//         //     std::cout << "=================================================" << std::endl;
//         //     std::cout << "rLeftHandSideMatrix: " << rLeftHandSideMatrix << std::endl;
//         //     std::cout << "nodal_displacement: " << nodal_displacement << std::endl;
//         //     std::cout << "local_reaction: " << local_reaction << std::endl;
//         // }

//         // Calculate nodal forces
//         Vector nodal_force = ZeroVector(3);
//         for (unsigned int i = 0; i < number_of_nodes; i++)
//         {
//             for (unsigned int j = 0; j < dimension; j++)
//             {
//                 nodal_force[j] = local_reaction[block_size * i + j];
//             }

//             // Check whether there nodes are active and associated to material point elements
//             const double& nodal_mass = rGeom[i].FastGetSolutionStepValue(NODAL_MASS, 0);

//             // if((this->GetId() == 8554 ||this->GetId() == 8560 )){
//             //     std::cout << "NODE: " << i << std::endl;
//             //     std::cout << "nodal_mass: " << nodal_mass << std::endl;
//             //     std::cout << "nodal_force: " << nodal_force << std::endl;
//             //     std::cout << "=================================================" << std::endl;
//             // }
//             if (nodal_mass > std::numeric_limits<double>::epsilon())
//             {
//                 rGeom[i].SetLock();
//                 rGeom[i].FastGetSolutionStepValue(REACTION) += nodal_force;
//                 rGeom[i].UnSetLock();
//             }

//         }
//     }
// }

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

    if (Is(CONTACT))
    {
        // Apply only in the normal direction
        array_1d<double, 3 > & unit_normal_vector = this->GetValue(MPC_NORMAL);
        ParticleMechanicsMathUtilities<double>::Normalize(unit_normal_vector);
        const double normal_force = MathUtils<double>::Dot(MPC_Force,unit_normal_vector);

        // if((this->GetId() == 8554 ||this->GetId() == 8560 )){
        //     std::cout << "=================================================" << std::endl;
        //     std::cout << "CONDITION ID: " << this->GetId() << std::endl;
        //     std::cout << "MPC_Imposed_Displacement" << this->GetValue(MPC_IMPOSED_DISPLACEMENT) << std::endl;
        //     std::cout << "MPC_Force: " << MPC_Force << std::endl;
        //     std::cout << "unit_normal_vector: " << unit_normal_vector << std::endl;
        //     std::cout << "normal_force: " << normal_force << std::endl;
        //     std::cout << "=================================================" << std::endl;
        // }

        // This check is done to avoid sticking forces
        if (normal_force > 0.0)
            MPC_Force = -1.0 * normal_force * unit_normal_vector;
        else
            MPC_Force = ZeroVector(3);
    }

    // Set Contact Force
    this->SetValue(MPC_CONTACT_FORCE, MPC_Force);

    // 2. COMPUTING FORCE FROM MOMENTUM
    // const double& delta_time = rCurrentProcessInfo[DELTA_TIME];
    // const array_1d<double, 3 > & prescribed_vel = this->GetValue(MPC_VELOCITY);

    // Initiated equivalent mass and kinematic variables
    // double equivalent_MPC_mass = 0.0;
    // // array_1d<double, 3 > new_MPC_velocity = ZeroVector(3);
    // array_1d<double, 3 > old_MPC_velocity = ZeroVector(3);
    // // array_1d<double, 3 > new_MPC_acceleration = ZeroVector(3);
    // // array_1d<double, 3 > old_MPC_acceleration = ZeroVector(3);

    // // Check whether background mesh is active: all background nodes should contain nodal_mass
    // bool mesh_active = true;
    // for (unsigned int i = 0; i < number_of_nodes; i++)
    // {
    //     // Check whether there is material point inside the node
    //     const double& nodal_mass = rGeom[i].FastGetSolutionStepValue(NODAL_MASS, 0);
    //     if (nodal_mass < std::numeric_limits<double>::epsilon())
    //         mesh_active = false;
    // }

    // Only apply force when mesh is active
    // if (mesh_active)
    // {
    //     // Interpolate equivalent mass and kinematic variables from connectivity
    //     for (unsigned int i = 0; i < number_of_nodes; i++)
    //     {
    //         const double& nodal_mass = rGeom[i].FastGetSolutionStepValue(NODAL_MASS, 0);
    //         const array_1d<double, 3 > & nodal_prev_velocity     = rGeom[i].FastGetSolutionStepValue(VELOCITY, 1);
    //         // const array_1d<double, 3 > & nodal_curr_acceleration = rGeom[i].FastGetSolutionStepValue(ACCELERATION, 0);
    //         // const array_1d<double, 3 > & nodal_prev_acceleration = rGeom[i].FastGetSolutionStepValue(ACCELERATION, 1);

    //         old_MPC_velocity += Variables.N[i] * nodal_prev_velocity;
    //         // new_MPC_acceleration += Variables.N[i] * nodal_curr_acceleration;
    //         // old_MPC_acceleration += Variables.N[i] * nodal_prev_acceleration;

    //         equivalent_MPC_mass += Variables.N[i] * Variables.N[i] / nodal_mass;
    //     }

    //     // Obtain MPC_mass and new_MPC_velocity
    //     const double MPC_mass = 1.0 / equivalent_MPC_mass;
    //     // new_MPC_velocity = old_MPC_velocity + 0.5 * delta_time * (new_MPC_acceleration + old_MPC_acceleration);

    //     // Compute force
    //     array_1d<double, 3 > MPC_Force = ZeroVector(3);
    //     MPC_Force = MPC_mass / delta_time * (old_MPC_velocity - prescribed_vel);

    //     if(this->GetId() == 1){
    //         std::cout << "=================================================" << std::endl;
    //         std::cout << "CONDITION ID: " << this->GetId() << std::endl;
    //         std::cout << "MPC_Force_1: " << MPC_Force << std::endl;
    //         std::cout << "MPC_mass: " << MPC_mass << std::endl;
    //         // std::cout << "new_MPC_velocity: " << new_MPC_velocity << std::endl;
    //         // std::cout << "new_MPC_acceleration: " << new_MPC_acceleration << std::endl;
    //         std::cout << "old_MPC_velocity: " << old_MPC_velocity << std::endl;
    //         // std::cout << "old_MPC_acceleration: " << old_MPC_acceleration << std::endl;
    //         std::cout << "prescribed_vel: " << prescribed_vel << std::endl;
    //     }

    //     // Apply only in the normal direction
    //     array_1d<double, 3 > & unit_normal_vector = this->GetValue(MPC_NORMAL);
    //     ParticleMechanicsMathUtilities<double>::Normalize(unit_normal_vector);
    //     const double normal_force = MathUtils<double>::Dot(MPC_Force,unit_normal_vector);
    //     if(this->GetId() == 1)
    //         std::cout << "normal_force: " << normal_force << std::endl;

    //     // Check velocity normal_velocity
    //     const double normal_velocity = MathUtils<double>::Dot(old_MPC_velocity,unit_normal_vector);
    //     if(this->GetId() == 1)
    //         std::cout << "normal_velocity: " << normal_velocity << std::endl;

    //     // This check is done to avoid sticking forces
    //     if (normal_force < 0.0)
    //         MPC_Force = normal_force * unit_normal_vector;
    //     else
    //         MPC_Force = ZeroVector(3);

    //     // Set Contact Force
    //     this->SetValue(MPC_CONTACT_FORCE, MPC_Force);

    //     if(this->GetId() == 1){
    //         std::cout << "MPC_FORCE_NORMAL: " << MPC_Force << std::endl;
    //         std::cout << "=================================================" << std::endl;
    //     }
    // }
}


} // Namespace Kratos


