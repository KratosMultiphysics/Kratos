//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Veronika Singer
//


// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_conditions/particle_based_conditions/mpm_particle_fixdof_dirichlet_condition.h"
#include "includes/kratos_flags.h"
#include "utilities/math_utils.h"
#include "custom_utilities/particle_mechanics_math_utilities.h"
#include "includes/checks.h"

namespace Kratos
{
//******************************* CONSTRUCTOR ****************************************
//************************************************************************************

MPMParticleFixDofDirichletCondition::MPMParticleFixDofDirichletCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : MPMParticleBaseDirichletCondition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

MPMParticleFixDofDirichletCondition::MPMParticleFixDofDirichletCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : MPMParticleBaseDirichletCondition( NewId, pGeometry, pProperties )
{
}

//********************************* CREATE *******************************************
//************************************************************************************

Condition::Pointer MPMParticleFixDofDirichletCondition::Create(IndexType NewId,GeometryType::Pointer pGeometry,PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<MPMParticleFixDofDirichletCondition>(NewId, pGeometry, pProperties);
}

//************************************************************************************
//************************************************************************************

Condition::Pointer MPMParticleFixDofDirichletCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<MPMParticleFixDofDirichletCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

//******************************* DESTRUCTOR *****************************************
//************************************************************************************

MPMParticleFixDofDirichletCondition::~MPMParticleFixDofDirichletCondition()
{
}

//************************************************************************************
//************************************************************************************

void MPMParticleFixDofDirichletCondition::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    MPMParticleBaseDirichletCondition::InitializeSolutionStep( rCurrentProcessInfo );

    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if (this->GetValue(FIX_DOF))
    {
        for (unsigned int i = 0; i < number_of_nodes; ++i)
        {
            r_geometry[i].SetLock();
            r_geometry[i].pGetDof(DISPLACEMENT_X)->FixDof();
            r_geometry[i].pGetDof(DISPLACEMENT_Y)->FixDof();
            if(dimension == 3)
                r_geometry[i].pGetDof(DISPLACEMENT_Z)->FixDof();
            r_geometry[i].UnSetLock();
        }

    }

    // Additional treatment for slip conditions
    if (Is(SLIP))
    {
        //GeometryType& r_geometry = GetGeometry();
        const unsigned int number_of_nodes = r_geometry.PointsNumber();
        const array_1d<double,3> & xg_c = this->GetValue(MPC_COORD);
        GeneralVariables Variables;

        // Calculating shape function
        Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg_c);

        // Normal Vector
        const array_1d<double,3> & unit_normal_vector = this->GetValue(MPC_NORMAL);

        // Here MPC contribution of normal vector are added
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            r_geometry[i].SetLock();
            r_geometry[i].Set(SLIP);
            r_geometry[i].FastGetSolutionStepValue(IS_STRUCTURE) = 2.0;
            r_geometry[i].FastGetSolutionStepValue(NORMAL) += Variables.N[i] * unit_normal_vector;
            r_geometry[i].UnSetLock();
        }
    }
}

//************************************************************************************
//************************************************************************************

void MPMParticleFixDofDirichletCondition::CalculateAll(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo,
    bool CalculateStiffnessMatrixFlag,
    bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY
    GeneralVariables Variables;
    const array_1d<double, 3 > & xg_c = this->GetValue(MPC_COORD);
    const array_1d<double, 3 > & imposed_displacement = this->GetValue (MPC_IMPOSED_DISPLACEMENT);
    Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg_c);
    Variables.CurrentDisp = this->CalculateCurrentDisp(Variables.CurrentDisp, rCurrentProcessInfo);

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int block_size = this->GetBlockSize();

    // Resizing as needed the LHS
    const unsigned int matrix_size = number_of_nodes * block_size;

    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != matrix_size )
        {
            rLeftHandSideMatrix.resize( matrix_size, matrix_size, false );
        }

        noalias( rLeftHandSideMatrix ) = ZeroMatrix(matrix_size,matrix_size); //resetting LHS
    }

    // Resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size( ) != matrix_size )
        {
            rRightHandSideVector.resize( matrix_size, false );
        }

        noalias( rRightHandSideVector ) = ZeroVector( matrix_size ); //resetting RHS
    }


    // Check contact: Check contact penetration: if <0 apply constraint, otherwise no
    this->SetValue(FIX_DOF,true);
    if (Is(CONTACT))
    {
        // NOTE: the unit_normal_vector is assumed always pointing outside the boundary
        array_1d<double, 3 > & unit_normal_vector = this->GetValue(MPC_NORMAL);
        ParticleMechanicsMathUtilities<double>::Normalize(unit_normal_vector);
        array_1d<double, 3 > field_displacement = ZeroVector(3);
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            if (Variables.N[i] > std::numeric_limits<double>::epsilon() )
            {
                for ( unsigned int j = 0; j < dimension; j++)
                {
                    field_displacement[j] += Variables.N[i] * Variables.CurrentDisp(i,j);
                }
            }
        }

        const double penetration = MathUtils<double>::Dot((field_displacement - imposed_displacement), unit_normal_vector);

        // If penetrates, apply constraint, otherwise no
        if (penetration >= 0.0001)
        {
            this->SetValue(FIX_DOF,false);
        }

    }


    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void MPMParticleFixDofDirichletCondition::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    MPMParticleBaseDirichletCondition::FinalizeSolutionStep(rCurrentProcessInfo);

    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    for (unsigned int i = 0; i < number_of_nodes; ++i)
        {
            r_geometry[i].SetLock();
            r_geometry[i].pGetDof(DISPLACEMENT_X)->FreeDof();
            r_geometry[i].pGetDof(DISPLACEMENT_Y)->FreeDof();
            if(dimension == 3)
                r_geometry[i].pGetDof(DISPLACEMENT_Z)->FreeDof();
            r_geometry[i].UnSetLock();
        }

    // Additional treatment for slip conditions
    if (Is(SLIP))
    {
        // Here MPC normal vector and IS_STRUCTURE are reset
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            r_geometry[i].SetLock();
            r_geometry[i].Reset(SLIP);
            r_geometry[i].FastGetSolutionStepValue(IS_STRUCTURE) = 0.0;
            r_geometry[i].FastGetSolutionStepValue(NORMAL).clear();
            r_geometry[i].UnSetLock();
        }
    }

    KRATOS_CATCH( "" )
}

int MPMParticleFixDofDirichletCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
    MPMParticleBaseDirichletCondition::Check(rCurrentProcessInfo);

    // Verify that the dofs exist
    for (const auto& r_node : this->GetGeometry().Points())
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NORMAL,r_node)

    return 0;
}

} // Namespace Kratos


