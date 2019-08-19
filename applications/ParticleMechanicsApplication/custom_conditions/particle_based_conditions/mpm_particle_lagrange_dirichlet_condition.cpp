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
#include "custom_conditions/particle_based_conditions/mpm_particle_lagrange_dirichlet_condition.h"
#include "includes/kratos_flags.h"
#include "utilities/math_utils.h"
#include "custom_utilities/particle_mechanics_math_utilities.h"
#include "includes/checks.h"

namespace Kratos
{
//******************************* CONSTRUCTOR ****************************************
//************************************************************************************

MPMParticleLagrangeDirichletCondition::MPMParticleLagrangeDirichletCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : MPMParticleBaseDirichletCondition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

MPMParticleLagrangeDirichletCondition::MPMParticleLagrangeDirichletCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : MPMParticleBaseDirichletCondition( NewId, pGeometry, pProperties )
{
}

//********************************* CREATE *******************************************
//************************************************************************************

Condition::Pointer MPMParticleLagrangeDirichletCondition::Create(IndexType NewId,GeometryType::Pointer pGeometry,PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<MPMParticleLagrangeDirichletCondition>(NewId, pGeometry, pProperties);
}

//************************************************************************************
//************************************************************************************

Condition::Pointer MPMParticleLagrangeDirichletCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<MPMParticleLagrangeDirichletCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

//******************************* DESTRUCTOR *****************************************
//************************************************************************************

MPMParticleLagrangeDirichletCondition::~MPMParticleLagrangeDirichletCondition()
{
}

//************************************************************************************
//************************************************************************************

void MPMParticleLagrangeDirichletCondition::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    MPMParticleBaseDirichletCondition::InitializeSolutionStep( rCurrentProcessInfo );

    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        array_1d<double, 3 > & lagrange_multiplier  = r_geometry[i].FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            lagrange_multiplier[j] = 0;
        }

        r_geometry[i].SetValue(VECTOR_LAGRANGE_MULTIPLIER,lagrange_multiplier);
    }



}

void MPMParticleLagrangeDirichletCondition::CalculateAll(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo,
    bool CalculateStiffnessMatrixFlag,
    bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int block_size = this->GetBlockSize();
    const unsigned int matrix_size = number_of_nodes * dimension * 2;

    // Resizing as needed the LHS
    if ( CalculateStiffnessMatrixFlag == true )
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

    // Get imposed displacement and normal vector
    const array_1d<double, 3 > & xg_c = this->GetValue(MPC_COORD);
    const array_1d<double, 3 > & imposed_displacement = this->GetValue (MPC_IMPOSED_DISPLACEMENT);

    // Prepare variables
    GeneralVariables Variables;
    const double augmentation_factor = this->GetValue(SCALAR_LAGRANGE_MULTIPLIER);

    // Calculating shape function
    Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg_c);
    Variables.CurrentDisp = this->CalculateCurrentDisp(Variables.CurrentDisp, rCurrentProcessInfo);

    // Check contact: Check contact penetration: if <0 apply constraint, otherwise no
    bool apply_constraints = true;
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
        if (penetration >= 0.0)
        {
            apply_constraints = false;
        }
    }

    if (apply_constraints)
    {
        Matrix lagrange_matrix = ZeroMatrix(matrix_size, matrix_size);

         // Loop over Lagrange Multipliers
        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            // Loop over shape functions of displacements
            for (unsigned int j = 0; j < number_of_nodes; j++)
            {
                const double NN = Variables.N[j] * Variables.N[i];
                const unsigned int ibase = i * dimension + dimension * number_of_nodes;
                const unsigned int jbase = j * dimension;

                // Matrix in following shape:
                // |0 H^T|
                // |H 0  |
                // Lambda in X
                for (unsigned int k = 0; k < dimension; k++)
                {
                    lagrange_matrix(ibase+k, jbase+k) = NN;
                    lagrange_matrix(jbase+k, ibase+k) = NN;
                }



                if (augmentation_factor > 0.0)
                {
                    //Stabilization
                    //if (lagrange_matrix(dimension * j, ibase) <= 0.01)
                        lagrange_matrix(ibase, ibase) = -augmentation_factor;

                    //if (lagrange_matrix(dimension * j + 1, ibase + 1) <= 0.01)
                        lagrange_matrix(ibase + 1, ibase + 1) = -augmentation_factor;

                    if (dimension == 3){
                        if (lagrange_matrix(dimension * j + 2, ibase + 2) <= 0.001)
                            lagrange_matrix(ibase +2, ibase +2) = -augmentation_factor;
                    }
                }
            }
        }

        lagrange_matrix  *= this->GetIntegrationWeight();

        // Calculate LHS Matrix and RHS Vector
        if ( CalculateStiffnessMatrixFlag == true )
        {
            rLeftHandSideMatrix = lagrange_matrix;
        }

        if ( CalculateResidualVectorFlag == true )
        {
            Vector gap_function = ZeroVector(matrix_size);
            for (unsigned int i = 0; i < number_of_nodes; i++)
            {
                const array_1d<double, 3> disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
                const int index = dimension * i;

                const array_1d<double, 3> lagrange_multiplier = GetGeometry()[i].FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);
                const int lagrange_index = dimension * i + dimension * number_of_nodes;

                for (unsigned int j = 0; j < dimension; j++){
                    gap_function[index+j]          = disp[j] - imposed_displacement[j];
                    gap_function[lagrange_index+j] = lagrange_multiplier[j];
                }
            }

            noalias(rRightHandSideVector) -= prod(lagrange_matrix, gap_function);
        }

    }
    else{
        // To improve stability: use identity matrix to avoid nonzero diagonal LHS matrix
        if ( CalculateStiffnessMatrixFlag == true )
        {
            noalias(rLeftHandSideMatrix) = IdentityMatrix(matrix_size);
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

int MPMParticleLagrangeDirichletCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
    MPMParticleBaseDirichletCondition::Check(rCurrentProcessInfo);
    const unsigned int number_of_nodes = GetGeometry().size();
    // Verify that the dofs exist
    for ( IndexType i = 0; i < number_of_nodes; i++ ) {
        const NodeType &rnode = this->GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VECTOR_LAGRANGE_MULTIPLIER,rnode)

        KRATOS_CHECK_DOF_IN_NODE(VECTOR_LAGRANGE_MULTIPLIER_X,rnode)
        KRATOS_CHECK_DOF_IN_NODE(VECTOR_LAGRANGE_MULTIPLIER_Y,rnode)
        KRATOS_CHECK_DOF_IN_NODE(VECTOR_LAGRANGE_MULTIPLIER_Z,rnode)
    }

    return 0;
}

void MPMParticleLagrangeDirichletCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    if (rResult.size() != dimension * number_of_nodes * 2)
    {
        rResult.resize(dimension * number_of_nodes * 2 ,false);
    }


    if(dimension == 2)
    {
        for (unsigned int i = 0; i < number_of_nodes; ++i)
        {
            unsigned int index = i * dimension;
            rResult[index    ] = r_geometry[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_geometry[i].GetDof(DISPLACEMENT_Y).EquationId();


            index = i * dimension + number_of_nodes * dimension;
            rResult[index    ] = r_geometry[i].GetDof(VECTOR_LAGRANGE_MULTIPLIER_X).EquationId();
            rResult[index + 1] = r_geometry[i].GetDof(VECTOR_LAGRANGE_MULTIPLIER_Y).EquationId();
        }
    }
    else
    {
        for (unsigned int i = 0; i < number_of_nodes; ++i)
        {
            unsigned int index = i * dimension;
            rResult[index    ] = r_geometry[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_geometry[i].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = r_geometry[i].GetDof(DISPLACEMENT_Z).EquationId();

            index = i * dimension + number_of_nodes * dimension;
            rResult[index    ] = r_geometry[i].GetDof(VECTOR_LAGRANGE_MULTIPLIER_X).EquationId();
            rResult[index + 1] = r_geometry[i].GetDof(VECTOR_LAGRANGE_MULTIPLIER_Y).EquationId();
            rResult[index + 2] = r_geometry[i].GetDof(VECTOR_LAGRANGE_MULTIPLIER_Z).EquationId();
        }
    }


    KRATOS_CATCH("")
}

void MPMParticleLagrangeDirichletCondition::GetDofList(
    DofsVectorType& rElementalDofList,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension =  r_geometry.WorkingSpaceDimension();
    rElementalDofList.resize(0);
    rElementalDofList.reserve(dimension * number_of_nodes *2);

    GeneralVariables Variables;
    const array_1d<double,3> & xg_c = this->GetValue(MPC_COORD);
    // Calculating shape function
    Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg_c);
    if(dimension == 2)
    {
        for (unsigned int i = 0; i < number_of_nodes; ++i)
        {
            rElementalDofList.push_back( r_geometry[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back( r_geometry[i].pGetDof(DISPLACEMENT_Y));

        }
        for (unsigned int i = 0; i < number_of_nodes; ++i)
        {
            rElementalDofList.push_back( r_geometry[i].pGetDof(VECTOR_LAGRANGE_MULTIPLIER_X));
            rElementalDofList.push_back( r_geometry[i].pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Y));
        }
    }
    else
    {
        for (unsigned int i = 0; i < number_of_nodes; ++i)
        {
            rElementalDofList.push_back( r_geometry[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back( r_geometry[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back( r_geometry[i].pGetDof(DISPLACEMENT_Z));
        }
        for (unsigned int i = 0; i < number_of_nodes; ++i)
        {
            rElementalDofList.push_back( r_geometry[i].pGetDof(VECTOR_LAGRANGE_MULTIPLIER_X));
            rElementalDofList.push_back( r_geometry[i].pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Y));
            rElementalDofList.push_back( r_geometry[i].pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Z));
        }
    }

    KRATOS_CATCH("")
}

void MPMParticleLagrangeDirichletCondition::GetValuesVector(
    Vector& rValues,
    int Step
    )
{
    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    const unsigned int mat_size = number_of_nodes * dimension * 2;

    if (rValues.size() != mat_size)
    {
        rValues.resize(mat_size, false);
    }
    rValues = ZeroVector(mat_size);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        const array_1d<double, 3 > & Displacement = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
        const array_1d<double, 3 > & lagrange_multiplier = r_geometry[i].FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER, Step);
        unsigned int index = i * dimension;
        const unsigned int lagrange_index = number_of_nodes * dimension;
        for(unsigned int k = 0; k < dimension; ++k)
        {
            rValues[index + k] = Displacement[k];
            rValues[lagrange_index + k] = lagrange_multiplier[k];
        }
    }
}

void MPMParticleLagrangeDirichletCondition::GetFirstDerivativesVector(
    Vector& rValues,
    int Step
    )
{
    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    const unsigned int mat_size = number_of_nodes * dimension * 2;

    if (rValues.size() != mat_size)
    {
        rValues.resize(mat_size, false);
    }
    rValues = ZeroVector(mat_size);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        const array_1d<double, 3 > & Velocity = r_geometry[i].FastGetSolutionStepValue(VELOCITY, Step);
        const unsigned int index = i * dimension;
        for(unsigned int k = 0; k < dimension; ++k)
        {
            rValues[index + k] = Velocity[k];
        }
    }
}

void MPMParticleLagrangeDirichletCondition::GetSecondDerivativesVector(
    Vector& rValues,
    int Step
    )
{
    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    const unsigned int mat_size = number_of_nodes * dimension * 2;

    if (rValues.size() != mat_size)
    {
        rValues.resize(mat_size, false);
    }
    rValues = ZeroVector(mat_size);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        const array_1d<double, 3 > & Acceleration = r_geometry[i].FastGetSolutionStepValue(ACCELERATION, Step);
        const unsigned int index = i * dimension;
        for(unsigned int k = 0; k < dimension; ++k)
        {
            rValues[index + k] = Acceleration[k];
        }
    }
}




} // Namespace Kratos


