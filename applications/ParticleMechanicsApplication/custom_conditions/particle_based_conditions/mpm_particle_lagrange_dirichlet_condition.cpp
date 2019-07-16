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

    // Additional treatment for slip conditions
    if (Is(SLIP))
    {
        GeometryType& r_geometry = GetGeometry();
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

    // Resizing as needed the LHS
    const unsigned int matrix_size = number_of_nodes * block_size * 2;

    //if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    //{
    if ( rLeftHandSideMatrix.size1() != matrix_size )
    {
        rLeftHandSideMatrix.resize( matrix_size, matrix_size, false );
    }

    noalias( rLeftHandSideMatrix ) = ZeroMatrix(matrix_size,matrix_size); //resetting LHS
    //}

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
    const double augmention_factor = this->GetValue(SCALAR_LAGRANGE_MULTIPLIER);

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
        for (unsigned int i = 0; i < number_of_nodes; i++) //loop over Lagrange Multipliers
        {
            for (unsigned int j = 0; j < number_of_nodes; j++) // lopp over shape functions of displacements
            {
                double NN = Variables.N[j] * Variables.N[i];

                const unsigned int ibase = i * dimension + dimension * number_of_nodes;
                const unsigned int jbase = j * dimension;

                // Matrix in following shape:
                // |0 H^T|
                // |H 0  |
                //lambda in X
                rLeftHandSideMatrix(ibase, jbase)         = NN;
                rLeftHandSideMatrix(ibase + 1, jbase + 1) = NN;
                if (dimension == 3){
                    rLeftHandSideMatrix(ibase + 2, jbase + 2) = NN;
                }


                rLeftHandSideMatrix(jbase, ibase)         = NN;
                rLeftHandSideMatrix(jbase + 1, ibase + 1) = NN;
                if (dimension == 3){
                    rLeftHandSideMatrix(jbase + 2, ibase + 2) = NN;
                }

                if (augmention_factor > 0.0)
                {
                    //Depending on the solver if needed.
                    if (rLeftHandSideMatrix(dimension * j, ibase) <= 0.000001)
                        rLeftHandSideMatrix(ibase, ibase) = augmention_factor;
                    if (rLeftHandSideMatrix(dimension * j + 1, ibase + 1) <= 0.000001)
                        rLeftHandSideMatrix(ibase + 1, ibase + 1) = augmention_factor;
                    if (dimension == 3){
                        if (rLeftHandSideMatrix(dimension * j + 2, ibase + 2) <= 0.000001)
                            rLeftHandSideMatrix(ibase +2, ibase +2) = augmention_factor;
                    }
                }
            }
        }

        // Calculate LHS Matrix and RHS Vector
        //if ( CalculateStiffnessMatrixFlag == true )
        //{
        rLeftHandSideMatrix  *= this->GetIntegrationWeight();
        //}


        if ( CalculateResidualVectorFlag == true )
        {
            Vector u(matrix_size);
            for (unsigned int i = 0; i < number_of_nodes; i++)
            {
                const array_1d<double, 3> disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
                int index = dimension * i;
                u[index]     = disp[0] - imposed_displacement[0];
                u[index + 1] = disp[1] - imposed_displacement[1];
                if (dimension == 3){
                    u[index + 2] = disp[2] - imposed_displacement[2];
                }
            }
            for (unsigned int i = 0; i < number_of_nodes; i++)
            {
                const array_1d<double, 3> lagrange_multiplier = GetGeometry()[i].FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);
                int index = dimension * i + dimension * number_of_nodes;
                u[index]     = lagrange_multiplier[0];
                u[index + 1] = lagrange_multiplier[1];
                if (dimension == 3){
                    u[index + 2] = lagrange_multiplier[2];
                }

            }

            noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, u);
        }


    }
    else{
        /*// To improve stability: use identity matrix to avoid nonzero diagonal LHS matrix
        if ( CalculateStiffnessMatrixFlag == true )
        {
            noalias(rLeftHandSideMatrix) = IdentityMatrix(matrix_size);
        }*/
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void MPMParticleLagrangeDirichletCondition::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    MPMParticleBaseDirichletCondition::FinalizeSolutionStep(rCurrentProcessInfo);

    // Additional treatment for slip conditions
    if (Is(SLIP))
    {
        GeometryType& r_geometry = GetGeometry();
        const unsigned int number_of_nodes = r_geometry.PointsNumber();

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

int MPMParticleLagrangeDirichletCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
    MPMParticleBaseDirichletCondition::Check(rCurrentProcessInfo);

    // Verify that the dofs exist
    for (const auto& r_node : this->GetGeometry().Points())
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NORMAL,r_node)

    return 0;
}

void MPMParticleLagrangeDirichletCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    GeometryType& rGeom = GetGeometry();
    const unsigned int number_of_nodes = rGeom.size();
    const unsigned int dim = rGeom.WorkingSpaceDimension();
    if (rResult.size() != dim * number_of_nodes * 2)
    {
        rResult.resize(dim * number_of_nodes * 2 ,false);
    }


    if(dim == 2)
    {
        for (unsigned int i = 0; i < number_of_nodes; ++i)
        {
            unsigned int index = i * 2;
            rResult[index    ] = rGeom[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = rGeom[i].GetDof(DISPLACEMENT_Y).EquationId();

            index = i * 2 + number_of_nodes * dim;
            rResult[index    ] = rGeom[i].GetDof(VECTOR_LAGRANGE_MULTIPLIER_X).EquationId();
            rResult[index + 1] = rGeom[i].GetDof(VECTOR_LAGRANGE_MULTIPLIER_Y).EquationId();
        }
    }
    else
    {
        for (unsigned int i = 0; i < number_of_nodes; ++i)
        {
            unsigned int index = i * 3;
            rResult[index    ] = rGeom[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = rGeom[i].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = rGeom[i].GetDof(DISPLACEMENT_Z).EquationId();

            index = i * 3 + number_of_nodes * dim;
            rResult[index    ] = rGeom[i].GetDof(VECTOR_LAGRANGE_MULTIPLIER_X).EquationId();
            rResult[index + 1] = rGeom[i].GetDof(VECTOR_LAGRANGE_MULTIPLIER_Y).EquationId();
            rResult[index + 2] = rGeom[i].GetDof(VECTOR_LAGRANGE_MULTIPLIER_Z).EquationId();
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

    GeometryType& rGeom = GetGeometry();
    const unsigned int number_of_nodes = rGeom.size();
    const unsigned int dim =  rGeom.WorkingSpaceDimension();
    rElementalDofList.resize(0);
    rElementalDofList.reserve(dim * number_of_nodes *2);

    GeneralVariables Variables;
    const array_1d<double,3> & xg_c = this->GetValue(MPC_COORD);
    // Calculating shape function
    Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg_c);
    array_1d<bool,3> zero_shape_function;
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        zero_shape_function[i]=false;
    }
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        if (Variables.N[i] <= 0.00001 && rGeom[i].FastGetSolutionStepValue(NODAL_MASS, 0) <= 0.00001)
        {
            zero_shape_function[i]=true;
        }
    }
    if(dim == 2)
    {
        for (unsigned int i = 0; i < number_of_nodes; ++i)
        {
            if (zero_shape_function[i]==true){
                rGeom[i].pGetDof(DISPLACEMENT_X)->FixDof();
                rGeom[i].pGetDof(DISPLACEMENT_Y)->FixDof();
            }
            rElementalDofList.push_back( rGeom[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back( rGeom[i].pGetDof(DISPLACEMENT_Y));

        }
        for (unsigned int i = 0; i < number_of_nodes; ++i)
        {
            if (zero_shape_function[i]==true){
                rGeom[i].pGetDof(VECTOR_LAGRANGE_MULTIPLIER_X)->FixDof();
                rGeom[i].pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Y)->FixDof();
            }
            rElementalDofList.push_back( rGeom[i].pGetDof(VECTOR_LAGRANGE_MULTIPLIER_X));
            rElementalDofList.push_back( rGeom[i].pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Y));
        }
    }
    else
    {
        for (unsigned int i = 0; i < number_of_nodes; ++i)
        {
            if (zero_shape_function[i]==true){
                rGeom[i].pGetDof(DISPLACEMENT_X)->FixDof();
                rGeom[i].pGetDof(DISPLACEMENT_Y)->FixDof();
                rGeom[i].pGetDof(DISPLACEMENT_Z)->FixDof();
            }
            rElementalDofList.push_back( rGeom[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back( rGeom[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back( rGeom[i].pGetDof(DISPLACEMENT_Z));
        }
        for (unsigned int i = 0; i < number_of_nodes; ++i)
        {
            if (zero_shape_function[i]==true){
                rGeom[i].pGetDof(VECTOR_LAGRANGE_MULTIPLIER_X)->FixDof();
                rGeom[i].pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Y)->FixDof();
                rGeom[i].pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Z)->FixDof();
            }
            rElementalDofList.push_back( rGeom[i].pGetDof(VECTOR_LAGRANGE_MULTIPLIER_X));
            rElementalDofList.push_back( rGeom[i].pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Y));
            rElementalDofList.push_back( rGeom[i].pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Z));
        }
    }

    KRATOS_CATCH("")
}

void MPMParticleLagrangeDirichletCondition::GetValuesVector(
    Vector& rValues,
    int Step
    )
{
    GeometryType& rGeom = GetGeometry();
    const unsigned int number_of_nodes = rGeom.size();
    const unsigned int dim = rGeom.WorkingSpaceDimension();
    const unsigned int mat_size = number_of_nodes * dim * 2;

    if (rValues.size() != mat_size)
    {
        rValues.resize(mat_size, false);
    }
    rValues = ZeroVector(mat_size);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        const array_1d<double, 3 > & Displacement = rGeom[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
        unsigned int index = i * dim;
        for(unsigned int k = 0; k < dim; ++k)
        {
            rValues[index + k] = Displacement[k];
        }
    }
}

void MPMParticleLagrangeDirichletCondition::GetFirstDerivativesVector(
    Vector& rValues,
    int Step
    )
{
    GeometryType& rGeom = GetGeometry();
    const unsigned int number_of_nodes = rGeom.size();
    const unsigned int dim = rGeom.WorkingSpaceDimension();
    const unsigned int mat_size = number_of_nodes * dim * 2;

    if (rValues.size() != mat_size)
    {
        rValues.resize(mat_size, false);
    }
    rValues = ZeroVector(mat_size);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        const array_1d<double, 3 > & Velocity = rGeom[i].FastGetSolutionStepValue(VELOCITY, Step);
        const unsigned int index = i * dim;
        for(unsigned int k = 0; k < dim; ++k)
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
    GeometryType& rGeom = GetGeometry();
    const unsigned int number_of_nodes = rGeom.size();
    const unsigned int dim = rGeom.WorkingSpaceDimension();
    const unsigned int mat_size = number_of_nodes * dim * 2;

    if (rValues.size() != mat_size)
    {
        rValues.resize(mat_size, false);
    }
    rValues = ZeroVector(mat_size);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        const array_1d<double, 3 > & Acceleration = rGeom[i].FastGetSolutionStepValue(ACCELERATION, Step);
        const unsigned int index = i * dim;
        for(unsigned int k = 0; k < dim; ++k)
        {
            rValues[index + k] = Acceleration[k];
        }
    }
}




} // Namespace Kratos


