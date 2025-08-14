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
#include "custom_utilities/mpm_math_utilities.h"

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

void MPMParticleLagrangeDirichletCondition::InitializeSolutionStep( const ProcessInfo& rCurrentProcessInfo )
{
    MPMParticleBaseDirichletCondition::InitializeSolutionStep( rCurrentProcessInfo );

    GeometryType& r_geometry = GetGeometry();

    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    
    // Prepare variables
    GeneralVariables Variables;
    MPMShapeFunctionPointValues(Variables.N);

    auto pLagrangeNode = r_geometry.GetGeometryParent(0).GetValue(MPC_LAGRANGE_NODE);
    array_1d<double, 3 > & r_lagrange_multiplier  = pLagrangeNode->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);

    for ( unsigned int j = 0; j < dimension; j++ )
    {
        #pragma omp atomic
        r_lagrange_multiplier[j] *= 0.0;
    }

    pLagrangeNode->SetLock();
    pLagrangeNode->FastGetSolutionStepValue(NODAL_AREA) += this->GetIntegrationWeight();
    pLagrangeNode->UnSetLock();

    // The calculation of the normal is required for slip and contact conditions
    if (Is(SLIP))
    {
        pLagrangeNode->SetLock();
        pLagrangeNode->Set(SLIP);
        pLagrangeNode->FastGetSolutionStepValue(NORMAL) += this->GetIntegrationWeight() * m_normal;
        pLagrangeNode->UnSetLock();
        

        // Here MPC contribution of normal vector are added
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            r_geometry[i].SetLock();
            r_geometry[i].Set(SLIP);
            r_geometry[i].SetValue(PARTICLE_BASED_SLIP, true);
            r_geometry[i].FastGetSolutionStepValue(NORMAL) += Variables.N[i] * m_normal * this->GetIntegrationWeight();
            r_geometry[i].UnSetLock();
        }
    }

}

//************************************************************************************
//************************************************************************************

void MPMParticleLagrangeDirichletCondition::CalculateAll(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    bool CalculateStiffnessMatrixFlag,
    bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY

    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int block_size = this->GetBlockSize();

    // Resizing as needed the LHS
    const unsigned int matrix_size = number_of_nodes * dimension + dimension;

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
    
    // Prepare variables
    GeneralVariables Variables;

    // Calculating shape function
    MPMShapeFunctionPointValues(Variables.N);

    bool apply_constraints = true;

    auto mp_counter = r_geometry.GetGeometryParent(0).GetValue(MP_COUNTER);
    if (mp_counter < 1 )
        this->Reset(ACTIVE);  

    if (Is(CONTACT) && this->Is(ACTIVE))
    {  
        auto pLagrangeNode = r_geometry.GetGeometryParent(0).GetValue(MPC_LAGRANGE_NODE);
        array_1d<double, 3>& r_lagrange_multiplier = pLagrangeNode->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER); 

        const double normal_force = MathUtils<double>::Dot(r_lagrange_multiplier, m_normal);
       
        if (normal_force > 0.0)
            apply_constraints=false;
     
    }
   
    if (apply_constraints && this->Is(ACTIVE))    
    {
        Matrix lagrange_matrix = ZeroMatrix(matrix_size, matrix_size);
        
        for (unsigned int i = 0; i < number_of_nodes; i++)
        {       
            const unsigned int ibase = dimension * number_of_nodes; 
            for (unsigned int k = 0; k < dimension; k++)
            {
                lagrange_matrix(i* dimension+k, ibase+k) = Variables.N[i];
                lagrange_matrix(ibase+k, i*dimension + k) = Variables.N[i];

                // perturbed Lagrangian
                if (m_penalty>0)
                    lagrange_matrix(ibase + k, ibase + k) = -1/m_penalty;

                    
            }
        }
         
        // Calculate LHS Matrix and RHS Vector
        if ( CalculateStiffnessMatrixFlag == true )
        {    
            rLeftHandSideMatrix = lagrange_matrix * this->GetIntegrationWeight();
        }

        if ( CalculateResidualVectorFlag == true )
        {
            //last rows of RHS
            Vector gap_function = ZeroVector(matrix_size);
            Vector right_hand_side = ZeroVector(matrix_size);
            for (unsigned int i = 0; i < number_of_nodes; i++)
            {
                    const array_1d<double, 3>& r_displacement = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT);
                    const int index = dimension * i;

                    for (unsigned int j = 0; j < dimension; j++)
                    {
                        gap_function[index+j] = r_displacement[j] ;
                    }
            }
            right_hand_side = prod(lagrange_matrix, gap_function);
            
            //first rows of RHS
            gap_function = ZeroVector(matrix_size);
            auto pLagrangeNode = r_geometry.GetGeometryParent(0).GetValue(MPC_LAGRANGE_NODE);
            const array_1d<double, 3>& r_lagrange_multiplier = pLagrangeNode->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);
           
            for (unsigned int j = 0; j < dimension; j++){
                gap_function[dimension * number_of_nodes+j] = r_lagrange_multiplier[j];
            }
            
            right_hand_side += prod(lagrange_matrix, gap_function);
        
            //add imposed displacement
            for (unsigned int j = 0; j < dimension; j++){
                    right_hand_side[dimension * number_of_nodes+j] -= m_imposed_displacement[j];
                }

            right_hand_side *= this->GetIntegrationWeight();
            noalias(rRightHandSideVector) = -right_hand_side;
        }
        
        if (Is(SLIP)){
            auto pLagrangeNode = r_geometry.GetGeometryParent(0).GetValue(MPC_LAGRANGE_NODE);
            // normalize normal at boundary particle node
            MPMMathUtilities<double>::Normalize(pLagrangeNode->FastGetSolutionStepValue(NORMAL));
            
            // rotate to normal-tangential frame
            if (CalculateStiffnessMatrixFlag == true){
                GetRotationTool().Rotate(rLeftHandSideMatrix, rRightHandSideVector, GetGeometry());
            } else {
                MatrixType temp = ZeroMatrix(matrix_size,matrix_size);
                GetRotationTool().Rotate(temp, rRightHandSideVector, GetGeometry());
            }

            if (CalculateStiffnessMatrixFlag == true) {
                for (unsigned int i = 0; i < matrix_size; ++i) {
                    for (unsigned int j = 0; j < matrix_size; ++j) {
                        // erase tangential DoFs
                        if (j > number_of_nodes * dimension)
                            rLeftHandSideMatrix(i, j) = 0.0;
                       
                        if (i > number_of_nodes * dimension)
                            rLeftHandSideMatrix(i, j) = 0.0;
                    }
                }

                for (unsigned int k = 0; k < dimension-1; k++){
                    rLeftHandSideMatrix(number_of_nodes * dimension + k + 1, number_of_nodes * dimension + k + 1) = 1.0;
                }
            }

            if (CalculateResidualVectorFlag == true) {
                for (unsigned int j = 0; j < matrix_size; j++) {
                    if (j % block_size != 0) // tangential DoF
                        rRightHandSideVector[j] = 0.0;
                }
            }

            // rotate back to global frame
            if (CalculateStiffnessMatrixFlag == true){
                GetRotationTool().RevertRotate(rLeftHandSideMatrix, rRightHandSideVector, GetGeometry());
            } else {
                MatrixType temp = ZeroMatrix(matrix_size,matrix_size);
                GetRotationTool().RevertRotate(temp, rRightHandSideVector, GetGeometry());
            }

        }
    }
    
    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
void MPMParticleLagrangeDirichletCondition::CalculateNodalReactions(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
        
    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = GetGeometry().size();

    GeneralVariables Variables;
    array_1d<double, 3 > mpc_force = ZeroVector(3);
    auto pLagrangeNode = r_geometry.GetGeometryParent(0).GetValue(MPC_LAGRANGE_NODE);

    mpc_force = - pLagrangeNode->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);

    // Calculating shape function
    MPMShapeFunctionPointValues(Variables.N);

    // Calculate nodal forces
    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
            r_geometry[i].SetLock();
            r_geometry[i].FastGetSolutionStepValue(REACTION) += mpc_force * Variables.N[i] * this->GetIntegrationWeight();
            r_geometry[i].UnSetLock();
    }
    KRATOS_CATCH( "" )
}

void MPMParticleLagrangeDirichletCondition::FinalizeSolutionStep( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    MPMParticleBaseDirichletCondition::FinalizeSolutionStep(rCurrentProcessInfo);

    // Additional treatment for slip conditions
    if (Is(SLIP))
    {
        GeometryType& r_geometry = GetGeometry();
        auto pLagrangeNode = r_geometry.GetGeometryParent(0).GetValue(MPC_LAGRANGE_NODE);
        pLagrangeNode->SetLock();
        pLagrangeNode->Reset(SLIP);
        pLagrangeNode->FastGetSolutionStepValue(NORMAL).clear();
        pLagrangeNode->UnSetLock();

        const unsigned int number_of_nodes = r_geometry.PointsNumber();

        // Here MPC normal vector and PARTICLE_BASED_SLIP are reset
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            r_geometry[i].SetLock();
            r_geometry[i].Reset(SLIP);
            r_geometry[i].SetValue(PARTICLE_BASED_SLIP, false);
            r_geometry[i].FastGetSolutionStepValue(NORMAL).clear();
            r_geometry[i].UnSetLock();
        }
    }

    KRATOS_CATCH( "" )
}


void MPMParticleLagrangeDirichletCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY

    const GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size() ;
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    if (rResult.size() != dimension * number_of_nodes + dimension)
    {
        rResult.resize(dimension * number_of_nodes + dimension ,false);
    }


    for (unsigned int i = 0; i < number_of_nodes; ++i)
    {
        unsigned int index = i * dimension;
        rResult[index    ] = r_geometry[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index + 1] = r_geometry[i].GetDof(DISPLACEMENT_Y).EquationId();
        if(dimension == 3)
            rResult[index + 2] = r_geometry[i].GetDof(DISPLACEMENT_Z).EquationId();

    }

    unsigned int index = number_of_nodes * dimension;
    auto pLagrangeNode = r_geometry.GetGeometryParent(0).GetValue(MPC_LAGRANGE_NODE);

    rResult[index    ] = pLagrangeNode->GetDof(VECTOR_LAGRANGE_MULTIPLIER_X).EquationId();
    rResult[index + 1] = pLagrangeNode->GetDof(VECTOR_LAGRANGE_MULTIPLIER_Y).EquationId();
    if(dimension == 3)
        rResult[index + 2] = pLagrangeNode->GetDof(VECTOR_LAGRANGE_MULTIPLIER_Z).EquationId();


    KRATOS_CATCH("")
}

void MPMParticleLagrangeDirichletCondition::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY

    const GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension =  r_geometry.WorkingSpaceDimension();
    rElementalDofList.resize(0);
    rElementalDofList.reserve(dimension * number_of_nodes + dimension);

    for (unsigned int i = 0; i < number_of_nodes; ++i)
    {
        rElementalDofList.push_back( r_geometry[i].pGetDof(DISPLACEMENT_X));
        rElementalDofList.push_back( r_geometry[i].pGetDof(DISPLACEMENT_Y));
        if(dimension == 3)
            rElementalDofList.push_back( r_geometry[i].pGetDof(DISPLACEMENT_Z));

    }
    auto pLagrangeNode = r_geometry.GetGeometryParent(0).GetValue(MPC_LAGRANGE_NODE);

    rElementalDofList.push_back(pLagrangeNode->pGetDof(VECTOR_LAGRANGE_MULTIPLIER_X));
    rElementalDofList.push_back(pLagrangeNode->pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Y));
    if(dimension == 3)
        rElementalDofList.push_back(pLagrangeNode->pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Z));



    KRATOS_CATCH("")
}

void MPMParticleLagrangeDirichletCondition::GetValuesVector(
    Vector& rValues,
    int Step
    ) const
{
    const GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    const unsigned int mat_size = number_of_nodes * dimension + dimension;

    if (rValues.size() != mat_size)
    {
        rValues.resize(mat_size, false);
    }
    rValues = ZeroVector(mat_size);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        const array_1d<double, 3 > & r_displacement = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT, Step);

        unsigned int index = i * dimension;

        for(unsigned int k = 0; k < dimension; ++k)
        {
            rValues[index + k] = r_displacement[k];
        }
    }
    auto pLagrangeNode = r_geometry.GetGeometryParent(0).GetValue(MPC_LAGRANGE_NODE);
    const array_1d<double, 3 > & r_lagrange_multiplier = pLagrangeNode->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER, Step);
    const unsigned int lagrange_index = number_of_nodes * dimension;
    for(unsigned int k = 0; k < dimension; ++k)
        {
            rValues[lagrange_index + k] = r_lagrange_multiplier[k];
        }
}

void MPMParticleLagrangeDirichletCondition::GetFirstDerivativesVector(
    Vector& rValues,
    int Step
    ) const
{
    const GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    const unsigned int mat_size = number_of_nodes * dimension + dimension;

    if (rValues.size() != mat_size)
    {
        rValues.resize(mat_size, false);
    }
    rValues = ZeroVector(mat_size);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        const array_1d<double, 3 > & r_velocity = r_geometry[i].FastGetSolutionStepValue(VELOCITY, Step);
        const unsigned int index = i * dimension;
        for(unsigned int k = 0; k < dimension; ++k)
        {
            rValues[index + k] = r_velocity[k];
        }
    }
}

void MPMParticleLagrangeDirichletCondition::GetSecondDerivativesVector(
    Vector& rValues,
    int Step
    ) const
{
    const GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    const unsigned int mat_size = number_of_nodes * dimension + dimension;

    if (rValues.size() != mat_size)
    {
        rValues.resize(mat_size, false);
    }
    rValues = ZeroVector(mat_size);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        const array_1d<double, 3 > & r_acceleration = r_geometry[i].FastGetSolutionStepValue(ACCELERATION, Step);
        const unsigned int index = i * dimension;
        for(unsigned int k = 0; k < dimension; ++k)
        {
            rValues[index + k] = r_acceleration[k];
        }
    }
}


void MPMParticleLagrangeDirichletCondition::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);

    if (rVariable == PENALTY_FACTOR) {
        rValues[0] = m_penalty;
    }
    else {
        MPMParticleBaseDirichletCondition::CalculateOnIntegrationPoints(
            rVariable, rValues, rCurrentProcessInfo);
    }
}


void MPMParticleLagrangeDirichletCondition::SetValuesOnIntegrationPoints(const Variable<double>& rVariable,
    const std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR_IF(rValues.size() > 1)
        << "Only 1 value per integration point allowed! Passed values vector size: "
        << rValues.size() << std::endl;
    
    if (rVariable == PENALTY_FACTOR) {
        m_penalty = rValues[0];
    }
    else {
        MPMParticleBaseDirichletCondition::SetValuesOnIntegrationPoints(
            rVariable, rValues, rCurrentProcessInfo);
    }
}


} // Namespace Kratos