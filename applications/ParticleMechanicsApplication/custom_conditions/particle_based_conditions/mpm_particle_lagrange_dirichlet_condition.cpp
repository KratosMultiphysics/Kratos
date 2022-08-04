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

void MPMParticleLagrangeDirichletCondition::InitializeSolutionStep( const ProcessInfo& rCurrentProcessInfo )
{
    MPMParticleBaseDirichletCondition::InitializeSolutionStep( rCurrentProcessInfo );

    GeometryType& r_geometry = GetGeometry();

    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    
    // Prepare variables
    GeneralVariables Variables;

    // Calculating shape function --> the normal vectors should not be affected by the small cut!
    MPMParticleBaseCondition::MPMShapeFunctionPointValues(Variables.N);

    // Get NODAL_AREA from MPC_Area
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        if (r_geometry[i].SolutionStepsDataHas(NODAL_AREA))
        {
            r_geometry[i].SetLock();
            r_geometry[i].FastGetSolutionStepValue(NODAL_AREA, 0) += Variables.N[i] * this->GetIntegrationWeight();
            r_geometry[i].UnSetLock();
        }
        else break;
    }

    auto pBoundaryParticle = r_geometry.GetGeometryParent(0).GetValue(MPC_LAGRANGE_NODE);
    array_1d<double, 3 > & r_lagrange_multiplier  = pBoundaryParticle->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);

    for ( unsigned int j = 0; j < dimension; j++ )
    {
        r_lagrange_multiplier[0] *= 0.0;
    }

    // Additional treatment for slip conditions
    if (Is(SLIP))
    {
        pBoundaryParticle->SetLock();
        pBoundaryParticle->Set(SLIP);
        pBoundaryParticle->FastGetSolutionStepValue(NORMAL) = m_unit_normal;
        pBoundaryParticle->UnSetLock();


        // Here MPC contribution of normal vector are added
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            r_geometry[i].SetLock();
            r_geometry[i].Set(SLIP);
            // distinguish between penalty and lagrange in bossak scheme
            r_geometry[i].FastGetSolutionStepValue(IS_STRUCTURE) = 3.0;
            r_geometry[i].FastGetSolutionStepValue(NORMAL) += Variables.N[i] * m_unit_normal;
            r_geometry[i].UnSetLock();
        }
    }
    

}

void MPMParticleLagrangeDirichletCondition::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    
    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();

    // At the beginning of NonLinearIteration, REACTION has to be reset to zero
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        r_geometry[i].SetLock();
        r_geometry[i].FastGetSolutionStepValue(REACTION).clear();
        r_geometry[i].UnSetLock();
    }

    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    auto pBoundaryParticle = r_geometry.GetGeometryParent(0).GetValue(MPC_LAGRANGE_NODE);

    array_1d<double, 3 > & r_lagrange_multiplier  = pBoundaryParticle->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);

    for ( unsigned int j = 0; j < dimension; j++ )
    {
        r_lagrange_multiplier[0] *= 0.0;
    }

}

void MPMParticleLagrangeDirichletCondition::MPMShapeFunctionPointValues( Vector& rResult ) const
{
    KRATOS_TRY

    MPMParticleBaseCondition::MPMShapeFunctionPointValues(rResult);

    // Additional check to modify zero shape function values
    // Lagrange Condition is more sensitive for small cut then Penalty condition
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const GeometryType& r_geometry = GetGeometry();

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

    // Nodes with zero mass are not connected to the body--> zero shape function result in zero line and columns in stiffness matrix
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        if (r_geometry[i].FastGetSolutionStepValue(NODAL_MASS, 0) <= std::numeric_limits<double>::epsilon()){
            rResult[i]=0.0;
        }
    }

    

    KRATOS_CATCH( "" )
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
    Variables.CurrentDisp = this->CalculateCurrentDisp(Variables.CurrentDisp, rCurrentProcessInfo);


    bool apply_constraints = true;
    double necessary_nodes=2;

    // avoid overconstrained systems for unstructured triangles
    const GeometryData::KratosGeometryType geo_type = GetGeometry().GetGeometryParent(0).GetGeometryType();
    if (geo_type == GeometryData::Kratos_Triangle2D3)
        necessary_nodes=3;
    else if (geo_type == GeometryData::Kratos_Tetrahedra3D4 || geo_type == GeometryData::Kratos_Hexahedra3D8)
        necessary_nodes=4;

    if (Is(CONTACT))
    {   
        // NOTE: the unit_normal_vector is assumed always pointing outside the boundary
        array_1d<double, 3 > field_displacement = ZeroVector(3);
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            for ( unsigned int j = 0; j < dimension; j++)
            {
                field_displacement[j] += Variables.N[i] * Variables.CurrentDisp(i,j);
            }
        }
            
        const double penetration = MathUtils<double>::Dot((field_displacement - m_imposed_displacement), m_unit_normal);
        // If penetrates, apply constraint, otherwise no
        if (penetration > 0.0)
        {
            apply_constraints = false;
            
        }
        necessary_nodes=number_of_nodes;
    }
    
    int counter = 0;
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        if (r_geometry[i].FastGetSolutionStepValue(NODAL_MASS, 0) >= std::numeric_limits<double>::epsilon()){
            counter +=1;
        }
    }

    // avoid singular matrices --> 2 nodes of the element should be conntected to the body to statisfy the inf-sub condition
    if (counter < necessary_nodes)
        apply_constraints = false;
    

    if (apply_constraints)
    {
        Matrix lagrange_matrix = ZeroMatrix(matrix_size, matrix_size);

        // Loop over shape functions of displacements
        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            const unsigned int ibase = dimension * number_of_nodes;
            for (unsigned int k = 0; k < dimension; k++)
            {
                lagrange_matrix(i* dimension+k, ibase+k) = Variables.N[i];
                lagrange_matrix(ibase+k, i*dimension + k) = Variables.N[i];
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
                    gap_function[index+j] = r_displacement[j] - m_imposed_displacement[j] ;
                }


            }
            right_hand_side = prod(lagrange_matrix, gap_function);

            //first rows of RHS
            gap_function = ZeroVector(matrix_size);
            auto pBoundaryParticle = r_geometry.GetGeometryParent(0).GetValue(MPC_LAGRANGE_NODE);
            const array_1d<double, 3>& r_lagrange_multiplier = pBoundaryParticle->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);

            for (unsigned int j = 0; j < dimension; j++){
                gap_function[dimension * number_of_nodes+j] = r_lagrange_multiplier[j];
            }
            
            right_hand_side += prod(lagrange_matrix, gap_function);

            // //add imposed displacement
            // for (unsigned int j = 0; j < dimension; j++){
            //     right_hand_side[dimension * number_of_nodes+j] -= m_imposed_displacement[j];
            // }

            right_hand_side *= this->GetIntegrationWeight();
            noalias(rRightHandSideVector) = -right_hand_side;
            
        }

    }
    

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
void MPMParticleLagrangeDirichletCondition::FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    MPMParticleBaseDirichletCondition::FinalizeNonLinearIteration(rCurrentProcessInfo);
    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = GetGeometry().size();

    GeneralVariables Variables;
    array_1d<double, 3 > mpc_force = ZeroVector(3);
    auto pBoundaryParticle = r_geometry.GetGeometryParent(0).GetValue(MPC_LAGRANGE_NODE);

    mpc_force = pBoundaryParticle->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);

    // Calculating shape function
    MPMShapeFunctionPointValues(Variables.N);

    // Calculate nodal forces
    Vector nodal_force = ZeroVector(3);
    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        
        // Check whether there nodes are active and associated to material point elements
        const double& nodal_mass = r_geometry[i].FastGetSolutionStepValue(NODAL_MASS, 0);
        if (nodal_mass > std::numeric_limits<double>::epsilon())
        {
            r_geometry[i].SetLock();
            r_geometry[i].FastGetSolutionStepValue(REACTION) += mpc_force * Variables.N[i] * this->GetIntegrationWeight();
            r_geometry[i].UnSetLock();
        }

    }
    KRATOS_CATCH( "" )
}

void MPMParticleLagrangeDirichletCondition::FinalizeSolutionStep( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    MPMParticleBaseDirichletCondition::FinalizeSolutionStep(rCurrentProcessInfo);

    this->CalculateContactForce(rCurrentProcessInfo);
    
    // Additional treatment for slip conditions
    if (Is(SLIP))
    {
        GeometryType& r_geometry = GetGeometry();
        
        auto pBoundaryParticle = r_geometry.GetGeometryParent(0).GetValue(MPC_LAGRANGE_NODE);
        pBoundaryParticle->SetLock();
        pBoundaryParticle->Reset(SLIP);
        pBoundaryParticle->FastGetSolutionStepValue(NORMAL).clear();
        pBoundaryParticle->UnSetLock();

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

void MPMParticleLagrangeDirichletCondition::CalculateContactForce( const ProcessInfo& rCurrentProcessInfo )
{
    
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();

    GeneralVariables Variables;
    Variables.N = row(GetGeometry().ShapeFunctionsValues(), 0);

    array_1d<double, 3 > mpc_force = ZeroVector(3);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        if (Variables.N[i] > std::numeric_limits<double>::epsilon() )
        {
            auto r_geometry = GetGeometry();
            
            const double& nodal_mass = r_geometry[i].FastGetSolutionStepValue(NODAL_MASS, 0);
            const double nodal_area  = r_geometry[i].FastGetSolutionStepValue(NODAL_AREA, 0);
            const Vector nodal_force = r_geometry[i].FastGetSolutionStepValue(REACTION);

            if (nodal_mass > std::numeric_limits<double>::epsilon() && nodal_area > std::numeric_limits<double>::epsilon())
            {
                mpc_force += Variables.N[i] * nodal_force * this->GetIntegrationWeight() / nodal_area;
            }
        }
    }

    // Apply in the normal contact direction and allow releasing motion
    if (Is(CONTACT))
    {
        // Apply only in the normal direction
        const double normal_force = MathUtils<double>::Dot(mpc_force, m_unit_normal);

        // This check is done to avoid sticking forces
        if (normal_force > 0.0)
            mpc_force = 1.0 * normal_force * m_unit_normal;
        else
            mpc_force = ZeroVector(3);
    }
    
    m_contact_force = mpc_force;

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
    auto pBoundaryParticle = r_geometry.GetGeometryParent(0).GetValue(MPC_LAGRANGE_NODE);

    rResult[index    ] = pBoundaryParticle->GetDof(VECTOR_LAGRANGE_MULTIPLIER_X).EquationId();
    rResult[index + 1] = pBoundaryParticle->GetDof(VECTOR_LAGRANGE_MULTIPLIER_Y).EquationId();
    if(dimension == 3)
        rResult[index + 2] = pBoundaryParticle->GetDof(VECTOR_LAGRANGE_MULTIPLIER_Z).EquationId();


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
    auto pBoundaryParticle = r_geometry.GetGeometryParent(0).GetValue(MPC_LAGRANGE_NODE);

    rElementalDofList.push_back(pBoundaryParticle->pGetDof(VECTOR_LAGRANGE_MULTIPLIER_X));
    rElementalDofList.push_back(pBoundaryParticle->pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Y));
    if(dimension == 3)
        rElementalDofList.push_back(pBoundaryParticle->pGetDof(VECTOR_LAGRANGE_MULTIPLIER_Z));



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
    auto pBoundaryParticle = r_geometry.GetGeometryParent(0).GetValue(MPC_LAGRANGE_NODE);
    const array_1d<double, 3 > & r_lagrange_multiplier = pBoundaryParticle->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER, Step);
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

void MPMParticleLagrangeDirichletCondition::CalculateOnIntegrationPoints(
    const Variable<int>& rVariable,
    std::vector<int>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);

    if (rVariable == MPC_CORRESPONDING_CONDITION_ID) {
        rValues[0] = m_corresponding_condition_id;
    }
    else if (rVariable == MPC_CORRESPONDING_NODE_ID) {
        rValues[0] = m_corresponding_node_id;
    }
    else if (rVariable == MPC_COUNTER) {
        rValues[0] = m_counter;
    }
    else
    {
        KRATOS_ERROR << "Variable " << rVariable << " is called in CalculateOnIntegrationPoints, but is not implemented." << std::endl;
    }
}

void MPMParticleLagrangeDirichletCondition::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);

    else {
        MPMParticleBaseDirichletCondition::CalculateOnIntegrationPoints(
            rVariable, rValues, rCurrentProcessInfo);
    }
}

void MPMParticleLagrangeDirichletCondition::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
    std::vector<array_1d<double, 3 > >& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);

    if (rVariable == MPC_NORMAL) {
        rValues[0] = m_unit_normal;
    }
    else if (rVariable == MPC_CONTACT_FORCE) {
        this->CalculateContactForce(rCurrentProcessInfo);
        rValues[0] = m_contact_force;
    }
    else {
        MPMParticleBaseDirichletCondition::CalculateOnIntegrationPoints(
            rVariable, rValues, rCurrentProcessInfo);
    }
}

void MPMParticleLagrangeDirichletCondition::SetValuesOnIntegrationPoints(
    const Variable<int>& rVariable,
    const std::vector<int>& rValues,
    const ProcessInfo& rCurrentProcessInfo) {
    
    KRATOS_ERROR_IF(rValues.size() > 1)
        << "Only 1 value per integration point allowed! Passed values vector size: "
        << rValues.size() << std::endl;

    if (rVariable == MPC_CORRESPONDING_CONDITION_ID) {
        m_corresponding_condition_id = rValues[0];
    }
    else if (rVariable == MPC_CORRESPONDING_NODE_ID) {
        m_corresponding_node_id = rValues[0];
    }
    else if (rVariable == MPC_COUNTER) {
        m_counter = rValues[0];
    }
    else
    {
        KRATOS_ERROR << "Variable " << rVariable << " is called in SetValuesIntegrationPoints, but is not implemented." << std::endl;
    }
}


void MPMParticleLagrangeDirichletCondition::SetValuesOnIntegrationPoints(const Variable<double>& rVariable,
    const std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR_IF(rValues.size() > 1)
        << "Only 1 value per integration point allowed! Passed values vector size: "
        << rValues.size() << std::endl;

    else {
        MPMParticleBaseDirichletCondition::SetValuesOnIntegrationPoints(
            rVariable, rValues, rCurrentProcessInfo);
    }
}

void MPMParticleLagrangeDirichletCondition::SetValuesOnIntegrationPoints(
    const Variable<array_1d<double, 3 > >& rVariable,
    const std::vector<array_1d<double, 3 > >& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR_IF(rValues.size() > 1)
        << "Only 1 value per integration point allowed! Passed values vector size: "
        << rValues.size() << std::endl;

    if (rVariable == MPC_NORMAL) {
        m_unit_normal = rValues[0];
        ParticleMechanicsMathUtilities<double>::Normalize(m_unit_normal);
    }
    else if (rVariable == MPC_CONTACT_FORCE) {
        m_contact_force = rValues[0];
    }
    else {
        MPMParticleBaseDirichletCondition::SetValuesOnIntegrationPoints(
            rVariable, rValues, rCurrentProcessInfo);
    }
}

int MPMParticleLagrangeDirichletCondition::Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    MPMParticleBaseDirichletCondition::Check(rCurrentProcessInfo);

    // Verify that the dofs exist
    for (const auto& r_node : this->GetGeometry().Points())
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NORMAL,r_node)

    return 0;
}

} // Namespace Kratos

