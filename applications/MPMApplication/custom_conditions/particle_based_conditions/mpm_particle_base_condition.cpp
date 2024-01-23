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
#include "custom_conditions/particle_based_conditions/mpm_particle_base_condition.h"
#include "includes/checks.h"
#include "custom_utilities/mpm_math_utilities.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************

void MPMParticleBaseCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY

    const GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    if (rResult.size() != dimension * number_of_nodes)
    {
        rResult.resize(dimension*number_of_nodes,false);
    }

    const unsigned int pos = r_geometry[0].GetDofPosition(DISPLACEMENT_X);

    for (unsigned int i = 0; i < number_of_nodes; ++i)
    {
        const unsigned int index = i * dimension;
        rResult[index    ] = r_geometry[i].GetDof(DISPLACEMENT_X,pos    ).EquationId();
        rResult[index + 1] = r_geometry[i].GetDof(DISPLACEMENT_Y,pos + 1).EquationId();
        if ( dimension == 3 )
            rResult[index + 2] = r_geometry[i].GetDof( DISPLACEMENT_Z,pos + 2 ).EquationId();
    }

    KRATOS_CATCH("")
}

//***********************************************************************
//***********************************************************************
void MPMParticleBaseCondition::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY

    const GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension =  r_geometry.WorkingSpaceDimension();
    rElementalDofList.resize(0);
    rElementalDofList.reserve(dimension * number_of_nodes);

    for (unsigned int i = 0; i < number_of_nodes; ++i)
    {
        rElementalDofList.push_back( r_geometry[i].pGetDof(DISPLACEMENT_X));
        rElementalDofList.push_back( r_geometry[i].pGetDof(DISPLACEMENT_Y));
        if ( dimension == 3 ){
            rElementalDofList.push_back( r_geometry[i].pGetDof( DISPLACEMENT_Z ) );
        }

    }

    KRATOS_CATCH("")
}

//***********************************************************************
//***********************************************************************

void MPMParticleBaseCondition::GetValuesVector(
    Vector& rValues,
    int Step
    ) const
{
    const GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    const unsigned int matrix_size = number_of_nodes * dimension;

    if (rValues.size() != matrix_size)
    {
        rValues.resize(matrix_size, false);
    }

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        const array_1d<double, 3 > & Displacement = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
        unsigned int index = i * dimension;
        for(unsigned int k = 0; k < dimension; ++k)
        {
            rValues[index + k] = Displacement[k];
        }
    }
}

//***********************************************************************
//***********************************************************************

void MPMParticleBaseCondition::GetFirstDerivativesVector(
    Vector& rValues,
    int Step
    ) const
{
    const GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    const unsigned int matrix_size = number_of_nodes * dimension;

    if (rValues.size() != matrix_size)
    {
        rValues.resize(matrix_size, false);
    }

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        const array_1d<double, 3 > & Velocity = r_geometry[i].FastGetSolutionStepValue(VELOCITY, Step);
        const unsigned int index = i * dimension;
        for(unsigned int k = 0; k<dimension; ++k)
        {
            rValues[index + k] = Velocity[k];
        }
    }
}

//***********************************************************************
//***********************************************************************

void MPMParticleBaseCondition::GetSecondDerivativesVector(
    Vector& rValues,
    int Step
    ) const
{
    const GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    const unsigned int matrix_size = number_of_nodes * dimension;

    if (rValues.size() != matrix_size)
    {
        rValues.resize(matrix_size, false);
    }

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

//************************************************************************************
//************************************************************************************

void MPMParticleBaseCondition::CalculateRightHandSide( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{
    // Calculation flags
    const bool CalculateStiffnessMatrixFlag = false;
    const bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
}

//************************************************************************************
//************************************************************************************
void MPMParticleBaseCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{
    //calculation flags
    const bool CalculateStiffnessMatrixFlag = true;
    const bool CalculateResidualVectorFlag = true;

    CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
}

//***********************************************************************
//***********************************************************************

void MPMParticleBaseCondition::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if(rMassMatrix.size1() != 0)
    {
        rMassMatrix.resize(0, 0, false);
    }
}

//***********************************************************************
//***********************************************************************

void MPMParticleBaseCondition::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if(rDampingMatrix.size1() != 0)
    {
        rDampingMatrix.resize(0, 0, false);
    }
}

//***********************************************************************
//***********************************************************************

void MPMParticleBaseCondition::CalculateAll(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    bool CalculateStiffnessMatrixFlag,
    bool CalculateResidualVectorFlag
    )
{
    KRATOS_ERROR << "You are calling the CalculateAll from the base class for loads" << std::endl;
}

//***********************************************************************
//***********************************************************************

int MPMParticleBaseCondition::Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    // Base check
    Condition::Check(rCurrentProcessInfo);

    // Check that the condition's nodes contain all required SolutionStepData and Degrees of freedom
    for (const auto& r_node : this->GetGeometry().Points()) {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,r_node)

        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, r_node)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, r_node)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, r_node)
    }

    return 0;
}

//***********************************************************************
//***********************************************************************
/**
   * Shape function values in given point. This method calculate the shape function
   * vector in given point.
   *
 */
void MPMParticleBaseCondition::MPMShapeFunctionPointValues(Vector& rResult) const
{
    KRATOS_TRY

    rResult = row(GetGeometry().ShapeFunctionsValues(), 0);

    KRATOS_CATCH( "" )
}

//*************************COMPUTE CURRENT DISPLACEMENT*******************************
//************************************************************************************
/*
This function convert the computed nodal displacement into matrix of (number_of_nodes, dimension)
*/
Matrix& MPMParticleBaseCondition::CalculateCurrentDisp(Matrix & rCurrentDisp, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();

    rCurrentDisp = ZeroMatrix(number_of_nodes, dimension);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        const array_1d<double, 3 > & current_displacement  = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT);

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            rCurrentDisp(i,j) = current_displacement[j];
        }
    }

    return rCurrentDisp;

    KRATOS_CATCH( "" )
}

void MPMParticleBaseCondition::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);

    if (rVariable == MPC_AREA) {
        rValues[0] = m_area;
    }
    else {
        KRATOS_ERROR << "Variable " << rVariable << " is called in CalculateOnIntegrationPoints, but is not implemented." << std::endl;
    }
}

void MPMParticleBaseCondition::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
    std::vector<array_1d<double, 3 > >& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);

    if (rVariable == MP_COORD || rVariable == MPC_COORD) {
        rValues[0] = m_xg;
    }
    else if (rVariable == MPC_DISPLACEMENT) {
        rValues[0] = m_displacement;
    }
    else if (rVariable == MPC_VELOCITY) {
        rValues[0] = m_velocity;
    }
    else if (rVariable == MPC_ACCELERATION) {
        rValues[0] = m_acceleration;
    }
    else if (rVariable == MPC_NORMAL ) {
        rValues[0] = m_normal;
    }
    else {
        KRATOS_ERROR << "Variable " << rVariable << " is called in CalculateOnIntegrationPoints, but is not implemented." << std::endl;
    }
}

void MPMParticleBaseCondition::SetValuesOnIntegrationPoints(
    const Variable<double>& rVariable,
    const std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo) {
    KRATOS_ERROR_IF(rValues.size() > 1)
        << "Only 1 value per integration point allowed! Passed values vector size: "
        << rValues.size() << std::endl;

    if (rVariable == MPC_AREA) {
        m_area = rValues[0];
    }
    else {
        KRATOS_ERROR << "Variable " << rVariable << " is called in SetValuesOnIntegrationPoints, but is not implemented." << std::endl;
    }
}

void MPMParticleBaseCondition::SetValuesOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
    const std::vector<array_1d<double, 3 > >& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR_IF(rValues.size() > 1)
        << "Only 1 value per integration point allowed! Passed values vector size: "
        << rValues.size() << std::endl;

    if (rVariable == MP_COORD || rVariable == MPC_COORD) {
        m_xg = rValues[0];
    }
    else if (rVariable == MPC_DISPLACEMENT) {
        m_displacement = rValues[0];
    }
    else if (rVariable == MPC_VELOCITY) {
        m_velocity = rValues[0];
    }
    else if (rVariable == MPC_ACCELERATION) {
        m_acceleration = rValues[0];
    }
    else if (rVariable == MPC_NORMAL) {
        m_normal = rValues[0];
        MPMMathUtilities<double>::Normalize(m_normal);
    }
    else {
        KRATOS_ERROR << "Variable " << rVariable << " is called in SetValuesOnIntegrationPoints, but is not implemented." << std::endl;
    }
}

//***********************************************************************
//***********************************************************************

double MPMParticleBaseCondition::GetIntegrationWeight()
{
    return m_area;
}

} // Namespace Kratos


