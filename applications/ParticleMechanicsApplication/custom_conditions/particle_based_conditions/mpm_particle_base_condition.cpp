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

namespace Kratos
{

//************************************************************************************
//************************************************************************************

void MPMParticleBaseCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    if (rResult.size() != dimension * number_of_nodes)
    {
        rResult.resize(dimension*number_of_nodes,false);
    }

    const unsigned int pos = r_geometry[0].GetDofPosition(DISPLACEMENT_X);

    if(dimension == 2)
    {
        for (unsigned int i = 0; i < number_of_nodes; ++i)
        {
            const unsigned int index = i * 2;
            rResult[index    ] = r_geometry[i].GetDof(DISPLACEMENT_X,pos    ).EquationId();
            rResult[index + 1] = r_geometry[i].GetDof(DISPLACEMENT_Y,pos + 1).EquationId();
        }
    }
    else
    {
        for (unsigned int i = 0; i < number_of_nodes; ++i)
        {
            const unsigned int index = i * 3;
            rResult[index    ] = r_geometry[i].GetDof(DISPLACEMENT_X,pos    ).EquationId();
            rResult[index + 1] = r_geometry[i].GetDof(DISPLACEMENT_Y,pos + 1).EquationId();
            rResult[index + 2] = r_geometry[i].GetDof(DISPLACEMENT_Z,pos + 2).EquationId();
        }
    }
    KRATOS_CATCH("")
}

//***********************************************************************
//***********************************************************************
void MPMParticleBaseCondition::GetDofList(
    DofsVectorType& rElementalDofList,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension =  r_geometry.WorkingSpaceDimension();
    rElementalDofList.resize(0);
    rElementalDofList.reserve(dimension * number_of_nodes);

    if(dimension == 2)
    {
        for (unsigned int i = 0; i < number_of_nodes; ++i)
        {
            rElementalDofList.push_back( r_geometry[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back( r_geometry[i].pGetDof(DISPLACEMENT_Y));
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
    }
    KRATOS_CATCH("")
}

//***********************************************************************
//***********************************************************************

void MPMParticleBaseCondition::GetValuesVector(
    Vector& rValues,
    int Step
    )
{
    GeometryType& r_geometry = GetGeometry();
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
    )
{
    GeometryType& r_geometry = GetGeometry();
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
    )
{
    GeometryType& r_geometry = GetGeometry();
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

void MPMParticleBaseCondition::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    // Calculation flags
    const bool CalculateStiffnessMatrixFlag = false;
    const bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
}

//************************************************************************************
//************************************************************************************
void MPMParticleBaseCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
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
    ProcessInfo& rCurrentProcessInfo
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
    ProcessInfo& rCurrentProcessInfo
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
    ProcessInfo& rCurrentProcessInfo,
    bool CalculateStiffnessMatrixFlag,
    bool CalculateResidualVectorFlag
    )
{
    KRATOS_ERROR << "You are calling the CalculateAll from the base class for loads" << std::endl;
}

//***********************************************************************
//***********************************************************************

int MPMParticleBaseCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
    // Base check
    Condition::Check(rCurrentProcessInfo);

    // Verify variable exists
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)

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
   * @param rPoint point which shape function values have to
   * be calculated in it.
   *
   * @return Vector of double which is shape function vector \f$ N \f$ in given point.
   *
 */
Vector& MPMParticleBaseCondition::MPMShapeFunctionPointValues(Vector& rResult, const array_1d<double,3>& rPoint)
{
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    rResult.resize(number_of_nodes, false);

    // Get local point coordinate
    array_1d<double,3> rPointLocal = ZeroVector(3);
    rPointLocal = GetGeometry().PointLocalCoordinates(rPointLocal, rPoint);

    if (dimension == 2)
    {
        // Get Shape functions: N depending on number of nodes
        switch (number_of_nodes)
        {
            case 3:
                rResult[0] = 1 - rPointLocal[0] - rPointLocal[1] ;
                rResult[1] = rPointLocal[0];
                rResult[2] = rPointLocal[1];
                break;
            case 4:
                rResult[0] = 0.25 * (1 - rPointLocal[0]) * (1 - rPointLocal[1]) ;
                rResult[1] = 0.25 * (1 + rPointLocal[0]) * (1 - rPointLocal[1]) ;
                rResult[2] = 0.25 * (1 + rPointLocal[0]) * (1 + rPointLocal[1]) ;
                rResult[3] = 0.25 * (1 - rPointLocal[0]) * (1 + rPointLocal[1]) ;
                break;
        }
    }
    else if (dimension == 3)
    {
        // Get Shape functions: N depending on number of nodes
        switch (number_of_nodes)
        {
            case 4:
                rResult[0] =  1.0-(rPointLocal[0]+rPointLocal[1]+rPointLocal[2]) ;
                rResult[1] = rPointLocal[0];
                rResult[2] = rPointLocal[1];
                rResult[3] = rPointLocal[2];
                break;
            case 8:
                rResult[0] = 0.125 * (1 - rPointLocal[0]) * (1 - rPointLocal[1]) * (1 - rPointLocal[2]) ;
                rResult[1] = 0.125 * (1 + rPointLocal[0]) * (1 - rPointLocal[1]) * (1 - rPointLocal[2]) ;
                rResult[2] = 0.125 * (1 + rPointLocal[0]) * (1 + rPointLocal[1]) * (1 - rPointLocal[2]) ;
                rResult[3] = 0.125 * (1 - rPointLocal[0]) * (1 + rPointLocal[1]) * (1 - rPointLocal[2]) ;
                rResult[4] = 0.125 * (1 - rPointLocal[0]) * (1 - rPointLocal[1]) * (1 + rPointLocal[2]) ;
                rResult[5] = 0.125 * (1 + rPointLocal[0]) * (1 - rPointLocal[1]) * (1 + rPointLocal[2]) ;
                rResult[6] = 0.125 * (1 + rPointLocal[0]) * (1 + rPointLocal[1]) * (1 + rPointLocal[2]) ;
                rResult[7] = 0.125 * (1 - rPointLocal[0]) * (1 + rPointLocal[1]) * (1 + rPointLocal[2]) ;
                break;
        }
    }

    return rResult;

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
    std::vector<double>& rValues,
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
    std::vector<array_1d<double, 3 > > rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR_IF(rValues.size() > 1)
        << "Only 1 value per integration point allowed! Passed values vector size: "
        << rValues.size() << std::endl;

    if (rVariable == MP_COORD || rVariable == MPC_COORD) {
        m_xg = rValues[0];
    }
    else if (rVariable == MPC_VELOCITY) {
        m_velocity = rValues[0];
    }
    else if (rVariable == MPC_ACCELERATION) {
        m_acceleration = rValues[0];
    }
    else if (rVariable == MPC_NORMAL) {
        m_normal = rValues[0];
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


