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
#include "custom_conditions/particle_based_conditions/mpm_particle_base_dirichlet_condition.h"

namespace Kratos
{
//************************************************************************************
//************************************************************************************

void MPMParticleBaseDirichletCondition::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    /* NOTE:
    In the InitializeSolutionStep of each time step the nodal initial conditions are evaluated.
    This function is called by the base scheme class.*/
    // Here MPC_DISPLACEMENT is updated in terms of velocity and acceleration is added
    array_1d<double,3>& MPC_Displacement = this->GetValue(MPC_DISPLACEMENT);
    array_1d<double,3>& MPC_Velocity = this->GetValue(MPC_VELOCITY);
    const array_1d<double,3>& MPC_Acceleration = this->GetValue(MPC_ACCELERATION);
    const double& delta_time = rCurrentProcessInfo[DELTA_TIME];

    MPC_Displacement += (MPC_Velocity * delta_time) + (0.5 * MPC_Acceleration * delta_time * delta_time);
    MPC_Velocity += (MPC_Acceleration * delta_time);
}

//************************************************************************************
//************************************************************************************

void MPMParticleBaseDirichletCondition::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // Update the MPC Position
    const array_1d<double,3> & xg_c = this->GetValue(MPC_COORD);
    array_1d<double,3> & displacement = this->GetValue(MPC_DISPLACEMENT);
    const array_1d<double,3> & new_xg_c = xg_c + displacement ;
    this->SetValue(MPC_COORD,new_xg_c);

    // Set displacement to zero (NOTE: to use incremental displacement, use MPC_VELOCITY)
    displacement.clear();

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void MPMParticleBaseDirichletCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    GeometryType& rGeom = GetGeometry();
    const unsigned int NumberOfNodes = rGeom.size();
    const unsigned int dim = rGeom.WorkingSpaceDimension();
    if (rResult.size() != dim * NumberOfNodes)
    {
        rResult.resize(dim*NumberOfNodes,false);
    }

    const unsigned int pos = rGeom[0].GetDofPosition(DISPLACEMENT_X);

    if(dim == 2)
    {
        for (unsigned int i = 0; i < NumberOfNodes; ++i)
        {
            const unsigned int index = i * 2;
            rResult[index    ] = rGeom[i].GetDof(DISPLACEMENT_X,pos    ).EquationId();
            rResult[index + 1] = rGeom[i].GetDof(DISPLACEMENT_Y,pos + 1).EquationId();
        }
    }
    else
    {
        for (unsigned int i = 0; i < NumberOfNodes; ++i)
        {
            const unsigned int index = i * 3;
            rResult[index    ] = rGeom[i].GetDof(DISPLACEMENT_X,pos    ).EquationId();
            rResult[index + 1] = rGeom[i].GetDof(DISPLACEMENT_Y,pos + 1).EquationId();
            rResult[index + 2] = rGeom[i].GetDof(DISPLACEMENT_Z,pos + 2).EquationId();
        }
    }
    KRATOS_CATCH("")
}

//***********************************************************************
//***********************************************************************
void MPMParticleBaseDirichletCondition::GetDofList(
    DofsVectorType& rElementalDofList,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    GeometryType& rGeom = GetGeometry();
    const unsigned int NumberOfNodes = rGeom.size();
    const unsigned int dim =  rGeom.WorkingSpaceDimension();
    rElementalDofList.resize(0);
    rElementalDofList.reserve(dim * NumberOfNodes);

    if(dim == 2)
    {
        for (unsigned int i = 0; i < NumberOfNodes; ++i)
        {
            rElementalDofList.push_back( rGeom[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back( rGeom[i].pGetDof(DISPLACEMENT_Y));
        }
    }
    else
    {
        for (unsigned int i = 0; i < NumberOfNodes; ++i)
        {
            rElementalDofList.push_back( rGeom[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back( rGeom[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back( rGeom[i].pGetDof(DISPLACEMENT_Z));
        }
    }
    KRATOS_CATCH("")
}

//***********************************************************************
//***********************************************************************

void MPMParticleBaseDirichletCondition::GetValuesVector(
    Vector& rValues,
    int Step
    )
{
    GeometryType& rGeom = GetGeometry();
    const unsigned int NumberOfNodes = rGeom.size();
    const unsigned int dim = rGeom.WorkingSpaceDimension();
    const unsigned int MatSize = NumberOfNodes * dim;

    if (rValues.size() != MatSize)
    {
        rValues.resize(MatSize, false);
    }

    for (unsigned int i = 0; i < NumberOfNodes; i++)
    {
        const array_1d<double, 3 > & Displacement = rGeom[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
        unsigned int index = i * dim;
        for(unsigned int k = 0; k < dim; ++k)
        {
            rValues[index + k] = Displacement[k];
        }
    }
}

//***********************************************************************
//***********************************************************************

void MPMParticleBaseDirichletCondition::GetFirstDerivativesVector(
    Vector& rValues,
    int Step
    )
{
    GeometryType& rGeom = GetGeometry();
    const unsigned int NumberOfNodes = rGeom.size();
    const unsigned int dim = rGeom.WorkingSpaceDimension();
    const unsigned int MatSize = NumberOfNodes * dim;

    if (rValues.size() != MatSize)
    {
        rValues.resize(MatSize, false);
    }

    for (unsigned int i = 0; i < NumberOfNodes; i++)
    {
        const array_1d<double, 3 > & Velocity = rGeom[i].FastGetSolutionStepValue(VELOCITY, Step);
        const unsigned int index = i * dim;
        for(unsigned int k = 0; k<dim; ++k)
        {
            rValues[index + k] = Velocity[k];
        }
    }
}

//***********************************************************************
//***********************************************************************

void MPMParticleBaseDirichletCondition::GetSecondDerivativesVector(
    Vector& rValues,
    int Step
    )
{
    GeometryType& rGeom = GetGeometry();
    const unsigned int NumberOfNodes = rGeom.size();
    const unsigned int dim = rGeom.WorkingSpaceDimension();
    const unsigned int MatSize = NumberOfNodes * dim;

    if (rValues.size() != MatSize)
    {
        rValues.resize(MatSize, false);
    }

    for (unsigned int i = 0; i < NumberOfNodes; i++)
    {
        const array_1d<double, 3 > & Acceleration = rGeom[i].FastGetSolutionStepValue(ACCELERATION, Step);
        const unsigned int index = i * dim;
        for(unsigned int k = 0; k < dim; ++k)
        {
            rValues[index + k] = Acceleration[k];
        }
    }
}

//************************************************************************************
//************************************************************************************

void MPMParticleBaseDirichletCondition::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    // Calculation flags
    const bool CalculateStiffnessMatrixFlag = false;
    const bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
}

//************************************************************************************
//************************************************************************************
void MPMParticleBaseDirichletCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    //calculation flags
    const bool CalculateStiffnessMatrixFlag = true;
    const bool CalculateResidualVectorFlag = true;

    CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
}

//***********************************************************************
//***********************************************************************

void MPMParticleBaseDirichletCondition::CalculateMassMatrix(
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

void MPMParticleBaseDirichletCondition::CalculateDampingMatrix(
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

void MPMParticleBaseDirichletCondition::CalculateAll(
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

int MPMParticleBaseDirichletCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
    if ( DISPLACEMENT.Key() == 0 )
    {
        KRATOS_ERROR <<  "DISPLACEMENT has Key zero! (check if the application is correctly registered" << std::endl;
    }

    //verify that the dofs exist
    for ( unsigned int i = 0; i < this->GetGeometry().size(); i++ )
    {
        if ( this->GetGeometry()[i].SolutionStepsDataHas( DISPLACEMENT ) == false )
        {
            KRATOS_ERROR << "missing variable DISPLACEMENT on node " << this->GetGeometry()[i].Id() << std::endl;
        }

        if ( this->GetGeometry()[i].HasDofFor( DISPLACEMENT_X ) == false ||
                this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Y ) == false ||
                this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Z ) == false )
        {
            KRATOS_ERROR << "missing one of the dofs for the variable DISPLACEMENT on node " << GetGeometry()[i].Id() << " of condition " << Id() << std::endl;
        }
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
Vector& MPMParticleBaseDirichletCondition::MPMShapeFunctionPointValues(Vector& rResult, const array_1d<double,3>& rPoint)
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
Matrix& MPMParticleBaseDirichletCondition::CalculateCurrentDisp(Matrix & rCurrentDisp, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    GeometryType& rGeom = GetGeometry();
    const unsigned int number_of_nodes = rGeom.PointsNumber();
    const unsigned int dimension = rGeom.WorkingSpaceDimension();

    rCurrentDisp = ZeroMatrix(number_of_nodes, dimension);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        const array_1d<double, 3 > & current_displacement  = rGeom[i].FastGetSolutionStepValue(DISPLACEMENT);

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            rCurrentDisp(i,j) = current_displacement[j];
        }
    }

    return rCurrentDisp;

    KRATOS_CATCH( "" )
}

//***********************************************************************
//***********************************************************************

double MPMParticleBaseDirichletCondition::GetIntegrationWeight()
{
    return this->GetValue(MPC_AREA);
}

} // Namespace Kratos


