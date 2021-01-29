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
#include "custom_conditions/grid_based_conditions/mpm_grid_base_load_condition.h"
#include "includes/checks.h"

namespace Kratos
{
    //************************************************************************************
    //************************************************************************************

    void MPMGridBaseLoadCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const
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
    void MPMGridBaseLoadCondition::GetDofList(
        DofsVectorType& ElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
        ) const
    {
        KRATOS_TRY

        const GeometryType& r_geometry = GetGeometry();
        const unsigned int number_of_nodes = r_geometry.size();
        const unsigned int dimension =  r_geometry.WorkingSpaceDimension();
        ElementalDofList.resize(0);
        ElementalDofList.reserve(dimension * number_of_nodes);

        if(dimension == 2)
        {
            for (unsigned int i = 0; i < number_of_nodes; ++i)
            {
                ElementalDofList.push_back( r_geometry[i].pGetDof(DISPLACEMENT_X));
                ElementalDofList.push_back( r_geometry[i].pGetDof(DISPLACEMENT_Y));
            }
        }
        else
        {
            for (unsigned int i = 0; i < number_of_nodes; ++i)
            {
                ElementalDofList.push_back( r_geometry[i].pGetDof(DISPLACEMENT_X));
                ElementalDofList.push_back( r_geometry[i].pGetDof(DISPLACEMENT_Y));
                ElementalDofList.push_back( r_geometry[i].pGetDof(DISPLACEMENT_Z));
            }
        }
        KRATOS_CATCH("")
    }

    //***********************************************************************
    //***********************************************************************

    void MPMGridBaseLoadCondition::GetValuesVector(
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
            const array_1d<double, 3 > & r_displacement = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
            unsigned int index = i * dimension;
            for(unsigned int k = 0; k < dimension; ++k)
            {
                rValues[index + k] = r_displacement[k];
            }
        }
    }

    //***********************************************************************
    //***********************************************************************

    void MPMGridBaseLoadCondition::GetFirstDerivativesVector(
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
            const array_1d<double, 3 > & r_velocity = r_geometry[i].FastGetSolutionStepValue(VELOCITY, Step);
            const unsigned int index = i * dimension;
            for(unsigned int k = 0; k<dimension; ++k)
            {
                rValues[index + k] = r_velocity[k];
            }
        }
    }

    //***********************************************************************
    //***********************************************************************

    void MPMGridBaseLoadCondition::GetSecondDerivativesVector(
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
            const array_1d<double, 3 > & r_acceleration = r_geometry[i].FastGetSolutionStepValue(ACCELERATION, Step);
            const unsigned int index = i * dimension;
            for(unsigned int k = 0; k < dimension; ++k)
            {
                rValues[index + k] = r_acceleration[k];
            }
        }
    }

    //************************************************************************************
    //************************************************************************************

    void MPMGridBaseLoadCondition::CalculateRightHandSide( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
    {
        // Calculation flags
        const bool CalculateStiffnessMatrixFlag = false;
        const bool CalculateResidualVectorFlag = true;
        MatrixType temp = Matrix();

        CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
    }

    //************************************************************************************
    //************************************************************************************
    void MPMGridBaseLoadCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
    {
        //calculation flags
        const bool CalculateStiffnessMatrixFlag = true;
        const bool CalculateResidualVectorFlag = true;

        CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
    }

    //***********************************************************************
    //***********************************************************************

    void MPMGridBaseLoadCondition::CalculateMassMatrix(
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

    void MPMGridBaseLoadCondition::CalculateDampingMatrix(
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

    void MPMGridBaseLoadCondition::CalculateAll(
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

    int MPMGridBaseLoadCondition::Check( const ProcessInfo& rCurrentProcessInfo ) const
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

    double MPMGridBaseLoadCondition::GetIntegrationWeight(
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
        const unsigned int PointNumber,
        const double detJ
        )
    {
        return IntegrationPoints[PointNumber].Weight() * detJ;
    }

    //***********************************************************************************
    //***********************************************************************************

    void MPMGridBaseLoadCondition::AddExplicitContribution(const VectorType& rRHS,
        const Variable<VectorType>& rRHSVariable,
        const Variable<array_1d<double,3> >& rDestinationVariable,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        GeometryType& r_geometry = GetGeometry();
        const unsigned int number_of_nodes = r_geometry.PointsNumber();
        const unsigned int dimension       = r_geometry.WorkingSpaceDimension();

        if( rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL )
        {

            for(SizeType i=0; i< number_of_nodes; i++)
            {
                SizeType index = dimension * i;

                if (r_geometry[i].SolutionStepsDataHas(FORCE_RESIDUAL))
                {
                    array_1d<double, 3 > & r_force_residual = r_geometry[i].FastGetSolutionStepValue(FORCE_RESIDUAL);
                    for(SizeType j=0; j<dimension; j++)
                    {
                        #pragma omp atomic
                        r_force_residual[j] += rRHS[index + j];
                    }
                }
            }
        }

    KRATOS_CATCH( "" )
    }

} // Namespace Kratos


