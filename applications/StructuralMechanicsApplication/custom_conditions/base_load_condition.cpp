// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "custom_conditions/base_load_condition.h"

namespace Kratos
{
    void BaseLoadCondition::Initialize()
    {
        // TODO: Add somethig if necessary
    }

    //************************************************************************************
    //************************************************************************************

    void BaseLoadCondition::InitializeSolutionStep( ProcessInfo& CurrentProcessInfo )
    {
        // TODO: Add somethig if necessary
    }

    //************************************************************************************
    //************************************************************************************

    void BaseLoadCondition::InitializeNonLinearIteration( ProcessInfo& CurrentProcessInfo )
    {
        // TODO: Add somethig if necessary
    }

    //************************************************************************************
    //************************************************************************************

    void BaseLoadCondition::FinalizeNonLinearIteration( ProcessInfo& CurrentProcessInfo )
    {
        // TODO: Add somethig if necessary
    }

    //************************************************************************************
    //************************************************************************************

    void BaseLoadCondition::FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo )
    {
        // TODO: Add somethig if necessary
    }

    //************************************************************************************
    //************************************************************************************

    void BaseLoadCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        const std::size_t number_of_nodes = GetGeometry().size();
        const std::size_t dim = GetGeometry().WorkingSpaceDimension();
        const std::size_t block_size = this->GetBlockSize();
        if (rResult.size() != dim * number_of_nodes) {
            rResult.resize(number_of_nodes * block_size, false);
        }

        const std::size_t pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

        if(dim == 2) {
            for (std::size_t i = 0; i < number_of_nodes; ++i) {
                const std::size_t index = i * block_size;
                rResult[index    ] = GetGeometry()[i].GetDof(DISPLACEMENT_X,pos    ).EquationId();
                rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y,pos + 1).EquationId();
                if (this->HasRotDof())
                    rResult[index + 2] = GetGeometry()[i].GetDof(ROTATION_Z,pos + 2).EquationId();
            }
        } else {
            for (std::size_t i = 0; i < number_of_nodes; ++i) {
                const std::size_t index = i * block_size;
                rResult[index    ] = GetGeometry()[i].GetDof(DISPLACEMENT_X,pos    ).EquationId();
                rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y,pos + 1).EquationId();
                rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z,pos + 2).EquationId();
            }
        }
        KRATOS_CATCH("")
    }

    //***********************************************************************
    //***********************************************************************
    void BaseLoadCondition::GetDofList(
        DofsVectorType& ElementalDofList,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_TRY

        const std::size_t number_of_nodes = GetGeometry().size();
        const std::size_t dim =  GetGeometry().WorkingSpaceDimension();
        const std::size_t block_size = this->GetBlockSize();
        ElementalDofList.resize(0);
        ElementalDofList.reserve(number_of_nodes * block_size);

        if(dim == 2) {
            for (std::size_t i = 0; i < number_of_nodes; ++i) {
                ElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_X));
                ElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
                if (this->HasRotDof())
                    ElementalDofList.push_back( GetGeometry()[i].pGetDof(ROTATION_Z));
            }
        } else {
            for (std::size_t i = 0; i < number_of_nodes; ++i) {
                ElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_X));
                ElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
                ElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
            }
        }
        KRATOS_CATCH("")
    }

    //***********************************************************************
    //***********************************************************************

    void BaseLoadCondition::GetValuesVector(
        Vector& rValues,
        int Step
        )
    {
        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        const unsigned int MatSize = NumberOfNodes * dim;

        if (rValues.size() != MatSize)
        {
            rValues.resize(MatSize, false);
        }

        for (unsigned int i = 0; i < NumberOfNodes; i++)
        {
            const array_1d<double, 3 > & Displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
            unsigned int index = i * dim;
            for(unsigned int k = 0; k < dim; ++k)
            {
                rValues[index + k] = Displacement[k];
            }
        }
    }

    //***********************************************************************
    //***********************************************************************

    void BaseLoadCondition::GetFirstDerivativesVector(
        Vector& rValues,
        int Step
        )
    {
        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        const unsigned int MatSize = NumberOfNodes * dim;

        if (rValues.size() != MatSize)
        {
            rValues.resize(MatSize, false);
        }

        for (unsigned int i = 0; i < NumberOfNodes; i++)
        {
            const array_1d<double, 3 > & Velocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
            const unsigned int index = i * dim;
            for(unsigned int k = 0; k<dim; ++k)
            {
                rValues[index + k] = Velocity[k];
            }
        }
    }

    //***********************************************************************
    //***********************************************************************

    void BaseLoadCondition::GetSecondDerivativesVector(
        Vector& rValues,
        int Step
        )
    {
        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        const unsigned int MatSize = NumberOfNodes * dim;

        if (rValues.size() != MatSize)
        {
            rValues.resize(MatSize, false);
        }

        for (unsigned int i = 0; i < NumberOfNodes; i++)
        {
            const array_1d<double, 3 > & Acceleration = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
            const unsigned int index = i * dim;
            for(unsigned int k = 0; k < dim; ++k)
            {
                rValues[index + k] = Acceleration[k];
            }
        }
    }

    //************************************************************************************
    //************************************************************************************

    void BaseLoadCondition::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
    {
        // Calculation flags
        const bool CalculateStiffnessMatrixFlag = false;
        const bool CalculateResidualVectorFlag = true;
        MatrixType temp = Matrix();

        CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
    }

    //************************************************************************************
    //************************************************************************************
    void BaseLoadCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
    {
        //calculation flags
        const bool CalculateStiffnessMatrixFlag = true;
        const bool CalculateResidualVectorFlag = true;

        CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
    }

    //***********************************************************************
    //***********************************************************************

    void BaseLoadCondition::CalculateMassMatrix(
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

    void BaseLoadCondition::CalculateDampingMatrix(
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

    void BaseLoadCondition::CalculateAll(
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

    int BaseLoadCondition::Check( const ProcessInfo& rCurrentProcessInfo )
    {
        // Base check
        Condition::Check(rCurrentProcessInfo);
            
        // Verify variable exists
        KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)

        // Check that the condition's nodes contain all required SolutionStepData and Degrees of freedom
        const std::size_t number_of_nodes = this->GetGeometry().size();
        for ( std::size_t i = 0; i < number_of_nodes; i++ ) {
            NodeType &rnode = this->GetGeometry()[i];
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,rnode)

            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, rnode)
            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, rnode)
            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, rnode)
        }

        return 0;
    }

    //***********************************************************************
    //***********************************************************************

    double BaseLoadCondition::GetIntegrationWeight(
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
        const unsigned int PointNumber,
        const double detJ
        )
    {
        return IntegrationPoints[PointNumber].Weight() * detJ;
    }

    //***********************************************************************************
    //***********************************************************************************

    void BaseLoadCondition::AddExplicitContribution(const VectorType& rRHS,
        const Variable<VectorType>& rRHSVariable,
        Variable<array_1d<double,3> >& rDestinationVariable,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const unsigned int number_of_nodes = GetGeometry().PointsNumber();
        const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

        if( rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL )
        {

            for(SizeType i=0; i< number_of_nodes; i++)
            {
                SizeType index = dimension * i;

                array_1d<double, 3 > &ForceResidual = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);
                for(SizeType j=0; j<dimension; j++)
                {
                    #pragma omp atomic
                    ForceResidual[j] += rRHS[index + j];
                }
            }
        }

    KRATOS_CATCH( "" )
    }

} // Namespace Kratos


