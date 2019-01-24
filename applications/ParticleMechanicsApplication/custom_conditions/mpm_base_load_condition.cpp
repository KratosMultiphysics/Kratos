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
#include "custom_conditions/mpm_base_load_condition.h"

namespace Kratos
{
    void MPMBaseLoadCondition::Initialize()
    {
        // TODO: Add somethig if necessary
    }

    //************************************************************************************
    //************************************************************************************

    void MPMBaseLoadCondition::InitializeSolutionStep( ProcessInfo& CurrentProcessInfo )
    {
        // TODO: Add somethig if necessary
    }

    //************************************************************************************
    //************************************************************************************

    void MPMBaseLoadCondition::InitializeNonLinearIteration( ProcessInfo& CurrentProcessInfo )
    {
        // TODO: Add somethig if necessary
    }

    //************************************************************************************
    //************************************************************************************

    void MPMBaseLoadCondition::FinalizeNonLinearIteration( ProcessInfo& CurrentProcessInfo )
    {
        // TODO: Add somethig if necessary
    }

    //************************************************************************************
    //************************************************************************************

    void MPMBaseLoadCondition::FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo )
    {
        // TODO: Add somethig if necessary
    }

    //************************************************************************************
    //************************************************************************************

    void MPMBaseLoadCondition::EquationIdVector(
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
    void MPMBaseLoadCondition::GetDofList(
        DofsVectorType& ElementalDofList,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_TRY

        GeometryType& rGeom = GetGeometry();
        const unsigned int NumberOfNodes = rGeom.size();
        const unsigned int dim =  rGeom.WorkingSpaceDimension();
        ElementalDofList.resize(0);
        ElementalDofList.reserve(dim * NumberOfNodes);

        if(dim == 2)
        {
            for (unsigned int i = 0; i < NumberOfNodes; ++i)
            {
                ElementalDofList.push_back( rGeom[i].pGetDof(DISPLACEMENT_X));
                ElementalDofList.push_back( rGeom[i].pGetDof(DISPLACEMENT_Y));
            }
        }
        else
        {
            for (unsigned int i = 0; i < NumberOfNodes; ++i)
            {
                ElementalDofList.push_back( rGeom[i].pGetDof(DISPLACEMENT_X));
                ElementalDofList.push_back( rGeom[i].pGetDof(DISPLACEMENT_Y));
                ElementalDofList.push_back( rGeom[i].pGetDof(DISPLACEMENT_Z));
            }
        }
        KRATOS_CATCH("")
    }

    //***********************************************************************
    //***********************************************************************

    void MPMBaseLoadCondition::GetValuesVector(
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

    void MPMBaseLoadCondition::GetFirstDerivativesVector(
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

    void MPMBaseLoadCondition::GetSecondDerivativesVector(
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

    void MPMBaseLoadCondition::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
    {
        // Calculation flags
        const bool CalculateStiffnessMatrixFlag = false;
        const bool CalculateResidualVectorFlag = true;
        MatrixType temp = Matrix();

        CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
    }

    //************************************************************************************
    //************************************************************************************
    void MPMBaseLoadCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
    {
        //calculation flags
        const bool CalculateStiffnessMatrixFlag = true;
        const bool CalculateResidualVectorFlag = true;

        CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
    }

    //***********************************************************************
    //***********************************************************************

    void MPMBaseLoadCondition::CalculateMassMatrix(
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

    void MPMBaseLoadCondition::CalculateDampingMatrix(
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

    void MPMBaseLoadCondition::CalculateAll(
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

    int MPMBaseLoadCondition::Check( const ProcessInfo& rCurrentProcessInfo )
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

    double MPMBaseLoadCondition::GetIntegrationWeight(
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
        const unsigned int PointNumber,
        const double detJ
        )
    {
        return IntegrationPoints[PointNumber].Weight() * detJ;
    }

    //***********************************************************************************
    //***********************************************************************************

    void MPMBaseLoadCondition::AddExplicitContribution(const VectorType& rRHS,
        const Variable<VectorType>& rRHSVariable,
        Variable<array_1d<double,3> >& rDestinationVariable,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        GeometryType& rGeom = GetGeometry();
        const unsigned int number_of_nodes = rGeom.PointsNumber();
        const unsigned int dimension       = rGeom.WorkingSpaceDimension();

        if( rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL )
        {

            for(SizeType i=0; i< number_of_nodes; i++)
            {
                SizeType index = dimension * i;

                if (rGeom[i].SolutionStepsDataHas(FORCE_RESIDUAL))
                {
                    array_1d<double, 3 > &ForceResidual = rGeom[i].FastGetSolutionStepValue(FORCE_RESIDUAL);
                    for(SizeType j=0; j<dimension; j++)
                    {
                        #pragma omp atomic
                        ForceResidual[j] += rRHS[index + j];
                    }
                }
            }
        }

    KRATOS_CATCH( "" )
    }

} // Namespace Kratos


