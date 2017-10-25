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
        
        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        const bool bHasRotDof = this->HasRotDof();

        unsigned int LocalSize = NumberOfNodes*dim;
        if (bHasRotDof) LocalSize *= 2;

        if (rResult.size() != LocalSize)
        {
            rResult.resize(LocalSize,false);
        }

        const unsigned int pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

        if(dim == 2)
        {
            for (unsigned int i = 0; i < NumberOfNodes; ++i)
            {
                unsigned int index = i * 2;
                if(bHasRotDof) index *= 2;

                rResult[index    ] = GetGeometry()[i].GetDof(DISPLACEMENT_X,pos    ).EquationId();
                rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y,pos + 1).EquationId();

                if (bHasRotDof)
                {
                    rResult[index + 2] = GetGeometry()[i].GetDof(ROTATION_X,pos + 2).EquationId();
                    rResult[index + 3] = GetGeometry()[i].GetDof(ROTATION_Y,pos + 3).EquationId();
                }

            }
        }
        else
        {
            for (unsigned int i = 0; i < NumberOfNodes; ++i)
            {
                unsigned int index = i * 3;
                if(bHasRotDof) index *= 2;

                rResult[index    ] = GetGeometry()[i].GetDof(DISPLACEMENT_X,pos    ).EquationId();
                rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y,pos + 1).EquationId();
                rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z,pos + 2).EquationId();

                if (bHasRotDof)
                {
                    rResult[index + 3] = GetGeometry()[i].GetDof(ROTATION_X,pos + 3).EquationId();
                    rResult[index + 4] = GetGeometry()[i].GetDof(ROTATION_Y,pos + 4).EquationId();
                    rResult[index + 5] = GetGeometry()[i].GetDof(ROTATION_Z,pos + 5).EquationId();
                }

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
        
        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int dim =  GetGeometry().WorkingSpaceDimension();
        const bool bHasRotDof = this->HasRotDof();
        
        unsigned int LocalSize = NumberOfNodes*dim;
        if (bHasRotDof) LocalSize *= 2;

        ElementalDofList.resize(0);
        ElementalDofList.reserve(LocalSize);

        if(dim == 2)
        {
            for (unsigned int i = 0; i < NumberOfNodes; ++i)
            {
                ElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_X));
                ElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_Y));

                if (bHasRotDof)
                {
                    ElementalDofList.push_back( GetGeometry()[i].pGetDof(ROTATION_X));
                    ElementalDofList.push_back( GetGeometry()[i].pGetDof(ROTATION_Y));
                }

            }
        }
        else
        {
            for (unsigned int i = 0; i < NumberOfNodes; ++i)
            {
                ElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_X));
                ElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
                ElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_Z));

                if (bHasRotDof)
                {
                    ElementalDofList.push_back( GetGeometry()[i].pGetDof(ROTATION_X));
                    ElementalDofList.push_back( GetGeometry()[i].pGetDof(ROTATION_Y));
                    ElementalDofList.push_back( GetGeometry()[i].pGetDof(ROTATION_Z));
                }
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
        const bool bHasRotDof = this->HasRotDof();
        
        unsigned int MatSize = NumberOfNodes*dim;
        if (bHasRotDof) MatSize *= 2;

        
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

            if (bHasRotDof)
            {
                const array_1d<double, 3 > & Rotation = GetGeometry()[i].FastGetSolutionStepValue(ROTATION, Step);
                for(unsigned int k = 0; k < dim; ++k)
                {
                    rValues[index + k + dim] = Rotation[k];
                }
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
        const bool bHasRotDof = this->HasRotDof();
        
        unsigned int MatSize = NumberOfNodes*dim;
        if (bHasRotDof) MatSize *= 2;
        
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

            if (bHasRotDof)
            {
                const array_1d<double, 3 > & AngularVelocity = GetGeometry()[i].FastGetSolutionStepValue(ANGULAR_VELOCITY, Step);
                for(unsigned int k = 0; k < dim; ++k)
                {
                    rValues[index + k + dim] = AngularVelocity[k];
                }
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
        const bool bHasRotDof = this->HasRotDof();
        
        unsigned int MatSize = NumberOfNodes*dim;
        if (bHasRotDof) MatSize *= 2;
        
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

            if (bHasRotDof)
            {
                const array_1d<double, 3 > & AngularAcceleration = GetGeometry()[i].FastGetSolutionStepValue(ANGULAR_ACCELERATION, Step);
                for(unsigned int k = 0; k < dim; ++k)
                {
                    rValues[index + k] = AngularAcceleration[k];
                }
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

            if (this->HasRotDof())
            {
                if ( this->GetGeometry()[i].HasDofFor( ROTATION_X ) == false ||
                this->GetGeometry()[i].HasDofFor( ROTATION_Y ) == false ||
                this->GetGeometry()[i].HasDofFor( ROTATION_Z ) == false )
                    {
                     KRATOS_ERROR << "missing one of the dofs for the variable ROTATION on node " << GetGeometry()[i].Id() << " of condition " << Id() << std::endl;
                    }                
            }


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

} // Namespace Kratos


