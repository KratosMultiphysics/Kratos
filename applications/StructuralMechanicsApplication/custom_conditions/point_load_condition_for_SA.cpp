// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder
//

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_conditions/point_load_condition_for_SA.h"
#include "utilities/math_utils.h"
#include "utilities/integration_utilities.h"

namespace Kratos
{
    //******************************* CONSTRUCTOR ****************************************
    //************************************************************************************
    
    PointLoadConditionForSA::PointLoadConditionForSA( IndexType NewId, GeometryType::Pointer pGeometry )
        : BaseLoadCondition( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************
    
    PointLoadConditionForSA::PointLoadConditionForSA( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
        : BaseLoadCondition( NewId, pGeometry, pProperties )
    {
    }

    //********************************* CREATE *******************************************
    //************************************************************************************
    
    Condition::Pointer PointLoadConditionForSA::Create(IndexType NewId,GeometryType::Pointer pGeom,PropertiesType::Pointer pProperties) const
    {
        return boost::make_shared<PointLoadConditionForSA>(NewId, pGeom, pProperties);
    }

    //************************************************************************************
    //************************************************************************************
    
    Condition::Pointer PointLoadConditionForSA::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
    {
        return boost::make_shared<PointLoadConditionForSA>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
    }

    //******************************* DESTRUCTOR *****************************************
    //************************************************************************************
    
    PointLoadConditionForSA::~PointLoadConditionForSA()
    {
    }

    //************************************************************************************
    //************************************************************************************

    void PointLoadConditionForSA::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY
  
        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        if (rResult.size() != dim * NumberOfNodes)
        {
            rResult.resize(dim*NumberOfNodes,false);
        }

        const unsigned int pos = this->GetGeometry()[0].GetDofPosition(ADJOINT_DISPLACEMENT_X);

        if(dim == 2)
        {
            for (unsigned int i = 0; i < NumberOfNodes; ++i)
            {
                const unsigned int index = i * 2;
                rResult[index    ] = GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_X,pos    ).EquationId();
                rResult[index + 1] = GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_Y,pos + 1).EquationId();
            }
        }
        else
        {
            for (unsigned int i = 0; i < NumberOfNodes; ++i)
            {
                const unsigned int index = i * 3;
                rResult[index    ] = GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_X,pos    ).EquationId();
                rResult[index + 1] = GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_Y,pos + 1).EquationId();
                rResult[index + 2] = GetGeometry()[i].GetDof(ADJOINT_DISPLACEMENT_Z,pos + 2).EquationId();
            }
        }
        KRATOS_CATCH("")
    }
    
    //***********************************************************************
    //***********************************************************************
    void PointLoadConditionForSA::GetDofList(
        DofsVectorType& ElementalDofList,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_TRY

        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int dim =  GetGeometry().WorkingSpaceDimension();
        ElementalDofList.resize(0);
        ElementalDofList.reserve(dim * NumberOfNodes);

        if(dim == 2)
        {
            for (unsigned int i = 0; i < NumberOfNodes; ++i)
            {
                ElementalDofList.push_back( GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_X));
                ElementalDofList.push_back( GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_Y));
            }
        }
        else
        {
            for (unsigned int i = 0; i < NumberOfNodes; ++i)
            {
                ElementalDofList.push_back( GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_X));
                ElementalDofList.push_back( GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_Y));
                ElementalDofList.push_back( GetGeometry()[i].pGetDof(ADJOINT_DISPLACEMENT_Z));
            }
        }

        KRATOS_CATCH("")
    }
    
    //***********************************************************************
    //***********************************************************************
    
    void PointLoadConditionForSA::GetValuesVector(Vector& rValues, int Step)
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
            const array_1d<double, 3 > & Displacement = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_DISPLACEMENT, Step);
            unsigned int index = i * dim;
            for(unsigned int k = 0; k < dim; ++k)
            {
                rValues[index + k] = Displacement[k];
            }
        }
    }

    //************************************************************************************
    //************************************************************************************

    void PointLoadConditionForSA::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
				       ProcessInfo& rCurrentProcessInfo) // remove this maybe to base load condition
    {
        KRATOS_TRY

        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int Dimension = GetGeometry().WorkingSpaceDimension();

        // Resizing as needed the LHS
        const unsigned int MatSize = NumberOfNodes * Dimension;

    
        if ( rLeftHandSideMatrix.size1() != MatSize )
        {
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );
        }

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize );

        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************

    void PointLoadConditionForSA::CalculateAll( 
        MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag 
        )
    {
        KRATOS_TRY
        
        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int Dimension = GetGeometry().WorkingSpaceDimension();

        // Resizing as needed the LHS
        const unsigned int MatSize = NumberOfNodes * Dimension;

        if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
        {
            if ( rLeftHandSideMatrix.size1() != MatSize )
            {
                rLeftHandSideMatrix.resize( MatSize, MatSize, false );
            }

            noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
        }

        //resizing as needed the RHS
        if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
        {
            if ( rRightHandSideVector.size( ) != MatSize )
            {
                rRightHandSideVector.resize( MatSize, false );
            }

            noalias( rRightHandSideVector ) = ZeroVector( MatSize ); //resetting RHS
        }

        // Vector with a loading applied to the condition
        array_1d<double, 3 > PointLoad = ZeroVector(3);
        if( this->Has( POINT_LOAD ) )
        {
            noalias(PointLoad) = this->GetValue( POINT_LOAD );
        }

        for (unsigned int ii = 0; ii < NumberOfNodes; ++ii)
        {
            const unsigned int base = ii*Dimension;
            
            if( GetGeometry()[ii].SolutionStepsDataHas( POINT_LOAD ) )
            {
                noalias(PointLoad) += GetGeometry()[ii].FastGetSolutionStepValue( POINT_LOAD );
            }
            for(unsigned int k = 0; k < Dimension; ++k)
            {
                rRightHandSideVector[base + k] += GetPointLoadIntegrationWeight() * PointLoad[k];
            }
        }

        KRATOS_CATCH( "" )
    }
    
    //************************************************************************************
    //************************************************************************************
    
    double PointLoadConditionForSA::GetPointLoadIntegrationWeight()
    {
        return 1.0;
    }

    //************************************************************************************
    //************************************************************************************

    void PointLoadConditionForSA::CalculateSensitivityMatrix(const Variable<double>& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) 
    {
        KRATOS_TRY

        std::cout << ("I was in CalculateSensitivityMatrix for scalar variables!!!") << std::endl; //----->change this

        KRATOS_CATCH( "" )

    }

    //************************************************************************************
    //************************************************************************************

    void PointLoadConditionForSA::CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) 
    {
        KRATOS_TRY

        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int Dimension = GetGeometry().WorkingSpaceDimension();
        const unsigned int MatSize = NumberOfNodes * Dimension;

        if( rDesignVariable == POINT_LOAD)
        {
            Vector RHS;
            Matrix dummy_LHS;  
            ProcessInfo testProcessInfo = rCurrentProcessInfo;  
            
            // alternative 1: derive w.r.t. load modulus
            /*rOutput = ZeroMatrix(1,MatSize);

            this->CalculateAll(dummy_LHS,RHS,testProcessInfo,false,true);

            double norm2_load = norm_2(RHS);

            for(unsigned int i = 0; i < RHS.size();++i)
            {
                if( abs(RHS[i]) > 1e-12)
                {
                    RHS[i] /= norm2_load;
                }
                rOutput(0, i) = RHS[i];
               // std::cout << ("Pseudo load = ") << RHS[i]  << std::endl;
            }*/

            // alternative 2: treat each load component as induvidual design variable
            rOutput = ZeroMatrix(MatSize,MatSize);

            this->CalculateAll(dummy_LHS,RHS,testProcessInfo,false,true);

            int k = 0;
            for(unsigned int i = 0; i < RHS.size();++i)
            {
                if( abs(RHS[i]) > 1e-12)
                {
                    rOutput(k, i) = 1.0; // or is it sign dependent?
                    
                }
                k++;
            }

        }
        else if( rDesignVariable == SHAPE_SENSITIVITY )
        {
            rOutput.resize(MatSize , MatSize);

            for (unsigned int i = 0; i < MatSize; ++i)
                for (unsigned int j = 0; j < MatSize; ++j)
                    rOutput(i,j) = 0.0;

        }
        else
        {
            KRATOS_THROW_ERROR(std::invalid_argument, "The chosen design variable is not provided by this condition!", "");
        }

        KRATOS_CATCH( "" )
    }

     //***********************************************************************
    //***********************************************************************
    
    int PointLoadConditionForSA::Check( const ProcessInfo& rCurrentProcessInfo )
    {
        if ( ADJOINT_DISPLACEMENT.Key() == 0 )
        {
            KRATOS_ERROR <<  "DISPLACEMENT has Key zero! (check if the application is correctly registered" << std::endl;
        }
         if ( DISPLACEMENT.Key() == 0 )
        {
            KRATOS_ERROR <<  "DISPLACEMENT has Key zero! (check if the application is correctly registered" << std::endl;
        }

        //verify that the dofs exist
        for ( unsigned int i = 0; i < this->GetGeometry().size(); i++ )
        {
            if ( this->GetGeometry()[i].SolutionStepsDataHas( ADJOINT_DISPLACEMENT ) == false )
            {
                KRATOS_ERROR << "missing variable DISPLACEMENT on node " << this->GetGeometry()[i].Id() << std::endl;
            }

            if ( this->GetGeometry()[i].HasDofFor( ADJOINT_DISPLACEMENT_X ) == false ||
                 this->GetGeometry()[i].HasDofFor( ADJOINT_DISPLACEMENT_Y ) == false ||
                 this->GetGeometry()[i].HasDofFor( ADJOINT_DISPLACEMENT_Z ) == false )
            {
                KRATOS_ERROR << "missing one of the dofs for the variable ADJOINT_DISPLACEMENT on node " << GetGeometry()[i].Id() << " of condition " << Id() << std::endl;
            }
        }
        
        return 0;
    }
    
} // Namespace Kratos


