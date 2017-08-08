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


    /**
     * Calculate the transposed gradient of the condition's residual w.r.t. design variable.
     */
    void PointLoadConditionForSA::CalculateSensitivityMatrix(const Variable<double>& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) 
    {
        KRATOS_TRY

        std::cout << ("I was in CalculateSensitivityMatrix for scalar variables!!!") << std::endl; //----->change this

        KRATOS_CATCH( "" )

    }

    /**
     * Calculate the transposed gradient of the condition's residual w.r.t. design variable.
     */
   void PointLoadConditionForSA::CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) 
    {
        KRATOS_TRY

        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int Dimension = GetGeometry().WorkingSpaceDimension();
        const unsigned int MatSize = NumberOfNodes * Dimension;

        if( rDesignVariable == POINT_LOAD )
        {
            Vector RHS;
            Matrix dummy_LHS;  
            ProcessInfo testProcessInfo = rCurrentProcessInfo;  

            // alternative 1: derive w.r.t. load modulus
            rOutput = ZeroMatrix(1,MatSize);

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
            }

            // alternative 2: treat each load component as induvidual design variable
            /*rOutput = ZeroMatrix(MatSize,MatSize);

            this->CalculateAll(dummy_LHS,RHS,testProcessInfo,false,true);

            double norm2_load = norm_2(RHS);

            int k = 0;
            for(unsigned int i = 0; i < RHS.size();++i)
            {
                if( abs(RHS[i]) > 1e-12)
                {
                    rOutput(k++, i) = 1.0; // or is it sign dependent?
                   
                }
                k++;
            }*/

        }
        else if( rDesignVariable == ADJOINT_NODE_COORD )
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
    
} // Namespace Kratos


