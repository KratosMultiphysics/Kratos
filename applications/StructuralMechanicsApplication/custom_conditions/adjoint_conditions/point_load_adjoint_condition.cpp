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
#include "custom_conditions/adjoint_conditions/point_load_adjoint_condition.h"
#include "utilities/math_utils.h"
#include "utilities/integration_utilities.h"

namespace Kratos
{
    //******************************* CONSTRUCTOR ****************************************
    //************************************************************************************
    
    PointLoadAdjointCondition::PointLoadAdjointCondition( IndexType NewId, GeometryType::Pointer pGeometry )
        : PointLoadCondition( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************
    
    PointLoadAdjointCondition::PointLoadAdjointCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
        : PointLoadCondition( NewId, pGeometry, pProperties )
    {
    }

    //********************************* CREATE *******************************************
    //************************************************************************************
    
    Condition::Pointer PointLoadAdjointCondition::Create(IndexType NewId,GeometryType::Pointer pGeom,PropertiesType::Pointer pProperties) const
    {
        return boost::make_shared<PointLoadAdjointCondition>(NewId, pGeom, pProperties);
    }

    //************************************************************************************
    //************************************************************************************
    
    Condition::Pointer PointLoadAdjointCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
    {
        return boost::make_shared<PointLoadAdjointCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
    }

    //******************************* DESTRUCTOR *****************************************
    //************************************************************************************
    
    PointLoadAdjointCondition::~PointLoadAdjointCondition()
    {
    }

      //************************************************************************************
    //************************************************************************************

    void PointLoadAdjointCondition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
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
    void PointLoadAdjointCondition::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo)
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
    
    void PointLoadAdjointCondition::GetValuesVector(Vector& rValues, int Step)
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

    void PointLoadAdjointCondition::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
				       ProcessInfo& rCurrentProcessInfo) 
    {
        KRATOS_TRY

        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int Dimension = GetGeometry().WorkingSpaceDimension();
        const unsigned int MatSize = NumberOfNodes * Dimension;

        if ( rLeftHandSideMatrix.size1() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );
    
        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize );

        KRATOS_CATCH( "" )
    }

    //***********************************************************************
    //***********************************************************************

     //************************************************************************************
    //************************************************************************************

    void PointLoadAdjointCondition::CalculateSensitivityMatrix(const Variable<double>& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) 
    {
        KRATOS_TRY

        std::cout << ("I was in CalculateSensitivityMatrix for scalar variables!!!") << std::endl; //----->change this

        KRATOS_CATCH( "" )

    }

    //************************************************************************************
    //************************************************************************************

    void PointLoadAdjointCondition::CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) 
    {
        KRATOS_TRY

        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int Dimension = GetGeometry().WorkingSpaceDimension();
        const unsigned int MatSize = NumberOfNodes * Dimension;

        rOutput = ZeroMatrix(MatSize,MatSize);

        if( rDesignVariable == POINT_LOAD )
        {
            Vector RHS;
            Matrix dummy_LHS;  
            ProcessInfo testProcessInfo = rCurrentProcessInfo;  

            PointLoadCondition::CalculateAll(dummy_LHS, RHS, testProcessInfo, false, true);

            int k = 0;
            for(unsigned int i = 0; i < RHS.size();++i)
            {
                if( abs(RHS[i]) > 1e-12)
                    rOutput(k, i) = 1.0; // or is it sign dependent?  
                k++;
            }
        }
        else if( rDesignVariable == SHAPE_SENSITIVITY )
        {
            rOutput.clear();
        }
        else
            rOutput.clear();
    
        KRATOS_CATCH( "" )
    }

    //***********************************************************************
    //***********************************************************************
    
    int PointLoadAdjointCondition::Check( const ProcessInfo& rCurrentProcessInfo )
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

    std::string PointLoadAdjointCondition::Info() const
    {
        std::string condition_name = "PointLoadAdjointCondition";
        const unsigned int dim =  GetGeometry().WorkingSpaceDimension();
        if(dim == 2)
            condition_name = "PointLoadAdjointCondition2D1N";
        else if(dim == 3) 
            condition_name = "PointLoadAdjointCondition3D1N";

		return condition_name;
    }
    
} // Namespace Kratos


