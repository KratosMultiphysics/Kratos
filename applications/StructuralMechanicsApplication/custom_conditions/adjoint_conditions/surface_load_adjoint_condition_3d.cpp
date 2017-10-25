// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:   Martin Fusseder
//

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_conditions/adjoint_conditions/surface_load_adjoint_condition_3d.h"
#include "utilities/math_utils.h"
#include "utilities/integration_utilities.h"

namespace Kratos
{
    //******************************* CONSTRUCTOR ****************************************
    //************************************************************************************
    
    SurfaceLoadAdjointCondition3D::SurfaceLoadAdjointCondition3D()
    {
    }

    //***********************************************************************************
    //***********************************************************************************

    SurfaceLoadAdjointCondition3D::SurfaceLoadAdjointCondition3D(
        IndexType NewId, 
        GeometryType::Pointer pGeometry
        )
        : SurfaceLoadCondition3D(NewId, pGeometry)
    {
    }

    //***********************************************************************************
    //***********************************************************************************

    SurfaceLoadAdjointCondition3D::SurfaceLoadAdjointCondition3D(
        IndexType NewId, 
        GeometryType::Pointer pGeometry, 
        PropertiesType::Pointer pProperties
        )
        : SurfaceLoadCondition3D(NewId, pGeometry, pProperties)
    {
    }

    //********************************* CREATE *******************************************
    //************************************************************************************

    Condition::Pointer SurfaceLoadAdjointCondition3D::Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const
    {
        return boost::make_shared<SurfaceLoadAdjointCondition3D>(NewId, pGeom, pProperties);
    }

    //***********************************************************************************
    //***********************************************************************************
    
    Condition::Pointer SurfaceLoadAdjointCondition3D::Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const
    {
        return boost::make_shared<SurfaceLoadAdjointCondition3D>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

       //************************************************************************************
    //************************************************************************************

    void SurfaceLoadAdjointCondition3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
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
    void SurfaceLoadAdjointCondition3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo)
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
    
    void SurfaceLoadAdjointCondition3D::GetValuesVector(Vector& rValues, int Step)
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

    void SurfaceLoadAdjointCondition3D::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
				       ProcessInfo& rCurrentProcessInfo) 
    {
        KRATOS_TRY

        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int Dimension = GetGeometry().WorkingSpaceDimension();
        const unsigned int MatSize = NumberOfNodes * Dimension;

        if ( rLeftHandSideMatrix.size1() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );
    
        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize );
        // return zero vector because surface load has no influence on adjoint problem
        KRATOS_CATCH( "" )
    }

    //***********************************************************************
    //***********************************************************************

     //************************************************************************************
    //************************************************************************************

    void SurfaceLoadAdjointCondition3D::CalculateSensitivityMatrix(const Variable<double>& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) 
    {
        KRATOS_TRY


        KRATOS_CATCH( "" )

    }

    //************************************************************************************
    //************************************************************************************

    void SurfaceLoadAdjointCondition3D::CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) 
    {
        KRATOS_TRY

    
        KRATOS_CATCH( "" )
    }

    //***********************************************************************
    //***********************************************************************
    
    int SurfaceLoadAdjointCondition3D::Check( const ProcessInfo& rCurrentProcessInfo )
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


    //******************************* DESTRUCTOR *****************************************
    //************************************************************************************

    SurfaceLoadAdjointCondition3D::~SurfaceLoadAdjointCondition3D()
    {
    }

    //***********************************************************************************
    //***********************************************************************************

} // Namespace Kratos.
