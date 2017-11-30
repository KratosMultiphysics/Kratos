// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Philipp Bucher
//

// System includes


// External includes


// Project includes
#include "custom_conditions/point_moment_condition.h"

namespace Kratos
{
    //******************************* CONSTRUCTOR ****************************************
    //************************************************************************************
    
    PointMomentCondition::PointMomentCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : BaseLoadCondition( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************

    PointMomentCondition::PointMomentCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
        : BaseLoadCondition( NewId, pGeometry, pProperties )
    {
    }

    //********************************* CREATE *******************************************
    //************************************************************************************

    Condition::Pointer PointMomentCondition::Create(IndexType NewId,GeometryType::Pointer pGeom,PropertiesType::Pointer pProperties) const
    {
        return boost::make_shared<PointMomentCondition>(NewId, pGeom, pProperties);
    }

    //************************************************************************************
    //************************************************************************************

    Condition::Pointer PointMomentCondition::Create( IndexType NewId, NodesArrayType const& rThisNodes,  PropertiesType::Pointer pProperties ) const
    {
        return boost::make_shared<PointMomentCondition>( NewId, GetGeometry().Create( rThisNodes ), pProperties );
    }

    //******************************* DESTRUCTOR *****************************************
    //************************************************************************************

    PointMomentCondition::~PointMomentCondition()
    {
    }
    
    //************************************************************************************
    //************************************************************************************


    void PointMomentCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY
        
        const unsigned int num_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        if (rResult.size() != dim * num_nodes) rResult.resize(dim*num_nodes,false);

        const unsigned int pos = this->GetGeometry()[0].GetDofPosition(ROTATION_X);

        if(dim == 2)
        {
            for (unsigned int i = 0; i < num_nodes; ++i)
            {
                const unsigned int index = i * 2;
                rResult[index    ] = GetGeometry()[i].GetDof(ROTATION_X,pos    ).EquationId();
                rResult[index + 1] = GetGeometry()[i].GetDof(ROTATION_Y,pos + 1).EquationId();
            }
        }
        else
        {
            for (unsigned int i = 0; i < num_nodes; ++i)
            {
                const unsigned int index = i * 3;
                rResult[index    ] = GetGeometry()[i].GetDof(ROTATION_X,pos    ).EquationId();
                rResult[index + 1] = GetGeometry()[i].GetDof(ROTATION_Y,pos + 1).EquationId();
                rResult[index + 2] = GetGeometry()[i].GetDof(ROTATION_Z,pos + 2).EquationId();
            }
        }
        KRATOS_CATCH("")
    }
    
    //***********************************************************************
    //***********************************************************************
    void PointMomentCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_TRY
        
        const unsigned int num_nodes = GetGeometry().size();
        const unsigned int dim =  GetGeometry().WorkingSpaceDimension();
        rElementalDofList.resize(0);
        rElementalDofList.reserve(dim * num_nodes);

        if(dim == 2)
        {
            for (unsigned int i = 0; i < num_nodes; ++i)
            {
                rElementalDofList.push_back( GetGeometry()[i].pGetDof(ROTATION_X));
                rElementalDofList.push_back( GetGeometry()[i].pGetDof(ROTATION_Y));
            }
        }
        else
        {
            for (unsigned int i = 0; i < num_nodes; ++i)
            {
                rElementalDofList.push_back( GetGeometry()[i].pGetDof(ROTATION_X));
                rElementalDofList.push_back( GetGeometry()[i].pGetDof(ROTATION_Y));
                rElementalDofList.push_back( GetGeometry()[i].pGetDof(ROTATION_Z));
            }
        }
        KRATOS_CATCH("")
    }
    
    //***********************************************************************
    //***********************************************************************
    
    void PointMomentCondition::GetValuesVector(
        Vector& rValues,
        int Step
        )
    {
        const unsigned int num_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        const unsigned int mat_size = num_nodes * dim;
        
        if (rValues.size() != mat_size)
        {
            rValues.resize(mat_size, false);
        }
        
        for (unsigned int i = 0; i < num_nodes; i++)
        {
            const array_1d<double, 3 > & r_rotation = GetGeometry()[i].FastGetSolutionStepValue(ROTATION, Step);
            unsigned int index = i * dim;
            for(unsigned int k = 0; k < dim; ++k)
            {
                rValues[index + k] = r_rotation[k];
            }
        }
    }
    
    //***********************************************************************
    //***********************************************************************
    
    void PointMomentCondition::GetFirstDerivativesVector(
        Vector& rValues,
        int Step 
        )
    {
        const unsigned int num_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        const unsigned int mat_size = num_nodes * dim;
        
        if (rValues.size() != mat_size)
        {
            rValues.resize(mat_size, false);
        }
        
        for (unsigned int i = 0; i < num_nodes; i++)
        {
            const array_1d<double, 3 > & r_angular_vel = GetGeometry()[i].FastGetSolutionStepValue(ANGULAR_VELOCITY, Step);
            const unsigned int index = i * dim;
            for(unsigned int k = 0; k<dim; ++k)
            {
                rValues[index + k] = r_angular_vel[k];
            }
        }
    }
    
    //***********************************************************************
    //***********************************************************************
    
    void PointMomentCondition::GetSecondDerivativesVector(
        Vector& rValues,
        int Step 
        )
    {
        const unsigned int num_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        const unsigned int mat_size = num_nodes * dim;
        
        if (rValues.size() != mat_size)
        {
            rValues.resize(mat_size, false);
        }
        
        for (unsigned int i = 0; i < num_nodes; i++)
        {
            const array_1d<double, 3 > & r_angular_acc = GetGeometry()[i].FastGetSolutionStepValue(ANGULAR_ACCELERATION, Step);
            const unsigned int index = i * dim;
            for(unsigned int k = 0; k < dim; ++k)
            {
                rValues[index + k] = r_angular_acc[k];
            }
        }
    }
    
    //***********************************************************************
    //***********************************************************************

    void PointMomentCondition::CalculateAll( 
        MatrixType& rLeftHandSideMatrix, 
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag 
        )
    {
        KRATOS_TRY
        
        const unsigned int num_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();

        // Resizing as needed the LHS
        const unsigned int mat_size = num_nodes * dim;

        if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
        {
            if ( rLeftHandSideMatrix.size1() != mat_size )
            {
                rLeftHandSideMatrix.resize( mat_size, mat_size, false );
            }

            noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS
        }

        //resizing as needed the RHS
        if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
        {
            if ( rRightHandSideVector.size( ) != mat_size )
            {
                rRightHandSideVector.resize( mat_size, false );
            }

            noalias( rRightHandSideVector ) = ZeroVector( mat_size ); //resetting RHS
        }

        // Vector with a loading applied to the condition
        array_1d<double, 3 > point_moment = ZeroVector(3);
        if( this->Has( POINT_MOMENT ) )
        {
            noalias(point_moment) = this->GetValue( POINT_MOMENT );
        }

        for (unsigned int ii = 0; ii < num_nodes; ++ii)
        {
            const unsigned int base = ii*dim;
            
            if( GetGeometry()[ii].SolutionStepsDataHas( POINT_MOMENT ) )
            {
                noalias(point_moment) += GetGeometry()[ii].FastGetSolutionStepValue( POINT_MOMENT );
            }
            
            for(unsigned int k = 0; k < dim; ++k)
            {
                rRightHandSideVector[base + k] += GetPointMomentIntegrationWeight() * point_moment[k];
            }
        }

        KRATOS_CATCH( "" )
    }
    
    //************************************************************************************
    //************************************************************************************
    
    double PointMomentCondition::GetPointMomentIntegrationWeight()
    {
        return 1.0;
    }
    
    //***********************************************************************
    //***********************************************************************
    
    int PointMomentCondition::Check( const ProcessInfo& rCurrentProcessInfo )
    {
        if ( ROTATION.Key() == 0 )
        {
            KRATOS_ERROR <<  "ROTATION has Key zero! (check if the application is correctly registered" << std::endl;
        }

        //verify that the dofs exist
        for ( unsigned int i = 0; i < this->GetGeometry().size(); i++ )
        {
            if ( this->GetGeometry()[i].SolutionStepsDataHas( ROTATION ) == false )
            {
                KRATOS_ERROR << "missing variable ROTATION on node " << this->GetGeometry()[i].Id() << std::endl;
            }

            if ( this->GetGeometry()[i].HasDofFor( ROTATION_X ) == false ||
                 this->GetGeometry()[i].HasDofFor( ROTATION_Y ) == false ||
                 this->GetGeometry()[i].HasDofFor( ROTATION_Z ) == false )
            {
                KRATOS_ERROR << "missing one of the dofs for the variable ROTATION on node " << GetGeometry()[i].Id() << " of condition " << Id() << std::endl;
            }
        }
        
        return 0;
    }

} // Namespace Kratos


