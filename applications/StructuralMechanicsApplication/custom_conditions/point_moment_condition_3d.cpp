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
#include "custom_conditions/point_moment_condition_3d.h"

namespace Kratos
{
    //******************************* CONSTRUCTOR ****************************************
    //************************************************************************************
    
    PointMomentCondition3D::PointMomentCondition3D( IndexType NewId, GeometryType::Pointer pGeometry )
    : BaseLoadCondition( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************

    PointMomentCondition3D::PointMomentCondition3D( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
        : BaseLoadCondition( NewId, pGeometry, pProperties )
    {
    }

    //********************************* CREATE *******************************************
    //************************************************************************************

    Condition::Pointer PointMomentCondition3D::Create(IndexType NewId,GeometryType::Pointer pGeom,PropertiesType::Pointer pProperties) const
    {
        return boost::make_shared<PointMomentCondition3D>(NewId, pGeom, pProperties);
    }

    //************************************************************************************
    //************************************************************************************

    Condition::Pointer PointMomentCondition3D::Create( IndexType NewId, NodesArrayType const& rThisNodes,  PropertiesType::Pointer pProperties ) const
    {
        return boost::make_shared<PointMomentCondition3D>( NewId, GetGeometry().Create( rThisNodes ), pProperties );
    }

    //******************************* DESTRUCTOR *****************************************
    //************************************************************************************

    PointMomentCondition3D::~PointMomentCondition3D()
    {
    }
    
    //************************************************************************************
    //************************************************************************************


    void PointMomentCondition3D::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY
        
        if (rResult.size() != 3) rResult.resize(3,false);
        rResult[0] = GetGeometry()[0].GetDof(ROTATION_X).EquationId();
        rResult[1] = GetGeometry()[0].GetDof(ROTATION_Y).EquationId();
        rResult[2] = GetGeometry()[0].GetDof(ROTATION_Z).EquationId();
        
        KRATOS_CATCH("")
    }
    
    //***********************************************************************
    //***********************************************************************
    void PointMomentCondition3D::GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_TRY
        
        if (rElementalDofList.size() != 3) rElementalDofList.resize(3,false);
        rElementalDofList[0] = ( GetGeometry()[0].pGetDof(ROTATION_X));
        rElementalDofList[1] = ( GetGeometry()[0].pGetDof(ROTATION_Y));
        rElementalDofList[2] = ( GetGeometry()[0].pGetDof(ROTATION_Z));
        
        KRATOS_CATCH("")
    }
    
    //***********************************************************************
    //***********************************************************************
    
    void PointMomentCondition3D::GetValuesVector(
        Vector& rValues,
        int Step
        )
    {
        const array_1d<double, 3 > & r_rotation = GetGeometry()[0].FastGetSolutionStepValue(ROTATION, Step);

        if (rValues.size() != 3) rValues.resize(3, false);
        rValues[0] = r_rotation[0];
        rValues[1] = r_rotation[1];
        rValues[2] = r_rotation[2];
    }
    
    //***********************************************************************
    //***********************************************************************
    
    void PointMomentCondition3D::GetFirstDerivativesVector(
        Vector& rValues,
        int Step 
        )
    {
        const array_1d<double, 3 > & r_angular_vel = GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY, Step);

        if (rValues.size() != 3) rValues.resize(3, false);
        rValues[0] = r_angular_vel[0];
        rValues[1] = r_angular_vel[1];
        rValues[2] = r_angular_vel[2];
    }
    
    //***********************************************************************
    //***********************************************************************
    
    void PointMomentCondition3D::GetSecondDerivativesVector(
        Vector& rValues,
        int Step 
        )
    {
        const array_1d<double, 3 > & r_angular_acc = GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_ACCELERATION, Step);

        if (rValues.size() != 3) rValues.resize(3, false);
        rValues[0] = r_angular_acc[0];
        rValues[1] = r_angular_acc[1];
        rValues[2] = r_angular_acc[2];
    }
    
    //***********************************************************************
    //***********************************************************************

    void PointMomentCondition3D::CalculateAll( 
        MatrixType& rLeftHandSideMatrix, 
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag 
        )
    {
        KRATOS_TRY
        
        const unsigned int dim = 3;

        // Resizing as needed the LHS
        const unsigned int mat_size = dim;

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
        
        if( GetGeometry()[0].SolutionStepsDataHas( POINT_MOMENT ) )
        {
            noalias(point_moment) += GetGeometry()[0].FastGetSolutionStepValue( POINT_MOMENT );
        }
        
        for(unsigned int k = 0; k < dim; ++k)
        {
            rRightHandSideVector[k] += GetPointMomentIntegrationWeight() * point_moment[k];
        }

        KRATOS_CATCH( "" )
    }
    
    //************************************************************************************
    //************************************************************************************
    
    double PointMomentCondition3D::GetPointMomentIntegrationWeight()
    {
        return 1.0;
    }
    
    //***********************************************************************
    //***********************************************************************
    
    int PointMomentCondition3D::Check( const ProcessInfo& rCurrentProcessInfo )
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


