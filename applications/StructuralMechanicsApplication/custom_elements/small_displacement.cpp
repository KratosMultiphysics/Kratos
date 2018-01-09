// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrándiz
//

// System includes

// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/small_displacement.h"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
    SmallDisplacement::SmallDisplacement( IndexType NewId, GeometryType::Pointer pGeometry )
            : BaseSolidElement( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************

    SmallDisplacement::SmallDisplacement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
            : BaseSolidElement( NewId, pGeometry, pProperties )
    {
    }

    Element::Pointer SmallDisplacement::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
    {
        return Element::Pointer( new SmallDisplacement( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
    }

    SmallDisplacement::~SmallDisplacement()
    {
    }

    //************************************************************************************
    //************************************************************************************

    void SmallDisplacement::CalculateAll( 
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag 
        )
    {
        KRATOS_TRY;
        
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        const unsigned int strain_size = GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();

        KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
        ConstitutiveVariables this_constitutive_variables(strain_size);

        // Resizing as needed the LHS
        const unsigned int mat_size = number_of_nodes * dimension;

        if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
        {
            if ( rLeftHandSideMatrix.size1() != mat_size )
            {
                rLeftHandSideMatrix.resize( mat_size, mat_size, false );
            }

            noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS
        }

        // Resizing as needed the RHS
        if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
        {
            if ( rRightHandSideVector.size() != mat_size )
            {
                rRightHandSideVector.resize( mat_size, false );
            }

            rRightHandSideVector = ZeroVector( mat_size ); //resetting RHS
        }

        // Reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(  );
        
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        
        // If strain has to be computed inside of the constitutive law with PK2
        Values.SetStrainVector(this_constitutive_variables.StrainVector); //this is the input  parameter

        // Displacements vector
        Vector displacements;
        GetValuesVector(displacements);
        
        for ( unsigned int point_number = 0; point_number < integration_points.size(); point_number++ )
        {
            // Contribution to external forces
            const Vector body_force = this->GetBodyForce(integration_points, point_number);
            
            // Compute element kinematics B, F, DN_DX ...
            CalculateKinematicVariables(this_kinematic_variables, point_number, integration_points);
            
            // Compute material reponse
            CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, GetStressMeasure(), displacements);
            
            // Calculating weights for integration on the reference configuration
            double int_to_reference_weight = GetIntegrationWeight(integration_points, point_number, this_kinematic_variables.detJ0); 

            if ( dimension == 2 && GetProperties().Has( THICKNESS )) 
            {
                int_to_reference_weight *= GetProperties()[THICKNESS];
            }

            if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
            {
                // Contributions to stiffness matrix calculated on the reference config
                this->CalculateAndAddKm( rLeftHandSideMatrix, this_kinematic_variables.B, this_constitutive_variables.D, int_to_reference_weight );
            }

            if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
            {
                this->CalculateAndAddResidualVector(rRightHandSideVector, this_kinematic_variables, rCurrentProcessInfo, body_force, this_constitutive_variables.StressVector, int_to_reference_weight);
            }
        }
        
        KRATOS_CATCH( "" )
    }
    
    //************************************************************************************
    //************************************************************************************
    
    void SmallDisplacement::CalculateKinematicVariables(
        KinematicVariables& rThisKinematicVariables, 
        const unsigned int PointNumber,
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints
        )
    {        
        // Shape functions
        rThisKinematicVariables.N = GetGeometry().ShapeFunctionsValues(rThisKinematicVariables.N, IntegrationPoints[PointNumber].Coordinates());
        
        rThisKinematicVariables.detJ0 = CalculateDerivativesOnReferenceConfiguration(rThisKinematicVariables.J0, rThisKinematicVariables.InvJ0, rThisKinematicVariables.DN_DX, PointNumber, GetGeometry().GetDefaultIntegrationMethod()); 
        
        if (rThisKinematicVariables.detJ0 < 0.0)
        {
            KRATOS_ERROR << "WARNING:: ELEMENT ID: " << this->Id() << " INVERTED. DETJ0: " << rThisKinematicVariables.detJ0 << std::endl;
        }
        
        // Compute B
        CalculateB( rThisKinematicVariables.B, rThisKinematicVariables.DN_DX, IntegrationPoints, PointNumber );
    }
    
    //************************************************************************************
    //************************************************************************************
    
    void SmallDisplacement::CalculateConstitutiveVariables(
        KinematicVariables& rThisKinematicVariables, 
        ConstitutiveVariables& rThisConstitutiveVariables, 
        ConstitutiveLaw::Parameters& rValues,
        const unsigned int PointNumber,
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure,
        const Vector Displacements
        )
    {        
        // Compute strain
        noalias(rThisConstitutiveVariables.StrainVector) = prod(rThisKinematicVariables.B, Displacements);

        // Compute equivalent F
        rThisKinematicVariables.F = ComputeEquivalentF(rThisConstitutiveVariables.StrainVector);

        // Here we essentially set the input parameters
        rThisKinematicVariables.detF = MathUtils<double>::Det(rThisKinematicVariables.F);
        rValues.SetDeterminantF(rThisKinematicVariables.detF); //assuming the determinant is computed somewhere else
        rValues.SetDeformationGradientF(rThisKinematicVariables.F); //F computed somewhere else
        
        // Here we set the space on which the results shall be written
        rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.D); //assuming the determinant is computed somewhere else
        rValues.SetStressVector(rThisConstitutiveVariables.StressVector); //F computed somewhere else
        
        // Actually do the computations in the ConstitutiveLaw    
        mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(rValues, ThisStressMeasure); //here the calculations are actually done 
    }

    //************************************************************************************
    //************************************************************************************

    void SmallDisplacement::CalculateB(
        Matrix& rB,
        const Matrix& DN_DX,
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
        const unsigned int PointNumber
        )
    {
        KRATOS_TRY;
        
        const unsigned int number_of_nodes = GetGeometry().PointsNumber();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        rB.clear();
        
        if(dimension == 2)
        {
            for ( unsigned int i = 0; i < number_of_nodes; ++i )
            {
                rB( 0, i*2     ) = DN_DX( i, 0 );
                rB( 1, i*2 + 1 ) = DN_DX( i, 1 );
                rB( 2, i*2     ) = DN_DX( i, 1 );
                rB( 2, i*2 + 1 ) = DN_DX( i, 0 );
            }
        }
        else if(dimension == 3)
        {
            for ( unsigned int i = 0; i < number_of_nodes; ++i )
            {
                rB( 0, i*3     ) = DN_DX( i, 0 );
                rB( 1, i*3 + 1 ) = DN_DX( i, 1 );
                rB( 2, i*3 + 2 ) = DN_DX( i, 2 );
                rB( 3, i*3     ) = DN_DX( i, 1 );
                rB( 3, i*3 + 1 ) = DN_DX( i, 0 );
                rB( 4, i*3 + 1 ) = DN_DX( i, 2 );
                rB( 4, i*3 + 2 ) = DN_DX( i, 1 );
                rB( 5, i*3     ) = DN_DX( i, 2 );
                rB( 5, i*3 + 2 ) = DN_DX( i, 0 );
            }
        }

        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************

    Matrix SmallDisplacement::ComputeEquivalentF(const Vector& rStrainTensor)
    {
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        Matrix F(dim,dim);
        
        if(dim == 2)
        {
            F(0,0) = 1.0+rStrainTensor(0);  
            F(0,1) = 0.5*rStrainTensor(2); 
            F(1,0) = 0.5*rStrainTensor(2);   
            F(1,1) = 1.0+rStrainTensor(1);
        }
        else
        {
            F(0,0) = 1.0+rStrainTensor(0);     
            F(0,1) = 0.5*rStrainTensor(3); 
            F(0,2) = 0.5*rStrainTensor(5);
            F(1,0) = 0.5*rStrainTensor(3);   
            F(1,1) = 1.0+rStrainTensor(1);   
            F(1,2) = 0.5*rStrainTensor(4);
            F(2,0) = 0.5*rStrainTensor(5);   
            F(2,1) = 0.5*rStrainTensor(4); 
            F(2,2) = 1.0+rStrainTensor(2);
        }
        
        return F;
    }

    //************************************************************************************
    //************************************************************************************
    
    int  SmallDisplacement::Check( const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        const unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();

        // Verify that the variables are correctly initialized
        if ( VELOCITY.Key() == 0 )
        {
            KRATOS_ERROR <<"VELOCITY has Key zero! (check if the application is correctly registered" << std::endl;
        }

        if ( DISPLACEMENT.Key() == 0 )
        {
            KRATOS_ERROR <<"DISPLACEMENT has Key zero! (check if the application is correctly registered" << std::endl;
        }

        if ( ACCELERATION.Key() == 0 )
        {
            KRATOS_ERROR <<"ACCELERATION has Key zero! (check if the application is correctly registered" << std::endl;
        }

        if ( DENSITY.Key() == 0 )
        {
            KRATOS_ERROR <<"DENSITY has Key zero! (check if the application is correctly registered" << std::endl;
        }

        if ( VOLUME_ACCELERATION.Key() == 0 )
        {
            KRATOS_ERROR <<"VOLUME_ACCELERATION has Key zero! (check if the application is correctly registered" << std::endl;
        }

        if ( THICKNESS.Key() == 0 )
        {
            KRATOS_ERROR <<"THICKNESS has Key zero! (check if the application is correctly registered" << std::endl;
        }

        // Verify that the dofs exist
        for ( unsigned int i = 0; i < this->GetGeometry().size(); i++ )
        {
            if ( this->GetGeometry()[i].SolutionStepsDataHas( DISPLACEMENT ) == false )
            {
                KRATOS_ERROR <<"Missing variable DISPLACEMENT on node " << this->GetGeometry()[i].Id() << std::endl;
            }

            if ( this->GetGeometry()[i].HasDofFor( DISPLACEMENT_X ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Y ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Z ) == false )
            {
                KRATOS_ERROR <<"Missing one of the dofs for the variable DISPLACEMENT on node " << GetGeometry()[i].Id() << std::endl;
            }
        }

        // Verify that the constitutive law exists
        if ( this->GetProperties().Has( CONSTITUTIVE_LAW ) == false )
        {
            KRATOS_ERROR << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;
        }


        // Verify that the constitutive law has the correct dimension
        if ( dimension == 2 )
        {
//             if ( this->GetProperties().Has( THICKNESS ) == false ) // NOTE: Not mandatory
//             {
//                 KRATOS_ERROR << "THICKNESS not provided for element " << this->Id() << std::endl;
//             }

            if ( this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize() < 3 || this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize() > 4 )
            {
                KRATOS_ERROR << "Wrong constitutive law used. This is a 2D element! expected strain size is 3 or 4 (el id = ) " << this->Id() << std::endl;
            }
        }
        else
        {
            if ( this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize() != 6 )
            {
                KRATOS_ERROR << "Wrong constitutive law used. This is a 3D element! expected strain size is 6 (el id = ) " << this->Id() << std::endl;
            }
        }

        //check constitutive law
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        {
            return mConstitutiveLawVector[i]->Check( GetProperties(), GetGeometry(), rCurrentProcessInfo );
        }

        //check if it is in the XY plane for 2D case

        return 0;

        KRATOS_CATCH( "" );
    }

    //************************************************************************************
    //************************************************************************************
    
    void SmallDisplacement::save( Serializer& rSerializer ) const
    {
        rSerializer.save( "Name", "SmallDisplacement" );
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseSolidElement );
    }
    
    //************************************************************************************
    //************************************************************************************
    
    void SmallDisplacement::load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseSolidElement );
    }

} // Namespace Kratos


