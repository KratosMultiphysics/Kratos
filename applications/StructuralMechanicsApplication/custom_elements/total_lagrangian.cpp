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
#include "includes/define.h"
#include "custom_elements/total_lagrangian.h"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/structural_mechanics_math_utilities.hpp"

namespace Kratos
{

    TotalLagrangian::TotalLagrangian( IndexType NewId, GeometryType::Pointer pGeometry )
            : BaseSolidElement( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************

    TotalLagrangian::TotalLagrangian( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
            : BaseSolidElement( NewId, pGeometry, pProperties )
    {
    }

    Element::Pointer TotalLagrangian::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
    {
        return Element::Pointer( new TotalLagrangian( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
    }

    TotalLagrangian::~TotalLagrangian()
    {
    }

    //************************************************************************************
    //************************************************************************************

    void TotalLagrangian::CalculateAll( 
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

        // Reading integration points
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(  );

        // Auxiliary terms
        Vector body_force = ZeroVector(3);
        
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);
        
        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRAIN, false);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        
        // If strain has to be computed inside of the constitutive law with PK2
        Values.SetStrainVector(this_kinematic_variables.StrainVector); //this is the input  parameter

        for ( unsigned int point_number = 0; point_number < integration_points.size(); point_number++ )
        {
            UpdateKinematics(this_kinematic_variables, Values, point_number, integration_points);

            // Calculating weights for integration on the reference configuration
            double IntToReferenceWeight = GetIntegrationWeight(integration_points, point_number, this_kinematic_variables.detJ0); 

            if ( dimension == 2 && GetProperties().Has( THICKNESS )) 
            {
                IntToReferenceWeight *= GetProperties()[THICKNESS];
            }

            if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
            {
                // Contributions to stiffness matrix calculated on the reference config
                /* Material stiffness matrix */
                typedef Matrix temp_type;
                noalias( rLeftHandSideMatrix ) += IntToReferenceWeight *prod( trans( this_kinematic_variables.B ), temp_type(prod(this_kinematic_variables.D, this_kinematic_variables.B)));
                /* Geometric stiffness matrix */
                CalculateAndAddKg( rLeftHandSideMatrix, this_kinematic_variables.DN_DX, this_kinematic_variables.StressVector, IntToReferenceWeight );
            }

            if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
            {
                // Contribution to external forces
                if (GetProperties().Has( VOLUME_ACCELERATION ) == true)
                {
                    body_force += GetProperties()[VOLUME_ACCELERATION];
                }
                if( GetGeometry()[0].SolutionStepsDataHas(VOLUME_ACCELERATION) )
                {
                    body_force += GetGeometry()[0].FastGetSolutionStepValue(VOLUME_ACCELERATION);
                }

                // Operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
                CalculateAndAdd_ExtForceContribution( this_kinematic_variables.N, rCurrentProcessInfo, body_force, rRightHandSideVector, IntToReferenceWeight );

                // Operation performed: rRightHandSideVector -= IntForce*IntToReferenceWeight
                noalias( rRightHandSideVector ) -= IntToReferenceWeight * prod( trans( this_kinematic_variables.B ), this_kinematic_variables.StressVector );
            }
        }

        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************
    
    void TotalLagrangian::UpdateKinematics(
        KinematicVariables& rThisKinematicVariables, 
        ConstitutiveLaw::Parameters& rValues,
        const unsigned int PointNumber,
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
        const Vector Displacements
        )
    {
        // Shape functions
        rThisKinematicVariables.N = GetGeometry().ShapeFunctionsValues(rThisKinematicVariables.N, IntegrationPoints[PointNumber].Coordinates());

        // Calculating actual jacobian
        Matrix J;
        J = GetGeometry().Jacobian( J, PointNumber );
    
        rThisKinematicVariables.detJ0 = CalculateDerivativesOnReference(rThisKinematicVariables.J0, rThisKinematicVariables.InvJ0, rThisKinematicVariables.DN_DX, PointNumber, GetGeometry().GetDefaultIntegrationMethod());

        // Deformation gradient
        noalias( rThisKinematicVariables.F ) = prod( J, rThisKinematicVariables.InvJ0 );
        
        // Axisymmetric case
        if (rThisKinematicVariables.StrainVector.size() == 4)
        {
            rThisKinematicVariables.F.resize(3, 3); // We keep the old values
            for (unsigned int index = 0; index < 1; index++)
            {
                rThisKinematicVariables.F(index, 2) = 0.0;
                rThisKinematicVariables.F(2, index) = 0.0;
            }

            rThisKinematicVariables.N = GetGeometry().ShapeFunctionsValues( rThisKinematicVariables.N, IntegrationPoints[PointNumber].Coordinates() );
            const double current_radius = StructuralMechanicsMathUtilities::CalculateRadius(rThisKinematicVariables.N, GetGeometry(), Current);
            const double initial_radius = StructuralMechanicsMathUtilities::CalculateRadius(rThisKinematicVariables.N, GetGeometry(), Initial);
            rThisKinematicVariables.F(2, 2) = current_radius/initial_radius;
        }
        
        // Here we essentially set the input parameters
        rThisKinematicVariables.detF = MathUtils<double>::Det(rThisKinematicVariables.F);
        rValues.SetDeterminantF(rThisKinematicVariables.detF); //assuming the determinant is computed somewhere else
        rValues.SetDeformationGradientF(rThisKinematicVariables.F); //F computed somewhere else
        
        // Here we set the space on which the results shall be written
        rValues.SetConstitutiveMatrix(rThisKinematicVariables.D); //assuming the determinant is computed somewhere else
        rValues.SetStressVector(rThisKinematicVariables.StressVector); //F computed somewhere else
        
        // Actually do the computations in the ConstitutiveLaw    
        mConstitutiveLawVector[PointNumber]->CalculateMaterialResponsePK2(rValues); //here the calculations are actually done 

        // Calculating operator B
        CalculateB( rThisKinematicVariables.B, rThisKinematicVariables.F, rThisKinematicVariables.DN_DX, rThisKinematicVariables.StrainVector.size(), IntegrationPoints, PointNumber );
    }

    //************************************************************************************
    //************************************************************************************

    inline void TotalLagrangian::CalculateAndAdd_ExtForceContribution(
        const Vector& N,
        const ProcessInfo& CurrentProcessInfo,
        Vector& body_force,
        VectorType& rRightHandSideVector,
        double weight
    )
    {
        KRATOS_TRY
        
        const unsigned int number_of_nodes = GetGeometry().PointsNumber();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            const unsigned int index = dimension * i;

            for ( unsigned int j = 0; j < dimension; j++ ) 
            {
                rRightHandSideVector[index + j] += weight * N[i] * body_force[j];
            }
        }

        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************

    void TotalLagrangian::CalculateAndAddKg(
        MatrixType& K,
        Matrix& DN_DX,
        Vector& stress_vector,
        double weight )
    {
        KRATOS_TRY
        unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        Matrix StressTensor = MathUtils<double>::StressVectorToTensor( stress_vector );
        Matrix ReducedKg = prod( DN_DX, weight * Matrix( prod( StressTensor, trans( DN_DX ) ) ) ); //to be optimized
        MathUtils<double>::ExpandAndAddReducedMatrix( K, ReducedKg, dimension );

        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************

    void TotalLagrangian::CalculateB(
        Matrix& B,
        Matrix& F,
        Matrix& DN_DX,
        unsigned int strain_size,
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
        const unsigned int PointNumber
        )
    {
        KRATOS_TRY
        
        const unsigned int number_of_nodes = GetGeometry().PointsNumber();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        // For axisymmetric case
        Vector N;
        double Radius = 0.0f;
        
        if ( strain_size == 4 )
        {
            N = GetGeometry().ShapeFunctionsValues( N, IntegrationPoints[PointNumber].Coordinates() );
            Radius = StructuralMechanicsMathUtilities::CalculateRadius(N, GetGeometry());
        }
        
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            const unsigned int index = dimension * i;

            if ( strain_size == 3 )
            {
                B( 0, index + 0 ) = F( 0, 0 ) * DN_DX( i, 0 );
                B( 0, index + 1 ) = F( 1, 0 ) * DN_DX( i, 0 );
                B( 1, index + 0 ) = F( 0, 1 ) * DN_DX( i, 1 );
                B( 1, index + 1 ) = F( 1, 1 ) * DN_DX( i, 1 );
                B( 2, index + 0 ) = F( 0, 0 ) * DN_DX( i, 1 ) + F( 0, 1 ) * DN_DX( i, 0 );
                B( 2, index + 1 ) = F( 1, 0 ) * DN_DX( i, 1 ) + F( 1, 1 ) * DN_DX( i, 0 );
            }
            else if ( strain_size == 4 )
            {
                B( 0, index + 0 ) = F( 0, 0 ) * DN_DX( i, 0 );
                B( 0, index + 1 ) = F( 1, 0 ) * DN_DX( i, 0 );
                B( 1, index + 1 ) = F( 0, 1 ) * DN_DX( i, 1 );
                B( 1, index + 1 ) = F( 1, 1 ) * DN_DX( i, 1 );
                B( 2, index + 0 ) = N[i]/Radius;
                B( 3, index + 0 ) = F( 0, 0 ) * DN_DX( i, 1 ) + F( 0, 1 ) * DN_DX( i, 0 );
                B( 3, index + 1 ) = F( 1, 0 ) * DN_DX( i, 1 ) + F( 1, 1 ) * DN_DX( i, 0 );
            }
            else
            {
                B( 0, index + 0 ) = F( 0, 0 ) * DN_DX( i, 0 );
                B( 0, index + 1 ) = F( 1, 0 ) * DN_DX( i, 0 );
                B( 0, index + 2 ) = F( 2, 0 ) * DN_DX( i, 0 );
                B( 1, index + 0 ) = F( 0, 1 ) * DN_DX( i, 1 );
                B( 1, index + 1 ) = F( 1, 1 ) * DN_DX( i, 1 );
                B( 1, index + 2 ) = F( 2, 1 ) * DN_DX( i, 1 );
                B( 2, index + 0 ) = F( 0, 2 ) * DN_DX( i, 2 );
                B( 2, index + 1 ) = F( 1, 2 ) * DN_DX( i, 2 );
                B( 2, index + 2 ) = F( 2, 2 ) * DN_DX( i, 2 );
                B( 3, index + 0 ) = F( 0, 0 ) * DN_DX( i, 1 ) + F( 0, 1 ) * DN_DX( i, 0 );
                B( 3, index + 1 ) = F( 1, 0 ) * DN_DX( i, 1 ) + F( 1, 1 ) * DN_DX( i, 0 );
                B( 3, index + 2 ) = F( 2, 0 ) * DN_DX( i, 1 ) + F( 2, 1 ) * DN_DX( i, 0 );
                B( 4, index + 0 ) = F( 0, 1 ) * DN_DX( i, 2 ) + F( 0, 2 ) * DN_DX( i, 1 );
                B( 4, index + 1 ) = F( 1, 1 ) * DN_DX( i, 2 ) + F( 1, 2 ) * DN_DX( i, 1 );
                B( 4, index + 2 ) = F( 2, 1 ) * DN_DX( i, 2 ) + F( 2, 2 ) * DN_DX( i, 1 );
                B( 5, index + 0 ) = F( 0, 2 ) * DN_DX( i, 0 ) + F( 0, 0 ) * DN_DX( i, 2 );
                B( 5, index + 1 ) = F( 1, 2 ) * DN_DX( i, 0 ) + F( 1, 0 ) * DN_DX( i, 2 );
                B( 5, index + 2 ) = F( 2, 2 ) * DN_DX( i, 0 ) + F( 2, 0 ) * DN_DX( i, 2 );
            }
        }

        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************

    int  TotalLagrangian::Check( const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        const unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();

        // Verify that the variables are correctly initialized

        if ( VELOCITY.Key() == 0 )
        {
            KRATOS_ERROR << "VELOCITY has Key zero! (check if the application is correctly registered"<< std::endl;
        }

        if ( DISPLACEMENT.Key() == 0 )
        {
            KRATOS_ERROR << "DISPLACEMENT has Key zero! (check if the application is correctly registered"<< std::endl;
        }

        if ( ACCELERATION.Key() == 0 )
        {
            KRATOS_ERROR << "ACCELERATION has Key zero! (check if the application is correctly registered"<< std::endl;
        }

        if ( DENSITY.Key() == 0 )
        {
            KRATOS_ERROR << "DENSITY has Key zero! (check if the application is correctly registered"<< std::endl;
        }

        if ( VOLUME_ACCELERATION.Key() == 0 )
        {
            KRATOS_ERROR << "VOLUME_ACCELERATION has Key zero! (check if the application is correctly registered"<< std::endl;
        }

        if ( THICKNESS.Key() == 0 )
        {
            KRATOS_ERROR << "THICKNESS has Key zero! (check if the application is correctly registered"<< std::endl;
        }

        //verify that the dofs exist
        for ( unsigned int i = 0; i < this->GetGeometry().size(); i++ )
        {
            if ( this->GetGeometry()[i].SolutionStepsDataHas( DISPLACEMENT ) == false )
            {
                KRATOS_ERROR << "Missing variable DISPLACEMENT on node " << this->GetGeometry()[i].Id() << std::endl;
            }

            if ( this->GetGeometry()[i].HasDofFor( DISPLACEMENT_X ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Y ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Z ) == false )
            {
                KRATOS_ERROR << "missing one of the dofs for the variable DISPLACEMENT on node " << GetGeometry()[i].Id() << std::endl;
            }
        }

        //verify that the constitutive law exists
        if ( this->GetProperties().Has( CONSTITUTIVE_LAW ) == false )
        {
            KRATOS_ERROR << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;
        }

        //verify that the constitutive law has the correct dimension
        if ( dimension == 2 )
        {
//             if ( this->GetProperties().Has( THICKNESS ) == false ) // NOTE: Not mandatory
//             {
//                 KRATOS_ERROR << "THICKNESS not provided for element " << this->Id() << std::endl;
//             }

            if ( this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize() < 3 || this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize() > 4)
            {
                KRATOS_ERROR << "Wrong constitutive law used. This is a 2D element! expected strain size is 3 or 4 (el id = ) " << this->Id() << std::endl;
            }
        }
        else
        {
            if ( this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize() != 6 )
            {
                KRATOS_ERROR << "Wrong constitutive law used. This is a 3D element! expected strain size is 6 (el id = ) "<<  this->Id() << std::endl;
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

    void TotalLagrangian::save( Serializer& rSerializer ) const
    {
        rSerializer.save( "Name", "TotalLagrangian" );
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseSolidElement );
    }
    
    //************************************************************************************
    //************************************************************************************
    
    void TotalLagrangian::load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseSolidElement );
    }

} // Namespace Kratos


