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
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        const unsigned int strain_size = GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();

        Matrix B( strain_size, number_of_nodes * dim );
        Matrix F( dim, dim );
        Matrix D( strain_size, strain_size );
        Matrix C( dim, dim );
        Vector strain_vector( strain_size );
        Vector stress_vector( strain_size );
        Matrix DN_DX( number_of_nodes, dim );
        Matrix J0(dim,dim), InvJ0(dim,dim);

        // Resizing as needed the LHS
        const unsigned int mat_size = number_of_nodes * dim;

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

        const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients(  );

        const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues(  );

        // Calculating actual jacobian
        GeometryType::JacobiansType J;

        GetGeometry().Jacobian( J );

        // Auxiliary terms
        Vector body_force = ZeroVector(3);
        
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRAIN, false);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        
        // if strain has to be computed inside of the constitutive law with PK2
        //rValues.SetDeformationGradientF(rVariables.F); //in this case F is the whole deformation gradient

        Values.SetStrainVector(strain_vector); //this is the input  parameter
        Values.SetStressVector(stress_vector); //this is the output parameter

        for ( unsigned int point_number = 0; point_number < integration_points.size(); point_number++ )
        {
            const double detJ0 = CalculateDerivativesOnReference(J0, InvJ0, DN_DX, point_number, GetGeometry().GetDefaultIntegrationMethod());
            
            // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
            noalias( DN_DX ) = prod( DN_De[point_number], InvJ0 );

            // Deformation gradient
            noalias( F ) = prod( J[point_number], InvJ0 );
            
            // Axisymmetric case
            if (strain_size == 4)
            {
                F.resize(3, 3); // We keep the old values
                for (unsigned int index = 0; index < 1; index++)
                {
                    F(index, 2) = 0.0;
                    F(2, index) = 0.0;
                }
                Vector N;
                N = GetGeometry().ShapeFunctionsValues( N, integration_points[point_number].Coordinates() );
                const double current_radius = StructuralMechanicsMathUtilities::CalculateRadius(N, GetGeometry(), Current);
                const double initial_radius = StructuralMechanicsMathUtilities::CalculateRadius(N, GetGeometry(), Initial);
                F(2, 2) = current_radius/initial_radius;
            }
            
            // Here we essentially set the input parameters
            const double detF = MathUtils<double>::Det(F);
            Values.SetDeterminantF(detF); //assuming the determinant is computed somewhere else
            Values.SetDeformationGradientF(F); //F computed somewhere else
            
            // Here we set the space on which the results shall be written
            Values.SetConstitutiveMatrix(D); //assuming the determinant is computed somewhere else
            Values.SetStressVector(stress_vector); //F computed somewhere else
            
            // Actually do the computations in the ConstitutiveLaw    
            mConstitutiveLawVector[point_number]->CalculateMaterialResponsePK2(Values); //here the calculations are actually done 

            // Calculating operator B
            CalculateB( B, F, DN_DX, strain_vector.size(), integration_points, point_number );

            // Calculating weights for integration on the reference configuration
            double IntToReferenceWeight = GetIntegrationWeight(integration_points, point_number, detJ0); 

            if ( dim == 2 && GetProperties().Has( THICKNESS )) 
            {
                IntToReferenceWeight *= GetProperties()[THICKNESS];
            }

            if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
            {
                // Contributions to stiffness matrix calculated on the reference config
                /* Material stiffness matrix */
                typedef Matrix temp_type;
                noalias( rLeftHandSideMatrix ) += IntToReferenceWeight *prod( trans( B ), temp_type(prod(D, B)));
                /* Geometric stiffness matrix */
                CalculateAndAddKg( rLeftHandSideMatrix, DN_DX, stress_vector, IntToReferenceWeight );
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
                CalculateAndAdd_ExtForceContribution( row( Ncontainer, point_number ), rCurrentProcessInfo, body_force, rRightHandSideVector, IntToReferenceWeight );

                // Operation performed: rRightHandSideVector -= IntForce*IntToReferenceWeight
                noalias( rRightHandSideVector ) -= IntToReferenceWeight * prod( trans( B ), stress_vector );
            }
        }

        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************

    void TotalLagrangian::CalculateOnIntegrationPoints( 
        const Variable<Matrix >& rVariable, 
        std::vector< Matrix >& Output, 
        const ProcessInfo& rCurrentProcessInfo
        ) 
    { 
        KRATOS_TRY 
 
        const unsigned int number_of_nodes = GetGeometry().size(); 
        const unsigned int dim = GetGeometry().WorkingSpaceDimension(); 
        const unsigned int strain_size = GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();
 
        Matrix F( dim, dim ); 
        Matrix D( strain_size, strain_size ); 
        Vector strain_vector( strain_size ); 
        Vector stress_vector( strain_size ); 
        Matrix DN_DX( number_of_nodes, dim ); 
 
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo); 
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRAIN, false); 
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS); 
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false); 
        Values.SetStrainVector(strain_vector); //this is the input  parameter 
        Values.SetStressVector(stress_vector); //this is the output parameter 
        Matrix J0(dim,dim), InvJ0(dim,dim); 
         
        //reading integration points and local gradients 
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(  ); 
 
//         const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues(  ); 
 
        //calculating actual jacobian 
        GeometryType::JacobiansType J; 
 
        J = GetGeometry().Jacobian( J ); 
 
        if ( Output.size() != integration_points.size() ) 
        {
            Output.resize( integration_points.size() ); 
        }
 
        for ( unsigned int point_number = 0; point_number < integration_points.size(); point_number++ ) 
        { 
            CalculateDerivativesOnReference(J0, InvJ0, DN_DX, point_number, GetGeometry().GetDefaultIntegrationMethod()); 
 
            // Deformation gradient 
            noalias( F ) = prod( J[point_number], InvJ0 ); 
 
            Matrix Plasticstrain_vector( GetGeometry().size(), GetGeometry().WorkingSpaceDimension() ); 
 
            if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR ) 
            { 
                // Strain calculation // TODO: MOVE THIS TO THE CL!!!!!!
                Matrix C = prod( trans( F ), F ); 
                CalculateStrain( C, strain_vector ); 
                if ( Output[point_number].size2() != strain_vector.size() )
                {
                    Output[point_number].resize( 1, strain_vector.size(), false ); 
                }
 
                for ( unsigned int ii = 0; ii < strain_vector.size(); ii++ )
                {
                    Output[point_number]( 0, ii ) = strain_vector[ii]; 
                }
            } 
            else if ( rVariable == PK2_STRESS_TENSOR ) 
            { 
                if ( Output[point_number].size2() != stress_vector.size() ) 
                {
                    Output[point_number].resize( 1, stress_vector.size(), false ); 
                }
                                         
                //here we essentially set the input parameters 
                const double detF = MathUtils<double>::Det(F); 
                Values.SetDeterminantF(detF); //assuming the determinant is computed somewhere else 
                Values.SetDeformationGradientF(F); //F computed somewhere else 
                 
                //actually do the computations in the ConstitutiveLaw     
                mConstitutiveLawVector[point_number]->CalculateMaterialResponsePK2(Values); //here the calculations are actually done  
 
                for ( unsigned int ii = 0; ii < stress_vector.size(); ii++ ) 
                { 
                    Output[point_number]( 0, ii ) = stress_vector[ii]; 
                } 
            } 
//             else if ( rVariable == GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR ) 
//             { 
//                 double size = strain_vector.size(); 
//                 Matrix Plasticstrain_vector( 1, size ); 
//  
//                 if ( Output[point_number].size2() != strain_vector.size() ) 
//                     Output[point_number].resize( 1, size, false ); 
//  
//                 mConstitutiveLawVector[point_number]->GetValue( GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR, Plasticstrain_vector ); 
//  
//                 Output[point_number] = Plasticstrain_vector; 
//             } 
        } 
 
        KRATOS_CATCH( "" ) 
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

    void TotalLagrangian::CalculateStrain(
        const Matrix& C,
        Vector& rStrainTensor 
        )
    {
        KRATOS_TRY
        
        const unsigned int strain_size = GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();
        
        if ( strain_size == 3 )
        {
            if ( rStrainTensor.size() != 3 ) rStrainTensor.resize( 3, false );

            rStrainTensor[0] = 0.5 * ( C( 0, 0 ) - 1.00 );
            rStrainTensor[1] = 0.5 * ( C( 1, 1 ) - 1.00 );
            rStrainTensor[2] = C( 0, 1 );
        }
        else if ( strain_size == 4 )
        {
            rStrainTensor[0] = 0.5 * ( C( 0, 0 ) - 1.00 );
            rStrainTensor[1] = 0.5 * ( C( 1, 1 ) - 1.00 );
            rStrainTensor[2] = 0.5 * ( C( 2, 2 ) - 1.00 );
            rStrainTensor[3] = C( 0, 1 );
        }
        else if ( strain_size == 6 )
        {
            if ( rStrainTensor.size() != 6 ) rStrainTensor.resize( 6, false );

            rStrainTensor[0] = 0.5 * ( C( 0, 0 ) - 1.00 );
            rStrainTensor[1] = 0.5 * ( C( 1, 1 ) - 1.00 );
            rStrainTensor[2] = 0.5 * ( C( 2, 2 ) - 1.00 );
            rStrainTensor[3] = C( 0, 1 ); // xy
            rStrainTensor[4] = C( 1, 2 ); // yz
            rStrainTensor[5] = C( 0, 2 ); // xz
        }

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


