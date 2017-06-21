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
#include "custom_elements/kinematic_linear.h"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
    KinematicLinear::KinematicLinear( IndexType NewId, GeometryType::Pointer pGeometry )
            : BaseSolidElement( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************

    KinematicLinear::KinematicLinear( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
            : BaseSolidElement( NewId, pGeometry, pProperties )
    {
    }

    Element::Pointer KinematicLinear::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
    {
        return Element::Pointer( new KinematicLinear( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
    }

    KinematicLinear::~KinematicLinear()
    {
    }

    //************************************************************************************
    //************************************************************************************

    void KinematicLinear::CalculateAll( 
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
        Matrix D( strain_size, strain_size );
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

        const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues(  );

        // Auxiliary terms
        Vector body_force = ZeroVector(3);
        
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRAIN, false);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        
        // If strain has to be computed inside of the constitutive law with PK2
        //rValues.SetDeformationGradientF(rVariables.F); //in this case F is the whole deformation gradient

        Values.SetStrainVector(strain_vector); //this is the input  parameter
        Values.SetStressVector(stress_vector); //this is the output parameter
        
        Vector displacements;
        GetValuesVector(displacements);

        for ( unsigned int point_number = 0; point_number < integration_points.size(); point_number++ )
        {
            const double detJ0 = CalculateDerivativesOnReference(J0, InvJ0, DN_DX, point_number, GetGeometry().GetDefaultIntegrationMethod()); 
            
            //Compute B and strain
            CalculateB( B, DN_DX, integration_points, point_number );
            noalias(strain_vector) = prod(B,displacements);

            // Compute equivalent F
            Matrix F = ComputeEquivalentF(strain_vector);

            // Here we essentially set the input parameters
            const double detF = MathUtils<double>::Det(F);
            Values.SetDeterminantF(detF); //assuming the determinant is computed somewhere else
            Values.SetDeformationGradientF(F); //F computed somewhere else
            
            // Here we set the space on which the results shall be written
            Values.SetConstitutiveMatrix(D); //assuming the determinant is computed somewhere else
            Values.SetStressVector(stress_vector); //F computed somewhere else
            
            // Actually do the computations in the ConstitutiveLaw    
            mConstitutiveLawVector[point_number]->CalculateMaterialResponsePK2(Values); //here the calculations are actually done 

            // Calculating weights for integration on the reference configuration
            double int_to_reference_weight = GetIntegrationWeight(integration_points, point_number, detJ0); 

            if ( dim == 2 && GetProperties().Has( THICKNESS )) 
            {
                int_to_reference_weight *= GetProperties()[THICKNESS];
            }

            if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
            {
                // Contributions to stiffness matrix calculated on the reference config
                typedef Matrix temp_type;
                noalias( rLeftHandSideMatrix ) += int_to_reference_weight * prod( trans( B ), temp_type(prod(D, B)));
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

                // Operation performed: rRightHandSideVector += ExtForce*int_to_reference_weight
                CalculateAndAdd_ExtForceContribution( row( Ncontainer, point_number ), rCurrentProcessInfo, body_force, rRightHandSideVector, int_to_reference_weight );

                // Operation performed: rRightHandSideVector -= IntForce*int_to_reference_weight
                noalias( rRightHandSideVector ) -= int_to_reference_weight * prod( trans( B ), stress_vector );
            }
        }
        
        KRATOS_CATCH( "" )
    }
    
    //************************************************************************************
    //************************************************************************************

    void KinematicLinear::CalculateOnIntegrationPoints( 
        const Variable<Matrix >& rVariable, 
        std::vector< Matrix >& rOutput, 
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
        Matrix B( strain_size, number_of_nodes * dim );
        Matrix J0(dim,dim), InvJ0(dim,dim);

        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRAIN, false);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
        Values.SetStrainVector(strain_vector); //this is the input  parameter
        Values.SetStressVector(stress_vector); //this is the output parameter
        
        //reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(  );

        if ( rOutput.size() != integration_points.size() )
        {
            rOutput.resize( integration_points.size() );
        }
        
        Vector displacements;
        GetValuesVector(displacements);

        for ( unsigned int point_number = 0; point_number < integration_points.size(); point_number++ )
        {
            CalculateDerivativesOnReference(J0, InvJ0, DN_DX, point_number, GetGeometry().GetDefaultIntegrationMethod()); 
            
            //Compute B and strain // TODO: MOVE THIS TO THE CL!!!!!!
            CalculateB( B, DN_DX, integration_points, point_number );
            noalias(strain_vector) = prod(B,displacements);
            F = ComputeEquivalentF(strain_vector);
            
            if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR ) 
            {
                if ( rOutput[point_number].size2() != strain_vector.size() )
                {
                    rOutput[point_number].resize( 1, strain_vector.size(), false );
                }

                for ( unsigned int ii = 0; ii < strain_vector.size(); ii++ )
                {
                    rOutput[point_number]( 0, ii ) = strain_vector[ii];
                }
            }
            else if ( rVariable == PK2_STRESS_TENSOR )
            {
                if ( rOutput[point_number].size2() != stress_vector.size() )
                {
                    rOutput[point_number].resize( 1, stress_vector.size(), false );
                }
                                        
                // Here we essentially set the input parameters
                const double detF = MathUtils<double>::Det(F); 
                Values.SetDeterminantF(detF); //assuming the determinant is computed somewhere else
                Values.SetDeformationGradientF(F); //F computed somewhere else
                
                //actually do the computations in the ConstitutiveLaw    
                mConstitutiveLawVector[point_number]->CalculateMaterialResponsePK2(Values); //here the calculations are actually done 

                for ( unsigned int ii = 0; ii < stress_vector.size(); ii++ )
                {
                    rOutput[point_number]( 0, ii ) = stress_vector[ii];
                }
            }
//             else if ( rVariable == GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR )
//             {
//                 double size = strain_vector.size();
//                 Matrix Plasticstrain_vector( 1, size );
// 
//                 if ( rOutput[point_number].size2() != strain_vector.size() )
//                     rOutput[point_number].resize( 1, size, false );
// 
//                 mConstitutiveLawVector[point_number]->GetValue( GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR, Plasticstrain_vector );
// 
//                 rOutput[point_number] = Plasticstrain_vector;
//             }
        }

        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************

    inline void KinematicLinear::CalculateAndAdd_ExtForceContribution(
        const Vector& N,
        const ProcessInfo& CurrentProcessInfo,
        Vector& body_force,
        VectorType& rRightHandSideVector,
        const double weight
        )
    {
        KRATOS_TRY;
        
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

    void KinematicLinear::CalculateB(
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
                rB( 0, i*2 ) = DN_DX( i, 0 );
                rB( 1, i*2 + 1 ) = DN_DX( i, 1 );
                rB( 2, i*2 ) = DN_DX( i, 1 );
                rB( 2, i*2 + 1 ) = DN_DX( i, 0 );
            }
        }
        else if(dimension == 3)
        {
            for ( unsigned int i = 0; i < number_of_nodes; ++i )
            {
                rB( 0, i*3 ) = DN_DX( i, 0 );
                rB( 1, i*3 + 1 ) = DN_DX( i, 1 );
                rB( 2, i*3 + 2 ) = DN_DX( i, 2 );
                rB( 3, i*3 ) = DN_DX( i, 1 );
                rB( 3, i*3 + 1 ) = DN_DX( i, 0 );
                rB( 4, i*3 + 1 ) = DN_DX( i, 2 );
                rB( 4, i*3 + 2 ) = DN_DX( i, 1 );
                rB( 5, i*3 ) = DN_DX( i, 2 );
                rB( 5, i*3 + 2 ) = DN_DX( i, 0 );
            }
        }

        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************

    Matrix KinematicLinear::ComputeEquivalentF(const Vector& rStrainTensor)
    {
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        Matrix F(dim,dim);
        
        if(dim == 2)
        {
            F(0,0) = 1+rStrainTensor(0);  
            F(0,1) = 0.5*rStrainTensor(2); 
            F(1,0) = 0.5*rStrainTensor(2);   
            F(1,1) = 1+rStrainTensor(1);
        }
        else
        {
            F(0,0) = 1+rStrainTensor(0);     
            F(0,1) = 0.5*rStrainTensor(3); 
            F(0,2) = 0.5*rStrainTensor(5);
            F(1,0) = 0.5*rStrainTensor(3);   
            F(1,1) = 1+rStrainTensor(1);   
            F(1,2) = 0.5*rStrainTensor(4);
            F(2,0) = 0.5*rStrainTensor(5);   
            F(2,1) = 0.5*rStrainTensor(4); 
            F(2,2) = 1+rStrainTensor(2);
        }
        
        return F;
    }

    //************************************************************************************
    //************************************************************************************
    
    int  KinematicLinear::Check( const ProcessInfo& rCurrentProcessInfo )
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
    
    void KinematicLinear::save( Serializer& rSerializer ) const
    {
        rSerializer.save( "Name", "KinematicLinear" );
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseSolidElement );
    }
    
    //************************************************************************************
    //************************************************************************************
    
    void KinematicLinear::load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseSolidElement );
    }

} // Namespace Kratos


