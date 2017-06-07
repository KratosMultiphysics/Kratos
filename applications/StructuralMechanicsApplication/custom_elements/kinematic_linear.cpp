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

    void KinematicLinear::Initialize()
    {
        KRATOS_TRY

        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = GetGeometry().IntegrationPoints(  );

        //Constitutive Law initialisation
        if ( mConstitutiveLawVector.size() != IntegrationPoints.size() )
        {
            mConstitutiveLawVector.resize( IntegrationPoints.size() );
        }

        InitializeMaterial();

        KRATOS_CATCH( "" )
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
        
        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        const unsigned int StrainSize = ( dim == 2 ) ? 3 : 6;

        Matrix B( StrainSize, NumberOfNodes * dim );
        Matrix D( StrainSize, StrainSize );
        Vector StrainVector( StrainSize );
        Vector StressVector( StrainSize );
        Matrix DN_DX( NumberOfNodes, dim );
        Matrix J0(dim,dim), InvJ0(dim,dim);

        // Resizing as needed the LHS
        const unsigned int MatSize = NumberOfNodes * dim;

        if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
        {
            if ( rLeftHandSideMatrix.size1() != MatSize )
            {
                rLeftHandSideMatrix.resize( MatSize, MatSize, false );
            }

            noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
        }

        // Resizing as needed the RHS
        if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
        {
            if ( rRightHandSideVector.size() != MatSize )
            {
                rRightHandSideVector.resize( MatSize, false );
            }

            rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS
        }

        // Reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = GetGeometry().IntegrationPoints(  );

        const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues(  );

        // Auxiliary terms
        Vector BodyForce = ZeroVector(3);
        
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRAIN, false);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        
        // If strain has to be computed inside of the constitutive law with PK2
        //rValues.SetDeformationGradientF(rVariables.F); //in this case F is the whole deformation gradient

        Values.SetStrainVector(StrainVector); //this is the input  parameter
        Values.SetStressVector(StressVector); //this is the output parameter
        
        Vector displacements;
        GetValuesVector(displacements);

        for ( unsigned int PointNumber = 0; PointNumber < IntegrationPoints.size(); PointNumber++ )
        {
            const Matrix& DN_De = GetGeometry().ShapeFunctionsLocalGradients()[PointNumber];
            double detJ0;
            CalculateDerivativesOnReference(J0, InvJ0, DN_DX, detJ0, DN_De);
            
            //Compute B and strain
            CalculateB( B, DN_DX );
            noalias(StrainVector) = prod(B,displacements);

            // Compute equivalent F
            Matrix F = ComputeEquivalentF(StrainVector);

            // Here we essentially set the input parameters
            const double detF = MathUtils<double>::Det(F);
            Values.SetDeterminantF(detF); //assuming the determinant is computed somewhere else
            Values.SetDeformationGradientF(F); //F computed somewhere else
            
            // Here we set the space on which the results shall be written
            Values.SetConstitutiveMatrix(D); //assuming the determinant is computed somewhere else
            Values.SetStressVector(StressVector); //F computed somewhere else
            
            // Actually do the computations in the ConstitutiveLaw    
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponsePK2(Values); //here the calculations are actually done 

            // Calculating weights for integration on the reference configuration
            double IntToReferenceWeight = GetIntegrationWeight(IntegrationPoints, PointNumber, detJ0); 

            if ( dim == 2 && GetProperties().Has( THICKNESS )) 
            {
                IntToReferenceWeight *= GetProperties()[THICKNESS];
            }

            if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
            {
                // Contributions to stiffness matrix calculated on the reference config
                typedef Matrix temp_type;
                noalias( rLeftHandSideMatrix ) += IntToReferenceWeight * prod( trans( B ), temp_type(prod(D, B)));
            }

            if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
            {
                // Contribution to external forces
                if (GetProperties().Has( VOLUME_ACCELERATION ) == true)
                {
                    BodyForce += GetProperties()[VOLUME_ACCELERATION];
                }

                // Operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
                CalculateAndAdd_ExtForceContribution( row( Ncontainer, PointNumber ), rCurrentProcessInfo, BodyForce, rRightHandSideVector, IntToReferenceWeight );

                // Operation performed: rRightHandSideVector -= IntForce*IntToReferenceWeight
                noalias( rRightHandSideVector ) -= IntToReferenceWeight * prod( trans( B ), StressVector );
            }
        }
        
        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************

    void KinematicLinear::InitializeSolutionStep( ProcessInfo& CurrentProcessInfo )
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        {
            mConstitutiveLawVector[i]->InitializeSolutionStep( GetProperties(),
                    GetGeometry(), row( GetGeometry().ShapeFunctionsValues(  ), i ),
                    CurrentProcessInfo );
        }
    }
    
    //************************************************************************************
    //************************************************************************************

    void KinematicLinear::InitializeNonLinearIteration( ProcessInfo& CurrentProcessInfo )
    {
        // TODO: Add somethig if necessary
    }
    
    //************************************************************************************
    //************************************************************************************

    void KinematicLinear::FinalizeNonLinearIteration( ProcessInfo& CurrentProcessInfo )
    {
        // TODO: Add somethig if necessary
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

        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        const unsigned int StrainSize = dim == 2 ? 3 : 6;

        Matrix F( dim, dim );
        Matrix D( StrainSize, StrainSize );
        Vector StrainVector( StrainSize );
        Vector StressVector( StrainSize );
        Matrix DN_DX( NumberOfNodes, dim );
        Matrix B( StrainSize, NumberOfNodes * dim );
        Matrix J0(dim,dim), InvJ0(dim,dim);

        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRAIN, false);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
        Values.SetStrainVector(StrainVector); //this is the input  parameter
        Values.SetStressVector(StressVector); //this is the output parameter
        
        //reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = GetGeometry().IntegrationPoints(  );

        if ( rOutput.size() != IntegrationPoints.size() )
        {
            rOutput.resize( IntegrationPoints.size() );
        }
        
        Vector displacements;
        GetValuesVector(displacements);

        for ( unsigned int PointNumber = 0; PointNumber < IntegrationPoints.size(); PointNumber++ )
        {
            const Matrix& DN_De = GetGeometry().ShapeFunctionsLocalGradients()[PointNumber];
            double detJ0;
            CalculateDerivativesOnReference(J0, InvJ0, DN_DX, detJ0, DN_De);
            
            //Compute B and strain
            CalculateB( B, DN_DX );
            noalias(StrainVector) = prod(B,displacements);
            F = ComputeEquivalentF(StrainVector);
            
            if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
            {
                if ( rOutput[PointNumber].size2() != StrainVector.size() )
                {
                    rOutput[PointNumber].resize( 1, StrainVector.size(), false );
                }

                for ( unsigned int ii = 0; ii < StrainVector.size(); ii++ )
                {
                    rOutput[PointNumber]( 0, ii ) = StrainVector[ii];
                }
            }
            else if ( rVariable == PK2_STRESS_TENSOR )
            {
                if ( rOutput[PointNumber].size2() != StressVector.size() )
                {
                    rOutput[PointNumber].resize( 1, StressVector.size(), false );
                }
                                        
                // Here we essentially set the input parameters
                const double detF = MathUtils<double>::Det(F);
                Values.SetDeterminantF(detF); //assuming the determinant is computed somewhere else
                Values.SetDeformationGradientF(F); //F computed somewhere else
                
                //actually do the computations in the ConstitutiveLaw    
                mConstitutiveLawVector[PointNumber]->CalculateMaterialResponsePK2(Values); //here the calculations are actually done 

                for ( unsigned int ii = 0; ii < StressVector.size(); ii++ )
                {
                    rOutput[PointNumber]( 0, ii ) = StressVector[ii];
                }
            }
//             else if ( rVariable == GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR )
//             {
//                 double size = StrainVector.size();
//                 Matrix PlasticStrainVector( 1, size );
// 
//                 if ( rOutput[PointNumber].size2() != StrainVector.size() )
//                     rOutput[PointNumber].resize( 1, size, false );
// 
//                 mConstitutiveLawVector[PointNumber]->GetValue( GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR, PlasticStrainVector );
// 
//                 rOutput[PointNumber] = PlasticStrainVector;
//             }
        }

        KRATOS_CATCH( "" )
    }
    
    //************************************************************************************
    //************************************************************************************

    void KinematicLinear::FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo )
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        {
            mConstitutiveLawVector[i]->FinalizeSolutionStep( GetProperties(),
                    GetGeometry(),
                    row( GetGeometry().ShapeFunctionsValues(  ), i ),
                    CurrentProcessInfo );
        }
    }

    //************************************************************************************
    //************************************************************************************

    void KinematicLinear::InitializeMaterial()
    {
        KRATOS_TRY;

        if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
        {
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            {
                mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
                mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), GetGeometry(),
                        row( GetGeometry().ShapeFunctionsValues(  ), i ) );
            }
        }
        else
        {
            KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;
        }
        
        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************
    
    void KinematicLinear::ResetConstitutiveLaw()
    {
        KRATOS_TRY

        if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
        {
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
                mConstitutiveLawVector[i]->ResetMaterial( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues(  ), i ) );
        }

        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************

    inline void KinematicLinear::CalculateAndAdd_ExtForceContribution(
        const Vector& N,
        const ProcessInfo& CurrentProcessInfo,
        Vector& BodyForce,
        VectorType& rRightHandSideVector,
        const double weight
        )
    {
        KRATOS_TRY;
        
        const unsigned int NumberOfNodes = GetGeometry().PointsNumber();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        
        for ( unsigned int i = 0; i < NumberOfNodes; i++ )
        {
            const unsigned int index = dimension * i;

            for ( unsigned int j = 0; j < dimension; j++ ) 
            {
                rRightHandSideVector[index + j] += weight * N[i] * BodyForce[j];
            }
        }

        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************

    void KinematicLinear::CalculateB(
        Matrix& B,
        const Matrix& DN_DX 
        )
    {
        KRATOS_TRY;
        
        const unsigned int NumberOfNodes = GetGeometry().PointsNumber();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();

        B.clear();
        
        if(dim == 2)
        {
            for ( unsigned int i = 0; i < NumberOfNodes; ++i )
            {
                B( 0, i*2 ) = DN_DX( i, 0 );
                B( 1, i*2 + 1 ) = DN_DX( i, 1 );
                B( 2, i*2 ) = DN_DX( i, 1 );
                B( 2, i*2 + 1 ) = DN_DX( i, 0 );
            }
        }
        else if(dim == 3)
        {
            for ( unsigned int i = 0; i < NumberOfNodes; ++i )
            {
                B( 0, i*3 ) = DN_DX( i, 0 );
                B( 1, i*3 + 1 ) = DN_DX( i, 1 );
                B( 2, i*3 + 2 ) = DN_DX( i, 2 );
                B( 3, i*3 ) = DN_DX( i, 1 );
                B( 3, i*3 + 1 ) = DN_DX( i, 0 );
                B( 4, i*3 + 1 ) = DN_DX( i, 2 );
                B( 4, i*3 + 2 ) = DN_DX( i, 1 );
                B( 5, i*3 ) = DN_DX( i, 2 );
                B( 5, i*3 + 2 ) = DN_DX( i, 0 );
            }
        }

        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************

    Matrix KinematicLinear::ComputeEquivalentF(const Vector& StrainVector)
    {
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        Matrix F(dim,dim);
        
        if(dim == 2)
        {
            F(0,0) = 1+StrainVector(0);  
            F(0,1) = 0.5*StrainVector(2); 
            F(1,0) = 0.5*StrainVector(2);   
            F(1,1) = 1+StrainVector(1);
        }
        else
        {
            F(0,0) = 1+StrainVector(0);     
            F(0,1) = 0.5*StrainVector(3); 
            F(0,2) = 0.5*StrainVector(5);
            F(1,0) = 0.5*StrainVector(3);   
            F(1,1) = 1+StrainVector(1);   
            F(1,2) = 0.5*StrainVector(4);
            F(2,0) = 0.5*StrainVector(5);   
            F(2,1) = 0.5*StrainVector(4); 
            F(2,2) = 1+StrainVector(2);
        }
        
        return F;
    }
    
    //***********************************************************************
    //***********************************************************************

    double KinematicLinear::GetIntegrationWeight(
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
        const unsigned int PointNumber,
        const double detJ
        )
    {
        return IntegrationPoints[PointNumber].Weight() * detJ;
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

            if ( this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize() != 3 )
            {
                KRATOS_ERROR << "Wrong constitutive law used. This is a 2D element! expected strain size is 3 (el id = ) " << this->Id() << std::endl;
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


