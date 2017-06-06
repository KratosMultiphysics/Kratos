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
            : Element( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }

//************************************************************************************
//************************************************************************************

    KinematicLinear::KinematicLinear( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
            : Element( NewId, pGeometry, pProperties )
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

        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(  );

        //Constitutive Law initialisation
        if ( mConstitutiveLawVector.size() != integration_points.size() )
        {
            mConstitutiveLawVector.resize( integration_points.size() );
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
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag 
        )
    {
        KRATOS_TRY
        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        const unsigned int StrainSize = ( dim == 2 ) ? 3 : 6;

        Matrix B( StrainSize, NumberOfNodes * dim );
        Matrix D( StrainSize, StrainSize );
        Vector StrainVector( StrainSize );
        Vector StressVector( StrainSize );
        Matrix DN_DX( NumberOfNodes, dim );
        Matrix J0(dim,dim), InvJ0(dim,dim);

        //resizing as needed the LHS
        const unsigned int MatSize = NumberOfNodes * dim;

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
            if ( rRightHandSideVector.size() != MatSize )
            {
                rRightHandSideVector.resize( MatSize, false );
            }

            rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS
        }

        //reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(  );

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

        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            double detJ0;
            ComputeDerivatives(J0, InvJ0, DN_DX, detJ0, PointNumber);
            
            //Compute B and strain
            CalculateB( B, DN_DX );
            noalias(StrainVector) = prod(B,displacements);

            // Compute equivalent F
            Matrix F = ComputeEquivalentF(StrainVector);

            //here we essentially set the input parameters
            double detF = MathUtils<double>::Det(F);
            Values.SetDeterminantF(detF); //assuming the determinant is computed somewhere else
            Values.SetDeformationGradientF(F); //F computed somewhere else
            
            // Here we set the space on which the results shall be written
            Values.SetConstitutiveMatrix(D); //assuming the determinant is computed somewhere else
            Values.SetStressVector(StressVector); //F computed somewhere else
            
            // Actually do the computations in the ConstitutiveLaw    
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponsePK2(Values); //here the calculations are actually done 

            //calculating weights for integration on the reference configuration
            double IntToReferenceWeight = integration_points[PointNumber].Weight() * detJ0;

            if ( dim == 2 && GetProperties().Has( THICKNESS )) 
            {
                IntToReferenceWeight *= GetProperties()[THICKNESS];
            }

            if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
            {
                // Contributions to stiffness matrix calculated on the reference config
                Matrix tmp =  IntToReferenceWeight *  prod( D, B );
                noalias( rLeftHandSideMatrix ) += prod( trans( B ), tmp ); //to be optimized to remove the temporary
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

    void KinematicLinear::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = false;
        bool CalculateResidualVectorFlag = true;
        MatrixType temp = Matrix();

        CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
    }

//************************************************************************************
//************************************************************************************

    void KinematicLinear::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = true;
        bool CalculateResidualVectorFlag = true;

        CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
        


    }

//************************************************************************************
//************************************************************************************

    double KinematicLinear::CalculateIntegrationWeight( double GaussPointWeight, double DetJ0 )
    {
        //to permorm the integration over the reference domain we need to include
        // the thickness in 2D
        unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        double weight = GaussPointWeight;

        weight *= DetJ0;

        if ( dimension == 2 && GetProperties().Has( THICKNESS )) 
        {
            weight *= GetProperties()[THICKNESS];
        }

        return weight;
    }

////************************************************************************************
////************************************************************************************

    void KinematicLinear::InitializeSolutionStep( ProcessInfo& CurrentProcessInfo )
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            mConstitutiveLawVector[i]->InitializeSolutionStep( GetProperties(),
                    GetGeometry(), row( GetGeometry().ShapeFunctionsValues(  ), i ),
                    CurrentProcessInfo );
    }

////************************************************************************************
////************************************************************************************

    void KinematicLinear::FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo )
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            mConstitutiveLawVector[i]->FinalizeSolutionStep( GetProperties(),
                    GetGeometry(),
                    row( GetGeometry().ShapeFunctionsValues(  ), i ),
                    CurrentProcessInfo );
    }

//************************************************************************************
//************************************************************************************

    void KinematicLinear::InitializeMaterial()
    {
        KRATOS_TRY

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
            KRATOS_ERROR << "a constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;
            KRATOS_CATCH( "" )
        }

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
        double weight
    )
    {
        KRATOS_TRY
        
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
        const Matrix& DN_DX )
    {
        KRATOS_TRY
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

    void KinematicLinear::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo )
    {
        int NumberOfNodes = GetGeometry().size();
        int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int dim2 = NumberOfNodes * dim;

        if ( rResult.size() != dim2 )
            rResult.resize( dim2, false );

        for ( int i = 0; i < NumberOfNodes; i++ )
        {
            int index = i * dim;
            rResult[index] = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();

            if ( dim == 3 )
                rResult[index + 2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
        }

    }

//************************************************************************************
//************************************************************************************

    void KinematicLinear::GetDofList( DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo )
    {
        ElementalDofList.resize( 0 );

        for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
        {
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );

            if ( GetGeometry().WorkingSpaceDimension() == 3 )
            {
                ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
            }
        }
    }

//************************************************************************************
//************************************************************************************

    void KinematicLinear::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        //lumped
        unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        unsigned int NumberOfNodes = GetGeometry().size();
        unsigned int MatSize = dimension * NumberOfNodes;

        if ( rMassMatrix.size1() != MatSize )
            rMassMatrix.resize( MatSize, MatSize, false );

        rMassMatrix = ZeroMatrix( MatSize, MatSize );
        
        
        
        Matrix DN_DX( NumberOfNodes, dimension );
        Matrix J0(dimension,dimension), InvJ0(dimension,dimension);
        
        //reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(  );

        double initial_area = 0.0;
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            double detJ0;
            ComputeDerivatives(J0, InvJ0, DN_DX, detJ0, PointNumber);
            initial_area += detJ0;
        }

        double TotalMass = initial_area * GetProperties()[DENSITY];

        if ( dimension == 2 && GetProperties().Has( THICKNESS )) 
        {
            TotalMass *= GetProperties()[THICKNESS];
        }

        Vector LumpFact;

        LumpFact = GetGeometry().LumpingFactors( LumpFact );

        for ( unsigned int i = 0; i < NumberOfNodes; i++ )
        {
            double temp = LumpFact[i] * TotalMass;

            for ( unsigned int j = 0; j < dimension; j++ )
            {
                unsigned int index = i * dimension + j;
                rMassMatrix( index, index ) = temp;
            }
        }

        KRATOS_CATCH( "" )
    }

//************************************************************************************
//************************************************************************************

    void KinematicLinear::CalculateDampingMatrix( MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY
        unsigned int NumberOfNodes = GetGeometry().size();
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        //resizing as needed the LHS
        unsigned int MatSize = NumberOfNodes * dim;

        if ( rDampingMatrix.size1() != MatSize )
            rDampingMatrix.resize( MatSize, MatSize, false );

        noalias( rDampingMatrix ) = ZeroMatrix( MatSize, MatSize );

        KRATOS_CATCH( "" )
    }

//************************************************************************************
//************************************************************************************

    void KinematicLinear::CalculateOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& Output, const ProcessInfo& rCurrentProcessInfo )
    {
        if ( Output.size() != GetGeometry().IntegrationPoints(  ).size() )
            Output.resize( GetGeometry().IntegrationPoints(  ).size() );

        for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
            Output[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, Output[ii] );
    }

//************************************************************************************
//************************************************************************************

    void KinematicLinear::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& Output, const ProcessInfo& rCurrentProcessInfo )
    {
        unsigned int StrainSize;

        if ( GetGeometry().WorkingSpaceDimension() == 2 ) StrainSize = 3;
        else StrainSize = 6;

        Vector StrainVector( StrainSize );

        if ( rVariable == INSITU_STRESS )
        {
            for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
            {
                if ( Output[ii].size() != StrainVector.size() )
                    Output[ii].resize( StrainVector.size(), false );

                Output[ii] = mConstitutiveLawVector[ii]->GetValue( INSITU_STRESS, Output[ii] );
            }
        }
        else
        {
            if ( Output.size() != GetGeometry().IntegrationPoints(  ).size() )
                Output.resize( GetGeometry().IntegrationPoints(  ).size() );

            for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
                Output[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, Output[ii] );
        }

    }

//************************************************************************************
//************************************************************************************

    void KinematicLinear::CalculateOnIntegrationPoints( const Variable<Matrix >& rVariable, std::vector< Matrix >& Output, const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int StrainSize;

        if ( dim == 2 )
            StrainSize = 3;
        else
            StrainSize = 6;

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
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(  );

        if ( Output.size() != integration_points.size() )
            Output.resize( integration_points.size() );
        
        Vector displacements;
        GetValuesVector(displacements);

        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            double detJ0;
            ComputeDerivatives(J0, InvJ0, DN_DX, detJ0, PointNumber);
            
            //Compute B and strain
            
            CalculateB( B, DN_DX );
            noalias(StrainVector) = prod(B,displacements);
            F = ComputeEquivalentF(StrainVector);
            
            if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
            {
                if ( Output[PointNumber].size2() != StrainVector.size() )
                    Output[PointNumber].resize( 1, StrainVector.size(), false );

                for ( unsigned int ii = 0; ii < StrainVector.size(); ii++ )
                    Output[PointNumber]( 0, ii ) = StrainVector[ii];
            }
            else if ( rVariable == PK2_STRESS_TENSOR )
            {
                if ( Output[PointNumber].size2() != StressVector.size() )
                    Output[PointNumber].resize( 1, StressVector.size(), false );
                                        
                //here we essentially set the input parameters
                double detF = MathUtils<double>::Det(F);
                Values.SetDeterminantF(detF); //assuming the determinant is computed somewhere else
                Values.SetDeformationGradientF(F); //F computed somewhere else
                
                //actually do the computations in the ConstitutiveLaw    
                mConstitutiveLawVector[PointNumber]->CalculateMaterialResponsePK2(Values); //here the calculations are actually done 

                for ( unsigned int ii = 0; ii < StressVector.size(); ii++ )
                {
                    Output[PointNumber]( 0, ii ) = StressVector[ii];
                }
            }
//             else if ( rVariable == GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR )
//             {
//                 double size = StrainVector.size();
//                 Matrix PlasticStrainVector( 1, size );
// 
//                 if ( Output[PointNumber].size2() != StrainVector.size() )
//                     Output[PointNumber].resize( 1, size, false );
// 
//                 mConstitutiveLawVector[PointNumber]->GetValue( GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR, PlasticStrainVector );
// 
//                 Output[PointNumber] = PlasticStrainVector;
//             }
        }

        KRATOS_CATCH( "" )
    }

//************************************************************************************
//************************************************************************************

    void KinematicLinear::SetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {

        for ( unsigned int PointNumber = 0; PointNumber < GetGeometry().IntegrationPoints(  ).size(); PointNumber++ )
        {
            mConstitutiveLawVector[PointNumber]->SetValue( rVariable,
                    rValues[PointNumber], rCurrentProcessInfo );
        }

    }


//************************************************************************************
//************************************************************************************

    void KinematicLinear::SetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        for ( unsigned int PointNumber = 0; PointNumber < GetGeometry().IntegrationPoints(  ).size(); PointNumber++ )
        {
            mConstitutiveLawVector[PointNumber]->SetValue( rVariable,
                    rValues[PointNumber], rCurrentProcessInfo );
        }

    }

//************************************************************************************
//************************************************************************************

    void KinematicLinear::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
            std::vector<double>& rValues,
            const ProcessInfo& rCurrentProcessInfo )
    {
        if ( rValues.size() != GetGeometry().IntegrationPoints(  ).size() )
            rValues.resize( GetGeometry().IntegrationPoints(  ).size(), false );

        for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
            rValues[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rValues[ii] );
    }


//************************************************************************************
//************************************************************************************

    void KinematicLinear::GetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        const unsigned int& size = GetGeometry().IntegrationPoints(  ).size();

        if ( rValues.size() != size )
            rValues.resize( size );

        //TODO: decide which is the correct one
  /*      if ( rVariable == INSITU_STRESS || rVariable == PRESTRESS )
        {
            for ( unsigned int PointNumber = 0;
                    PointNumber < GetGeometry().IntegrationPoints(  ).size();
                    PointNumber++ )
            {
                rValues[PointNumber] =
                    mConstitutiveLawVector[PointNumber]->GetValue( PRESTRESS, rValues[PointNumber] );
            }
        }
        if ( rVariable == PLASTIC_STRAIN_VECTOR )
        {
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            {
                if ( rValues[i].size() != 6 )
                    rValues[i].resize( 6 );
                noalias( rValues[i] ) = mConstitutiveLawVector[i]->GetValue( PLASTIC_STRAIN_VECTOR, rValues[i] );
            }
        }
  */      if ( rVariable == STRESSES )
        {
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            {
                if ( rValues[i].size() != 6 )
                    rValues[i].resize( 6 );
                noalias( rValues[i] ) = mConstitutiveLawVector[i]->GetValue( STRESSES, rValues[i] );
            }
        }
        if ( rVariable == MATERIAL_PARAMETERS )
        {
            for ( unsigned int PointNumber = 0;
                    PointNumber < GetGeometry().IntegrationPoints(  ).size(); PointNumber++ )
            {
                rValues[PointNumber] =
                    mConstitutiveLawVector[PointNumber]->GetValue( MATERIAL_PARAMETERS, rValues[PointNumber] );
            }
        }

        if ( rVariable == INTERNAL_VARIABLES )
        {
            for ( unsigned int PointNumber = 0;
                    PointNumber < GetGeometry().IntegrationPoints(  ).size();
                    PointNumber++ )
            {
                rValues[PointNumber] =
                    mConstitutiveLawVector[PointNumber]->GetValue( INTERNAL_VARIABLES, rValues[PointNumber] );

            }
        }


    }

//************************************************************************************
//************************************************************************************

    void KinematicLinear::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable,
            std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
        {
            CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
        }

        if ( rVariable == PK2_STRESS_TENSOR )
        {
            CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
        }

//         if ( rVariable == GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR )
//         {
//             CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
//         }

    }

//************************************************************************************
//************************************************************************************

    void KinematicLinear::GetValuesVector( Vector& values, int Step )
    {
        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        const unsigned int MatSize = NumberOfNodes * dim;

        if ( values.size() != MatSize ) values.resize( MatSize, false );

        for ( unsigned int i = 0; i < NumberOfNodes; i++ )
        {
            const unsigned int index = i * dim;
            
            const auto& d  = GetGeometry()[i].FastGetSolutionStepValue( DISPLACEMENT, Step );
            
            for(unsigned int k=0; k<dim;k++)
                values[index+k] = d[k];
        }
    }


//************************************************************************************
//************************************************************************************

    void KinematicLinear::GetFirstDerivativesVector( Vector& values, int Step )
    {
        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int MatSize = NumberOfNodes * dim;

        if ( values.size() != MatSize ) values.resize( MatSize, false );

        for ( unsigned int i = 0; i < NumberOfNodes; i++ )
        {
            const unsigned int index = i * dim;
            
            const auto& d  = GetGeometry()[i].FastGetSolutionStepValue( VELOCITY, Step );
            
            for(unsigned int k=0; k<dim;k++)
                values[index+k] = d[k];
        }
    }

//************************************************************************************
//************************************************************************************

    void KinematicLinear::GetSecondDerivativesVector( Vector& values, int Step )
    {
        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int MatSize = NumberOfNodes * dim;

        if ( values.size() != MatSize ) values.resize( MatSize, false );
        for ( unsigned int i = 0; i < NumberOfNodes; i++ )
        {
            const unsigned int index = i * dim;
            
            const auto& d  = GetGeometry()[i].FastGetSolutionStepValue( ACCELERATION, Step );
            
            for(unsigned int k=0; k<dim;k++)
                values[index+k] = d[k];
        }    }

//************************************************************************************
//************************************************************************************

    void KinematicLinear::Calculate( const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo )
    {

        //HERE IT IS ESTIMATING DELTA_TIME for explicit elements
//         double lamda = 1.00; // parametro que depende del tipo de problema y del elemento pag 308 libro dinamica de Barbat
//         double c1 = 0.00; //sqrt(GetProperties()[YOUNG_MODULUS]/GetProperties()[DENSITY]); velocidad del sonido en el medio
//         double c2 = 0.00; // norma de la velocidad actual dentro del elemento
//         double c = 0.00;
//         double wmax = 0.00;
//         Vector Values( GetGeometry().IntegrationPoints(  ).size() );
//         Vector Velocities;
// 
//         GetFirstDerivativesVector( Velocities, 0 );
// 
//         if ( rVariable == DELTA_TIME )
//         {
//             for ( unsigned int PointNumber = 0;
//                     PointNumber < GetGeometry().IntegrationPoints(  ).size();
//                     PointNumber++ )
//             {
//                 mConstitutiveLawVector[PointNumber]-> GetValue( DELTA_TIME, c1 );
//                 Values[PointNumber] = c1;
//             }
//         }
// 
//         c1 = ( *std::max_element( Values.begin(), Values.end() ) );
// 
//         c2 = norm_2( Velocities );
// 
//         c = ( c1 > c2 ) ? c1 : c2;
// 
// 
//         double le = GetGeometry().Length();
//         //KRATOS_WATCH(le)
// 
//         /// maxima frecuencia de un elemento
//         wmax = ( lamda * c ) / le;
//         Output = 2.0 / wmax;
//         //KRATOS_WATCH(Output)

    }

//************************************************************************************
//************************************************************************************
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int  KinematicLinear::Check( const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();



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

    void  KinematicLinear::ComputeDerivatives(Matrix& J0, 
                             Matrix& InvJ0, 
                             Matrix& DN_DX, 
                             double& detJ0, 
                             const unsigned int PointNumber)
    {
        auto& DN_De = GetGeometry().ShapeFunctionsLocalGradients()[PointNumber];
        
        J0.clear();
        for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
        {
            const auto& coords = GetGeometry()[i].GetInitialPosition(); //NOTE: here we refer to the original, undeformed position!!
            for(unsigned int k=0; k<GetGeometry().WorkingSpaceDimension(); k++)
            {
                for(unsigned int m=0; m<GetGeometry().LocalSpaceDimension(); m++)
                {
                    J0(k,m) += coords[k]*DN_De(i,m);
                }
            }
        }
        
        MathUtils<double>::InvertMatrix( J0, InvJ0, detJ0 );
        
        noalias( DN_DX ) = prod( DN_De, InvJ0);
    }

    void KinematicLinear::save( Serializer& rSerializer ) const
    {
        rSerializer.save( "Name", "KinematicLinear" );
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element );
    }

    void KinematicLinear::load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element );
    }





} // Namespace Kratos


