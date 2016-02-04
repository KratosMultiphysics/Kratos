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
#include "custom_elements/membrane_element.hpp"
#include "structural_mechanics_application_variables.h"


namespace Kratos
{


// Constructor
MembraneElement::MembraneElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : Element( NewId, pGeometry )
{
}

// Constructor
MembraneElement::MembraneElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : Element( NewId, pGeometry, pProperties )

{
}

//***********************************************************************************
//***********************************************************************************

Element::Pointer MembraneElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties ) const

{
    return Element::Pointer( new MembraneElement( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

//***********************************************************************************
//***********************************************************************************
// Destructor
MembraneElement::~MembraneElement()
{
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo )

{
    KRATOS_TRY
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int dim = number_of_nodes * 3;

    if ( rResult.size() != dim )
        rResult.resize( dim );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int index = i * 3;
        rResult[index]   = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[index+1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();
        rResult[index+2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
    }

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::GetDofList(
    DofsVectorType& ElementalDofList,
    ProcessInfo& rCurrentProcessInfo )

{
    ElementalDofList.resize( 0 );

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
    {
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
    }
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::Initialize()

{
    KRATOS_TRY
    //KRATOS_WATCH( "INITIALIZE ELEMENT" )
    //getting all "Actual" info from the geometry
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();

    //resizing jacobian inverses containers
    mDetJ0.resize( integration_points.size() );

    GeometryType::JacobiansType J0;
    J0 = GetGeometry().Jacobian( J0 );
    mTotalDomainInitialSize = 0.00;

    mStrainsVector.resize( integration_points.size() );
    mStressesVector.resize( integration_points.size() );
    mCauchyStressesVector.resize( integration_points.size() ); //VM

    mV1.resize( integration_points.size() );
    mV2.resize( integration_points.size() );
    mG_Vector.resize( integration_points.size(), ZeroMatrix( 2, 2 ) );

    array_1d<double, 3> g0e;
    array_1d<double, 3> g0n;
    array_1d<double, 3> V1;
    array_1d<double, 3> V2;
    array_1d<double, 3> V3;
    array_1d<double, 3> N;

    // Initialize Variables
    mdensity = GetProperties()[DENSITY];
    mThickness0 = GetProperties()[THICKNESS];
    //KRATOS_WATCH( mdensity )
    //KRATOS_WATCH( mThickness0 )
    mThickness.resize( integration_points.size(), 0.00 );

    // Initialize Material

    if ( mConstitutiveLawVector.size() == 0 )
    {
        mConstitutiveLawVector.resize( integration_points.size() );

        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        {
            /*    ConstitutiveLaw<Node<3> >::Pointer material = ConstitutiveLaw<Node<3> >::Pointer( new Isotropic2D() );
                mConstitutiveLawVector[i] = material;*/


            mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();

            mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues(), i ) );
        }
    }

    //calculating the inverse J0
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        //getting informations for integration
        double IntegrationWeight = integration_points[PointNumber].Weight();

        //calculating and storing inverse of the jacobian and the parameters needed
        g0e[0] = J0[PointNumber]( 0, 0 );
        g0n[0] = J0[PointNumber]( 0, 1 );
        g0e[1] = J0[PointNumber]( 1, 0 );
        g0n[1] = J0[PointNumber]( 1, 1 );
        g0e[2] = J0[PointNumber]( 2, 0 );
        g0n[2] = J0[PointNumber]( 2, 1 );

        //calculate base vectors
        CrossProduct( V3, g0e, g0n );
        N = V3;
        N /= norm_2( V3 );
        V1 = g0e;
        V1 /= norm_2( V1 );
        CrossProduct( V2, N, V1 );

        //saving the initial local base vectors
        mV1[PointNumber] = V1;
        mV2[PointNumber] = V2;

        // Calculation of Matrix G (a sort of inverse jacobian)
        double J11 = norm_2( g0e );
        double J12 = inner_prod( g0e, g0n ) / norm_2( g0e );
        double J22 = norm_2( V3 ) / norm_2( g0e );

        Matrix G( 2, 2 );
        G( 0, 0 ) = 1 / J11;
        G( 0, 1 ) = -J12 / ( J11 * J22 );
        G( 1, 0 ) = 0.00;
        G( 1, 1 ) = 1 / J22;

        //saving the G matrix for the point number
        noalias( mG_Vector[PointNumber] ) = G;

        //Calculate the reduced mass matrix
        mDetJ0[PointNumber] = norm_2( V3 );

        //calculating the total area
        mTotalDomainInitialSize += mDetJ0[PointNumber] * IntegrationWeight;
    }

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo )

{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo )

{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;

    CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::CalculateOnIntegrationPoints(
    const Variable<Matrix>& rVariable,
    std::vector<Matrix>& Output,
    const ProcessInfo& rCurrentProcessInfo )

{

    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
    //const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues();

    //calculating actual jacobian
    GeometryType::JacobiansType J;
    J = GetGeometry().Jacobian( J );

    //auxiliary terms
    boost::numeric::ublas::bounded_matrix<double, 2, 2> j;
    boost::numeric::ublas::bounded_matrix<double, 2, 2> g;
    boost::numeric::ublas::bounded_matrix<double, 2, 2> C;
    array_1d<double, 3> ge;
    array_1d<double, 3> gn;
    array_1d<double, 3> v3;

    Vector StrainVector( 3 );
    Vector StressVector( 3 );


    if ( Output.size() != integration_points.size() )
        Output.resize( integration_points.size() );

    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        //double IntegrationWeight = GetGeometry().IntegrationPoints()[PointNumber].Weight();

        ge[0] = J[PointNumber]( 0, 0 );
        gn[0] = J[PointNumber]( 0, 1 );
        ge[1] = J[PointNumber]( 1, 0 );
        gn[1] = J[PointNumber]( 1, 1 );
        ge[2] = J[PointNumber]( 2, 0 );
        gn[2] = J[PointNumber]( 2, 1 );

        CrossProduct( v3, ge, gn );
        CalculateJ( j, ge, gn, v3 );

        // Calculation of matrix g = jtrans*j;
        noalias( g ) = prod( trans( j ), j );

        // calculation of the Right Cauchy-Green Tensor C = Gtrans*g*G
        boost::numeric::ublas::bounded_matrix<double, 2, 2> tmp;
        tmp = prod( g, mG_Vector[PointNumber] );
        noalias( C ) = prod( trans( mG_Vector[PointNumber] ), tmp );

        // Calculation of the StrainVector
        CalculateStrain( StrainVector, C );

        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRAIN, false);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        // if strain has to be computed inside of the constitutive law with PK2
        //rValues.SetDeformationGradientF(rVariables.F); //in this case F is the whole deformation gradient

        Values.SetStrainVector(StrainVector); //this is the input  parameter
        Values.SetStressVector(StressVector); //this is the output parameter


        mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(Values,ConstitutiveLaw::StressMeasure_PK2 );


        noalias( mStressesVector[PointNumber] ) = ZeroVector( 6 );
        Calculate_GlobalStressVector( mStressesVector[PointNumber], StressVector, mV1[PointNumber], mV2[PointNumber] ); //saving the stress vector


        if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
        {
//             if ( Output[PointNumber].size2() != StrainVector.size() )
//                 Output[PointNumber].resize( 1, StrainVector.size() );
//
//             for ( unsigned int ii = 0; ii < StrainVector.size(); ii++ )
//                 Output[PointNumber]( 0, ii ) = StrainVector[ii];

            //BIG CHAPUZA! need to have this to look like a stress
            Vector strain_as_stress = StrainVector;
            strain_as_stress[2] *= 0.5;
            array_1d<double,6> global_strain = ZeroVector(6);
            Calculate_GlobalStressVector( global_strain, strain_as_stress, mV1[PointNumber], mV2[PointNumber] );

            Matrix StrainMatrix = MathUtils<double>::StrainVectorToTensor(global_strain);
            Output[PointNumber] = StrainMatrix;
        }
        else if ( rVariable == PK2_STRESS_TENSOR )
        {
//             if ( Output[PointNumber].size2() != 6 )
//                 Output[PointNumber].resize( 1, 6 );
//
//             for ( unsigned int ii = 0; ii < 6; ii++ )
//                 Output[PointNumber]( 0, ii ) = mStressesVector[PointNumber][ii];

            Matrix StressMatrix = MathUtils<double>::StressVectorToTensor(mStressesVector[PointNumber]);
            Output[PointNumber] = StressMatrix;
        }
        else if(rVariable==CAUCHY_STRESS_TENSOR)  // to compute Cauchy_Stress
        {
//             if(Output[PointNumber].size2() != 6)
//                 Output[PointNumber].resize(1,6);

            boost::numeric::ublas::bounded_matrix<double, 2, 2> F;
            noalias(F)=tmp; //VM
            Vector CauchyStressVector = StressVector;
            double detF = MathUtils<double>::Det(F);

            //calculate base vectors in the current configuration
            array_1d<double, 3> v1,v2,n;
            CrossProduct( v3, ge, gn );
            n = v3;
            n /= norm_2( v3 );
            v1 = ge;
            v1 /= norm_2( v1 );
            CrossProduct( v2, n, v1 );


            mConstitutiveLawVector[PointNumber]->TransformPK2Stresses(CauchyStressVector,F,detF,ConstitutiveLaw::StressMeasure_Cauchy);

            noalias(mCauchyStressesVector[PointNumber])= ZeroVector(6);
            Calculate_GlobalStressVector(mCauchyStressesVector[PointNumber], CauchyStressVector, v1, v2);   //saving the stress vector
            // Calculate_GlobalStressVector(mCauchyStressesVector[PointNumber], CauchyStressVector, mV1[PointNumber], mV2[PointNumber]);	//saving the stress vector

            // for(unsigned int ii = 0; ii<6; ii++)
            //   Output[PointNumber](0,ii) = mCauchyStressesVector[PointNumber][ii];

            Matrix StressMatrix = MathUtils<double>::StressVectorToTensor(mCauchyStressesVector[PointNumber]);
            Output[PointNumber] = StressMatrix;

        }
        else
        {
            Output[PointNumber] = mConstitutiveLawVector[PointNumber]->GetValue( rVariable, Output[PointNumber] );
        }


    }

}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    ProcessInfo& rCurrentProcessInfo )

{
    KRATOS_TRY

    //rMassMatrix.resize(0,0);
    // LUMPED MASS MATRIX
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes * 3;

    if ( rMassMatrix.size1() != MatSize )
        rMassMatrix.resize( MatSize, MatSize );

    rMassMatrix = ZeroMatrix( MatSize, MatSize );

    double TotalMass = mTotalDomainInitialSize * GetProperties()[THICKNESS] * GetProperties()[DENSITY];

    Vector LumpFact;

    LumpFact = GetGeometry().LumpingFactors( LumpFact );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        double temp = LumpFact[i] * TotalMass;

        for ( unsigned int j = 0; j < 3; j++ )
        {
            unsigned int index = i * 3 + j;
            rMassMatrix( index, index ) = temp;
        }
    }

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    ProcessInfo& rCurrentProcessInfo )

{
    KRATOS_TRY

    // LUMPED DAMPING MATRIX
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes * 3;

    if ( rDampingMatrix.size1() != MatSize )
        rDampingMatrix.resize( MatSize, MatSize );

    rDampingMatrix = ZeroMatrix( MatSize, MatSize );

    double TotalMass = mTotalDomainInitialSize * GetProperties()[THICKNESS] * GetProperties()[DENSITY];

    Vector LumpFact;

    LumpFact = GetGeometry().LumpingFactors( LumpFact );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        double temp = LumpFact[i] * TotalMass;

        for ( unsigned int j = 0; j < 3; j++ )
        {
            unsigned int index = i * 3 + j;
            rDampingMatrix( index, index ) = temp;
        }
    }

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::FinalizeSolutionStep(
    ProcessInfo& rCurrentProcessInfo )
{
    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
    {
//
//            ConstitutiveLaw::Parameters Values (GetGeometry(),GetProperties(),rCurrentProcessInfo);
//            Values.GetOptions().Set (ConstitutiveLaw::COMPUTE_STRAIN, false);
//            Values.GetOptions().Set (ConstitutiveLaw::COMPUTE_STRESS);
//            Values.GetOptions().Set (ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
//            Matrix dummy = ZeroMatrix ( 0, 0 );
//            Vector StrainVector = mStrainsVector[i];
//            Values.SetStrainVector (StrainVector); //this has to be the input parameter
//            Values.SetStressVector (StressVector);
//            Values.SetConstitutiveMatrix (dummy);
//            Values.SetShapeFunctionsValues ( row ( GetGeometry().ShapeFunctionsValues(), PointNumber ) );
//
//            mConstitutiveLawVector[PointNumber]->FinalizeMaterialResponse (Values,ConstitutiveLaw::StressMeasure_PK2 );

        mConstitutiveLawVector[i]->FinalizeSolutionStep( GetProperties(),
                GetGeometry(),
                row( GetGeometry().ShapeFunctionsValues(), i ),
                rCurrentProcessInfo );
    }
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::GetValuesVector(
    Vector& values,
    int Step )

{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int MatSize = number_of_nodes * 3;

    if ( values.size() != MatSize )
        values.resize( MatSize );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        const array_1d<double, 3>& disp = GetGeometry()[i].FastGetSolutionStepValue( DISPLACEMENT, Step );
        unsigned int index = i * 3;
        values[index]   = disp[0];
        values[index+1] = disp[1];
        values[index+2] = disp[2];
    }
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::GetFirstDerivativesVector(
    Vector& values,
    int Step )

{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int MatSize = number_of_nodes * 3;

    if ( values.size() != MatSize )
        values.resize( MatSize );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        const array_1d<double, 3>& vel = GetGeometry()[i].FastGetSolutionStepValue( VELOCITY, Step );
        unsigned int index = i * 3;
        values[index]   = vel[0];
        values[index+1] = vel[1];
        values[index+2] = vel[2];
    }

}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::GetSecondDerivativesVector(
    Vector& values,
    int Step )

{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int MatSize = number_of_nodes * 3;

    if ( values.size() != MatSize )
        values.resize( MatSize );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        const array_1d<double, 3>& acc = GetGeometry()[i].FastGetSolutionStepValue( ACCELERATION, Step );
        unsigned int index = i * 3;
        values[index]   = acc[0];
        values[index+1] = acc[1];
        values[index+2] = acc[2];
    }
}

//***********************************************************************************
//***********************************************************************************
// --------- //
//  PRIVATE  //
// --------- //
//***********************************************************************************
//***********************************************************************************

void MembraneElement::CalculateAndAddKm(
    Matrix& K,
    Matrix& B,
    Matrix& D,
    double weight )

{
    KRATOS_TRY

    unsigned int dim = B.size2();
    Matrix temp( 3, dim );
    noalias( temp ) = prod( D, B );
    temp *= weight;
    Matrix Km( dim, dim );
    noalias( Km ) = prod( trans( B ), temp );
    noalias( K ) += Km;

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::CalculateAndAddKg(
    Matrix& K,
    boost::numeric::ublas::bounded_matrix<double, 3, 3>& Q,
    const Matrix& DN_De,
    Vector& StressVector,
    double weight )

{
    KRATOS_TRY

    unsigned int number_of_nodes = GetGeometry().size();
    Vector s( 3 );
    noalias( s ) = prod( trans( Q ), StressVector );

    double s11 = s[0];
    double s22 = s[1];
    double s12 = s[2];

    Matrix Kloc( number_of_nodes, number_of_nodes );
    Vector a = ZeroVector( number_of_nodes );
    Vector b = ZeroVector( number_of_nodes );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        a[i] = DN_De( i, 0 );
        b[i] = DN_De( i, 1 );
    }

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
            Kloc( i, j ) = a[i] * a[j] * s11 + b[i] * b[j] * s22 + ( b[i] * a[j] + a[i] * b[j] ) * s12;
        }
    }

    Kloc *= weight;

    ExpandReducedMatrix( K, Kloc );

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::CalculateAndSubKp(
    Matrix& K,
    array_1d<double, 3>& ge,
    array_1d<double, 3>& gn,
    const Matrix& DN_De,
    const Vector& N,
    double pressure,
    double weight )

{
    KRATOS_TRY

    boost::numeric::ublas::bounded_matrix<double, 3, 3> Kij;
    boost::numeric::ublas::bounded_matrix<double, 3, 3> Cross_ge;
    boost::numeric::ublas::bounded_matrix<double, 3, 3> Cross_gn;
    double coeff;
    unsigned int number_of_nodes = GetGeometry().size();

    MakeCrossMatrix( Cross_ge, ge );
    MakeCrossMatrix( Cross_gn, gn );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int RowIndex = i * 3;

        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
            int ColIndex = j * 3;

            coeff = pressure * N[i] * DN_De( j, 1 ) * weight;
            noalias( Kij )  = coeff * Cross_ge;

            coeff = pressure * N[i] * DN_De( j, 0 ) * weight;

            noalias( Kij ) -= coeff * Cross_gn;
//Kij *= -1;
            //TAKE CARE: the load correction matrix should be SUBTRACTED not added
            SubtractMatrix( K, Kij, RowIndex, ColIndex );
        }
    }

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************
void MembraneElement::ClearNodalForces()
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        if( GetGeometry()[i].SolutionStepsDataHas(EXTERNAL_FORCE) && GetGeometry()[i].SolutionStepsDataHas(INTERNAL_FORCE) )
        {

            array_1d<double, 3 > & ExternalForce = GetGeometry()[i].FastGetSolutionStepValue(EXTERNAL_FORCE);
            array_1d<double, 3 > & InternalForce = GetGeometry()[i].FastGetSolutionStepValue(INTERNAL_FORCE);

            GetGeometry()[i].SetLock();
            ExternalForce.clear();
            InternalForce.clear();
            GetGeometry()[i].UnSetLock();

        }

    }

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::AddExplicitContribution(const VectorType& rRHSVector,
        const Variable<VectorType>& rRHSVariable,
        Variable<array_1d<double,3> >& rDestinationVariable,
        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    if( rRHSVariable == EXTERNAL_FORCES_VECTOR && rDestinationVariable == EXTERNAL_FORCE )
    {

        for(unsigned int i=0; i< number_of_nodes; i++)
        {
            int index = dimension * i;

            GetGeometry()[i].SetLock();

            array_1d<double, 3 > &ExternalForce = GetGeometry()[i].FastGetSolutionStepValue(EXTERNAL_FORCE);
            for(unsigned int j=0; j<dimension; j++)
            {
                ExternalForce[j] += rRHSVector[index + j];
            }

            GetGeometry()[i].UnSetLock();
        }
    }

    if( rRHSVariable == INTERNAL_FORCES_VECTOR && rDestinationVariable == INTERNAL_FORCE )
    {

        for(unsigned int i=0; i< number_of_nodes; i++)
        {
            int index = dimension * i;

            GetGeometry()[i].SetLock();

            array_1d<double, 3 > &InternalForce = GetGeometry()[i].FastGetSolutionStepValue(INTERNAL_FORCE);
            for(unsigned int j=0; j<dimension; j++)
            {
                InternalForce[j] += rRHSVector[index + j];
            }

            GetGeometry()[i].UnSetLock();
        }
    }


    if( rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL )
    {

        for(unsigned int i=0; i< number_of_nodes; i++)
        {
            int index = dimension * i;

            GetGeometry()[i].SetLock();

            array_1d<double, 3 > &ForceResidual = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);
            for(unsigned int j=0; j<dimension; j++)
            {
                ForceResidual[j] += rRHSVector[index + j];
            }

            GetGeometry()[i].UnSetLock();
        }
    }

    KRATOS_CATCH( "" )

}


//***********************************************************************************
//***********************************************************************************




void MembraneElement::MakeCrossMatrix(
    boost::numeric::ublas::bounded_matrix<double, 3, 3>& M,
    array_1d<double, 3>& U )

{
    M( 0, 0 ) =  0.00;
    M( 0, 1 ) = -U[2];
    M( 0, 2 ) =  U[1];
    M( 1, 0 ) =  U[2];
    M( 1, 1 ) =  0.00;
    M( 1, 2 ) = -U[0];
    M( 2, 0 ) = -U[1];
    M( 2, 1 ) =  U[0];
    M( 2, 2 ) =  0.00;
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::CrossProduct(
    array_1d<double, 3>& cross,
    array_1d<double, 3>& a,
    array_1d<double, 3>& b )

{
    cross[0] = a[1] * b[2] - a[2] * b[1];
    cross[1] = a[2] * b[0] - a[0] * b[2];
    cross[2] = a[0] * b[1] - a[1] * b[0];
}

//***********************************************************************************
//***********************************************************************************

void  MembraneElement::ExpandReducedMatrix(
    Matrix& Destination,
    Matrix& ReducedMatrix )

{
    KRATOS_TRY

    unsigned int size = ReducedMatrix.size2();

    for ( unsigned int i = 0; i < size; i++ )
    {
        int rowindex = i * 3;

        for ( unsigned int j = 0; j < size; j++ )
        {
            unsigned int colindex = j * 3;

            for ( unsigned int ii = 0; ii < 3; ii++ )
                Destination( rowindex + ii, colindex + ii ) += ReducedMatrix( i, j );
        }
    }

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void  MembraneElement::SubtractMatrix(
    MatrixType& Destination,
    boost::numeric::ublas::bounded_matrix<double, 3, 3>& InputMatrix,
    int InitialRow,
    int InitialCol )

{
    KRATOS_TRY

    for ( unsigned int i = 0; i < 3; i++ )
        for ( unsigned int j = 0; j < 3; j++ )
            Destination( InitialRow + i, InitialCol + j ) -= InputMatrix( i, j );

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::CalculateQ(
    boost::numeric::ublas::bounded_matrix<double, 3, 3>& Q,
    Matrix& mG )

{
    KRATOS_TRY

    Q( 0, 0 ) = pow( mG( 0, 0 ), 2 );
    Q( 1, 0 ) = pow( mG( 0, 1 ), 2 );
    Q( 1, 1 ) = pow( mG( 1, 1 ), 2 );
    Q( 1, 2 ) = mG( 0, 1 ) * mG( 1, 1 );
    Q( 2, 0 ) = 2.00 * mG( 0, 0 ) * mG( 0, 1 );
    Q( 2, 2 ) = mG( 0, 0 ) * mG( 1, 1 );

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::CalculateB(
    Matrix& B,
    boost::numeric::ublas::bounded_matrix<double, 3, 3>& Q,
    const Matrix& DN_De,
    array_1d<double, 3>& ge,
    array_1d<double, 3>& gn )

{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    Matrix b( 3, number_of_nodes*3 );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = 3 * i;

        //first line
        b( 0, index )   = DN_De( i, 0 ) * ge[0];
        b( 0, index + 1 ) = DN_De( i, 0 ) * ge[1];
        b( 0, index + 2 ) = DN_De( i, 0 ) * ge[2];

        //second line
        b( 1, index )   = DN_De( i, 1 ) * gn[0];
        b( 1, index + 1 ) = DN_De( i, 1 ) * gn[1];
        b( 1, index + 2 ) = DN_De( i, 1 ) * gn[2];

        //third line
        b( 2, index )   = DN_De( i, 1 ) * ge[0] + DN_De( i, 0 ) * gn[0];
        b( 2, index + 1 ) = DN_De( i, 1 ) * ge[1] + DN_De( i, 0 ) * gn[1];
        b( 2, index + 2 ) = DN_De( i, 1 ) * ge[2] + DN_De( i, 0 ) * gn[2];
    }

    B = prod( Q, b );

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::CalculateJ(
    boost::numeric::ublas::bounded_matrix<double, 2, 2>& j,
    array_1d<double, 3>& ge,
    array_1d<double, 3>& gn,
    array_1d<double, 3>& v3 )

{
    double Norm_v3 = norm_2( v3 );
    double Norm_ge = norm_2( ge );

    double j11 = Norm_ge;
    double j12 = inner_prod( ge, gn ) / Norm_ge;
    double j22 = Norm_v3 / Norm_ge;

    j( 0, 0 ) = j11;
    j( 0, 1 ) = j12;
    j( 1, 0 ) = 0.00;
    j( 1, 1 ) = j22;
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::CalculateStrain(
    Vector& StrainVector,
    boost::numeric::ublas::bounded_matrix<double, 2, 2>& C )

{
    KRATOS_TRY

    StrainVector[0] = 0.5 * ( C( 0, 0 ) - 1.00 );
    StrainVector[1] = 0.5 * ( C( 1, 1 ) - 1.00 );
    StrainVector[2] = C( 0, 1 );

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::CalculateAndAdd_BodyForce(
    const Vector& N,
    const ProcessInfo& rCurrentProcessInfo,
    array_1d<double, 3>& BodyForce,
    VectorType& rRightHandSideVector,
    double weight )

{
    KRATOS_TRY

    unsigned int number_of_nodes = this->GetGeometry().size();
    const double density = this->GetProperties()[DENSITY];


    noalias(BodyForce) = ZeroVector(3);
    for(unsigned int i=0; i<number_of_nodes; i++)
        BodyForce += N[i]*this->GetGeometry()[i].FastGetSolutionStepValue(VOLUME_ACCELERATION);
    BodyForce *= density;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int index = 3 * i;

        for ( unsigned int j = 0; j < 3; j++ )
            rRightHandSideVector[index+j] += weight * N[i] * BodyForce[j];
    }

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::CalculateAndAdd_PressureForce(
    VectorType& residualvector,
    const Vector& N,
    const array_1d<double, 3>& v3,
    double pressure,
    double weight,
    const ProcessInfo& rCurrentProcessInfo )

{
    KRATOS_TRY

    unsigned int number_of_nodes = GetGeometry().size();

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int index = 3 * i;
        double coeff = pressure * N[i] * weight;
        residualvector[index]   += coeff * v3[0];
        residualvector[index+1] += coeff * v3[1];
        residualvector[index+2] += coeff * v3[2];
    }

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    bool CalculateStiffnessMatrixFlag,
    bool CalculateResidualVectorFlag )

{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes * 3;

    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRAIN, false);
    Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
    Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);


    Matrix B( 3, MatSize );
    Vector StrainVector( 3 );
    Vector StressVector( 3 );
    boost::numeric::ublas::bounded_matrix<double, 2, 2>  C = ZeroMatrix( 2, 2 );
    Matrix D = ZeroMatrix( 3, 3 );
    boost::numeric::ublas::bounded_matrix<double, 3, 3>  Q = ZeroMatrix( 3, 3 );

    //resizing as needed the LHS
    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }

    //resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size() != MatSize )
            rRightHandSideVector.resize( MatSize );

        rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS
    }

    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();

    const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer = GetGeometry().ShapeFunctionsLocalGradients();

    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues();

    //calculating actual jacobian
    GeometryType::JacobiansType J;

    J = GetGeometry().Jacobian( J );

    //auxiliary terms
    array_1d<double, 3> BodyForce;
    array_1d<double, 3> ge;
    array_1d<double, 3> gn;
    array_1d<double, 3> v3;
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        double IntegrationWeight = GetGeometry().IntegrationPoints()[PointNumber].Weight();

        ge[0] = J[PointNumber]( 0, 0 );
        gn[0] = J[PointNumber]( 0, 1 );
        ge[1] = J[PointNumber]( 1, 0 );
        gn[1] = J[PointNumber]( 1, 1 );
        ge[2] = J[PointNumber]( 2, 0 );
        gn[2] = J[PointNumber]( 2, 1 );

        CrossProduct( v3, ge, gn );
        boost::numeric::ublas::bounded_matrix<double, 2, 2> j;
        CalculateJ( j, ge, gn, v3 );

        // calculation of matrix g = jtrans*j;
        boost::numeric::ublas::bounded_matrix<double, 2, 2> g;
        noalias( g ) = prod( trans( j ), j );

        // calculation of the Right Cauchy-Green Tensor C = Gtrans*g*G
        boost::numeric::ublas::bounded_matrix<double, 2, 2> tmp;
        tmp = prod( g, mG_Vector[PointNumber] );
        noalias( C ) = prod( trans( mG_Vector[PointNumber] ), tmp );

        // calculation of the StrainVector
        CalculateStrain( StrainVector, C );
        mStrainsVector[PointNumber] = StrainVector; //saving the strain vector

        // if strain has to be computed inside of the constitutive law with PK2
        //rValues.SetDeformationGradientF(rVariables.F); //in this case F is the whole deformation gradient

        Values.SetStrainVector(StrainVector); //this is the input parameter
        Values.SetStressVector(StressVector); //this is an ouput parameter
        Values.SetConstitutiveMatrix(D);      //this is an ouput parameter

        mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(Values,ConstitutiveLaw::StressMeasure_PK2 );

        noalias( mStressesVector[PointNumber] ) = ZeroVector( 6 );
        Calculate_GlobalStressVector( mStressesVector[PointNumber], StressVector, mV1[PointNumber], mV2[PointNumber] ); //saving the stress vector

        CalculateQ( Q, mG_Vector[PointNumber] );
        CalculateB( B, Q, DN_DeContainer[PointNumber], ge, gn );

        // integration on the REFERENCE CONFIGURATION
        double DetJ0 = mDetJ0[PointNumber];
        double IntToReferenceWeight = IntegrationWeight * DetJ0 * mThickness0;

        // LEFT HAND SIDE MATRIX
        if ( CalculateStiffnessMatrixFlag == true )
        {
            //adding contributions to the stiffness matrix
            CalculateAndAddKm( rLeftHandSideMatrix, B, D, IntToReferenceWeight );
            CalculateAndAddKg( rLeftHandSideMatrix, Q, DN_DeContainer[PointNumber], StressVector, IntToReferenceWeight );
        }

//                      if(this->Id() == 51) //TODO: remove this! it is just for debugging purposes
//         {
//             KRATOS_WATCH(IntToReferenceWeight)
//             KRATOS_WATCH(DetJ0);
//             KRATOS_WATCH(mThickness0);
//             KRATOS_WATCH(StrainVector)
//             KRATOS_WATCH(StressVector)
//         }

        // RIGHT HAND SIDE VECTOR
        if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
        {
            //contribution to external forces
            //BodyForce = GetProperties()[BODY_FORCE];
            // operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
            CalculateAndAdd_BodyForce( row( Ncontainer, PointNumber ), rCurrentProcessInfo, BodyForce, rRightHandSideVector, IntToReferenceWeight );
            // operation performed: rRightHandSideVector -= Weight*IntForce
            noalias( rRightHandSideVector ) -= IntToReferenceWeight * prod( trans( B ), StressVector );
        }
    }
    // if(this->Id() == 51) //TODO: remove this! it is just for debugging purposes
    //   {
    //         KRATOS_WATCH(rLeftHandSideMatrix)
    //         KRATOS_WATCH(rRightHandSideVector)
    //         Vector displacements;
    //         this->GetValuesVector(displacements,0);
    //         KRATOS_WATCH( displacements );

    //         KRATOS_WATCH("coordintates")
    //         for(unsigned int i=0; i<GetGeometry().size(); i++)
    //             std::cout << " " << GetGeometry()[i].Id() << " " << GetGeometry()[i].Coordinates() << std::endl;

    //         KRATOS_WATCH("initial coordintates")
    //         for(unsigned int i=0; i<GetGeometry().size(); i++)
    //             std::cout << " " << GetGeometry()[i].Id() << " " << GetGeometry()[i].GetInitialPosition() << std::endl;
    //     }

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void MembraneElement::Calculate_GlobalStressVector(
    array_1d<double, 6>& GlobalVector,
    Vector& LocalStressVector,
    array_1d<double, 3>& v1,
    array_1d<double, 3>& v2 )

{
    KRATOS_TRY

    array_1d<double, 6> temp;

    //adding the component S11
    noalias( temp )  = VoigtTensorComponents( v1, v1 );
    temp *= LocalStressVector[0];
    noalias( GlobalVector ) += temp;

    //adding the component S22
    noalias( temp )  = VoigtTensorComponents( v2, v2 );
    temp *= LocalStressVector[1];
    noalias( GlobalVector ) += temp;

    //adding the component S12 (& S21)
    noalias( temp )  = VoigtTensorComponents( v1, v2 );
    noalias( temp ) += VoigtTensorComponents( v2, v1 );
    temp *= LocalStressVector[2];
    noalias( GlobalVector ) += temp;

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

//auxiliary function needed in the calculation of output stresses
inline array_1d<double, 6> MembraneElement::VoigtTensorComponents(
    array_1d<double, 3>& a,
    array_1d<double, 3>& b )

{
    array_1d<double, 6> v;

    v[0] = a[0] * b[0];
    v[1] = a[1] * b[1];
    v[2] = a[2] * b[2];
    v[3] = a[0] * b[1];
    v[4] = a[1] * b[2];
    v[5] = a[0] * b[2];

    return v;
}

//************************************************************************************
//************************************************************************************
void MembraneElement::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable,
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
    // VM
    if(rVariable==CAUCHY_STRESS_TENSOR)
    {
        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
    }
    // VM
}

int  MembraneElement::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //verify that the variables are correctly initialized

    if ( VELOCITY.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "VELOCITY has Key zero! (check if the application is correctly registered", "" )

        if ( DISPLACEMENT.Key() == 0 )
            KRATOS_THROW_ERROR( std::invalid_argument, "DISPLACEMENT has Key zero! (check if the application is correctly registered", "" )

            if ( ACCELERATION.Key() == 0 )
                KRATOS_THROW_ERROR( std::invalid_argument, "ACCELERATION has Key zero! (check if the application is correctly registered", "" )

                if ( DENSITY.Key() == 0 )
                    KRATOS_THROW_ERROR( std::invalid_argument, "DENSITY has Key zero! (check if the application is correctly registered", "" )

                    if ( VOLUME_ACCELERATION.Key() == 0 )
                        KRATOS_THROW_ERROR( std::invalid_argument, "BODY_FORCE has Key zero! (check if the application is correctly registered", "" )

                        if ( THICKNESS.Key() == 0 )
                            KRATOS_THROW_ERROR( std::invalid_argument, "THICKNESS has Key zero! (check if the application is correctly registered", "" )

                            //verify that the dofs exist
                            for ( unsigned int i = 0; i < this->GetGeometry().size(); i++ )
                            {
                                if ( this->GetGeometry()[i].SolutionStepsDataHas( DISPLACEMENT ) == false )
                                    KRATOS_THROW_ERROR( std::invalid_argument, "missing variable DISPLACEMENT on node ", this->GetGeometry()[i].Id() )

                                    if ( this->GetGeometry()[i].HasDofFor( DISPLACEMENT_X ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Y ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Z ) == false )
                                        KRATOS_THROW_ERROR( std::invalid_argument, "missing one of the dofs for the variable DISPLACEMENT on node ", GetGeometry()[i].Id() )
                                    }

    //verify that the constitutive law exists
    if ( this->GetProperties().Has( CONSTITUTIVE_LAW ) == false )
    {
        KRATOS_THROW_ERROR( std::logic_error, "constitutive law not provided for property ", this->GetProperties().Id() )
    }

    //verify that the constitutive law has the correct dimension
    if ( this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize() != 3 )
        KRATOS_THROW_ERROR( std::logic_error, "wrong constitutive law used. This is a 3D element with expected strain size is 3 (el id = ) ", this->Id() )

        //check constitutive law
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        {
            mConstitutiveLawVector[i]->Check( GetProperties(), GetGeometry(), rCurrentProcessInfo );

            ConstitutiveLaw::Features LawFeatures;
            mConstitutiveLawVector[i]->GetLawFeatures(LawFeatures);

            if(LawFeatures.mOptions.IsNot(ConstitutiveLaw::PLANE_STRESS_LAW))
                KRATOS_THROW_ERROR( std::logic_error,"Constitutive law is compatible only with a plane stress 2D law for membrane element with Id", this->Id())

                if(LawFeatures.mOptions.IsNot(ConstitutiveLaw::INFINITESIMAL_STRAINS))
                    KRATOS_THROW_ERROR( std::logic_error,"Constitutive law is compatible only with a law using infinitessimal strains for membrane element with Id", this->Id())

                    if(LawFeatures.mStrainSize != 3) KRATOS_THROW_ERROR( std::logic_error,"Constitutive law expects a strain size different from 3 for membrane element with Id", this->Id() )
                    }

    return 0;

    KRATOS_CATCH( "" );
}


//***********************************************************************************
//***********************************************************************************

} // Namespace Kratos.
