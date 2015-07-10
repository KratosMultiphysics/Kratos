/*
==============================================================================
KratosStructuralApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last modified by:    $Author: virginia $
//   Date:                $Date: 2009-01-23 14:39:59 $
//   Revision:            $Revision: 1.27 $
//
//


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/linear_element.h"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "structural_application.h"

//#include <omp.h>

namespace Kratos
{
//************************************************************************************
//************************************************************************************
LinearElement::LinearElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : Element( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
LinearElement::LinearElement( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : Element( NewId, pGeometry, pProperties )
{
//         const unsigned int dim = GetGeometry().WorkingSpaceDimension();
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();


}

Element::Pointer LinearElement::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new LinearElement( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

LinearElement::~LinearElement()
{
}

std::string LinearElement::Info() const
{
    std::stringstream buffer;
    buffer << "Linear Element" << std::endl;
    return buffer.str();
}


//************************************************************************************
//************************************************************************************
void LinearElement::Initialize()
{
    KRATOS_TRY

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    //resizing jacobian inverses containers
    mInvJ0.resize( integration_points.size() );
    mDetJ0.resize( integration_points.size(), false );


    GeometryType::JacobiansType J0;
    J0 = GetGeometry().Jacobian( J0, mThisIntegrationMethod );
    mTotalDomainInitialSize = 0.00;

    //Constitutive Law initialisation

    if ( mConstitutiveLawVector.size() != integration_points.size() )
    {
        mConstitutiveLawVector.resize( integration_points.size() );
        //InitializeMaterial();
    }

    InitializeMaterial();

    //calculating the inverse J0

    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        //getting informations for integration
        double IntegrationWeight = integration_points[PointNumber].Weight();

        //calculating and storing inverse of the jacobian and the parameters needed
        MathUtils<double>::InvertMatrix( J0[PointNumber], mInvJ0[PointNumber], mDetJ0[PointNumber] );

        //calculating the total area
        mTotalDomainInitialSize += mDetJ0[PointNumber] * IntegrationWeight;
    }


    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
void LinearElement::CalculateAll( MatrixType& rLeftHandSideMatrix,
                                  VectorType& rRightHandSideVector,
                                  ProcessInfo& rCurrentProcessInfo,
                                  bool CalculateStiffnessMatrixFlag,
                                  bool CalculateResidualVectorFlag )
{
    KRATOS_TRY
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();

    unsigned int StrainSize;

    if ( dim == 2 )
        StrainSize = 3;
    else
        StrainSize = 6;

    Matrix B( StrainSize, number_of_nodes*dim );

    Matrix F( dim, dim );

    Matrix D( StrainSize, StrainSize );

    Matrix C( dim, dim );

    Vector StrainVector( StrainSize );

    Vector StressVector( StrainSize );

    Matrix DN_DX( number_of_nodes, dim );


    //resizing as needed the LHS
    unsigned int MatSize = number_of_nodes * dim;

    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }

    //resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size() != MatSize )
            rRightHandSideVector.resize( MatSize, false );

        rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS
    }

    if ( DN_DX.size1() != GetGeometry().size() )
    {
        DN_DX.resize( GetGeometry().size(), dim, false );
    }

    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );

    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

    //calculating actual jacobian
    GeometryType::JacobiansType J;

    GetGeometry().Jacobian( J );

    //KRATOS_WATCH(J)

    //auxiliary terms
    Vector BodyForce;

    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        //Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
//    DN_DX = ZeroMatrix(DN_DX.size1(),DN_DX.size2());
        noalias( DN_DX ) = prod( DN_De[PointNumber], mInvJ0[PointNumber] );
        noalias( F ) = prod( J[PointNumber], mInvJ0[PointNumber] );
        identity_matrix<double> unit( dim );

        //strain calculation
        noalias( C ) =  F - unit;
        CalculateStrain( C, StrainVector );
        Comprobate_State_Vector( StrainVector );


        //material update (considering the level of strain achieved)
//    mConstitutiveLawVector[PointNumber]->UpdateMaterial( StrainVector,
//      GetProperties(),
//      GetGeometry(),
//      row(Ncontainer,PointNumber),
//      rCurrentProcessInfo );
        //Calculation of stress
        //KRATOS_WATCH( Id() )
        //KRATOS_WATCH(PointNumber)
        //mConstitutiveLawVector[PointNumber]->CalculateStress(StrainVector,StressVector);
        mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(
            StrainVector,
            ZeroMatrix( 1 ),
            StressVector,
            D,
            rCurrentProcessInfo,
            GetProperties(),
            GetGeometry(),
            row( Ncontainer, PointNumber ),
            true,
            CalculateStiffnessMatrixFlag,
            true );

        //calculating operator B
        CalculateB( B, DN_DX, StrainVector.size() );

        //calculating weights for integration on the reference configuration
        double IntToReferenceWeight = integration_points[PointNumber].Weight() * mDetJ0[PointNumber];

        if ( dim == 2 ) IntToReferenceWeight *= GetProperties()[THICKNESS];

        if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
        {
            //mConstitutiveLawVector[PointNumber]->CalculateConstitutiveMatrix(StrainVector,D);
            //mConstitutiveLawVector[PointNumber]->CalculateStressAndTangentMatrix(StressVector,StrainVector,D);
            //contributions to stiffness matrix calculated on the reference config
            noalias( rLeftHandSideMatrix ) += prod( trans( B ), ( IntToReferenceWeight ) * Matrix( prod( D, B ) ) ); //to be optimized to remove the temporary
            //KRATOS_WATCH(D)

            CalculateAndAddKg( rLeftHandSideMatrix, DN_DX, StressVector, IntToReferenceWeight );
        }

        if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
        {
            //contribution to external forces
            BodyForce = GetProperties()[BODY_FORCE];

            // operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
            CalculateAndAdd_ExtForceContribution( row( Ncontainer, PointNumber ), rCurrentProcessInfo, BodyForce, rRightHandSideVector, IntToReferenceWeight );

            // operation performed: rRightHandSideVector -= IntForce*IntToReferenceWeight
            noalias( rRightHandSideVector ) -= IntToReferenceWeight * prod( trans( B ), StressVector );
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
void LinearElement::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
}

//************************************************************************************
//************************************************************************************
void LinearElement::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;

    CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );


}

//************************************************************************************
//************************************************************************************
double LinearElement::CalculateIntegrationWeight( double GaussPointWeight, double DetJ0 )
{
    //to permorm the integration over the reference domain we need to include
    // the thickness in 2D
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    double weight = GaussPointWeight;

    weight *= DetJ0;

    if ( dimension == 2 ) weight *= GetProperties()[THICKNESS];

    return weight;
}

////************************************************************************************
////************************************************************************************

void LinearElement::InitializeSolutionStep( ProcessInfo& CurrentProcessInfo )
{
    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        mConstitutiveLawVector[i]->InitializeSolutionStep( GetProperties(),
                GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ),
                CurrentProcessInfo );
}

////************************************************************************************
////************************************************************************************
void LinearElement::FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo )
{
//         std::cout << "in TL: calling FinalizeSolutionStep" << std::endl;
    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        mConstitutiveLawVector[i]->FinalizeSolutionStep( GetProperties(),
                GetGeometry(),
                row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ),
                CurrentProcessInfo );
}

//************************************************************************************
//************************************************************************************
void LinearElement::InitializeMaterial()
{
    KRATOS_TRY
    
    if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        {
            mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
            // const std::string type = mConstitutiveLawVector[i]->Info();
            // KRATOS_WATCH(Id())
            // KRATOS_WATCH(Info())
            // KRATOS_WATCH(type)
            //KRATOS_WATCH(mConstitutiveLawVector[i])
        }
    }
    else
        KRATOS_THROW_ERROR( std::logic_error, "a constitutive law needs to be specified for the element with ID ", this->Id() )
        KRATOS_CATCH( "" )
    }

//************************************************************************************
//************************************************************************************
inline void LinearElement::CalculateAndAdd_ExtForceContribution(
    const Vector& N,
    const ProcessInfo& CurrentProcessInfo,
    Vector& BodyForce,
    VectorType& rRightHandSideVector,
    double weight
)
{
    KRATOS_TRY
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int index = dimension * i;

        for ( unsigned int j = 0; j < dimension; j++ )  rRightHandSideVector[index+j] += weight * N[i] * BodyForce[j];
    }

    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************
void LinearElement::CalculateAndAddKg(
    MatrixType& K,
    Matrix& DN_DX,
    Vector& StressVector,
    double weight )
{
    KRATOS_TRY
    // unsigned int dimension = mpReferenceGeometry->WorkingSpaceDimension();
    //Matrix<double> StressTensor = MathUtils<double>::StressVectorToTensor(StressVector);
    //Matrix<double> ReducedKg(DN_Dx.RowsNumber(),DN_Dx.RowsNumber());
    //Matrix<double>::MatMulAndAdd_B_D_Btrans(ReducedKg,weight,DN_Dx,StressTensor);
    //MathUtils<double>::ExpandAndAddReducedMatrix(K,ReducedKg,dimension);

    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    Matrix StressTensor = MathUtils<double>::StressVectorToTensor( StressVector );
    Matrix ReducedKg = prod( DN_DX, weight * Matrix( prod( StressTensor, trans( DN_DX ) ) ) ); //to be optimized
    MathUtils<double>::ExpandAndAddReducedMatrix( K, ReducedKg, dimension );

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
void LinearElement::CalculateStrain(
    const Matrix& C,
    Vector& StrainVector )
{
    KRATOS_TRY
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    Matrix Aux( dimension, dimension );
    noalias( Aux ) = 0.50 * ( C + trans( C ) );


    if ( dimension == 2 )
    {
        if ( StrainVector.size() != 3 ) StrainVector.resize( 3, false );

        StrainVector[0] = Aux( 0, 0 );

        StrainVector[1] = Aux( 1, 1 );

        StrainVector[2] = 2.00 * Aux( 0, 1 );
    }

    if ( dimension == 3 )
    {
        if ( StrainVector.size() != 6 ) StrainVector.resize( 6, false );

        StrainVector[0] = Aux( 0, 0 );

        StrainVector[1] = Aux( 1, 1 );

        StrainVector[2] = Aux( 2, 2 );

        StrainVector[3] = 2.00 * Aux( 0, 1 );

        StrainVector[4] = 2.00 * Aux( 1, 2 );

        StrainVector[5] = 2.00 * Aux( 0, 2 );
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
void LinearElement::CalculateB(
    Matrix& B,
    //Matrix& F,
    Matrix& DN_DX,
    unsigned int StrainSize )
{
    KRATOS_TRY
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //unsigned int dim2 = number_of_nodes*dimension;
    //if(B.size1() != StrainSize || B.size2()!=dim2)
    // B.resize(StrainSize,dim2);
    //Matrix Bi;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = dimension * i;

        if ( dimension == 2 )
        {
            identity_matrix<double> F( 2 );
            B( 0, index + 0 ) = F( 0, 0 ) * DN_DX( i, 0 );
            B( 0, index + 1 ) = F( 1, 0 ) * DN_DX( i, 0 );
            B( 1, index + 0 ) = F( 0, 1 ) * DN_DX( i, 1 );
            B( 1, index + 1 ) = F( 1, 1 ) * DN_DX( i, 1 );
            B( 2, index + 0 ) = F( 0, 0 ) * DN_DX( i, 1 ) + F( 0, 1 ) * DN_DX( i, 0 );
            B( 2, index + 1 ) = F( 1, 0 ) * DN_DX( i, 1 ) + F( 1, 1 ) * DN_DX( i, 0 );
        }
        else
        {
            identity_matrix<double> F( 3 );
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

        //CalculateBi(Bi,F,DN_DX,i);
        //MathUtils<double>::WriteMatrix(B,Bi,0,index);
    }

    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************
void LinearElement::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo )
{
    int number_of_nodes = GetGeometry().size();
    int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int dim2 = number_of_nodes * dim;

    if ( rResult.size() != dim2 )
        rResult.resize( dim2, false );

    for ( int i = 0; i < number_of_nodes; i++ )
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
void LinearElement::GetDofList( DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo )
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
void LinearElement::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //lumped
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int NumberOfNodes = GetGeometry().size();
    unsigned int MatSize = dimension * NumberOfNodes;

    if ( rMassMatrix.size1() != MatSize )
        rMassMatrix.resize( MatSize, MatSize, false );

    rMassMatrix = ZeroMatrix( MatSize, MatSize );

    double TotalMass = mTotalDomainInitialSize * GetProperties()[DENSITY];

    if ( dimension == 2 ) TotalMass *= GetProperties()[THICKNESS];

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
void LinearElement::CalculateDampingMatrix( MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
    unsigned int MatSize = number_of_nodes * dim;

    if ( rDampingMatrix.size1() != MatSize )
        rDampingMatrix.resize( MatSize, MatSize, false );

    noalias( rDampingMatrix ) = ZeroMatrix( MatSize, MatSize );

    KRATOS_CATCH( "" )
}

LinearElement::IntegrationMethod LinearElement::GetIntegrationMethod() const
{
    return mThisIntegrationMethod;
}

//************************************************************************************
//************************************************************************************
void LinearElement::CalculateOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& Output, const ProcessInfo& rCurrentProcessInfo )
{
    if ( Output.size() != GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size() )
        Output.resize( GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size() );

    for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
        Output[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, Output[ii] );
}

//************************************************************************************
//************************************************************************************
void LinearElement::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& Output, const ProcessInfo& rCurrentProcessInfo )
{
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();

    unsigned int StrainSize;

    if ( dim == 2 )
        StrainSize = 3;
    else
        StrainSize = 6;

    Vector StrainVector( StrainSize );

    if ( Output.size() != GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size() )
        Output.resize( GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size() );

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
        for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
            Output[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, Output[ii] );
    }
}

//************************************************************************************
//************************************************************************************
void LinearElement::CalculateOnIntegrationPoints( const Variable<Matrix >& rVariable, std::vector< Matrix >& Output, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();

    unsigned int StrainSize;

    if ( dim == 2 )
        StrainSize = 3;
    else
        StrainSize = 6;

    Matrix F( dim, dim );

    Matrix D( StrainSize, StrainSize );

    Matrix C( dim, dim );

    Vector StrainVector( StrainSize );

    Vector StressVector( StrainSize );

    Matrix DN_DX( number_of_nodes, dim );

    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    //calculating actual jacobian
    GeometryType::JacobiansType J;

    J = GetGeometry().Jacobian( J );

    if ( Output.size() != integration_points.size() )
        Output.resize( integration_points.size() );

    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {

        //deformation gradient
        noalias( F ) = prod( J[PointNumber], mInvJ0[PointNumber] );
        //KRATOS_WATCH(Id());
        //KRATOS_WATCH(F);
        //KRATOS_WATCH("*************************");



        //strain calculation
        identity_matrix<double> unit( dim );
        noalias( C ) =  F - unit;
        CalculateStrain( C, StrainVector );
        Comprobate_State_Vector( StrainVector );

        if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
        {
            if ( Output[PointNumber].size2() != StrainVector.size() )
                Output[PointNumber].resize( 1, StrainVector.size(), false );

            for ( unsigned int ii = 0; ii < StrainVector.size(); ii++ )
                Output[PointNumber]( 0, ii ) = StrainVector[ii];
        }
        else if ( rVariable == PK2_STRESS_TENSOR )
        {
            if ( Output[PointNumber].size2() != StrainVector.size() )
                Output[PointNumber].resize( 1, StrainVector.size(), false );

            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(
                StrainVector,
                F,
                StressVector,
                D,
                rCurrentProcessInfo,
                GetProperties(),
                GetGeometry(),
                row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), PointNumber ),
                true,
                0,
                true );

            for ( unsigned int ii = 0; ii < StrainVector.size(); ii++ )
            {
                Output[PointNumber]( 0, ii ) = StressVector[ii];
            }
        }

        else if ( rVariable == GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR )
        {
            Matrix PlasticStrainVector( 1, StrainVector.size() );

            mConstitutiveLawVector[PointNumber]->GetValue( GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR, PlasticStrainVector );

            Output[PointNumber] = PlasticStrainVector;
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
void LinearElement::SetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    for ( unsigned int PointNumber = 0; PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(); PointNumber++ )
    {
        mConstitutiveLawVector[PointNumber]->SetValue( rVariable ,
                rValues[PointNumber], rCurrentProcessInfo );
    }

}


//************************************************************************************
//************************************************************************************
void LinearElement::SetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    for ( unsigned int PointNumber = 0; PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(); PointNumber++ )
    {
        mConstitutiveLawVector[PointNumber]->SetValue( rVariable ,
                rValues[PointNumber], rCurrentProcessInfo );
    }

}

//************************************************************************************
//************************************************************************************
void LinearElement::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo )
{
    if ( rValues.size() != GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size() )
        rValues.resize( GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(), false );

    for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
        rValues[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rValues[ii] );
}


//************************************************************************************
//************************************************************************************

void LinearElement::GetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
    if ( rVariable == INSITU_STRESS )
    {
        for ( unsigned int PointNumber = 0;
                PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size();
                PointNumber++ )
        {
            rValues[PointNumber] =
                mConstitutiveLawVector[PointNumber]->GetValue( INSITU_STRESS, rValues[PointNumber] );
        }
    }

    if ( rVariable == MATERIAL_PARAMETERS )
    {
        for ( unsigned int PointNumber = 0;
                PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(); PointNumber++ )
        {
            rValues[PointNumber] =
                mConstitutiveLawVector[PointNumber]->GetValue( MATERIAL_PARAMETERS, rValues[PointNumber] );
        }
    }

    if ( rVariable == INTERNAL_VARIABLES )
    {
        for ( unsigned int PointNumber = 0;
                PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size();
                PointNumber++ )
        {
            rValues[PointNumber] =
                mConstitutiveLawVector[PointNumber]->GetValue( INTERNAL_VARIABLES, rValues[PointNumber] );

        }
    }


}

//************************************************************************************
//************************************************************************************
void LinearElement::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable,
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

    if ( rVariable == GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR )
    {
        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
    }

}

//************************************************************************************
//************************************************************************************
void LinearElement::GetValuesVector( Vector& values, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes * dim;

    if ( values.size() != MatSize ) values.resize( MatSize, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dim;
        values[index] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
        values[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );

        if ( dim == 3 )
            values[index + 2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );
    }
}


//************************************************************************************
//************************************************************************************
void LinearElement::GetFirstDerivativesVector( Vector& values, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes * dim;

    if ( values.size() != MatSize )   values.resize( MatSize, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dim;
        values[index] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
        values[index + 1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );

        if ( dim == 3 )
            values[index + 2] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Z, Step );
    }
}

//************************************************************************************
//************************************************************************************
void LinearElement::GetSecondDerivativesVector( Vector& values, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes * dim;

    if ( values.size() != MatSize ) values.resize( MatSize, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dim;
        values[index] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
        values[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );

        if ( dim == 3 )
            values[index + 2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );
    }
}

//************************************************************************************
//************************************************************************************

void LinearElement::Calculate( const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo )
{

    double c =  0.00; //sqrt(GetProperties()[YOUNG_MODULUS]/GetProperties()[DENSITY]);
    Vector Values( GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size() );
    //KRATOS_WATCH(Values.size())
    //KRATOS_WATCH(GetGeometry().IntegrationPoints(mThisIntegrationMethod).size())

    if ( rVariable == DELTA_TIME )
    {
        for ( unsigned int PointNumber = 0;
                PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size();
                PointNumber++ )
        {
            mConstitutiveLawVector[PointNumber]-> GetValue( DELTA_TIME, c );
            Values[PointNumber] = c;
        }
    }

    c = ( *std::max_element( Values.begin(), Values.end() ) );

    //KRATOS_WATCH(Id())
    //KRATOS_WATCH(c)
    double le =  GetGeometry().Length();
    //KRATOS_WATCH(le)

    Output          =  le / c;
    //KRATOS_WATCH(Output)

}

//************************************************************************************
//************************************************************************************

void  LinearElement::Comprobate_State_Vector( Vector& Result )
{
    for ( unsigned int i = 0.00; i < Result.size(); i++ )
    {
        if ( fabs( Result( i ) ) < 1E-9 )
        {
            Result( i ) = 0.00;
        }
    }
}


} // Namespace Kratos


