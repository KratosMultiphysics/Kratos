/*
== == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==
MengmengApplication
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

== == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==
*/
//
//   Project Name:        Kratos
//   Last modified by:    $Author: mengmeng $
//   Date:                $Date: 2010-01-06 10:23:55 $
//   Revision:            $Revision: 1.2 $
//
//
/// Default dimension: u[mm], -> pL[kPa], T[°C],
/// In Problem type Soil: The units have been converted.
/// Rescaling is needed in the GID when create the model geometry and assign unit-given parameters: 1 m --> 1000 mm

// System includes
#include <math.h>

// 20130828
//  Assume: vL[Grad_p, t], comment DvL_Dp and DphiM_Dp
// 

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/solid.h"
#include "utilities/math_utils.h"
//#include "includes/constitutive_law.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/tetrahedra_3d_4.h"
#include "freezing_soil_application.h"
#include "freezing_soil.h"
#include "structural_application/custom_utilities/sd_math_utils.h"

namespace Kratos
{

Solid::Solid( IndexType NewId, GeometryType::Pointer pGeometry )
        : Element( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}


Solid::Solid( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
        : Element( NewId, pGeometry, pProperties )
{
    mNodesDispMin = 1; // Disp: displacement 
    mNodesDispMax = GetGeometry().size();

    if ( mNodesDispMax == 8 || mNodesDispMax == 27 || mNodesDispMax == 20 )
    { 
        mThisIntegrationMethod = GeometryData::GI_GAUSS_3;//methods for hexahedra elements
    } 
    //select Integration Method for 10 and 4 node tetrahedron
    else if (mNodesDispMax== 4 || mNodesDispMax== 10 )
    { 
        mThisIntegrationMethod = GeometryData::GI_GAUSS_2;//methods for tetrahedra elements
    }
    else
        KRATOS_THROW_ERROR( std::logic_error, "This element matches only with a linear hexahedra/tetrahedra (8/3) or quadratic hexahedra/tetrahedra (20,27/10 ) geometry" , *this ); 
}

Element::Pointer Solid::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new Solid( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

Solid::~Solid()
{
}


void Solid::ResetConstitutiveLaw()
{
    KRATOS_TRY
    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        mConstitutiveLawVector[i]->ResetMaterial( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
    KRATOS_CATCH( "" )
}


void Solid::InitializeSolutionStep( ProcessInfo& CurrentProcessInfo )
{
}

////****************************************************************
//************************************************************************************
//************************************************************************************

// Degree of freedom: displacment: u[m], water pressure p[Pa]
void Solid::Initialize()
{
    KRATOS_TRY

    mNodesNumberDisp = mNodesDispMax - mNodesDispMin + 1;
    mDimension = GetGeometry().WorkingSpaceDimension();
    mMatSizeU = mNodesNumberDisp * mDimension; 
    mNumberU = 1; 
    mMatSize = mMatSizeU * mNumberU; 
    mScaleU = GetProperties()[SCALE_U]; 
    mTol = 1.0e-8; 

//     mMaterialParameters = GetProperties()[ELEMENT_PARAMETERS];  
//     mUnitRatio = GetProperties()[SCALE];				// Unit convertor: from [m] to [mm]
    mUnitRatio = 1.0;				// Unit convertor: from [m] to [mm]
    mrho = GetProperties()[DENSITY] * pow( mUnitRatio, 3.0 ); 	// [ kg/m^3 ] 
//     mK = GetProperties()[BULK_MODULUS] / mUnitRatio; 			// [ Pa = kg/(m*s^2) ] 
//     mG = GetProperties()[SHEAR_MODULUS] / mUnitRatio; 			// [ Pa = kg/(m*s^2) ]
    mK = 10e+6  / mUnitRatio; 			// [ Pa = kg/(m*s^2) ] 
    mG = 3e+6  / mUnitRatio; 			// [ Pa = kg/(m*s^2) ]
    mGravity = GetProperties()[GRAVITY] * mUnitRatio; 			// [m/s^2] 
 
    // Needed for assigning nonzero Dirchlet bc: eg. Ti=5 or pLi= 1000;
    for ( unsigned int i = ( mNodesDispMin - 1 ) ; i < mNodesDispMax ; i++ )
    {
        GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_NULL ) = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT );
        GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_EINS ) = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT );
        GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_DT ) = ZeroVector( 3 );
        GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_NULL_DT ) = ZeroVector( 3 );
        GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_EINS_DT ) = ZeroVector( 3 );
        GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_OLD ) = ZeroVector( 3 );
        GetGeometry()[i].GetSolutionStepValue( ACCELERATION ) = ZeroVector( 3 );
        GetGeometry()[i].GetSolutionStepValue( ACCELERATION_NULL ) = ZeroVector( 3 );
        GetGeometry()[i].GetSolutionStepValue( ACCELERATION_EINS ) = ZeroVector( 3 );
    } 

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );
    GeometryType::JacobiansType J0( integration_points.size() );
    J0 = GetGeometry().Jacobian( J0, mThisIntegrationMethod );
    mInvJ0.resize( integration_points.size() );
    mDetJ0.resize( integration_points.size(), false );
    noalias( mDetJ0 ) = ZeroVector( integration_points.size() );
 
    for ( unsigned int GaussPoint = 0; GaussPoint < integration_points.size(); GaussPoint++ )
    { 
        // InvJ0
        mInvJ0[GaussPoint].resize( mDimension, mDimension, false );
        noalias( mInvJ0[GaussPoint] ) = ZeroMatrix( mDimension, mDimension );

        //calculating the inverse J0
        MathUtils<double>::InvertMatrix( J0[GaussPoint], mInvJ0[GaussPoint], mDetJ0[GaussPoint] );
    }

    // initialize material
    if ( mConstitutiveLawVector.size() != integration_points.size() )
        mConstitutiveLawVector.resize( integration_points.size() );
    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
    {
        mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
        mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
    }
 

    KRATOS_CATCH( "" )
}

void Solid::CalculateAll( MatrixType& rLeftHandSideMatrix,
                                 VectorType& rRightHandSideVector, ProcessInfo& rProcessInfo,
                                 bool CalculateStiffnessMatrixFlag, bool CalculateResidualVectorFlag )
{
    KRATOS_TRY
    //resizing as needed the RHS

    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size() != mMatSize )
            rRightHandSideVector.resize( mMatSize );
        noalias( rRightHandSideVector ) = ZeroVector( mMatSize );   //resetting RHS
    }

    //resizing as needed the LHS
    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != mMatSize )
            rLeftHandSideMatrix.resize( mMatSize, mMatSize );
        noalias( rLeftHandSideMatrix ) = ZeroMatrix( mMatSize, mMatSize );   //resetting LHS
    }

    //reading integration points and local Gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );
    const GeometryType::ShapeFunctionsGradientsType& DNu_Dxi = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );

    const Matrix& Nu_container = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ); 

    Vector Help_R_1( mMatSizeU ); 
    Matrix Help_K_UU( mMatSizeU, mMatSizeU );  

    // ++++++ RightHandSideVector: R ++++++
    if ( CalculateResidualVectorFlag == true )
    {
        noalias( Help_R_1 ) = ZeroVector( mMatSizeU ); 
    }

    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        noalias( Help_K_UU ) = ZeroMatrix( mMatSizeU, mMatSizeU ); 
    }

    double weight, DetJ = 0.0 ;
    Vector Nu( mNodesNumberDisp );
    Matrix DNu_DX( mNodesNumberDisp, mDimension );

    double Div_dot_u; 
    Matrix Grad_u( mDimension, mDimension );

    // adding for calling constitutive law
    Vector InputVector( 6 );
    Vector strainVector( 6 );

    Vector stressVector( 6 );
    noalias( stressVector ) = ZeroVector( 6 );
    Matrix CtanEff( 6, 6 );
    noalias( CtanEff ) = ZeroMatrix( 6, 6 );
     

    // START GAUSS INTEGRATION:
    for ( unsigned int GaussPoint = 0; GaussPoint < integration_points.size(); GaussPoint++ )
    {
        weight = integration_points[GaussPoint].Weight();
        DetJ = mDetJ0[GaussPoint];

        noalias( DNu_DX ) = prod( DNu_Dxi[GaussPoint], mInvJ0[GaussPoint] );
        noalias( Nu ) = row( Nu_container, GaussPoint );

        noalias( Grad_u ) = GetGradu( DNu_DX );
        Div_dot_u = GetDivdotu( DNu_DX );

        // calling constitutive law for computing effective stress,
        noalias( strainVector ) = GetstrainVector( Grad_u );   // Compressive ! 

        if ( GetValue( PLASTIC_FLAG ) == 0 ) // ELATIC 
        {
            /// METHOD 1: CtanEff[ K[t], G[t] ]
            noalias( CtanEff ) = GetElasticTangent(); 
            noalias( stressVector ) = GetstressVector( Grad_u );  
	    
	    if ( GetValue( KRATOS_WATCH_FLAG ) == 1 && GaussPoint == 1 )
	    {
		std::cout << "1---- u= "<<Getu(Nu)<<",\t CtanEff= " << CtanEff << std::endl; 
	    }
        }
        else // PLATIC 
        {
            /// METHOD 2: mBBM model
            for ( unsigned int i = 0; i < 6; i++ )
                InputVector[i] = strainVector[i]; 

            bool onlyGP1 = false;
	    if ( GetValue( KRATOS_WATCH_FLAG ) == 1 && GaussPoint == 1 )
	    {
                onlyGP1 = true;
		std::cout << "1---- u= "<<Getu(Nu)<<",\t Div_dot_u= " << Div_dot_u << std::endl; 
	    }
	      
            // INPUT: strainVector (Tensile-Positive !!);  OUTPUT: stressVector [kPa] (Tensile-Positive !!) , CtanEff [kPa]
            mConstitutiveLawVector[GaussPoint]->CalculateMaterialResponse( InputVector, ZeroMatrix( 1, 1 ), stressVector, CtanEff, rProcessInfo, GetProperties(), GetGeometry(), Nu, true, 1, onlyGP1 );
        } 
        
        if ( CalculateResidualVectorFlag == true ) 
            AddInternalForcesToRHS1( Help_R_1, weight, DetJ, Grad_u, Div_dot_u, Nu, DNu_DX, stressVector ); 
        if ( CalculateStiffnessMatrixFlag == true ) 
            CalculateStiffnessMatrixUU( Help_K_UU, weight, DetJ, Grad_u, Div_dot_u, Nu, DNu_DX, CtanEff );  

    }// END GAUSS INTEGRATION.

    // Assemble Kt and Ri : { ux_n1, uy_n1, uz_n1, ux_n2, uy_n2, uz_n2, ... n_n1, p_n1, t_n1, n_n2, p_n2, t_n2, ...}
    if ( CalculateStiffnessMatrixFlag == true )
    {
        for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
            for ( unsigned int k = 0; k < mDimension; k++ )
            {
                for ( unsigned int j = ( mNodesDispMin - 1 ); j < mNodesDispMax; j++ )
                    for ( unsigned int l = 0; l < mDimension; l++ )
                        rLeftHandSideMatrix( i*mDimension + k, j*mDimension + l ) += Help_K_UU( i * mDimension + k, j * mDimension + l );
            }
    }

    if ( CalculateResidualVectorFlag == true )
    {
        for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
            for ( unsigned int k = 0; k < mDimension; k++ )
                rRightHandSideVector( i*mDimension + k ) -= Help_R_1( i * mDimension + k );
    }

if ( GetValue( KRATOS_WATCH_FLAG ) == 1 )
  KRATOS_WATCH(Help_R_1)
  
    KRATOS_CATCH( "" )
}


void Solid::DampMatrix( MatrixType& rDampMatrix, ProcessInfo& rProcessInfo )
{
    KRATOS_TRY
 
    KRATOS_CATCH( "" )
}
///-//////////////////////////////////////////////////////////////////////////////////////
///-//////////////////////////////////////////////////////////////////////////////////////
/// R1 ***************************
void Solid::AddInternalForcesToRHS1( Vector& Help_R_1, double weight, double DetJ, Matrix& Grad_u, double Div_dot_u, const Vector& Nu, const Matrix& DNu_DX, Vector& stressVector )
{
    Matrix Bu( 6, mNodesNumberDisp*mDimension );
    noalias( Bu ) = GetBu( DNu_DX );
 
	
    for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
        for ( unsigned int m = 0; m < mDimension; m++ )
        {
            Help_R_1( i*mDimension + m ) += Nu( i ) * (
                                                mrho * mGravity( m )
                                            ) * weight * DetJ * mScaleU;
   
            for ( unsigned int k = 0; k < 6; k++ )
                Help_R_1( i*mDimension + m ) -= Bu( k, i * mDimension + m ) * (
                                                    stressVector( k )
                                                ) * weight * DetJ * mScaleU;
        }
}

void Solid::CalculateStiffnessMatrixUU( Matrix& Help_K_UU, double weight, double DetJ, Matrix& Grad_u, double Div_dot_u, const Vector& Nu, const Matrix& DNu_DX, Matrix& CtanEff )
{
    Matrix Bu( 6, mNodesNumberDisp*mDimension );
    noalias( Bu ) = GetBu( DNu_DX );

    for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
        for ( unsigned int m = 0; m < mDimension; m++ )
            for ( unsigned int j = ( mNodesDispMin - 1 ); j < mNodesDispMax; j++ )
                for ( unsigned int n = 0; n < mDimension; n++ )
                {  
                    for ( unsigned int k = 0; k < 6; k++ )
                        for ( unsigned int l = 0; l < 6; l++ )
                        {
                            Help_K_UU( i*mDimension + m, j*mDimension + n ) -= Bu( k, i * mDimension + m ) * (
                                        CtanEff( k, l )
                                    ) * Bu( l, j * mDimension + n ) * weight * DetJ * mScaleU;
                        }

                }
}

///-//////////////////////////////////////////////////////////////////////////////////////
///  CONSTITUTIVE QUANTITES
///-//////////////////////////////////////////////////////////////////////////////////////
/// KnoneckerDelta
double Solid::KnoneckerDelta( int i, int j )
{
    if ( i == j )
        return 1.0;
    else
        return 0.0;
}

/// Bu operator
Matrix Solid::GetBu( const Matrix& DNu_DX )
{
    Matrix result( 6, mNodesNumberDisp*mDimension );
    noalias( result ) = ZeroMatrix( 6, mNodesNumberDisp * mDimension );

    for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
    {
        result( 0, i*3 ) = DNu_DX( i, 0 );
        result( 1, i*3 + 1 ) = DNu_DX( i, 1 );
        result( 2, i*3 + 2 ) = DNu_DX( i, 2 );
        result( 3, i*3 ) = DNu_DX( i, 1 );
        result( 3, i*3 + 1 ) = DNu_DX( i, 0 );
        result( 4, i*3 + 1 ) = DNu_DX( i, 2 );
        result( 4, i*3 + 2 ) = DNu_DX( i, 1 );
        result( 5, i*3 ) = DNu_DX( i, 2 );
        result( 5, i*3 + 2 ) = DNu_DX( i, 0 );
    }
    return result;
}

///--------------------------------------------------------------------------
/// Displacement
///--------------------------------------------------------------------------
Vector Solid::Getu( const Vector& Nu )
{
    Vector result( mDimension );
    noalias( result ) = ZeroVector( mDimension );

    for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ ) 
    {
	result( 0 ) += GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X ) * Nu( i );
	result( 1 ) += GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y ) * Nu( i );
	result( 2 ) += GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z ) * Nu( i );
    }
  		
    return result;
}
Matrix Solid::GetGradu( const Matrix& DNu_DX )
{
    Matrix result( mDimension, mDimension );
    noalias( result ) = ZeroMatrix( mDimension, mDimension );

    for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
        for ( unsigned int k = 0; k < mDimension; k++ )
        {
            result( 0, k ) += GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X ) * DNu_DX( i, k );
            result( 1, k ) += GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y ) * DNu_DX( i, k );
            result( 2, k ) += GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z ) * DNu_DX( i, k );
        }

    return result;
}

double Solid::GetDivdotu( const Matrix& DNu_DX )
{
    double result = 0.0;

    for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
        for ( unsigned int k = 0; k < mDimension; k++ )
            result  += GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_DT )[k] * DNu_DX( i, k );

    return result;
}  

/// ++++++++++++++ Group c1 +++++++++++++++++++
/// +++++++++++++++++++++++++++++++++++++++++++
/// elastic modulus: CtanEff
Matrix Solid::GetElasticTangent()
{
    Matrix result( 6, 6 );
    noalias( result ) = ZeroMatrix( 6, 6 );
  
    double c1 = mK + 4.0 * mG / 3.0;
    double c2 = mK - 2.0 * mG / 3.0;
    for ( unsigned int i = 0; i < 3; i++ )
    {
	result( i, i ) = c1;
	result( i + 3, i + 3 ) = mG;

	for ( unsigned int j = 0; j < 3; j++ )
	    if ( i != j )
		result( i, j ) = c2;
    } 
    
    return result;
}
/// strain: linear strain tensor
Vector Solid::GetstrainVector( Matrix Grad_u )
{
    Vector result( 6 );
    noalias( result ) = ZeroVector( 6 );
    result( 0 ) = Grad_u( 0, 0 );
    result( 1 ) = Grad_u( 1, 1 );
    result( 2 ) = Grad_u( 2, 2 );
    result( 3 ) = Grad_u( 0, 1 ) + Grad_u( 1, 0 );
    result( 4 ) = Grad_u( 1, 2 ) + Grad_u( 2, 1 );
    result( 5 ) = Grad_u( 0, 2 ) + Grad_u( 2, 0 );

    return result;
}

/// strainVol: volumetric strain
double Solid::GetstrainVol( Matrix Grad_u )
{
    return  Grad_u( 0, 0 ) + Grad_u( 1, 1 ) + Grad_u( 2, 2 );
}


/// -------------- stress components ---------------
/// ------------------------------------------------ 
Vector Solid::GetstressVector( Matrix Grad_u )
{
    Vector result( 6 );
    noalias( result ) = ZeroVector( 6 ); 
    
    noalias( result ) = prod( GetElasticTangent(), GetstrainVector( Grad_u ) ); 
    return result;
}

double Solid::Getp ( Vector& stressVector )
{
    return ( stressVector[0] + stressVector[1] + stressVector[2] ) / 3.0;
}

double Solid::Getq ( Vector& stressVector )
{
//     return sqrt ( 1.5 ) * MathUtils<double>::Norm ( GetDevStress ( stress ) );
    return sqrt( ( pow(stressVector[0]-stressVector[1],2.0) +  pow(stressVector[1]-stressVector[2],2.0) + pow(stressVector[2]-stressVector[0],2.0) )/2.0 + 3.0*( stressVector[3]*stressVector[3] + stressVector[4]*stressVector[4] + stressVector[5]*stressVector[5] ) );
}
 
/// SOLUTION

Solid::IntegrationMethod Solid::GetIntegrationMethod()
{
    return mThisIntegrationMethod;
}

void Solid::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rProcessInfo )
{
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll( temp, rRightHandSideVector, rProcessInfo, CalculateStiffnessMatrixFlag,  CalculateResidualVectorFlag );
}


void Solid::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rProcessInfo )
{
    //calculation flamgS
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;

    CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
}

void Solid::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& ProcessInfo )
{
    if ( rResult.size() != mMatSize )
        rResult.resize( mMatSize );

    for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
    {
        rResult[i*mDimension]   = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[i*mDimension+1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();
        rResult[i*mDimension+2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
    } 
}


void Solid::GetDofList( DofsVectorType& ElementalDofList, ProcessInfo& ProcessInfo )
{
    ElementalDofList.resize( 0 );

    for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
    {
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
    } 
}


void Solid::GetValuesVector( Vector& values, int Step )
{
    if ( values.size() != mMatSize )
        values.resize( mMatSize );

    for ( unsigned int i = ( mNodesDispMin - 1 );i < mNodesDispMax;i++ )
    {
        values( i*mDimension )   = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
        values( i*mDimension + 1 ) = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );
        values( i*mDimension + 2 ) = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );
    } 
}

///-////////////////////////////////////
/////
///- //////////////////////////////////////////////////////
void Solid::GetValueOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rProcessInfo )
{
	KRATOS_TRY
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );
    if(rValues.size() != integration_points.size())
	    rValues.resize( integration_points.size() );
 
    const GeometryType::ShapeFunctionsGradientsType& DNu_Dxi = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod ); 
    const Matrix& Nu_container = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ); 

    Vector Nu( mNodesNumberDisp );
    Matrix DNu_DX( mNodesNumberDisp, mDimension );
    
    Vector u( mDimension ); // solid displacement 
    Matrix Grad_u( mDimension, mDimension ); 
    Vector dot_u( mDimension ); 
    double Div_dot_u; 
    Vector stressVector( mDimension ); 

    /////////////////////////////////////////////////////////////////////////
    //// Integration in space sum_(beta= 0)^(number of quadrature points)
    /////////////////////////////////////////////////////////////////////////

   for ( unsigned int GaussPoint = 0; GaussPoint < integration_points.size(); GaussPoint++ )
    {
        // in element
        noalias( DNu_DX ) = prod( DNu_Dxi[GaussPoint], mInvJ0[GaussPoint] ); 
        noalias( Nu ) = row( Nu_container, GaussPoint ); 

        noalias( u ) = Getu( Nu );
        noalias( Grad_u ) = GetGradu( DNu_DX );
        Div_dot_u = GetDivdotu( DNu_DX ); 
        noalias( stressVector ) = GetstressVector( Grad_u );
	
        if ( rVariable == LINEAR_STRAIN )
            rValues[GaussPoint] = GetstrainVol( Grad_u );
        else if ( rVariable == EQUIVALENT_VOLUMETRIC_STRESS )
            rValues[GaussPoint] = Getp ( stressVector );
        else if ( rVariable == EQUIVALENT_DEVIATORIC_STRESS )
            rValues[GaussPoint] = Getq ( stressVector );
        else // in constitutive law
            rValues[GaussPoint] = mConstitutiveLawVector[GaussPoint]->GetValue( rVariable, rValues[GaussPoint] );
    }

    KRATOS_CATCH( "" )
}
void Solid::GetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rProcessInfo )
{
    KRATOS_TRY

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );
    if(rValues.size() != integration_points.size())
	    rValues.resize( integration_points.size() );

    for ( unsigned int GaussPoint = 0; GaussPoint < integration_points.size(); GaussPoint++ )
    {  
	if ( rVariable == INTERNAL_VARIABLES ) 
	{
	    if ( rValues[GaussPoint].size() != 6 )
		rValues[GaussPoint].resize( 6 );
	    noalias( rValues[GaussPoint] )= mConstitutiveLawVector[GaussPoint]->GetValue( INTERNAL_VARIABLES, rValues[GaussPoint] ); 
	}

	//To Plot Stresses
	if ( rVariable == INSITU_STRESS )
	{
	    if ( rValues[GaussPoint].size() != 6 )
		rValues[GaussPoint].resize( 6 );
	    noalias( rValues[GaussPoint] )= mConstitutiveLawVector[GaussPoint]->GetValue( INSITU_STRESS, rValues[GaussPoint] );  
	}
    }

    KRATOS_CATCH( "" )
}
void Solid::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rProcessInfo ) {}
void Solid::SetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rProcessInfo ) {}
void Solid::SetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rProcessInfo ) {}
void Solid::CalculateOnIntegrationPoints( const Variable<double>& rVariable, Vector& Output, const ProcessInfo& rProcessInfo ) {}
int Solid::Check( const Kratos::ProcessInfo& rProcessInfo )
{
    return 0;
}


void Solid::Interpolate( const Variable<double>& rVariable, const ProcessInfo& rProcessInfo )
{
    for ( unsigned int k = 0; k < 4; k++ )
    {
        GetGeometry()[8+k].GetSolutionStepValue( rVariable ) = ( GetGeometry()[k%4].GetSolutionStepValue( rVariable ) + GetGeometry()[( k+1 )%4].GetSolutionStepValue( rVariable ) ) / 2.0;
        GetGeometry()[12+k].GetSolutionStepValue( rVariable ) = ( GetGeometry()[4+k%4].GetSolutionStepValue( rVariable ) + GetGeometry()[4+( k+1 )%4].GetSolutionStepValue( rVariable ) ) / 2.0;
        GetGeometry()[16+k].GetSolutionStepValue( rVariable ) = ( GetGeometry()[k%4].GetSolutionStepValue( rVariable ) + GetGeometry()[4+k%4].GetSolutionStepValue( rVariable ) ) / 2.0;
    }
}

void Solid::Interpolate( const Variable<Kratos::array_1d<double, 3> >& rVariable, const ProcessInfo& rProcessInfo )
{
    for ( unsigned int k = 0; k < 4; k++ )
    {
        GetGeometry()[8+k].GetSolutionStepValue( rVariable ) = ( GetGeometry()[k%4].GetSolutionStepValue( rVariable ) + GetGeometry()[( k+1 )%4].GetSolutionStepValue( rVariable ) ) / 2.0;
        GetGeometry()[12+k].GetSolutionStepValue( rVariable ) = ( GetGeometry()[4+k%4].GetSolutionStepValue( rVariable ) + GetGeometry()[4+( k+1 )%4].GetSolutionStepValue( rVariable ) ) / 2.0;
        GetGeometry()[16+k].GetSolutionStepValue( rVariable ) = ( GetGeometry()[k%4].GetSolutionStepValue( rVariable ) + GetGeometry()[4+k%4].GetSolutionStepValue( rVariable ) ) / 2.0;
    }
}

/// - WriteNodalResults
void Solid::FinalizeSolutionStep( ProcessInfo& ProcessInfo )
{  
    int nodesNumberOther = 8, nodesOtherMax= 8, nodesOtherMin=1;
    Matrix DNo_DX( nodesNumberOther, mDimension );
    Matrix DNo_Dxi( nodesNumberOther, mDimension );
    Matrix nodesLocalCoords( nodesNumberOther, mDimension );

    if ( GetGeometry().size() == 8 || GetGeometry().size() == 20 )
    {
        for ( unsigned int i = ( nodesOtherMin - 1 ) ; i < nodesOtherMax ;i++ )
            for ( unsigned int m = 0; m < mDimension; m++ )
            {
                nodesLocalCoords( i, m ) = 1.0;
                if ( m >= i )
                    nodesLocalCoords( i, m ) = -1.0;
            }
        nodesLocalCoords( 3, 0 ) = -1.0;
        nodesLocalCoords( 3, 2 ) = -1.0;
        nodesLocalCoords( 4, 0 ) = -1.0;
        nodesLocalCoords( 4, 1 ) = -1.0;
        nodesLocalCoords( 5, 1 ) = -1.0;
        nodesLocalCoords( 7, 0 ) = -1.0;
    }
    else if ( GetGeometry().size() == 4 )
    {
        nodesLocalCoords(0,0)=0.0;
        nodesLocalCoords(0,1)=0.0;
        nodesLocalCoords(0,2)=0.0;
        nodesLocalCoords(1,0)=1.0;
        nodesLocalCoords(1,1)=0.0;
        nodesLocalCoords(1,2)=0.0;
        nodesLocalCoords(2,0)=0.0;
        nodesLocalCoords(2,1)=1.0;
        nodesLocalCoords(2,2)=0.0;
        nodesLocalCoords(3,0)=0.0;
        nodesLocalCoords(3,1)=0.0;
        nodesLocalCoords(3,2)=1.0;
    }
 
    Matrix Grad_u( mDimension, mDimension );
    Vector stressVector ( 6 );

    for ( unsigned int i = ( nodesOtherMin - 1 ) ; i < nodesOtherMax ;i++ )
    {
        Vector local_coords( 3 );
        for ( unsigned int m = 0; m < mDimension; m++ )
            local_coords( m ) = nodesLocalCoords( i, m );

        if ( GetGeometry().size() == 8 || GetGeometry().size() == 20 )
            noalias( DNo_Dxi ) = GetGeometry().ShapeFunctionsLocalGradients( DNo_Dxi, local_coords );
        else if ( GetGeometry().size() == 4 )
        {
// 		  noalias( DNo_Dxi ) = GetGeometry().CalculateShapeFunctionsLocalGradients( DNo_Dxi, local_coords );
            DNo_Dxi(0,0) = -1.0;
            DNo_Dxi(0,1) = -1.0;
            DNo_Dxi(0,2) = -1.0;
            DNo_Dxi(1,0) =  1.0;
            DNo_Dxi(1,1) =  0.0;
            DNo_Dxi(1,2) =  0.0;
            DNo_Dxi(2,0) =  0.0;
            DNo_Dxi(2,1) =  1.0;
            DNo_Dxi(2,2) =  0.0;
            DNo_Dxi(3,0) =  0.0;
            DNo_Dxi(3,1) =  0.0;
            DNo_Dxi(3,2) =  1.0;
        }
        noalias( DNo_DX ) = prod( DNo_Dxi, mInvJ0[0] );
        noalias( Grad_u ) = GetGradu( DNo_DX );
        noalias( stressVector ) = GetstressVector ( Grad_u );
  
        GetGeometry()[i].GetSolutionStepValue( LINEAR_STRAIN ) = GetstrainVol( Grad_u ); 
//         GetGeometry()[i].GetSolutionStepValue( STRESS_VECTOR ) = stressVector; 
        GetGeometry()[i].GetSolutionStepValue( EQUIVALENT_VOLUMETRIC_STRESS ) =  Getp( stressVector );
        GetGeometry()[i].GetSolutionStepValue( EQUIVALENT_DEVIATORIC_STRESS ) =  Getq( stressVector );  
    }

    if ( GetGeometry().size() == 20 )
    { 
        Interpolate( LINEAR_STRAIN, ProcessInfo ); 
        Interpolate( EQUIVALENT_VOLUMETRIC_STRESS, ProcessInfo );
        Interpolate( EQUIVALENT_DEVIATORIC_STRESS, ProcessInfo );
    }

}
} // Namespace Kratos
