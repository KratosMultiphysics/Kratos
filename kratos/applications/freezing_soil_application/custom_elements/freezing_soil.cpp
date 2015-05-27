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

// 20130828
//  Assume: vL[Grad_p, t], comment DvL_Dp and DphiM_Dp
// 

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/freezing_soil.h"
#include "utilities/math_utils.h"
//#include "includes/constitutive_law.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/tetrahedra_3d_4.h"
#include "freezing_soil_application.h"
#include <math.h>
#include "freezing_soil.h"
#include "../applications/structural_application/custom_utilities/sd_math_utils.h"

namespace Kratos
{

FreezingSoil::FreezingSoil( IndexType NewId, GeometryType::Pointer pGeometry )
        : Element( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}


FreezingSoil::FreezingSoil( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
        : Element( NewId, pGeometry, pProperties )
{
    mNodesDispMin = 1; // Disp: displacement
    mNodesOtherMin = 1; // Other: ice volume fraction, water pressure, or temperature
    mNodesDispMax = GetGeometry().size();

    if ( mNodesDispMax == 8 || mNodesDispMax == 27 || mNodesDispMax == 20 )
    {
        mNodesOtherMax = 8;
        mThisIntegrationMethod = GeometryData::GI_GAUSS_3;//methods for hexahedra elements
        mThisGeometryOther = Geometry< Node<3> >::Pointer( new Hexahedra3D8 <Node<3> > (
                                 GetGeometry()( 0 ), GetGeometry()( 1 ), GetGeometry()( 2 ), GetGeometry()( 3 ),
                                 GetGeometry()( 4 ), GetGeometry()( 5 ), GetGeometry()( 6 ), GetGeometry()( 7 ) ) );          
// 	std::cout << "FreezingSoil 3D"<<mNodesDispMax<<"N"<< std::endl;

    } 
    //select Integration Method for 10 and 4 node tetrahedron
    else if (mNodesDispMax== 4 || mNodesDispMax== 10 )
    { 
        mNodesOtherMax = 4;
        mThisIntegrationMethod = GeometryData::GI_GAUSS_2;//methods for tetrahedra elements
        mThisGeometryOther = Geometry< Node<3> >::Pointer( new Tetrahedra3D4 <Node<3> > (
                                 GetGeometry()( 0 ), GetGeometry()( 1 ), GetGeometry()( 2 ), GetGeometry()( 3 ) ) );
// 	std::cout << "FreezingSoil 3D"<<mNodesDispMax<<"N"<< std::endl;
				 
    }
    else
        KRATOS_THROW_ERROR( std::logic_error, "This element matches only with a linear hexahedra/tetrahedra (8/3) or quadratic hexahedra/tetrahedra (20,27/10 ) geometry" , *this ); 
}

Element::Pointer FreezingSoil::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new FreezingSoil( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

FreezingSoil::~FreezingSoil()
{
}


void FreezingSoil::ResetConstitutiveLaw()
{
    KRATOS_TRY
    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        mConstitutiveLawVector[i]->ResetMaterial( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
    KRATOS_CATCH( "" )
}


void FreezingSoil::InitializeSolutionStep( ProcessInfo& CurrentProcessInfo )
{
}

////****************************************************************
//************************************************************************************
//************************************************************************************

// Degree of freedom: displacment: u[m], water pressure p[Pa]
void FreezingSoil::Initialize()
{
    KRATOS_TRY

    mNodesNumberDisp = mNodesDispMax - mNodesDispMin + 1;
    mNodesNumberOther = mNodesOtherMax - mNodesOtherMin + 1;
    mDimension = GetGeometry().WorkingSpaceDimension();
    mMatSizeU = mNodesNumberDisp * mDimension;
    mMatSizeO = mNodesNumberOther * 1;
    mNumberU = 1;
    mNumberO = 2;
    mMatSize = mMatSizeU * mNumberU + mMatSizeO * mNumberO;
    mAddIndexU = mMatSizeU * mNumberU;
    mScaleU = GetProperties()[SCALE_U];
    mScaleP = GetProperties()[SCALE_O];
    mTol = 1.0e-8; 

    mMaterialParameters = GetProperties()[ELEMENT_PARAMETERS]; 
// mMaterialParameters : [18](2700,1000,917,5e+10,2.2e+09,8.6e+09,3.75e+10,1e+08,3.4e+09,874,4190,2095,3,0.6,2.2,0,0,0) 
    mUnitRatio = GetProperties()[SCALE];				// Unit convertor: from [m] to [mm]
    mrhoS0 = mMaterialParameters[0] / pow( mUnitRatio, 3.0 ); 	// [ kg/m^3 ]
    mrhoL0 = mMaterialParameters[1] / pow( mUnitRatio, 3.0 );		// [ kg/m^3 ]
    mrhoC0 = mMaterialParameters[2] / pow( mUnitRatio, 3.0 ); 	// [ kg/m^3 ]
    mkS = mMaterialParameters[3] / mUnitRatio; 			// [ Pa = kg/(m*s^2) ]
    mkL = mMaterialParameters[4] / mUnitRatio; 			// [ Pa = kg/(m*s^2) ]
    mkC = mMaterialParameters[5] / mUnitRatio; 			// [ Pa = kg/(m*s^2) ]
    mgS = mMaterialParameters[6] / mUnitRatio; 			// [ Pa = kg/(m*s^2) ]
    mgL = mMaterialParameters[7] / mUnitRatio; 			// [ Pa = kg/(m*s^2) ]
    mgC = mMaterialParameters[8] / mUnitRatio; 			// [ Pa = kg/(m*s^2) ]
    mcS = mMaterialParameters[9] * pow( mUnitRatio, 2.0 ); 		// [ J/(kg*K) = m^2/(s^2*K) ]
    mcL = mMaterialParameters[10] * pow( mUnitRatio, 2.0 ); 		// [ J/(kg*K) = m^2/(s^2*K) ]
    mcC = mMaterialParameters[11] * pow( mUnitRatio, 2.0 ); 		// [ J/(kg*K) = m^2/(s^2*K) ]
    mlambdaS = mMaterialParameters[12] * mUnitRatio; 		// [ W/(m*K) = kg*m/(s^3*K) ]
    mlambdaL = mMaterialParameters[13] * mUnitRatio; 		// [ W/(m*K) = kg*m/(s^3*K) ]
    mlambdaC = mMaterialParameters[14] * mUnitRatio; 		// [ W/(m*K) = kg*m/(s^3*K) ]
    malphaS = mMaterialParameters[15]; 				// [ 1/K ]
    malphaL = mMaterialParameters[16]; 				// [ 1/K ]
    malphaC = mMaterialParameters[17]; 				// [ 1/K ] : depends on temperature:  (54-0.18(-t))x1e-6; see "Frozen Ground" page 53
    mTf = mMaterialParameters[18]; 				// [ K ]
    mSf = mMaterialParameters[19] / mUnitRatio; 			// [ Pa = kg/(m*s^2) ]
    mGravity = GetProperties()[GRAVITY] * mUnitRatio; 			// [m/s^2]
    mGravityDisp = mGravity; 						// [m/s^2]
    mn0 = GetProperties()[POROSITY]; 					// [ - ]
    mtstar = GetProperties()[FIRST_SATURATION_PARAM]; 			// [ K ]
    mm = GetProperties()[SECOND_SATURATION_PARAM]; 			// [ - ]
    mkappa0 = GetProperties()[PERMEABILITY_WATER] * mUnitRatio; 	// [ m/s ]

    if ( GetValue( KRATOS_WATCH_FLAG ) == 1 )
    {
        KRATOS_WATCH( mMaterialParameters );
        KRATOS_WATCH( mUnitRatio );
        KRATOS_WATCH( mn0 );
    }
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

    for ( unsigned int i = ( mNodesOtherMin - 1 ) ; i < mNodesOtherMax ; i++ )
    {
        GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE_NULL ) = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE );
        GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE_EINS ) = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE );
        GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE_DT ) = 0;
        GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE_NULL_DT ) = 0;
        GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE_EINS_DT ) = 0;
        GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE_ACCELERATION ) = 0;
        GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE_NULL_ACCELERATION ) = 0;
        GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE_EINS_ACCELERATION ) = 0;

        GetGeometry()[i].GetSolutionStepValue( TEMPERATURE_NULL ) = GetGeometry()[i].GetSolutionStepValue( TEMPERATURE );
        GetGeometry()[i].GetSolutionStepValue( TEMPERATURE_EINS ) = GetGeometry()[i].GetSolutionStepValue( TEMPERATURE );
        GetGeometry()[i].GetSolutionStepValue( TEMPERATURE_DT ) = 0;
        GetGeometry()[i].GetSolutionStepValue( TEMPERATURE_NULL_DT ) = 0;
        GetGeometry()[i].GetSolutionStepValue( TEMPERATURE_EINS_DT ) = 0;
        GetGeometry()[i].GetSolutionStepValue( TEMPERATURE_ACCELERATION ) = 0;
        GetGeometry()[i].GetSolutionStepValue( TEMPERATURE_NULL_ACCELERATION ) = 0;
        GetGeometry()[i].GetSolutionStepValue( TEMPERATURE_EINS_ACCELERATION ) = 0;
    }

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );
    GeometryType::JacobiansType J0( integration_points.size() );
    J0 = GetGeometry().Jacobian( J0, mThisIntegrationMethod );
    mInvJ0.resize( integration_points.size() );
    mDetJ0.resize( integration_points.size(), false );
    noalias( mDetJ0 ) = ZeroVector( integration_points.size() );

    //initial nodal temperature
    mT0e.resize( mNodesNumberOther );
    noalias( mT0e ) = ZeroVector( mNodesNumberOther );
    for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
        mT0e[i] = GetGeometry()[i].GetSolutionStepValue( TEMPERATURE_NULL );

    //initial temperature t0 in GaussPoints and InvJ0
    mT0.resize( integration_points.size() );
    noalias( mT0 ) = ZeroVector( integration_points.size() );
    const Matrix& No_container = mThisGeometryOther->ShapeFunctionsValues( mThisIntegrationMethod ); // for mT0
    for ( unsigned int GaussPoint = 0; GaussPoint < integration_points.size(); GaussPoint++ )
    {
        //initial temperature t0
        Vector No( mNodesNumberOther );
        noalias( No ) = row( No_container, GaussPoint );
        for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
            mT0[GaussPoint] += mT0e[i] * No( i );

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

    // initialized prestress
    mPrestressAssigned = false;
    mp0e.resize( mNodesNumberOther );
    noalias( mp0e ) = ZeroVector( mNodesNumberOther );
    mp0.resize( integration_points.size() );
    noalias( mp0 ) = ZeroVector( integration_points.size() );
    mStrainVol0.resize( integration_points.size() );
    noalias( mStrainVol0 ) = ZeroVector( integration_points.size() );

    KRATOS_CATCH( "" )
}

void FreezingSoil::CalculateAll( MatrixType& rLeftHandSideMatrix,
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
    const GeometryType::ShapeFunctionsGradientsType& DNo_Dxi = mThisGeometryOther->ShapeFunctionsLocalGradients( mThisIntegrationMethod );

    const Matrix& Nu_container = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );
    const Matrix& No_container = mThisGeometryOther->ShapeFunctionsValues( mThisIntegrationMethod );

    Vector Help_R_1( mMatSizeU );
    Vector Help_R_2( mMatSizeO );
    Vector Help_R_3( mMatSizeO );
    Matrix Help_K_UU( mMatSizeU, mMatSizeU );
    Matrix Help_K_UP( mMatSizeU, mMatSizeO );
    Matrix Help_K_UT( mMatSizeU, mMatSizeO );
    Matrix Help_K_PU( mMatSizeO, mMatSizeU );
    Matrix Help_K_PP( mMatSizeO, mMatSizeO );
    Matrix Help_K_PT( mMatSizeO, mMatSizeO );
    Matrix Help_K_TU( mMatSizeO, mMatSizeU );
    Matrix Help_K_TP( mMatSizeO, mMatSizeO );
    Matrix Help_K_TT( mMatSizeO, mMatSizeO );

    // ++++++ RightHandSideVector: R ++++++
    if ( CalculateResidualVectorFlag == true )
    {
        noalias( Help_R_1 ) = ZeroVector( mMatSizeU );
        noalias( Help_R_2 ) = ZeroVector( mMatSizeO );
        noalias( Help_R_3 ) = ZeroVector( mMatSizeO );
    }

    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        noalias( Help_K_UU ) = ZeroMatrix( mMatSizeU, mMatSizeU );
        noalias( Help_K_UP ) = ZeroMatrix( mMatSizeU, mMatSizeO );
        noalias( Help_K_UT ) = ZeroMatrix( mMatSizeU, mMatSizeO );
        noalias( Help_K_PU ) = ZeroMatrix( mMatSizeO, mMatSizeU );
        noalias( Help_K_PP ) = ZeroMatrix( mMatSizeO, mMatSizeO );
        noalias( Help_K_PT ) = ZeroMatrix( mMatSizeO, mMatSizeO );
        noalias( Help_K_TU ) = ZeroMatrix( mMatSizeO, mMatSizeU );
        noalias( Help_K_TP ) = ZeroMatrix( mMatSizeO, mMatSizeO );
        noalias( Help_K_TT ) = ZeroMatrix( mMatSizeO, mMatSizeO );
    }

    double weight, DetJ = 0.0 ;
    Vector Nu( mNodesNumberDisp ), No( mNodesNumberOther );
    Matrix DNu_DX( mNodesNumberDisp, mDimension ), DNo_DX( mNodesNumberOther, mDimension );

    double Div_dot_u, p, dot_p, t, dot_t;
    Vector Grad_p( mDimension ), Grad_t ( mDimension );
    Matrix Grad_u( mDimension, mDimension );

    // adding for calling constitutive law
    Vector InputVector( 9 );
    Vector strainVector( 6 );

    Vector stressVectorEff( 6 );
    noalias( stressVectorEff ) = ZeroVector( 6 );
    Vector DstressVectorEff_Dt( 6 );
    noalias( DstressVectorEff_Dt ) = ZeroVector( 6 );
    Matrix CtanEff( 6, 6 );
    noalias( CtanEff ) = ZeroMatrix( 6, 6 );
    Matrix dCtanEff_dt( 6, 6 );
    noalias( dCtanEff_dt ) = ZeroMatrix( 6, 6 );
     

    // START GAUSS INTEGRATION:
    for ( unsigned int GaussPoint = 0; GaussPoint < integration_points.size(); GaussPoint++ )
    {
        weight = integration_points[GaussPoint].Weight();
        DetJ = mDetJ0[GaussPoint];

        noalias( DNu_DX ) = prod( DNu_Dxi[GaussPoint], mInvJ0[GaussPoint] );
        noalias( DNo_DX ) = prod( DNo_Dxi[GaussPoint], mInvJ0[GaussPoint] );
        noalias( Nu ) = row( Nu_container, GaussPoint );
        noalias( No ) = row( No_container, GaussPoint );

        noalias( Grad_u ) = GetGradu( DNu_DX );
        Div_dot_u = GetDivdotu( DNu_DX );

        p = Getp( No );
        dot_p = Getdotp( No );
        noalias( Grad_p ) = GetGradp( DNo_DX );

        t = Gett( No );
        dot_t = Getdott( No );
        noalias( Grad_t ) = GetGradt( DNo_DX );

	
        // calling constitutive law for computing effective stress,
        noalias( strainVector ) = GetstrainVector( Grad_u );   // Compressive ! 

        if ( GetValue( PLASTIC_FLAG ) == 0 ) // ELATIC 
        {
            /// METHOD 1: CtanEff[ K[t], G[t] ]
            noalias( CtanEff ) = GetElasticTangent( t, 0 );
            noalias( dCtanEff_dt ) = GetElasticTangent( t, 1 );
            noalias( stressVectorEff ) = prod( CtanEff, strainVector );
            noalias( DstressVectorEff_Dt ) = prod( dCtanEff_dt, strainVector );
        }
        else // PLATIC 
        {
            /// METHOD 2: mBBM model
            for ( unsigned int i = 0; i < 6; i++ )
                InputVector[i] = strainVector[i];
            for ( unsigned int i = 0; i < 3; i++ )
                InputVector[i] += mStrainVol0[GaussPoint];
            InputVector[6] = p;
            InputVector[7] = t;
//             InputVector[8] = GetPorosity( Grad_u, p, t, mT0[GaussPoint], 0 );

            bool onlyGP1 = false;
	    if ( GetValue( KRATOS_WATCH_FLAG ) == 1 && GaussPoint == 1 )
	    {
                onlyGP1 = true;
		std::cout << "1---- u= "<<Getu(Nu)<<",\t Div_dot_u= " << Div_dot_u << std::endl;
		double pC = 0; 
		if ( t < 0 )    pC = p + mSf * ( -t );
		std::cout << "2---- pL= " << p << " (pC= " << pC << "),\t dot_pL= " << dot_p << ",\t Grad_pL= " << Grad_p << std::endl;
		std::cout << "3---- t= " << t << ",\t dot_t= " << dot_t << ",\t Grad_t= " << Grad_t << ",\t t0= " << mT0[GaussPoint] << std::endl;
	    }
	      
            // INPUT: strainVector (Tensile-Positive !!);  OUTPUT: stressVector [kPa] (Tensile-Positive !!) , CtanEff [kPa]
            mConstitutiveLawVector[GaussPoint]->CalculateMaterialResponse( InputVector, ZeroMatrix( 1, 1 ), stressVectorEff, CtanEff, rProcessInfo, GetProperties(), GetGeometry(), Nu, true, 1, onlyGP1 );
        } 
        
        if ( CalculateResidualVectorFlag == true )
        {
            //Calculation of spatial load vector
            AddInternalForcesToRHS1( Help_R_1, weight, DetJ, mp0[GaussPoint], mT0[GaussPoint], Grad_u, Grad_p, Grad_t, p, t, Div_dot_u, dot_p, dot_t, Nu, DNu_DX, stressVectorEff );
            AddInternalForcesToRHS2( Help_R_2, weight, DetJ, mp0[GaussPoint], mT0[GaussPoint], Grad_u, Grad_p, Grad_t, p, t, Div_dot_u, dot_p, dot_t, No, DNo_DX );
            AddInternalForcesToRHS3( Help_R_3, weight, DetJ, mp0[GaussPoint], mT0[GaussPoint], Grad_u, Grad_p, Grad_t, p, t, Div_dot_u, dot_p, dot_t, No, DNo_DX );
        }

        if ( CalculateStiffnessMatrixFlag == true )
        {
            //Calculation of spatial Stiffnes and Mass Matrix
            CalculateStiffnessMatrixUU( Help_K_UU, weight, DetJ, mp0[GaussPoint], mT0[GaussPoint], Grad_u, Grad_p, Grad_t, p, t, Div_dot_u, dot_p, dot_t, Nu, DNu_DX, CtanEff );
            CalculateStiffnessMatrixUP( Help_K_UP, weight, DetJ, mp0[GaussPoint], mT0[GaussPoint], Grad_u, Grad_p, Grad_t, p, t, Div_dot_u, dot_p, dot_t, Nu, DNu_DX, No );
            CalculateStiffnessMatrixUT( Help_K_UT, weight, DetJ, mp0[GaussPoint], mT0[GaussPoint], Grad_u, Grad_p, Grad_t, p, t, Div_dot_u, dot_p, dot_t, Nu, DNu_DX, No, DstressVectorEff_Dt );
            CalculateStiffnessMatrixPU( Help_K_PU, weight, DetJ, mp0[GaussPoint], mT0[GaussPoint], Grad_u, Grad_p, Grad_t, p, t, Div_dot_u, dot_p, dot_t, No, DNu_DX );
            CalculateStiffnessMatrixPP( Help_K_PP, weight, DetJ, mp0[GaussPoint], mT0[GaussPoint], Grad_u, Grad_p, Grad_t, p, t, Div_dot_u, dot_p, dot_t, No, DNo_DX );
            CalculateStiffnessMatrixPT( Help_K_PT, weight, DetJ, mp0[GaussPoint], mT0[GaussPoint], Grad_u, Grad_p, Grad_t, p, t, Div_dot_u, dot_p, dot_t, No, DNo_DX );
            CalculateStiffnessMatrixTU( Help_K_TU, weight, DetJ, mp0[GaussPoint], mT0[GaussPoint], Grad_u, Grad_p, Grad_t, p, t, Div_dot_u, dot_p, dot_t, No, DNu_DX );
            CalculateStiffnessMatrixTP( Help_K_TP, weight, DetJ, mp0[GaussPoint], mT0[GaussPoint], Grad_u, Grad_p, Grad_t, p, t, Div_dot_u, dot_p, dot_t, No, DNo_DX );
            CalculateStiffnessMatrixTT( Help_K_TT, weight, DetJ, mp0[GaussPoint], mT0[GaussPoint], Grad_u, Grad_p, Grad_t, p, t, Div_dot_u, dot_p, dot_t, No, DNo_DX );
        }

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

                for ( unsigned int j = ( mNodesOtherMin - 1 ); j < mNodesOtherMax; j++ )
                {
                    rLeftHandSideMatrix( i*mDimension + k, mAddIndexU + j*mNumberO ) += Help_K_UP( i * mDimension + k, j );
                    rLeftHandSideMatrix( i*mDimension + k, mAddIndexU + j*mNumberO + 1 ) += Help_K_UT( i * mDimension + k, j );
                }
            }

        for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
        {
            for ( unsigned int j = ( mNodesDispMin - 1 ); j < mNodesDispMax; j++ )
                for ( unsigned int l = 0; l < mDimension; l++ )
                {
                    rLeftHandSideMatrix( mAddIndexU + i*mNumberO, j*mDimension + l ) += Help_K_PU( i, j * mDimension + l );
                    rLeftHandSideMatrix( mAddIndexU + i*mNumberO + 1, j*mDimension + l ) += Help_K_TU( i, j * mDimension + l );
                }

            for ( unsigned int j = ( mNodesOtherMin - 1 ); j < mNodesOtherMax; j++ )
            {
                rLeftHandSideMatrix( mAddIndexU + i*mNumberO, mAddIndexU + j*mNumberO ) += Help_K_PP( i, j );
                rLeftHandSideMatrix( mAddIndexU + i*mNumberO, mAddIndexU + j*mNumberO + 1 ) += Help_K_PT( i, j );
                rLeftHandSideMatrix( mAddIndexU + i*mNumberO + 1, mAddIndexU + j*mNumberO ) += Help_K_TP( i, j );
                rLeftHandSideMatrix( mAddIndexU + i*mNumberO + 1, mAddIndexU + j*mNumberO + 1 ) += Help_K_TT( i, j );
            }
        }
    }

    if ( CalculateResidualVectorFlag == true )
    {
        for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
            for ( unsigned int k = 0; k < mDimension; k++ )
                rRightHandSideVector( i*mDimension + k ) -= Help_R_1( i * mDimension + k );

        for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
        {
            rRightHandSideVector( mAddIndexU + i*mNumberO ) -= Help_R_2( i );
            rRightHandSideVector( mAddIndexU + i*mNumberO + 1 ) -= Help_R_3( i );
        }
    }

if(GetValue(KRATOS_WATCH_FLAG)==1)
{
// KRATOS_WATCH(Help_R_2);
// KRATOS_WATCH(Help_K_PP);
}


    KRATOS_CATCH( "" )
}


void FreezingSoil::DampMatrix( MatrixType& rDampMatrix, ProcessInfo& rProcessInfo )
{
    KRATOS_TRY

    if ( rDampMatrix.size1() != mMatSize )
        rDampMatrix.resize( mMatSize, mMatSize );
    noalias( rDampMatrix ) = ZeroMatrix( mMatSize, mMatSize );

    //reading integration points and local Gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );
    const GeometryType::ShapeFunctionsGradientsType& DNu_Dxi = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );
    const GeometryType::ShapeFunctionsGradientsType& DNo_Dxi = mThisGeometryOther->ShapeFunctionsLocalGradients( mThisIntegrationMethod );
    const Matrix& Nu_container = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );
    const Matrix& No_container = mThisGeometryOther->ShapeFunctionsValues( mThisIntegrationMethod );

    Matrix Help_D_PU( mMatSizeO, mMatSizeU );
    Matrix Help_D_PP( mMatSizeO, mMatSizeO );
    Matrix Help_D_PT( mMatSizeO, mMatSizeO );
    Matrix Help_D_TU( mMatSizeO, mMatSizeU );
    Matrix Help_D_TP( mMatSizeO, mMatSizeO );
    Matrix Help_D_TT( mMatSizeO, mMatSizeO );
    noalias( Help_D_PU ) = ZeroMatrix( mMatSizeO, mMatSizeU );
    noalias( Help_D_PP ) = ZeroMatrix( mMatSizeO, mMatSizeO );
    noalias( Help_D_PT ) = ZeroMatrix( mMatSizeO, mMatSizeO );
    noalias( Help_D_TU ) = ZeroMatrix( mMatSizeO, mMatSizeU );
    noalias( Help_D_TP ) = ZeroMatrix( mMatSizeO, mMatSizeO );
    noalias( Help_D_TT ) = ZeroMatrix( mMatSizeO, mMatSizeO );

    double weight, DetJ = 0.0;
    Vector Nu( mNodesNumberDisp ), No( mNodesNumberOther );
    Matrix DNu_DX( mNodesNumberDisp, mDimension ), DNo_DX( mNodesNumberOther, mDimension );

    double Div_dot_u, p, dot_p, t, dot_t;
    Vector Grad_p( mDimension ), Grad_t ( mDimension );
    Matrix Grad_u( mDimension, mDimension );

    for ( unsigned int GaussPoint = 0; GaussPoint < integration_points.size(); GaussPoint++ )
    {
        weight = integration_points[GaussPoint].Weight();
        DetJ = mDetJ0[GaussPoint];

        noalias( DNu_DX ) = prod( DNu_Dxi[GaussPoint], mInvJ0[GaussPoint] );
        noalias( DNo_DX ) = prod( DNo_Dxi[GaussPoint], mInvJ0[GaussPoint] );
        noalias( Nu ) = row( Nu_container, GaussPoint );
        noalias( No ) = row( No_container, GaussPoint );

        noalias( Grad_u ) = GetGradu( DNu_DX );
        Div_dot_u = GetDivdotu( DNu_DX );

        p = Getp( No );
        dot_p = Getdotp( No );
        noalias( Grad_p ) = GetGradp( DNo_DX );

        t = Gett( No );
        dot_t = Getdott( No );
        noalias( Grad_t ) = GetGradt( DNo_DX );

	/// for fully frozen soil, no pL --- 20130821, C_AGF
	if( GetLiquidSaturation(t,0) < 0.01)
	  p = 0.0; 
	
        //Calculation of spatial damping matrix
        CalculateDampingMatrixPU( Help_D_PU, weight, DetJ, mp0[GaussPoint], mT0[GaussPoint], Grad_u, Grad_p, Grad_t, p, t, Div_dot_u, dot_p, dot_t, No, DNu_DX );
        CalculateDampingMatrixPP( Help_D_PP, weight, DetJ, mp0[GaussPoint], mT0[GaussPoint], Grad_u, Grad_p, Grad_t, p, t, Div_dot_u, dot_p, dot_t, No, DNo_DX );
        CalculateDampingMatrixPT( Help_D_PT, weight, DetJ, mp0[GaussPoint], mT0[GaussPoint], Grad_u, Grad_p, Grad_t, p, t, Div_dot_u, dot_p, dot_t, No, DNo_DX );
        CalculateDampingMatrixTU( Help_D_TU, weight, DetJ, mp0[GaussPoint], mT0[GaussPoint], Grad_u, Grad_p, Grad_t, p, t, Div_dot_u, dot_p, dot_t, No, DNu_DX );
        CalculateDampingMatrixTP( Help_D_TP, weight, DetJ, mp0[GaussPoint], mT0[GaussPoint], Grad_u, Grad_p, Grad_t, p, t, Div_dot_u, dot_p, dot_t, No, DNo_DX );
        CalculateDampingMatrixTT( Help_D_TT, weight, DetJ, mp0[GaussPoint], mT0[GaussPoint], Grad_u, Grad_p, Grad_t, p, t, Div_dot_u, dot_p, dot_t, No, DNo_DX );
    }

    for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
    {
        for ( unsigned int j = ( mNodesDispMin - 1 ); j < mNodesDispMax; j++ )
            for ( unsigned int l = 0; l < mDimension; l++ )
            {
                rDampMatrix( mAddIndexU + i*mNumberO, j*mDimension + l ) += Help_D_PU( i, j * mDimension + l );
                rDampMatrix( mAddIndexU + i*mNumberO + 1, j*mDimension + l ) += Help_D_TU( i, j * mDimension + l );
            }

        for ( unsigned int j = ( mNodesOtherMin - 1 ); j < mNodesOtherMax; j++ )
        {
            rDampMatrix( mAddIndexU + i*mNumberO, mAddIndexU + j*mNumberO ) += Help_D_PP( i, j );
            rDampMatrix( mAddIndexU + i*mNumberO, mAddIndexU + j*mNumberO + 1 ) += Help_D_PT( i, j );
            rDampMatrix( mAddIndexU + i*mNumberO + 1, mAddIndexU + j*mNumberO ) += Help_D_TP( i, j );
            rDampMatrix( mAddIndexU + i*mNumberO + 1, mAddIndexU + j*mNumberO + 1 ) += Help_D_TT( i, j );
        }
    }
    KRATOS_CATCH( "" )
}
///-//////////////////////////////////////////////////////////////////////////////////////
///-//////////////////////////////////////////////////////////////////////////////////////
/// R1 ***************************
void FreezingSoil::AddInternalForcesToRHS1( Vector& Help_R_1, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& Nu, const Matrix& DNu_DX, Vector& stressVectorEff )
{
    Matrix Bu( 6, mNodesNumberDisp*mDimension );
    noalias( Bu ) = GetBu( DNu_DX );

    double MLC = GetWaterAndIceMass( Grad_u, p, t, t0, 0 );
    double XL = GetLiquidSaturation( t, 0 );
    double K = GetBulkModulus( t, 0 );
    double b = GetBiotCoefficient( t, 0 );
  
    double pC = p;
    if ( t < 0.0 )
            pC = p + mSf * ( -t );
	
    for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
        for ( unsigned int m = 0; m < mDimension; m++ )
        {
            Help_R_1( i*mDimension + m ) += Nu( i ) * (
                                                (( 1.0 - mn0 ) * mrhoS0 + MLC ) * mGravityDisp( m )
                                            ) * weight * DetJ * mScaleU;

            Help_R_1( i*mDimension + m ) += DNu_DX( i, m ) * (
//                                                 b * ( p - p0 + ( 1.0 - XL ) * mSf * ( -t ) ) // COUSSY 
                                                 pC // NEW GENS 
                                            ) * weight * DetJ * mScaleU;

//             Help_R_1( i*mDimension + m ) += DNu_DX( i, m ) * (
//                                                 3.0 * malphaS * K * ( t - t0 ) 
// //                                                 3.0 * (malphaS * XL ) *  K * ( t - t0 ) // T=Tf: alphaS, T>Tf or T<<Tf ->0
//                                             ) * weight * DetJ * mScaleU;

            for ( unsigned int k = 0; k < 6; k++ )
                Help_R_1( i*mDimension + m ) -= Bu( k, i * mDimension + m ) * (
                                                    stressVectorEff( k )
                                                ) * weight * DetJ * mScaleU;
        }
}

void FreezingSoil::CalculateStiffnessMatrixUU( Matrix& Help_K_UU, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& Nu, const Matrix& DNu_DX, Matrix& CtanEff )
{
    Matrix Bu( 6, mNodesNumberDisp*mDimension );
    noalias( Bu ) = GetBu( DNu_DX );
    double DMLC_De = GetWaterAndIceMass( Grad_u, p, t, t0, 1 );

    for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
        for ( unsigned int m = 0; m < mDimension; m++ )
            for ( unsigned int j = ( mNodesDispMin - 1 ); j < mNodesDispMax; j++ )
                for ( unsigned int n = 0; n < mDimension; n++ )
                {
                    Help_K_UU( i*mDimension + m, j*mDimension + n ) += Nu( i ) * (
                                DMLC_De * mGravityDisp( m )
                            ) * DNu_DX( j, n ) * weight * DetJ * mScaleU;

                    for ( unsigned int k = 0; k < 6; k++ )
                        for ( unsigned int l = 0; l < 6; l++ )
                        {
                            Help_K_UU( i*mDimension + m, j*mDimension + n ) -= Bu( k, i * mDimension + m ) * (
                                        CtanEff( k, l )
                                    ) * Bu( l, j * mDimension + n ) * weight * DetJ * mScaleU;
                        }

                }
}

void FreezingSoil::CalculateStiffnessMatrixUP( Matrix& Help_K_UP, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& Nu, const Matrix& DNu_DX, const Vector& No )
{
    double DMLC_Dp = GetWaterAndIceMass( Grad_u, p, t, t0, 2 );
    double b = GetBiotCoefficient( t, 0 );
    double XL = GetLiquidSaturation( t, 0 );

    for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
        for ( unsigned int m = 0; m < mDimension; m++ )
            for ( unsigned int j = ( mNodesOtherMin - 1 ); j < mNodesOtherMax; j++ )
            {
                Help_K_UP( i*mDimension + m, j ) += Nu( i ) * (
                                                        DMLC_Dp * mGravityDisp( m )
                                                    ) * No( j ) * weight * DetJ * mScaleU;

                Help_K_UP( i*mDimension + m, j ) += DNu_DX( i , m ) * (
//                                                         b  // COUSSY
							1.0 // NEW GENS
                                                    ) * No( j ) * weight * DetJ * mScaleU;
            }
}

void FreezingSoil::CalculateStiffnessMatrixUT( Matrix& Help_K_UT, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& Nu, const Matrix& DNu_DX, const Vector& No, Vector& DstressVectorEff_Dt )
{
    Matrix Bu( 6, mNodesNumberDisp*mDimension );
    noalias( Bu ) = GetBu( DNu_DX );
    double DMLC_Dt = GetWaterAndIceMass( Grad_u, p, t, t0, 3 );
    double b = GetBiotCoefficient( t, 0 );
    double db_dt = GetBiotCoefficient( t, 1 );
    double K = GetBulkModulus( t, 0 );
    double dK_dt = GetBulkModulus( t, 1 );
    double XL = GetLiquidSaturation( t, 0 );
    double dXL_dt = GetLiquidSaturation( t, 1 );  
     

    double dpC_dt = 0.0;
    if ( t < 0.0 )
            dpC_dt = - mSf ;

    for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
        for ( unsigned int m = 0; m < mDimension; m++ )
            for ( unsigned int j = ( mNodesOtherMin - 1 ); j < mNodesOtherMax; j++ )
            {
                Help_K_UT( i*mDimension + m, j ) += Nu( i ) * (
                                                        DMLC_Dt * mGravityDisp( m )
                                                    ) * No( j ) * weight * DetJ * mScaleU;

                Help_K_UT( i*mDimension + m, j ) += DNu_DX( i , m ) * (
//                                                         db_dt * ( p - p0 + ( 1.0 - XL ) * mSf * ( -t ) )  - b * mSf * ( dXL_dt * ( -t ) + ( 1.0 - XL ) ) // COUSSY
                                                        dpC_dt  // NEW GENS
                                                    ) * No( j ) * weight * DetJ * mScaleU;

//                 Help_K_UT( i*mDimension + m, j ) += DNu_DX( i , m ) * (
//                                                         3.0 * malphaS * ( dK_dt * ( t - t0 ) + K )  
//                                                     ) * No( j ) * weight * DetJ * mScaleU;

                for ( unsigned int k = 0; k < 6; k++ )
                    Help_K_UT( i*mDimension + m, j ) -= Bu( k, i * mDimension + m ) * (
                                                            DstressVectorEff_Dt( k )
                                                        ) * No( j ) * weight * DetJ * mScaleU;
            }
}

///-//////////////////////////////////////////////////////////////////////////////////////
/// R2 **************************
void FreezingSoil::AddInternalForcesToRHS2( Vector& Help_R_2, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& No, const Matrix& DNo_DX )
{
    double rhoL = GetWaterDensity( p, t, 0 );
    double DrhoL_Dp = GetWaterDensity( p, t, 1 );
    double DrhoL_Dt = GetWaterDensity( p, t, 2 );

    double DMLC_De = GetWaterAndIceMass( Grad_u, p, t, t0, 1 );
    double DMLC_Dp = GetWaterAndIceMass( Grad_u, p, t, t0, 2 );
    double DMLC_Dt = GetWaterAndIceMass( Grad_u, p, t, t0, 3 );

    Vector vL( mDimension );
    noalias( vL ) = GetDarcyFlow( Grad_p, p, t, 0 );

    for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
    {
        Help_R_2( i ) += No( i ) * (
                             1.0 / rhoL * DMLC_De * Div_dot_u
                         ) * weight * DetJ * mScaleP;

        Help_R_2( i ) += No( i ) * (
                             1.0 / rhoL * DMLC_Dp * dot_p
                         ) * weight * DetJ * mScaleP;

        Help_R_2( i ) += No( i ) * (
                             1.0 / rhoL * DMLC_Dt * dot_t
                         ) * weight * DetJ * mScaleP;


        for ( unsigned int k = 0; k < mDimension; k++ )
        {
            Help_R_2( i ) += No( i ) * (
                                 1.0 / rhoL * ( DrhoL_Dp * Grad_p( k ) + DrhoL_Dt * Grad_t( k ) ) * vL( k )
                             ) * weight * DetJ * mScaleP;

            Help_R_2( i ) -= DNo_DX( i, k ) * (
                                 vL( k )
                             ) * weight * DetJ * mScaleP;
        }
    }
}

void FreezingSoil::CalculateStiffnessMatrixPU( Matrix& Help_K_PU, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& No, const Matrix& DNu_DX )
{
    double rhoL = GetWaterDensity( p, t, 0 );
    double D2MLC_De2 = GetWaterAndIceMass( Grad_u, p, t, t0, 11 );
    double D2MLC_DeDp = GetWaterAndIceMass( Grad_u, p, t, t0, 12 );
    double D2MLC_DeDt = GetWaterAndIceMass( Grad_u, p, t, t0, 13 );

    for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
        for ( unsigned int j = ( mNodesDispMin - 1 ); j < mNodesDispMax; j++ )
            for ( unsigned int n = 0; n < mDimension; n++ )
            {
                Help_K_PU( i, j*mDimension + n ) += No( i ) * (
                                                        1.0 / rhoL * D2MLC_De2 * Div_dot_u
                                                    ) * DNu_DX( j, n ) * weight * DetJ * mScaleU;

                Help_K_PU( i, j*mDimension + n ) += No( i ) * (
                                                        1.0 / rhoL * D2MLC_DeDp * dot_p
                                                    ) * DNu_DX( j, n ) * weight * DetJ * mScaleU;

                Help_K_PU( i, j*mDimension + n ) += No( i ) * (
                                                        1.0 / rhoL * D2MLC_DeDt * dot_t
                                                    ) * DNu_DX( j, n ) * weight * DetJ * mScaleU;
            }
}

void FreezingSoil::CalculateStiffnessMatrixPP( Matrix& Help_K_PP, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& No, const Matrix& DNo_DX )
{
    double rhoL = GetWaterDensity( p, t, 0 );
    double DrhoL_Dp = GetWaterDensity( p, t, 1 );
    double DrhoL_Dt = GetWaterDensity( p, t, 2 );

    double DMLC_De = GetWaterAndIceMass( Grad_u, p, t, t0, 1 );
    double DMLC_Dp = GetWaterAndIceMass( Grad_u, p, t, t0, 2 );
    double DMLC_Dt = GetWaterAndIceMass( Grad_u, p, t, t0, 3 );
    double D2MLC_DeDp = GetWaterAndIceMass( Grad_u, p, t, t0, 12 );
    double D2MLC_Dp2 = GetWaterAndIceMass( Grad_u, p, t, t0, 22 );
    double D2MLC_DpDt = GetWaterAndIceMass( Grad_u, p, t, t0, 23 );

    Vector vL( mDimension );
    noalias( vL ) = GetDarcyFlow( Grad_p, p, t, 0 );
    double DvL_DGradp = GetDarcyFlow( Grad_p, p, t, 1 )[0];
//     Vector DvL_Dp( mDimension );
//     noalias( DvL_Dp ) = GetDarcyFlow( Grad_p, p, t, 2 );

    for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
        for ( unsigned int j = ( mNodesOtherMin - 1 ); j < mNodesOtherMax; j++ )
        {
            Help_K_PP( i, j ) += No( i ) * (
                                     ( 1.0 / rhoL * D2MLC_DeDp - 1.0 / pow( rhoL, 2.0 ) * DrhoL_Dp * DMLC_De )
                                     * Div_dot_u
                                 ) * No( j ) * weight * DetJ * mScaleP;

            Help_K_PP( i, j ) += No( i ) * (
                                     ( 1.0 / rhoL * D2MLC_Dp2 - 1.0 / pow( rhoL, 2.0 ) * DrhoL_Dp * DMLC_Dp )
                                     * dot_p
                                 ) * No( j ) * weight * DetJ * mScaleP;

            Help_K_PP( i, j ) += No( i ) * (
                                     ( 1.0 / rhoL * D2MLC_DpDt - 1.0 / pow( rhoL, 2.0 ) * DrhoL_Dp * DMLC_Dt )
                                     * dot_t
                                 ) * No( j ) * weight * DetJ * mScaleP;

            for ( unsigned int k = 0; k < mDimension; k++ )
            {
                Help_K_PP( i, j ) += No( i ) * (
                                         ( DrhoL_Dp * Grad_p( k ) + DrhoL_Dt * Grad_t( k ) )
                                         * ( /*1.0 / rhoL * DvL_Dp( k ) */- 1.0 / pow( rhoL, 2.0 ) * DrhoL_Dp * vL( k ) )
                                     ) * No( j ) * weight * DetJ * mScaleP;

                Help_K_PP( i, j ) += No( i ) * (
                                         1.0 / rhoL * ( DrhoL_Dp * vL( k ) + ( DrhoL_Dp * Grad_p( k ) + DrhoL_Dt * Grad_t( k ) ) * DvL_DGradp )
                                     ) * DNo_DX( j , k ) * weight * DetJ * mScaleP;


//                 Help_K_PP( i, j ) -= DNo_DX( i , k ) * (
//                                          DvL_Dp( k )
//                                      ) *  No( j ) * weight * DetJ * mScaleP;

                Help_K_PP( i, j ) -= DNo_DX( i , k ) * (
                                         DvL_DGradp
                                     ) * DNo_DX( j , k ) * weight * DetJ * mScaleP;
            }
        }
}

void FreezingSoil::CalculateStiffnessMatrixPT( Matrix& Help_K_PT, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& No, const Matrix& DNo_DX )
{
    double rhoL = GetWaterDensity( p, t, 0 );
    double DrhoL_Dp = GetWaterDensity( p, t, 1 );
    double DrhoL_Dt = GetWaterDensity( p, t, 2 );

    double DMLC_De = GetWaterAndIceMass( Grad_u, p, t, t0, 1 );
    double DMLC_Dp = GetWaterAndIceMass( Grad_u, p, t, t0, 2 );
    double DMLC_Dt = GetWaterAndIceMass( Grad_u, p, t, t0, 3 );
    double D2MLC_DeDt = GetWaterAndIceMass( Grad_u, p, t, t0, 13 );
    double D2MLC_DpDt = GetWaterAndIceMass( Grad_u, p, t, t0, 23 );
    double D2MLC_Dt2 = GetWaterAndIceMass( Grad_u, p, t, t0, 33 );

    Vector vL( mDimension );
    noalias( vL ) = GetDarcyFlow( Grad_p, p, t, 0 );
    Vector DvL_Dt( mDimension );
    noalias( DvL_Dt ) = GetDarcyFlow( Grad_p, p, t, 3 );

    for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
        for ( unsigned int j = ( mNodesOtherMin - 1 ); j < mNodesOtherMax; j++ )
        {
            Help_K_PT( i, j ) += No( i ) * (
                                     ( 1.0 / rhoL * D2MLC_DeDt - 1.0 / pow( rhoL, 2.0 ) * DrhoL_Dt * DMLC_De )
                                     * Div_dot_u
                                 ) * No( j ) * weight * DetJ * mScaleP;

            Help_K_PT( i, j ) += No( i ) * (
                                     ( 1.0 / rhoL * D2MLC_DpDt - 1.0 / pow( rhoL, 2.0 ) * DrhoL_Dt * DMLC_Dp )
                                     * dot_p
                                 ) * No( j ) * weight * DetJ * mScaleP;

            Help_K_PT( i, j ) += No( i ) * (
                                     ( 1.0 / rhoL * D2MLC_Dt2 - 1.0 / pow( rhoL, 2.0 ) * DrhoL_Dt * DMLC_Dt )
                                     * dot_t
                                 ) * No( j ) * weight * DetJ * mScaleP;

            for ( unsigned int k = 0; k < mDimension; k++ )
            {

                Help_K_PT( i, j ) += No( i ) * (
                                         ( DrhoL_Dp * Grad_p( k ) + DrhoL_Dt * Grad_t( k ) )
                                         * ( 1.0 / rhoL * DvL_Dt( k ) - 1.0 / pow( rhoL, 2.0 ) * DrhoL_Dt * vL( k ) )
                                     ) * No( j ) * weight * DetJ * mScaleP;

                Help_K_PT( i, j ) += No( i ) * (
                                         1.0 / rhoL * DrhoL_Dt * vL( k )
                                     ) * DNo_DX( j , k ) * weight * DetJ * mScaleP;

                Help_K_PT( i, j ) -= DNo_DX( i , k ) * (
                                         DvL_Dt( k )
                                     ) * No( j ) * weight * DetJ * mScaleP;
            }
        }
}

// Damping
void FreezingSoil::CalculateDampingMatrixPU( Matrix& Help_D_PU, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& No, const Matrix& DNu_DX )
{
    double rhoL = GetWaterDensity( p, t, 0 );
    double DMLC_De = GetWaterAndIceMass( Grad_u, p, t, t0, 1 );

    for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
        for ( unsigned int j = ( mNodesDispMin - 1 ); j < mNodesDispMax; j++ )
            for ( unsigned int n = 0; n < mDimension; n++ )
                Help_D_PU( i, j*mDimension + n ) += No( i ) * (
                                                        1.0 / rhoL * DMLC_De
                                                    ) * DNu_DX( j, n ) * weight * DetJ * mScaleU;
}

void FreezingSoil::CalculateDampingMatrixPP( Matrix& Help_D_PP, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& No, const Matrix& DNo_DX )
{
    double rhoL = GetWaterDensity( p, t, 0 );
    double DMLC_Dp = GetWaterAndIceMass( Grad_u, p, t, t0, 2 );

    for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
        for ( unsigned int j = ( mNodesOtherMin - 1 ); j < mNodesOtherMax; j++ )
            Help_D_PP( i, j ) += No( i ) * (
                                     1.0 / rhoL * DMLC_Dp
                                 ) * No( j ) * weight * DetJ * mScaleP;
}

void FreezingSoil::CalculateDampingMatrixPT( Matrix& Help_D_PT, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& No, const Matrix& DNo_DX )
{
    double rhoL = GetWaterDensity( p, t, 0 );
    double DMLC_Dt = GetWaterAndIceMass( Grad_u, p, t, t0, 3 );

    for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
        for ( unsigned int j = ( mNodesOtherMin - 1 ); j < mNodesOtherMax; j++ )
            Help_D_PT( i, j ) += No( i ) * (
                                     1.0 / rhoL * DMLC_Dt
                                 ) * No( j ) * weight * DetJ * mScaleP;
}

///-//////////////////////////////////////////////////////////////////////////////////////
/// R3 **************************
void FreezingSoil::AddInternalForcesToRHS3( Vector& Help_R_3, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& No, const Matrix& DNo_DX )
{
    double DSS_De = GetSkeletonEntropy( Grad_u, p, t, 1 );
    double DSS_Dp = GetSkeletonEntropy( Grad_u, p, t, 2 );
    double DSS_Dt = GetSkeletonEntropy( Grad_u, p, t, 3 );
    double DsL_Dp = GetWaterEntropy( t, 1 );
    double DsL_Dt = GetWaterEntropy( t, 2 );
    double rhoL = GetWaterDensity( p, t, 0 );
    double ML = GetWaterMass( Grad_u, p, t, t0, 0 );
    double PhiM = GetMechanicalDissipation( Grad_p, p, t, 0 )[0];

    double dsCL = 0.0;
    double DsC_Dp = 0.0;
    double DsC_Dt = 0.0;
    double MC = 0.0;
    double DMC_De = 0.0;
    double DMC_Dp = 0.0;
    double DMC_Dt = 0.0;
    if ( t < 0.0 ) //SAVE COMPUTING TIME
    {
        dsCL = GetEntropyChange( p, t, 0 );
        DsC_Dp = GetIceEntropy( t, 1 );
        DsC_Dt = GetIceEntropy( t, 2 );
        MC = GetIceMass( Grad_u, p, t, t0, 0 );
        DMC_De = GetIceMass( Grad_u, p, t, t0, 1 );
        DMC_Dp = GetIceMass( Grad_u, p, t, t0, 2 );
        DMC_Dt = GetIceMass( Grad_u, p, t, t0, 3 );
    }

    Vector vL( mDimension );
    noalias( vL ) = GetDarcyFlow( Grad_p, p, t,  0 );
    Vector q( mDimension );
    noalias( q ) = GetHeatFlow( Grad_t, t, 0 );

    for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
    {
        Help_R_3( i ) += No( i ) * (
                             ( t + mTf ) * ( DSS_De + DMC_De * dsCL ) * Div_dot_u
                         ) * weight * DetJ;

        Help_R_3( i ) += No( i ) * (
                             ( t + mTf ) * ( DSS_Dp + DMC_Dp * dsCL + ML * DsL_Dp + MC * DsC_Dp ) * dot_p
                         ) * weight * DetJ;

        Help_R_3( i ) += No( i ) * (
                             ( t + mTf ) * ( DSS_Dt + DMC_Dt * dsCL + ML * DsL_Dt + MC * DsC_Dt ) * dot_t
                         ) * weight * DetJ;

        //mechanical dissipaton 2011-04-12
        Help_R_3( i ) -= No( i ) * (
                             PhiM
                         ) * weight * DetJ;

        for ( unsigned int k = 0; k < mDimension; k++ )
        {
            //heat convection 2011-05-31
            Help_R_3( i ) += No( i ) * (
                                 ( t + mTf ) * DsL_Dp * Grad_p( k ) * rhoL * vL( k )
                             ) * weight * DetJ;

            //heat convection 2011-05-31
            Help_R_3( i ) += No( i ) * (
                                 ( t + mTf ) * DsL_Dt * Grad_t ( k ) * rhoL * vL( k )
                             ) * weight * DetJ;

            //heat conduction
            Help_R_3( i ) -= DNo_DX( i, k ) * (
                                 q( k )
                             ) * weight * DetJ;
        }
    }
}


void FreezingSoil::CalculateStiffnessMatrixTU( Matrix& Help_K_TU, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& No, const Matrix& DNu_DX )
{
    double D2SS_DeDt = GetSkeletonEntropy( Grad_u, p, t, 13 );
    double DML_De = GetWaterMass( Grad_u, p, t, t0, 1 );
    double DsL_Dp = GetWaterEntropy( t, 1 );
    double DsL_Dt = GetWaterEntropy( t, 2 );

    double dsCL = 0.0;
    double DsC_Dp = 0.0;
    double DsC_Dt = 0.0;
    double DMC_De = 0.0;
    double D2MC_De2 = 0.0;
    double D2MC_DeDp = 0.0;
    double D2MC_DeDt = 0.0;
    if ( t < 0.0 ) //SAVE COMPUTING TIME
    {
        dsCL = GetEntropyChange( p, t, 0 );
        DsC_Dp = GetIceEntropy( t, 1 );
        DsC_Dt = GetIceEntropy( t, 2 );
        DMC_De = GetIceMass( Grad_u, p, t, t0, 1 );
        D2MC_De2 = GetIceMass( Grad_u, p, t, t0, 11 );
        D2MC_DeDp = GetIceMass( Grad_u, p, t, t0, 12 );
        D2MC_DeDt = GetIceMass( Grad_u, p, t, t0, 13 );
    }

    for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
        for ( unsigned int j = ( mNodesDispMin - 1 ); j < mNodesDispMax; j++ )
            for ( unsigned int n = 0; n < mDimension; n++ )
            {
                Help_K_TU( i, j*mDimension + n ) += No( i ) * (
                                                        ( t + mTf ) * ( D2MC_De2 * dsCL ) * Div_dot_u
                                                    ) * DNu_DX( j, n ) * weight * DetJ;

                Help_K_TU( i, j*mDimension + n ) += No( i ) * (
                                                        ( t + mTf ) * ( D2MC_DeDp * dsCL +  DML_De * DsL_Dp + DMC_De * DsC_Dp ) * dot_p
                                                    ) * DNu_DX( j, n ) * weight * DetJ;

                Help_K_TU( i, j*mDimension + n ) += No( i ) * (
                                                        ( t + mTf ) * ( D2SS_DeDt + D2MC_DeDt * dsCL + DML_De * DsL_Dt + DMC_De * DsC_Dt ) * dot_t
                                                    ) * DNu_DX( j, n ) * weight * DetJ;
            }
}

void FreezingSoil::CalculateStiffnessMatrixTP( Matrix& Help_K_TP, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& No, const Matrix& DNo_DX )
{
    double D2SS_DpDt = GetSkeletonEntropy( Grad_u, p, t, 23 );
    double DsL_Dp = GetWaterEntropy( t, 1 );
    double DsL_Dt = GetWaterEntropy( t, 2 );
    double rhoL = GetWaterDensity( p, t, 0 );
    double DrhoL_Dp = GetWaterDensity( p, t, 1 );
    double DML_Dp = GetWaterMass( Grad_u, p, t, t0, 2 );

    double dsCL = 0.0;
    double DdsCL_Dp = 0.0;
    double DsC_Dp = 0.0;
    double DsC_Dt = 0.0;
    double DMC_De = 0.0;
    double DMC_Dp = 0.0;
    double DMC_Dt = 0.0;
    double D2MC_DeDp = 0.0;
    double D2MC_DpDt = 0.0;
    double D2MC_Dp2 = 0.0;
    if ( t < 0.0 ) //SAVE COMPUTING TIME
    {
        dsCL = GetEntropyChange( p, t, 0 );
        DdsCL_Dp = GetEntropyChange( p, t, 1 );
        DsC_Dp = GetIceEntropy( t, 1 );
        DsC_Dt = GetIceEntropy( t, 2 );
        DMC_De = GetIceMass( Grad_u, p, t, t0, 1 );
        DMC_Dp = GetIceMass( Grad_u, p, t, t0, 2 );
        DMC_Dt = GetIceMass( Grad_u, p, t, t0, 3 );
        D2MC_DeDp = GetIceMass( Grad_u, p, t, t0, 12 );
        D2MC_DpDt = GetIceMass( Grad_u, p, t, t0, 23 );
        D2MC_Dp2 = GetIceMass( Grad_u, p, t, t0, 22 );
    }

    Vector vL( mDimension );
    noalias( vL ) = GetDarcyFlow( Grad_p, p, t,  0 );
    double DvL_DGradp = GetDarcyFlow( Grad_p, p, t, 1 )[0];
//     Vector DvL_Dp( mDimension );
//     noalias( DvL_Dp ) = GetDarcyFlow( Grad_p, p, t, 2 );

    Vector DPhiM_DGradp( mDimension );
    noalias( DPhiM_DGradp ) = GetMechanicalDissipation( Grad_p, p, t, 1 );
//     double DPhiM_Dp = GetMechanicalDissipation( Grad_p, p, t, 2 )[0];

    for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
        for ( unsigned int j = ( mNodesOtherMin - 1 ); j < mNodesOtherMax; j++ )
        {
            Help_K_TP( i, j ) += No( i ) * (
                                     ( t + mTf ) * ( D2MC_DeDp * dsCL + DMC_De * DdsCL_Dp ) * Div_dot_u
                                 ) * No( j ) * weight * DetJ;

            Help_K_TP( i, j ) += No( i ) * (
                                     ( t + mTf ) * ( D2MC_Dp2 * dsCL + DMC_Dp * DdsCL_Dp + DML_Dp * DsL_Dp + DMC_Dp * DsC_Dp ) * dot_p
                                 ) * No( j ) * weight * DetJ;

            Help_K_TP( i, j ) += No( i ) * (
                                     ( t + mTf ) * ( D2SS_DpDt + D2MC_DpDt * dsCL + DMC_Dt * DdsCL_Dp + DML_Dp * DsL_Dt + DMC_Dp * DsC_Dt ) * dot_t
                                 ) * No( j ) * weight * DetJ;

//             Help_K_TP( i, j ) -= No( i ) * (
//                                      DPhiM_Dp
//                                  ) * No( j ) * weight * DetJ;

            for ( unsigned int k = 0; k < mDimension; k++ )
            {
                Help_K_TP( i, j ) -= No( i ) * (
                                         DPhiM_DGradp( k )
                                     ) * DNo_DX( j , k ) * weight * DetJ;

		// old
//                 Help_K_TP( i, j ) += No( i ) * (
//                                          ( t + mTf ) * DsL_Dp * mrhoL0 * ( vL( k ) + Grad_p( k ) * DvL_DGradp )
//                                      ) * DNo_DX( j, k ) * weight * DetJ;

//                 Help_K_TP( i, j ) += No( i ) * (
//                                          ( t + mTf ) * DsL_Dt * Grad_t ( k ) * mrhoL0 * DvL_DGradp
//                                      ) * DNo_DX( j , k ) * weight * DetJ;
                // new
                Help_K_TP( i, j ) += No( i ) * (
                                         ( t + mTf ) * ( DsL_Dp * Grad_p( k ) + DsL_Dt * Grad_t( k ) )
                                         * ( DrhoL_Dp * vL( k ) /*+ rhoL * DvL_Dp( k ) */ )
                                     ) *  No( j ) * weight * DetJ;

                // new
                Help_K_TP( i, j ) += No( i ) * (
                                         ( t + mTf ) * DsL_Dp * rhoL * vL( k )
                                     ) * DNo_DX( j, k ) * weight * DetJ;

                // new
                Help_K_TP( i, j ) += No( i ) * (
                                         ( t + mTf )  * ( DsL_Dp * Grad_p( k ) + DsL_Dt * Grad_t( k ) )
                                         * rhoL * DvL_DGradp
                                     ) * DNo_DX( j , k ) * weight * DetJ;
            }
        }
}

void FreezingSoil::CalculateStiffnessMatrixTT( Matrix& Help_K_TT, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& No, const Matrix& DNo_DX )
{
    double DSS_De = GetSkeletonEntropy( Grad_u, p, t, 1 );
    double DSS_Dp = GetSkeletonEntropy( Grad_u, p, t, 2 );
    double DSS_Dt = GetSkeletonEntropy( Grad_u, p, t, 3 );
    double D2SS_DeDt = GetSkeletonEntropy( Grad_u, p, t, 13 );
    double D2SS_DpDt = GetSkeletonEntropy( Grad_u, p, t, 23 );
    double D2SS_Dt2 = GetSkeletonEntropy( Grad_u, p, t, 33 );
    double DsL_Dp = GetWaterEntropy( t, 1 );
    double DsL_Dt = GetWaterEntropy( t, 2 );
    double rhoL = GetWaterDensity( p, t, 0 );
    double DrhoL_Dt = GetWaterDensity( p, t, 2 );
    double ML = GetWaterMass( Grad_u, p, t, t0, 0 );
    double DML_Dt = GetWaterMass( Grad_u, p, t, t0, 3 );

    double dsCL = 0.0;
    double DdsCL_Dt =  0.0;
    double DsC_Dp =  0.0;
    double DsC_Dt =  0.0;
    double MC = 0.0;
    double DMC_De = 0.0;
    double DMC_Dp = 0.0;
    double DMC_Dt = 0.0;
    double D2MC_DeDt = 0.0;
    double D2MC_DpDt = 0.0;
    double D2MC_Dt2 = 0.0;
    if ( t < 0.0 ) //SAVE COMPUTING TIME
    {
        dsCL = GetEntropyChange( p, t, 0 );
        DdsCL_Dt = GetEntropyChange( p, t, 2 );
        DsC_Dp = GetIceEntropy( t, 1 );
        DsC_Dt = GetIceEntropy( t, 2 );
        MC = GetIceMass( Grad_u, p, t, t0, 0 );
        DMC_De = GetIceMass( Grad_u, p, t, t0, 1 );
        DMC_Dp = GetIceMass( Grad_u, p, t, t0, 2 );
        DMC_Dt = GetIceMass( Grad_u, p, t, t0, 3 );
        D2MC_DeDt = GetIceMass( Grad_u, p, t, t0, 13 );
        D2MC_DpDt = GetIceMass( Grad_u, p, t, t0, 23 );
        D2MC_Dt2 = GetIceMass( Grad_u, p, t, t0, 33 );
    }

    double Dq_DGradt = GetHeatFlow( Grad_t, t, 1 )[0];
    double DPhiM_Dt = GetMechanicalDissipation( Grad_p, p, t, 2 )[0];

    Vector vL( mDimension );
    noalias( vL ) = GetDarcyFlow( Grad_p, p, t,  0 );
    Vector DvL_Dt( mDimension );
    noalias( DvL_Dt ) = GetDarcyFlow( Grad_p, p, t,  2 );
    Vector q( mDimension );
    noalias( q ) = GetHeatFlow( Grad_t, t, 0 );
    Vector Dq_Dt( mDimension );
    noalias( Dq_Dt ) = GetHeatFlow( Grad_t, t, 2 );

    for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
        for ( unsigned int j = ( mNodesOtherMin - 1 ); j < mNodesOtherMax; j++ )
        {
            Help_K_TT( i, j ) += No( i ) * (
                                     ( DSS_De + DMC_De * dsCL ) * Div_dot_u
                                 ) * No( j ) * weight * DetJ;

            Help_K_TT( i, j ) += No( i ) * (
                                     ( t + mTf ) * ( D2SS_DeDt + D2MC_DeDt * dsCL + DMC_De * DdsCL_Dt ) * Div_dot_u
                                 ) * No( j ) * weight * DetJ;

            Help_K_TT( i, j ) += No( i ) * (
                                     ( DSS_Dp + DMC_Dp * dsCL + ML * DsL_Dp + MC * DsC_Dp ) * dot_p
                                 ) * No( j ) * weight * DetJ;

            Help_K_TT( i, j ) += No( i ) * (
                                     ( t + mTf ) * ( D2SS_DpDt + D2MC_DpDt * dsCL + DMC_Dp * DdsCL_Dt + DML_Dt * DsL_Dp + DMC_Dt * DsC_Dp ) * dot_p
                                 ) * No( j ) * weight * DetJ;

            Help_K_TT( i, j ) += No( i ) * (
                                     ( DSS_Dt + DMC_Dt * dsCL + ML * DsL_Dt + MC * DsC_Dt ) * dot_t
                                 ) * No( j ) * weight * DetJ;

            Help_K_TT( i, j ) += No( i ) * (
                                     ( t + mTf ) * ( D2SS_Dt2 + D2MC_Dt2 * dsCL + DMC_Dt * DdsCL_Dt + DML_Dt * DsL_Dt  + DMC_Dt * DsC_Dt ) * dot_t
                                 ) * No( j ) * weight * DetJ;

            Help_K_TT( i, j ) -= No( i ) * (
                                     DPhiM_Dt
                                 ) * No( j ) * weight * DetJ;

            for ( unsigned int k = 0; k < mDimension; k++ )
            {
	      // old
// //                 Help_K_TT( i, j )  += No( i ) * (
// //                                           DsL_Dp * Grad_p( k ) * mrhoL0 * vL( k )
// //                                       ) * No( j ) * weight * DetJ;
// //
// //                 Help_K_TT( i, j )  += No( i ) * (
// //                                           ( t + mTf ) * DsL_Dp * Grad_p( k ) * mrhoL0 * DvL_Dt( k )
// //                                       ) * No( j ) * weight * DetJ;
// //
// //                 Help_K_TT( i, j )  += No( i ) * (
// //                                           DsL_Dt * Grad_t ( k ) * mrhoL0 * vL( k )
// //                                       ) * No( j ) * weight * DetJ;
// //
// //                 Help_K_TT( i, j )  += No( i ) * (
// //                                           ( t + mTf ) * DsL_Dt * mrhoL0 * vL( k )
// //                                       ) * DNo_DX( j , k ) * weight * DetJ;
// //
// //                 Help_K_TT( i, j ) += No( i ) * (
// //                                          ( t + mTf ) * DsL_Dt * Grad_t ( k ) * mrhoL0 * DvL_Dt( k )
// //                                      ) * No( j ) * weight * DetJ;
                // new
                Help_K_TT( i, j )  += No( i ) * (
                                          ( DsL_Dp * Grad_p( k ) +  DsL_Dt * Grad_t ( k ) )
                                          * rhoL   * vL( k )
                                      ) * No( j ) * weight * DetJ;
                // new
                Help_K_TT( i, j )  += No( i ) * (
                                          ( t + mTf ) * ( DsL_Dp * Grad_p( k ) +  DsL_Dt * Grad_t ( k ) )
                                          * ( DrhoL_Dt  * vL( k )  + rhoL * DvL_Dt( k ) )
                                      ) * No( j ) * weight * DetJ;
                // new
                Help_K_TT( i, j )  += No( i ) * (
                                          ( t + mTf ) * DsL_Dt * rhoL * vL( k )
                                      ) * DNo_DX( j , k ) * weight * DetJ;

                Help_K_TT( i, j ) -= DNo_DX( i, k ) * (
                                         Dq_Dt( k )
                                     ) * No( j ) * weight * DetJ;

                Help_K_TT( i, j ) -= DNo_DX( i, k ) * (
                                         Dq_DGradt
                                     ) * DNo_DX( j , k ) * weight * DetJ;

            }
        }
}


void FreezingSoil::CalculateDampingMatrixTU( Matrix& Help_D_TU, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& No, const Matrix& DNu_DX )
{
    double DSS_De = GetSkeletonEntropy( Grad_u, p, t, 1 );

    double dsCL = 0.0;
    double DMC_De = 0.0;
    if ( t < 0.0 ) //SAVE COMPUTING TIME
    {
        dsCL = GetEntropyChange( p, t, 0 );
        DMC_De = GetIceMass( Grad_u, p, t, t0, 1 );
    }

    for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
        for ( unsigned int j = ( mNodesDispMin - 1 ); j < mNodesDispMax; j++ )
            for ( unsigned int n = 0; n < mDimension; n++ )
                Help_D_TU( i, j*mDimension + n ) += No( i ) * (
                                                        ( t + mTf ) * ( DSS_De + DMC_De * dsCL )
                                                    ) * DNu_DX( j, n ) * weight * DetJ;
}

void FreezingSoil::CalculateDampingMatrixTP( Matrix& Help_D_TP, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& No, const Matrix& DNo_DX )
{
    double DSS_Dp = GetSkeletonEntropy( Grad_u, p, t, 2 );
    double DsL_Dp = GetWaterEntropy( t, 1 );
    double ML = GetWaterMass( Grad_u, p, t, t0, 0 );

    double dsCL = 0.0;
    double DsC_Dp = 0.0;
    double MC = 0.0;
    double DMC_Dp = 0.0;
    if ( t < 0.0 ) //SAVE COMPUTING TIME
    {
        dsCL = GetEntropyChange( p, t, 0 );
        DsC_Dp = GetIceEntropy( t, 1 );
        MC = GetIceMass( Grad_u, p, t, t0, 0 );
        DMC_Dp = GetIceMass( Grad_u, p, t, t0, 2 );
    }
    for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
        for ( unsigned int j = ( mNodesOtherMin - 1 ); j < mNodesOtherMax; j++ )
            Help_D_TP( i, j ) += No( i ) * (
                                     ( t + mTf ) * ( DSS_Dp + DMC_Dp * dsCL + ML * DsL_Dp + MC * DsC_Dp )
                                 ) * No( j ) * weight * DetJ;

}

void FreezingSoil::CalculateDampingMatrixTT( Matrix& Help_D_TT, double weight, double DetJ, double p0, double t0, Matrix& Grad_u, Vector& Grad_p, Vector& Grad_t, double p, double t, double Div_dot_u, double dot_p, double dot_t, const Vector& No, const Matrix& DNo_DX )
{

    double DSS_Dt = GetSkeletonEntropy( Grad_u, p, t, 3 );
    double ML = GetWaterMass( Grad_u, p, t, t0, 0 );
    double DsL_Dt = GetWaterEntropy( t, 2 );

    double dsCL = 0.0;
    double DsC_Dt = 0.0;
    double MC = 0.0;
    double DMC_Dt = 0.0;
    if ( t < 0.0 ) //SAVE COMPUTING TIME
    {
        dsCL = GetEntropyChange( p, t, 0 );
        DsC_Dt = GetIceEntropy( t, 2 );
        MC = GetIceMass( Grad_u, p, t, t0, 0 );
        DMC_Dt = GetIceMass( Grad_u, p, t, t0, 3 );
    }

    for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
        for ( unsigned int j = ( mNodesOtherMin - 1 ); j < mNodesOtherMax; j++ )
            Help_D_TT( i, j ) += No( i ) * (
                                     ( t + mTf ) * ( DSS_Dt + DMC_Dt * dsCL + ML * DsL_Dt + MC * DsC_Dt )
                                 ) * No( j ) * weight * DetJ;
}

///-//////////////////////////////////////////////////////////////////////////////////////
///  CONSTITUTIVE QUANTITES
///-//////////////////////////////////////////////////////////////////////////////////////
/// KnoneckerDelta
double FreezingSoil::KnoneckerDelta( int i, int j )
{
    if ( i == j )
        return 1.0;
    else
        return 0.0;
}

/// Bu operator
Matrix FreezingSoil::GetBu( const Matrix& DNu_DX )
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
Vector FreezingSoil::Getu( const Vector& Nu )
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
Matrix FreezingSoil::GetGradu( const Matrix& DNu_DX )
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

double FreezingSoil::GetDivdotu( const Matrix& DNu_DX )
{
    double result = 0.0;

    for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
        for ( unsigned int k = 0; k < mDimension; k++ )
            result  += GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_DT )[k] * DNu_DX( i, k );

    return result;
}

///--------------------------------------------------------------------------
/// Water pressure
///--------------------------------------------------------------------------
double FreezingSoil::Getp( const Vector& No )
{
    double result = 0.0;

    for ( unsigned int i = mNodesOtherMin - 1 ; i < mNodesOtherMax ;i++ )
        result += GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE ) * No( i );

    return result;
}

double FreezingSoil::Getdotp( const Vector& No )
{
    double result = 0.0;

    for ( unsigned int i = mNodesOtherMin - 1 ; i < mNodesOtherMax ;i++ )
        result += GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE_DT ) * No( i );

    return result;
}

Vector FreezingSoil::GetGradp( const Matrix& DNo_DX )
{
    Vector result( mDimension );
    noalias( result ) = ZeroVector( mDimension );

    for ( unsigned int i = mNodesOtherMin - 1 ; i < mNodesOtherMax ;i++ )
        for ( unsigned int k = 0; k < mDimension; k++ )
            result[k] += GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE ) * DNo_DX( i, k );

    return result;
}

///--------------------------------------------------------------------------
/// Temperature
///--------------------------------------------------------------------------

double FreezingSoil::Gett( const Vector& No )
{
    double result = 0.0;

    for ( unsigned int i = mNodesOtherMin - 1 ; i < mNodesOtherMax ;i++ )
        result += GetGeometry()[i].GetSolutionStepValue( TEMPERATURE ) * No( i );

    return result;
}

double FreezingSoil::Getdott( const Vector& No )
{
    double result = 0.0;

    for ( unsigned int i = mNodesOtherMin - 1 ; i < mNodesOtherMax ;i++ )
        result += GetGeometry()[i].GetSolutionStepValue( TEMPERATURE_DT ) * No( i );

    return result;
}

Vector FreezingSoil::GetGradt( const Matrix& DNo_DX )
{
    Vector result( mDimension );
    noalias( result ) = ZeroVector( mDimension );

    for ( unsigned int i = mNodesOtherMin - 1 ; i < mNodesOtherMax ;i++ )
        for ( unsigned int k = 0; k < mDimension; k++ )
            result[k] += GetGeometry()[i].GetSolutionStepValue( TEMPERATURE ) * DNo_DX( i, k );

    return result;
}


/// ++++++++++++++ Group c1 +++++++++++++++++++
/// +++++++++++++++++++++++++++++++++++++++++++
/// elastic modulus: CtanEff
Matrix FreezingSoil::GetElasticTangent( double t, int index )
{
    Matrix result( 6, 6 );
    noalias( result ) = ZeroMatrix( 6, 6 );

    switch ( index )
    {
        case 0: // CtanEff
            {
                double K = GetBulkModulus( t, 0 );
                double G = GetShearModulus( t, 0 );
                double c1 = K + 4.0 * G / 3.0;
                double c2 = K - 2.0 * G / 3.0;
                for ( unsigned int i = 0; i < 3; i++ )
                {
                    result( i, i ) = c1;
                    result( i + 3, i + 3 ) = G;

                    for ( unsigned int j = 0; j < 3; j++ )
                        if ( i != j )
                            result( i, j ) = c2;
                }
            }
            break;
        case 1: // dCtanEff_dt
            {
                double dK_dt = GetBulkModulus( t, 1 );
                double dG_dt = GetShearModulus( t, 1 );
                double dc1_dt = dK_dt + 4.0 * dG_dt / 3.0;
                double dc2_dt = dK_dt - 2.0 * dG_dt / 3.0;
                for ( unsigned int i = 0; i < 3; i++ )
                {
                    result( i, i ) = dc1_dt;
                    result( i + 3, i + 3 ) = dG_dt;

                    for ( unsigned int j = 0; j < 3; j++ )
                        if ( i != j )
                            result( i, j ) = dc2_dt;
                }
            }
            break;
    }

    return result;
}
/// strain: linear strain tensor
Vector FreezingSoil::GetstrainVector( Matrix Grad_u )
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
double FreezingSoil::GetstrainVol( Matrix Grad_u )
{
    return  Grad_u( 0, 0 ) + Grad_u( 1, 1 ) + Grad_u( 2, 2 );
}

/// liquid saturation
double FreezingSoil::GetLiquidSaturation( double t, int index )
{
    double result = 0.0;
    switch ( index )
    {
        case 0: // XL
            {
                if ( t < 0.0 )
                    result =  pow( 1.0 + pow( -t / mtstar, 1.0 / ( 1.0 - mm ) ), -mm );
                else
                    result = 1.0;
            }
            break;
        case 1: // dXL_dt
            {
                if ( t < 0.0 )
                    result =  mm * pow( pow( -t / mtstar, 1.0 / ( 1.0 - mm ) ) + 1.0, -mm - 1.0 ) * pow( -t / mtstar, 1.0 / ( 1.0 - mm ) )
                              / (( mm - 1.0 ) * t );
            }
            break;
        case 2: // d2XL_dt2
            {
                if ( t < 0.0 )
                    result =  mm / pow(( mm - 1 ) * t, 2.0 ) * pow( pow( -t / mtstar, 1 / ( 1 - mm ) ) + 1, -mm - 2 ) * ( pow( -t / mtstar, 1 / ( 1 - mm ) ) - mm ) * pow( -t / mtstar, 1 / ( 1 - mm ) );
            }
            break;
    }

    return result;
}
/// K: effective drained bulk modulus
double FreezingSoil::GetBulkModulus( double t, int index )  // TODO add ice to the skeleton (drained bulk modulus of skeleton, elininate liquid phase!!!)
{
    double result = 0.0;

    /// Method 3 :Linear interpolation + Mori-Tanaka,
    /// M3-A: Matrix: Solid, Inclusion: Void/Ice:  Stiffer!!
//     double KSG= (4 *mgS *mkS *(1- mn0))/(4 *mgS + 3* mkS* mn0);
//     double KSC= (3*mkC*mkS + 4*mgS*(mkS + mkC*mn0 - mkS*mn0))/(4*mgS + 3* (mkC - mkC* mn0 + mkS* mn0));
    /// M3-B: Matrix: Liquid/Ice, Inclusion:Solid/Solid
    double kG = mgL;
    double gG = mgL;
    double KGS = ( 3 * kG * mkS + 4 * gG * ( mkS + kG * mn0 - mkS * mn0 ) ) / ( 4 * gG + 3 * ( kG - kG * mn0 + mkS * mn0 ) );
    double KCS = ( 3 * mkC * mkS + 4 * mgC * ( mkS + mkC * mn0 - mkS * mn0 ) ) / ( 4 * mgC + 3 * ( mkC - mkC * mn0 + mkS * mn0 ) );


    switch ( index )
    {
        case 0: // K
            {
                double XL = GetLiquidSaturation( t, 0 );
//                 return XL * KSG + ( 1.0 - XL ) * KSC; // M3-A
                return XL * KGS + ( 1.0 - XL ) * KCS; // M3-B

//                 result = 4.0 * mgS * mkS * ( 1 - mn0 ) / ( 3.0 * mn0 * mkS + 4.0 * mgS ); //M0: Coussy

//                 result= mkS*(1-mn0) + mkC*mn0*(1-XL);  // M3: linear interpolation and skeleton (Solid+Ice)

// 		double k0= mgL;
// 		double g0= k0;
//                 result = ( 3.0 * k0 * mkS + 4.0 * g0 * ( mkS + k0 * mn0 - mkS * mn0 ) ) / ( 4.0 * g0 + 3.0 * ( k0 - k0 * mn0 + mkS * mn0 ) ); //M4: Matrix:void, Inclusion;sand --> instable!!
            }
            break;
        case 1: // dK_dt
            {
                double dXL_dt = GetLiquidSaturation( t, 1 );
//                 return dXL_dt * ( KSG - KSC );  // M3-A
                return dXL_dt * ( KGS - KCS ); // M3-B
//                 result= mkC*mn0*(-dXL_dt);  // M3: linear interpolation and skeleton (Solid+Ice)
            }
            break;
        case 2: // d2K_dt2
            {
                double d2XL_dt2 = GetLiquidSaturation( t, 2 );
//                 return d2XL_dt2 * ( KSG - KSC ); //  M3-A
                return d2XL_dt2 * ( KGS - KCS ); // M3-B
//                 result= mkC*mn0*(-d2XL_dt2);  // M3: linear interpolation and skeleton (Solid+Ice)
            }
            break;
    }
    return result;
}

/// G: effective shear modulus
double FreezingSoil::GetShearModulus( double t, int index )  // TODO add ice to the skeleton (drained bulk modulus of skeleton, elininate liquid phase!!!)
{
    double result = 0.0;
//     double beta = 1.2 * ( mkC + 2.0 * mgC ) / ( 3.0 * mkC + 4.0 * mgC );
//     double GSL = mgL * ( 1.0 + ( 1.0 - mn0 ) * ( 1.0 - mgL / mgS ) / ( mgL / mgS + beta * mn0 * ( 1.0 - mgL / mgS ) ) );
//     double GSC = mgC * ( 1.0 + ( 1.0 - mn0 ) * ( 1.0 - mgC / mgS ) / ( mgC / mgS + beta * mn0 * ( 1.0 - mgC / mgS ) ) );

    // Method 3 :Linear interpolation + Mori-Tanaka,
    // M3-A: Matrix: Solid, Inclusion: Void/Ice:  Stiffer!!
//     double GSG = -((mgS* (8* mgS + 9* mkS)* (-1 + mn0))/( 8* mgS + 9 *mkS + 6 *(2 *mgS + mkS) *mn0));
//     double GSC =(mgS *(-mgS *(8* mgS + 9 *mkS)* (-1 + mn0) +  mgC* (12 *mgS + 6* mkS + 8 *mgS* mn0 + 9* mkS *mn0)))/(-6* mgC* (2 *mgS + mkS) *(-1 + mn0) + mgS* (8* mgS + 9* mkS + 6* (2 *mgS + mkS) * mn0)) ;

    // M3-B: Matrix: Liquid/Ice, Inclusion:Solid/Solid
    double kG = mgL;
    double gG = mgL;
    double GGS = -(( 5 * gG * ( 4 * gG + 3 * kG ) * mgS +   gG * ( 8 * gG + 9 * kG ) * ( gG - mgS ) * mn0 ) / ( -5 * gG * ( 4 * gG + 3 * kG ) +   6 * ( 2 * gG + kG ) * ( gG - mgS ) * mn0 ) );
    double GCS = -(( 5 * mgC * mgS * ( 4 * mgC + 3 * mkC ) +   mgC * ( mgC - mgS ) * ( 8 * mgC + 9 * mkC ) * mn0 ) / ( -5 * mgC * ( 4 * mgC + 3 * mkC ) +   6 * ( mgC - mgS ) * ( 2 * mgC + mkC ) * mn0 ) );

    switch ( index )
    {
        case 0: // G
            {
                double XL = GetLiquidSaturation( t, 0 );
// 		return  XL * GSG + ( 1.0 - XL ) * GSC; // M3-A
                return  XL * GGS + ( 1.0 - XL ) * GCS; // M3-B

//                 result =  mgS * ( 8.0 * mgS + 9.0 * mkS ) * ( 1 - mn0 ) / ( 6.0 * mn0 * ( 2.0 * mgS + mkS ) + 8.0 * mgS + 9.0 * mkS );  //M0: Coussy
//                 result= mgS*(1-mn0) + mgC*mn0*(1-XL);  //
            }
            break;
        case 1: // dG_dt
            {
                double dXL_dt = GetLiquidSaturation( t, 1 );
// 		return  dXL_dt * ( GSG - GSC );  // M3-A
                return  dXL_dt * ( GGS - GCS ); // M3-B
//                 result= mgC*mn0*(-dXL_dt);  //
            }
            break;
    }
    return result;
}
/// b: Biot coeffecient : M3: b= 1-KSCL/KSC
double FreezingSoil::GetBiotCoefficient( double t, int index )
{
    double result = 0.0;

//     double KCS = ( 3 * mkC * mkS + 4 * mgC * ( mkS + mkC * mn0 - mkS * mn0 ) ) / ( 4 * mgC + 3 * ( mkC - mkC * mn0 + mkS * mn0 ) );

//     switch ( index )
//     {
//         case 0: // b
//             {
//                 double K = GetBulkModulus( t, 0 );
//                 result =  1.0 - K / mkS;
//             }
//             break;
//         case 1: // db_dt
//             {
//                 double dK_dt = GetBulkModulus( t, 1 );
//                 result =  -dK_dt / mkS;
//             }
//             break;
//         case 2: // d2b_dt2
//             {
//                 double d2K_dt2 = GetBulkModulus( t, 2 );
//                 result = -d2K_dt2 / mkS;
//             }
//             break;
//     }
    switch ( index )
    {
        case 0: // b
            result =  1.0;
    }
    return result;
}

/// n: current porosity
double FreezingSoil::GetPorosity( Matrix Grad_u, double p, double t, double t0, int index )
{
    double result = 0.0;

//     double XL = GetLiquidSaturation( t, 0 );
//     double dXL_dt = GetLiquidSaturation( t, 1 );
//     double b = GetBiotCoefficient( t, 0 );
//     double db_dt = GetBiotCoefficient( t, 1 );
    double strainVol = GetstrainVol( Grad_u );

    // Method 1: Coussy 2011
//     switch ( index )
//     {
//         case 0: // n
//             result = mn0 + b * strainVol + ( b - mn0 ) / mkS * ( p + ( 1.0 - XL ) * mSf * ( -t ) ) - 3.0 * malphaS * ( b - mn0 ) * ( t - t0 );
//             break;
//         case 1: // Dn_De
//             result = b;
//             break;
//         case 2: // Dn_Dp
//             result = ( b - mn0 ) / mkS ;
//             break;
//         case 3: // Dn_Dt
//             result = db_dt * strainVol
//                      + db_dt / mkS * ( p + ( 1.0 - XL ) * mSf * ( -t ) )
//                      - ( b - mn0 ) / mkS * mSf * ( dXL_dt * ( -t ) + ( 1.0 - XL ) )
//                      - 3.0 * malphaS * ( db_dt * ( t - t0 )  + ( b - mn0 ) );
//             break;
//         case 13: // D2n_DeDt
//             result = db_dt;
//             break;
//         case 23: // D2n_DpDt
//             result = db_dt / mkS;
//             break;
//         case 33: // D2n_Dt2
//             {
//                 double d2XL_dt2 = GetLiquidSaturation( t, 2 );
//                 double d2b_dt2 = GetBiotCoefficient( t, 2 );
//                 result = d2b_dt2 * strainVol
//                          + d2b_dt2 / mkS * ( p + ( 1.0 - XL ) * mSf * ( -t ) )
//                          - 2.0 * db_dt / mkS * mSf * ( dXL_dt * ( -t ) + ( 1.0 - XL ) )
//                          - ( b - mn0 ) / mkS * mSf * ( d2XL_dt2 * ( -t ) - 2.0 * dXL_dt )
//                          - 3.0 * malphaS * ( d2b_dt2 * ( t - t0 )  + 2.0 * db_dt );
//             }
//     }
    // Method 2: Nagel 2009
    switch ( index )
    {
        case 0: // n
            result = 1.0 - ( 1.0 - mn0 ) * exp( - strainVol );
            break;
        case 1: // Dn_De
            result = ( 1.0 - mn0 ) * exp( - strainVol );
            break;
        case 11: // D2n_De2
            result = - ( 1.0 - mn0 ) * exp( - strainVol );
            break;
    }

    return result;
}


/// rhoL: current water density
double FreezingSoil::GetWaterDensity( double p, double t, int index )
{
    double result = 0.0;

//     double alphaL = 0.0;
// //     if ( t < 0.0 )
// 	      alphaL = malphaL;
    switch ( index )
    {
        case 0: // rhoL
            result =  mrhoL0  * ( 1.0 - 3.0 * malphaL * ( t )  + p / mkL );
            break;
        case 1: // DrhoL_Dp
            result =  mrhoL0 / mkL;
            break;
        case 2: // DrhoL_Dt
            result =  -3.0 * malphaL * mrhoL0;
            break;
    }  
    
    /// Boyle Marriott's law
//     double M_w = /*GetProperties()[MOLAR_MASS_WATER];*/ 24.042835;
// //     double R = /*GetProperties()[GAS_CONSTANT];*/  8.314472 * pow( mUnitRatio, 2.0 );
//     double R = /*GetProperties()[GAS_CONSTANT];*/  8.314472 * pow( mUnitRatio, 2.0 ); // 20131023--differenece bt water and gas
// 
//     switch ( index )
//     {
//         case 0: // rhoL
//             result = mrhoL0 + p * M_w / ( R * ( mTf + t ) ) ;
//             break;
//         case 1: // DrhoL_Dp
//             result =  M_w / ( R * ( mTf + t ) );
//             break;
//         case 2: // DrhoL_Dt
//             result = - M_w / R * pow( mTf + t, -2 ) * p;
//             break;
//     }

    return result;
}

/// rhoC: current ice density
double FreezingSoil::GetIceDensity( double p, double t, int index )
{
    double result = 0.0;
    switch ( index )
    {
        case 0: // rhoC
            result =  mrhoC0  * ( 1.0 - 3.0 * malphaC * t + ( p + mSf * ( -t ) ) / mkC ); //approximation
            break;
        case 1: // DrhoC_Dp
            result =  mrhoC0 / mkC;
            break;
        case 2: // DrhoC_Dt
            result =  -mrhoC0 * ( 3.0 * malphaC + mSf / mkC );
            break;
    }
    return result;
}

/// ML: water mass
double FreezingSoil::GetWaterMass( Matrix Grad_u, double p, double t, double t0, int index )
{
    double result = 0.0;
    double n = GetPorosity( Grad_u, p, t, t0, 0 );
    double rhoL = GetWaterDensity( p, t, 0 );
    double XL = GetLiquidSaturation( t, 0 );

    switch ( index )
    {
        case 0: // ML
            result = rhoL * n * XL;
            break;
        case 1: // DML_De
            {
                double Dn_De = GetPorosity( Grad_u, p, t, t0, 1 );
                result = rhoL * Dn_De * XL;
            }
            break;
        case 2: // DML_Dp
            {
                double Dn_Dp = GetPorosity( Grad_u, p, t, t0, 2 );
                double DrhoL_Dp = GetWaterDensity( p, t, 1 );
                result = ( DrhoL_Dp * n + rhoL * Dn_Dp ) * XL;
            }
            break;
        case 3: // DML_Dt
            {
                double Dn_Dt = GetPorosity( Grad_u, p, t, t0, 3 );
                double dXL_dt = GetLiquidSaturation( t, 1 );
                double DrhoL_Dt = GetWaterDensity( p, t, 2 );
                result = DrhoL_Dt * n * XL + rhoL * ( Dn_Dt * XL + n * dXL_dt );
            }
            break;
    }
    return result;
}

/// MC: ice mass
double FreezingSoil::GetIceMass( Matrix Grad_u, double p, double t, double t0, int index )
{
    double result = 0.0;
    if ( t < 0.0 )
    {
        double rhoC = GetIceDensity( p, t, 0 );
        double DrhoC_Dp = GetIceDensity( p, t, 1 );
        double DrhoC_Dt = GetIceDensity( p, t, 2 );
        double n = GetPorosity( Grad_u, p, t, t0, 0 );
        double Dn_De = GetPorosity( Grad_u, p, t, t0, 1 );
//         double Dn_Dp = GetPorosity( Grad_u, p, t, t0, 2 );
//         double Dn_Dt = GetPorosity( Grad_u, p, t, t0, 3 );
        double XL = GetLiquidSaturation( t, 0 );
        double dXL_dt = GetLiquidSaturation( t, 1 );
        switch ( index )
        {
            case 0: // MC
                result = rhoC * n * ( 1.0 - XL );
                break;
            case 1: // DMC_De
                result = rhoC * Dn_De * ( 1.0 - XL );
                break;
            case 2: // DMC_Dp
                result = ( DrhoC_Dp * n /*+ rhoC * Dn_Dp*/ ) * ( 1.0 - XL );
                break;
            case 3: // DMC_Dt
                result = DrhoC_Dt * n * ( 1.0 - XL ) + rhoC * ( /*Dn_Dt * ( 1.0 - XL )*/ - n * dXL_dt );
                break;
            case 11: // D2MC_De2
                {
                    double D2n_De2 = GetPorosity( Grad_u, p, t, t0, 11 );
                    result = rhoC * D2n_De2 * ( 1.0 - XL );
                }
                break;
            case 12: // D2MC_DeDp
                result = DrhoC_Dp * Dn_De * ( 1.0 - XL );
                break;
            case 13: // D2MC_DeDt
                {
//                     double D2n_DeDt = GetPorosity( Grad_u, p, t, t0, 13 );
                    result = DrhoC_Dt * Dn_De * ( 1.0 - XL ) + rhoC * (/* D2n_DeDt * ( 1.0 - XL ) */- Dn_De * dXL_dt );
                }
                break;
            case 23: // D2MC_DpDt
                {
//                     double D2n_DpDt = GetPorosity( Grad_u, p, t, t0, 23 );
                    result =/* (  DrhoC_Dp * Dn_Dt +DrhoC_Dt * Dn_Dp + rhoC * D2n_DpDt  ) * ( 1.0 - XL ) */- dXL_dt * ( DrhoC_Dp * n /*+ rhoC * Dn_Dp*/ );
                }
                break;
            case 22: // D2MC_Dp2
//                 result = ( 2.0 * DrhoC_Dp * Dn_Dp ) * ( 1.0 - XL );
                break;
            case 33: // D2MC_Dt2
                {
//                     double D2n_Dt2 = GetPorosity( Grad_u, p, t, t0, 33 );
                    double d2XL_dt2 = GetLiquidSaturation( t, 2 );
                    result = /*2.0 * DrhoC_Dt * Dn_Dt * ( 1.0 - XL )*/
                             - 2.0 * DrhoC_Dt * n * dXL_dt
//                              - 2.0 * rhoC * Dn_Dt * dXL_dt
                             + rhoC * ( /*D2n_Dt2 * ( 1.0 - XL ) */- n * d2XL_dt2 );
                }
                break;
        }
    }
    return result;
}

/// mass of ice and water
double FreezingSoil::GetWaterAndIceMass( Matrix Grad_u, double p, double t, double t0, int index )
{
    double result = 0.0;
    double n = GetPorosity( Grad_u, p, t, t0, 0 );
    double Dn_De = GetPorosity( Grad_u, p, t, t0, 1 );
    double Dn_Dp = GetPorosity( Grad_u, p, t, t0, 2 );
    double Dn_Dt = GetPorosity( Grad_u, p, t, t0, 3 );
    double D2n_De2 = GetPorosity( Grad_u, p, t, t0, 11 );

    double rhoL = GetWaterDensity( p, t, 0 );
    double DrhoL_Dp = GetWaterDensity( p, t, 1 );
    double DrhoL_Dt = GetWaterDensity( p, t, 2 );
    double rhoC = GetIceDensity( p, t, 0 );
    double DrhoC_Dp = GetIceDensity( p, t, 1 );
    double DrhoC_Dt = GetIceDensity( p, t, 2 );
    double XL = GetLiquidSaturation( t, 0 );
    double dXL_dt = GetLiquidSaturation( t, 1 );

    switch ( index )
    {
        case 0: // MLC
            result = n * ( rhoL * XL + rhoC * ( 1.0 - XL ) );
            break;
        case 1: // DMLC_De
            result = Dn_De * ( rhoL * XL + rhoC * ( 1.0 - XL ) );
            break;
        case 2: // DMLC_Dp
            result = Dn_Dp * ( rhoL * XL + rhoC * ( 1.0 - XL ) ) + n * ( DrhoL_Dp * XL + DrhoC_Dp * ( 1.0 - XL ) );
            break;
        case 3: // DMLC_Dt
            result = Dn_Dt * ( rhoL * XL + rhoC * ( 1.0 - XL ) ) + n * ( DrhoL_Dt * XL + DrhoC_Dt * ( 1.0 - XL ) + ( rhoL - rhoC ) * dXL_dt );
            break;
        case 11: // D2MLC_De2
            result = D2n_De2 * ( rhoL * XL + rhoC * ( 1.0 - XL ) );
            break;
        case 12: // D2MLC_DeDp
            result = Dn_De * ( DrhoL_Dp * XL + DrhoC_Dp * ( 1.0 - XL ) );
            break;
        case 13: // D2MLC_DeDt
            {
                double D2n_DeDt = GetPorosity( Grad_u, p, t, t0, 13 );
                result = D2n_DeDt * ( rhoL * XL + rhoC * ( 1.0 - XL ) ) + Dn_De * ( DrhoL_Dt * XL + DrhoC_Dt * ( 1.0 - XL ) + ( rhoL - rhoC ) * dXL_dt );
            }
            break;
        case 23: // D2MLC_DpDt
            {
                double D2n_DpDt = GetPorosity( Grad_u, p, t, t0, 23 );
                result = D2n_DpDt * ( rhoL * XL + rhoC * ( 1.0 - XL ) ) + Dn_Dp * ( DrhoL_Dt * XL + DrhoC_Dt * ( 1.0 - XL ) + ( rhoL - rhoC ) * dXL_dt ) + Dn_Dt * ( DrhoL_Dp * XL + DrhoC_Dp * ( 1.0 - XL ) ) + n * ( DrhoL_Dp - DrhoC_Dp ) * dXL_dt;
            }
            break;
        case 22: // D2MLC_Dp2
            result = 2.0 * Dn_Dp * ( DrhoL_Dp * XL + DrhoC_Dp * ( 1.0 - XL ) );
            break;
        case 33: // D2MLC_Dt2
            {
                double D2n_Dt2 = GetPorosity( Grad_u, p, t, t0, 33 );
                double d2XL_dt2 = GetLiquidSaturation( t, 2 );
                result = D2n_Dt2 * ( rhoL * XL + rhoC * ( 1.0 - XL ) ) + 2.0 * Dn_Dt * ( DrhoL_Dt * XL + DrhoC_Dt * ( 1.0 - XL ) + ( rhoL - rhoC ) * dXL_dt ) + n * ( 2.0 * ( DrhoL_Dt - DrhoC_Dt ) * dXL_dt + ( rhoL - rhoC ) * d2XL_dt2 );
            }
            break;
    }
    return result;
}

/// ++++++++++++++ Group c2 +++++++++++++++++++
/// +++++++++++++++++++++++++++++++++++++++++++
/// vL: liquid Darcy flow
Vector FreezingSoil::GetDarcyFlow( Vector Grad_p, double p, double t, int index )
{
    Vector result( mDimension );
    noalias( result ) = ZeroVector( mDimension );
    double lambda = 1.0;
    double dlambda_dXL = 0.0;
    double g = 9.81 * mUnitRatio;
    
    /// too steep for small mm --> convergence problem
    if ( t < 0.0 )
    {
        double XL = GetLiquidSaturation( t, 0 );
        lambda = sqrt( XL ) * pow( 1.0 - pow( 1.0 - pow( XL, 1.0 / mm ), mm ), 2.0 );     // Luckner et al.,1989
        dlambda_dXL = 0.5 * ( 1 - pow( 1 - pow( XL, 1 / mm ), mm ) ) * (( 1 - pow( 1 - pow( XL, 1 / mm ), mm ) ) / sqrt( XL ) + 4.0 * pow( XL, 1 / mm - 0.5 ) * pow( 1 - pow( XL, 1 / mm ), mm - 1 ) );
    }

    /// approximation --> good in convergence
//     double dlambda_dt = 0.0;
//     if ( t < 0.0 )
//     { 
// 	mtstar = 0.5*mtstar;
//         lambda =  GetLiquidSaturation( t, 0 ); 
//         dlambda_dt = GetLiquidSaturation( t, 1 );
// 	mtstar = 2.0*mtstar;
//     }
    
    switch ( index )
    {
        case 0: // vL
            for ( unsigned int k = 0; k < mDimension; k++ )
                result[k] =  mkappa0 * lambda / ( mrhoL0 * g ) * ( -Grad_p( k ) + mrhoL0 * mGravity( k ) );
            break;
        case 1: // DvL_DGradp : double
            result[0] = - mkappa0 * lambda / ( mrhoL0 * g );
            break;
        case 2: // DvL_Dp :
            {
            }
            break;
        case 3: // DvL_Dt
            {  
		double dXL_dt = GetLiquidSaturation( t, 1 );
                for ( unsigned int k = 0; k < mDimension; k++ )
//                     result[k] = ( mkappa0 / ( mrhoL0 * g ) * dlambda_dt ) * ( -Grad_p( k ) + mrhoL0 * mGravity( k ) ) ;
                    result[k] = ( mkappa0 / ( mrhoL0 * g ) * dlambda_dXL * dXL_dt ) * ( -Grad_p( k ) + mrhoL0 * mGravity( k ) ) ;
            }
            break;
    } 
    
    /// original  
//     double rhoL = GetWaterDensity( p, t, 0 );
//     if ( t < 0.0 )
//     {
//         double XL = GetLiquidSaturation( t, 0 );
//         lambda = sqrt( XL ) * pow( 1.0 - pow( 1.0 - pow( XL, 1.0 / mm ), mm ), 2.0 );     // Luckner et al.,1989
//         dlambda_dXL = 0.5 * ( 1 - pow( 1 - pow( XL, 1 / mm ), mm ) ) * (( 1 - pow( 1 - pow( XL, 1 / mm ), mm ) ) / sqrt( XL ) + 4.0 * pow( XL, 1 / mm - 0.5 ) * pow( 1 - pow( XL, 1 / mm ), mm - 1 ) );
//     }
// 
//     switch ( index )
//     {
//         case 0: // vL
//             for ( unsigned int k = 0; k < mDimension; k++ )
//                 result[k] =  mkappa0 * lambda / ( rhoL * g ) * ( -Grad_p( k ) + rhoL * mGravity( k ) );
//             break;
//         case 1: // DvL_DGradp : double
//             result[0] = - mkappa0 * lambda / ( rhoL * g );
//             break;
//         case 2: // DvL_Dp :
//             {
//                 double DrhoL_Dp = GetWaterDensity( p, t, 1 ); 
//                 for ( unsigned int k = 0; k < mDimension; k++ )
//                     result[k] =  mkappa0 * lambda / ( pow( rhoL, 2.0 ) * g ) * DrhoL_Dp * Grad_p( k );
//             }
//             break;
//         case 3: // DvL_Dt
//             {
//                 double DrhoL_Dt = GetWaterDensity( p, t, 2 ); 
//                 double dXL_dt = GetLiquidSaturation( t, 1 );
//                 for ( unsigned int k = 0; k < mDimension; k++ )
//                     result[k] = ( mkappa0 / ( rhoL * g ) * dlambda_dXL * dXL_dt - mkappa0 * lambda / ( pow( rhoL, 2.0 ) * g ) * DrhoL_Dt ) 
// 					* ( -Grad_p( k ) + rhoL * mGravity( k ) )
//                                 + mkappa0 * lambda / ( rhoL * g ) * DrhoL_Dt * mGravity( k );
//             }
//             break;
//     } 

    return result;
}

/// PhiM: mechanical dissipation
Vector FreezingSoil::GetMechanicalDissipation( Vector Grad_p, double p, double t, int index )
{
    Vector result( mDimension );
    noalias( result ) = ZeroVector( mDimension );

//     double rhoL = GetWaterDensity( p, t, 0 );
//     Vector vL( mDimension );
//     noalias( vL ) = GetDarcyFlow( Grad_p, p, t, 0 );
// 
//     switch ( index )
//     {
//         case 0: // PhiM : double
//             {
//                 for ( unsigned int k = 0; k < mDimension; k++ )
//                     result[0] += vL( k ) * ( -Grad_p( k ) + mrhoL0 * mGravity( k ) );
//             }
//             break;
//         case 1: // DPhiM_DGradp
//             {
//                 double DvL_DGradp = GetDarcyFlow( Grad_p, p, t, 1 )[0];
//                 for ( unsigned int k = 0; k < mDimension; k++ )
//                     result[k] = DvL_DGradp * ( -Grad_p( k ) + mrhoL0 * mGravity( k ) ) - vL( k )  ;
//             }
//             break;
// 	case 2: // DPhiM_Dp: double
//             {
// //                 Vector DvL_Dp( mDimension );
// //                 noalias( DvL_Dp ) = GetDarcyFlow( Grad_p, p, t, 2 );
// //                 double DrhoL_Dp = GetWaterDensity( p, t, 1 );
// //                 for ( unsigned int k = 0; k < mDimension; k++ )
// //                     result[0] += DvL_Dp( k ) * ( -Grad_p( k ) + mrhoL0 * mGravity( k ) ) /*+ vL( k ) * DrhoL_Dp * mGravity( k )*/;
//             }
//             break;
//         case 3: // DPhiM_Dt : double
//             {
//                 Vector DvL_Dt( mDimension );
//                 noalias( DvL_Dt ) = GetDarcyFlow( Grad_p, p, t, 3 );
//                 double DrhoL_Dt = GetWaterDensity( p, t, 2 );
//                 for ( unsigned int k = 0; k < mDimension; k++ )
//                     result[0] += DvL_Dt( k ) * ( -Grad_p( k ) + mrhoL0 * mGravity( k ) ) /*+ vL( k ) * DrhoL_Dt * mGravity( k )*/;
//             }
//             break;
//     }
    return result;
}


/// ++++++++++++++ Group c3 +++++++++++++++++++
/// +++++++++++++++++++++++++++++++++++++++++++
/// SS: entropy of skeleton
double FreezingSoil::GetSkeletonEntropy( Matrix Grad_u, double p, double t, int index )
{
    double result = 0.0;
    double K = GetBulkModulus( t, 0 );
    double dK_dt = GetBulkModulus( t, 1 );
    double b = GetBiotCoefficient( t, 0 );
    double db_dt = GetBiotCoefficient( t, 1 );
    double XL = GetLiquidSaturation( t, 0 );
    double dXL_dt = GetLiquidSaturation( t, 1 );
    double strainVol = GetstrainVol( Grad_u );

    // M3
//     double KM2= 4.0*mgS*mkS*( 1 - mn0 ) / ( 3.0*mn0*mkS + 4.0*mgS ); //M2: Coussy
//     double b = 1- KM2/mkS;
//     double dK_dt = 0.0;
//     double db_dt = 0;

    switch ( index )
    {
        case 1: // DSS_De
            result = 3.0 * malphaS * K; 
            break;
        case 2: // DSS_Dp
            result = - 3.0 * malphaS * ( b - mn0 ); 
            break;
        case 3: // DSS_Dt
            result = 3.0 * malphaS * dK_dt * strainVol
                     - 3.0 * malphaS * db_dt * ( p + ( 1.0 - XL ) * mSf * ( -t ) )
                     + 3.0 * malphaS * ( b - mn0 ) * mSf * ( dXL_dt * ( -t ) + ( 1.0 - XL ) )
                     + mrhoS0 * ( 1.0 - mn0 ) * mcS / mTf;
//                      + mrhoS0 * ( 1.0 - mn0 ) * mcS / (mTf + t); //NEW 20130827 
            break;
        case 13: // D2SSDeDt
            result =  3.0 * malphaS * dK_dt;
            break;
        case 23: // D2SSDpDt
            result = - 3.0 * malphaS * db_dt ;
            break;
        case 33: // D2SSDt2
            {
                double d2XL_dt2 = GetLiquidSaturation( t, 2 );
                double d2K_dt2 = GetBulkModulus( t, 2 );
                double d2b_dt2 = GetBiotCoefficient( t, 2 );
                result =  3.0 * malphaS * d2K_dt2 * strainVol
                          - 3.0 * malphaS * d2b_dt2 * ( p + ( 1.0 - XL ) * mSf * ( -t ) )
                          + 6.0 * malphaS * db_dt * mSf * ( dXL_dt * ( -t ) + ( 1.0 - XL ) )
                          + 3.0 * malphaS * ( b - mn0 ) * mSf * ( d2XL_dt2 * ( -t ) - 2.0 * dXL_dt );
            }
            break;
    }
    return result;
}
/// sL: specific water entropy
double FreezingSoil::GetWaterEntropy( double t, int index )
{
    double result = 0.0;

    double alphaL = 0.0;
//     if ( t < 0.0 )
        alphaL = malphaL;
    switch ( index )
    {
        case 1: // DsL_Dp
            result = - 3.0 * alphaL / mrhoL0;
            break;
        case 2: // DsL_Dt
            result = mcL / mTf;
//             result = mcL / (mTf + t); //NEW 20130827
            break;
    }
    return result;
}

/// sC: specific entropy
double FreezingSoil::GetIceEntropy( double t, int index )
{
    double result = 0.0;
    switch ( index )
    {
        case 1: // DsL_Dp
            result = -3.0 * malphaC / mrhoC0;
            break;
        case 2: // DsL_Dt
            result = 3.0 * malphaC * mSf / mrhoC0 + mcC / mTf;
//             result = 3.0 * malphaC * mSf / mrhoC0 + mcC / (mTf + t); //NEW 20130827
            break;
    }
    return result;
}

/// ds = sC-sL : entropy differenece
double FreezingSoil::GetEntropyChange( double p, double t, int index )
{
    double result = 0.0;

    double alphaL = 0.0;
//     if ( t < 0.0 )
        alphaL = malphaL;
    switch ( index )
    {
        case 0: // dsCL
            result = ( mcC - mcL ) * t / mTf + 3.0 * ( alphaL / mrhoL0 - malphaC / mrhoC0 ) * p + 3.0 * malphaC * mSf * t / mrhoC0;
//             result = ( mcC - mcL ) * log( (mTf + t) / mTf ) + 3.0 * ( alphaL / mrhoL0 - malphaC / mrhoC0 ) * p + 3.0 * malphaC * mSf * t / mrhoC0; //NEW 20130827
            break;
        case 1: // DdsCL_Dp
            result = 3.0 * ( alphaL / mrhoL0 - malphaC / mrhoC0 );
            break;
        case 2: // DdsCL_Dt
            result = ( mcC - mcL ) / mTf + 3.0 * malphaC * mSf / mrhoC0;
//             result = ( mcC - mcL ) / (mTf+t) + 3.0 * malphaC * mSf / mrhoC0; //NEW 20130827
    }
    return result;
}

/// +++++++++++++++++++++++++++++++++++++++++++
/// q: heat flux
Vector FreezingSoil::GetHeatFlow( Vector Grad_t, double t, int index )
{
    Vector result( mDimension );
    noalias( result ) = ZeroVector( mDimension );
    double XL = GetLiquidSaturation( t, 0 );
    double lambdaTot = pow( mlambdaS, 1.0 - mn0 ) * pow( mlambdaL, mn0 * XL ) * pow( mlambdaC, mn0 * ( 1.0 - XL ) );    //Zhu

    switch ( index )
    {
        case 0: // q
            for ( unsigned int k = 0; k < mDimension; k++ )
                result[k] = -lambdaTot * Grad_t ( k );
            break;
        case 1: // Dq_DGradt : double
            result[0] = -lambdaTot;
            break;
        case 2: // Dq_Dt
            {
                double dXL_dt = GetLiquidSaturation( t, 1 );
                double dlambdaTot_dXL = mn0 * lambdaTot * ( log( mlambdaL ) - log( mlambdaC ) );
                for ( unsigned int k = 0; k < mDimension; k++ )
                    result[k] = -dlambdaTot_dXL * dXL_dt * Grad_t ( k );
            }
    }
    return result;
}
/// SOLUTION

FreezingSoil::IntegrationMethod FreezingSoil::GetIntegrationMethod()
{
    return mThisIntegrationMethod;
}

void FreezingSoil::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rProcessInfo )
{
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll( temp, rRightHandSideVector, rProcessInfo, CalculateStiffnessMatrixFlag,  CalculateResidualVectorFlag );
}


void FreezingSoil::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rProcessInfo )
{
    //calculation flamgS
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;

    CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
}

void FreezingSoil::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& ProcessInfo )
{
    if ( rResult.size() != mMatSize )
        rResult.resize( mMatSize );

    for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
    {
        rResult[i*mDimension]   = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[i*mDimension+1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();
        rResult[i*mDimension+2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
    }

    for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
    {
        rResult[mAddIndexU+i*mNumberO]   = GetGeometry()[i].GetDof( WATER_PRESSURE ).EquationId();
        rResult[mAddIndexU+i*mNumberO+1] = GetGeometry()[i].GetDof( TEMPERATURE ).EquationId();
    }
}


void FreezingSoil::GetDofList( DofsVectorType& ElementalDofList, ProcessInfo& ProcessInfo )
{
    ElementalDofList.resize( 0 );

    for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
    {
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
    }

    for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
    {
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( WATER_PRESSURE ) );
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( TEMPERATURE ) );
    }
}


void FreezingSoil::GetValuesVector( Vector& values, int Step )
{
    if ( values.size() != mMatSize )
        values.resize( mMatSize );

    for ( unsigned int i = ( mNodesDispMin - 1 );i < mNodesDispMax;i++ )
    {
        values( i*mDimension )   = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
        values( i*mDimension + 1 ) = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );
        values( i*mDimension + 2 ) = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );
    }

    for ( unsigned int i = ( mNodesOtherMin - 1 );i < mNodesOtherMax;i++ )
    {
        values( mAddIndexU + i*mNumberO )   = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE, Step );
        values( mAddIndexU + i*mNumberO + 1 ) = GetGeometry()[i].GetSolutionStepValue( TEMPERATURE, Step );
    }
}

///-////////////////////////////////////
/////
///- //////////////////////////////////////////////////////
void FreezingSoil::GetValueOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rProcessInfo )
{
	KRATOS_TRY
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );
    if(rValues.size() != integration_points.size())
	    rValues.resize( integration_points.size() );
 
    const GeometryType::ShapeFunctionsGradientsType& DNu_Dxi = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );
    const GeometryType::ShapeFunctionsGradientsType& DNo_Dxi = mThisGeometryOther->ShapeFunctionsLocalGradients( mThisIntegrationMethod );
    const Matrix& Nu_container = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );
    const Matrix& No_container = mThisGeometryOther->ShapeFunctionsValues( mThisIntegrationMethod );

    Vector Nu( mNodesNumberDisp ), No( mNodesNumberOther );
    Matrix DNu_DX( mNodesNumberDisp, mDimension ), DNo_DX( mNodesNumberOther, mDimension );
    
    Vector u( mDimension ); // solid displacement
    double p; // relative water pressure
    double t; // temperature in Celsius degree

    Matrix Grad_u( mDimension, mDimension );
    Vector Grad_p( mDimension );
    Vector Grad_t( mDimension );

    Vector dot_u( mDimension );
    double dot_p;
    double dot_t;

    double Div_dot_u;
    double XL, pC = 0;

    /////////////////////////////////////////////////////////////////////////
    //// Integration in space sum_(beta= 0)^(number of quadrature points)
    /////////////////////////////////////////////////////////////////////////

   for ( unsigned int GaussPoint = 0; GaussPoint < integration_points.size(); GaussPoint++ )
    {
        // in element
        noalias( DNu_DX ) = prod( DNu_Dxi[GaussPoint], mInvJ0[GaussPoint] );
        noalias( DNo_DX ) = prod( DNo_Dxi[GaussPoint], mInvJ0[GaussPoint] );
        noalias( Nu ) = row( Nu_container, GaussPoint );
        noalias( No ) = row( No_container, GaussPoint );

        noalias( u ) = Getu( Nu );
        noalias( Grad_u ) = GetGradu( DNu_DX );
        Div_dot_u = GetDivdotu( DNu_DX );

        p = Getp( No );
        dot_p = Getdotp( No );
        noalias( Grad_p ) = GetGradp( DNo_DX );

        t = Gett( No );
        dot_t = Getdott( No );
        noalias( Grad_t ) = GetGradt( DNo_DX );


        XL = GetLiquidSaturation( t, 0 );
        if ( t < 0.0 )
            pC = p + mSf * ( -t );

//         if ( rVariable == DISPLACEMENT )
// 	{
// 	    if ( rValues[GaussPoint].size() != mDimension )
// 		rValues[GaussPoint].resize( mDimension );
// 	    noalias( rValues[GaussPoint] ) = u;
//          }
//         if ( rVariable == WATER_FLOW )
// 	{
// 	    if ( rValues[GaussPoint].size() != mDimension )
// 		rValues[GaussPoint].resize( mDimension );
// 	    noalias( rValues[GaussPoint] ) = GetDarcyFlow( Grad_p, p, t, 0 ) * 24.0 * 3600.0 / mUnitRatio;
//          }
//
//         if ( rVariable == HEAT_FLOW )
// 	{
// 	    if ( rValues[GaussPoint].size() != mDimension )
// 		rValues[GaussPoint].resize( mDimension );
// 	    noalias( rValues[GaussPoint] ) = GetHeatFlow( Grad_t, t, 0 ) ;
//          }

        if ( rVariable == POROSITY )
            rValues[GaussPoint] = GetPorosity( Grad_u, p, t, mT0[GaussPoint], 0 );

        else if ( rVariable == LINEAR_STRAIN )
            rValues[GaussPoint] = GetstrainVol( Grad_u );

        else if ( rVariable == ICE_SATURATION )
            rValues[GaussPoint] = 1.0 - XL;

        else if ( rVariable == ICE_VOLUME_FRACTION )
            rValues[GaussPoint] = GetPorosity( Grad_u, p, t, mT0[GaussPoint], 0 ) * ( 1.0 - XL );

        else if ( rVariable == WATER_DENSITY )
            rValues[GaussPoint] = GetWaterDensity( p, t, 0 ) * pow( mUnitRatio, 3.0 ); // kg/m^3

        else if ( rVariable == WATER_MASS )
            rValues[GaussPoint] = XL * GetPorosity( Grad_u, p, t, mT0[GaussPoint], 0 ) * GetWaterDensity( p, t, 0 ) * pow( mUnitRatio, 3.0 ); //  kg/m^3

        else if ( rVariable == ICE_MASS )
            rValues[GaussPoint] = ( 1 - XL ) * GetPorosity( Grad_u, p, t, mT0[GaussPoint], 0 ) * GetIceDensity( p, t, 0 ) * pow( mUnitRatio, 3.0 ); //  kg/m^3

        else if ( rVariable == WATER_PRESSURE )
            rValues[GaussPoint] = p;

        else if ( rVariable == ICE_PRESSURE )
            rValues[GaussPoint] = pC;

        else  if ( rVariable == EFFECTIVE_STRESS )
            rValues[GaussPoint] = XL * p + ( 1 - XL ) * pC ; // kPa

        else // in constitutive law
            rValues[GaussPoint] = mConstitutiveLawVector[GaussPoint]->GetValue( rVariable, rValues[GaussPoint] );
    }

    KRATOS_CATCH( "" )
}
void FreezingSoil::GetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rProcessInfo )
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
void FreezingSoil::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rProcessInfo ) {}
void FreezingSoil::SetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rProcessInfo ) {}
void FreezingSoil::SetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rProcessInfo ) {}
void FreezingSoil::CalculateOnIntegrationPoints( const Variable<double>& rVariable, Vector& Output, const ProcessInfo& rProcessInfo ) {}
int FreezingSoil::Check( const Kratos::ProcessInfo& rProcessInfo )
{
    return 0;
}


void FreezingSoil::Interpolate( const Variable<double>& rVariable, const ProcessInfo& rProcessInfo )
{
    for ( unsigned int k = 0; k < 4; k++ )
    {
        GetGeometry()[8+k].GetSolutionStepValue( rVariable ) = ( GetGeometry()[k%4].GetSolutionStepValue( rVariable ) + GetGeometry()[( k+1 )%4].GetSolutionStepValue( rVariable ) ) / 2.0;
        GetGeometry()[12+k].GetSolutionStepValue( rVariable ) = ( GetGeometry()[4+k%4].GetSolutionStepValue( rVariable ) + GetGeometry()[4+( k+1 )%4].GetSolutionStepValue( rVariable ) ) / 2.0;
        GetGeometry()[16+k].GetSolutionStepValue( rVariable ) = ( GetGeometry()[k%4].GetSolutionStepValue( rVariable ) + GetGeometry()[4+k%4].GetSolutionStepValue( rVariable ) ) / 2.0;
    }
}

void FreezingSoil::Interpolate( const Variable<Kratos::array_1d<double, 3> >& rVariable, const ProcessInfo& rProcessInfo )
{
    for ( unsigned int k = 0; k < 4; k++ )
    {
        GetGeometry()[8+k].GetSolutionStepValue( rVariable ) = ( GetGeometry()[k%4].GetSolutionStepValue( rVariable ) + GetGeometry()[( k+1 )%4].GetSolutionStepValue( rVariable ) ) / 2.0;
        GetGeometry()[12+k].GetSolutionStepValue( rVariable ) = ( GetGeometry()[4+k%4].GetSolutionStepValue( rVariable ) + GetGeometry()[4+( k+1 )%4].GetSolutionStepValue( rVariable ) ) / 2.0;
        GetGeometry()[16+k].GetSolutionStepValue( rVariable ) = ( GetGeometry()[k%4].GetSolutionStepValue( rVariable ) + GetGeometry()[4+k%4].GetSolutionStepValue( rVariable ) ) / 2.0;
    }
}

/// - WriteNodalResults
void FreezingSoil::FinalizeSolutionStep( ProcessInfo& ProcessInfo )
{
    if ( GetGeometry().size() == 8 || GetGeometry().size() == 20 )
    {
	  for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
	      mConstitutiveLawVector[i]->FinalizeSolutionStep( GetProperties(),
		      GetGeometry(),
		      row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ),
		      ProcessInfo );

	  // Added for reset prestress, prestrain and displacement
	  if ( GetValue( ASSIGN_PRESTRESS_FLAG ) == 1 && !mPrestressAssigned ) // Assign only once!
	  {
	      // Assign prestress on Gauss points
	      const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );
	      const Matrix& No_container = mThisGeometryOther->ShapeFunctionsValues( mThisIntegrationMethod ); // for mT0
	      const GeometryType::ShapeFunctionsGradientsType& DNu_Dxi = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );
	      Vector No( mNodesNumberOther );
	      Matrix DNu_DX( mNodesNumberDisp, mDimension );
	      Matrix Grad_u( mDimension, mDimension );

	      // Assign nodal prestress
	      for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
		  mp0e[i] = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE );

	      for ( unsigned int GaussPoint = 0; GaussPoint < integration_points.size(); GaussPoint++ )
	      {
		  // pre liquid pressure to deal with inital water pressure
		  noalias( No ) = row( No_container, GaussPoint );
		  for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
		      mp0[GaussPoint] += mp0e[i] * No( i );

		  // pre volumetric strain to deal with initial hydrostatic pressure
		  noalias( DNu_DX ) = prod( DNu_Dxi[GaussPoint], mInvJ0[GaussPoint] );
		  noalias( Grad_u ) = GetGradu( DNu_DX );
		  mStrainVol0[GaussPoint] = GetstrainVol( Grad_u ) / 3.0;
		  
      //  mConstitutiveLawVector[GaussPoint]->SetValue( PRESTRESS, -1, ProcessInfo );
	      }

	      // Reset displacement
	      for ( unsigned int i = ( mNodesDispMin - 1 ) ; i < mNodesDispMax ; i++ )
	      {
		  GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_NULL ) = ZeroVector( 3 );
		  GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_EINS ) = ZeroVector( 3 );
		  GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_DT ) = ZeroVector( 3 );
		  GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_NULL_DT ) = ZeroVector( 3 );
		  GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_EINS_DT ) = ZeroVector( 3 );
		  GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_OLD ) = ZeroVector( 3 );
		  GetGeometry()[i].GetSolutionStepValue( ACCELERATION ) = ZeroVector( 3 );
		  GetGeometry()[i].GetSolutionStepValue( ACCELERATION_NULL ) = ZeroVector( 3 );
		  GetGeometry()[i].GetSolutionStepValue( ACCELERATION_EINS ) = ZeroVector( 3 );
	      }
	      //Cancel Gravity
	      noalias( mGravityDisp ) = ZeroVector( mDimension );

	      //Disable the reset of prestress
	      mPrestressAssigned = true;
	  }


	  Matrix DNo_DX( mNodesNumberOther, mDimension );
	  Matrix DNo_Dxi( mNodesNumberOther, mDimension );
	  Matrix nodesLocalCoords( mNodesNumberOther, mDimension );

	  for ( unsigned int i = ( mNodesOtherMin - 1 ) ; i < mNodesOtherMax ;i++ )
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

	  double p, t, pC = 0.0, XL;
	  Vector Grad_p( mDimension ), Grad_t( mDimension );
	  Matrix Grad_u( mDimension, mDimension );

	  for ( unsigned int i = ( mNodesOtherMin - 1 ) ; i < mNodesOtherMax ;i++ )
	  {
	      Vector local_coords( 3 );
	      for ( unsigned int m = 0; m < mDimension; m++ )
		  local_coords( m ) = nodesLocalCoords( i, m );

	      noalias( DNo_Dxi ) = GetGeometry().ShapeFunctionsLocalGradients( DNo_Dxi, local_coords );
	      noalias( DNo_DX ) = prod( DNo_Dxi, mInvJ0[0] );

	      p = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE );
	      t = GetGeometry()[i].GetSolutionStepValue( TEMPERATURE );
	      if ( t < 0.0 )
		  pC = p + mSf * ( -t );
	      noalias( Grad_p ) = GetGradp( DNo_DX );
	      noalias( Grad_t ) = GetGradt( DNo_DX );
	      noalias( Grad_u ) = GetGradu( DNo_DX );

	      XL = GetLiquidSaturation( t, 0 );

	      GetGeometry()[i].GetSolutionStepValue( POROSITY ) = GetPorosity( Grad_u, p, t, mT0e[i], 0 );
	      GetGeometry()[i].GetSolutionStepValue( ICE_SATURATION ) = 1.0 - XL;
	      GetGeometry()[i].GetSolutionStepValue( ICE_VOLUME_FRACTION ) = GetPorosity( Grad_u, p, t, mT0e[i], 0 ) * ( 1.0 - XL );
	      GetGeometry()[i].GetSolutionStepValue( LINEAR_STRAIN ) = GetstrainVol( Grad_u );

	      GetGeometry()[i].GetSolutionStepValue( WATER_DENSITY ) = GetWaterDensity( p, t, 0 ) * pow( mUnitRatio, 3.0 ); // kg/m^3
	      GetGeometry()[i].GetSolutionStepValue( ICE_DENSITY ) = GetIceDensity( p, t, 0 ) * pow( mUnitRatio, 3.0 ); // kg/m^3
	      GetGeometry()[i].GetSolutionStepValue( WATER_FLOW ) =  GetDarcyFlow( Grad_p, p, t, 0 ) * 24.0 * 3600.0 / mUnitRatio; // m/day
	      GetGeometry()[i].GetSolutionStepValue( HEAT_FLOW ) = GetHeatFlow( Grad_t, t, 0 ) ;
	      GetGeometry()[i].GetSolutionStepValue( ICE_PRESSURE ) = pC; // kPa
	      GetGeometry()[i].GetSolutionStepValue( WATER_MASS ) = XL * GetPorosity( Grad_u, p, t, mT0e[i], 0 ) * GetWaterDensity( p, t, 0 ) * pow( mUnitRatio, 3.0 ); //  kg/m^3
	      GetGeometry()[i].GetSolutionStepValue( ICE_MASS ) = ( 1 - XL ) * GetPorosity( Grad_u, p, t, mT0e[i], 0 ) * GetIceDensity( p, t, 0 ) * pow( mUnitRatio, 3.0 ); //  kg/m^3

	      GetGeometry()[i].GetSolutionStepValue( EFFECTIVE_STRESS ) =  XL * p + ( 1 - XL ) * pC ; // kPa
      //         GetGeometry()[i].GetSolutionStepValue( FRICTION_COEFFICIENT ) = GetBiotCoefficient( t, 0 );   // pSkeleton
      //         GetGeometry()[i].GetSolutionStepValue( BULK_MODULUS ) = GetBulkModulus( t, 0 );   // kPa

	  }
	  
	  if ( GetGeometry().size() == 20 )
	  {
	      Interpolate( WATER_PRESSURE, ProcessInfo );
	      Interpolate( TEMPERATURE, ProcessInfo );
	      Interpolate( DISPLACEMENT, ProcessInfo );

	      Interpolate( POROSITY, ProcessInfo );
	      Interpolate( ICE_SATURATION, ProcessInfo );
	      Interpolate( ICE_VOLUME_FRACTION, ProcessInfo );
	      Interpolate( LINEAR_STRAIN, ProcessInfo );
      //
	      Interpolate( WATER_DENSITY, ProcessInfo );
	      Interpolate( ICE_DENSITY, ProcessInfo );
	      Interpolate( WATER_FLOW, ProcessInfo );
	      Interpolate( HEAT_FLOW, ProcessInfo );
	      Interpolate( ICE_PRESSURE, ProcessInfo );
	      Interpolate( WATER_MASS, ProcessInfo );
	      Interpolate( ICE_MASS, ProcessInfo );

      //         Interpolate( FRICTION_COEFFICIENT, ProcessInfo );
      //         Interpolate( BULK_MODULUS, ProcessInfo );
	      Interpolate( EFFECTIVE_STRESS, ProcessInfo );
	  }
    }

}
} // Namespace Kratos
