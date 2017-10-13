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
- Ruhr-University DNo_DXchum, Institute for Structural Mechanics, Germany


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


// System includes
#include <math.h>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/soil_3phase.h"
#include "utilities/math_utils.h"
//#include "includes/constitutive_law.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/tetrahedra_3d_4.h"
#include "freezing_soil_application.h"

namespace Kratos
{
///done
Soil3Phase::Soil3Phase( IndexType NewId, GeometryType::Pointer pGeometry )
        : Element( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

///done
Soil3Phase::Soil3Phase( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
        : Element( NewId, pGeometry, pProperties )
{
    //setting up the nodal degrees of freedom
    //DOFs at the end of time step
    //All calculations are made on the general midpoint alpha
    //Variables DOF_ALPHA are updated in the scheme####
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

    } 
    //select Integration Method for 10 and 4 node tetrahedron
    else if (mNodesDispMax== 4 || mNodesDispMax== 10 )
    { 
        mNodesOtherMax = 4;
        mThisIntegrationMethod = GeometryData::GI_GAUSS_2;//methods for tetrahedra elements
        mThisGeometryOther = Geometry< Node<3> >::Pointer( new Tetrahedra3D4 <Node<3> > (
                                 GetGeometry()( 0 ), GetGeometry()( 1 ), GetGeometry()( 2 ), GetGeometry()( 3 ) ) ); 
				 
    }
    else
        KRATOS_THROW_ERROR( std::logic_error, "This element matches only with a linear hexahedra/tetrahedra (8/3) or quadratic hexahedra/tetrahedra (20,27/10 ) geometry" , *this ); 

}

///done
Element::Pointer Soil3Phase::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new Soil3Phase( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

///done
Soil3Phase::~Soil3Phase()
{
}


void Soil3Phase::ResetConstitutiveLaw()
{
    KRATOS_TRY
    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        mConstitutiveLawVector[i]->ResetMaterial( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************
///done
void Soil3Phase::Initialize()
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
    
    mMaterialParameters = GetProperties()[ELEMENT_PARAMETERS];
    mUnitRatio = GetProperties()[SCALE];			// Unit convertor: from [m] to [mm]
    mrhoS0 = mMaterialParameters[0] / pow(mUnitRatio,3.0); 	// [ kg/m^3 ]
    mrhoL0 = mMaterialParameters[1] / pow(mUnitRatio,3.0);		// [ kg/m^3 ]
    mrhoC0 = mMaterialParameters[2] / pow(mUnitRatio,3.0); 	// [ kg/m^3 ]
    mkS = mMaterialParameters[3] / mUnitRatio; 			// [ Pa = kg/(m*s^2) ]
    mkL = mMaterialParameters[4] / mUnitRatio; 			// [ Pa = kg/(m*s^2) ]
    mkC = mMaterialParameters[5] / mUnitRatio; 			// [ Pa = kg/(m*s^2) ]
    mgS = mMaterialParameters[6] / mUnitRatio; 			// [ Pa = kg/(m*s^2) ]
    mgL = mMaterialParameters[7] / mUnitRatio; 			// [ Pa = kg/(m*s^2) ]
    mgC = mMaterialParameters[8] / mUnitRatio; 			// [ Pa = kg/(m*s^2) ]
    mcS = mMaterialParameters[9] * pow(mUnitRatio,2.0); 		// [ J/(kg*K) = m^2/(s^2*K) ]
    mcL = mMaterialParameters[10] * pow(mUnitRatio,2.0); 		// [ J/(kg*K) = m^2/(s^2*K) ]
    mcC = mMaterialParameters[11] * pow(mUnitRatio,2.0); 		// [ J/(kg*K) = m^2/(s^2*K) ]
    mlambdaS = mMaterialParameters[12] * mUnitRatio; 		// [ W/(m*K) = kg*m/(s^3*K) ]
    mlambdaL = mMaterialParameters[13] * mUnitRatio; 		// [ W/(m*K) = kg*m/(s^3*K) ]
    mlambdaC = mMaterialParameters[14] * mUnitRatio; 		// [ W/(m*K) = kg*m/(s^3*K) ]
    malphaS = mMaterialParameters[15]; 				// [ 1/K ]
    malphaL = mMaterialParameters[16]; 				// [ 1/K ]
    malphaC = mMaterialParameters[17]; 				// [ 1/K ] : depends on temperature:  (54-0.18(-t))x1e-6; see "Frozen Ground" page 53
    mTf = mMaterialParameters[18]; 				// [ K ]
    mSf = mMaterialParameters[19] / mUnitRatio; 			// [ Pa = kg/(m*s^2) ]
    mGravity = GetProperties()[GRAVITY] * mUnitRatio; 			// [m/s^2] 
    mn0 = GetProperties()[POROSITY]; 					// [ - ]
    mtstar = GetProperties()[FIRST_SATURATION_PARAM]; 			// [ K ]
    mm = GetProperties()[SECOND_SATURATION_PARAM]; 			// [ - ]
    mkappa0 = GetProperties()[PERMEABILITY_WATER] * mUnitRatio; 	// [ m/s ] 
 

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
        mT0e[i] = ( GetGeometry()[i] ).GetSolutionStepValue( TEMPERATURE_NULL );
    
    //initial temperature t0 in GaussPoints and InvJ0
    mT0.resize( integration_points.size() );
    noalias( mT0 ) = ZeroVector( integration_points.size() );       
    const Matrix& No_container = mThisGeometryOther->ShapeFunctionsValues( mThisIntegrationMethod ); // for mT0
    for ( unsigned int GaussPoint = 0; GaussPoint < integration_points.size(); GaussPoint++ )
    { 
	//initial temperature t0
	Vector No(mNodesNumberOther);
	noalias(No)= row(No_container,GaussPoint);
	for(unsigned int i = (mNodesOtherMin-1); i < mNodesOtherMax; i++)
	    mT0[GaussPoint] += mT0e[i]*No(i); 

	// InvJ0
	mInvJ0[GaussPoint].resize(mDimension,mDimension,false);
	noalias(mInvJ0[GaussPoint])= ZeroMatrix(mDimension,mDimension);
	
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

///done
void Soil3Phase::CalculateAll( MatrixType& rLeftHandSideMatrix,
                               VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo,
                               bool CalculateStiffnessMatrixFlag, bool CalculateResidualVectorFlag )
{
    KRATOS_TRY
    //resizing as needed the RHS

    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size() != mMatSize )
            rRightHandSideVector.resize( mMatSize );

        noalias( rRightHandSideVector ) = ZeroVector( mMatSize ); //resetting RHS
    }

    //resizing as needed the LHS

    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != mMatSize )
            rLeftHandSideMatrix.resize( mMatSize, mMatSize );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( mMatSize, mMatSize ); //resetting LHS
    }
 
    //reading integration points and local Gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );
    const GeometryType::ShapeFunctionsGradientsType& DNu_De = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );
    const GeometryType::ShapeFunctionsGradientsType& DNo_De = mThisGeometryOther->ShapeFunctionsLocalGradients( mThisIntegrationMethod );
    const Matrix& Nu_container = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );
    const Matrix& No_container = mThisGeometryOther->ShapeFunctionsValues( mThisIntegrationMethod );

    Vector Help_R_1( mMatSizeU );
    Vector Help_R_2( mMatSizeO );
    Vector Help_R_3( mMatSizeO );

    Vector Help_R_1_x( mMatSizeU );
    Vector Help_R_2_x( mMatSizeO );
    Vector Help_R_3_x( mMatSizeO );
    
//     Matrix Help_K_TT( mMatSizeO, mMatSizeO ); 
//     noalias( Help_K_TT ) = ZeroMatrix( mMatSizeO, mMatSizeO );

    Matrix DNu_DX( mNodesNumberDisp, mDimension ); // DNu_DX= DNu_DX
    Matrix DNo_DX( mNodesNumberOther, mDimension ); // DNo_DX= DNo_DX
    double weight;
    Vector Nu( mNodesNumberDisp );
    Vector No( mNodesNumberOther );
    double DetJ = 0.0; 

    // nodal solution of primary variable set
    Vector ue = ZeroVector( mMatSize );
    Vector ue_x = ZeroVector( mMatSize );
    Vector dot_ue = ZeroVector( mMatSize );
    Vector epsilon = ZeroVector( mMatSize );

    for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
    {
        ue[i*mDimension] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X );
        ue[i*mDimension+1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y );
        ue[i*mDimension+2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z );

        dot_ue[i*mDimension] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_DT )[0];
        dot_ue[i*mDimension+1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_DT )[1];
        dot_ue[i*mDimension+2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_DT )[2];
    }

    for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
    {
        ue[mAddIndexU+i*mNumberO] = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE );
        ue[mAddIndexU+i*mNumberO+1] = GetGeometry()[i].GetSolutionStepValue( TEMPERATURE );

        dot_ue[mAddIndexU+i*mNumberO] = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE_DT );
        dot_ue[mAddIndexU+i*mNumberO+1] = GetGeometry()[i].GetSolutionStepValue( TEMPERATURE_DT );
    }

    // vector ue: { ux_n1, uy_n1, uz_n1, ux_n2, uy_n2, uz_n2, ... n_n1, p_n1, t_n1, n_n2, p_n2, t_n2, ...}
    Vector u( mDimension ); // somlambdaCd displacement
    double p; // relative water pressure
    double t; // temperature in Celsius degree

    Vector u_x( mDimension );
    double p_x;
    double t_x;

    Matrix Grad_u( mDimension, mDimension );
    Vector Grad_p( mDimension );
    Vector Grad_t( mDimension );

    Matrix Grad_u_x( mDimension, mDimension );
    Vector Grad_p_x( mDimension );
    Vector Grad_t_x( mDimension );

    Vector dot_u( mDimension );
    double dot_p;
    double dot_t;

    double Div_dot_u; 

    // adding for calling constitutive law
    Vector InputVector( 9 );
    Vector strainVector( 6 ); 
    Vector stressVectorEff( 6 );
    noalias( stressVectorEff ) = ZeroVector( 6 );  
    Matrix CtanEff( 6, 6 );
    noalias( CtanEff ) = ZeroMatrix( 6, 6 );
    
    // ++++++ RightHandSideVector: R ++++++
    if ( CalculateResidualVectorFlag == true )
    {
        noalias( Help_R_1 ) = ZeroVector( mMatSizeU );
        noalias( Help_R_2 ) = ZeroVector( mMatSizeO );
        noalias( Help_R_3 ) = ZeroVector( mMatSizeO );

	
        // START GAUSS INTEGRATION: 
        for ( unsigned int GaussPoint = 0; GaussPoint < integration_points.size(); GaussPoint++ )
        {
            // setting DNo_DX, DNo_DX, Nu, No for each integration point
            noalias( DNu_DX ) = prod( DNu_De[GaussPoint], mInvJ0[GaussPoint] );
            noalias( DNo_DX ) = prod( DNo_De[GaussPoint], mInvJ0[GaussPoint] );

            noalias( Nu ) = row( Nu_container, GaussPoint );
            noalias( No ) = row( No_container, GaussPoint );

            weight = integration_points[GaussPoint].Weight();
            DetJ = mDetJ0[GaussPoint]; 

            u = ZeroVector( mDimension );
            p = 0.0; // relative water pressure
            t = 0.0; // temperature in Celsius degree

            Grad_u = ZeroMatrix( mDimension, mDimension );
            Grad_p = ZeroVector( mDimension );
            Grad_t = ZeroVector( mDimension );

            dot_u = ZeroVector( mDimension );
            dot_p = 0.0;
            dot_t = 0.0;

            Div_dot_u = 0.0;

            for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
            {
                for ( unsigned int k = 0; k < mDimension; k++ )
                {
                    u( k ) += ue( i * mDimension + k ) * Nu( i );
                    dot_u( k ) += dot_ue( i * mDimension + k ) * Nu( i );
                    Div_dot_u += dot_ue( i * mDimension + k ) * DNu_DX( i, k );

                    for ( unsigned int l = 0; l < mDimension; l++ )
                        Grad_u( k, l ) += ue( i * mDimension + k ) * DNu_DX( i, l );
                }
            }

            for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
            {
                p += ue( mAddIndexU + i * mNumberO ) * No( i );
                t += ue( mAddIndexU + i * mNumberO + 1 ) * No( i );

                dot_p += dot_ue( mAddIndexU + i * mNumberO ) * No( i );
                dot_t += dot_ue( mAddIndexU + i * mNumberO + 1 ) * No( i );

                for ( unsigned int l = 0; l < mDimension; l++ )
                {
                    Grad_p( l ) += ue( mAddIndexU + i * mNumberO ) * DNo_DX( i, l );
                    Grad_t( l ) += ue( mAddIndexU + i * mNumberO + 1 ) * DNo_DX( i, l );
                }
            }
            
	    if ( GetValue( PLASTIC_FLAG ) != 0 ) // PLATIC  --> ATTENTION 20131030: NOT WORKING PROPERLY, NEED TO BE CORRECTED!!
	    { 
		noalias( strainVector ) = GetstrainVector( Grad_u );   // Compressive ! 
		/// METHOD 2: mBBM model
		for ( unsigned int i = 0; i < 6; i++ )
		    InputVector[i] = strainVector[i]; 
		InputVector[6] = p;
		InputVector[7] = t;
// 		InputVector[8] = GetPorosity( Grad_u, p, t, mT0[GaussPoint], 0 );

// 		bool onlyGP1 = false;
// 		if ( GaussPoint == 1 && GetValue( KRATOS_WATCH_FLAG ) == 1 )
// 		    onlyGP1 = true;

		// INPUT: strainVector (Tensile-Positive !!);  OUTPUT: stressVector [kPa] (Tensile-Positive !!) , CtanEff [kPa]
		mConstitutiveLawVector[GaussPoint]->CalculateMaterialResponse( InputVector, ZeroMatrix( 1, 1 ), stressVectorEff, CtanEff, rCurrentProcessInfo, GetProperties(), GetGeometry(), Nu, true, 0, true ); //last term: update p0
	    }
       
if(GetValue(KRATOS_WATCH_FLAG)==1 && GaussPoint == 1)
{
std::cout << "1---- Div_dot_u= "<<Div_dot_u<<",\t Grad_u= "<<Grad_u<< std::endl;
std::cout << "2---- p= "<<p<<",\t dot_p= "<<dot_p<<",\t Grad_p= "<<Grad_p<< std::endl;
std::cout << "3---- t0= "<<mT0[GaussPoint]<<",\t t= "<<t<<",\t dot_t= "<<dot_t<<",\t Grad_t= "<<Grad_t<< std::endl;
KRATOS_WATCH(stressVectorEff);
}
	    
            //Calculation of spatial load vector
            AddInternalForcesToRHS1( Help_R_1, Nu, DNu_DX, weight, DetJ, mT0[GaussPoint], u, p, t, Grad_u, Grad_p, Grad_t, dot_u, dot_p, dot_t, Div_dot_u, stressVectorEff );
            AddInternalForcesToRHS2( Help_R_2, No, DNo_DX, weight, DetJ, mT0[GaussPoint], u, p, t, Grad_u, Grad_p, Grad_t, dot_u, dot_p, dot_t, Div_dot_u );
            AddInternalForcesToRHS3( Help_R_3, No, DNo_DX, weight, DetJ, mT0[GaussPoint], u, p, t, Grad_u, Grad_p, Grad_t, dot_u, dot_p, dot_t, Div_dot_u );
        }// END GAUSS INTEGRATION.

if(GetValue(KRATOS_WATCH_FLAG)==1)
{
// KRATOS_WATCH(Help_R_1);
// KRATOS_WATCH(Help_R_2);
// KRATOS_WATCH(Help_R_3);
}

        // Assemble the RightHandSideVector R
        for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
            for ( unsigned int k = 0; k < mDimension; k++ )
                rRightHandSideVector( i*mDimension + k ) -= Help_R_1( i * mDimension + k );

        for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
        {
            rRightHandSideVector( mAddIndexU + i*mNumberO ) -= Help_R_2( i );
            rRightHandSideVector( mAddIndexU + i*mNumberO + 1 ) -= Help_R_3( i );
        }
    }

    // ++++++ LeftHandSideMatrix: K_col ++++++
    if ( CalculateStiffnessMatrixFlag == true )
    {
        double eps = GetProperties()[DP_EPSILON];

        for ( unsigned int i = 0; i < mMatSize; i++ )
        {
            epsilon[i] = ue[i] * eps;

            if ( fabs( epsilon[i] ) < eps )
                epsilon[i] = eps;

            ue_x[i] = ue[i] + epsilon[i];
        }

        // START LOOP K_col
        Vector ue_x_col( mMatSize );

        for ( unsigned int col = 0; col < mMatSize; col++ )
        {
            noalias( Help_R_1_x ) = ZeroVector( mMatSizeU );
            noalias( Help_R_2_x ) = ZeroVector( mMatSizeO );
            noalias( Help_R_3_x ) = ZeroVector( mMatSizeO );
            ue_x_col =  ZeroVector( mMatSize );

            for ( unsigned int i = 0; i < mMatSize; i++ )
            {
                if ( i == col )
                    ue_x_col[i] = ue_x[i];
                else
                    ue_x_col[i] = ue[i];
            }

            // START GAUSS INTEGRATION:

            for ( unsigned int GaussPoint = 0; GaussPoint < integration_points.size(); GaussPoint++ )
            {
                // setting DNo_DX, DNo_DX, Nu, No for each integration point
                noalias( DNu_DX ) = prod( DNu_De[GaussPoint], mInvJ0[GaussPoint] );
                noalias( DNo_DX ) = prod( DNo_De[GaussPoint], mInvJ0[GaussPoint] );

                noalias( Nu ) = row( Nu_container, GaussPoint );
                noalias( No ) = row( No_container, GaussPoint );

                weight = integration_points[GaussPoint].Weight();
                DetJ = mDetJ0[GaussPoint]; 

                u_x = ZeroVector( mDimension );
                p_x = 0.0;
                t_x = 0.0;

                Grad_u_x = ZeroMatrix( mDimension, mDimension );
                Grad_p_x = ZeroVector( mDimension );
                Grad_t_x = ZeroVector( mDimension );

                dot_u = ZeroVector( mDimension );
                dot_p = 0.0;
                dot_t = 0.0;

                Div_dot_u = 0.0;

                for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
                {
                    for ( unsigned int k = 0; k < mDimension; k++ )
                    {
                        u_x( k ) += ue_x_col( i * mDimension + k ) * Nu( i );
                        dot_u( k ) += dot_ue( i * mDimension + k ) * Nu( i );
                        Div_dot_u += dot_ue( i * mDimension + k ) * DNu_DX( i, k );

                        for ( unsigned int l = 0; l < mDimension; l++ )
                        {
                            Grad_u_x( k, l ) += ue_x_col( i * mDimension + k ) * DNu_DX( i, l );
                        }
                    }
                }

                for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
                {
                    p_x += ue_x_col( mAddIndexU + i * mNumberO ) * No( i );
                    t_x += ue_x_col( mAddIndexU + i * mNumberO + 1 ) * No( i ); 

                    dot_p += dot_ue( mAddIndexU + i * mNumberO ) * No( i );
                    dot_t += dot_ue( mAddIndexU + i * mNumberO + 1 ) * No( i );

                    for ( unsigned int l = 0; l < mDimension; l++ )
                    {
                        Grad_p_x( l ) += ue_x_col( mAddIndexU + i * mNumberO ) * DNo_DX( i, l );
                        Grad_t_x( l ) += ue_x_col( mAddIndexU + i * mNumberO + 1 ) * DNo_DX( i, l );
                    }
                }

		if ( GetValue( PLASTIC_FLAG ) != 0 ) // ELATIC 
		{ 
		    noalias( strainVector ) = GetstrainVector( Grad_u_x );   // Compressive ! 
		    /// METHOD 2: mBBM model
		    for ( unsigned int i = 0; i < 6; i++ )
			InputVector[i] = strainVector[i]; 
		    InputVector[6] = p_x;
		    InputVector[7] = t_x;
    // 		InputVector[8] = GetPorosity( Grad_u, p, t, mT0[GaussPoint], 0 );
 
		    // INPUT: strainVector (Tensile-Positive !!);  OUTPUT: stressVector [kPa] (Tensile-Positive !!) , CtanEff [kPa]
		    mConstitutiveLawVector[GaussPoint]->CalculateMaterialResponse( InputVector, ZeroMatrix( 1, 1 ), stressVectorEff, CtanEff, rCurrentProcessInfo, GetProperties(), GetGeometry(), Nu, true, 0, false );
		}
		
                AddInternalForcesToRHS1( Help_R_1_x, Nu, DNu_DX, weight, DetJ, mT0[GaussPoint], u_x, p_x, t_x, Grad_u_x, Grad_p_x, Grad_t_x, dot_u, dot_p, dot_t, Div_dot_u, stressVectorEff);
                AddInternalForcesToRHS2( Help_R_2_x, No, DNo_DX, weight, DetJ, mT0[GaussPoint], u_x, p_x, t_x, Grad_u_x, Grad_p_x, Grad_t_x, dot_u, dot_p, dot_t, Div_dot_u );
                AddInternalForcesToRHS3( Help_R_3_x, No, DNo_DX, weight, DetJ, mT0[GaussPoint], u_x, p_x, t_x, Grad_u_x, Grad_p_x, Grad_t_x, dot_u, dot_p, dot_t, Div_dot_u );
            }// END GAUSS INTEGRATION.

            for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
                for ( unsigned int k = 0; k < mDimension; k++ )
                    rLeftHandSideMatrix( i*mDimension + k, col ) = ( Help_R_1_x[i*mDimension+k] - Help_R_1[i*mDimension+k] ) / epsilon[col];
		
            for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
            {
                rLeftHandSideMatrix( mAddIndexU + i*mNumberO, col ) = ( Help_R_2_x[i] - Help_R_2[i] ) / epsilon[col];
                rLeftHandSideMatrix( mAddIndexU + i*mNumberO + 1, col ) = ( Help_R_3_x[i] - Help_R_3[i] ) / epsilon[col];
            }

        }// END LOOP K_col
        
if ( GetValue( KRATOS_WATCH_FLAG ) == 1 )
{
//     Matrix Help_K_UU( mMatSizeU, mMatSizeU );
//     Matrix Help_K_UP( mMatSizeU, mMatSizeO );
//     Matrix Help_K_UT( mMatSizeU, mMatSizeO );
//     Matrix Help_K_PU( mMatSizeO, mMatSizeU );
//     Matrix Help_K_PP( mMatSizeO, mMatSizeO );
//     Matrix Help_K_PT( mMatSizeO, mMatSizeO );
//     Matrix Help_K_TU( mMatSizeO, mMatSizeU );
//     Matrix Help_K_TP( mMatSizeO, mMatSizeO );
//     Matrix Help_K_TT( mMatSizeO, mMatSizeO );
//     noalias( Help_K_UU ) = ZeroMatrix( mMatSizeU, mMatSizeU );
//     noalias( Help_K_UP ) = ZeroMatrix( mMatSizeU, mMatSizeO );
//     noalias( Help_K_UT ) = ZeroMatrix( mMatSizeU, mMatSizeO );
//     noalias( Help_K_PU ) = ZeroMatrix( mMatSizeO, mMatSizeU );
//     noalias( Help_K_PP ) = ZeroMatrix( mMatSizeO, mMatSizeO );
//     noalias( Help_K_PT ) = ZeroMatrix( mMatSizeO, mMatSizeO );
//     noalias( Help_K_TU ) = ZeroMatrix( mMatSizeO, mMatSizeU );
//     noalias( Help_K_TP ) = ZeroMatrix( mMatSizeO, mMatSizeO );
//     noalias( Help_K_TT ) = ZeroMatrix( mMatSizeO, mMatSizeO );
// 
//     for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
//         for ( unsigned int k = 0; k < mDimension; k++ )
//         {
//             for ( unsigned int j = ( mNodesDispMin - 1 ); j < mNodesDispMax; j++ )
//                 for ( unsigned int l = 0; l < mDimension; l++ )
//                     Help_K_UU( i * mDimension + k, j * mDimension + l ) = rLeftHandSideMatrix( i * mDimension + k, j * mDimension + l );
// 
//             for ( unsigned int j = ( mNodesOtherMin - 1 ); j < mNodesOtherMax; j++ )
//             {
//                 Help_K_UP( i * mDimension + k, j ) = rLeftHandSideMatrix( i * mDimension + k, mAddIndexU + j * mNumberO ) ;
//                 Help_K_UT( i * mDimension + k, j ) = rLeftHandSideMatrix( i * mDimension + k, mAddIndexU + j * mNumberO + 1 )  ;
//             }
//         }
// 
//     for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
//     {
//         for ( unsigned int j = ( mNodesDispMin - 1 ); j < mNodesDispMax; j++ )
//             for ( unsigned int l = 0; l < mDimension; l++ )
//             {
//                 Help_K_PU( i, j * mDimension + l ) = rLeftHandSideMatrix( mAddIndexU + i * mNumberO, j * mDimension + l );
//                 Help_K_TU( i, j * mDimension + l ) = rLeftHandSideMatrix( mAddIndexU + i * mNumberO + 1, j * mDimension + l );
//             }
// 
//         for ( unsigned int j = ( mNodesOtherMin - 1 ); j < mNodesOtherMax; j++ )
//         {
//             Help_K_PP( i, j ) = rLeftHandSideMatrix( mAddIndexU + i * mNumberO, mAddIndexU + j * mNumberO );
//             Help_K_PT( i, j ) = rLeftHandSideMatrix( mAddIndexU + i * mNumberO, mAddIndexU + j * mNumberO + 1 );
//             Help_K_TP( i, j ) = rLeftHandSideMatrix( mAddIndexU + i * mNumberO + 1, mAddIndexU + j * mNumberO );
//             Help_K_TT( i, j ) = rLeftHandSideMatrix( mAddIndexU + i * mNumberO + 1, mAddIndexU + j * mNumberO + 1 );
//         }
//     }
// 
//     KRATOS_WATCH( Help_K_UU );
// //     KRATOS_WATCH( Help_K_UP );
//     KRATOS_WATCH( Help_K_UT );
//     
// //     KRATOS_WATCH( Help_K_PU );
// //     KRATOS_WATCH( Help_K_PP );
// //     KRATOS_WATCH( Help_K_PT ); 
//     
//     KRATOS_WATCH( Help_K_TU );
// //     KRATOS_WATCH( Help_K_TP );
//     KRATOS_WATCH( Help_K_TT ); 
}
    }
    
	
    KRATOS_CATCH( "" )
}

///done
void Soil3Phase::DampMatrix( MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    if ( rDampMatrix.size1() != mMatSize )
        rDampMatrix.resize( mMatSize, mMatSize );

    noalias( rDampMatrix ) = ZeroMatrix( mMatSize, mMatSize );
 
    //reading integration points and local Gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    const GeometryType::ShapeFunctionsGradientsType& DNu_De = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );

    const GeometryType::ShapeFunctionsGradientsType& DNo_De = mThisGeometryOther->ShapeFunctionsLocalGradients( mThisIntegrationMethod );

    const Matrix& Nu_container = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

    const Matrix& No_container = mThisGeometryOther->ShapeFunctionsValues( mThisIntegrationMethod );

    Vector Help_R_1( mMatSizeU );
    Vector Help_R_2( mMatSizeO );
    Vector Help_R_3( mMatSizeO );

    Vector Help_R_1_x( mMatSizeU );
    Vector Help_R_2_x( mMatSizeO );
    Vector Help_R_3_x( mMatSizeO );

    Matrix Help_D_TT( mMatSizeO, mMatSizeO ); 
    noalias( Help_D_TT ) = ZeroMatrix( mMatSizeO, mMatSizeO );
    
    Matrix DNu_DX( mNodesNumberDisp, mDimension );
    Matrix DNo_DX( mNodesNumberOther, mDimension );
    Vector Nu( mNodesNumberDisp );
    Vector No( mNodesNumberOther );

    double weight; 
    double DetJ = 0.0; 

    // nodal solution of primary variable set
    Vector ue = ZeroVector( mMatSize );
    Vector dot_ue = ZeroVector( mMatSize );
    Vector dot_ue_x = ZeroVector( mMatSize );
    Vector dot_epsilon = ZeroVector( mMatSize );

    for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
    {
        ue[i*mDimension] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X );
        ue[i*mDimension+1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y );
        ue[i*mDimension+2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z );

        dot_ue[i*mDimension] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_DT )[0];
        dot_ue[i*mDimension+1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_DT )[1];
        dot_ue[i*mDimension+2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_DT )[2];
    }

    for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
    {
        ue[mAddIndexU+i*mNumberO] = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE );
        ue[mAddIndexU+i*mNumberO+1] = GetGeometry()[i].GetSolutionStepValue( TEMPERATURE );

        dot_ue[mAddIndexU+i*mNumberO] = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE_DT );
        dot_ue[mAddIndexU+i*mNumberO+1] = GetGeometry()[i].GetSolutionStepValue( TEMPERATURE_DT );
    }

    Vector u( mDimension ); 
    double p;
    double t;

    Matrix Grad_u( mDimension, mDimension );
    Vector Grad_p( mDimension );
    Vector Grad_t( mDimension );

    Vector dot_u( mDimension );
    double dot_p;
    double dot_t;

    Vector dot_u_x( mDimension );
    double dot_p_x;
    double dot_t_x;

    double Div_dot_u;
    double Div_dot_u_x;

    // ++++++ RightHandSideVector: R ++++++
    noalias( Help_R_1 ) = ZeroVector( mMatSizeU );
    noalias( Help_R_2 ) = ZeroVector( mMatSizeO );
    noalias( Help_R_3 ) = ZeroVector( mMatSizeO );

    // START GAUSS INTEGRATION:

    for ( unsigned int GaussPoint = 0; GaussPoint < integration_points.size(); GaussPoint++ )
    {
        // setting DNo_DX, DNo_DX, Nu, No for each integration point
        noalias( DNu_DX ) = prod( DNu_De[GaussPoint], mInvJ0[GaussPoint] );
        noalias( DNo_DX ) = prod( DNo_De[GaussPoint], mInvJ0[GaussPoint] );

        noalias( Nu ) = row( Nu_container, GaussPoint );
        noalias( No ) = row( No_container, GaussPoint );

        weight = integration_points[GaussPoint].Weight();
        DetJ = mDetJ0[GaussPoint]; 

        u = ZeroVector( mDimension );
        p = 0.0; // relative water pressure
        t = 0.0; // temperature in Celsius degree

        Grad_u = ZeroMatrix( mDimension, mDimension );
        Grad_p = ZeroVector( mDimension );
        Grad_t = ZeroVector( mDimension );

        dot_u = ZeroVector( mDimension );
        dot_p = 0.0;
        dot_t = 0.0;

        Div_dot_u = 0.0;

        for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
        {
            for ( unsigned int k = 0; k < mDimension; k++ )
            {
                u( k ) += ue( i * mDimension + k ) * Nu( i );
                dot_u( k ) += dot_ue( i * mDimension + k ) * Nu( i );
                Div_dot_u += dot_ue( i * mDimension + k ) * DNu_DX( i, k );

                for ( unsigned int l = 0; l < mDimension; l++ )
                    Grad_u( k, l ) += ue( i * mDimension + k ) * DNu_DX( i, l );
            }
        }

        for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
        {
            p += ue( mAddIndexU + i * mNumberO ) * No( i );
            t += ue( mAddIndexU + i * mNumberO + 1 ) * No( i ); 

            dot_p += dot_ue( mAddIndexU + i * mNumberO ) * No( i );
            dot_t += dot_ue( mAddIndexU + i * mNumberO + 1 ) * No( i );

            for ( unsigned int l = 0; l < mDimension; l++ )
            {
                Grad_p( l ) += ue( mAddIndexU + i * mNumberO ) * DNo_DX( i, l );
                Grad_t( l ) += ue( mAddIndexU + i * mNumberO + 1 ) * DNo_DX( i, l );
            }
        }

        //Calculation of spatial load vector
//         AddInternalForcesToRHS1( Help_R_1, Nu, DNu_DX, weight, DetJ, mT0[GaussPoint], u, p, t, Grad_u, Grad_p, Grad_t, dot_u, dot_p, dot_t, Div_dot_u );
        AddInternalForcesToRHS2( Help_R_2, No, DNo_DX, weight, DetJ, mT0[GaussPoint], u, p, t, Grad_u, Grad_p, Grad_t, dot_u, dot_p, dot_t, Div_dot_u );
        AddInternalForcesToRHS3( Help_R_3, No, DNo_DX, weight, DetJ, mT0[GaussPoint], u, p, t, Grad_u, Grad_p, Grad_t, dot_u, dot_p, dot_t, Div_dot_u );
    }// END GAUSS INTEGRATION.

    // ++++++ DampingMatrix: D_col ++++++
    double eps = GetProperties()[DP_EPSILON];

    for ( unsigned int i = 0; i < mMatSize; i++ )
    {
        dot_epsilon[i] = dot_ue[i] * eps;

        if ( fabs( dot_epsilon[i] ) < eps )
            dot_epsilon[i] = eps;

        dot_ue_x[i] = dot_ue[i] + dot_epsilon[i];
    }

    // START LOOP K_col
    Vector dot_ue_x_col( mMatSize );

    for ( unsigned int col = 0; col < mMatSize; col++ )
    {
        noalias( Help_R_1_x ) = ZeroVector( mMatSizeU );
        noalias( Help_R_2_x ) = ZeroVector( mMatSizeO );
        noalias( Help_R_3_x ) = ZeroVector( mMatSizeO );

        dot_ue_x_col =  ZeroVector( mMatSize );

        for ( unsigned int i = 0; i < mMatSize; i++ )
        {
            if ( i == col )
                dot_ue_x_col[i] = dot_ue_x[i];
            else
                dot_ue_x_col[i] = dot_ue[i];
        }

        // START GAUSS INTEGRATION:

        for ( unsigned int GaussPoint = 0; GaussPoint < integration_points.size(); GaussPoint++ )
        {
            // setting DNo_DX, DNo_DX, Nu, No for each integration point
            noalias( DNu_DX ) = prod( DNu_De[GaussPoint], mInvJ0[GaussPoint] );
            noalias( DNo_DX ) = prod( DNo_De[GaussPoint], mInvJ0[GaussPoint] );
            noalias( Nu ) = row( Nu_container, GaussPoint );
            noalias( No ) = row( No_container, GaussPoint );

            weight = integration_points[GaussPoint].Weight();
            DetJ = mDetJ0[GaussPoint]; 

            u = ZeroVector( mDimension );
            p = 0.0;
            t = 0.0;

            Grad_u = ZeroMatrix( mDimension, mDimension );
            Grad_p = ZeroVector( mDimension );
            Grad_t = ZeroVector( mDimension );

            dot_u_x = ZeroVector( mDimension );
            dot_p_x = 0.0;
            dot_t_x = 0.0;

            Div_dot_u_x = 0.0;

            for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
            {
                for ( unsigned int k = 0; k < mDimension; k++ )
                {
                    u( k ) += ue( i * mDimension + k ) * Nu( i );
                    dot_u_x( k ) += dot_ue_x_col( i * mDimension + k ) * Nu( i );
                    Div_dot_u_x += dot_ue_x_col( i * mDimension + k ) * DNu_DX( i, k );

                    for ( unsigned int l = 0; l < mDimension; l++ )
                        Grad_u( k, l ) += ue( i * mDimension + k ) * DNu_DX( i, l );
                }
            }

            for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
            {
                p += ue( mAddIndexU + i * mNumberO ) * No( i );
                t += ue( mAddIndexU + i * mNumberO + 1 ) * No( i );

                dot_p_x += dot_ue_x_col( mAddIndexU + i * mNumberO ) * No( i );
                dot_t_x += dot_ue_x_col( mAddIndexU + i * mNumberO + 1 ) * No( i );

                for ( unsigned int l = 0; l < mDimension; l++ )
                {
                    Grad_p( l ) += ue( mAddIndexU + i * mNumberO ) * DNo_DX( i, l );
                    Grad_t( l ) += ue( mAddIndexU + i * mNumberO + 1 ) * DNo_DX( i, l );
                }
            }

//             AddInternalForcesToRHS1( Help_R_1_x, Nu, DNu_DX, weight, DetJ, mT0[GaussPoint], u, p, t, Grad_u, Grad_p, Grad_t, dot_u_x, dot_p_x, dot_t_x, Div_dot_u_x );
            AddInternalForcesToRHS2( Help_R_2_x, No, DNo_DX, weight, DetJ, mT0[GaussPoint], u, p, t, Grad_u, Grad_p, Grad_t, dot_u_x, dot_p_x, dot_t_x, Div_dot_u_x );
            AddInternalForcesToRHS3( Help_R_3_x, No, DNo_DX, weight, DetJ, mT0[GaussPoint], u, p, t, Grad_u, Grad_p, Grad_t, dot_u_x, dot_p_x, dot_t_x, Div_dot_u_x );
	    
        }// END GAUSS INTEGRATION.

        for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
            for ( unsigned int k = 0; k < mDimension; k++ )
//                 rDampMatrix( i*mDimension + k, col ) = ( Help_R_1_x[i*mDimension+k] - Help_R_1[i*mDimension+k] ) / dot_epsilon[col];
                rDampMatrix( i*mDimension + k, col ) = 0.0;

        for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
        {
            rDampMatrix( mAddIndexU + i*mNumberO, col ) = ( Help_R_2_x[i] - Help_R_2[i] ) / dot_epsilon[col];
            rDampMatrix( mAddIndexU + i*mNumberO + 1, col ) = ( Help_R_3_x[i] - Help_R_3[i] ) / dot_epsilon[col];
        }  
    }// END LOOP K_col
 

if ( GetValue( KRATOS_WATCH_FLAG ) == 1 )
{ 
//     Matrix Help_D_PU( mMatSizeO, mMatSizeU );
//     Matrix Help_D_PP( mMatSizeO, mMatSizeO );
//     Matrix Help_D_PT( mMatSizeO, mMatSizeO );
//     Matrix Help_D_TU( mMatSizeO, mMatSizeU );
//     Matrix Help_D_TP( mMatSizeO, mMatSizeO );
//     Matrix Help_D_TT( mMatSizeO, mMatSizeO ); 
//     noalias( Help_D_PU ) = ZeroMatrix( mMatSizeO, mMatSizeU );
//     noalias( Help_D_PP ) = ZeroMatrix( mMatSizeO, mMatSizeO );
//     noalias( Help_D_PT ) = ZeroMatrix( mMatSizeO, mMatSizeO );
//     noalias( Help_D_TU ) = ZeroMatrix( mMatSizeO, mMatSizeU );
//     noalias( Help_D_TP ) = ZeroMatrix( mMatSizeO, mMatSizeO );
//     noalias( Help_D_TT ) = ZeroMatrix( mMatSizeO, mMatSizeO ); 
// 
//     for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
//     {
//         for ( unsigned int j = ( mNodesDispMin - 1 ); j < mNodesDispMax; j++ )
//             for ( unsigned int l = 0; l < mDimension; l++ )
//             {
//                 Help_D_PU( i, j * mDimension + l ) = rDampMatrix( mAddIndexU + i * mNumberO, j * mDimension + l );
//                 Help_D_TU( i, j * mDimension + l ) = rDampMatrix( mAddIndexU + i * mNumberO + 1, j * mDimension + l );
//             }
// 
//         for ( unsigned int j = ( mNodesOtherMin - 1 ); j < mNodesOtherMax; j++ )
//         {
//             Help_D_PP( i, j ) = rDampMatrix( mAddIndexU + i * mNumberO, mAddIndexU + j * mNumberO );
//             Help_D_PT( i, j ) = rDampMatrix( mAddIndexU + i * mNumberO, mAddIndexU + j * mNumberO + 1 );
//             Help_D_TP( i, j ) = rDampMatrix( mAddIndexU + i * mNumberO + 1, mAddIndexU + j * mNumberO );
//             Help_D_TT( i, j ) = rDampMatrix( mAddIndexU + i * mNumberO + 1, mAddIndexU + j * mNumberO + 1 );
//         }
//     }
 
//     KRATOS_WATCH( Help_D_PU );
//     KRATOS_WATCH( Help_D_PP );
//     KRATOS_WATCH( Help_D_PT ); 
    
//     KRATOS_WATCH( Help_D_TU );
//     KRATOS_WATCH( Help_D_TP );
//     KRATOS_WATCH( Help_D_TT ); 
}

    KRATOS_CATCH( "" )
}

///done
Soil3Phase::IntegrationMethod Soil3Phase::GetIntegrationMethod()
{
    return mThisIntegrationMethod;
}

///done
void Soil3Phase::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag,  CalculateResidualVectorFlag );
}

///done
void Soil3Phase::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    //calculation flamgS
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;

    CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );

}

///done
void Soil3Phase::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo )
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

///done
void Soil3Phase::GetDofList( DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo )
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

///done
void Soil3Phase::GetValuesVector( Vector& values, int Step )
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

/// - Initialize Material
void Soil3Phase::InitializeMaterial() {}

void Soil3Phase::AssembleTimeSpaceRHSFromSubVectors( VectorType& rRightHandSideVector, const Vector& R_1, const Vector& R_2, const Vector& R_3 )
{
    KRATOS_TRY

    for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
        for ( unsigned int m = 0; m < mDimension; m++ )
            rRightHandSideVector( i*mDimension + m ) -= R_1( i * mDimension + m );

    for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
    {
        rRightHandSideVector( mMatSizeU + i ) -= R_2( i );
        rRightHandSideVector( mMatSizeU + mMatSizeO + i ) -= R_3( i );
    }

    KRATOS_CATCH( "" )
}

///-//////////////////////////////////////////////////////////////////////////////////////
///-//////////////////////////////////////////////////////////////////////////////////////
/// R1 **************************
void Soil3Phase::AddInternalForcesToRHS1( Vector& Help_R_1, const Vector& Nu, const Matrix& DNu_DX, double weight, double DetJ, double t0, Vector u, double p, double t, Matrix Grad_u, Vector Grad_p, Vector Grad_t, Vector dot_u, double dot_p, double dot_t, double Div_dot_u, Vector& stressVectorEff )
{ 
	Matrix stressEff( mDimension, mDimension );
	noalias( stressEff ) = ZeroMatrix( mDimension, mDimension );
	      
        if ( GetValue( PLASTIC_FLAG ) == 0 ) // ELATIC   
	      noalias( stressEff ) = GetstressEff( Grad_u, p, t );  
    
	double MLC = GetMLC( Grad_u, p, t, t0 );
	double XL = GetXL( t );
	double b = Getb( t );
	double K = GetK( t ); 
	double n = Getn( Grad_u, p, t, t0 );
 
//     double pC = 0.0;
//     double pS = -(stressEff(0,0)+stressEff(1,1)+stressEff(2,2) )/3.0;
//     if(t<0.0)
//     {
//       pC = p + mSf*(-t);
// //       if(pC>pS)
// // 	pC = pS;
//     }
    
    for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
        for ( unsigned int m = 0; m < mDimension; m++ )
        {
            Help_R_1( i*mDimension + m ) += Nu( i ) * (
                                                (( 1 - n) * mrhoS0 + MLC ) * mGravity( m )
                                            ) * weight * DetJ * mScaleU;

            Help_R_1( i*mDimension + m ) += DNu_DX( i, m ) * (
//                                                 b * ( XL * p + ( 1 - XL ) * pC )
                                                b * ( p + ( 1 - XL ) * mSf * ( -t ) ) 
                                            ) * weight * DetJ * mScaleU;

            Help_R_1( i*mDimension + m ) += DNu_DX( i, m ) * (
                                                3 * malphaS * K * ( t - t0 )
                                            ) * weight * DetJ * mScaleU;

	      if ( GetValue( PLASTIC_FLAG ) == 0 ) // ELATIC 
	      {
		  for ( unsigned int k = 0; k < mDimension; k++ )
		      Help_R_1( i*mDimension + m ) -= DNu_DX( i, k ) * (
							  stressEff( k, m )
						      ) * weight * DetJ * mScaleU;
	      }
	      else
	      {
		  for ( unsigned int k = 0; k < 6; k++ )
		      Help_R_1( i*mDimension + m ) -= DNu_DX( k, i * mDimension + m ) * (
							  stressVectorEff( k )
						      ) * weight * DetJ * mScaleU;
	      }
        }
}

///-//////////////////////////////////////////////////////////////////////////////////////
/// R2 **************************
void Soil3Phase::AddInternalForcesToRHS2( Vector& Help_R_2, const Vector& No, const Matrix& DNo_DX, double weight, double DetJ, double t0, Vector u, double p, double t, Matrix Grad_u, Vector Grad_p, Vector Grad_t, Vector dot_u, double dot_p, double dot_t, double Div_dot_u )
{

// START: Method 1: condensed form : failed for hydral conduction !! wL = rhoL * vL
    double DMLC_De = GetDMLCDe( Grad_u, p, t );
    double DMLC_Dp = GetDMLCDp( Grad_u, p, t, t0 );
    double DMLC_Dt = GetDMLCDt( Grad_u, p, t, t0 );
    
    double rhoL= GetrhoL(p, t);
    double DrhoL_Dp= GetDrhoLDp( t );
    double DrhoL_Dt= GetDrhoLDt( p, t );

    Vector vL( mDimension );
    noalias( vL ) = GetvL( Grad_p, p, t );
    
//     Vector wL( mDimension );
//     noalias( wL ) = GetwL( Grad_p, t );

    for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
    {
        Help_R_2( i ) += No( i ) * (
                             DMLC_De / rhoL * Div_dot_u
                         ) * weight * DetJ * mScaleP;

        Help_R_2( i ) += No( i ) * (
                             DMLC_Dp / rhoL * dot_p
                         ) * weight * DetJ * mScaleP;

        Help_R_2( i ) += No( i ) * (
                             DMLC_Dt / rhoL * dot_t
                         ) * weight * DetJ * mScaleP;

        for ( unsigned int k = 0; k < mDimension; k++ )
        {
            // %%% R2---4
            Help_R_2( i ) += No( i )* (
                               1.0 / rhoL * ( DrhoL_Dp * Grad_p( k ) + DrhoL_Dt * Grad_t(k) )
                             ) *vL(k)
                             * weight * DetJ * mScaleP;
			     
            // %%% R2---5
            Help_R_2( i ) -= DNo_DX( i, k ) * (
                                 vL( k )
                             ) * weight * DetJ * mScaleP;
        }
    }
// // END

// START: Method 2: distributed form : divided by rhoL, GOOD!
//     double n= Getn( Grad_u, p, t, t0 );
//     double Dn_De= GetDnDe( Grad_u, t );
//     double Dn_Dp= GetDnDp(t);
//     double Dn_Dt= GetDnDt(Grad_u, p, t, t0);
//     double XL= GetXL(t);
//     double dXL_dt= GetdXLdt(t);
//     double rhoL= GetrhoL(p, t);
//     double rhoC= GetrhoC(p, t);
//     double DrhoL_Dp= GetDrhoLDp( t );
//     double DrhoL_Dt= GetDrhoLDt( p, t );
//     double DrhoC_Dp= GetDrhoCDp();
//     double DrhoC_Dt= GetDrhoCDt();

//     Vector vL( mDimension );
//     noalias( vL ) = GetvL( Grad_p, p, t );
// 
//     for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
//     {
//         // %%% R2---1
//         Help_R_2( i ) += No( i ) * (
//                              ( XL + rhoC / rhoL * (1.0-XL) ) * Dn_De
//                          ) * Div_dot_u
//                          * weight * DetJ * mScaleP;
//         // %%% R2---2
//         Help_R_2( i ) += No( i ) * (
//                              ( XL + rhoC / rhoL * (1.0-XL) ) * Dn_Dp //leads to vibration in pL profile, even problem in convergence with small ks, eg e+6.
//                              + n / rhoL * ( XL * DrhoL_Dp + (1.0-XL) * DrhoC_Dp )
//                          ) * dot_p
//                          * weight * DetJ * mScaleP;
//         // %%% R2---3
//         Help_R_2( i ) += No( i ) * (
//                              ( XL + rhoC / rhoL * (1.0-XL) ) * Dn_Dt
//                              + n / rhoL * ( XL * DrhoL_Dt + (1.0-XL) * DrhoC_Dt )
//                              + n * dXL_dt * (1.0-rhoC/rhoL)
//                          ) * dot_t
//                          * weight * DetJ * mScaleP;
// 
//         for ( unsigned int k = 0; k < mDimension; k++ )
//         {
//             // %%% R2---4
//             Help_R_2( i ) += No( i )* (
//                                1.0 / rhoL * ( DrhoL_Dp * Grad_p( k ) + DrhoL_Dt * Grad_t(k) )
//                              ) *vL(k)
//                              * weight * DetJ * mScaleP;
//             // %%% R2---5
//             Help_R_2( i ) -= DNo_DX( i, k ) * (
//                                  vL( k )
//                              ) * weight * DetJ * mScaleP;
//         }
//     }
// END

}

///-//////////////////////////////////////////////////////////////////////////////////////
/// R3 **************************
void Soil3Phase::AddInternalForcesToRHS3( Vector& Help_R_3, const Vector& No, const Matrix& DNo_DX, double weight, double DetJ, double t0, Vector u, double p, double t, Matrix Grad_u, Vector Grad_p, Vector Grad_t, Vector dot_u, double dot_p, double dot_t, double Div_dot_u )
{
    double rhoL= GetrhoL(p, t);
    
    double DSS_De = GetDSSDe( t );
    double DSS_Dp = GetDSSDp( t );
    double DSS_Dt = GetDSSDt( Grad_u, p, t );
    double DsL_Dp = GetDsLDp();
    double DsL_Dt = GetDsLDt();
    double ML = GetML( Grad_u, p, t, t0 );
    double PhiM = GetPhiM( Grad_p, p, t );

    double dsCL = 0.0;
    double DsC_Dp = 0.0;
    double DsC_Dt = 0.0;
    double MC = 0.0;
    double DMC_De = 0.0;
    double DMC_Dp = 0.0;
    double DMC_Dt = 0.0;
    if ( t < 0.0 ) //SAVE COMPUTING TIME
    {
        dsCL = GetdsCL( p, t );
        DsC_Dp = GetDsCDp();
        DsC_Dt = GetDsCDt();
        MC = GetMC( Grad_u, p, t, t0 );
        DMC_De = GetDMCDe( Grad_u, p, t );
        DMC_Dp = GetDMCDp( Grad_u, p, t, t0 );
        DMC_Dt = GetDMCDt( Grad_u, p, t, t0 );
    }
 
    Vector vL( mDimension );
    noalias( vL ) = GetvL( Grad_p, p, t );
    Vector q( mDimension );
    noalias( q ) = Getq( Grad_t, t );

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

///-//////////////////////////////////////////////////////////////////////////////////////
///  CONSTITUTIVE QUANTITES
///-//////////////////////////////////////////////////////////////////////////////////////
/// +++++++++++++++++++++++++++++++++++++++++++
/// ++++++++++++++ Group c1 +++++++++++++++++++
/// +++++++++++++++++++++++++++++++++++++++++++
/// solid grain pressure
// double Soil3Phase::GetpS( Matrix Grad_u, double t )
// {
//     double K = GetK( t );  
//     double strainVol = GetstrainVol( Grad_u );
// 
//     return K*strainVol;
// }
// /// grain pressure
// double Soil3Phase::GetpC( Matrix Grad_u, double t )
// {
//     double result= 0.0;
//     if ( t < 0.0 ) 
//     {
//       result = p + mSf*(-t);
//       double pS = GetpS(Grad_u, t);
//       
//       if(result >= pS);
//       result = pS;
//     }  
//     return result;
// }


/// XL: liquid saturation
double Soil3Phase::GetXL( double t )
{
    double result = 1.0;

    if ( t < 0.0 ) 
        result = pow( 1.0 + pow( -t / mtstar, 1.0 / ( 1.0 - mm ) ), -mm ); 

    return result;
}

double Soil3Phase::GetdXLdt( double t )
{
    double result = 0.0;

    if ( t < 0.0 ) 
        result = mm * pow( pow( -t / mtstar, 1.0 / ( 1.0 - mm ) ) + 1.0, -mm - 1.0 ) * pow( -t / mtstar, 1.0 / ( 1.0 - mm ) )
                 / (( mm - 1.0 ) * t ); 

    return result;
}

/// stress: Cauchy stress tensor
Matrix Soil3Phase::GetstressEff( Matrix Grad_u, double p, double t)
{
    Matrix result( mDimension, mDimension );
    noalias( result ) = ZeroMatrix( mDimension, mDimension );
    Matrix strain( mDimension, mDimension );
    noalias( strain ) = Getstrain( Grad_u );
    double K = GetK( t );
    double G = GetG( t );
    double strainVol = GetstrainVol( Grad_u );
    
    for ( unsigned int i = 0; i < mDimension; i++ )
        for ( unsigned int j = 0; j < mDimension; j++ )
            result( i, j ) = ( K - 2.0 * G / 3.0 ) * strainVol * KnoneckerDelta( i, j )
                             + 2.0 * G * strain( i, j );
    return result;
}

/// KnoneckerDelta
double Soil3Phase::KnoneckerDelta( int i, int j )
{
    if ( i == j )
        return 1.0;
    else
        return 0.0;
}

/// strain: linear strain tensor
Matrix Soil3Phase::Getstrain( Matrix Grad_u )
{
    Matrix result( mDimension, mDimension );
    noalias( result ) = ZeroMatrix( mDimension, mDimension );

    for ( unsigned int i = 0; i < mDimension; i++ )
        for ( unsigned int j = 0; j < mDimension; j++ )
            result( i, j ) = 0.5 * ( Grad_u( i, j ) + Grad_u( j, i ) );

    return result;
}

/// strainVol: volumetric strain
double Soil3Phase::GetstrainVol( Matrix Grad_u )
{
    double result = 0.0;
    Matrix strain( mDimension, mDimension );
    noalias( strain ) = Getstrain( Grad_u );

    for ( unsigned int i = 0; i < mDimension; i++ )
        result += strain( i, i );

    return result;
}

/// K: effective bulk modulus
double Soil3Phase::GetK( double t )
{
 
//     double e0 = mn0 / ( 1.0 - mn0 ); 	
//     return 1000.0 * ( 1.0 + e0 )  / mMaterialParameters[2]; 
    
    
//     double result = 0.0;
//     
//     //Lackener + zhou
//     double XL = GetXL( t );
//     double alpha = 3.0 * mkC / ( 3.0 * mkC + 4.0 * mgC );
//     double KSL = mkL * ( 1.0 + ( 1.0 - mn0) * ( 1.0 - mkL / mkS ) / ( mkL / mkS + alpha * mn0* ( 1.0 - mkL / mkS ) ) );
//     double KSC = mkC * ( 1.0 + ( 1.0 - mn0) * ( 1.0 - mkC / mkS ) / ( mkC / mkS + alpha * mn0* ( 1.0 - mkC / mkS ) ) );
// 
//     result = XL * KSL + ( 1.0 - XL ) * KSC;
//     return result;
    
    return 4.0*mgS*mkS*( 1 - mn0 ) / ( 3.0*mn0*mkS + 4.0*mgS );


//     double kG = mgL;
//     double gG = mgL;
//     double KGS = ( 3 * kG * mkS + 4 * gG * ( mkS + kG * mn0 - mkS * mn0 ) ) / ( 4 * gG + 3 * ( kG - kG * mn0 + mkS * mn0 ) );
//     double KCS = ( 3 * mkC * mkS + 4 * mgC * ( mkS + mkC * mn0 - mkS * mn0 ) ) / ( 4 * mgC + 3 * ( mkC - mkC * mn0 + mkS * mn0 ) );
//     
//     return XL * KGS + ( 1.0 - XL ) * KCS;
    
}

double Soil3Phase::GetdKdt( double t )
{
    double result = 0.0;

//     // Lackener + zhou
//     double dXL_dt = GetdXLdt( t );
//     double alpha = 3.0 * mkC / ( 3.0 * mkC + 4.0 * mgC );
//     double KSL = mkL * ( 1.0 + ( 1.0 - mn0) * ( 1.0 - mkL / mkS ) / ( mkL / mkS + alpha * mn0* ( 1.0 - mkL / mkS ) ) );
//     double KSC = mkC * ( 1.0 + ( 1.0 - mn0) * ( 1.0 - mkC / mkS ) / ( mkC / mkS + alpha * mn0* ( 1.0 - mkC / mkS ) ) );
// 
//     result = dXL_dt * ( KSL - KSC );
    return result;
    
//     double kG = mgL;
//     double gG = mgL;
//     double KGS = ( 3 * kG * mkS + 4 * gG * ( mkS + kG * mn0 - mkS * mn0 ) ) / ( 4 * gG + 3 * ( kG - kG * mn0 + mkS * mn0 ) );
//     double KCS = ( 3 * mkC * mkS + 4 * mgC * ( mkS + mkC * mn0 - mkS * mn0 ) ) / ( 4 * mgC + 3 * ( mkC - mkC * mn0 + mkS * mn0 ) );
//     
//     return dXL_dt * ( KGS - KCS );
    
}

/// G: effective shear modulus
double Soil3Phase::GetG( double t )
{
//     double nu = GetProperties()[POISSON_RATIO];	
//     double K = GetK(t);
//     return 1.5 * K * ( 1.0 - 2.0 * nu ) / ( 1.0 + nu );
    
//     double result = 0.0;
//     // Lackener + zhou
//     double XL = GetXL( t );
//     double beta = 1.2 * ( mkC + 2.0 * mgC ) / ( 3.0 * mkC + 4.0 * mgC );
//     double GSL = mgL * ( 1.0 + ( 1.0 - mn0) * ( 1.0 - mgL / mgS ) / ( mgL / mgS + beta * mn0* ( 1.0 - mgL / mgS ) ) );
//     double GSC = mgC * ( 1.0 + ( 1.0 - mn0) * ( 1.0 - mgC / mgS ) / ( mgC / mgS + beta * mn0* ( 1.0 - mgC / mgS ) ) );
// 
//     result = XL * GSL + ( 1.0 - XL ) * GSC;
//     return result;
    
    return  mgS*( 8.0*mgS + 9.0*mkS )*( 1 - mn0 ) / ( 6.0*mn0*( 2.0*mgS + mkS ) + 8.0*mgS + 9.0*mkS );

//     double kG = mgL;
//     double gG = mgL;
//     double GGS = -(( 5 * gG * ( 4 * gG + 3 * kG ) * mgS +   gG * ( 8 * gG + 9 * kG ) * ( gG - mgS ) * mn0 ) / ( -5 * gG * ( 4 * gG + 3 * kG ) +   6 * ( 2 * gG + kG ) * ( gG - mgS ) * mn0 ) );
//     double GCS = -(( 5 * mgC * mgS * ( 4 * mgC + 3 * mkC ) +   mgC * ( mgC - mgS ) * ( 8 * mgC + 9 * mkC ) * mn0 ) / ( -5 * mgC * ( 4 * mgC + 3 * mkC ) +   6 * ( mgC - mgS ) * ( 2 * mgC + mkC ) * mn0 ) );
//     
//     return  XL * GGS + ( 1.0 - XL ) * GCS;
}

/// b: Biot coeffecient
double Soil3Phase::Getb( double t )
{
//     double K = GetK( t ); 
//     return 1.0 - K / mkS; 
//     double XL = GetXL( t );
//     return 1.0 - K / (XL*mkS+(1-XL)*mKS);
   return 1.0;
}

double Soil3Phase::Getdbdt( double t )
{
//     double dK_dt = GetdKdt( t ); 
//     return -dK_dt / mkS;
    return 0.0;
}

/// +++++++++++++++++++++++++++++++++++++++++++
/// ++++++++++++++ Group c2 +++++++++++++++++++
/// +++++++++++++++++++++++++++++++++++++++++++


/// mass of ice and water
double Soil3Phase::GetMLC( Matrix Grad_u, double p, double t, double t0 )
{
    double result = 0.0;
    double rhoL = GetrhoL( p, t );
    double rhoC = GetrhoC( p, t );
    double n = Getn( Grad_u, p, t, t0 );
    double XL = GetXL( t );

    result = n * ( rhoL * XL + rhoC * ( 1 - XL ) );
    return result;
}

double Soil3Phase::GetDMLCDe( Matrix Grad_u, double p, double t )
{
    double result = 0.0;
    double rhoL = GetrhoL( p, t );
    double rhoC = GetrhoC( p, t );
    double Dn_De = GetDnDe( Grad_u, t );
    double XL = GetXL( t );

    result = Dn_De * ( rhoL * XL + rhoC * ( 1.0 - XL ) );
    return result;
}

double Soil3Phase::GetDMLCDp( Matrix Grad_u, double p, double t, double t0 )
{
    double result = 0.0;
    double rhoL = GetrhoL( p, t );
    double DrhoL_Dp = GetDrhoLDp( t );
    double rhoC = GetrhoC( p, t );
    double DrhoC_Dp = GetDrhoCDp();
    double n = Getn( Grad_u, p, t, t0 );
    double Dn_Dp = GetDnDp( t );
    double XL = GetXL( t );

    result = Dn_Dp * ( rhoL * XL + rhoC * ( 1 - XL ) ) + n * ( DrhoL_Dp * XL + DrhoC_Dp * ( 1 - XL ) );
    return result;
}

double Soil3Phase::GetDMLCDt( Matrix Grad_u, double p, double t, double t0 )
{
    double result = 0.0;
    double rhoL = GetrhoL( p, t );
    double DrhoL_Dt = GetDrhoLDt( p, t );
    double rhoC = GetrhoC( p, t );
    double DrhoC_Dt = GetDrhoCDt();
    double n = Getn( Grad_u, p, t, t0 );
    double Dn_Dt = GetDnDt( Grad_u, p, t, t0 );
    double XL = GetXL( t );
    double dXL_dt = GetdXLdt( t );

    result = Dn_Dt * ( rhoL * XL + rhoC * ( 1 - XL ) ) + n * ( DrhoL_Dt * XL + DrhoC_Dt * ( 1 - XL ) + ( rhoL - rhoC ) * dXL_dt );
    return result;
}


/// n: current porosity
double Soil3Phase::Getn( Matrix Grad_u, double p, double t, double t0 )
{
    double result = mn0;
//     double XL = GetXL( t );
//     double b = Getb( t );
    double strainVol = GetstrainVol( Grad_u );

//     result = mn0 + b * strainVol + ( b - mn0) / mkS * ( p + ( 1 - XL ) * mSf * ( -t ) ) - 3.0 * malphaS * ( b - mn0) * ( t - t0 );
//     result = mn0 + b * strainVol + ( b - mn0 * XL) / mkS * ( p ) - 3.0 * (XL* malphaS+ (1-XL)*malphaC) * ( b - mn0 * XL) * ( t - t0 );
    result = 1.0 - (1.0- mn0)*exp(-strainVol); 
    return result;
}

// double Soil3Phase::GetDnDe( double t )
double Soil3Phase::GetDnDe( Matrix Grad_u, double t )
{
//     return Getb( t );
    double strainVol = GetstrainVol( Grad_u ); 
    return  (1.0- mn0)*exp(-strainVol);
}

double Soil3Phase::GetDnDp( double t ) //WRONG!
{
//     double b = Getb( t ); 
//     return ( b - mn0) / mkS ;
// //     return ( b - mn0*XL) / mkS ;
    return  0.0;
}

double Soil3Phase::GetDnDt( Matrix Grad_u, double p, double t, double t0 )
{
    double result = 0.0;
//     double XL = GetXL( t );
//     double dXL_dt = GetdXLdt( t );
//     double b = Getb( t );
//     double db_dt = Getdbdt( t );
//     double strainVol = GetstrainVol( Grad_u );
// 
//     result = db_dt * strainVol
//              + db_dt / mkS * ( p + ( 1 - XL ) * mSf * ( -t ) )
//              - ( b - mn0) / mkS * mSf * ( dXL_dt * ( -t ) + ( 1 - XL ) )
//              - 3.0 * malphaS * ( db_dt * ( t - t0 )  + ( b - mn0) );
// //     result = db_dt * strainVol
// //              + (db_dt - mn0*dXL_dt)/ mkS * ( p ) 
// //              - 3.0 * dXL_dt* (malphaS - malphaC) * ( b - mn0 * XL) * ( t - t0 )
// //              - 3.0 * (XL* malphaS+ (1-XL)*malphaC) * ( (db_dt-mn0*dXL_dt) * ( t - t0 )  + ( b - mn0*XL) );
    return result;
}


/// rhoL: current water density
double Soil3Phase::GetrhoL( double p, double t )
{ 
    return mrhoL0 * ( 1.0 - 3.0 * malphaL * ( t ) + p / mkL );
     
//     double M_w = /*GetProperties()[MOLAR_MASS_WATER];*/ 24.042835;
//     double R = /*GetProperties()[GAS_CONSTANT];*/  8.314472 * pow( 1000, 2.0 );
//   
//     return mrhoL0 + p * M_w / ( R * ( 273.15 + t ) ) ; 
    
}

double Soil3Phase::GetDrhoLDp( double t )
{
    return mrhoL0 / mkL;
    
//     double M_w = /*GetProperties()[MOLAR_MASS_WATER];*/ 24.042835;
//     double R = /*GetProperties()[GAS_CONSTANT];*/  8.314472 * pow( 1000, 2.0 );
//     return M_w / ( R * ( 273.15 + t ) );
}

double Soil3Phase::GetDrhoLDt( double p, double t )
{  
    return -3*malphaL*mrhoL0;
    
//     double M_w = /*GetProperties()[MOLAR_MASS_WATER];*/ 24.042835;
//     double R = /*GetProperties()[GAS_CONSTANT];*/  8.314472 * pow( 1000, 2.0 );
// 
//     return (-1)*M_w/R*pow(273.15+t,-2)*p;
}
/// rhoC: current ice density
double Soil3Phase::GetrhoC( double p, double t )
{ 
    double pC=  0.0;
    if(t < 0.0)
      pC= mSf * ( -t );
//     double strainVol = GetstrainVol( Grad_u );
//     double K = GetK(t);
//     double pEff = K*strainVol;
//     if( pC > = pEff)
//     {
// 	pC= pEff;
// 	std::cout<<"pC >= pEff: RESET pC = pEff = "<<pC<<std::endl;
//     }
    
    return mrhoC0 * ( 1.0 - 3.0 * malphaC * t + ( pC ) / mkC ); //approximation
}

double Soil3Phase::GetDrhoCDp()
{  
    return mrhoC0 / mkC;
}

double Soil3Phase::GetDrhoCDt()
{ 
    return -mrhoC0*( 3.0*malphaC + mSf / mkC );
}
 
/// ML: water mass
double Soil3Phase::GetML( Matrix Grad_u, double p, double t, double t0 )
{
    double result = 0.0;
    double rhoL = GetrhoL( p, t );
    double n = Getn( Grad_u, p, t, t0 );
    double XL = GetXL( t );

    result = rhoL * n * XL;
    return result;
}

double Soil3Phase::GetDMLDe( Matrix Grad_u, double p, double t )
{
    double result = 0.0;
    double rhoL = GetrhoL( p, t );
    double Dn_De = GetDnDe( Grad_u, t );
    double XL = GetXL( t );

    result = rhoL * Dn_De * XL;
    return result;
}

double Soil3Phase::GetDMLDp( Matrix Grad_u, double p, double t, double t0 )
{
    double result = 0.0;
    double rhoL = GetrhoL( p, t );
    double DrhoL_Dp = GetDrhoLDp( t );
    double n = Getn( Grad_u, p, t, t0 );
    double Dn_Dp = GetDnDp( t );
    double XL = GetXL( t );

    result = ( DrhoL_Dp * n + rhoL * Dn_Dp ) * XL;
    return result;
}

double Soil3Phase::GetDMLDt( Matrix Grad_u, double p, double t, double t0 )
{
    double result = 0.0;
    double rhoL = GetrhoL( p, t );
    double DrhoL_Dt = GetDrhoLDt( p, t );
    double n = Getn( Grad_u, p, t, t0 );
    double Dn_Dt = GetDnDt( Grad_u, p, t, t0 );
    double XL = GetXL( t );
    double dXL_dt = GetdXLdt( t );

    result = DrhoL_Dt * n * XL + rhoL * ( Dn_Dt * XL + n * dXL_dt );
    return result;
}

/// MC: ice mass
double Soil3Phase::GetMC( Matrix Grad_u, double p, double t, double t0 )
{
    double result = 0.0;
    double rhoC = GetrhoC( p, t );
    double n = Getn( Grad_u, p, t, t0 );
    double XL = GetXL( t );

    result = rhoC * n * ( 1.0 - XL );
    return result;
}

double Soil3Phase::GetDMCDe( Matrix Grad_u, double p, double t )
{
    double result = 0.0;
    double rhoC = GetrhoC( p, t );
    double Dn_De = GetDnDe( Grad_u, t );
    double XL = GetXL( t );

    result = rhoC * Dn_De * ( 1.0 - XL );
    return result;
}

double Soil3Phase::GetDMCDp( Matrix Grad_u, double p, double t, double t0 )
{
    double result = 0.0;
    double rhoC = GetrhoC( p, t );
    double DrhoC_Dp = GetDrhoCDp();
    double n = Getn( Grad_u, p, t, t0 );
    double Dn_Dp = GetDnDp( t );
    double XL = GetXL( t );

    result = ( DrhoC_Dp * n + rhoC * Dn_Dp ) * ( 1.0 - XL );
    return result;
}

double Soil3Phase::GetDMCDt( Matrix Grad_u, double p, double t, double t0 )
{
    double result = 0.0;
    double rhoC = GetrhoC( p, t );
    double DrhoC_Dt = GetDrhoCDt();
    double n = Getn( Grad_u, p, t, t0 );
    double Dn_Dt = GetDnDt( Grad_u, p, t, t0 );
    double XL = GetXL( t );
    double dXL_dt = GetdXLdt( t );

    result = DrhoC_Dt * n * ( 1.0 - XL ) + rhoC * ( Dn_Dt * ( 1.0 - XL ) - n * dXL_dt );
    return result;
}

/// wL: liquid water flow
// Vector Soil3Phase::GetwL( Vector Grad_p, double t )
// {
//     Vector result( mDimension );
//     noalias( result ) = ZeroVector( mDimension ); 
//     double XL = GetXL( t );
//     double lambda = 1.0;
//     double delta = 1.0; 
// 
//     if ( t < 0.0 )
//     { 
//         lambda = sqrt( XL ) * pow( 1.0 - pow( 1.0 - pow( XL, 1.0 / mm ), mm ), 2.0 ); // Luckner et al.,1989
//         double a = 1.5963e-2;
//         double b = 509.53;
//         double c = 123.15;
//         delta = a * exp( b / ( c + t ) );
//     }
// 
//     for ( unsigned int k = 0; k < mDimension; k++ )
//         result( k ) = mrhoL0 * mkappa0 * lambda / ( metaL * delta ) * ( -Grad_p( k ) + mrhoL0 * mGravity( k ) );
// 
//     return result;
// }

/// vL: darcy water flow
Vector Soil3Phase::GetvL( Vector Grad_p, double p, double t )
{
    Vector result( mDimension );
    noalias( result ) = ZeroVector( mDimension ); 
    double lambda = 1.0; 
    double rhoL = GetrhoL( p, t );
    double g = 9.81 * mUnitRatio;

    /// too steep for small mm --> convergence problem
//     if ( t < 0.0 )
//     { 
// 	double XL = GetXL( t );
//         lambda = sqrt( XL ) * pow( 1.0 - pow( 1.0 - pow( XL, 1.0 / mm ), mm ), 2.0 ); // Luckner et al.,1989 
//     } 
    /// approximation --> good in convergence 
    if ( t < 0.0 )
    { 
	mtstar = 0.5*mtstar;
        lambda =  GetXL( t );  
	mtstar = 2.0*mtstar;
    }

    for ( unsigned int k = 0; k < mDimension; k++ ) 
	  result[k] =  mkappa0 * lambda / ( rhoL * g ) * ( -Grad_p( k ) + rhoL * mGravity( k ) ); 

    return result;
}

/// PhiM: mechanical dissipation
double Soil3Phase::GetPhiM( Vector Grad_p, double p, double t )
{
    double result = 0.0;
    
// //     double XL = GetXL( t );
// //     double lambda = 1.0;
// //     double delta = 1.0; 
// //     
// //     if ( t < 0.0 )
// //     {
// //         lambda = sqrt( XL ) * pow( 1.0 - pow( 1.0 - pow( XL, 1.0 / mm ), mm ), 2.0 ); // Luckner et al.,1989
// // 
// //         double a = 1.5963e-2;
// //         double b = 509.53;
// //         double c = 123.15;
// //         delta = a * exp( b / ( c + t ) );
// //     }
// // 
// //     for ( unsigned int k = 0; k < mDimension; k++ )
// //         result += ( -Grad_p( k ) + mrhoL0 * mGravity( k ) ) * ( -Grad_p( k ) + mrhoL0 * mGravity( k ) );
// // 
// //     result = mkappa0 * lambda / ( metaL * delta ) * result;
//     
//     double rhoL = GetrhoL( p, t );
//     Vector vL( mDimension );
//     noalias( vL ) = GetvL( Grad_p, p, t );
// 
//     for ( unsigned int k = 0; k < mDimension; k++ )
// 	result += vL( k ) * ( -Grad_p( k ) + rhoL * mGravity( k ) );

    return result;
}
/// +++++++++++++++++++++++++++++++++++++++++++
/// ++++++++++++++ Group c3 +++++++++++++++++++
/// +++++++++++++++++++++++++++++++++++++++++++


/// SS: entropy of skeleton
double Soil3Phase::GetDSSDe( double t )
{
    double K = GetK( t );
    return 3.0 * malphaS * K;
//     double XL = GetXL( t );
//     return 3.0 * (XL*malphaS+(1-XL)*malphaC) * K;
}

double Soil3Phase::GetDSSDp( double t )
{ 
    double b = Getb( t ); 
    return - 3.0 * malphaS * ( b - mn0);
//     return - 3.0 * (XL*malphaS+(1-XL)*malphaC) * ( b - mn0) ;
}

double Soil3Phase::GetDSSDt( Matrix Grad_u, double p, double t )
{
    double result = 0.0; 
    double XL = GetXL( t );
    double dXL_dt = GetdXLdt( t );
    double dK_dt = GetdKdt( t );
    double b = Getb( t );
    double db_dt = Getdbdt( t );
    double strainVol = GetstrainVol( Grad_u );

    result = 3.0 * malphaS * dK_dt * strainVol
             - 3.0 * malphaS * db_dt * ( p + ( 1 - XL ) * mSf * ( -t ) )
             + 3.0 * malphaS * ( b - mn0) * mSf * ( dXL_dt * ( -t ) + ( 1 - XL ) )
             + mrhoS0 * ( 1.0 - mn0) * mcS / mTf;

    return result;
}

/// sL: specific water entropy
double Soil3Phase::GetDsLDp()
{  
    return - 3.0 * malphaL / mrhoL0;
}

double Soil3Phase::GetDsLDt()
{ 
    return mcL / mTf;
}
 
/// sC: specific entropy
double Soil3Phase::GetDsCDp()
{ 
    return -3.0 * malphaC / mrhoC0;
}

double Soil3Phase::GetDsCDt()
{ 
    return 3.0 * malphaC * mSf / mrhoC0 + mcC / mTf;
}

/// ds = sC-sL : entropy differenece
double Soil3Phase::GetdsCL( double p, double t )
{
    return ( mcC - mcL ) * t / mTf + 3.0 * ( malphaL / mrhoL0 - malphaC / mrhoC0 ) * p + 3.0 * malphaC * mSf * t / mrhoC0;
}
double Soil3Phase::GetDdsCLDp( )
{
    return 3.0*( malphaL / mrhoL0 - malphaC / mrhoC0 );
}
double Soil3Phase::GetDdsCLDt()
{ 
    return ( mcC - mcL ) / mTf + 3.0 * malphaC * mSf / mrhoC0;
}


/// +++++++++++++++++++++++++++++++++++++++++++
/// q: heat flux

Vector Soil3Phase::Getq( Vector Grad_t, double t )
{
    Vector result( mDimension );
    noalias( result ) = ZeroVector( mDimension );
    double XL = GetXL( t );
    double lambdaTot = pow( mlambdaS, 1.0 - mn0) * pow( mlambdaL, mn0* XL ) * pow( mlambdaC, mn0* ( 1.0 - XL ) ); //Zhu

    for ( unsigned int k = 0; k < mDimension; k++ )
        result( k ) = -lambdaTot * Grad_t( k );

    return result;

}


/// +++++++++++++++++++++++++++++++++++++++++++
/// ++++++++++ additional quantities ++++++++++
/// +++++++++++++++++++++++++++++++++++++++++++
/// +++++++++++++++++++++++++++++++++++++++++++
/// +++++++++++++++++++++++++++++++++++++++++++
/// - WriteNodalResults
void Soil3Phase::FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo )
{
  for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
      mConstitutiveLawVector[i]->FinalizeSolutionStep( GetProperties(),
	      GetGeometry(),
	      row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ),
	      CurrentProcessInfo );
	      
    Matrix DNo_DX( mNodesNumberOther, mDimension );
    Matrix DNo_DX_e( mNodesNumberOther, mDimension );
    Matrix nodesLocalCoords( mNodesNumberOther, mDimension );

    nodesLocalCoords( 0, 0 ) = -1.0;
    nodesLocalCoords( 0, 1 ) = -1.0;
    nodesLocalCoords( 0, 2 ) = -1.0;
    nodesLocalCoords( 1, 0 ) = 1.0;
    nodesLocalCoords( 1, 1 ) = -1.0;
    nodesLocalCoords( 1, 2 ) = -1.0;
    nodesLocalCoords( 2, 0 ) = 1.0;
    nodesLocalCoords( 2, 1 ) = 1.0;
    nodesLocalCoords( 2, 2 ) = -1.0;
    nodesLocalCoords( 3, 0 ) = -1.0;
    nodesLocalCoords( 3, 1 ) = 1.0;
    nodesLocalCoords( 3, 2 ) = -1.0;
    nodesLocalCoords( 4, 0 ) = -1.0;
    nodesLocalCoords( 4, 1 ) = -1.0;
    nodesLocalCoords( 4, 2 ) = 1.0;
    nodesLocalCoords( 5, 0 ) = 1.0;
    nodesLocalCoords( 5, 1 ) = -1.0;
    nodesLocalCoords( 5, 2 ) = 1.0;
    nodesLocalCoords( 6, 0 ) = 1.0;
    nodesLocalCoords( 6, 1 ) = 1.0;
    nodesLocalCoords( 6, 2 ) = 1.0;
    nodesLocalCoords( 7, 0 ) = -1.0;
    nodesLocalCoords( 7, 1 ) = 1.0;
    nodesLocalCoords( 7, 2 ) = 1.0;
    
    for ( unsigned int i = ( mNodesOtherMin - 1 ) ; i < mNodesOtherMax ;i++ )
    {
        Vector local_coords( 3 );
        local_coords( 0 ) = nodesLocalCoords( i, 0 );
        local_coords( 1 ) = nodesLocalCoords( i, 1 );
        local_coords( 2 ) = nodesLocalCoords( i, 2 );

        noalias( DNo_DX_e ) = GetGeometry().ShapeFunctionsLocalGradients( DNo_DX_e, local_coords );
        noalias( DNo_DX ) = prod( DNo_DX_e, mInvJ0[0] );

        Vector Grad_p( mDimension );
        noalias( Grad_p ) = ZeroVector( mDimension );
        Vector Grad_t( mDimension );
        noalias( Grad_t ) = ZeroVector( mDimension );
        Matrix Grad_u( mDimension, mDimension );
        noalias( Grad_u ) = ZeroMatrix( mDimension, mDimension );

        for ( unsigned int j = mNodesOtherMin - 1 ; j < mNodesOtherMax ;j++ )
            for ( unsigned int k = 0; k < mDimension; k++ )
            {
                Grad_p( k ) += GetGeometry()[j].GetSolutionStepValue( WATER_PRESSURE ) * DNo_DX( j, k );
                Grad_t( k ) += GetGeometry()[j].GetSolutionStepValue( TEMPERATURE ) * DNo_DX( j, k );
                Grad_u( k, 0 ) += GetGeometry()[j].GetSolutionStepValue( DISPLACEMENT_X ) * DNo_DX( j, k ); //Only valid for 8N
                Grad_u( k, 1 ) += GetGeometry()[j].GetSolutionStepValue( DISPLACEMENT_Y ) * DNo_DX( j, k );
                Grad_u( k, 2 ) += GetGeometry()[j].GetSolutionStepValue( DISPLACEMENT_Z ) * DNo_DX( j, k );
            }

        double p = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE );

        double t = GetGeometry()[i].GetSolutionStepValue( TEMPERATURE ); 
    double pC = 0.0;
//     double pS = -(GetstressEff( Grad_u, p, t )(0,0)+GetstressEff( Grad_u, p, t )(1,1)+GetstressEff( Grad_u, p, t )(2,2) )/3.0; 
    double pS =  -GetK( t )*GetstrainVol( Grad_u );
    if(t<0.0)
    {
      pC = p + mSf*(-t);
//       if(pC>pS)
// 	pC = pS;
    }
    
        GetGeometry()[i].GetSolutionStepValue( ICE_SATURATION ) = 1.0 - GetXL( t );
        GetGeometry()[i].GetSolutionStepValue( WATER_DENSITY ) = GetrhoL( p, t ) * pow(mUnitRatio, 3.0); // kg/m^3
        GetGeometry()[i].GetSolutionStepValue( ICE_DENSITY ) = GetrhoC( p, t ) * pow(mUnitRatio, 3.0); // kg/m^3
//         GetGeometry()[i].GetSolutionStepValue( ICE_PRESSURE ) = p + mSf*(-t);
        GetGeometry()[i].GetSolutionStepValue( ICE_PRESSURE ) = pC;
        GetGeometry()[i].GetSolutionStepValue( MECH_DISSIPATION ) = GetPhiM( Grad_p, p, t );

        GetGeometry()[i].GetSolutionStepValue( LINEAR_STRAIN ) = GetstrainVol( Grad_u );
        GetGeometry()[i].GetSolutionStepValue( POROSITY ) = Getn( Grad_u, p, t, mT0e[i] );
        GetGeometry()[i].GetSolutionStepValue( ICE_VOLUME_FRACTION ) = Getn( Grad_u, p, t, mT0e[i] ) * ( 1.0 - GetXL( t ) );
        GetGeometry()[i].GetSolutionStepValue( WATER_MASS ) = GetML( Grad_u, p, t, mT0e[i] );
        GetGeometry()[i].GetSolutionStepValue( ICE_MASS ) = GetMC( Grad_u, p, t, mT0e[i] );

        GetGeometry()[i].GetSolutionStepValue( EFFECTIVE_STRESS ) = -pS;
        GetGeometry()[i].GetSolutionStepValue( TOTAL_STRESS ) = -pS-Getb( t )*(GetXL( t )*p+(1-GetXL( t ))*pC);
        GetGeometry()[i].GetSolutionStepValue( FRICTION_COEFFICIENT ) = Getb( t );   // pSkeleton

        for ( unsigned int m = 0; m < mDimension; m++ )
        {
            GetGeometry()[i].GetSolutionStepValue( WATER_FLOW )( m ) = 24.0*3600*GetvL( Grad_p, p, t )( m ) / mUnitRatio;  // m/day
            GetGeometry()[i].GetSolutionStepValue( HEAT_FLOW )( m ) = Getq( Grad_t, t )( m );
        }
    }

    if ( GetGeometry().size() == 20 )
    {
        Interpolate( WATER_PRESSURE, CurrentProcessInfo );
        Interpolate( TEMPERATURE, CurrentProcessInfo );

        Interpolate( ICE_SATURATION, CurrentProcessInfo );
        Interpolate( WATER_DENSITY, CurrentProcessInfo );
        Interpolate( ICE_DENSITY, CurrentProcessInfo );
        Interpolate( ICE_PRESSURE, CurrentProcessInfo );
//         Interpolate( MECH_DISSIPATION, CurrentProcessInfo );

// 		Interpolate( LINEAR_STRAIN, CurrentProcessInfo );
// 		Interpolate( POROSITY, CurrentProcessInfo );
// 		Interpolate( ICE_VOLUME_FRACTION, CurrentProcessInfo );
// 		Interpolate( WATER_MASS, CurrentProcessInfo );
// 		Interpolate( ICE_MASS, CurrentProcessInfo );

        
    }
}

void Soil3Phase::Interpolate( const Variable<double>& rVariable, const ProcessInfo& rCurrentProcessInfo )
{
    GetGeometry()[8].GetSolutionStepValue( rVariable ) = ( GetGeometry()[0].GetSolutionStepValue( rVariable ) + GetGeometry()[1].GetSolutionStepValue( rVariable ) ) / 2.0;
    GetGeometry()[9].GetSolutionStepValue( rVariable ) = ( GetGeometry()[1].GetSolutionStepValue( rVariable ) + GetGeometry()[2].GetSolutionStepValue( rVariable ) ) / 2.0;
    GetGeometry()[10].GetSolutionStepValue( rVariable ) = ( GetGeometry()[2].GetSolutionStepValue( rVariable ) + GetGeometry()[3].GetSolutionStepValue( rVariable ) ) / 2.0;
    GetGeometry()[11].GetSolutionStepValue( rVariable ) = ( GetGeometry()[3].GetSolutionStepValue( rVariable ) + GetGeometry()[0].GetSolutionStepValue( rVariable ) ) / 2.0;
    GetGeometry()[12].GetSolutionStepValue( rVariable ) = ( GetGeometry()[4].GetSolutionStepValue( rVariable ) + GetGeometry()[5].GetSolutionStepValue( rVariable ) ) / 2.0;
    GetGeometry()[13].GetSolutionStepValue( rVariable ) = ( GetGeometry()[5].GetSolutionStepValue( rVariable ) + GetGeometry()[6].GetSolutionStepValue( rVariable ) ) / 2.0;
    GetGeometry()[14].GetSolutionStepValue( rVariable ) = ( GetGeometry()[6].GetSolutionStepValue( rVariable ) + GetGeometry()[7].GetSolutionStepValue( rVariable ) ) / 2.0;
    GetGeometry()[15].GetSolutionStepValue( rVariable ) = ( GetGeometry()[7].GetSolutionStepValue( rVariable ) + GetGeometry()[4].GetSolutionStepValue( rVariable ) ) / 2.0;
    GetGeometry()[16].GetSolutionStepValue( rVariable ) = ( GetGeometry()[0].GetSolutionStepValue( rVariable ) + GetGeometry()[4].GetSolutionStepValue( rVariable ) ) / 2.0;
    GetGeometry()[17].GetSolutionStepValue( rVariable ) = ( GetGeometry()[1].GetSolutionStepValue( rVariable ) + GetGeometry()[5].GetSolutionStepValue( rVariable ) ) / 2.0;
    GetGeometry()[18].GetSolutionStepValue( rVariable ) = ( GetGeometry()[2].GetSolutionStepValue( rVariable ) + GetGeometry()[6].GetSolutionStepValue( rVariable ) ) / 2.0;
    GetGeometry()[19].GetSolutionStepValue( rVariable ) = ( GetGeometry()[3].GetSolutionStepValue( rVariable ) + GetGeometry()[7].GetSolutionStepValue( rVariable ) ) / 2.0;
}

///-////////////////////////////////////
/////
///- //////////////////////////////////////////////////////
void Soil3Phase::GetValueOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo ) 
{
    KRATOS_TRY
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );
    if(rValues.size() != integration_points.size())
	    rValues.resize( integration_points.size() ); 
    
   for ( unsigned int GaussPoint = 0; GaussPoint < integration_points.size(); GaussPoint++ )
    { 
        rValues[GaussPoint] = mConstitutiveLawVector[GaussPoint]->GetValue( rVariable, rValues[GaussPoint] );
    }

    KRATOS_CATCH( "" )

}

void Soil3Phase::GetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo ) 
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

void Soil3Phase::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo ) {}

void Soil3Phase::SetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo ) {}

void Soil3Phase::SetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo ) {}

void Soil3Phase::CalculateOnIntegrationPoints( const Variable<double>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo ) {}

int Soil3Phase::Check( const Kratos::ProcessInfo& rCurrentProcessInfo )
{
    return 0;
}

/// strain: linear strain tensor
Vector Soil3Phase::GetstrainVector( Matrix Grad_u ) //new
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
}
