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


// System includes
#include <cmath>


// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/unfrozen_soil.h"
#include "utilities/math_utils.h"
//#include "includes/constitutive_law.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/tetrahedra_3d_4.h"
#include "freezing_soil_application.h"
#include "freezing_soil.h"
#include "structural_application/custom_utilities/sd_math_utils.h"

namespace Kratos
{

UnfrozenSoil::UnfrozenSoil ( IndexType NewId, GeometryType::Pointer pGeometry )
        : Element ( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}


UnfrozenSoil::UnfrozenSoil ( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
        : Element ( NewId, pGeometry, pProperties )
{  
    mNodesDispMin = 1; // Disp: displacement
    mNodesOtherMin = 1; // Other: ice volume fraction, water pressure, or temperature
        
    if ( GetGeometry().size() == 27 || GetGeometry().size() == 20 || GetGeometry().size() == 8 )
    {
        mNodesOtherMax = 8;
        mThisGeometryOther = Geometry< Node<3> >::Pointer ( new Hexahedra3D8 <Node<3> > (
                                 GetGeometry() ( 0 ), GetGeometry() ( 1 ), GetGeometry() ( 2 ), GetGeometry() ( 3 ),
                                 GetGeometry() ( 4 ), GetGeometry() ( 5 ), GetGeometry() ( 6 ), GetGeometry() ( 7 ) ) ); 
        mThisIntegrationMethod = GeometryData::GI_GAUSS_3;

        if ( GetGeometry().size() == 8 )
            mNodesDispMax = 8;
        else if ( GetGeometry().size() == 20 )
            mNodesDispMax = 20;
        else
            mNodesDispMax = 27;
    }
    else if ( GetGeometry().size() == 10 )
    {
        mNodesDispMax = 10;
        mNodesOtherMax = 4;
        mThisGeometryOther = Geometry< Node<3> >::Pointer ( new Tetrahedra3D4 <Node<3> > (
                                 GetGeometry() ( 0 ), GetGeometry() ( 1 ), GetGeometry() ( 2 ), GetGeometry() ( 3 ) ) );
        mThisIntegrationMethod = GeometryData::GI_GAUSS_5;
    }
    else
        KRATOS_THROW_ERROR ( std::logic_error, "This element matches only with a linear hexahedra (8) or quadratic hexahedra (20 or 27) geometry" , *this );
}

Element::Pointer UnfrozenSoil::Create ( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer ( new UnfrozenSoil ( NewId, GetGeometry().Create ( ThisNodes ), pProperties ) );
}

UnfrozenSoil::~UnfrozenSoil()
{
}


void UnfrozenSoil::ResetConstitutiveLaw()
{
    KRATOS_TRY
    for( unsigned int i=0; i<mConstitutiveLawVector.size(); i++ )
        mConstitutiveLawVector[i]->ResetMaterial(GetProperties(), GetGeometry(), row(GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod), i) );
    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

// Degree of freedom: displacment: u[m], water pressure p[Pa]
void UnfrozenSoil::Initialize()
{
    KRATOS_TRY

    mNodesNumberDisp = mNodesDispMax - mNodesDispMin + 1;
    mNodesNumberOther = mNodesOtherMax - mNodesOtherMin + 1;
    mDimension = GetGeometry().WorkingSpaceDimension();
    mMatSizeU = mNodesNumberDisp * mDimension;
    mMatSizeO = mNodesNumberOther * 1;
    mNumberU = 1;
    mNumberO = 1;
    mMatSize = mMatSizeU * mNumberU + mMatSizeO * mNumberO;
    mAddIndexU = mMatSizeU * mNumberU;
    mScaleU = GetProperties()[SCALE_U];
    mScaleP = GetProperties()[SCALE_O];

    mGravity = GetProperties()[GRAVITY];
    mn0 = GetProperties()[POROSITY];
//     mtstar = GetProperties()[FIRST_SATURATION_PARAM];
//     mm = GetProperties()[SECOND_SATURATION_PARAM];
    mkappa0 = GetProperties()[PERMEABILITY_WATER];
    mMaterialParameters = GetProperties()[MATERIAL_PARAMETERS];
    int addN = 7;
//     mrhoS0 = mMaterialParameters[addN+0];
    mrhoL0 = mMaterialParameters[addN+1];
//     mrhoC0 = mMaterialParameters[addN+2];
    mkS = mMaterialParameters[addN+3];
    mkL = mMaterialParameters[addN+4];
//     mkC = mMaterialParameters[addN+5];
//     mgS = mMaterialParameters[addN+6];
//     mgL = mMaterialParameters[addN+7];
//     mgC = mMaterialParameters[addN+8];
//     mcS = mMaterialParameters[addN+9];
//     mcL = mMaterialParameters[addN+10];
//     mcC = mMaterialParameters[addN+11];
//     mlambdaS = mMaterialParameters[addN+12];
//     mlambdaL = mMaterialParameters[addN+13];
//     mlambdaC = mMaterialParameters[addN+14];
//     malphaS = mMaterialParameters[addN+15];
//     malphaL = mMaterialParameters[addN+16];
//     malphaC = mMaterialParameters[addN+17];
//     mTf = mMaterialParameters[addN+18];
//     mSf = mMaterialParameters[addN+19];
    metaL = mMaterialParameters[addN+20];
    
if ( GetValue(KRATOS_WATCH_FLAG)==1 ) 
  KRATOS_WATCH(mMaterialParameters);

    // Needed for assigning nonzero Dirchlet bc: eg. Ti=5 or pLi= 1000;
    for ( unsigned int i = ( mNodesDispMin - 1 ) ; i < mNodesDispMax ; i++ )
    {
        GetGeometry()[i].GetSolutionStepValue ( DISPLACEMENT_NULL ) = GetGeometry()[i].GetSolutionStepValue ( DISPLACEMENT );
        GetGeometry()[i].GetSolutionStepValue ( DISPLACEMENT_EINS ) = GetGeometry()[i].GetSolutionStepValue ( DISPLACEMENT );
        GetGeometry()[i].GetSolutionStepValue ( DISPLACEMENT_DT ) = ZeroVector ( 3 );
        GetGeometry()[i].GetSolutionStepValue ( DISPLACEMENT_NULL_DT ) = ZeroVector ( 3 );
        GetGeometry()[i].GetSolutionStepValue ( DISPLACEMENT_EINS_DT ) = ZeroVector ( 3 );
        GetGeometry()[i].GetSolutionStepValue ( DISPLACEMENT_OLD ) = ZeroVector ( 3 );
        GetGeometry()[i].GetSolutionStepValue ( ACCELERATION ) = ZeroVector ( 3 );
        GetGeometry()[i].GetSolutionStepValue ( ACCELERATION_NULL ) = ZeroVector ( 3 );
        GetGeometry()[i].GetSolutionStepValue ( ACCELERATION_EINS ) = ZeroVector ( 3 );
    }

    for ( unsigned int i = ( mNodesOtherMin - 1 ) ; i < mNodesOtherMax ; i++ )
    {
        GetGeometry()[i].GetSolutionStepValue ( WATER_PRESSURE_NULL ) = GetGeometry()[i].GetSolutionStepValue ( WATER_PRESSURE );
        GetGeometry()[i].GetSolutionStepValue ( WATER_PRESSURE_EINS ) = GetGeometry()[i].GetSolutionStepValue ( WATER_PRESSURE );
        GetGeometry()[i].GetSolutionStepValue ( WATER_PRESSURE_DT ) = 0;
        GetGeometry()[i].GetSolutionStepValue ( WATER_PRESSURE_NULL_DT ) = 0;
        GetGeometry()[i].GetSolutionStepValue ( WATER_PRESSURE_EINS_DT ) = 0;
        GetGeometry()[i].GetSolutionStepValue ( WATER_PRESSURE_ACCELERATION ) = 0;
        GetGeometry()[i].GetSolutionStepValue ( WATER_PRESSURE_NULL_ACCELERATION ) = 0;
        GetGeometry()[i].GetSolutionStepValue ( WATER_PRESSURE_EINS_ACCELERATION ) = 0;
 
    }

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints ( mThisIntegrationMethod );
    mInvJ0.resize ( integration_points.size() );
    mDetJ0.resize ( integration_points.size(), false );
    noalias ( mDetJ0 ) = ZeroVector ( integration_points.size() );
    GeometryType::JacobiansType J0 ( integration_points.size() );
    J0 = GetGeometry().Jacobian ( J0, mThisIntegrationMethod );
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        MathUtils<double>::InvertMatrix ( J0[PointNumber], mInvJ0[PointNumber], mDetJ0[PointNumber] );
 
    // initialize material
    if ( mConstitutiveLawVector.size() != integration_points.size() )
        mConstitutiveLawVector.resize ( integration_points.size() );
    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
    {
        mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
        mConstitutiveLawVector[i]->InitializeMaterial ( GetProperties(), GetGeometry(), row ( GetGeometry().ShapeFunctionsValues ( mThisIntegrationMethod ), i ) );
    } 
    KRATOS_CATCH ( "" )
} 

void UnfrozenSoil::CalculateAll ( MatrixType& rLeftHandSideMatrix,
                                  VectorType& rRightHandSideVector, ProcessInfo& rProcessInfo,
                                  bool CalculateStiffnessMatrixFlag, bool CalculateResidualVectorFlag )
{
    KRATOS_TRY
    //resizing as needed the RHS

    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size() != mMatSize )
            rRightHandSideVector.resize ( mMatSize );
        noalias ( rRightHandSideVector ) = ZeroVector ( mMatSize ); //resetting RHS
    }

    //resizing as needed the LHS
    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != mMatSize )
            rLeftHandSideMatrix.resize ( mMatSize, mMatSize );
        noalias ( rLeftHandSideMatrix ) = ZeroMatrix ( mMatSize, mMatSize ); //resetting LHS
    }

    //reading integration points and local Gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints ( mThisIntegrationMethod );
    const GeometryType::ShapeFunctionsGradientsType& DNu_Dxi = GetGeometry().ShapeFunctionsLocalGradients ( mThisIntegrationMethod );
    const GeometryType::ShapeFunctionsGradientsType& DNo_Dxi = mThisGeometryOther->ShapeFunctionsLocalGradients ( mThisIntegrationMethod );

    const Matrix& Nu_container = GetGeometry().ShapeFunctionsValues ( mThisIntegrationMethod );
    const Matrix& No_container = mThisGeometryOther->ShapeFunctionsValues ( mThisIntegrationMethod );

    Vector Help_R_1 ( mMatSizeU );
    Vector Help_R_2 ( mMatSizeO ); 
    Matrix Help_K_UU ( mMatSizeU, mMatSizeU );
    Matrix Help_K_UP ( mMatSizeU, mMatSizeO ); 
    Matrix Help_K_PU ( mMatSizeO, mMatSizeU );
    Matrix Help_K_PP ( mMatSizeO, mMatSizeO ); 

    // ++++++ RightHandSideVector: R ++++++
    if ( CalculateResidualVectorFlag == true )
    {
        noalias ( Help_R_1 ) = ZeroVector ( mMatSizeU );
        noalias ( Help_R_2 ) = ZeroVector ( mMatSizeO ); 
    }

    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        noalias ( Help_K_UU ) = ZeroMatrix ( mMatSizeU, mMatSizeU );
        noalias ( Help_K_UP ) = ZeroMatrix ( mMatSizeU, mMatSizeO ); 
        noalias ( Help_K_PU ) = ZeroMatrix ( mMatSizeO, mMatSizeU );
        noalias ( Help_K_PP ) = ZeroMatrix ( mMatSizeO, mMatSizeO ); 
    }

    double weight, DetJ = 0.0;
    Vector Nu ( mNodesNumberDisp ), No ( mNodesNumberOther );
    Matrix DNu_DX ( mNodesNumberDisp, mDimension ), DNo_DX ( mNodesNumberOther, mDimension );

    double Div_dot_u, p, dot_p;
    Vector Grad_p ( mDimension );
    Matrix Grad_u ( mDimension, mDimension );

    // adding for calling constitutive law
    Vector strainVector ( 6 );
    Vector stressVectorEff ( 6 );
    Matrix CtanEff ( 6, 6 );
    
    Matrix InternalVariables ( 2, 2 ); 
    
    // START GAUSS INTEGRATION:
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        weight = integration_points[PointNumber].Weight();
        DetJ = mDetJ0[PointNumber]; 

        noalias ( DNu_DX ) = prod ( DNu_Dxi[PointNumber], mInvJ0[PointNumber] );
        noalias ( DNo_DX ) = prod ( DNo_Dxi[PointNumber], mInvJ0[PointNumber] );
        noalias ( Nu ) = row ( Nu_container, PointNumber );
        noalias ( No ) = row ( No_container, PointNumber );

        noalias ( Grad_u ) = GetGradu ( DNu_DX );
        Div_dot_u = GetDivdotu ( DNu_DX );

        p = Getp ( No );
        dot_p = Getdotp ( No );
        noalias ( Grad_p ) = GetGradp ( DNo_DX );
 
        // calling constitutive law for computing effective stress,
        noalias ( strainVector ) = GetstrainVector ( Grad_u ); // Compressive !
        noalias ( stressVectorEff ) = ZeroVector ( 6 );
        noalias ( CtanEff ) = ZeroMatrix ( 6, 6 );
	
	 //Row 1 -- INPUT: n0;    Row 2 -- OUTPUT: K and G
	noalias ( InternalVariables ) = ZeroMatrix (2, 2); 
	InternalVariables(0,0) = GetPorosity ( Grad_u, p, 0 );
 
	bool onlyPN1 = false;
//	if ( PointNumber < 1 && GetValue(KRATOS_WATCH_FLAG)==1 )
	if ( GetValue(KRATOS_WATCH_FLAG)==1 )
	  onlyPN1 = true;
	
        // INPUT: strainVector (Tensile-Positive !!);  OUTPUT: stressVector [kPa] (Tensile-Positive !!) , CtanEff [kPa]
        mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse ( strainVector, InternalVariables, stressVectorEff, CtanEff, rProcessInfo, GetProperties(), GetGeometry(), Nu, true, 1, onlyPN1 ); 
	
// 	mK = InternalVariables(1,0);
// 	mG = InternalVariables(1,1);
	 
        if ( CalculateResidualVectorFlag == true )
        {
            //Calculation of spatial load vector
            AddInternalForcesToRHS1 ( Help_R_1, weight, DetJ, Grad_u, Grad_p, p, Div_dot_u, dot_p, Nu, DNu_DX, stressVectorEff );
            AddInternalForcesToRHS2 ( Help_R_2, weight, DetJ, Grad_u, Grad_p, p, Div_dot_u, dot_p, No, DNo_DX ); 
        }

        if ( CalculateStiffnessMatrixFlag == true )
        {
            //Calculation of spatial Stiffnes and Mass Matrix
            CalculateStiffnessMatrixUU ( Help_K_UU, weight, DetJ, Grad_u, Grad_p, p, Div_dot_u, dot_p, Nu, DNu_DX, CtanEff );
            CalculateStiffnessMatrixUP ( Help_K_UP, weight, DetJ, Grad_u, Grad_p, p, Div_dot_u, dot_p, Nu, DNu_DX, No ); 
            CalculateStiffnessMatrixPU ( Help_K_PU, weight, DetJ, Grad_u, Grad_p, p, Div_dot_u, dot_p, No, DNu_DX );
            CalculateStiffnessMatrixPP ( Help_K_PP, weight, DetJ, Grad_u, Grad_p, p, Div_dot_u, dot_p, No, DNo_DX ); 
        }

    }// END GAUSS INTEGRATION.

// if (GetValue(KRATOS_WATCH_FLAG)==1 )
// {
//   KRATOS_WATCH(Help_R_1);
//   KRATOS_WATCH(Help_K_UU);
// }

    // Assemble Kt and Ri : { ux_n1, uy_n1, uz_n1, ux_n2, uy_n2, uz_n2, ... n_n1, p_n1, t_n1, n_n2, p_n2, t_n2, ...}
    if ( CalculateStiffnessMatrixFlag == true )
    {
        for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
            for ( unsigned int k = 0; k < mDimension; k++ )
            {
                for ( unsigned int j = ( mNodesDispMin - 1 ); j < mNodesDispMax; j++ )
                    for ( unsigned int l = 0; l < mDimension; l++ )
                        rLeftHandSideMatrix ( i*mDimension + k, j*mDimension + l ) += Help_K_UU ( i * mDimension + k, j * mDimension + l );

                for ( unsigned int j = ( mNodesOtherMin - 1 ); j < mNodesOtherMax; j++ )
                {                                                                                                             
                    rLeftHandSideMatrix ( i*mDimension + k, mAddIndexU + j*mNumberO ) += Help_K_UP ( i * mDimension + k, j );
                }
            }

        for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
        {
            for ( unsigned int j = ( mNodesDispMin - 1 ); j < mNodesDispMax; j++ )
                for ( unsigned int l = 0; l < mDimension; l++ ) 
                    rLeftHandSideMatrix ( mAddIndexU + i*mNumberO, j*mDimension + l ) += Help_K_PU ( i, j * mDimension + l );  

            for ( unsigned int j = ( mNodesOtherMin - 1 ); j < mNodesOtherMax; j++ )
            {
                rLeftHandSideMatrix ( mAddIndexU + i*mNumberO, mAddIndexU + j*mNumberO ) += Help_K_PP ( i, j ); 
            }
        }
    }

    if ( CalculateResidualVectorFlag == true )
    {
        for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
            for ( unsigned int k = 0; k < mDimension; k++ )
                rRightHandSideVector ( i*mDimension + k ) -= Help_R_1 ( i * mDimension + k );

        for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ ) 
            rRightHandSideVector ( mAddIndexU + i*mNumberO ) -= Help_R_2 ( i );  
    }


    KRATOS_CATCH ( "" )
}


void UnfrozenSoil::DampMatrix ( MatrixType& rDampMatrix, ProcessInfo& rProcessInfo )
{
    KRATOS_TRY

    if ( rDampMatrix.size1() != mMatSize )
        rDampMatrix.resize ( mMatSize, mMatSize );
    noalias ( rDampMatrix ) = ZeroMatrix ( mMatSize, mMatSize );

    //reading integration points and local Gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints ( mThisIntegrationMethod );
    const GeometryType::ShapeFunctionsGradientsType& DNu_Dxi = GetGeometry().ShapeFunctionsLocalGradients ( mThisIntegrationMethod );
    const GeometryType::ShapeFunctionsGradientsType& DNo_Dxi = mThisGeometryOther->ShapeFunctionsLocalGradients ( mThisIntegrationMethod );
    const Matrix& Nu_container = GetGeometry().ShapeFunctionsValues ( mThisIntegrationMethod );
    const Matrix& No_container = mThisGeometryOther->ShapeFunctionsValues ( mThisIntegrationMethod );

    Matrix Help_D_PU ( mMatSizeO, mMatSizeU );
    Matrix Help_D_PP ( mMatSizeO, mMatSizeO ); 
    noalias ( Help_D_PU ) = ZeroMatrix ( mMatSizeO, mMatSizeU );
    noalias ( Help_D_PP ) = ZeroMatrix ( mMatSizeO, mMatSizeO ); 

    double weight, DetJ = 0.0;
    Vector Nu ( mNodesNumberDisp ), No ( mNodesNumberOther );
    Matrix DNu_DX ( mNodesNumberDisp, mDimension ), DNo_DX ( mNodesNumberOther, mDimension );

    double Div_dot_u, p, dot_p;
    Vector Grad_p ( mDimension );
    Matrix Grad_u ( mDimension, mDimension );

    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        weight = integration_points[PointNumber].Weight();
        DetJ = mDetJ0[PointNumber]; 

        noalias ( DNu_DX ) = prod ( DNu_Dxi[PointNumber], mInvJ0[PointNumber] );
        noalias ( DNo_DX ) = prod ( DNo_Dxi[PointNumber], mInvJ0[PointNumber] );
        noalias ( Nu ) = row ( Nu_container, PointNumber );
        noalias ( No ) = row ( No_container, PointNumber );

        noalias ( Grad_u ) = GetGradu ( DNu_DX );
        Div_dot_u = GetDivdotu ( DNu_DX );

        p = Getp ( No );
        dot_p = Getdotp ( No );
        noalias ( Grad_p ) = GetGradp ( DNo_DX );
 
        //Calculation of spatial damping matrix
        CalculateDampingMatrixPU ( Help_D_PU, weight, DetJ, Grad_u, Grad_p, p, Div_dot_u, dot_p, No, DNu_DX );
        CalculateDampingMatrixPP ( Help_D_PP, weight, DetJ, Grad_u, Grad_p, p, Div_dot_u, dot_p, No, DNo_DX ); 
    }

    for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
    {
        for ( unsigned int j = ( mNodesDispMin - 1 ); j < mNodesDispMax; j++ )
            for ( unsigned int l = 0; l < mDimension; l++ ) 
                rDampMatrix ( mAddIndexU + i*mNumberO, j*mDimension + l ) += Help_D_PU ( i, j * mDimension + l );  

        for ( unsigned int j = ( mNodesOtherMin - 1 ); j < mNodesOtherMax; j++ ) 
            rDampMatrix ( mAddIndexU + i*mNumberO, mAddIndexU + j*mNumberO ) += Help_D_PP ( i, j );  
    }
    KRATOS_CATCH ( "" )
}


UnfrozenSoil::IntegrationMethod UnfrozenSoil::GetIntegrationMethod()
{
    return mThisIntegrationMethod;
}

void UnfrozenSoil::CalculateRightHandSide ( VectorType& rRightHandSideVector, ProcessInfo& rProcessInfo )
{
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll ( temp, rRightHandSideVector, rProcessInfo, CalculateStiffnessMatrixFlag,  CalculateResidualVectorFlag );
}


void UnfrozenSoil::CalculateLocalSystem ( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rProcessInfo )
{
    //calculation flamgS
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;

    CalculateAll ( rLeftHandSideMatrix, rRightHandSideVector, rProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
}

void UnfrozenSoil::EquationIdVector ( EquationIdVectorType& rResult, ProcessInfo& ProcessInfo )
{
    if ( rResult.size() != mMatSize )
        rResult.resize ( mMatSize );

    for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
    {
        rResult[i*mDimension]   = GetGeometry()[i].GetDof ( DISPLACEMENT_X ).EquationId();
        rResult[i*mDimension+1] = GetGeometry()[i].GetDof ( DISPLACEMENT_Y ).EquationId();
        rResult[i*mDimension+2] = GetGeometry()[i].GetDof ( DISPLACEMENT_Z ).EquationId();
    }

    for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
    {
        rResult[mAddIndexU+i*mNumberO]   = GetGeometry()[i].GetDof ( WATER_PRESSURE ).EquationId(); 
    }
}


void UnfrozenSoil::GetDofList ( DofsVectorType& ElementalDofList, ProcessInfo& ProcessInfo )
{
    ElementalDofList.resize ( 0 );

    for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
    {
        ElementalDofList.push_back ( GetGeometry()[i].pGetDof ( DISPLACEMENT_X ) );
        ElementalDofList.push_back ( GetGeometry()[i].pGetDof ( DISPLACEMENT_Y ) );
        ElementalDofList.push_back ( GetGeometry()[i].pGetDof ( DISPLACEMENT_Z ) );
    }

    for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ ) 
        ElementalDofList.push_back ( GetGeometry()[i].pGetDof ( WATER_PRESSURE ) );  
}


void UnfrozenSoil::GetValuesVector ( Vector& values, int Step )
{
    if ( values.size() != mMatSize )
        values.resize ( mMatSize );

    for ( unsigned int i = ( mNodesDispMin - 1 );i < mNodesDispMax;i++ )
    {
        values ( i*mDimension )   = GetGeometry()[i].GetSolutionStepValue ( DISPLACEMENT_X, Step );
        values ( i*mDimension + 1 ) = GetGeometry()[i].GetSolutionStepValue ( DISPLACEMENT_Y, Step );
        values ( i*mDimension + 2 ) = GetGeometry()[i].GetSolutionStepValue ( DISPLACEMENT_Z, Step );
    }

    for ( unsigned int i = ( mNodesOtherMin - 1 );i < mNodesOtherMax;i++ ) 
        values ( mAddIndexU + i*mNumberO )   = GetGeometry()[i].GetSolutionStepValue ( WATER_PRESSURE, Step );  
}


///-//////////////////////////////////////////////////////////////////////////////////////
///-//////////////////////////////////////////////////////////////////////////////////////
/// R1 ***************************
void UnfrozenSoil::AddInternalForcesToRHS1 ( Vector& Help_R_1, double weight, double DetJ, Matrix& Grad_u, Vector& Grad_p, double p, double Div_dot_u, double dot_p, const Vector& Nu, const Matrix& DNu_DX, Vector& stressVectorEff )
{
    Matrix Bu ( 6, mNodesNumberDisp*mDimension );
    noalias ( Bu ) = GetBu ( DNu_DX );

    double ML = GetWaterMass ( Grad_u, p, 0 ); 
    double b = GetBiotCoefficient (); 
     
    for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
        for ( unsigned int m = 0; m < mDimension; m++ )
        {
            Help_R_1 ( i*mDimension + m ) += Nu ( i ) * (
                                                 ( ( 1.0 - mn0 ) * mrhoS0 + ML ) * mGravity ( m )
                                             ) * weight * DetJ * mScaleU;

            Help_R_1 ( i*mDimension + m ) += DNu_DX ( i, m ) * (
                                                 b * p
                                             ) * weight * DetJ * mScaleU;  

            for ( unsigned int k = 0; k < 6; k++ )
                Help_R_1 ( i*mDimension + m ) -= Bu ( k, i * mDimension + m ) * (
                                                     /*1000.0 * */stressVectorEff ( k )  // convert kPa to Pa 
                                                 ) * weight * DetJ * mScaleU; 
        }
}

void UnfrozenSoil::CalculateStiffnessMatrixUU ( Matrix& Help_K_UU, double weight, double DetJ, Matrix& Grad_u, Vector& Grad_p, double p, double Div_dot_u, double dot_p, const Vector& Nu, const Matrix& DNu_DX, Matrix& CtanEff )
{
    Matrix Bu ( 6, mNodesNumberDisp*mDimension );
    noalias ( Bu ) = GetBu ( DNu_DX );
    
    double DML_De = GetWaterMass ( Grad_u, p, 1);

    for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
        for ( unsigned int m = 0; m < mDimension; m++ )
            for ( unsigned int j = ( mNodesDispMin - 1 ); j < mNodesDispMax; j++ )
                for ( unsigned int n = 0; n < mDimension; n++ )
                {
                    Help_K_UU ( i*mDimension + m, j*mDimension + n ) += Nu ( i ) * (
                                DML_De * mGravity ( m )
                            ) * DNu_DX ( j, n ) * weight * DetJ * mScaleU;

                    for ( unsigned int k = 0; k < 6; k++ )
                        for ( unsigned int l = 0; l < 6; l++ )
                        {
                            Help_K_UU ( i*mDimension + m, j*mDimension + n ) -= Bu ( k, i * mDimension + m ) * (
                                       /* 1000.0 * */CtanEff ( k, l )  // convert kPa to Pa 
                                    ) * Bu ( l, j * mDimension + n ) * weight * DetJ * mScaleU;
                        }

                }
}

void UnfrozenSoil::CalculateStiffnessMatrixUP ( Matrix& Help_K_UP, double weight, double DetJ, Matrix& Grad_u, Vector& Grad_p, double p, double Div_dot_u, double dot_p, const Vector& Nu, const Matrix& DNu_DX, const Vector& No )
{
    double DML_Dp = GetWaterMass ( Grad_u, p, 2);
    double b = GetBiotCoefficient (); 

    for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
        for ( unsigned int m = 0; m < mDimension; m++ )
            for ( unsigned int j = ( mNodesOtherMin - 1 ); j < mNodesOtherMax; j++ )
            {
                Help_K_UP ( i*mDimension + m, j ) += Nu ( i ) * (
                                                         DML_Dp * mGravity ( m )
                                                     ) * No ( j ) * weight * DetJ * mScaleU;

                Help_K_UP ( i*mDimension + m, j ) += DNu_DX ( i , m ) * (
                                                         b 
                                                     ) * No ( j ) * weight * DetJ * mScaleU;
            }
} 

///-//////////////////////////////////////////////////////////////////////////////////////
/// R2 **************************
void UnfrozenSoil::AddInternalForcesToRHS2 ( Vector& Help_R_2, double weight, double DetJ, Matrix& Grad_u, Vector& Grad_p, double p, double Div_dot_u, double dot_p, const Vector& No, const Matrix& DNo_DX )
{
    double DML_De = GetWaterMass ( Grad_u, p, 1 );
    double DML_Dp = GetWaterMass ( Grad_u, p, 2 ); 

    Vector wL ( mDimension );
    noalias ( wL ) = GetWaterFlow ( Grad_p, 0 );

    for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
    {
        Help_R_2 ( i ) += No ( i ) * (
                              DML_De * Div_dot_u
                          ) * weight * DetJ * mScaleP;

        Help_R_2 ( i ) += No ( i ) * (
                              DML_Dp * dot_p
                          ) * weight * DetJ * mScaleP;

        for ( unsigned int k = 0; k < mDimension; k++ )
        {
            // %%% R2---5
            Help_R_2 ( i ) -= DNo_DX ( i, k ) * (
                                  wL ( k )
                              ) * weight * DetJ * mScaleP;
        }
    }
}

void UnfrozenSoil::CalculateStiffnessMatrixPU ( Matrix& Help_K_PU, double weight, double DetJ, Matrix& Grad_u, Vector& Grad_p, double p, double Div_dot_u, double dot_p, const Vector& No, const Matrix& DNu_DX )
{
    double D2ML_DeDp = GetWaterMass ( Grad_u, p, 3 ); 

    for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
        for ( unsigned int j = ( mNodesDispMin - 1 ); j < mNodesDispMax; j++ )
            for ( unsigned int n = 0; n < mDimension; n++ )
            {
                Help_K_PU ( i, j*mDimension + n ) += No ( i ) * (
                                                         D2ML_DeDp * dot_p
                                                     ) * DNu_DX ( j, n ) * weight * DetJ * mScaleU;
 
            }
}

void UnfrozenSoil::CalculateStiffnessMatrixPP ( Matrix& Help_K_PP, double weight, double DetJ, Matrix& Grad_u, Vector& Grad_p, double p, double Div_dot_u, double dot_p, const Vector& No, const Matrix& DNo_DX )
{
    double D2ML_DeDp = GetWaterMass ( Grad_u, p, 3 );
    double D2ML_Dp2 = GetWaterMass ( Grad_u, p, 4 ); 
    double DwL_DGradp = GetWaterFlow ( Grad_p, 1 )[0];

    for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
        for ( unsigned int j = ( mNodesOtherMin - 1 ); j < mNodesOtherMax; j++ )
        {
            Help_K_PP ( i, j ) += No ( i ) * (
                                      D2ML_DeDp * Div_dot_u
                                  ) * No ( j ) * weight * DetJ * mScaleP;

            Help_K_PP ( i, j ) += No ( i ) * (
                                      D2ML_Dp2 * dot_p
                                  ) * No ( j ) * weight * DetJ * mScaleP; 


            for ( unsigned int k = 0; k < mDimension; k++ )
                Help_K_PP ( i, j ) -= DNo_DX ( i , k ) * (
                                          DwL_DGradp
                                      ) * DNo_DX ( j , k ) * weight * DetJ * mScaleP;
        }
}
 
// Damping
void UnfrozenSoil::CalculateDampingMatrixPU ( Matrix& Help_D_PU, double weight, double DetJ, Matrix& Grad_u, Vector& Grad_p, double p, double Div_dot_u, double dot_p, const Vector& No, const Matrix& DNu_DX )
{
    double DML_De = GetWaterMass ( Grad_u, p, 1);

    for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
        for ( unsigned int j = ( mNodesDispMin - 1 ); j < mNodesDispMax; j++ )
            for ( unsigned int n = 0; n < mDimension; n++ )
                Help_D_PU ( i, j*mDimension + n ) += No ( i ) * (
                                                         DML_De
                                                     ) * DNu_DX ( j, n ) * weight * DetJ * mScaleU;
}

void UnfrozenSoil::CalculateDampingMatrixPP ( Matrix& Help_D_PP, double weight, double DetJ, Matrix& Grad_u, Vector& Grad_p, double p, double Div_dot_u, double dot_p, const Vector& No, const Matrix& DNo_DX )
{
    double DML_Dp = GetWaterMass ( Grad_u, p, 2);

    for ( unsigned int i = ( mNodesOtherMin - 1 ); i < mNodesOtherMax; i++ )
        for ( unsigned int j = ( mNodesOtherMin - 1 ); j < mNodesOtherMax; j++ )
            Help_D_PP ( i, j ) += No ( i ) * (
                                      DML_Dp
                                  ) * No ( j ) * weight * DetJ * mScaleP;
}
 
 
///-//////////////////////////////////////////////////////////////////////////////////////
///  CONSTITUTIVE QUANTITES
///-//////////////////////////////////////////////////////////////////////////////////////
/// KnoneckerDelta
double UnfrozenSoil::KnoneckerDelta ( int i, int j )
{
    if ( i == j )
        return 1.0;
    else
        return 0.0;
}

/// Bu operator
Matrix UnfrozenSoil::GetBu ( Matrix DNu_DX )
{
    Matrix result ( 6, mNodesNumberDisp*mDimension );
    noalias ( result ) = ZeroMatrix ( 6, mNodesNumberDisp * mDimension );

    for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
    {
        result ( 0, i*3 ) = DNu_DX ( i, 0 );
        result ( 1, i*3 + 1 ) = DNu_DX ( i, 1 );
        result ( 2, i*3 + 2 ) = DNu_DX ( i, 2 );
        result ( 3, i*3 ) = DNu_DX ( i, 1 );
        result ( 3, i*3 + 1 ) = DNu_DX ( i, 0 );
        result ( 4, i*3 + 1 ) = DNu_DX ( i, 2 );
        result ( 4, i*3 + 2 ) = DNu_DX ( i, 1 );
        result ( 5, i*3 ) = DNu_DX ( i, 2 );
        result ( 5, i*3 + 2 ) = DNu_DX ( i, 0 );
    }
    return result;
}

///--------------------------------------------------------------------------
/// Displacement
///--------------------------------------------------------------------------
Matrix UnfrozenSoil::GetGradu ( const Matrix& DNu_DX )
{
    Matrix result ( mDimension, mDimension );
    noalias ( result ) = ZeroMatrix ( mDimension, mDimension );

    for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
        for ( unsigned int k = 0; k < mDimension; k++ )
        {
            result ( 0, k ) += GetGeometry()[i].GetSolutionStepValue ( DISPLACEMENT_X ) * DNu_DX ( i, k );
            result ( 1, k ) += GetGeometry()[i].GetSolutionStepValue ( DISPLACEMENT_Y ) * DNu_DX ( i, k );
            result ( 2, k ) += GetGeometry()[i].GetSolutionStepValue ( DISPLACEMENT_Z ) * DNu_DX ( i, k );
        }

    return result;
}

double UnfrozenSoil::GetDivdotu ( const Matrix& DNu_DX )
{
    double result = 0.0;

    for ( unsigned int i = ( mNodesDispMin - 1 ); i < mNodesDispMax; i++ )
        for ( unsigned int k = 0; k < mDimension; k++ )
            result  += GetGeometry()[i].GetSolutionStepValue ( DISPLACEMENT_DT )[k] * DNu_DX ( i, k );

    return result;
}

///--------------------------------------------------------------------------
/// Water pressure
///--------------------------------------------------------------------------
double UnfrozenSoil::Getp ( const Vector& No )
{
    double result = 0.0;

    for ( unsigned int i = mNodesOtherMin - 1 ; i < mNodesOtherMax ;i++ )
        result += GetGeometry()[i].GetSolutionStepValue ( WATER_PRESSURE ) * No ( i );

    return result;
}

double UnfrozenSoil::Getdotp ( const Vector& No )
{
    double result = 0.0;

    for ( unsigned int i = mNodesOtherMin - 1 ; i < mNodesOtherMax ;i++ )
        result += GetGeometry()[i].GetSolutionStepValue ( WATER_PRESSURE_DT ) * No ( i );

    return result;
}

Vector UnfrozenSoil::GetGradp ( const Matrix& DNo_DX )
{
    Vector result ( mDimension );
    noalias ( result ) = ZeroVector ( mDimension );

    for ( unsigned int i = mNodesOtherMin - 1 ; i < mNodesOtherMax ;i++ )
        for ( unsigned int k = 0; k < mDimension; k++ )
            result[k] += GetGeometry()[i].GetSolutionStepValue ( WATER_PRESSURE ) * DNo_DX ( i, k );

    return result;
}
 
/// ++++++++++++++ Group c1 +++++++++++++++++++
/// +++++++++++++++++++++++++++++++++++++++++++
/// strain: linear strain tensor
Vector UnfrozenSoil::GetstrainVector ( Matrix Grad_u )
{
    Vector result ( 6 );
    noalias(result) = ZeroVector (6);
    result ( 0 ) = Grad_u ( 0, 0 );
    result ( 1 ) = Grad_u ( 1, 1 );
    result ( 2 ) = Grad_u ( 2, 2 );
    result ( 3 ) = Grad_u ( 0, 1 ) + Grad_u ( 1, 0 );
    result ( 4 ) = Grad_u ( 1, 2 ) + Grad_u ( 2, 1 );
    result ( 5 ) = Grad_u ( 0, 2 ) + Grad_u ( 2, 0 );

//     result ( 0 ) += 0.5*( Grad_u ( 0, 0 )*Grad_u ( 0, 0 ) +  Grad_u ( 1, 0 )*Grad_u ( 1, 0 )+  Grad_u ( 2, 0 )*Grad_u ( 2, 0 ) );
//     result ( 1 ) += 0.5*( Grad_u ( 0, 1 )*Grad_u ( 0, 1 ) +  Grad_u ( 1, 1 )*Grad_u ( 1, 1 )+  Grad_u ( 2, 1 )*Grad_u ( 2, 1 ) );
//     result ( 2 ) += 0.5*( Grad_u ( 0, 2 )*Grad_u ( 0, 2 ) +  Grad_u ( 1, 2 )*Grad_u ( 1, 2 )+  Grad_u ( 2, 2 )*Grad_u ( 2, 2 ) );
//     result ( 3 ) += Grad_u ( 0, 0 )*Grad_u ( 0, 1 ) +  Grad_u ( 1, 0 )*Grad_u ( 1, 1 )+  Grad_u ( 2, 0 )*Grad_u ( 2, 1 );
//     result ( 4 ) += Grad_u ( 0, 1 )*Grad_u ( 0, 2 ) +  Grad_u ( 1, 1 )*Grad_u ( 1, 2 )+  Grad_u ( 2, 1 )*Grad_u ( 2, 2 );
//     result ( 5 ) += Grad_u ( 0, 0 )*Grad_u ( 0, 2 ) +  Grad_u ( 1, 0 )*Grad_u ( 1, 2 )+  Grad_u ( 2, 0 )*Grad_u ( 2, 2 );
    return result;
}

/// strainVol: volumetric strain
double UnfrozenSoil::GetstrainVol ( Matrix Grad_u )
{
    return  Grad_u ( 0, 0 ) + Grad_u ( 1, 1 ) + Grad_u ( 2, 2 );
}
 
/// b: Biot coeffecient
double UnfrozenSoil::GetBiotCoefficient ()
{
//     double alpha = 3.0 * mkC / ( 3.0 * mkC + 4.0 * mgC );
//     double K = mkL * ( 1.0 + ( 1.0 - mn0 ) * ( 1.0 - mkL / mkS ) / ( mkL / mkS + alpha * mn0 * ( 1.0 - mkL / mkS ) ) );  
//  
//     return  1.0 - K / mkS;   
    return  1.0;   
}

/// n: current porosity
double UnfrozenSoil::GetPorosity ( Matrix Grad_u, double p, int index )
{
    double result = 0.0; 
    double b = GetBiotCoefficient (); 
    double strainVol = GetstrainVol ( Grad_u );

    switch ( index )
    {
    case 0: // n
        result = mn0 + b * strainVol + ( b - mn0 ) / mkS * p;
        break;
    case 1: // Dn_De
        result = b;
        break;
    case 2: // Dn_Dp
        result = ( b - mn0 ) / mkS ;
        break;
    }
    return result;
}


/// rhoL: current water density
// double UnfrozenSoil::GetrhoL ( double p )
double UnfrozenSoil::GetWaterDensity ( double p, int index )
{
    switch ( index )
    {
    case 0: // rhoL
        return mrhoL0 * ( 1.0 + p / mkL );
    case 1: // DrhoL_Dp
        return mrhoL0 / mkL; 
    }
}
 

/// ML: water mass
double UnfrozenSoil::GetWaterMass ( Matrix Grad_u, double p, int index )
{
    double n = GetPorosity ( Grad_u, p, 0 ) ;
    double rhoL = GetWaterDensity ( p, 0 );

    switch ( index )
    {
    case 0: // ML
        return rhoL * n;

    case 1: // DML_De
    {
        double Dn_De = GetPorosity ( Grad_u, p, 1 );
        return rhoL * Dn_De;
    }

    case 2: // DML_Dp
    {
        double Dn_Dp = GetPorosity ( Grad_u, p, 2 );
        double DrhoL_Dp = GetWaterDensity ( p, 1 );
        return DrhoL_Dp * n + rhoL * Dn_Dp ;
    } 
    
    case 3: // D2ML_DeDp
    {
        double Dn_De = GetPorosity ( Grad_u, p, 1 );
        double DrhoL_Dp = GetWaterDensity ( p, 1 );
        return DrhoL_Dp * Dn_De ;
    } 
    
    case 4: // D2ML_Dp2
    {
        double Dn_Dp = GetPorosity ( Grad_u, p, 2 );
        double DrhoL_Dp = GetWaterDensity ( p, 1 );
        return 2.0 * DrhoL_Dp * Dn_Dp ;
    } 
    }
} 
/// ++++++++++++++ Group c2 +++++++++++++++++++
/// +++++++++++++++++++++++++++++++++++++++++++
/// wL: liquid water flow
Vector UnfrozenSoil::GetWaterFlow ( Vector Grad_p, int index )
{
    Vector result ( mDimension );
    noalias ( result ) = ZeroVector ( mDimension );  
    
    switch ( index )
    {
    case 0: // wL
        for ( unsigned int k = 0; k < mDimension; k++ )
            result[k] =  mrhoL0 * mkappa0 / metaL * ( -Grad_p ( k ) + mrhoL0 * mGravity ( k ) );
        break;
    case 1: // DwL_DGradp : double
        result[0] = -mrhoL0 * mkappa0 / metaL;
        break; 
    }
    return result;
}
 

/// - WriteNodalResults
void UnfrozenSoil::FinalizeSolutionStep ( ProcessInfo& ProcessInfo )
{
    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        mConstitutiveLawVector[i]->FinalizeSolutionStep( GetProperties(),
                GetGeometry(),
                row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ),
                ProcessInfo );

    Matrix DNo_DX ( mNodesNumberOther, mDimension );
    Matrix DNo_Dxi ( mNodesNumberOther, mDimension );
    Matrix nodesLocalCoords ( mNodesNumberOther, mDimension );

    for ( unsigned int i = ( mNodesOtherMin - 1 ) ; i < mNodesOtherMax ;i++ )
        for ( unsigned int m = 0; m < mDimension; m++ )
        {
            nodesLocalCoords ( i, m ) = 1.0;
            if ( m >= i )
                nodesLocalCoords ( i, m ) = -1.0;
        }
    nodesLocalCoords ( 3, 0 ) = -1.0;
    nodesLocalCoords ( 3, 2 ) = -1.0;
    nodesLocalCoords ( 4, 0 ) = -1.0;
    nodesLocalCoords ( 4, 1 ) = -1.0;
    nodesLocalCoords ( 5, 1 ) = -1.0;
    nodesLocalCoords ( 7, 0 ) = -1.0;
 
    double p;
    Vector Grad_p(mDimension),Grad_t(mDimension);
    Matrix Grad_u(mDimension,mDimension);

    for ( unsigned int i = ( mNodesOtherMin - 1 ) ; i < mNodesOtherMax ;i++ )
    {
        Vector local_coords ( 3 );
        for ( unsigned int m = 0; m < mDimension; m++ )
            local_coords ( m ) = nodesLocalCoords ( i, m );

        noalias ( DNo_Dxi ) = GetGeometry().ShapeFunctionsLocalGradients ( DNo_Dxi, local_coords );
        noalias ( DNo_DX ) = prod ( DNo_Dxi, mInvJ0[0] );
 
	p= GetGeometry()[i].GetSolutionStepValue ( WATER_PRESSURE );
        noalias(Grad_p)= GetGradp(DNo_DX); 
        noalias(Grad_u)= GetGradu(DNo_DX);
 
        GetGeometry()[i].GetSolutionStepValue ( POROSITY ) = GetPorosity( Grad_u, p, 0 );
        GetGeometry()[i].GetSolutionStepValue ( LINEAR_STRAIN ) = GetstrainVol (Grad_u); 
        GetGeometry()[i].GetSolutionStepValue ( WATER_FLOW ) = GetWaterFlow ( Grad_p, 0 ); 
    }

    if ( GetGeometry().size() == 20 )
    {
        Interpolate ( POROSITY, ProcessInfo ); 
        Interpolate ( WATER_PRESSURE, ProcessInfo ); 
        Interpolate ( LINEAR_STRAIN, ProcessInfo );

        Interpolate ( DISPLACEMENT, ProcessInfo );
        Interpolate ( WATER_FLOW, ProcessInfo ); 
    }
}

void UnfrozenSoil::Interpolate ( const Variable<double>& rVariable, const ProcessInfo& rProcessInfo )
{
    for ( unsigned int k = 0; k < 4; k++ )
    {
        GetGeometry()[8+k].GetSolutionStepValue ( rVariable ) = ( GetGeometry()[k%4].GetSolutionStepValue ( rVariable ) + GetGeometry()[(k+1)%4].GetSolutionStepValue ( rVariable ) ) / 2.0;
        GetGeometry()[12+k].GetSolutionStepValue ( rVariable ) = ( GetGeometry()[4+k%4].GetSolutionStepValue ( rVariable ) + GetGeometry()[4+(k+1)%4].GetSolutionStepValue ( rVariable ) ) / 2.0;
        GetGeometry()[16+k].GetSolutionStepValue ( rVariable ) = ( GetGeometry()[k%4].GetSolutionStepValue ( rVariable ) + GetGeometry()[4+k%4].GetSolutionStepValue ( rVariable ) ) / 2.0;
    }
}

void UnfrozenSoil::Interpolate ( const Variable<Kratos::array_1d<double, 3> >& rVariable, const ProcessInfo& rProcessInfo )
{
    for ( unsigned int k = 0; k < 4; k++ )
    {
        GetGeometry()[8+k].GetSolutionStepValue ( rVariable ) = ( GetGeometry()[k%4].GetSolutionStepValue ( rVariable ) + GetGeometry()[(k+1)%4].GetSolutionStepValue ( rVariable ) ) / 2.0;
        GetGeometry()[12+k].GetSolutionStepValue ( rVariable ) = ( GetGeometry()[4+k%4].GetSolutionStepValue ( rVariable ) + GetGeometry()[4+(k+1)%4].GetSolutionStepValue ( rVariable ) ) / 2.0;
        GetGeometry()[16+k].GetSolutionStepValue ( rVariable ) = ( GetGeometry()[k%4].GetSolutionStepValue ( rVariable ) + GetGeometry()[4+k%4].GetSolutionStepValue ( rVariable ) ) / 2.0;
    }
}

///-////////////////////////////////////
/////
///- //////////////////////////////////////////////////////
void UnfrozenSoil::GetValueOnIntegrationPoints ( const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rProcessInfo )
{
    if ( rValues.size() != GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size() )
        rValues.resize( GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(), false );

    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
	rValues[i] = mConstitutiveLawVector[i]->GetValue( rVariable, rValues[i] ); 
	
}
void UnfrozenSoil::GetValueOnIntegrationPoints ( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rProcessInfo )
{
        const unsigned int& size = GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size();

        if ( rValues.size() != size )
            rValues.resize( size );
	
    //To Plot Internal variables
    if ( rVariable == INTERNAL_VARIABLES ) 
        for ( unsigned int i = 0; i < size; i++ )
        {
            if ( rValues[i].size() != 6 )
                rValues[i].resize( 6 );
            noalias( rValues[i] ) = mConstitutiveLawVector[i]->GetValue( INTERNAL_VARIABLES, rValues[i] );
        } 

    //To Plot Stresses
    if ( rVariable == INSITU_STRESS ) 
        for ( unsigned int i = 0; i < size; i++ )
        {
            if ( rValues[i].size() != 6 )
                rValues[i].resize( 6 );
            noalias( rValues[i] ) = mConstitutiveLawVector[i]->GetValue( INSITU_STRESS, rValues[i] );
        } 
        
}
void UnfrozenSoil::GetValueOnIntegrationPoints ( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rProcessInfo ) {}
void UnfrozenSoil::SetValueOnIntegrationPoints ( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rProcessInfo ) {}
void UnfrozenSoil::SetValueOnIntegrationPoints ( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rProcessInfo ) {}
void UnfrozenSoil::CalculateOnIntegrationPoints ( const Variable<double>& rVariable, Vector& Output, const ProcessInfo& rProcessInfo ) {}
int UnfrozenSoil::Check ( const Kratos::ProcessInfo& rProcessInfo )
{
    return 0;
}
} // Namespace Kratos
