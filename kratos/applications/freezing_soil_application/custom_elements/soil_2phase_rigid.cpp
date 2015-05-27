/*
==============================================================================
KratosR1StructuralApplication
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
//   Last modified by:    $Author: mengmeng $
//   Date:                $Date: 2009-02-23 16:02:32 $
//   Revision:            $Revision: 1.1 $
//
//


// System includes


// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/soil_2phase_rigid.h"
#include "utilities/math_utils.h"
//#include "includes/constitutive_law.h"
#include "geometries/hexahedra_3d_8.h"
#include "freezing_soil_application.h"

namespace Kratos
{
///done
Soil2PhaseRigid::Soil2PhaseRigid( IndexType NewId, GeometryType::Pointer pGeometry )
        : Element( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}
///done
Soil2PhaseRigid::Soil2PhaseRigid( IndexType NewId,
                                  GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
        : Element( NewId, pGeometry, pProperties )
{
    //setting up the nodal degrees of freedom
    //DOFs at the end of time step
    //All calculations are made on the general midpoint alpha
    //Variables DOF_ALPHA are updated in the scheme####

    if ( GetGeometry().size() == 27 || GetGeometry().size() == 20 || GetGeometry().size() == 8 )
    {
        if ( GetGeometry().size() == 8 )
        {
//             std::cout << "+++ Soil2PhaseRigid3D8N(NumericalTangent) +++" << std::endl;
            mNodesMin = 1; // Disp: displacement
            mNodesMax = 8;
            mThisIntegrationMethod = GeometryData::GI_GAUSS_3;
        }

        if ( GetGeometry().size() == 20 )
        {
//             std::cout << "+++ Soil2PhaseRigid3D20N(NumericalTangent) +++" << std::endl;
            mNodesMin = 1;
            mNodesMax = 20;
            mThisIntegrationMethod = GeometryData::GI_GAUSS_3;
        }

        if ( GetGeometry().size() == 27 )
        {
//             std::cout << "+++ Soil2PhaseRigid3D27N(NumericalTangent) +++" << std::endl;
            mNodesMin = 1; // Disp: displacement
            mNodesMax = 27;
            mThisIntegrationMethod = GeometryData::GI_GAUSS_3;
        }

        if ( GetGeometry().size() == 10 )
        {
//             std::cout << "+++ Soil2PhaseRigid3D10N(NumericalTangent) +++" << std::endl;
            mNodesMin = 1;
            mNodesMax = 10;
            mThisIntegrationMethod = GeometryData::GI_GAUSS_5;
        }
    }

    else
        KRATOS_THROW_ERROR( std::logic_error, "This element matches only with a hybrid hexaeder 8 geometry" , *this );
}
///done
Element::Pointer Soil2PhaseRigid::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new Soil2PhaseRigid( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}
///done
Soil2PhaseRigid::~Soil2PhaseRigid()
{
}
//************************************************************************************
//************************************************************************************
///done
void Soil2PhaseRigid::Initialize()
{
    KRATOS_TRY

    mNodesNumber = mNodesMax - mNodesMin + 1;
    mMatSize = mNodesNumber * 2;
    mDimension = GetGeometry().WorkingSpaceDimension();
    mPorosity =  GetProperties()[POROSITY]; 
    mScale = GetProperties()[SCALE_O];
    mMaterialParameters = GetProperties()[ELEMENT_PARAMETERS];
//     mUnitRatio = GetProperties()[SCALE_U]; // convert m to mm
//     mUnitRatio = 1000.0; // in mm 
    mUnitRatio = GetProperties()[SCALE];				// Unit convertor: from [m] to [mm]
//     mUnitRatio = 1.0; // in m

    for ( unsigned int i = ( mNodesMin - 1 ) ; i < mNodesMax ; i++ )
    {
        ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE_NULL ) = ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE );
        ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE_EINS ) = ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE );
        ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE_DT ) = 0;
        ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE_NULL_DT ) = 0;
        ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE_EINS_DT ) = 0;
        ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE_ACCELERATION ) = 0;
        ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE_NULL_ACCELERATION ) = 0;
        ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE_EINS_ACCELERATION ) = 0;

        ( GetGeometry()[i] ).GetSolutionStepValue( TEMPERATURE_NULL ) = ( GetGeometry()[i] ).GetSolutionStepValue( TEMPERATURE );
        ( GetGeometry()[i] ).GetSolutionStepValue( TEMPERATURE_EINS ) = ( GetGeometry()[i] ).GetSolutionStepValue( TEMPERATURE );
        ( GetGeometry()[i] ).GetSolutionStepValue( TEMPERATURE_DT ) = 0;
        ( GetGeometry()[i] ).GetSolutionStepValue( TEMPERATURE_NULL_DT ) = 0;
        ( GetGeometry()[i] ).GetSolutionStepValue( TEMPERATURE_EINS_DT ) = 0;
        ( GetGeometry()[i] ).GetSolutionStepValue( TEMPERATURE_ACCELERATION ) = 0;
        ( GetGeometry()[i] ).GetSolutionStepValue( TEMPERATURE_NULL_ACCELERATION ) = 0;
        ( GetGeometry()[i] ).GetSolutionStepValue( TEMPERATURE_EINS_ACCELERATION ) = 0;
    }

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );
    mInvJ0.resize( integration_points.size() );

    for ( unsigned int i = 0; i < integration_points.size(); i++ )
    {
        mInvJ0[i].resize( mDimension, mDimension );
        noalias( mInvJ0[i] ) = ZeroMatrix( mDimension, mDimension );
    }

    mDetJ0.resize( integration_points.size() );
    noalias( mDetJ0 ) = ZeroVector( integration_points.size() );


    GeometryType::JacobiansType J0( integration_points.size() );
    J0 = GetGeometry().Jacobian( J0, mThisIntegrationMethod );

    //calculating the inverse J0
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        //calculating and storing inverse of the jacobian and the parameters needed
        MathUtils<double>::InvertMatrix( J0[PointNumber], mInvJ0[PointNumber], mDetJ0[PointNumber] );
    }

    mInitialTemperature.resize( mNodesNumber );

    for ( unsigned int i = ( mNodesMin - 1 ) ; i < mNodesMax ;i++ )
    {
        mInitialTemperature( i ) = GetGeometry()[i].GetSolutionStepValue( TEMPERATURE );
    }

    KRATOS_CATCH( "" )
}

///done
void Soil2PhaseRigid::CalculateAll( MatrixType& rLeftHandSideMatrix,
                                    VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo,
                                    bool CalculateStiffnessMatrixFlag, bool CalculateResidualVectorFlag )
{
    KRATOS_TRY
    //resizing as needed the LHS
    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != mMatSize )
            rLeftHandSideMatrix.resize( mMatSize, mMatSize );
        noalias( rLeftHandSideMatrix ) = ZeroMatrix( mMatSize, mMatSize ); //resetting LHS
    }

    //resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size() != mMatSize )
            rRightHandSideVector.resize( mMatSize );
        noalias( rRightHandSideVector ) = ZeroVector( mMatSize ); //resetting RHS
    }
    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );
    const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );
    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

    Matrix Help_K_PP( mNodesNumber, mNodesNumber );
    Matrix Help_K_PT( mNodesNumber, mNodesNumber );
    Matrix Help_K_TP( mNodesNumber, mNodesNumber );
    Matrix Help_K_TT( mNodesNumber, mNodesNumber );

    Vector Help_R_P( mNodesNumber );
    Vector Help_R_T( mNodesNumber );

    Matrix DN_DX( mNodesNumber, mDimension );
    double weight;
    Vector N( mNodesNumber );

    double p;
    double t;
    double DetJ = 0.0;

    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        noalias( Help_K_PP ) = ZeroMatrix( mNodesNumber, mNodesNumber );
        noalias( Help_K_PT ) = ZeroMatrix( mNodesNumber, mNodesNumber );
        noalias( Help_K_TP ) = ZeroMatrix( mNodesNumber, mNodesNumber );
        noalias( Help_K_TT ) = ZeroMatrix( mNodesNumber, mNodesNumber );
    }
    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
        noalias( Help_R_P ) = ZeroVector( mNodesNumber );
        noalias( Help_R_T ) = ZeroVector( mNodesNumber );
    }

    Vector grad_t( mDimension );
double dot_t;//ADD

    /////////////////////////////////////////////////////////////////////////
    //// Integration in space sum_(beta=0)^(number of quadrature points)
    /////////////////////////////////////////////////////////////////////////
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        if ( DN_DX.size1() != ( mNodesNumber ) || DN_DX.size2() != mDimension )
            DN_DX.resize( mNodesNumber, mDimension );
        noalias( DN_DX ) = ZeroMatrix( mNodesNumber, mDimension );
        noalias( DN_DX ) = prod( DN_De[PointNumber], mInvJ0[PointNumber] );

        weight = integration_points[PointNumber].Weight();
        DetJ = mDetJ0[PointNumber];

        // Shape Functions on current spatial quadrature point
        if ( N.size() != mNodesNumber )
            N.resize( mNodesNumber );
        noalias( N ) = row( Ncontainer, PointNumber );

        p = Getp( N );
        t = Gett( N );
	
dot_t = Getdott( N );//ADD 
    noalias( grad_t ) = Getgradt( DN_DX );

        if ( GetValue( KRATOS_WATCH_FLAG ) == 1 && PointNumber == 8 )
        { 
            std::cout << "3---- t= " << t << ",\t dot_t= " << dot_t << ",\t grad_t= " << grad_t << std::endl;
        }
        
        if ( CalculateStiffnessMatrixFlag == true )
        {
            //Calculation of spatial Stiffnes and Mass Matrix
            CalculateStiffnessMatrixPP( Help_K_PP, DN_DX, N, weight, DetJ, p, t );
            CalculateStiffnessMatrixPT( Help_K_PT, DN_DX, N, weight, DetJ, p, t );
            CalculateStiffnessMatrixTP( Help_K_TP, DN_DX, N, weight, DetJ, p, t );
            CalculateStiffnessMatrixTT( Help_K_TT, DN_DX, N, weight, DetJ, p, t );
        }

        if ( CalculateResidualVectorFlag == true )
        {
            //Calculation of spatial Loadvector
            AddInternalForcesToRHSP( Help_R_P, DN_DX, N, weight, DetJ, p, t );
            AddInternalForcesToRHST( Help_R_T, DN_DX, N, weight, DetJ, p, t );
        }
        ///////////////////////////////////////////////////////////////////////
        // END Integration in space sum_(beta=0)^(number of quadrature points)
        ///////////////////////////////////////////////////////////////////////
    }

    if ( CalculateStiffnessMatrixFlag == true )
    {
        AssembleTimeSpaceStiffnessFromStiffSubMatrices( rLeftHandSideMatrix, Help_K_PP, Help_K_PT, Help_K_TP, Help_K_TT );
// 		KRATOS_WATCH(rLeftHandSideMatrix);

    }

    if ( CalculateResidualVectorFlag == true )
    {
        AssembleTimeSpaceRHSFromSubVectors( rRightHandSideVector, Help_R_P, Help_R_T );
    }

    KRATOS_CATCH( "" )
}

///done
void Soil2PhaseRigid::DampMatrix( MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

// 	mMatSize= mNodesNumber;

    //resizing as needed the damping matrix
    if ( rDampMatrix.size1() != mMatSize )
        rDampMatrix.resize( mMatSize, mMatSize );
    noalias( rDampMatrix ) = ZeroMatrix( mMatSize, mMatSize ); //resetting LHS

    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );
    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );
    const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );

    Matrix Help_D_PP( mNodesNumber, mNodesNumber );
    Matrix Help_D_PT( mNodesNumber, mNodesNumber );
    Matrix Help_D_TP( mNodesNumber, mNodesNumber );
    Matrix Help_D_TT( mNodesNumber, mNodesNumber );

    noalias( Help_D_PP ) = ZeroMatrix( mNodesNumber, mNodesNumber );
    noalias( Help_D_PT ) = ZeroMatrix( mNodesNumber, mNodesNumber );
    noalias( Help_D_TP ) = ZeroMatrix( mNodesNumber, mNodesNumber );
    noalias( Help_D_TT ) = ZeroMatrix( mNodesNumber, mNodesNumber );


    Matrix DN_DX( mNodesNumber, mDimension );
    double weight;
    Vector N( mNodesNumber );

    double p;
    double t;
    double DetJ = 0.0;
    
    /////////////////////////////////////////////////////////////////////////
    //// Integration in space sum_(beta=0)^(number of quadrature points)
    /////////////////////////////////////////////////////////////////////////
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        if ( DN_DX.size1() != ( mNodesNumber ) || DN_DX.size2() != mDimension )
            DN_DX.resize( mNodesNumber, mDimension );
        noalias( DN_DX ) = ZeroMatrix( mNodesNumber, mDimension );
        noalias( DN_DX ) = prod( DN_De[PointNumber], mInvJ0[PointNumber] );

        weight = integration_points[PointNumber].Weight();

        // Jacobian on current quadrature point
        DetJ = mDetJ0[PointNumber];

        // Shape Functions on current spatial quadrature point
        if ( N.size() != mNodesNumber )
            N.resize( mNodesNumber );
        noalias( N ) = row( Ncontainer, PointNumber );

        p = Getp( N );
        t = Gett( N ); 
	

        // Calculation of spatial damping and Mass Matrix
        CalculateDampingMatrixPP( Help_D_PP, N, weight, DetJ, p, t );
        CalculateDampingMatrixPT( Help_D_PT, N, weight, DetJ, p, t );
        CalculateDampingMatrixTP( Help_D_TP, N, weight, DetJ, p, t );
        CalculateDampingMatrixTT( Help_D_TT, DN_DX, N, weight, DetJ, p, t );
//		CalculateMassMatrix(HelpMassMatrix, N, weight,DetJ,density);

        ///////////////////////////////////////////////////////////////////////
        // END Integration in space sum_(beta=0)^(number of quadrature points)
        ///////////////////////////////////////////////////////////////////////

    }
    AssembleTimeSpaceStiffnessFromDampSubMatrices( rDampMatrix, Help_D_PP, Help_D_PT, Help_D_TP, Help_D_TT );
// 	KRATOS_WATCH(rDampMatrix);

    KRATOS_CATCH( "" )
}
///done
void Soil2PhaseRigid::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag,  CalculateResidualVectorFlag );
}
///done
void Soil2PhaseRigid::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;

    CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );

}

///done
void Soil2PhaseRigid::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo )
{
    unsigned int dim_variables = 2;
    unsigned int index;

    if ( rResult.size() != mMatSize )
        rResult.resize( mMatSize );

    for ( unsigned int i = ( mNodesMin - 1 ); i < mNodesMax; i++ )
    {
        index = i * dim_variables;
        rResult[index] = GetGeometry()[i].GetDof( WATER_PRESSURE ).EquationId();
        rResult[index+1] = GetGeometry()[i].GetDof( TEMPERATURE ).EquationId();
    }
}
///done
void Soil2PhaseRigid::GetDofList( DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo )
{
    ElementalDofList.resize( 0 );
    for ( unsigned int i = ( mNodesMin - 1 ); i < mNodesMax; i++ )
    {
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( WATER_PRESSURE ) );
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( TEMPERATURE ) );
    }
}
///done
void Soil2PhaseRigid::GetValuesVector( Vector& values, int Step )
{
    unsigned int dim_variables = 2;
    unsigned int index;

    if ( values.size() != mMatSize )
        values.resize( mMatSize );

    for ( unsigned int i = ( mNodesMin - 1 );i < mNodesMax;i++ )
    {
        index = i * dim_variables;
        values( index ) = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE, Step );
        values( index + 1 ) = GetGeometry()[i].GetSolutionStepValue( TEMPERATURE, Step );
    }
}

///done
void Soil2PhaseRigid::CalculateOnIntegrationPoints( const Variable<double >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );
    const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );
    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );


    if ( Output.size() != integration_points.size() )
        Output.resize( integration_points.size() );

    Vector N( mNodesNumber );
    Matrix DN_DX( mNodesNumber, mDimension );

    double p;
    double t;
    Vector vL( mDimension );
    Vector qsl( mDimension );


    /////////////////////////////////////////////////////////////////////////
    //// Integration in space sum_(beta=0)^(number of quadrature points)
    /////////////////////////////////////////////////////////////////////////
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        if ( N.size() != mNodesNumber )
            N.resize( mNodesNumber );
        noalias( N ) = row( Ncontainer, PointNumber );

        if ( DN_DX.size1() != ( mNodesNumber ) || DN_DX.size2() != mDimension )
            DN_DX.resize( mNodesNumber, mDimension );
        noalias( DN_DX ) = ZeroMatrix( mNodesNumber, mDimension );
        noalias( DN_DX ) = prod( DN_De[PointNumber], mInvJ0[PointNumber] );

        p = Getp( N );
        t = Gett( N );
        noalias( vL ) = GetvL( DN_DX, p, t );
        noalias( qsl ) = Getq( DN_DX );

        if ( rVariable == WATER_PRESSURE )
        {
            Output[PointNumber] = p;
        }
        if ( rVariable == TEMPERATURE )
        {
            Output[PointNumber] = t;
        }
    }

    KRATOS_CATCH( "" )
}
///done
Soil2PhaseRigid::IntegrationMethod Soil2PhaseRigid::GetIntegrationMethod()
{
    return mThisIntegrationMethod;
}

///-//////////////////////////////////////////////////////////////////////////////////////
///-//////////////////////////////////////////////////////////////////////////////////////
///  PRIVATE OPERATIONS
///-//////////////////////////////////////////////////////////////////////////////////////
///-//////////////////////////////////////////////////////////////////////////////////////

///+++++++++++++++++++++++++++++++++++++++
/// CALCULATE INTERNAL FORCE VECTORS
///+++++++++++++++++++++++++++++++++++++++
///done
void Soil2PhaseRigid::AddInternalForcesToRHSP( Vector& Help_R_P, const Matrix& DN_DX, const Vector& N, double weight, double DetJ, double p, double t )
{
    //p
    Vector grad_p( mDimension );
    noalias( grad_p ) = Getgradp( DN_DX );
    double dot_p = Getdotp( N );

    //t
    Vector grad_t( mDimension );
    noalias( grad_t ) = Getgradt( DN_DX );
    double dot_t = Getdott( N );

    //rhoL
    double rhoL = GetrhoL( p, t );
    double DrhoL_Dp = GetDrhoLDp( 1, p, t );
    double DrhoL_Dt = GetDrhoLDt( 1, p, t );

    //vL
    Vector vL( mDimension );
    noalias( vL ) = GetvL( DN_DX, p, t );

    for ( unsigned int prim = 0; prim < mNodesNumber; prim++ )
    {
        //R1 P1.
        Help_R_P( prim ) += N( prim )
                            * mPorosity / rhoL
                            * ( DrhoL_Dp * dot_p
                                + DrhoL_Dt * dot_t )
                            * weight * DetJ * mScale;


        for ( unsigned int gamma = 0; gamma < mDimension; gamma++ )
        {
            //P2.
            Help_R_P( prim ) += N( prim )
                                / rhoL
                                * ( DrhoL_Dp * grad_p( gamma )
                                    + DrhoL_Dt * grad_t( gamma ) )
                                * vL( gamma )
                                * weight * DetJ * mScale;

            //P3.
            Help_R_P( prim ) += ( -1 ) * DN_DX( prim, gamma )
                                * vL( gamma )
                                * weight * DetJ * mScale;

        }
    }
}
///
void Soil2PhaseRigid::AddInternalForcesToRHST( Vector& Help_R_T, const Matrix& DN_DX, const Vector& N, double weight, double DetJ, double p, double t )
{
    //p
    double p0 = 101325 / mUnitRatio; // kPa
    Vector grad_p( mDimension );
    noalias( grad_p ) = Getgradp( DN_DX );
    double dot_p = Getdotp( N );

    //t
    Vector grad_t( mDimension );
    noalias( grad_t ) = Getgradt( DN_DX );
    double dot_t = Getdott( N );

    //density
    double rhoS0 = 2300 / pow(mUnitRatio, 3.0); //mMaterialParameters[8+0];
    double rhoL = GetrhoL( p, t );
    double DrhoL_Dp = GetDrhoLDp( 1, p, t );
    double DrhoL_Dt = GetDrhoLDt( 1, p, t );

    //vL
    Vector vL( mDimension );
    noalias( vL ) = GetvL( DN_DX, p, t );

    //internalEnergy
    double DeS_Dt = GetDeSDt( 1 );
    double DeL_Dp = GetDeLDp( N, 1, p, t );
    double DeL_Dt = GetDeLDt( N, 1, p, t );

    //heatFlux
    Vector q( mDimension );
    noalias( q ) = Getq( DN_DX );

    //supply of momentum
    Vector hatpL( mDimension );
    noalias( hatpL ) = GethatpL( DN_DX, p, t );

// std::cout<<"t= "<<t<<", p= "<<p<<", vL= "<<vL(0)<<", DeL_Dt= "<<DeL_Dt<<", q= "<<q(0)<<", hatpL= "<<hatpL(0)<<std::endl;

    for ( unsigned int prim = 0; prim < mNodesNumber; prim++ )
    {
        //R2 T1.
        Help_R_T( prim ) += N( prim )
                            * ( 1 - mPorosity ) * rhoS0
                            * ( DeS_Dt * dot_t )
                            * weight * DetJ * mScale;
        //T2.
        Help_R_T( prim ) += N( prim )
                            * mPorosity * rhoL
                            * ( DeL_Dp * dot_p
                                + DeL_Dt * dot_t )
                            * weight * DetJ * mScale;
        //T4.
        Help_R_T( prim ) += ( -1 ) * N( prim )
                            * ( p0 + p ) / rhoL
                            * ( DrhoL_Dp * dot_p
                                + DrhoL_Dt * dot_t )
                            * weight * DetJ * mScale;


        for ( unsigned int gamma = 0; gamma < mDimension; gamma++ )
        {
            //T3.
            Help_R_T( prim ) += N( prim )
                                * rhoL
                                * ( DeL_Dp * grad_p( gamma )
                                    + DeL_Dt * grad_t( gamma ) )
                                * vL( gamma )
                                * weight * DetJ * mScale;
            //T5.
            Help_R_T( prim ) += ( -1 ) * N( prim )
                                * ( p0 + p ) / ( mPorosity * rhoL )
                                * ( DrhoL_Dp * grad_p( gamma )
                                    + DrhoL_Dt * grad_t( gamma ) )
                                * vL( gamma )
                                * weight * DetJ * mScale;
            //T6.
            Help_R_T( prim ) += ( -1 ) * DN_DX( prim, gamma )
                                * ( q( gamma ) )
                                * weight * DetJ * mScale;
            //T7
            Help_R_T( prim ) += N( prim )
                                / mPorosity
                                * hatpL( gamma )
                                * vL( gamma )
                                * weight * DetJ * mScale;

        }
    }
}
///+++++++++++++++++++++++++++++++++++++++
/// CALCULATE STIFFNESS MATRICES
///+++++++++++++++++++++++++++++++++++++++
///done
void Soil2PhaseRigid::CalculateStiffnessMatrixPP( Matrix& Help_K_PP, const Matrix& DN_DX, const Vector& N, double weight, double DetJ, double p, double t )
{

    //p
    Vector grad_p( mDimension );
    noalias( grad_p ) = Getgradp( DN_DX );
    double dot_p = Getdotp( N );

    //t
    Vector grad_t( mDimension );
    noalias( grad_t ) = Getgradt( DN_DX );
    double dot_t = Getdott( N );

    //density
    double rhoL = GetrhoL( p, t );
    double DrhoL_Dp = GetDrhoLDp( 1, p, t );
    double D2rhoL_Dp2 = GetDrhoLDp( 2, p, t );
    double DrhoL_Dt = GetDrhoLDt( 1, p, t );
    double D2rhoL_DpDt = GetD2rhoLDpDt( p, t );

    //vL
    Vector vL( mDimension );
    noalias( vL ) = GetvL( DN_DX, p, t );
    Vector DvL_Dp( mDimension );
    noalias( DvL_Dp ) = GetDvLDp( DN_DX, p, t );
    double DvL_Dgradp = GetDvLDgradp( p, t );

    //stiffness matrix
    for ( unsigned int prim = 0; prim < mNodesNumber; prim++ )
    {
        for ( unsigned int sec = 0; sec < mNodesNumber; sec++ )
        {
            //PP1.
            Help_K_PP( prim, sec ) += ( -1 ) * N( prim )
                                      * mPorosity / pow( rhoL, 2 )
                                      * DrhoL_Dp
                                      * ( DrhoL_Dp * dot_p
                                          + DrhoL_Dt * dot_t )
                                      * N( sec )
                                      * weight * DetJ * mScale;

            Help_K_PP( prim, sec ) += N( prim )
                                      * mPorosity / rhoL
                                      * ( D2rhoL_Dp2 * dot_p
                                          + D2rhoL_DpDt * dot_t )
                                      * N( sec )
                                      * weight * DetJ * mScale;


            for ( unsigned int gamma = 0; gamma < mDimension; gamma++ )
            {
                //PP2.
                Help_K_PP( prim, sec ) += ( -1 ) * N( prim )
                                          / pow( rhoL, 2 )
                                          * DrhoL_Dp
                                          * ( DrhoL_Dp * grad_p( gamma )
                                              + DrhoL_Dt * grad_t( gamma ) )
                                          * vL( gamma )
                                          * N( sec )
                                          * weight * DetJ * mScale;

                Help_K_PP( prim, sec ) += N( prim )
                                          / rhoL
                                          * ( D2rhoL_Dp2 * grad_p( gamma )
                                              + D2rhoL_DpDt * grad_t( gamma ) )
                                          * vL( gamma )
                                          * N( sec )
                                          * weight * DetJ * mScale;

                Help_K_PP( prim, sec ) += N( prim )
                                          / rhoL
                                          * ( DrhoL_Dp )
                                          * vL( gamma )
                                          * DN_DX( sec, gamma )
                                          * weight * DetJ * mScale;

                Help_K_PP( prim, sec ) += N( prim )
                                          / rhoL
                                          * ( DrhoL_Dp * grad_p( gamma )
                                              + DrhoL_Dt * grad_t( gamma ) )
                                          * ( DvL_Dp( gamma ) * N( sec )
                                              + DvL_Dgradp * DN_DX( sec, gamma ) )
                                          * weight * DetJ * mScale;
                //PP3.
                Help_K_PP( prim, sec ) += ( -1 ) * DN_DX( prim, gamma )
                                          * ( DvL_Dp( gamma ) * N( sec )
                                              + DvL_Dgradp * DN_DX( sec, gamma ) )
                                          * weight * DetJ * mScale;
            }
        }
    }
}
///done
void Soil2PhaseRigid::CalculateStiffnessMatrixPT( Matrix& Help_K_PT, const Matrix& DN_DX, const Vector& N, double weight, double DetJ, double p, double t )
{
    //p
    Vector grad_p( mDimension );
    noalias( grad_p ) = Getgradp( DN_DX );
    double dot_p = Getdotp( N );

    //t
    Vector grad_t( mDimension );
    noalias( grad_t ) = Getgradt( DN_DX );
    double dot_t = Getdott( N );

    //density
    double rhoL = GetrhoL( p, t );
    double DrhoL_Dp = GetDrhoLDp( 1, p, t );
    double DrhoL_Dt = GetDrhoLDt( 1, p, t );
    double D2rhoL_Dt2 = GetDrhoLDt( 2, p, t );
    double D2rhoL_DpDt = GetD2rhoLDpDt( p, t );

    //vL
    Vector vL( mDimension );
    noalias( vL ) = GetvL( DN_DX, p, t );
    Vector DvL_Dt( mDimension );
    noalias( DvL_Dt ) = GetDvLDt( DN_DX, p, t );


    //stiffness matrix
    for ( unsigned int prim = 0; prim < mNodesNumber; prim++ )
    {
        for ( unsigned int sec = 0; sec < mNodesNumber; sec++ )
        {
            //PT1.
            Help_K_PT( prim, sec ) += ( -1 ) * N( prim )
                                      * mPorosity / pow( rhoL, 2 )
                                      * DrhoL_Dt
                                      * ( DrhoL_Dp * dot_p
                                          + DrhoL_Dt * dot_t )
                                      * N( sec )
                                      * weight * DetJ * mScale;

            Help_K_PT( prim, sec ) += N( prim )
                                      * mPorosity / rhoL
                                      * ( D2rhoL_DpDt * dot_p
                                          + D2rhoL_Dt2 * dot_t )
                                      * N( sec )
                                      * weight * DetJ * mScale;


            for ( unsigned int gamma = 0; gamma < mDimension; gamma++ )
            {
                //PT2.
                Help_K_PT( prim, sec ) += ( -1 ) * N( prim )
                                          / pow( rhoL, 2 )
                                          * DrhoL_Dt
                                          * ( DrhoL_Dp * grad_p( gamma )
                                              + DrhoL_Dt * grad_t( gamma ) )
                                          * vL( gamma )
                                          * N( sec )
                                          * weight * DetJ * mScale;

                Help_K_PT( prim, sec ) += N( prim )
                                          / rhoL
                                          * ( D2rhoL_DpDt * grad_p( gamma )
                                              + D2rhoL_Dt2 * grad_t( gamma ) )
                                          * vL( gamma )
                                          * N( sec )
                                          * weight * DetJ * mScale;

                Help_K_PT( prim, sec ) += N( prim )
                                          / rhoL
                                          * ( DrhoL_Dt )
                                          * vL( gamma )
                                          * DN_DX( sec, gamma )
                                          * weight * DetJ * mScale;

                Help_K_PT( prim, sec ) += N( prim )
                                          / rhoL
                                          * ( DrhoL_Dp * grad_p( gamma )
                                              + DrhoL_Dt * grad_t( gamma ) )
                                          * DvL_Dt( gamma )
                                          * N( sec )
                                          * weight * DetJ * mScale;
                //PT3.
                Help_K_PT( prim, sec ) += ( -1 ) * DN_DX( prim, gamma )
                                          * DvL_Dt( gamma )
                                          * N( sec )
                                          * weight * DetJ * mScale;
            }
        }
    }


}
///done
void Soil2PhaseRigid::CalculateStiffnessMatrixTP( Matrix& Help_K_TP, const Matrix& DN_DX, const Vector& N, double weight, double DetJ, double p, double t )
{
    //p
    double p0 = 101325 / mUnitRatio; // kPa
    Vector grad_p( mDimension );
    noalias( grad_p ) = Getgradp( DN_DX );
    double dot_p = Getdotp( N );

    //t
    Vector grad_t( mDimension );
    noalias( grad_t ) = Getgradt( DN_DX );
    double dot_t = Getdott( N );

    //density
    double rhoL = GetrhoL( p, t );
    double DrhoL_Dp = GetDrhoLDp( 1, p, t );
    double D2rhoL_Dp2 = GetDrhoLDp( 2, p, t );
    double DrhoL_Dt = GetDrhoLDt( 1, p, t );
    double D2rhoL_DpDt = GetD2rhoLDpDt( p, t );

    //vL
    Vector vL( mDimension );
    noalias( vL ) = GetvL( DN_DX, p, t );
    double DvL_Dgradp = GetDvLDgradp( p, t );
    Vector DvL_Dp( mDimension );
    noalias( DvL_Dp ) = GetDvLDp( DN_DX, p, t );

    //internalEnergy
    double DeL_Dp = GetDeLDp( N, 1, p, t );
    double D2wEpsilon_Dp2 = GetDeLDp( N, 2, p, t );
    double DeL_Dt = GetDeLDt( N, 1, p, t );
    double D2wEpsilon_DpDt = GetD2eLDpDt( N, p, t );

    //heatFlux
// 	Vector q(mDimension);
// 	noalias(q)= Getq(DN_DX);

    //supply of momentum
    Vector hatpL( mDimension );
    noalias( hatpL ) = GethatpL( DN_DX, p, t );
    double DhatpL_DvL = GetDhatpLDvL();

    //stiffness matrix
    for ( unsigned int prim = 0; prim < mNodesNumber; prim++ )
    {
        for ( unsigned int sec = 0; sec < mNodesNumber; sec++ )
        {
            //TP2.
            Help_K_TP( prim, sec ) += N( prim )
                                      * mPorosity
                                      * DrhoL_Dp
                                      * ( DeL_Dp * dot_p
                                          + DeL_Dt * dot_t )
                                      * N( sec )
                                      * weight * DetJ * mScale;

            Help_K_TP( prim, sec ) += N( prim )
                                      * mPorosity * rhoL
                                      * ( D2wEpsilon_Dp2 * dot_p
                                          + D2wEpsilon_DpDt * dot_t )
                                      * N( sec )
                                      * weight * DetJ * mScale;
            //TP4.
            Help_K_TP( prim, sec ) += ( -1 ) * N( prim )
                                      / rhoL
                                      * ( DrhoL_Dp * dot_p
                                          + DrhoL_Dt * dot_t )
                                      * N( sec )
                                      * weight * DetJ * mScale;

            Help_K_TP( prim, sec ) += N( prim )
                                      * ( p0 + p ) / pow( rhoL, 2 )
                                      * DrhoL_Dp
                                      * ( DrhoL_Dp * dot_p
                                          + DrhoL_Dt * dot_t )
                                      * N( sec )
                                      * weight * DetJ * mScale;

            Help_K_TP( prim, sec ) += ( -1 ) * N( prim )
                                      * ( p0 + p ) / rhoL
                                      * ( D2rhoL_Dp2 * dot_p
                                          + D2rhoL_DpDt * dot_t )
                                      * N( sec )
                                      * weight * DetJ * mScale;


            for ( unsigned int gamma = 0; gamma < mDimension; gamma++ )
            {
                //TP3.
                Help_K_TP( prim, sec ) += N( prim )
                                          * DrhoL_Dp
                                          * ( DeL_Dp * grad_p( gamma )
                                              + DeL_Dt * grad_t( gamma ) )
                                          * vL( gamma )
                                          * N( sec )
                                          * weight * DetJ * mScale;

                Help_K_TP( prim, sec ) += N( prim )
                                          * rhoL
                                          * ( D2wEpsilon_Dp2 * grad_p( gamma )
                                              + D2wEpsilon_DpDt * grad_t( gamma ) )
                                          * vL( gamma )
                                          * N( sec )
                                          * weight * DetJ * mScale;

                Help_K_TP( prim, sec ) += N( prim )
                                          * rhoL
                                          * ( DeL_Dp )
                                          * vL( gamma )
                                          * DN_DX( sec, gamma )
                                          * weight * DetJ * mScale;

                Help_K_TP( prim, sec ) += N( prim )
                                          * rhoL
                                          * ( DeL_Dp * grad_p( gamma )
                                              + DeL_Dt * grad_t( gamma ) )
                                          * ( DvL_Dp( gamma ) * N( sec )
                                              + DvL_Dgradp * DN_DX( sec, gamma ) )
                                          * weight * DetJ * mScale;
                //TP5.
                Help_K_TP( prim, sec ) += ( -1 ) * N( prim )
                                          / ( mPorosity * rhoL )
                                          * ( DrhoL_Dp * grad_p( gamma )
                                              + DrhoL_Dt * grad_t( gamma ) )
                                          * vL( gamma )
                                          * N( sec )
                                          * weight * DetJ * mScale;

                Help_K_TP( prim, sec ) += N( prim )
                                          * ( p0 + p ) / ( mPorosity * pow( rhoL, 2 ) )
                                          * DrhoL_Dp
                                          * ( DrhoL_Dp * grad_p( gamma )
                                              + DrhoL_Dt * grad_t( gamma ) )
                                          * vL( gamma )
                                          * N( sec )
                                          * weight * DetJ * mScale;

                Help_K_TP( prim, sec ) += ( -1 ) * N( prim )
                                          * ( p0 + p ) / ( mPorosity * rhoL )
                                          * ( D2rhoL_Dp2 * grad_p( gamma )
                                              + D2rhoL_DpDt * grad_t( gamma ) )
                                          * vL( gamma )
                                          * N( sec )
                                          * weight * DetJ * mScale;

                Help_K_TP( prim, sec ) += ( -1 ) * N( prim )
                                          * ( p0 + p ) / ( mPorosity * rhoL )
                                          * ( DrhoL_Dp )
                                          * vL( gamma )
                                          * DN_DX( sec, gamma )
                                          * weight * DetJ * mScale;

                Help_K_TP( prim, sec ) += ( -1 ) * N( prim )
                                          * ( p0 + p ) / ( mPorosity * rhoL )
                                          * ( DrhoL_Dp * grad_p( gamma )
                                              + DrhoL_Dt * grad_t( gamma ) )
                                          * ( DvL_Dp( gamma ) * N( sec )
                                              + DvL_Dgradp * DN_DX( sec, gamma ) )
                                          * weight * DetJ * mScale;
                //TP7.
                Help_K_TP( prim, sec ) += N( prim )
                                          / mPorosity
                                          * ( DhatpL_DvL * vL( gamma ) + hatpL( gamma ) )
                                          * ( DvL_Dp( gamma ) * N( sec )
                                              + DvL_Dgradp * DN_DX( sec, gamma ) )
                                          * weight * DetJ * mScale;

            }
        }
    }
}
///done
void Soil2PhaseRigid::CalculateStiffnessMatrixTT( Matrix& Help_K_TT, const Matrix& DN_DX, const Vector& N, double weight, double DetJ, double p, double t )
{
    //p
    double p0 = 101325 / mUnitRatio; // kPa
    Vector grad_p( mDimension );
    noalias( grad_p ) = Getgradp( DN_DX );
    double dot_p = Getdotp( N );

    //t
    Vector grad_t( mDimension );
    noalias( grad_t ) = Getgradt( DN_DX );
    double dot_t = Getdott( N );

    //density
    double rhoS0 = 2300 / pow(mUnitRatio, 3.0); //mMaterialParameters[8+0];
    double rhoL = GetrhoL( p, t );
    double DrhoL_Dp = GetDrhoLDp( 1, p, t );
    double DrhoL_Dt = GetDrhoLDt( 1, p, t );
    double D2rhoL_Dt2 = GetDrhoLDt( 2, p, t );
    double D2rhoL_DpDt = GetD2rhoLDpDt( p, t );

    //vL
    Vector vL( mDimension );
    noalias( vL ) = GetvL( DN_DX, p, t );
    Vector DvL_Dt( mDimension );
    noalias( DvL_Dt ) = GetDvLDt( DN_DX, p, t );

    //internalEnergy
    double D2sEpsilon_Dt2 = GetDeSDt( 2 );
    double DeL_Dp = GetDeLDp( N, 1, p, t );
    double DeL_Dt = GetDeLDt( N, 1, p, t );
    double D2wEpsilon_Dt2 = GetDeLDt( N, 2, p, t );
    double D2wEpsilon_DpDt = GetD2eLDpDt( N, p, t );

    //heatFlux
    Vector q( mDimension );
    noalias( q ) = Getq( DN_DX );
    double Dq_Dgradt = GetDqDgradt();

    //supply of momentum
    Vector hatpL( mDimension );
    noalias( hatpL ) = GethatpL( DN_DX, p, t );
    double DhatpL_DvL = GetDhatpLDvL();

    //stiffness matrix
    for ( unsigned int prim = 0; prim < mNodesNumber; prim++ )
    {
        for ( unsigned int sec = 0; sec < mNodesNumber; sec++ )
        {
            //TT1.
            Help_K_TT( prim, sec ) += N( prim )
                                      * ( 1 - mPorosity )
                                      * rhoS0
                                      * ( D2sEpsilon_Dt2 * dot_t )
                                      * N( sec )
                                      * weight * DetJ * mScale;

            //TT2.
            Help_K_TT( prim, sec ) += N( prim )
                                      * mPorosity
                                      * DrhoL_Dt
                                      * ( DeL_Dp * dot_p
                                          + DeL_Dt * dot_t )
                                      * N( sec )
                                      * weight * DetJ * mScale;

            Help_K_TT( prim, sec ) += N( prim )
                                      * mPorosity * rhoL
                                      * ( D2wEpsilon_DpDt * dot_p
                                          + D2wEpsilon_Dt2 * dot_t )
                                      * N( sec )
                                      * weight * DetJ * mScale;
            //TT4.
            Help_K_TT( prim, sec ) += N( prim )
                                      * ( p0 + p ) / pow( rhoL, 2 )
                                      * DrhoL_Dt
                                      * ( DrhoL_Dp * dot_p
                                          + DrhoL_Dt * dot_t )
                                      * N( sec )
                                      * weight * DetJ * mScale;

            Help_K_TT( prim, sec ) += ( -1 ) * N( prim )
                                      * ( p0 + p ) / rhoL
                                      * ( D2rhoL_DpDt * dot_p
                                          + D2rhoL_Dt2 * dot_t )
                                      * N( sec )
                                      * weight * DetJ * mScale;


            for ( unsigned int gamma = 0; gamma < mDimension; gamma++ )
            {
                //TT3.
                Help_K_TT( prim, sec ) += N( prim )
                                          * DrhoL_Dt
                                          * ( DeL_Dp * grad_p( gamma )
                                              + DeL_Dt * grad_t( gamma ) )
                                          * vL( gamma )
                                          * N( sec )
                                          * weight * DetJ * mScale;

                Help_K_TT( prim, sec ) += N( prim )
                                          * rhoL
                                          * ( D2wEpsilon_DpDt * grad_p( gamma )
                                              + D2wEpsilon_Dt2 * grad_t( gamma ) )
                                          * vL( gamma )
                                          * N( sec )
                                          * weight * DetJ * mScale;

                Help_K_TT( prim, sec ) += N( prim )
                                          * rhoL
                                          * ( DeL_Dt )
                                          * vL( gamma )
                                          * DN_DX( sec, gamma )
                                          * weight * DetJ * mScale;

                Help_K_TT( prim, sec ) += N( prim )
                                          * rhoL
                                          * ( DeL_Dp * grad_p( gamma )
                                              + DeL_Dt * grad_t( gamma ) )
                                          * DvL_Dt( gamma )
                                          * N( sec )
                                          * weight * DetJ * mScale;
                //TT5.
                Help_K_TT( prim, sec ) += N( prim )
                                          * ( p0 + p ) / ( mPorosity * pow( rhoL, 2 ) )
                                          * DrhoL_Dt
                                          * ( DrhoL_Dp * grad_p( gamma )
                                              + DrhoL_Dt * grad_t( gamma ) )
                                          * vL( gamma )
                                          * N( sec )
                                          * weight * DetJ * mScale;

                Help_K_TT( prim, sec ) += ( -1 ) * N( prim )
                                          * ( p0 + p ) / ( mPorosity * rhoL )
                                          * ( D2rhoL_DpDt * grad_p( gamma )
                                              + D2rhoL_Dt2 * grad_t( gamma ) )
                                          * vL( gamma )
                                          * N( sec )
                                          * weight * DetJ * mScale;

                Help_K_TT( prim, sec ) += ( -1 ) * N( prim )
                                          * ( p0 + p ) / ( mPorosity * rhoL )
                                          * ( DrhoL_Dt )
                                          * vL( gamma )
                                          * DN_DX( sec, gamma )
                                          * weight * DetJ * mScale;

                Help_K_TT( prim, sec ) += ( -1 ) * N( prim )
                                          * ( p0 + p ) / ( mPorosity * rhoL )
                                          * ( DrhoL_Dp * grad_p( gamma )
                                              + DrhoL_Dt * grad_t( gamma ) )
                                          * DvL_Dt( gamma )
                                          * N( sec )
                                          * weight * DetJ * mScale;
                //TT6.
                Help_K_TT( prim, sec ) += ( -1 ) * DN_DX( prim, gamma )
                                          * ( Dq_Dgradt )
                                          * DN_DX( sec, gamma )
                                          * weight * DetJ * mScale;
                //TP7.
                Help_K_TT( prim, sec ) += N( prim )
                                          / mPorosity
                                          * ( DhatpL_DvL * vL( gamma ) + hatpL( gamma ) )
                                          * DvL_Dt( gamma )
                                          * N( sec )
                                          * weight * DetJ * mScale;
            }
        }
    }
}
///+++++++++++++++++++++++++++++++++++++++
/// CALCULATE DAMPING MATRICES
///+++++++++++++++++++++++++++++++++++++++
///done
void Soil2PhaseRigid::CalculateDampingMatrixPP( Matrix& Help_D_PP, const Vector& N, double weight, double DetJ, double p, double t )
{
    //rhoL
    double rhoL = GetrhoL( p, t );
    double DrhoL_Dp = GetDrhoLDp( 1, p, t );

    for ( unsigned int prim = 0; prim < mNodesNumber; prim++ )
    {
        for ( unsigned int sec = 0; sec < mNodesNumber; sec++ )
        {
            //PP1.
            Help_D_PP( prim, sec ) +=
                N( prim )
                * mPorosity / rhoL
                * DrhoL_Dp
                * N( sec )
                * weight * DetJ * mScale;
        }
    }
}
///
void Soil2PhaseRigid::CalculateDampingMatrixPT( Matrix& Help_D_PT, const Vector& N, double weight, double DetJ, double p, double t )
{
    //rhoL
    double rhoL = GetrhoL( p, t );
    double DrhoL_Dt = GetDrhoLDt( 1, p, t );

    for ( unsigned int prim = 0; prim < mNodesNumber; prim++ )
    {
        for ( unsigned int sec = 0; sec < mNodesNumber; sec++ )
        {
            //PT1.
            Help_D_PT( prim, sec ) +=
                N( prim )
                * mPorosity / rhoL
                * DrhoL_Dt
                * N( sec )
                * weight * DetJ * mScale;
        }
    }
}
///
void Soil2PhaseRigid::CalculateDampingMatrixTP( Matrix& Help_D_TP, const Vector& N, double weight, double DetJ, double p, double t )
{
    //water Pressure
    double p0 = 101325 / mUnitRatio;  // kPa

    //rhoL
    double rhoL = GetrhoL( p, t );
    double DrhoL_Dp = GetDrhoLDp( 1, p, t );

    //internalEnergy
    double DeL_Dp = GetDeLDp( N, 1, p, t );

    for ( unsigned int prim = 0; prim < mNodesNumber; prim++ )
    {
        for ( unsigned int sec = 0; sec < mNodesNumber; sec++ )
        {
            //TP2.
            Help_D_TP( prim, sec ) +=
                N( prim )
                * mPorosity * rhoL
                * DeL_Dp
                * N( sec )
                * weight * DetJ * mScale;
            //TP4.
            Help_D_TP( prim, sec ) +=
                ( -1 ) * N( prim )
                * ( p0 + p ) / rhoL
                * DrhoL_Dp
                * N( sec )
                * weight * DetJ * mScale;
        }
    }
}
///
void Soil2PhaseRigid::CalculateDampingMatrixTT( Matrix& Help_D_TT, const Matrix& DN_DX, const Vector& N, double weight, double DetJ, double p, double t )
{
    double p0 = 101325 / mUnitRatio;  // kPa

    //density
    double rhoS0 = 2300 / pow(mUnitRatio, 3.0); //mMaterialParameters[8+0];
    double rhoL = GetrhoL( p, t );
    double DrhoL_Dt = GetDrhoLDt( 1, p, t );

    //vL
    Vector vL( mDimension );
    noalias( vL ) = GetvL( DN_DX, p, t );

    //internalEnergy
    double DeS_Dt = GetDeSDt( 1 );
    double DeL_Dt = GetDeLDt( N, 1, p, t );

    for ( unsigned int prim = 0; prim < mNodesNumber; prim++ )
    {
        for ( unsigned int sec = 0; sec < mNodesNumber; sec++ )
        {
            //TT1.
            Help_D_TT( prim, sec ) +=
                N( prim )
                * ( 1 - mPorosity ) * rhoS0
                * DeS_Dt
                * N( sec )
                * weight * DetJ * mScale;
            //TT2.
            Help_D_TT( prim, sec ) +=
                N( prim )
                * mPorosity * rhoL
                * DeL_Dt
                * N( sec )
                * weight * DetJ * mScale;
            //TT4.
            Help_D_TT( prim, sec ) +=
                ( -1 ) * N( prim )
                * ( p0 + p ) / rhoL
                * DrhoL_Dt
                * N( sec )
                * weight * DetJ * mScale;
        }
    }
}

///+++++++++++++++++++++++++++++++++++++++
/// Assemble
///+++++++++++++++++++++++++++++++++++++++
///done
void Soil2PhaseRigid::AssembleTimeSpaceRHSFromSubVectors( VectorType& rRightHandSideVector, const Vector& R_P, const Vector& R_T )
{
    KRATOS_TRY

    unsigned int dim_variables = 2;
    unsigned int index_time;

    for ( unsigned int prim = 0; prim < mNodesNumber; prim++ )
    {
        index_time = prim * dim_variables;

        rRightHandSideVector( index_time ) -= R_P( prim );
    }

    for ( unsigned int prim = 0; prim < mNodesNumber; prim++ )
    {
        index_time = prim * dim_variables + 1;

        rRightHandSideVector( index_time ) -= R_T( prim );
    }
    KRATOS_CATCH( "" )

}
///done
void Soil2PhaseRigid::AssembleTimeSpaceStiffnessFromStiffSubMatrices( MatrixType& rLeftHandSideMatrix, const Matrix& K_PP, const Matrix& K_PT, const Matrix& K_TP, const Matrix& K_TT )
{
    KRATOS_TRY
    unsigned int dim_variables = 2;
    unsigned int index_time_prim;
    unsigned int index_time_sec;

    for ( unsigned int prim = 0; prim < mNodesNumber; prim++ )
    {
        index_time_prim = prim * dim_variables;

        for ( unsigned int sec = 0; sec < mNodesNumber; sec++ )
        {
            index_time_sec = sec * dim_variables;

            rLeftHandSideMatrix( index_time_prim, index_time_sec )
            += K_PP( prim, sec );
        }
    }

    for ( unsigned int prim = 0; prim < mNodesNumber; prim++ )
    {
        index_time_prim = prim * dim_variables;

        for ( unsigned int sec = 0; sec < mNodesNumber; sec++ )
        {
            index_time_sec = sec * dim_variables + 1;

            rLeftHandSideMatrix( index_time_prim, index_time_sec )
            += K_PT( prim, sec );
        }
    }

    for ( unsigned int prim = 0; prim < mNodesNumber; prim++ )
    {
        index_time_prim = prim * dim_variables + 1;

        for ( unsigned int sec = 0; sec < mNodesNumber; sec++ )
        {
            index_time_sec = sec * dim_variables;

            rLeftHandSideMatrix( index_time_prim, index_time_sec )
            += K_TP( prim, sec );
        }
    }

    for ( unsigned int prim = 0; prim < mNodesNumber; prim++ )
    {
        index_time_prim = prim * dim_variables + 1;

        for ( unsigned int sec = 0; sec < mNodesNumber; sec++ )
        {
            index_time_sec = sec * dim_variables + 1;

            rLeftHandSideMatrix( index_time_prim, index_time_sec )
            += K_TT( prim, sec );
        }
    }
    KRATOS_CATCH( "" )
}
///done
void Soil2PhaseRigid::AssembleTimeSpaceStiffnessFromDampSubMatrices( MatrixType& rLeftHandSideMatrix, const Matrix& D_PP, const Matrix& D_PT, const Matrix& D_TP, const Matrix& D_TT )
{
    unsigned int dim_variables = 2;
    unsigned int index_time_prim;
    unsigned int index_time_sec;

    for ( unsigned int prim = 0; prim < mNodesNumber; prim++ )
    {
        index_time_prim = prim * dim_variables;

        for ( unsigned int sec = 0; sec < mNodesNumber; sec++ )
        {
            index_time_sec = sec * dim_variables;

            rLeftHandSideMatrix( index_time_prim, index_time_sec )
            += D_PP( prim, sec );
        }
    }

    for ( unsigned int prim = 0; prim < mNodesNumber; prim++ )
    {
        index_time_prim = prim * dim_variables;

        for ( unsigned int sec = 0; sec < mNodesNumber; sec++ )
        {
            index_time_sec = sec * dim_variables + 1;

            rLeftHandSideMatrix( index_time_prim, index_time_sec )
            += D_PT( prim, sec );
        }
    }

    for ( unsigned int prim = 0; prim < mNodesNumber; prim++ )
    {
        index_time_prim = prim * dim_variables + 1;

        for ( unsigned int sec = 0; sec < mNodesNumber; sec++ )
        {
            index_time_sec = sec * dim_variables;

            rLeftHandSideMatrix( index_time_prim, index_time_sec )
            += D_TP( prim, sec );
        }
    }

    for ( unsigned int prim = 0; prim < mNodesNumber; prim++ )
    {
        index_time_prim = prim * dim_variables + 1;

        for ( unsigned int sec = 0; sec < mNodesNumber; sec++ )
        {
            index_time_sec = sec * dim_variables + 1;

            rLeftHandSideMatrix( index_time_prim, index_time_sec )
            += D_TT( prim, sec );
        }
    }
}

///-//////////////////////////////////////////////////////////////////////////////////////
///-//////////////////////////////////////////////////////////////////////////////////////
///  PRIMARY VARIABLES AND THEIR DERIVATIVES
///-//////////////////////////////////////////////////////////////////////////////////////
///-//////////////////////////////////////////////////////////////////////////////////////
///--------------------------------------------------------------------------
/// PRESSURE (water)
///--------------------------------------------------------------------------
///done
double Soil2PhaseRigid::Getp( const Vector& N )
{
    double result = 0.0;

    for ( unsigned int i = mNodesMin - 1 ; i < mNodesMax ;i++ )
    {
        result += GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE ) * N( i );
    }
    return result;
}
///done
Vector Soil2PhaseRigid::Getgradp( const Matrix& DN_DX )
{
    Vector result( mDimension );
    noalias( result ) = ZeroVector( mDimension );

    double p_e;

    for ( unsigned int i = mNodesMin - 1 ; i < mNodesMax ;i++ )
    {

        p_e = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE );

        for ( unsigned int k = 0; k < 3; k++ )
        {
            result( k ) += p_e * DN_DX( i, k );
        }
    }

    return result;
}
///done
double Soil2PhaseRigid::Getdotp( const Vector& N )
{
    double result = 0.0;

    for ( unsigned int i = mNodesMin - 1 ; i < mNodesMax ;i++ )
    {
        result += GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE_DT ) * N( i );
    }

    return result;
}

///--------------------------------------------------------------------------
/// Temperature
///--------------------------------------------------------------------------
///done
double Soil2PhaseRigid::Gett( const Vector& N )
{
    double result = 0.0;

    for ( unsigned int i = mNodesMin - 1 ; i < mNodesMax ;i++ )
    {
        result += GetGeometry()[i].GetSolutionStepValue( TEMPERATURE ) * N( i );
    }
    return result;
}
///done
Vector Soil2PhaseRigid::Getgradt( const Matrix& DN_DX )
{
    Vector result( mDimension );
    noalias( result ) = ZeroVector( mDimension );

    double t_e;

    for ( unsigned int i = mNodesMin - 1 ; i < mNodesMax ;i++ )
    {

        t_e = GetGeometry()[i].GetSolutionStepValue( TEMPERATURE );

        for ( unsigned int k = 0; k < 3; k++ )
        {
            result( k ) += t_e * DN_DX( i, k );
        }
    }

    return result;
}
///done
double Soil2PhaseRigid::Getdott( const Vector& N )
{
    double result = 0.0;

    for ( unsigned int i = mNodesMin - 1 ; i < mNodesMax ;i++ )
    {
        result += GetGeometry()[i].GetSolutionStepValue( TEMPERATURE_DT ) * N( i );
    }

    return result;
}

///--------------------------------------------------------------------------
/// Density (water)
///--------------------------------------------------------------------------
///done
double Soil2PhaseRigid::GetrhoL( double p, double t )
{
    double result = 0.0;

//         if ( GetValue( PLASTIC_FLAG ) == 0 )
//         {
// 	double c1= -14596.0 / mUnitRatio;
// 	double c2= 61494.0 / mUnitRatio;
// 	double c3= -97570.0 / mUnitRatio;
// 	double c4= 69059.0 / mUnitRatio;
// 	double c5= -17387.0 / mUnitRatio;
// // 	double p0 = 101325;
// 	double t0= 273.15;
// 	double tau= (t + t0)/t0;
// 	double E=2.15e+9 / mUnitRatio;
// 
// 	// real water density for t from 243K to 303K
// 	result= (c1*pow(tau,4) + c2*pow(tau,3) + c3*pow(tau,2) + c4*tau + c5)
// 		/(1.0 - p/E);
// 	}
// 	else
// 	{
    double T = 273.15 + t;
    double M_w = /*GetProperties()[MOLAR_MASS_WATER];*/ 24.042835;
    double R = /*GetProperties()[GAS_CONSTANT];*/  8.314472 * pow(mUnitRatio, 2.0);

    result = 1000.0 / pow(mUnitRatio, 3.0) + p * M_w / ( R * T ) ;

// 	}
    return result;
}
  
///done
double Soil2PhaseRigid::GetDrhoLDp( int n, double p, double t )
{
    double result = 0.0;

//         if ( GetValue( PLASTIC_FLAG ) == 0 )
//         {
// 	double c1= -14596.0 / mUnitRatio;
// 	double c2= 61494.0 / mUnitRatio;
// 	double c3= -97570.0 / mUnitRatio;
// 	double c4= 69059.0 / mUnitRatio;
// 	double c5= -17387.0 / mUnitRatio;
// 	double t0= 273.15;
// 	double tau= (t + t0)/t0;
// 	double E=2.15e+9 / mUnitRatio;
// 
//     switch ( n )
//     {
//         case 1: result= ( c1*pow( tau,4 ) + c2*pow( tau,3 ) + c3*pow( tau,2 ) + c4*tau + c5 )
//                             / pow(( 1.0 - p / E ), 2 )
//                             / E;
//             break;
// 
//         case 2: result= ( c1*pow( tau,4 ) + c2*pow( tau,3 ) + c3*pow( tau,2 ) + c4*tau + c5 )
//                             / pow(( 1.0 - p / E ), 3 )
//                             * 2 / pow( E, 2 );
//     }
// 	}
// 	else
// 	{
    double T = 273.15 + t;
    double M_w = /*GetProperties()[MOLAR_MASS_WATER];*/ 24.042835;
    double R = /*GetProperties()[GAS_CONSTANT];*/  8.314472 * pow(mUnitRatio, 2.0);
    
    if(n==1)
    result= M_w / ( R * T ) ;
// 	}
    return result;
}
///done
double Soil2PhaseRigid::GetDrhoLDt( int n, double p, double t )
{
    double result = 0.0;

//         if ( GetValue( PLASTIC_FLAG ) == 0 )
//         {
// 	double c1= -14596.0 / mUnitRatio;
// 	double c2= 61494.0 / mUnitRatio;
// 	double c3= -97570.0 / mUnitRatio;
// 	double c4= 69059.0 / mUnitRatio;
// 	double c5= -17387.0 / mUnitRatio;
// 	double t0= 273.15;
// 	double tau= (t + t0)/t0;
// 	double E=2.15e+9 / mUnitRatio;
// 
//     switch ( n )
//     {
//         case 1: result= ( 4*c1*pow( tau,3 ) + 3*c2*pow( tau,2 ) + 2*c3*tau + c4 )
//                             / ( 1.0 - p / E )
//                             / t0;
//             break;
// 
//         case 2: result= ( 12*c1*pow( tau,2 ) + 6*c2*tau + 2*c3 )
//                             / ( 1.0 - p / E )
//                             / pow( t0, 2 );
//             break;
// 
//         case 3: result= ( 24*c1*tau + 6*c2 )
//                             / ( 1.0 - p / E )
//                             / pow( t0, 3 );
//             break;
//     }
// 	}
//     else
//     {
 
	double T = 273.15 + t;
	double M_w = /*GetProperties()[MOLAR_MASS_WATER];*/ 24.042835;
	double R = /*GetProperties()[GAS_CONSTANT];*/  8.314472 * pow(mUnitRatio, 2.0);

	switch(n)
	{
		case 1: result= (-1)*M_w/R*pow(T,-2)*p;
			break;

		case 2: result= (+2)*M_w/R*pow(T,-3)*p;
			break;

		case 3: result= (-6)*M_w/R*pow(T,-4)*p;
			break;
	}
//     }
 
    return result;
}
///done
double Soil2PhaseRigid::GetD2rhoLDpDt( double p, double t )
{
    double result = 0.0;

        if ( GetValue( PLASTIC_FLAG ) == 0 )
        {
	double c1= -14596.0 / mUnitRatio;
	double c2= 61494.0 / mUnitRatio;
	double c3= -97570.0 / mUnitRatio;
	double c4= 69059.0 / mUnitRatio;
	double t0= 273.15;
	double tau= (t + t0)/t0;
	double E=2.15e+9 / mUnitRatio;

    result = ( 4 * c1 * pow( tau, 3 ) + 3 * c2 * pow( tau, 2 ) + 2 * c3 * tau + c4 )
             / pow(( 1.0 - p / E ), 2 )
             / E
             / t0;
	}
	else
	{

	double T = 273.15 + t;
	double M_w = /*GetProperties()[MOLAR_MASS_WATER];*/ 24.042835;
	double R = /*GetProperties()[GAS_CONSTANT];*/  8.314472 * pow(mUnitRatio, 2.0);
	
	result= (-1)*M_w/R*pow(T,-2) ;
	}
    return result;
}
///TODO
double Soil2PhaseRigid::GetD3rhoLDp2Dt( double p, double t )
{
    double result = 0.0;

        if ( GetValue( PLASTIC_FLAG ) == 0 )
        {
	double c1= -14596.0 / mUnitRatio;
	double c2= 61494.0 / mUnitRatio;
	double c3= -97570.0 / mUnitRatio;
	double c4= 69059.0 / mUnitRatio;
	double t0= 273.15;
	double tau= (t + t0)/t0;
	double E=2.15e+9 / mUnitRatio;


    result = ( 4 * c1 * pow( tau, 3 ) + 3 * c2 * pow( tau, 2 ) + 2 * c3 * tau + c4 )
             / pow(( 1.0 - p / E ), 3 )
             * 2 / pow( E, 2 )
             / t0; 
	}
	
    return result;
}
///TODO
double Soil2PhaseRigid::GetD3rhoLDpDt2( double p, double t )
{
    double result = 0.0;


        if ( GetValue( PLASTIC_FLAG ) == 0 )
        {
	double c1= -14596.0 / mUnitRatio;
	double c2= 61494.0 / mUnitRatio;
	double c3= -97570.0 / mUnitRatio;
	double c4= 69059.0 / mUnitRatio;
	double t0= 273.15;
	double tau= (t + t0)/t0;
	double E=2.15e+9 / mUnitRatio;

    result = ( 4 * c1 * pow( tau, 3 ) + 3 * c2 * pow( tau, 2 ) + 2 * c3 * tau + c4 )
             / pow(( 1.0 - p / E ), 2 )
             / E
             / pow( t0, 2 );
	}
	else
	{
	double T = 273.15 + t;
	double M_w = /*GetProperties()[MOLAR_MASS_WATER];*/ 24.042835;
	double R = /*GetProperties()[GAS_CONSTANT];*/  8.314472 * pow(mUnitRatio, 2.0); // J/mol/K or Kg m^2 s-2/mol/K
	
	result= (2)*M_w/R*pow(T,-3) ;
	}
    return result;
}
///--------------------------------------------------------------------------
/// Darcy velocity  (water)
///--------------------------------------------------------------------------
///done
Vector Soil2PhaseRigid::GetvL( const Matrix& DN_DX, double p, double t )
{
    Vector result( mDimension );
    noalias( result ) = ZeroVector( mDimension );

    Vector grad_p( mDimension );
    noalias( grad_p ) = Getgradp( DN_DX );

    Vector gravity( mDimension );
    noalias( gravity ) = GetProperties()[GRAVITY];

    double rhoL = GetrhoL( p, t );
    double kappa0 = 4.4e-6 * mUnitRatio; //GetProperties()[PERMEABILITY_WATER];
//     double kappa0 = GetProperties()[PERMEABILITY_WATER] / mMaterialParameters[8+20];

    for ( unsigned int i = 0; i < mDimension; i++ )
    {
		result(i)= - kappa0 / (rhoL * 9.81 * mUnitRatio) * (grad_p(i) - rhoL * gravity(i));
//         result( i ) = - kappa0 * ( grad_p( i ) - rhoL * gravity( i ) );
    }

    return result;
}
///done
Vector Soil2PhaseRigid::GetDvLDrhoL( const Matrix& DN_DX, double p, double t )
{
    Vector result( mDimension );
    noalias( result ) = ZeroVector( mDimension );

    Vector gravity( mDimension );
    noalias( gravity ) = GetProperties()[GRAVITY];

	Vector grad_p(mDimension);
	noalias(grad_p)= Getgradp(DN_DX);

	double rhoL = GetrhoL(p, t);
    double kappa0 = 4.4e-6 * mUnitRatio; //GetProperties()[PERMEABILITY_WATER];
//     double kappa0 = GetProperties()[PERMEABILITY_WATER] / mMaterialParameters[8+20];

    for ( unsigned int i = 0; i < mDimension; i++ )
    {
		result(i)= kappa0 / (9.81 * mUnitRatio * pow(rhoL,2)) *grad_p(i);
//         result( i ) =   kappa0  * gravity( i );
    }

    return result;

}
///done
double Soil2PhaseRigid::GetDvLDgradp( double p, double t )
{
    double result = 0.0;
    double rhoL = GetrhoL( p, t );
    double kappa0 = 4.4e-6 * mUnitRatio; //GetProperties()[PERMEABILITY_WATER];
//     double kappa0 = GetProperties()[PERMEABILITY_WATER] / mMaterialParameters[8+20]

	result= - kappa0 / (rhoL * 9.81 * mUnitRatio);
//     result = - kappa0;

    return result;

}
///done
Vector Soil2PhaseRigid::GetDvLDp( const Matrix& DN_DX, double p, double t )
{
    Vector result( mDimension );
    noalias( result ) = ZeroVector( mDimension );

    Vector DvL_DrhoL( mDimension );
    noalias( DvL_DrhoL ) = GetDvLDrhoL( DN_DX, p, t );

    double DrhoL_Dp = GetDrhoLDp( 1, p, t );

    for ( unsigned int i = 0; i < mDimension; i++ )
    {
        result( i ) = DvL_DrhoL( i ) * DrhoL_Dp;
    }

    return result;

}
///done
Vector Soil2PhaseRigid::GetDvLDt( const Matrix& DN_DX, double p, double t )
{
    Vector result( mDimension );
    noalias( result ) = ZeroVector( mDimension );

    Vector DvL_DrhoL( mDimension );
    noalias( DvL_DrhoL ) = GetDvLDrhoL( DN_DX, p, t );

    double DrhoL_Dt = GetDrhoLDt( 1, p, t );

    for ( unsigned int i = 0; i < mDimension; i++ )
    {
        result( i ) = DvL_DrhoL( i ) * DrhoL_Dt;
    }

    return result;
}


///--------------------------------------------------------------------------
/// Free Helmholtz energy (water)
///--------------------------------------------------------------------------
///done
double Soil2PhaseRigid::GetDpsiLDt( const Vector& N, int n, double p, double t )
{
    double result = 0.0;

    double T0 = 273.15;
    for ( unsigned int i = mNodesMin - 1 ; i < mNodesMax ;i++ )
    {
        T0  += mInitialTemperature( i ) * N( i );
    }


    double T = 273.15 + t;
    double alphaL = 1.8e-6; //mMaterialParameters[8+16];
    double kL = 20e+9 / mUnitRatio; //mMaterialParameters[8+4];
    double cL = 4100 * pow(mUnitRatio, 2.0); //mMaterialParameters[8+10];

    double rhoL = GetrhoL( p, t );
    double DrhoL_Dt = GetDrhoLDt( 1, p, t );
    double D2rhoL_Dt2 = GetDrhoLDt( 2, p, t );
    double D3rhoL_Dt3 = GetDrhoLDt( 3, p, t );

    switch ( n )
    {
        case 1: result= 3*alphaL*kL/mPorosity
                            * ( pow( rhoL, -2 ) * DrhoL_Dt * ( T - T0 )
                                - pow( rhoL, -1 ) )
                            - cL * log( T / T0 );
            break;

        case 2: result= 3*alphaL*kL/mPorosity
                            * ( -2 * pow( rhoL, -3 ) * pow( DrhoL_Dt, 2 ) * ( T - T0 )
                                + pow( rhoL, -2 ) * D2rhoL_Dt2 * ( T - T0 )
                                + 2 * pow( rhoL, -2 ) * DrhoL_Dt )
                            - cL * pow( T, -1 );
            break;

        case 3: result= 3*alphaL*kL/mPorosity
                            * ( 6 * pow( rhoL, -4 ) * pow( DrhoL_Dt, 3 ) * ( T - T0 )
                                - 6 * pow( rhoL, -3 ) * D2rhoL_Dt2 * DrhoL_Dt * ( T - T0 )
                                - 6 * pow( rhoL, -3 ) * pow( DrhoL_Dt, 2 )
                                + pow( rhoL, -2 ) * D3rhoL_Dt3 * ( T - T0 )
                                + 3 * pow( rhoL, -2 ) * D2rhoL_Dt2 )
                            + cL * pow( T, -2 );
            break;
    }

    return result;
}
///done
double Soil2PhaseRigid::GetDpsiLDp( const Vector& N, int n, double p, double t )
{
    double result = 0.0;

    double T0 = 273.15;
    for ( unsigned int i = mNodesMin - 1 ; i < mNodesMax ;i++ )
    {
        T0  += mInitialTemperature( i ) * N( i );
    }

    double T = 273.15 + t;
    double alphaL = 1.8e-6; //mMaterialParameters[8+16];
    double kL = 20e+9 / mUnitRatio; //mMaterialParameters[8+4];

    double rhoL = GetrhoL( p, t );
    double DrhoL_Dp = GetDrhoLDp( 1, p, t );
    double D2rhoL_Dp2 = GetDrhoLDp( 2, p, t );

    switch ( n )
    {
        case 1: result= 3*alphaL*kL/mPorosity
                            * ( T - T0 )
                            * pow( rhoL, -2 ) * DrhoL_Dp ;
            break;

        case 2: result= 3*alphaL*kL/mPorosity
                            * ( T - T0 )
                            * (( -2 ) * pow( rhoL, -3 ) * pow( DrhoL_Dp, 2 )
                               + pow( rhoL, -2 ) * D2rhoL_Dp2 );
            break;
    }

    return result;
}
///
double Soil2PhaseRigid::GetD2psiLDpDt( const Vector& N, double p, double t )
{
    double result = 0.0;

    double T0 = 273.15;
    for ( unsigned int i = mNodesMin - 1 ; i < mNodesMax ;i++ )
    {
        T0  += mInitialTemperature( i ) * N( i );
    }

    double T = 273.15 + t;
    double alphaL = 1.8e-6; //mMaterialParameters[8+16];
    double kL = 20e+9 / mUnitRatio; //mMaterialParameters[8+4];

    double rhoL = GetrhoL( p, t );
    double DrhoL_Dt = GetDrhoLDt( 1, p, t );
    double DrhoL_Dp = GetDrhoLDp( 1, p, t );
    double D2rhoL_DpDt = GetD2rhoLDpDt( p, t );

    result = 3 * alphaL * kL / mPorosity
             * ( pow( rhoL, -2 ) * DrhoL_Dp
                 - ( T - T0 ) * 2 * pow( rhoL, -3 ) * DrhoL_Dp * DrhoL_Dt
                 + ( T - T0 ) * pow( rhoL, -2 ) * D2rhoL_DpDt ) ;

    return result;
}
///
double Soil2PhaseRigid::GetD3psiLDp2Dt( const Vector& N, double p, double t )
{
    double result = 0.0;

    double T0 = 273.15;
    for ( unsigned int i = mNodesMin - 1 ; i < mNodesMax ;i++ )
    {
        T0  += mInitialTemperature( i ) * N( i );
    }

    double T = 273.15 + t;
    double alphaL = 1.8e-6; //mMaterialParameters[8+16];
    double kL = 20e+9 / mUnitRatio; //mMaterialParameters[8+4];

    double rhoL = GetrhoL( p, t );
    double DrhoL_Dt = GetDrhoLDt( 1, p, t );
    double DrhoL_Dp = GetDrhoLDp( 1, p, t );
    double D2rhoL_Dp2 = GetDrhoLDp( 2, p, t );
    double D2rhoL_DpDt = GetD2rhoLDpDt( p, t );
    double D3rhoL_Dp2Dt = GetD3rhoLDp2Dt( p, t );

    result = 3 * alphaL * kL / mPorosity
             * (( -2 ) * pow( rhoL, -3 ) * pow( DrhoL_Dp, 2 )
                + ( T - T0 ) * 6 * pow( rhoL, -4 ) * pow( DrhoL_Dp, 2 ) * DrhoL_Dt
                - ( T - T0 ) * 4 * pow( rhoL, -3 ) * DrhoL_Dp * D2rhoL_DpDt
                - ( T - T0 ) * 2 * pow( rhoL, -3 ) * D2rhoL_Dp2 * DrhoL_Dt
                + ( T - T0 ) * pow( rhoL, -2 ) * D3rhoL_Dp2Dt
                + pow( rhoL, -2 ) * D2rhoL_Dp2 ) ;
    return result;
}
///
double Soil2PhaseRigid::GetD3psiLDpDt2( const Vector& N, double p, double t )
{
    double result = 0.0;

    double T0 = 273.15 ;
    for ( unsigned int i = mNodesMin - 1 ; i < mNodesMax ;i++ )
    {
        T0  += mInitialTemperature( i ) * N( i );
    }

    double T = 273.15 + t;
    double alphaL = 1.8e-6; //mMaterialParameters[8+16];
    double kL = 20e+9 / mUnitRatio; //mMaterialParameters[8+4];

    double rhoL = GetrhoL( p, t );
    double DrhoL_Dt = GetDrhoLDt( 1, p, t );
    double D2rhoL_Dt2 = GetDrhoLDt( 2, p, t );
    double DrhoL_Dp = GetDrhoLDp( 1, p, t );
    double D2rhoL_DpDt = GetD2rhoLDpDt( p, t );
    double D3rhoL_DpDt2 = GetD3rhoLDpDt2( p, t );

    result = 3 * alphaL * kL / mPorosity
             * ( 6 * pow( rhoL, -4 ) * DrhoL_Dp * pow( DrhoL_Dt, 2 ) * ( T - T0 )
                 - 4 * pow( rhoL, -3 ) * DrhoL_Dp * DrhoL_Dt
                 - 2 * pow( rhoL, -3 ) * D2rhoL_Dt2 * DrhoL_Dp * ( T - T0 )
                 - 4 * pow( rhoL, -3 ) * DrhoL_Dt * D2rhoL_DpDt * ( T - T0 )
                 + 2 * pow( rhoL, -2 ) * D2rhoL_DpDt
                 + pow( rhoL, -2 ) * D3rhoL_DpDt2 * ( T - T0 ) ) ;
    return result;
}
///--------------------------------------------------------------------------
/// internal energy (solid)
///--------------------------------------------------------------------------
///done
double Soil2PhaseRigid::GetDeSDt( int n )
{
    double result = 0.0;
    double cs = 800 * pow(mUnitRatio, 2.0); //mMaterialParameters[8+9];

    switch ( n )
    {
        case 1: result= cs;
            break;

        case 2: result= 0;
            break;
    }

    return result;
}
///--------------------------------------------------------------------------
/// internal energy (water)
///--------------------------------------------------------------------------
///done
double Soil2PhaseRigid::GetDeLDt( const Vector& N, int n, double p, double t )
{
    double result = 0.0;
    double T = 273.15 + t;

    switch ( n )
    {
        case 1: result= - T * GetDpsiLDt( N, 2, p, t );
            break;

        case 2: result= - GetDpsiLDt( N, 2, p, t )
                            - T * GetDpsiLDt( N, 3, p, t );
            break;
    }

    return result;
}
///
double Soil2PhaseRigid::GetDeLDp( const Vector& N, int n, double p, double t )
{
    double result = 0.0;
    double T = 273.15 + t;

    switch ( n )
    {
        case 1: result=  GetDpsiLDp( N, 1, p, t )
                             - T * GetD2psiLDpDt( N, p, t );
            break;

        case 2: result=  GetDpsiLDp( N, 2, p, t )
                             - T * GetD3psiLDp2Dt( N, p, t );
            break;
    }

    return result;
}
///
double Soil2PhaseRigid::GetD2eLDpDt( const Vector& N, double p, double t )
{
    double result = 0.0;
    double T = 273.15 + t;

    result =  - T * GetD3psiLDpDt2( N, p, t );

    return result;
}
///--------------------------------------------------------------------------
/// heat flux (mixture)
///--------------------------------------------------------------------------
///done
Vector Soil2PhaseRigid::Getq( const Matrix& DN_DX )
{
    Vector result( mDimension );
    noalias( result ) = ZeroVector( mDimension );

    Vector grad_t( mDimension );
    noalias( grad_t ) = Getgradt( DN_DX );

    for ( unsigned int i = 0; i < mDimension; i++ )
    {
        result( i ) = GetDqDgradt() * grad_t( i );
    }

    return result;
}
///done
double Soil2PhaseRigid::GetDqDgradt()
{
    double lambdaS = 2.3 * mUnitRatio; //mMaterialParameters[8+12];
    double lambdaL = 0.58 * mUnitRatio; //mMaterialParameters[8+13];
    double v = 1 - mPorosity;
    double ve = pow( v, 1.0 / 3.0 ) ;
    double lambdaEff = ( lambdaL * lambdaS * ( 1 - ve + v ) + pow( lambdaL, 2.0 ) * ( ve - v ) ) / ( lambdaS * ( 1 - ve ) + lambdaL * ve );
// 	double alpha_grad_t_eff = (1.0-mPorosity)*lambdaS + mPorosity*lambdaL;

    return ( -1 )* lambdaEff; /* 1.11357*/
}
///--------------------------------------------------------------------------
/// momentum supply (water)
///--------------------------------------------------------------------------
///done
Vector Soil2PhaseRigid::GethatpL( const Matrix& DN_DX, double p, double t )
{
    Vector result( mDimension );
    noalias( result ) = ZeroVector( mDimension );

    Vector vL( mDimension );
    noalias( vL ) = GetvL( DN_DX, p, t );

    for ( unsigned int i = 0; i < mDimension; i++ )
    {
        result( i ) = GetDhatpLDvL() * vL( i );
    }

    return result;
}
///done
double Soil2PhaseRigid::GetDhatpLDvL()
{
    double lambda_w_v = 1.0e-3 / pow(mUnitRatio, 3.0); // kg/m^3/s
// 	double lambda_w_v=0.0;
    return -lambda_w_v / mPorosity;  // kg/m^3/s
}
////////////////////////////////////////////////////////
void Soil2PhaseRigid::FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo )
{
    Matrix DN_DX( mNodesNumber, mDimension );
    Matrix DN_DX_e( mNodesNumber, mDimension );
    Matrix nodesLocalCoords( mNodesNumber, mDimension );

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

    for ( unsigned int i = ( mNodesMin - 1 ) ; i < mNodesMax ;i++ )
    {
        Vector local_coords( 3 );
        local_coords( 0 ) = nodesLocalCoords( i, 0 );
        local_coords( 1 ) = nodesLocalCoords( i, 1 );
        local_coords( 2 ) = nodesLocalCoords( i, 2 );

        noalias( DN_DX_e ) = GetGeometry().ShapeFunctionsLocalGradients( DN_DX_e, local_coords );
        noalias( DN_DX ) = prod( DN_DX_e, mInvJ0[0] );

        ( GetGeometry()[i] ).GetSolutionStepValue( WATER_DENSITY ) =
            GetrhoL(( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE ), ( GetGeometry()[i] ).GetSolutionStepValue( TEMPERATURE ) ) * pow(mUnitRatio, 3.0);
//
        for ( unsigned int m = 0; m < mDimension; m++ )
        {
            ( GetGeometry()[i] ).GetSolutionStepValue( WATER_FLOW )( m ) = 3600 * 24.0 / mUnitRatio * GetvL( DN_DX, ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE ), ( GetGeometry()[i] ).GetSolutionStepValue( TEMPERATURE ) )( m );
            ( GetGeometry()[i] ).GetSolutionStepValue( HEAT_FLOW )( m ) = Getq( DN_DX )( m );
        }
    }
}
///
void Soil2PhaseRigid::GetValueOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo )
{}
///not included
void Soil2PhaseRigid::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo )
{}

///not included
void Soil2PhaseRigid::GetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
{}
///not included
void Soil2PhaseRigid::SetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo )
{}

///not included
void Soil2PhaseRigid::SetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
{}

int Soil2PhaseRigid::Check( const Kratos::ProcessInfo& rCurrentProcessInfo )
{
    return 0;
}



} // Namespace Kratos


