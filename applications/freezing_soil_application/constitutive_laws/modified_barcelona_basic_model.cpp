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
/* *********************************************************
 *
 *   Last Modified by:    $Author: janosch $
 *   Date:                $Date: 2009-01-14 17:14:12 $
 *   Revision:            $Revision: 1.13 $
 *
 * ***********************************************************/

// System includes
#include <iostream>
#include <math.h>
#include <cmath>

// External includes

// Project includes

#include "includes/define.h"
#include "constitutive_laws/modified_barcelona_basic_model.h"

#include "includes/constitutive_law.h"

#include "utilities/math_utils.h"
#include "includes/variables.h"
#include "includes/process_info.h"
#include "freezing_soil_application.h"
#include "includes/properties.h"
#include "structural_application/custom_utilities/sd_math_utils.h"
 

namespace Kratos
{

/**
 * TO BE TESTED!!!
 */
ModifiedBarcelonaBasicModel::ModifiedBarcelonaBasicModel() :
        ConstitutiveLaw()
{
}

/**
 * TO BE TESTED!!!
 */
ModifiedBarcelonaBasicModel::~ModifiedBarcelonaBasicModel()
{
}

void ModifiedBarcelonaBasicModel::InitializeMaterial( const Properties& props,
        const GeometryType& geom, const Vector& ShapeFunctionsValues )
{
//     std::cout << "+++ ModifiedBarcelonaBasicModel::InitializeMaterial in kPa+++" << std::endl;

    mOldStress.resize( 6 );
    noalias( mOldStress ) = ZeroVector( 6 );

    mCurrentStress.resize( 6 );
    noalias( mCurrentStress ) = mOldStress;

    mOldStrain.resize( 6 );
    noalias( mOldStrain ) = ZeroVector( 6 );

    mCurrentStrain.resize( 6 );
    mCurrentStrain = mOldStrain;

//     mn0 = props[POROSITY];
    mn0 = props[POROSITY];
    mtstar = props[FIRST_SATURATION_PARAM];
    mm = props[SECOND_SATURATION_PARAM]; 
    
    
    mMaterialParameters = props[MATERIAL_PARAMETERS];
    mUnitRatio = props[SCALE]; 				// Unit convertor: 1:[m], 1000:[mm]
//     mE = props[YOUNG_MODULUS]/ mUnitRatio; 	//Young's modulus (mE)
    me = mn0 / ( 1.0 - mn0 ); 		//initial void ratio (e_0)

    // 8 material constants for BBM
    // -----------------------------------------------------------
    mnu = props[POISSON_RATIO];		//poisson ratio (mue)
    mMSL = mMaterialParameters[0];	//slope of CSL (MSL) of unfrozen state
    mMSC = mMaterialParameters[8]; 	//slope of CSL (MSC) of fully frozen state
    mlambda = mMaterialParameters[1];	//slope of NCL (lambda) of unfrozen state
    mkappa = mMaterialParameters[2];	//slope of URL (kappa) of unfrozen state
    mphiCS = mMaterialParameters[3];	//Friction Angle (phi_cs [degrees])
    mnn = mMaterialParameters[4];	//shape parameter of yield function (n)
    mrr = mMaterialParameters[5];	//spacing ratio (r)
    // -----------------------------------------------------------
    mOldPreconsolidation = mMaterialParameters[6] / mUnitRatio; // preconsolidation (p0[kPa]!!)
    mpi = mMaterialParameters[7] / mUnitRatio; // initial compression (pi[kPa]!!)
//     mpi= 0.0;
    mCurrentPreconsolidation = mOldPreconsolidation; 
 
    mTol = 1.0e-6; 
    mpmin = 1000.0 / mUnitRatio;  // 1 kPa 
    mtheta = ( 1.0 + me ) / ( mlambda - mkappa );
    mSf = 1.2e+6 / mUnitRatio;
//     mK = mE / ( 3.0 * ( 1.0 - 2.0 * mnu ) );
//     mG = mE / ( 2.0 * ( 1.0 + mnu ) );
}

void ModifiedBarcelonaBasicModel::ResetMaterial( const Properties& props,
        const GeometryType& geom, const Vector& ShapeFunctionsValues )
{
    noalias( mOldStrain ) = ZeroVector( 6 );
    noalias( mOldStress ) = ZeroVector( 6 );
//     for ( unsigned int i = 0; i < 3; i++ )
// 	  mOldStress[i] = mPreHydrostaticPressure;

    mOldPreconsolidation = mMaterialParameters[6] / mUnitRatio; // preconsolidation (p0[kPa]!!)

}

void ModifiedBarcelonaBasicModel::InitializeSolutionStep(
    const Properties& props,
    const GeometryType& geom, //this is just to give the array of nodes
    const Vector& ShapeFunctionsValues,
    const ProcessInfo& CurrentProcessInfo )
{ 
}

void ModifiedBarcelonaBasicModel::InitializeNonLinearIteration( const Kratos::Properties& props, const Kratos::ConstitutiveLaw::GeometryType& geom, const Kratos::Vector& ShapeFunctionsValues, const Kratos::ProcessInfo& CurrentProcessInfo )
{
}

void ModifiedBarcelonaBasicModel::FinalizeSolutionStep(
    const Properties& props,
    const GeometryType& geom, //this is just to give the array of nodes
    const Vector& ShapeFunctionsValues,
    const ProcessInfo& CurrentProcessInfo )
{
//     std::cout << "+++ ModifiedBarcelonaBasicModel::FinalizeSolutionStep +++" << std::endl;
    mOldStress = mCurrentStress;
    mOldStrain = mCurrentStrain;
    mOldPreconsolidation = mCurrentPreconsolidation; 
}

void ModifiedBarcelonaBasicModel::CalculateMaterialResponse(
    const Vector& InputVector,		// INPUT: StrainVector(6)+WaterPressure(1)+Temperature(1)
    const Matrix& InternalVariables, 	// INPUT: Empty
    Vector& StressVector, 			// OUTPUT: StressVector(6)
    Matrix& AlgorithmicTangent, 		// OUTPUT: CtanEff(6,6)
    const ProcessInfo& CurrentProcessInfo,
    const Properties& props,
    const GeometryType& geom,
    const Vector& ShapeFunctionsValues,
    bool CalculateStresses,
    int CalculateTangent,
    bool ShowIterativeResults /*SaveInternalVariables*/ )
{
    /// Step 0. Initialized input variables
    //---- assign input variables from element constitutive call: compressive positive, freezing positive!
    for ( unsigned int i = 0; i < 6; i++ )
        mCurrentStrain[i] = -InputVector[i];  // convert compressive strain to be positive
    mpL = InputVector[6]; // liquid pressure: compressive positive
    mt = InputVector[7];  // temperature
//     double n = InputVector[8]; 	// current porosity 
//     me = n / ( 1.0 - n ); 		//current void ratio (e)
    if ( mt > 0.0 )  mt = 0.0;  // set all positive t -> 0 

    if ( ShowIterativeResults )
        std::cout << "INPUT ### " << "t= " << mt << ",\t pL= " << mpL << ",\t pi= " << mpi <<",\t strain= "<< mCurrentStrain << std::endl;
// 
//     //--- resize output variables
    StressVector.resize( 6 );
    AlgorithmicTangent.resize( 6, 6 );

     
    ///---- consider temperature-dependency of preconsolidation
    //++ Method 1; Unchanged
    double p0Tr = mOldPreconsolidation;
//     //++ Method 2; Zhou-LinearInterpolation
// //     double p0Tr = mOldPreconsolidation * ( XL + 2.0 * ( 1 - XL ) );
//     //++ Method 3; Gens-BBM
// //     double rs = 0.9;
// //     double beta = 0.2;
// //     double Sf = 1.2;
// //     double lambdat = mlambda*( (1-rs)*exp(-beta*Sf*t) + rs);
// //     double pRef = 100.0;
// //     double p0Tr = pRef * pow(mOldPreconsolidation/pRef, (mlambda-mkappa)/(lambdat-mkappa) );
//     mtheta = ( 1.0 + me ) / ( mlambda - mkappa ); // kappa(t?!)
// //         double mtheta = ( 1.0 + me ) / ( lambdat - mkappa ); // kappa(t?!)

    ///---- update temperature-dependent Mth and S
    UpdateMS(); 

    /// Step 1. Elastic predictor (Trial step)
    //---- elastic tangent modulus
    double pTrOld = Getp( mOldStress ) + mpi; 
    UpdateKG( pTrOld );
     
    Matrix Ctan( 6, 6 );
    noalias( Ctan ) = GetElasticTangent();

    ///---- trial stress
    Vector sigmaTr( 6 );
    noalias( sigmaTr ) = mOldStress + prod( Ctan, mCurrentStrain - mOldStrain );
    double pTr = Getp( sigmaTr ) + mpi;
    double qTr = Getq( sigmaTr ); 
    if ( pTr + mS < mpmin ) pTr = -mS + mpmin;
    if ( qTr < mTol ) qTr = 0.0;
  
    ///---- check yield function at trial state
    double f = GetYieldFunctionAndDerivatives( pTr, qTr, p0Tr, 0 );
    if ( ShowIterativeResults ) 
        std::cout << "### Trial-Step (p0= " << p0Tr << ") :\t sigmaTr= " << sigmaTr << ",\t f( pTr= " << pTr << ", qTr= " << qTr << " )= " << f << std::endl;
    
    if ( f < mTol )
    {
        mCurrentStress = sigmaTr;
        mCurrentPreconsolidation = p0Tr;
    }
    else
    {
        /// Step 2. Return Map Algorithm
        double pk = pTr;
        double qk = qTr;
        double p0k = p0Tr;

        double DEpsPk = 0.0;
        double DEpsQk = 0.0;

        Vector Rk( 2 );  //Residual
        noalias( Rk ) = GetResidual( pk, qk, p0k, DEpsPk, DEpsQk );
        Vector Dk( 2 );  //Increment
        noalias( Dk ) = ZeroVector( 2 );
        Matrix Kk( 2, 2 );  //Stiffness
        noalias( Kk ) = ZeroMatrix( 2, 2 );
        Matrix InvKk( 2, 2 );  // Inverse of Stiffness
        noalias( InvKk ) = ZeroMatrix( 2, 2 );

        double DetKk;
        int k = 0;
        int maxIteration = 30;

        if ( ShowIterativeResults )
            std::cout << "### Return Map: k= " << k << ": p0k= " << p0k << ", pk= " << pk << ", qk= " << qk << ", Rk= " << Rk << std::endl;

        while ( MathUtils<double>::Norm( Rk ) > mTol && k < maxIteration )
        {
            noalias( Kk ) = GetStiffness( pk, qk, p0k, DEpsPk, DEpsQk );
            MathUtils<double>::InvertMatrix( Kk, InvKk, DetKk );  //DetKk not needed;
            noalias( Dk ) = -prod( InvKk, Rk );
            DEpsPk += Dk[0];
            DEpsQk += Dk[1];
            pk = pTr - mK * DEpsPk ;
            qk = qTr - 3.0 * mG * DEpsQk;
            if ( pk + mS < mpmin ) pk = -mS + mpmin;
            if ( qk < mTol ) qk = 0.0;	 
            p0k = p0Tr * exp( mtheta * DEpsPk );
            UpdateKG(pk);

            noalias( Rk ) = GetResidual( pk, qk, p0k, DEpsPk, DEpsQk );
            k = k + 1;

            if ( ShowIterativeResults )
                std::cout << "###### k= " << k << ": norm(Rk)= " << MathUtils<double>::Norm( Rk ) << ", DEpsPk= " << DEpsPk << ", DEpsQk= " << DEpsQk << ", p0k= " << p0k << ", pk= " << pk << ", qk= " << qk << ", Rk= " << Rk << std::endl;
        }
        
// 	if ( pk + mS < 0.0 || qk < 0.0)
	if ( pk + mS < mTol )  //Modified on 2013-07-08
	{
	    if ( ShowIterativeResults )
		std::cout << "mBBM - ATTENTION: unbearable expansion! pk= " << pk << " < -S= " << -mS << "\t when qk= " << qk << " !! and t= " << mt << std::endl; 
	}
	

        ///---- consider temperature-dependency of preconsolidation
        mCurrentPreconsolidation = p0k;
//         mCurrentPreconsolidation = p0k / ( XL + 2.0 * ( 1 - XL ) );

        if ( qTr < mTol )
            noalias( mCurrentStress ) = ZeroVector( 6 );
        else
            noalias( mCurrentStress ) = qk / qTr * GetDevStress( sigmaTr );
        for ( unsigned int i = 0; i < 3; i++ )
            mCurrentStress [i] += pk - mpi;

        /// Step 3. Elasto-plastic Tangent stiffness
        noalias( Ctan ) = GetElastoPlasticTangent( pk, qk, p0k, sigmaTr, DEpsPk, DEpsQk );

        if ( ShowIterativeResults )
        {
            if ( k == maxIteration )
                std::cout << "mBBM - ATTENTION: MAXIMUM Return Map Iteration Reached !!!!########" << std::endl;
            std::cout << "### k= " << k << ": stress=" << mCurrentStress << std::endl;
            KRATOS_WATCH( Ctan );
        }

    }

    if ( CalculateStresses )
        StressVector = -mCurrentStress; // convert back compressive stress to be negative

    if ( CalculateTangent != 0 )
        AlgorithmicTangent = Ctan;


    return;
}


void ModifiedBarcelonaBasicModel::UpdateMS()
{ 
//     ///---- update temperature-dependent material parameters of freezing materials
    double XL =  pow( 1.0 + pow( - mt / mtstar, 1.0 / ( 1.0 - mm ) ), -mm );
    double cSC = ( 2.1e+6 / mUnitRatio ) * ( 18 * pow( mn0, 3.0 ) - 37 * pow( mn0, 2.0 ) + 20 * mn0 ); // in kPa [see Andersland_Ladanyi:95] for -12°C, peak cohesion 6 MPa
    mMmax = XL * mMSL + ( 1 - XL ) * mMSC;
    mc = ( 1 - XL ) * cSC;
//     mMth = GetMth ( sigmaTr );  // Lode angel dependency
    mMth = mMmax;
    mS = mc / mMth;
}
void ModifiedBarcelonaBasicModel::UpdateKG( const double p )
{ 
    double px = p;
    if ( px < mpmin )    px = mpmin;  
//     double volStrain = 3.0 * Getp (mCurrentStrain); // positive when compressive
//     mK = ( 1.0 + me ) * pTrOld * exp(volStrain) / mkappa;   // stiffer with more compression
    mK = ( 1.0 + me ) * px / mkappa;  
    mG = 1.5 * mK * ( 1.0 - 2.0 * mnu ) / ( 1.0 + mnu );
    return;
}
/// -------------- stress components ---------------
/// ------------------------------------------------
double ModifiedBarcelonaBasicModel::Getp( Vector& stress )
{
    return ( stress[0] + stress[1] + stress[2] ) / 3.0;
}

double ModifiedBarcelonaBasicModel::Getq( Vector& stress )
{
//     return sqrt ( 1.5 ) * MathUtils<double>::Norm ( GetDevStress ( stress ) );
    return sqrt(( pow( stress[0] - stress[1], 2.0 ) +  pow( stress[1] - stress[2], 2.0 ) + pow( stress[2] - stress[0], 2.0 ) ) / 2.0 + 3.0*( stress[3]*stress[3] + stress[4]*stress[4] + stress[5]*stress[5] ) );
}

Vector ModifiedBarcelonaBasicModel::GetDevStress( Vector& stress )
{
    Vector result( 6 );
    noalias( result ) = stress;
    double p = Getp( stress );

    for ( unsigned int i = 0; i < 3; i++ )
        result [i] -= p;

    return result;
}


/// ------------------- solve ----------------------
/// ------------------------------------------------
Matrix ModifiedBarcelonaBasicModel::GetElasticTangent()
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

Matrix ModifiedBarcelonaBasicModel::GetElastoPlasticTangent( double p, double q, double p0, Vector& stressTr, double DEpsP, double DEpsQ )
{
    Matrix result( 6, 6 );
    noalias( result ) = ZeroMatrix( 6, 6 );

    double Df_Dp = GetYieldFunctionAndDerivatives( p, q, p0, 1 );
    double Df_Dq = GetYieldFunctionAndDerivatives( p, q, p0, 2 );
    double Df_Dp0 = GetYieldFunctionAndDerivatives( p, q, p0, 3 );

    double Dg_Dp = GetPotentialFunctionDerivatives( p, q, 1 );
    double Dg_Dq = GetPotentialFunctionDerivatives( p, q, 2 );
    double D2g_Dp2 = GetPotentialFunctionDerivatives( p, q, 3 );
    double D2g_Dq2 = GetPotentialFunctionDerivatives( p, q, 4 );
    double D2g_DpDq = GetPotentialFunctionDerivatives( p, q, 5 );

    double A11 = Df_Dp * mK - Df_Dp0 * p0 * mtheta;
    double A12 = Df_Dq * 3.0 * mG;
    double A21 = Dg_Dq + mK * ( -DEpsP * D2g_DpDq + DEpsQ * D2g_Dp2 );
    double A22 = -Dg_Dp + 3.0 * mG * ( -DEpsP * D2g_Dq2 + DEpsQ * D2g_DpDq );

    double B11 = Df_Dp * mK;
    double B12 = Df_Dq * 2.0 * mG;
    double B21 = mK * ( -DEpsP * D2g_DpDq + DEpsQ * D2g_Dp2 );
    double B22 = 2.0 * mG * ( -DEpsP * D2g_Dq2 + DEpsQ * D2g_DpDq );

    double DetA = A11 * A22 - A21 * A12;

    double mpI = ( A22 * B11 - A12 * B21 ) / DetA;
    double mpN = ( A22 * B12 - A12 * B22 ) / DetA;
    double mqI = ( A11 * B21 - A21 * B11 ) / DetA;
    double mqN = ( A11 * B22 - A21 * B12 ) / DetA;

    double qTr = Getq( stressTr );

    Vector N2( 6 );
    if ( qTr < mTol )
        noalias( N2 ) = ZeroVector( 6 );
    else
        noalias( N2 ) = 1.5 * GetDevStress( stressTr ) / qTr;
    Matrix I4( 6, 6 );
    noalias( I4 ) = ZeroMatrix( 6, 6 );
    Vector I2( 6 );
    noalias( I2 ) = ZeroVector( 6 );

    for ( unsigned int i = 0; i < 3; i++ )
    {
        I4( i, i ) = 1.0;
        I4( i + 3, i + 3 ) = 0.5;
        I2 [i] = 1.0;
    }

    for ( unsigned int i = 0; i < 6; i++ )
        for ( unsigned int j = 0; j < 6; j++ )
            result( i, j ) = 2.0 * mG * q / qTr * I4( i , j )
                             + ( mK * ( 1.0 - mpI ) - mG / 1.5 * q / qTr ) * I2 [i] * I2 [j]
                             + 4.0 / 3.0 * mG * ( 1.0 - q / qTr - 1.5 * mqN ) * N2 [i] * N2 [j]
                             - mK * mpN * I2 [i] * N2 [j]
                             - 2.0 * mG * mqI * N2 [i] * I2 [j];

    return result;

}

Vector ModifiedBarcelonaBasicModel::GetResidual( double p, double q,  double p0, double DEpsP, double DEpsQ )
{
    Vector result( 2 );
    noalias( result ) = ZeroVector( 2 );
    double Dg_Dp = GetPotentialFunctionDerivatives( p, q, 1 );
    double Dg_Dq = GetPotentialFunctionDerivatives( p, q, 2 );

    result[0] = DEpsP * Dg_Dq - DEpsQ * Dg_Dp;
    result[1] = GetYieldFunctionAndDerivatives( p, q, p0, 0 );

    return result;
}


Matrix ModifiedBarcelonaBasicModel::GetStiffness( double p, double q,  double p0, double DEpsP, double DEpsQ )
{
    Matrix result( 2, 2 );
    noalias( result ) = ZeroMatrix( 2, 2 );

    double Df_Dp = GetYieldFunctionAndDerivatives( p, q, p0, 1 );
    double Df_Dq = GetYieldFunctionAndDerivatives( p, q, p0, 2 );
    double Df_Dp0 = GetYieldFunctionAndDerivatives( p, q, p0, 3 );

    double Dg_Dp = GetPotentialFunctionDerivatives( p, q, 1 );
    double Dg_Dq = GetPotentialFunctionDerivatives( p, q, 2 );
    double D2g_Dp2 = GetPotentialFunctionDerivatives( p, q, 3 );
    double D2g_Dq2 = GetPotentialFunctionDerivatives( p, q, 4 );
    double D2g_DpDq = GetPotentialFunctionDerivatives( p, q, 5 );

    result( 0, 0 ) = Dg_Dq - mK * ( DEpsP * D2g_DpDq  - DEpsQ * D2g_Dp2 );
    result( 0, 1 ) = -Dg_Dp + 3.0 * mG * ( DEpsQ * D2g_DpDq - DEpsP * D2g_Dq2 );
    result( 1, 0 ) = -mK * Df_Dp + p0 * mtheta * Df_Dp0;
    result( 1, 1 ) = -3.0 * mG * Df_Dq;

    return result;
}

/// ----------- functions and derivatives-----------
/// ------------------------------------------------
///DONE
double ModifiedBarcelonaBasicModel::GetYieldFunctionAndDerivatives( double p, double q, double p0, int index )
{
    // CASM: S and c/M are equal
    if ( p + mS < mpmin ) p = -mS + mpmin; // quadratic convergence speed lost! 
    if ( q < mTol )	q = 0.0;

    // Lai: S and c/M are not equal
//     double Bstar = 6, p0star = 12; //MPa
//     double BLai = Bstar * ( mMth * p0 + mc ) / ( mMth * p0star + mc );
//     double mS = p0 - 2.0 * ( BLai - mc ) / mMth;
    switch ( index )
    {
        case 0: // f
//             if ( p <= mpmin )
//                 return q;
//             else
            {

//         return pow ( q / mMth, mnn) + pow ( p, mnn) * log ( p / p0 ) / log ( mrr ); //CASM
                return pow( q / mMth, mnn ) + pow( p + mS, mnn ) * log(( p + mS ) / ( p0 + mS ) ) / log( mrr ); //CASM + Lai
            }

        case 1: // Df_Dp
//             if ( p <= mpmin )
//                 return 0.0;
//             else
//         return pow ( p, mnn- 1.0 ) * ( 1.0 + mnn* log ( p / p0 ) ) / log ( mrr ); //CASM
            return pow( p + mS, mnn - 1.0 ) * ( 1.0 + mnn* log(( p + mS ) / ( p0 + mS ) ) ) / log( mrr ); //CASM + Lai

        case 2: // Df_Dq
//             if ( p <= mpmin )
//                 return 1.0;
//             else
            return mnn* pow( q / mMth, mnn - 1.0 ) / mMth;

        case 3: // Df_Dp0
//             if ( p <= mpmin )
//                 return 0.0;
//             else
//         return -pow ( p, mnn) / ( p0 * log ( mrr ) );//CASM
            return -pow(( p + mS ), mnn ) / (( p0 + mS ) * log( mrr ) );//CASM + Lai
    }
    
    return 0.0;
}

double ModifiedBarcelonaBasicModel::GetPotentialFunctionDerivatives( double p, double q, int index )
{
    if ( p + mS < mTol )	p = -mS + mpmin; // quadratic convergence speed lost! 
    if ( q < mTol )		q = 0.0;

    double N1 = 3.0 + 2.0 * mMth;
    double D1 = 3.0 * ( p + mS ) + 2.0 * q;
    double N2 = 3.0 - mMth; 
    double D2 = 3.0 * ( p + mS ) - q;

    switch ( index )
    {
        case 1: // Dg_Dp
            return 3.0 * ( N1 / D1 - N2 / D2 );

        case 2: // Dg_Dq
            return 2.0 * N1 / D1 + N2 / D2;

        case 3: // D2g_Dp2
            return - 9.0 * N1 / ( D1 * D1 ) + 9.0 * N2 / ( D2 * D2 );

        case 4: // D2g_Dq2
            return - 4.0 * N1 / ( D1 * D1 ) + N2 / ( D2 * D2 ) ;

        case 5: // D2g_DpDq
            return - 6.0 * N1 / ( D1 * D1 ) - 3.0 * N2 / ( D2 * D2 ) ;
//                 return 0.0;
    }

    return 0.0;
}

double ModifiedBarcelonaBasicModel::GetMth( Vector& stress )
{
    double result;
    double alpha = ( 3.0 - sin( mphiCS * PI / 180.0 ) ) / ( 3.0 + sin( mphiCS * PI / 180.0 ) );
    double alpha4 = pow( alpha, 4.0 );
    Vector s( 6 );
    noalias( s ) = GetDevStress( stress );
//     double q = Getq ( stress );
//     double J3 = s[0] * s[1] * s[2] + 2.0 * s[3] * s[4] * s[5] - s[0] * s[4] * s[4] - s[1] * s[5] * s[5] - s[2] * s[3] * s[3];
//
//     double nom = 2.0 * alpha4;
//     double denom = 1.0 + alpha4 + (1.0 - alpha4) * 13.5 * J3 / pow(q, 3.0);

    double J2 = 0.5 * ( s[0] * s[0] + s[1] * s[1] + s[2] * s[2] ) + s[3] * s[3] + s[4] * s[4] + s[5] * s[5];
    double J3 = s[0] * s[1] * s[2] + 2.0 * s[3] * s[4] * s[5] - s[0] * s[4] * s[4] - s[1] * s[5] * s[5] - s[2] * s[3] * s[3];


    double sin3ThL = J3 / sqrt( J2 * J2 * J2 );
    sin3ThL *=  - 3.0 * sqrt( 3.0 ) / 2.0;

    double nom = 2.0 * alpha4;
    double denom = 1.0 + alpha4 + ( 1.0 - alpha4 ) * sin3ThL;

    result = mMmax * pow( nom / denom, 0.25 );

    return result;
//     return mMmax;
}


///

bool ModifiedBarcelonaBasicModel::Has( const Variable<double>& rThisVariable )
{
    if ( rThisVariable == PRECONSOLIDATION ) return true;

    if ( rThisVariable == EQUIVALENT_VOLUMETRIC_STRAIN ) return true;

    if ( rThisVariable == EQUIVALENT_DEVIATORIC_STRAIN ) return true;

    if ( rThisVariable == EQUIVALENT_VOLUMETRIC_STRESS ) return true;

    if ( rThisVariable == EQUIVALENT_DEVIATORIC_STRESS ) return true;

    if ( rThisVariable == LOG_EQUIVALENT_VOLUMETRIC_STRESS ) return true;

    if ( rThisVariable == TOTAL_STRESS ) return true;
    
    if ( rThisVariable == PRESTRESS ) return true;
    
    return false;
}

bool ModifiedBarcelonaBasicModel::Has( const Variable<Vector>& rThisVariable )
{
    if ( rThisVariable == INTERNAL_VARIABLES ) return true;

    if ( rThisVariable == INSITU_STRESS ) return true;

    return false;
}

bool ModifiedBarcelonaBasicModel::Has( const Variable<Matrix>& rThisVariable )
{
    return Kratos::ConstitutiveLaw::Has( rThisVariable );
}

bool ModifiedBarcelonaBasicModel::Has( const Kratos::Variable< Kratos::array_1d< double, 3 > >& rThisVariable )
{
    return Kratos::ConstitutiveLaw::Has( rThisVariable );
}

bool ModifiedBarcelonaBasicModel::Has( const Kratos::Variable< Kratos::array_1d< double, 6 > >& rThisVariable )
{
    return Kratos::ConstitutiveLaw::Has( rThisVariable );
}

double& ModifiedBarcelonaBasicModel::GetValue( const Variable<double>& rThisVariable,
        double& rValue )
{
    if ( rThisVariable == PRECONSOLIDATION )
        rValue = mCurrentPreconsolidation;

    if ( rThisVariable == EQUIVALENT_VOLUMETRIC_STRAIN ) //epsilon_p
        rValue = -( mCurrentStrain[0] + mCurrentStrain[1] + mCurrentStrain[2] );

    if ( rThisVariable == EQUIVALENT_DEVIATORIC_STRAIN ) //epsilon_q
    {
//         rValue = sqrt( 2.0*( pow(mCurrentStrain[0]-mCurrentStrain[1], 2.0) + pow(mCurrentStrain[1]-mCurrentStrain[2], 2.0) + pow(mCurrentStrain[2]-mCurrentStrain[0], 2.0) ) + 3.0*( pow(mCurrentStrain[3]/2.0,2.0)+pow(mCurrentStrain[4]/2.0,2.0)+pow(mCurrentStrain[5]/2.0,2.0) ) )/3.0;
        rValue = sqrt( 2.0 * ( pow( mCurrentStrain[0] - mCurrentStrain[1], 2.0 ) + pow( mCurrentStrain[1] - mCurrentStrain[2], 2.0 ) + pow( mCurrentStrain[2] - mCurrentStrain[0], 2.0 ) ) + 3.0 * ( pow( mCurrentStrain[3], 2.0 ) + pow( mCurrentStrain[4], 2.0 ) + pow( mCurrentStrain[5], 2.0 ) ) ) / 3.0;
    }

    if ( rThisVariable == EQUIVALENT_VOLUMETRIC_STRESS ) //p
        rValue = Getp( mCurrentStress );

//     if ( rThisVariable == NET_EQUIVALENT_VOLUMETRIC_STRESS ) //p
//         rValue = Getp( mCurrentStress );

    if ( rThisVariable == EQUIVALENT_DEVIATORIC_STRESS ) //q
        rValue = Getq( mCurrentStress );

    if ( rThisVariable == LOG_EQUIVALENT_VOLUMETRIC_STRESS ) //ln(p)
        rValue = log10( Getp( mCurrentStress ));

    if ( rThisVariable == TOTAL_STRESS ) //ptot
    {
	double XL =  pow( 1.0 + pow( -mt / mtstar, 1.0 / ( 1.0 - mm ) ), -mm );
        rValue = Getp( mCurrentStress ) + mpL + ( 1.0 - XL ) * mSf * ( -mt ) ;
    }
    
    if ( rThisVariable == PRESTRESS ) //pi
        rValue = mpi;
      
      
    return rValue;
}

Vector& ModifiedBarcelonaBasicModel::GetValue( const Variable<Vector>& rThisVariable,
        Vector& rValue )
{
    if ( rThisVariable == INTERNAL_VARIABLES )
        rValue = mCurrentStrain;

    if ( rThisVariable == INSITU_STRESS )
        rValue = mCurrentStress;

    return rValue;
}

Matrix& ModifiedBarcelonaBasicModel::GetValue( const Variable<Matrix>& rThisVariable,
        Matrix& rValue )
{
    return rValue;
}

Kratos::array_1d< double, 3 >& ModifiedBarcelonaBasicModel::GetValue( const Kratos::Variable< Kratos::array_1d< double, 3 > >& rVariable, Kratos::array_1d< double, 3 >& rValue )
{
    return rValue;
}

Kratos::array_1d< double, 6 >& ModifiedBarcelonaBasicModel::GetValue( const Kratos::Variable< Kratos::array_1d< double, 6 > >& rVariable, Kratos::array_1d< double, 6 >& rValue )
{
    return rValue;
}

void ModifiedBarcelonaBasicModel::SetValue( const Variable<double>& rThisVariable,
        const double& rValue, const ProcessInfo& rCurrentProcessInfo )
{
    if ( rThisVariable == PRECONSOLIDATION )
        mOldPreconsolidation = rValue;
    
//     if ( rThisVariable == PRESTRESS )
//       if(rValue < 0)
// 	mpi= Getp( mCurrentStress );
//       else
//         mpi = rValue;
}

void ModifiedBarcelonaBasicModel::SetValue( const Variable<Vector>& rThisVariable,
        const Vector& rValue, const ProcessInfo& rCurrentProcessInfo )
{
}

void ModifiedBarcelonaBasicModel::SetValue( const Variable<Matrix>& rThisVariable,
        const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo )
{
}

void ModifiedBarcelonaBasicModel::SetValue( const Variable<array_1d<double, 3> >& rThisVariable, const array_1d <double, 3 > & rValue,
        const ProcessInfo& rCurrentProcessInfo )
{
}

void ModifiedBarcelonaBasicModel::SetValue( const Kratos::Variable< Kratos::array_1d< double, 6 > >& rVariable, const Kratos::array_1d< double, 6 >& Value, const Kratos::ProcessInfo& rCurrentProcessInfo )
{
}

bool ModifiedBarcelonaBasicModel::ValidateInput( const Kratos::Properties& props )
{
    return false;
}

int ModifiedBarcelonaBasicModel::Check( const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo )
{
    if ( ValidateInput( props ) )
        return 0;

    return -1;
}

} // Namespace Kratos
