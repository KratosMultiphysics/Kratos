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

// External includes
#include<cmath>

// Project includes

#include "includes/define.h"
#include "constitutive_laws/clay_and_sand_model.h"

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
ClayAndSandModel::ClayAndSandModel() :
        ConstitutiveLaw()
{
}

/**
 * TO BE TESTED!!!
 */
ClayAndSandModel::~ClayAndSandModel()
{
}

void ClayAndSandModel::InitializeMaterial ( const Properties& props,
        const GeometryType& geom, const Vector& ShapeFunctionsValues )
{
//     std::cout << "+++ ClayAndSandModel::InitializeMaterial +++" << std::endl;
    mOldStress.resize ( 6 );
    noalias ( mOldStress ) = ZeroVector ( 6 );
    
//     mPreHydrostaticPressure = props[PRESSURE];
//     for ( unsigned int i = 0; i < 3; i++ )
//       mOldStress[i] = mPreHydrostaticPressure; 
    
    mCurrentStress.resize ( 6 );
    noalias ( mCurrentStress ) = mOldStress;
    
    mOldStrain.resize ( 6 );
    noalias ( mOldStrain ) = ZeroVector ( 6 );

    mCurrentStrain.resize ( 6 );
    mCurrentStrain = mOldStrain;

    mMaterialParameters = props[MATERIAL_PARAMETERS]; 
    double n0 =  props[POROSITY];
    me = n0 / ( 1.0 - n0 );         //initial void ratio (e_0) 
    mnu = props[POISSON_RATIO];        //poisson ratio (mue)
    mMmax = mMaterialParameters[0];    //slope of CSL (M)
    mlambda = mMaterialParameters[1];    //slope of NCL (lambda)
    mkappa = mMaterialParameters[2];    //slope of URL (kappa)
    mphiCS = mMaterialParameters[3];    //Friction Angle (phi_cs [degrees])
    mn = mMaterialParameters[4];    //shape parameter of yield function (n)
    mr = mMaterialParameters[5];    //spacing ratio (r)
    mOldPreconsolidation = mMaterialParameters[6]; //preconsolidation (p0[kPa]!!)
    mCurrentPreconsolidation = mOldPreconsolidation;
     
    mTol = 1.0e-10;
    mpmin = 1.0;
}

void ClayAndSandModel::ResetMaterial ( const Properties& props,
                                       const GeometryType& geom, const Vector& ShapeFunctionsValues )
{
    noalias ( mOldStrain ) = ZeroVector ( 6 );
    
    noalias ( mOldStress ) = ZeroVector ( 6 );  
//     for ( unsigned int i = 0; i < 3; i++ )
//       mOldStress[i] = mPreHydrostaticPressure; 
    
    mOldPreconsolidation = mMaterialParameters[6];
     
}

void ClayAndSandModel::InitializeSolutionStep (
    const Properties& props,
    const GeometryType& geom, //this is just to give the array of nodes
    const Vector& ShapeFunctionsValues,
    const ProcessInfo& CurrentProcessInfo )
{
}

void ClayAndSandModel::InitializeNonLinearIteration ( const Kratos::Properties& props, const Kratos::ConstitutiveLaw::GeometryType& geom, const Kratos::Vector& ShapeFunctionsValues, const Kratos::ProcessInfo& CurrentProcessInfo )
{
}

void ClayAndSandModel::FinalizeSolutionStep (
    const Properties& props,
    const GeometryType& geom, //this is just to give the array of nodes
    const Vector& ShapeFunctionsValues,
    const ProcessInfo& CurrentProcessInfo )
{
//     std::cout << "+++ ClayAndSandModel::FinalizeSolutionStep +++" << std::endl;
    mOldStress = mCurrentStress;
    mOldStrain = mCurrentStrain;
    mOldPreconsolidation = mCurrentPreconsolidation;
}

void ClayAndSandModel::CalculateMaterialResponse ( const Vector& StrainVector,
        const Matrix& DeformationGradient, Vector& StressVector,
        Matrix& AlgorithmicTangent, const ProcessInfo& CurrentProcessInfo,
        const Properties& props, const GeometryType& geom,
        const Vector& ShapeFunctionsValues, bool CalculateStresses,
        int CalculateTangent, bool SaveInternalVariables )
{
    mCurrentStrain = -StrainVector; // convert compressive strain to be positive
    StressVector.resize ( 6 );
    AlgorithmicTangent.resize ( 6, 6 );

    // 1. Elastic predictor (Trial step)
    // 2013-01-02 : update the porosoity from the soil element
//     double n =  DeformationGradient(0,0); 
//     me = n / ( 1.0 - n );  
    
    double pTrOld = Getp ( mOldStress );
    if ( pTrOld < mpmin ) pTrOld = mpmin;  
    mK = ( 1.0 + me ) * pTrOld  / mkappa;   // approximate pk by pTrOld -->  non-quadratic convergence ??
    mG = 1.5 * mK * ( 1.0 - 2.0 * mnu ) / ( 1.0 + mnu ); 
    
    // 2013-01-02 : return the elastic stiffness moduli to the soil element
//     DeformationGradient(1,0) = mK;
//     DeformationGradient(1,1) = mG;
    
    Matrix Ctan ( 6, 6 );
    noalias ( Ctan ) = GetElasticTangent();

    Vector sigmaTr ( 6 );
    noalias ( sigmaTr ) = mOldStress + prod ( Ctan, mCurrentStrain - mOldStrain );
    mMth = mMmax;
//     mMth = GetMth ( sigmaTr );
    double p0Tr = mOldPreconsolidation;
    double pTr = Getp ( sigmaTr );
    double qTr = Getq ( sigmaTr );
    if ( pTr < mpmin ) pTr = mpmin; 
    if ( qTr < mTol ) qTr = mTol;  
        
    double f = GetYieldFunctionAndDerivatives ( pTr, qTr, p0Tr, 0 );

if(SaveInternalVariables)
std::cout << "### Trial Step ( p0= "<<p0Tr<<" ) :\t pTr= "<<pTr<<",\t qTr= "<<qTr<<";\t f= "<<f<<";\t K= "<<mK<<";\t G= "<<mG<<";\t stress= "<<sigmaTr<< std::endl;
    
    if ( f < mTol )
    { 
        mCurrentStress = sigmaTr;
        mCurrentPreconsolidation = p0Tr;
    }
    else
    {
        // 2. Return Map Algorithm
        double pk = pTr;
        double qk = qTr;
        double p0k = p0Tr; 
    
        double DEpsPk = 0.0;
        double DEpsQk = 0.0;

        Vector Rk ( 2 ); //Residual
        noalias ( Rk ) = GetResidual ( pk, qk, p0k, DEpsPk, DEpsQk );
        Vector Dk ( 2 ); //Increment
        noalias ( Dk ) = ZeroVector ( 2 );
        Matrix Kk ( 2, 2 ); //Stiffness
        noalias ( Kk ) = ZeroMatrix ( 2, 2 );
        Matrix InvKk ( 2, 2 ); // Inverse of Stiffness
        noalias ( InvKk ) = ZeroMatrix ( 2, 2 );

        double DetKk;
        int k = 0;
        int maxIteration = 30;
        double theta = ( 1.0 + me ) / ( mlambda - mkappa );
 
        if(SaveInternalVariables)
        std::cout << "##### Return Map: k= " << k << ": p0k= " << p0k << ", pk= " << pk << ", qk= " << qk<< ", Rk= " << Rk << std::endl;

        while ( MathUtils<double>::Norm ( Rk ) > mTol && k < maxIteration )
        {
            noalias ( Kk ) = GetStiffness ( pk, qk, p0k, DEpsPk, DEpsQk );
            MathUtils<double>::InvertMatrix ( Kk, InvKk, DetKk ); //DetKk not needed; 
            noalias ( Dk ) = -prod ( InvKk, Rk );
            DEpsPk += Dk[0];
            DEpsQk += Dk[1];
            pk = pTr - mK * DEpsPk ;
            qk = qTr - 3.0 * mG * DEpsQk;
            if ( pk < mpmin ) pk = mpmin; 
            if ( qk < mTol ) qk = 0.0; 
            p0k = p0Tr * exp ( theta * DEpsPk );
        
            noalias ( Rk ) = GetResidual ( pk, qk, p0k, DEpsPk, DEpsQk ); 
            k = k + 1;
            if(SaveInternalVariables)
                std::cout << "++ k= " << k << ": norm(Rk)= " << MathUtils<double>::Norm ( Rk )<<", DEpsPk= " << DEpsPk << ", DEpsQk= " << DEpsQk << ", p0k= " << p0k << ", pk= " << pk << ", qk= " << qk<< ", Rk= " << Rk << std::endl;
        }
         
        if ( SaveInternalVariables && k == maxIteration )
        {
            std::cout << "########## MAXIMUM Iteration Reached!: sigmak= "<< mCurrentStress <<" ########" << std::endl;
            KRATOS_THROW_ERROR(std::logic_error, "", "")
        }

        mCurrentPreconsolidation = p0k;
    
        if ( qTr < mTol )
            noalias(mCurrentStress) = ZeroVector(6);
        else
            noalias(mCurrentStress) = qk / qTr * GetDevStress ( sigmaTr );
        for ( unsigned int i = 0; i < 3; i++ )
            mCurrentStress [i] += pk; 
        
        // 3. Elasto-plastic Tangent
        noalias ( Ctan ) = GetElastoPlasticTangent ( pk, qk, p0k, sigmaTr, DEpsPk, DEpsQk );
    }

    if ( CalculateStresses )
        noalias(StressVector) = -mCurrentStress; // convert back compressive stress to be negative

    if ( CalculateTangent != 0 )
        noalias(AlgorithmicTangent) = Ctan;

if(SaveInternalVariables)
{
//    KRATOS_WATCH(StressVector)
//    KRATOS_WATCH(AlgorithmicTangent);
}

    return;
}

/// -------------- stress components ---------------
/// ------------------------------------------------
double ClayAndSandModel::Getp ( Vector& stress )
{
    return ( stress[0] + stress[1] + stress[2] ) / 3.0;
}

double ClayAndSandModel::Getq ( Vector& stress )
{
//     return sqrt ( 1.5 ) * MathUtils<double>::Norm ( GetDevStress ( stress ) ); 
    return sqrt( ( pow(stress[0]-stress[1],2.0) +  pow(stress[1]-stress[2],2.0) + pow(stress[2]-stress[0],2.0) )/2.0 + 3.0*( stress[3]*stress[3] + stress[4]*stress[4] + stress[5]*stress[5] ) );
}

Vector ClayAndSandModel::GetDevStress ( Vector& stress )
{
    Vector result ( 6 );
    noalias ( result ) = stress;
    double p = Getp ( stress );

    for ( unsigned int i = 0; i < 3; i++ )
        result [i] -= p;

    return result;
}


/// ------------------- solve ----------------------
/// ------------------------------------------------
Matrix ClayAndSandModel::GetElasticTangent ()
{
    Matrix result ( 6, 6 );
    noalias ( result ) = ZeroMatrix ( 6, 6 );

    double c1 = mK + 4.0 * mG / 3.0;
    double c2 = mK - 2.0 * mG / 3.0;

    for ( unsigned int i = 0; i < 3; i++ )
    {
        result ( i, i ) = c1;
        result ( i + 3, i + 3 ) = mG;

        for ( unsigned int j = 0; j < 3; j++ )
            if ( i != j )
                result ( i, j ) = c2;
    }

    return result;
}

Matrix ClayAndSandModel::GetElastoPlasticTangent ( double p, double q, double p0, Vector& stressTr, double DEpsP, double DEpsQ )
{
    Matrix result ( 6, 6 );
    noalias ( result ) = ZeroMatrix ( 6, 6 ); 
    
    double Df_Dp = GetYieldFunctionAndDerivatives ( p, q, p0, 1 );
    double Df_Dq = GetYieldFunctionAndDerivatives ( p, q, p0, 2 );
    double Df_Dp0 = GetYieldFunctionAndDerivatives ( p, q, p0, 3 );

    double Dg_Dp = GetPotentialFunctionDerivatives ( p, q, 1 );
    double Dg_Dq = GetPotentialFunctionDerivatives ( p, q, 2 ); 
    double D2g_Dp2 = GetPotentialFunctionDerivatives ( p, q, 3 );
    double D2g_Dq2 = GetPotentialFunctionDerivatives ( p, q, 4 );
    double D2g_DpDq = GetPotentialFunctionDerivatives ( p, q, 5 );

    double A11 = Df_Dp * mK - Df_Dp0 * p0 * ( 1.0 + me ) / ( mlambda - mkappa );
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
    
    double qTr = Getq ( stressTr );

    Vector N2 ( 6 );
    if( qTr < mTol ) 
        noalias ( N2 ) = ZeroVector(6);
    else
        noalias ( N2 ) = 1.5 * GetDevStress ( stressTr ) / qTr;
    Matrix I4 ( 6, 6 );
    noalias ( I4 ) = ZeroMatrix ( 6, 6 );
    Vector I2 ( 6 );
    noalias ( I2 ) = ZeroVector ( 6 );

    for ( unsigned int i = 0; i < 3; i++ )
    {
        I4 ( i, i ) = 1.0;
        I4 ( i + 3, i + 3 ) = 0.5;
        I2 [i] = 1.0;
    }
    
//    std::cout << "q at GetElastoPlasticTangent: " << q << std::endl;
//    std::cout << "qTr at GetElastoPlasticTangent: " << qTr << std::endl;

    if(qTr > mTol)
        for ( unsigned int i = 0; i < 6; i++ )
            for ( unsigned int j = 0; j < 6; j++ )
                result ( i, j ) = 2.0 * mG * q / qTr * I4 ( i , j )
                                  + ( mK * ( 1.0 - mpI ) - mG / 1.5 * q / qTr ) * I2 [i] * I2 [j]
                                  + 4.0 / 3.0 * mG * ( 1.0 - q / qTr - 1.5 * mqN ) * N2 [i] * N2 [j]
                                  - mK * mpN * I2 [i] * N2 [j]
                                  - 2.0 * mG * mqI * N2 [i] * I2 [j];
    else
        for ( unsigned int i = 0; i < 6; i++ )
            for ( unsigned int j = 0; j < 6; j++ )
                result ( i, j ) = 2.0 * mG * (I4 ( i , j ) - 1.0 / 3 * I2[i] * I2[j])
                                  + mK * ( 1.0 - mpI ) * I2 [i] * I2 [j];
//                                  - mK * mpN * I2 [i] * N2 [j];

//    std::cout << "result at GetElastoPlasticTangent: " << result << std::endl;

    return result;

}

Vector ClayAndSandModel::GetResidual ( double p, double q,  double p0, double DEpsP, double DEpsQ )
{
    Vector result ( 2 );
    noalias ( result ) = ZeroVector ( 2 );
    double Dg_Dp = GetPotentialFunctionDerivatives ( p, q, 1 );
    double Dg_Dq = GetPotentialFunctionDerivatives ( p, q, 2 ); 

    result[0] = DEpsP * Dg_Dq - DEpsQ * Dg_Dp;
    result[1] = GetYieldFunctionAndDerivatives ( p, q, p0, 0 );

    return result;
}


Matrix ClayAndSandModel::GetStiffness ( double  p, double q,  double p0, double DEpsP, double DEpsQ )
{
    Matrix result ( 2, 2 );
    noalias ( result ) = ZeroMatrix ( 2, 2 ); 

    double Df_Dp = GetYieldFunctionAndDerivatives ( p, q, p0, 1 );
    double Df_Dq = GetYieldFunctionAndDerivatives ( p, q, p0, 2 );
    double Df_Dp0 = GetYieldFunctionAndDerivatives ( p, q, p0, 3 );

    double Dg_Dp = GetPotentialFunctionDerivatives ( p, q, 1 );
    double Dg_Dq = GetPotentialFunctionDerivatives ( p, q, 2 ); 
    double D2g_Dp2 = GetPotentialFunctionDerivatives ( p, q, 3 );
    double D2g_Dq2 = GetPotentialFunctionDerivatives ( p, q, 4 );
    double D2g_DpDq = GetPotentialFunctionDerivatives ( p, q, 5 );
    
    double theta =  ( 1.0 + me ) / ( mlambda - mkappa );

    result ( 0, 0 ) = Dg_Dq - mK * ( DEpsP * D2g_DpDq  - DEpsQ * D2g_Dp2 );
    result ( 0, 1 ) = -Dg_Dp + 3.0 * mG * ( DEpsQ * D2g_DpDq - DEpsP * D2g_Dq2 );
    result ( 1, 0 ) = -mK * Df_Dp + p0 * theta * Df_Dp0;
    result ( 1, 1 ) = -3.0 * mG * Df_Dq;

    return result;
}

/// ----------- functions and derivatives-----------
/// ------------------------------------------------
double ClayAndSandModel::GetYieldFunctionAndDerivatives ( double p, double q, double p0, int index )
{ 
    if ( p < mpmin )    p = mpmin;
    if ( q < mTol )    q = 0.0;
      
    switch ( index )
    {
        case 0: // f 
//             if ( p <= mpmin )
//                 return q;
//             else
                return pow ( q / mMth, mn ) + pow ( p, mn ) * log ( p / p0 ) / log ( mr );

        case 1: // Df_Dp
//             if ( p <= mpmin )
//                 return 0.0;
//             else
                return pow ( p, mn - 1.0 ) * ( 1.0 + mn * log ( p / p0 ) ) / log ( mr );

        case 2: // Df_Dq
//             if ( p <= mpmin )
//                 return 1.0;
//             else
                return mn * pow ( q / mMth, mn - 1.0 ) / mMth;

        case 3: // Df_Dp0
//             if ( p <= mpmin )
//                 return 0.0;
//             else
                return -pow ( p, mn ) / ( p0 * log ( mr ) ); 
    }
}

double ClayAndSandModel::GetPotentialFunctionDerivatives ( double  p, double q, int index )
{
    if ( p < mpmin )    p = mpmin;
    if ( q < mTol )    q = 0.0;
    
    double N1 = 3.0 + 2.0 * mMth;
    double D1 = 3.0 * p + 2.0 * q;
    double N2 = 3.0 - mMth;
    double D2 = 3.0 * p - q;

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

double ClayAndSandModel::GetMth ( Vector& stress )
{
    double result;
    double alpha = ( 3.0 - sin ( mphiCS*PI/180.0 ) ) / ( 3.0 + sin ( mphiCS*PI/180.0 ) );
    double alpha4 = pow ( alpha, 4.0 );
    Vector s ( 6 );
    noalias ( s ) = GetDevStress ( stress );
//     double q = Getq ( stress );
//     double J3 = s[0] * s[1] * s[2] + 2.0 * s[3] * s[4] * s[5] - s[0] * s[4] * s[4] - s[1] * s[5] * s[5] - s[2] * s[3] * s[3];
//     
//     double nom = 2.0 * alpha4; 
//     double denom = 1.0 + alpha4 + (1.0 - alpha4) * 13.5 * J3 / pow(q, 3.0);

    double J2 = 0.5 * (s[0] * s[0] + s[1] * s[1] + s[2] * s[2])+ s[3] * s[3] + s[4] * s[4] + s[5] * s[5];
    double J3 = s[0]*s[1]*s[2] + 2.0*s[3]*s[4]*s[5] - s[0]*s[4]*s[4] - s[1]*s[5]*s[5] - s[2]*s[3]*s[3];


    double sin3ThL = J3 / sqrt(J2*J2*J2);
    sin3ThL *=  - 3.0 * sqrt(3.0) / 2.0;
     
    double nom = 2.0 * alpha4; 
    double denom = 1.0 + alpha4 + (1.0 - alpha4) * sin3ThL;
    
    result = mMmax * pow ( nom / denom, 0.25 );

    return result;
//     return mMmax;
}

///

bool ClayAndSandModel::Has ( const Variable<double>& rThisVariable )
{
    if ( rThisVariable == PRECONSOLIDATION ) return true;

    if ( rThisVariable == EQUIVALENT_VOLUMETRIC_STRAIN ) return true;

    if ( rThisVariable == EQUIVALENT_DEVIATORIC_STRAIN ) return true;

    if ( rThisVariable == EQUIVALENT_VOLUMETRIC_STRESS ) return true;

    if ( rThisVariable == EQUIVALENT_DEVIATORIC_STRESS ) return true;
    
    if ( rThisVariable == LOG_EQUIVALENT_VOLUMETRIC_STRESS ) return true;

    return false;
}

bool ClayAndSandModel::Has ( const Variable<Vector>& rThisVariable )
{
    if ( rThisVariable == INTERNAL_VARIABLES ) return true;

    if ( rThisVariable == INSITU_STRESS ) return true;

    return false;
}

bool ClayAndSandModel::Has ( const Variable<Matrix>& rThisVariable )
{
    return Kratos::ConstitutiveLaw::Has ( rThisVariable );
}

bool ClayAndSandModel::Has ( const Kratos::Variable< Kratos::array_1d< double, 3 > >& rThisVariable )
{
    return Kratos::ConstitutiveLaw::Has ( rThisVariable );
}

bool ClayAndSandModel::Has ( const Kratos::Variable< Kratos::array_1d< double, 6 > >& rThisVariable )
{
    return Kratos::ConstitutiveLaw::Has ( rThisVariable );
}

double& ClayAndSandModel::GetValue ( const Variable<double>& rThisVariable,
                                     double& rValue )
{
    if ( rThisVariable == PRECONSOLIDATION )
        rValue = mCurrentPreconsolidation;

    if ( rThisVariable == EQUIVALENT_VOLUMETRIC_STRAIN ) //epsilon_p
        rValue = -( mCurrentStrain[0] + mCurrentStrain[1] + mCurrentStrain[2] ); 

    if ( rThisVariable == EQUIVALENT_DEVIATORIC_STRAIN ) //epsilon_q
    { 
//         rValue = sqrt( 2.0*( pow(mCurrentStrain[0]-mCurrentStrain[1], 2.0) + pow(mCurrentStrain[1]-mCurrentStrain[2], 2.0) + pow(mCurrentStrain[2]-mCurrentStrain[0], 2.0) ) + 3.0*( pow(mCurrentStrain[3]/2.0,2.0)+pow(mCurrentStrain[4]/2.0,2.0)+pow(mCurrentStrain[5]/2.0,2.0) ) )/3.0;
        rValue = sqrt( 2.0*( pow(mCurrentStrain[0]-mCurrentStrain[1], 2.0) + pow(mCurrentStrain[1]-mCurrentStrain[2], 2.0) + pow(mCurrentStrain[2]-mCurrentStrain[0], 2.0) ) + 3.0*( pow(mCurrentStrain[3],2.0)+pow(mCurrentStrain[4],2.0)+pow(mCurrentStrain[5],2.0) ) )/3.0;
    }

    if ( rThisVariable == EQUIVALENT_VOLUMETRIC_STRESS ) //p
        rValue = Getp ( mCurrentStress );

    if ( rThisVariable == EQUIVALENT_DEVIATORIC_STRESS ) //q
        rValue = Getq ( mCurrentStress );

    if ( rThisVariable == LOG_EQUIVALENT_VOLUMETRIC_STRESS ) //ln(p)
        rValue = log10( Getp ( mCurrentStress ) );

    return rValue;
}

Vector& ClayAndSandModel::GetValue ( const Variable<Vector>& rThisVariable,
                                     Vector& rValue )
{
    if ( rThisVariable == INTERNAL_VARIABLES )
        rValue = mCurrentStrain;

    if ( rThisVariable == INSITU_STRESS )
        rValue = mCurrentStress;

    return rValue;
}

Matrix& ClayAndSandModel::GetValue ( const Variable<Matrix>& rThisVariable,
                                     Matrix& rValue )
{
    return rValue;
}

Kratos::array_1d< double, 3 >& ClayAndSandModel::GetValue ( const Kratos::Variable< Kratos::array_1d< double, 3 > >& rVariable, Kratos::array_1d< double, 3 >& rValue )
{
    return rValue;
}

Kratos::array_1d< double, 6 >& ClayAndSandModel::GetValue ( const Kratos::Variable< Kratos::array_1d< double, 6 > >& rVariable, Kratos::array_1d< double, 6 >& rValue )
{
    return rValue;
}

void ClayAndSandModel::SetValue ( const Variable<double>& rThisVariable,
                                  const double& rValue, const ProcessInfo& rCurrentProcessInfo )
{
    if ( rThisVariable == PRECONSOLIDATION )
        mOldPreconsolidation = rValue;
}

void ClayAndSandModel::SetValue ( const Variable<Vector>& rThisVariable,
                                  const Vector& rValue, const ProcessInfo& rCurrentProcessInfo )
{
}

void ClayAndSandModel::SetValue ( const Variable<Matrix>& rThisVariable,
                                  const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo )
{
}

void ClayAndSandModel::SetValue ( const Variable<array_1d<double, 3> >& rThisVariable, const array_1d <double, 3 > & rValue,
                                  const ProcessInfo& rCurrentProcessInfo )
{
}

void ClayAndSandModel::SetValue ( const Kratos::Variable< Kratos::array_1d< double, 6 > >& rVariable, const Kratos::array_1d< double, 6 >& Value, const Kratos::ProcessInfo& rCurrentProcessInfo )
{
}

bool ClayAndSandModel::ValidateInput ( const Kratos::Properties& props )
{
    return false;
}

int ClayAndSandModel::Check ( const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo )
{
    if ( ValidateInput ( props ) )
        return 0;

    return -1;
}

} // Namespace Kratos
