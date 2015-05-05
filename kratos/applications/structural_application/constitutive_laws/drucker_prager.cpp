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
#include "constitutive_laws/drucker_prager.h"

#include "includes/constitutive_law.h"

#include "utilities/math_utils.h"
#include "includes/variables.h"
#include "includes/process_info.h"
#include "structural_application.h"
#include "custom_utilities/sd_math_utils.h"
#include "includes/properties.h"
#include <math.h>
#include "drucker_prager.h"

namespace Kratos
{

/**
 * TO BE TESTED!!!
 */
DruckerPrager::DruckerPrager() :
        ConstitutiveLaw()
{
}

/**
 * TO BE TESTED!!!
 */
DruckerPrager::~DruckerPrager()
{
}

bool DruckerPrager::Has( const Variable<double>& rThisVariable )
{
    if ( rThisVariable == PLASTICITY_INDICATOR ) return true;

    if ( rThisVariable == EQUIVALENT_VOLUMETRIC_PLASTIC_STRAIN ) return true;

    if ( rThisVariable == EQUIVALENT_DEVIATORIC_PLASTIC_STRAIN ) return true;

    return false;
}

bool DruckerPrager::Has( const Variable<Vector>& rThisVariable )
{
    if ( rThisVariable == INSITU_STRESS ) return true;

    if ( rThisVariable == PRESTRESS ) return true;

    if ( rThisVariable == STRESSES ) return true;

    if ( rThisVariable == PLASTIC_STRAIN_VECTOR ) return true;

    return false;
}

bool DruckerPrager::Has( const Variable<Matrix>& rThisVariable )
{
    return Kratos::ConstitutiveLaw::Has( rThisVariable );
}

bool DruckerPrager::Has( const Kratos::Variable< Kratos::array_1d< double, 3 > >& rThisVariable )
{
    return Kratos::ConstitutiveLaw::Has( rThisVariable );
}

bool DruckerPrager::Has( const Kratos::Variable< Kratos::array_1d< double, 6 > >& rThisVariable )
{
    return Kratos::ConstitutiveLaw::Has( rThisVariable );
}


/**
 * TODO: check definitions of equivalent strains!!!
 */
double& DruckerPrager::GetValue( const Variable<double>& rThisVariable,
                                 double& rValue )
{
    if ( rThisVariable == PLASTICITY_INDICATOR )
    {
        //rValue = norm_2( mCurrentPlasticStrains );
        rValue = mAlpha;
    }

    if ( rThisVariable == EQUIVALENT_VOLUMETRIC_PLASTIC_STRAIN )
    {
        rValue = ( mCurrentPlasticStrains[0] + mCurrentPlasticStrains[1] + mCurrentPlasticStrains[2] ) / 3.0;
    }

    if ( rThisVariable == EQUIVALENT_DEVIATORIC_PLASTIC_STRAIN )
    {
        Vector deviatoric_plastic_strains = mCurrentPlasticStrains;
        double volumetric_plastic_strains = (( mCurrentPlasticStrains[0] + mCurrentPlasticStrains[1] + mCurrentPlasticStrains[2] ) / 3.0 );

        for ( unsigned int i = 0; i < 3; i++ )
            deviatoric_plastic_strains[i] -= volumetric_plastic_strains;

        rValue = sqrt( 2.0 / 3.0 * inner_prod( deviatoric_plastic_strains, deviatoric_plastic_strains ) );
    }

    return( rValue );
}

Vector& DruckerPrager::GetValue( const Variable<Vector>& rThisVariable,
                                 Vector& rValue )
{
    if ( rThisVariable == INSITU_STRESS || rThisVariable == PRESTRESS )
    {
        rValue = mPrestress;
    }

    if ( rThisVariable == INTERNAL_VARIABLES )
    {
        //INTERNAL VARIABLES:
        // 0 - Plastic strain xx
        // 1 - Plastic strain yy
        // 2 - Plastic strain zz
        // 3 - Plastic strain xy
        // 4 - Plastic strain yz
        // 5 - Plastic strain xz
        // 6 - hardening alpha
        rValue = ZeroVector( 7 );

        for ( unsigned int i = 0; i < 6; i++ )
            rValue[i] = mCurrentPlasticStrains[i];

        rValue[6] = mAlpha;
    }

    if ( rThisVariable == STRESSES )
    {
        rValue = mCurrentStress;
    }

    if ( rThisVariable == PLASTIC_STRAIN_VECTOR )
    {
        rValue = mCurrentPlasticStrains;
    }

    return rValue;
}

Matrix& DruckerPrager::GetValue( const Variable<Matrix>& rThisVariable,
                                 Matrix& rValue )
{
    return rValue;
}

Kratos::array_1d< double, 3 >& DruckerPrager::GetValue( const Kratos::Variable< Kratos::array_1d< double, 3 > >& rVariable, Kratos::array_1d< double, 3 >& rValue )
{
    return rValue;
}

Kratos::array_1d< double, 6 >& DruckerPrager::GetValue( const Kratos::Variable< Kratos::array_1d< double, 6 > >& rVariable, Kratos::array_1d< double, 6 >& rValue )
{
    return rValue;
}

void DruckerPrager::SetValue( const Variable<double>& rThisVariable,
                              const double& rValue, const ProcessInfo& rCurrentProcessInfo )
{
    if ( rThisVariable == PRESTRESS_FACTOR )
        mPrestressFactor = rValue;
    if ( rThisVariable == YOUNG_MODULUS )
        mE = rValue;
    if ( rThisVariable == POISSON_RATIO )
        mNU = rValue;
    if ( rThisVariable == COHESION )
        mCohesion = rValue;
    if ( rThisVariable == ISOTROPIC_HARDENING_MODULUS )
        mHardening = rValue;
    if ( rThisVariable == INTERNAL_FRICTION_ANGLE )
    {
        double tan_phi = tan( rValue*PI/180 );
        mEta = 3.0 * tan_phi / ( sqrt( 9.0 + 12.0 * tan_phi * tan_phi ) );
        mXi = 3.0 / ( sqrt( 9.0 + 12.0 * tan_phi * tan_phi ) );
    }
}

void DruckerPrager::SetValue( const Variable<Vector>& rThisVariable,
                              const Vector& rValue, const ProcessInfo& rCurrentProcessInfo )
{
    if ( rThisVariable == PRESTRESS || rThisVariable == INSITU_STRESS )
    {
        noalias( mPrestress ) = rValue;
    }

    if ( rThisVariable == INTERNAL_VARIABLES )
    {
        //INTERNAL VARIABLES:
        // 0 - Plastic strain xx
        // 1 - Plastic strain yy
        // 2 - Plastic strain zz
        // 3 - Plastic strain xy
        // 4 - Plastic strain yz
        // 5 - Plastic strain xz
        // 6 - hardening alpha
        for ( unsigned int i = 0; i < 6; i++ )
            mCurrentPlasticStrains[i] = rValue[i];

        mAlpha = rValue[6];
    }
}

void DruckerPrager::SetValue( const Variable<Matrix>& rThisVariable,
                              const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo )
{
}

void DruckerPrager::SetValue( const Variable<array_1d<double, 3> >& rThisVariable, const array_1d <double, 3 > & rValue,
                              const ProcessInfo& rCurrentProcessInfo )
{
}


void DruckerPrager::SetValue( const Kratos::Variable< Kratos::array_1d< double, 6 > >& rVariable, const Kratos::array_1d< double, 6 >& Value, const Kratos::ProcessInfo& rCurrentProcessInfo )
{
}

/**
 * Input paramters for DruckerPrager:
 * YOUNG_MODULUS
 * POISSON_RATIO
 * COHESION
 * INTERNAL_FRICTION_ANGLE
 * ISOTROPIC_HARDENING_MODULUS
 */
bool DruckerPrager::ValidateInput( const Kratos::Properties& props )
{
    if ( props.Has( YOUNG_MODULUS ) && props.Has( POISSON_RATIO ) && props.Has( COHESION ) && props.Has( INTERNAL_FRICTION_ANGLE ) && props.Has( ISOTROPIC_HARDENING_MODULUS ) )
    {
        if ( props.GetValue( YOUNG_MODULUS ) > 0.0 && props.GetValue( POISSON_RATIO ) < 0.5 &&  props.GetValue( ISOTROPIC_HARDENING_MODULUS ) >= 0.0 && props.GetValue( INTERNAL_FRICTION_ANGLE ) > 0.0 && props.GetValue( INTERNAL_FRICTION_ANGLE ) < 90.0 && props.GetValue( COHESION ) >= 0.0 )
        {
            return true;
        }
    }

    return false;
}

int DruckerPrager::Check( const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo )
{
    if ( ValidateInput( props ) )
        return 0;
    return -1;
}


void DruckerPrager::InitializeMaterial( const Properties& props,
                                        const GeometryType& geom, const Vector& ShapeFunctionsValues )
{
    mOldStress.resize( 6 );
    mOldStress = ZeroVector( 6 );
    mCtangent.resize(6, 6);
    mCtangent = ZeroMatrix( 6, 6 );
    mCtangentInv.resize(6, 6);
    mCtangentInv = ZeroMatrix(6, 6);
    mOldStrain.resize( 6 );
    mOldStrain = ZeroVector( 6 );
    mCurrentStrain.resize( 6 );
    mCurrentStrain = ZeroVector( 6 );
//     mCurrentStrainInc.resize( 6 );
//     mCurrentStrainInc = ZeroVector( 6 );
    mCurrentStress.resize( 6 );
    mCurrentStress = ZeroVector( 6 );
    mCurrentPlasticStrains.resize( 6 );
    mCurrentPlasticStrains = ZeroVector( 6 );
    mOldSElasticStrain.resize( 6 );
    mOldSElasticStrain = ZeroVector( 6 );
    mCurrentElasticStrain.resize( 6 );
    mCurrentElasticStrain = ZeroVector( 6 );
    mOldPlasticStrains.resize( 6 );
    mOldPlasticStrains = ZeroVector( 6 );
    mPrestress.resize( 6 );
    mPrestress = ZeroVector( 6 );
    mPrestressFactor = 1.0;

//     mdGamma = 0.0;
    mAlpha = 0.0;
    mOldAlpha = 0.0;
//     mIsYielded = false;
//     mIsApex = false;

    mE = props[YOUNG_MODULUS];
    mNU = props[POISSON_RATIO];
    mCohesion = props[COHESION];
    mHardening = props[ISOTROPIC_HARDENING_MODULUS];

    double tan_phi = tan(( props[INTERNAL_FRICTION_ANGLE] ) * PI / 180 );
    mEta = 3.0 * tan_phi / ( sqrt( 9.0 + 12.0 * tan_phi * tan_phi ) );
    mXi = 3.0 / ( sqrt( 9.0 + 12.0 * tan_phi * tan_phi ) );

    
//     ResetMaterial( props, geom, ShapeFunctionsValues );

//     mG = mE / ( 2.0 * ( 1.0 + mNU ) );
//     mK = mE / ( 3.0 * ( 1.0 - 2.0 * mNU ) );

//     CalculateUnit4thSym3D();

//     CalculateUnit2nd3D();
    
    ResetMaterial( props, geom, ShapeFunctionsValues );

}

void DruckerPrager::ResetMaterial( const Properties& props,
                                   const GeometryType& geom, const Vector& ShapeFunctionsValues )
{
    mG = mE / ( 2.0 * ( 1.0 + mNU ) );
    mK = mE / ( 3.0 * ( 1.0 - 2.0 * mNU ) );
// 
    CalculateElasticMatrix( mCtangent, mE, mNU );
    mPrestressFactor = 1.0;
    mPrestress = ZeroVector( 6 );
// 
//     CalculateElasticMatrix( mCtangent, mE, mNU );
// 
//     CalculateUnit4thSym3D();
// 
//     CalculateUnit2nd3D();

}

void DruckerPrager::InitializeSolutionStep(
    const Properties& props,
    const GeometryType& geom, //this is just to give the array of nodes
    const Vector& ShapeFunctionsValues,
    const ProcessInfo& CurrentProcessInfo )
{
}

void DruckerPrager::InitializeNonLinearIteration( const Kratos::Properties& props, const Kratos::ConstitutiveLaw::GeometryType& geom, const Kratos::Vector& ShapeFunctionsValues, const Kratos::ProcessInfo& CurrentProcessInfo )
{
}


void DruckerPrager::FinalizeSolutionStep(
    const Properties& props,
    const GeometryType& geom, //this is just to give the array of nodes
    const Vector& ShapeFunctionsValues,
    const ProcessInfo& CurrentProcessInfo )
{

    mOldStress = mCurrentStress;
    mOldStrain = mCurrentStrain;
    mOldAlpha = mAlpha;
    mOldPlasticStrains = mCurrentPlasticStrains;
    mOldSElasticStrain = mCurrentElasticStrain;

}

void DruckerPrager::CalculateMaterialResponse( const Vector& StrainVector,
        const Matrix& DeformationGradient, Vector& StressVector,
        Matrix& AlgorithmicTangent, const ProcessInfo& CurrentProcessInfo,
        const Properties& props, const GeometryType& geom,
        const Vector& ShapeFunctionsValues, bool CalculateStresses,
        int CalculateTangent, bool SaveInternalVariables )
{    
    bool isYielded = false;
    bool isApex = false;
    double dGamma = 0.0;
    if ( CalculateStresses )
        CalculateStress( StrainVector, StressVector, isYielded, isApex, dGamma );

    if ( CalculateTangent != 0 )
        CalculateConstitutiveMatrix( StrainVector, AlgorithmicTangent, isYielded, isApex, dGamma );
}

/**
 * TO BE TESTED!!!
 */
void DruckerPrager::CalculateStress( const Vector& StrainVector, Vector& StressVector, bool& isYielded, bool& isApex, double& dGamma )
{
    mCurrentStrain = StrainVector;
    mCurrentElasticStrain = StrainVector;
    Vector CurrentStrainInc = StrainVector - mOldStrain;

    
//     StressVector = prod( mCtangent, StrainVector )  - mPrestressFactor * mPrestress;
    StressVector = prod( mCtangent, CurrentStrainInc);
    for ( unsigned int i = 0; i < 6; i++ )
        StressVector( i ) += mOldStress( i );

    mCurrentStress = StressVector - mPrestressFactor * mPrestress;

    isYielded = false;
    isApex = false;
    dGamma = 0.0;
    mAlpha = mOldAlpha;

    double pTr = ( mCurrentStress[0] + mCurrentStress[1] + mCurrentStress[2] ) / 3.0;

    Vector deviatoric_Stress_Vector( 6 );

    deviatoric_Stress_Vector[0] = mCurrentStress[0] - pTr;

    deviatoric_Stress_Vector[1] = mCurrentStress[1] - pTr;

    deviatoric_Stress_Vector[2] = mCurrentStress[2] - pTr;

    deviatoric_Stress_Vector[3] = mCurrentStress[3];

    deviatoric_Stress_Vector[4] = mCurrentStress[4];

    deviatoric_Stress_Vector[5] = mCurrentStress[5];

    double J2 = 0.5 * ( deviatoric_Stress_Vector[ 0 ]
                        * deviatoric_Stress_Vector[ 0 ] + deviatoric_Stress_Vector[ 1 ]
                        * deviatoric_Stress_Vector[ 1 ] + deviatoric_Stress_Vector[ 2 ]
                        * deviatoric_Stress_Vector[ 2 ] ) + deviatoric_Stress_Vector[ 3 ]
                * deviatoric_Stress_Vector[ 3 ] + deviatoric_Stress_Vector[ 4 ]
                * deviatoric_Stress_Vector[ 4 ] + deviatoric_Stress_Vector[ 5 ]
                * deviatoric_Stress_Vector[ 5 ];

    double sqrtJ2 = sqrt( J2 );

    double yield_function = sqrtJ2 + mEta * pTr - mXi * ( mCohesion + mHardening * mAlpha );

    double Tol = 1e-10;

    double dummy = 0.0;

    double dPlasticStrainVol = 0.0;

    double dDPlasticStrainVol = 0.0;

    double pTr_n = 0.0;

    double residual = 0.0;

    double df_residual = 0.0;

    if ( yield_function > Tol )
    {
// 	KRATOS_WATCH( " ############ ########## ######## Plasticity ########## ########## #  # # #  #" );

      
        //Plastic regione is activated
        isYielded = true;

        // returnmap is needed

        for ( int itr = 0; itr < 50; itr++ )
        {
            double df_dGamma  = -mG - mK * mEta * mEta - mXi * mXi * mHardening;

            double dDGamma = -yield_function / df_dGamma;
            dGamma += dDGamma;

            mAlpha += mXi * dDGamma;
            pTr_n = pTr - mK * mEta * dGamma;

            yield_function = sqrtJ2 - mG * dGamma + mEta * pTr_n - mXi
                             * ( mCohesion + mHardening * mAlpha );
                             
            if( itr > 48 )
            {
//                 mdGamma = 0.0;
//                 mAlpha = mOldAlpha;
//                 break;
//                 KRATOS_THROW_ERROR( std::logic_error, "return map did not converge", "" );
            }

            if ( fabs( yield_function ) < Tol )
            {
                if (( sqrtJ2 - mG * dGamma ) > 0.0 )
                {
                    if ( fabs( sqrtJ2 ) < 1e-6 )
                        dummy = 0.0;
                    else
                        dummy = 1.0 - mG * dGamma / sqrtJ2;

                    break;
                }
                else {
		  
// 		              KRATOS_WATCH( " ############ ########## ######## Apex ########## ########## #  # # #  #" );
//                     KRATOS_THROW_ERROR( std::logic_error, "Kratos was killed because a material point in DP-model reached the APEX region", "" );

                    isApex = true;
                    dPlasticStrainVol = 0.0;
                    mAlpha = mOldAlpha;
                    residual = (mXi / mEta) * (mCohesion + mHardening * mAlpha) - pTr;

                    for (int itr_apex = 0; itr_apex < 40; itr_apex++) {
                        df_residual = (mXi / mEta) * (mXi / mEta) * mHardening + mK;
                        dDPlasticStrainVol = -residual / df_residual;
                        dPlasticStrainVol += dDPlasticStrainVol;

                        // compute new residual
                        mAlpha = mOldAlpha + (mXi / mEta) * dPlasticStrainVol;
                        pTr_n = pTr - mK * dPlasticStrainVol;
                        residual = (mXi / mEta)
                                   * (mCohesion + mHardening * mAlpha)
                                   - pTr_n;
                        if (fabs(residual) < Tol) {
                            dGamma = dPlasticStrainVol / mEta;
                            dummy = 0.0;
                            break;
                        }
//                         if( itr_apex > 2 )
//                         KRATOS_WATCH("##### ITERATING TOO MUCH IN APEX #####");
                    }
//                     if( itr > 2 )
//                         KRATOS_WATCH("##### ITERATING TOO MUCH #####");
                }
            }

        }

        mCurrentStress[0] = dummy * deviatoric_Stress_Vector[0] + pTr_n;
        mCurrentStress[1] = dummy * deviatoric_Stress_Vector[1] + pTr_n;
        mCurrentStress[2] = dummy * deviatoric_Stress_Vector[2] + pTr_n;
        mCurrentStress[3] = dummy * deviatoric_Stress_Vector[3];
        mCurrentStress[4] = dummy * deviatoric_Stress_Vector[4];
        mCurrentStress[5] = dummy * deviatoric_Stress_Vector[5];

        StressVector = mCurrentStress;
        mCurrentStress += (mPrestressFactor * mPrestress);

        mCurrentElasticStrain[0] = dummy / (2.0 * mG)
                                   * deviatoric_Stress_Vector[0] + pTr_n / (mK*3.0);
        mCurrentElasticStrain[1] = dummy / (2.0 * mG)
                                   * deviatoric_Stress_Vector[1] + pTr_n / (mK*3.0);
        mCurrentElasticStrain[2] = dummy / (2.0 * mG)
                                   * deviatoric_Stress_Vector[2] + pTr_n / (mK*3.0);
        mCurrentElasticStrain[3] = dummy / (2.0 * mG)
                                   * deviatoric_Stress_Vector[3] * 2.0;
        mCurrentElasticStrain[4] = dummy / (2.0 * mG)
                                   * deviatoric_Stress_Vector[4] * 2.0;
        mCurrentElasticStrain[5] = dummy / (2.0 * mG)
                                   * deviatoric_Stress_Vector[5] * 2.0;

    }
    else
    {

        StressVector = mCurrentStress;
        mCurrentStress += (mPrestressFactor * mPrestress);
    }
    mCurrentPlasticStrains = mOldPlasticStrains + CurrentStrainInc - (mCurrentElasticStrain - mOldSElasticStrain);
    
//     mCurrentStress += mPrestressFactor * mPrestress;

//     mCurrentPlasticStrains = -prod(InverseC(mCtangentInv),StressVector);

//     for (int i = 0 ; i < 3 ; i++)
//     {
//         mCurrentPlasticStrains[i] += StrainVector[i];
//     }
//     for (int i = 3 ; i < 6 ; i++)
//     {
//         mCurrentPlasticStrains[i] = 2*mCurrentPlasticStrains[i] + StrainVector[i];
//     }
    
    //     mCurrentPlasticStrains = StrainVector - prod(InverseC(mCtangentInv),StressVector);
//     mCurrentPlasticStrains[1] = StrainVector[1] - prod(InverseC(mCtangentInv),StressVector)(1);
//     mCurrentPlasticStrains[2] = StrainVector[2] - prod(InverseC(mCtangentInv),StressVector)(2);
//     mCurrentPlasticStrains[3] = StrainVector[3] - prod(InverseC(mCtangentInv),StressVector)(3)*2.0;
//     mCurrentPlasticStrains[4] = StrainVector[4] - prod(InverseC(mCtangentInv),StressVector)(4)*2.0;
//     mCurrentPlasticStrains[5] = StrainVector[5] - prod(InverseC(mCtangentInv),StressVector)(5)*2.0;


}


/**
 * TO BE REVIEWED!!!
 */
void DruckerPrager::CalculateConstitutiveMatrix( const Vector& StrainVector,
        Matrix& rResult, bool& isYielded, bool& isApex, double& dGamma )
{

//     Vector Strain_trial( 6 );
//     Strain_trial = StrainVector;
    
//     if ( false )
    if ( isYielded )
    {
        if (isApex) {
//             KRATOS_WATCH("############################# Apex Apex Apex Apex Apex ####################");

            rResult = mCtangent;
        }
        else
        {
            Matrix ConsistentTangent = ZeroMatrix(6,6);
            Matrix unit4thSym3D = ZeroMatrix(6,6);
            Vector unit2nd3D = ZeroVector(6);
            for( unsigned int i=0; i<3; i++ )
            {
                unit4thSym3D(i,i) = 1.0;
                unit4thSym3D(i+3,i+3) = 0.5; 
                unit2nd3D[i] = 1.0;
            }

            double StrainVol = ( StrainVector[0] + StrainVector[1] + StrainVector[2] ) / 3.0;

            Vector StrainDev( 6 );
            StrainDev[0] = StrainVector[0] - StrainVol;
            StrainDev[1] = StrainVector[1] - StrainVol;
            StrainDev[2] = StrainVector[2] - StrainVol;
            StrainDev[3] = StrainVector[3] / 2.0;
            StrainDev[4] = StrainVector[4] / 2.0;
            StrainDev[5] = StrainVector[5] / 2.0;

            double StrainDevNorm = sqrt( StrainDev[0] * StrainDev[0] + StrainDev[1]
                                         * StrainDev[1] + StrainDev[2] * StrainDev[2] + 2.0
                                         * ( StrainDev[3] * StrainDev[3] + StrainDev[4] * StrainDev[4]
                                             + StrainDev[5] * StrainDev[5] ) );

            Vector UnitStrainDev( 6 );

            if ( fabs( StrainDevNorm ) > 1e-8 )
            {
                UnitStrainDev = StrainDev / StrainDevNorm;
            }
            else
                UnitStrainDev = ZeroVector( 6 );

            double A_aux = 1.0 / ( mG + mK * mEta * mEta + mXi * mXi * mHardening );
//             double A_aux = 1.43014e-09;
            

            double B_aux = 2.0 * mG * ( 1.0 - dGamma / ( sqrt( 2.0 ) * StrainDevNorm ) );

            double C_aux = 2.0 * mG * ( dGamma / ( sqrt( 2.0 ) * StrainDevNorm ) - mG * A_aux );
            
            double D_aux = sqrt( 2.0 ) * mG * mK * A_aux;

            double E_aux = mK * ( 1 - mK * mEta * mEta * A_aux );

            for ( int i = 0; i < 6; i++ )
                for ( int j = 0; j < 6; j++ )
                    ConsistentTangent( i, j ) = B_aux * unit4thSym3D(i,j) + C_aux
                                                 * UnitStrainDev[i] * UnitStrainDev[j] - D_aux * ( mEta
                                                         * UnitStrainDev[i] * unit2nd3D[j] + mEta * unit2nd3D[i]
                                                         * UnitStrainDev[j] ) + ( E_aux - B_aux / 3.0 )
                                                 * unit2nd3D[i] * unit2nd3D[j];
                                                 
            for ( int i = 0; i < 6; i++ )
                for ( int j = 0; j < 6; j++ )
                    if( fabs(ConsistentTangent(i,j)) < 1.0e-7 )
                        ConsistentTangent(i,j) = 0.0;

             
            rResult = ConsistentTangent;
        }
    }
    else
    {
        rResult = mCtangent;
    }

}

Matrix DruckerPrager::InverseC(Matrix& InvC) {
    if (InvC.size1() != 6 || InvC.size2() != 6)
        InvC.resize(6, 6);

    noalias(InvC) = ZeroMatrix(6, 6);

    double lambda = mNU * mE / ((1
                    + mNU) * (1 - 2 * mNU));
    double mu = mE / (2 * (1 + mNU));

    double a = (4 * mu * mu + 4 * mu * lambda) / (8 * mu * mu * mu + 12 * mu
               * mu * lambda);
    double b = -(2 * mu * lambda) / (8 * mu * mu * mu + 12 * mu * mu * lambda);

    InvC(0, 0) = a;
    InvC(0, 1) = b;
    InvC(0, 2) = b;
    InvC(1, 0) = b;
    InvC(1, 1) = a;
    InvC(1, 2) = b;
    InvC(2, 0) = b;
    InvC(2, 1) = b;
    InvC(2, 2) = a;
    InvC(3, 3) = 1 / (2 * mu);
    InvC(4, 4) = 1 / (2 * mu);
    InvC(5, 5) = 1 / (2 * mu);
    return InvC;
}

void DruckerPrager::CalculateElasticMatrix( Matrix& C, const double E,
        const double NU )
{
    //setting up material matrix
    double c1 = E / (( 1.0 + NU ) * ( 1.0 - 2.0 * NU ) );
    double c2 = c1 * ( 1.0 - NU );
    double c3 = c1 * NU;
    double c4 = c1 * 0.5 * ( 1.0 - 2.0 * NU );
    //filling material matrix
    C( 0, 0 ) = c2;
    C( 0, 1 ) = c3;
    C( 0, 2 ) = c3;
    C( 0, 3 ) = 0.0;
    C( 0, 4 ) = 0.0;
    C( 0, 5 ) = 0.0;
    C( 1, 0 ) = c3;
    C( 1, 1 ) = c2;
    C( 1, 2 ) = c3;
    C( 1, 3 ) = 0.0;
    C( 1, 4 ) = 0.0;
    C( 1, 5 ) = 0.0;
    C( 2, 0 ) = c3;
    C( 2, 1 ) = c3;
    C( 2, 2 ) = c2;
    C( 2, 3 ) = 0.0;
    C( 2, 4 ) = 0.0;
    C( 2, 5 ) = 0.0;
    C( 3, 0 ) = 0.0;
    C( 3, 1 ) = 0.0;
    C( 3, 2 ) = 0.0;
    C( 3, 3 ) = c4;
    C( 3, 4 ) = 0.0;
    C( 3, 5 ) = 0.0;
    C( 4, 0 ) = 0.0;
    C( 4, 1 ) = 0.0;
    C( 4, 2 ) = 0.0;
    C( 4, 3 ) = 0.0;
    C( 4, 4 ) = c4;
    C( 4, 5 ) = 0.0;
    C( 5, 0 ) = 0.0;
    C( 5, 1 ) = 0.0;
    C( 5, 2 ) = 0.0;
    C( 5, 3 ) = 0.0;
    C( 5, 4 ) = 0.0;
    C( 5, 5 ) = c4;

}

// void DruckerPrager::CalculateUnit4thSym3D()
// {
//     for ( unsigned int i = 0; i < 6; i++ )
//         for ( unsigned int j = 0; j < 6; j++ )
//             unit4thSym3D[i][j] = 0;
// 
//     unit4thSym3D[0][0] = 1.0;
// 
//     unit4thSym3D[1][1] = 1.0;
// 
//     unit4thSym3D[2][2] = 1.0;
// 
//     unit4thSym3D[3][3] = 0.5;
// 
//     unit4thSym3D[4][4] = 0.5;
// 
//     unit4thSym3D[5][5] = 0.5;
//
//}

//void DruckerPrager::CalculateUnit2nd3D()
//{
//    unit2nd3D[0] = 1.0;
//    unit2nd3D[1] = 1.0;
//    unit2nd3D[2] = 1.0;
//    unit2nd3D[3] = 0.0;
//    unit2nd3D[4] = 0.0;
//    unit2nd3D[5] = 0.0;
//
//}


} // Namespace Kratos