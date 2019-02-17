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
*   Last Modified by:    $Author: hbui $
*   Date:                $Date: 9 Nov 2015 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/

// System includes
#include <iostream>
#include <cmath>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/process_info.h"
#include "includes/properties.h"
#include "utilities/math_utils.h"
#include "utilities/kratos_log.h"
#include "constitutive_laws/cam_clay_3d.h"
#include "structural_application.h"

// #define DEBUG_CAM_CLAY
// #define LOG_CAM_CLAY

namespace Kratos
{

#ifdef BOOST_NO_CXX11_CONSTEXPR
const double CamClay3D::TOL = 1.0e-5;
#endif

const double CamClay3D::unit2nd3D[6] = {1.0, 1.0, 1.0, 0.0, 0.0, 0.0};

const double CamClay3D::unit4thSym3D[][6] = {
                {1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
                {0.0, 0.0, 1.0, 0.0, 0.0, 0.0},
                {0.0, 0.0, 0.0, 0.5, 0.0, 0.0},
                {0.0, 0.0, 0.0, 0.0, 0.5, 0.0},
                {0.0, 0.0, 0.0, 0.0, 0.0, 0.5}
            };

/**
 * TO BE TESTED!!!
 */
CamClay3D::CamClay3D()
    : ConstitutiveLaw()
{
}

/**
 * TO BE TESTED!!!
 */
CamClay3D::~CamClay3D()
{
}

bool CamClay3D::Has( const Variable<int>& rThisVariable )
{
    return false;
}

bool CamClay3D::Has( const Variable<double>& rThisVariable )
{
    if ( rThisVariable == PRESTRESS_FACTOR || rThisVariable == YOUNG_MODULUS || rThisVariable == POISSON_RATIO )
        return true;

    return false;
}

bool CamClay3D::Has( const Variable<Vector>& rThisVariable )
{
    if ( rThisVariable == INSITU_STRESS )
        return true;

    if ( rThisVariable == PRESTRESS )
        return true;

    if ( rThisVariable == STRESSES )
        return true;

    if ( rThisVariable == CURRENT_STRAIN_VECTOR )
        return true;

    return false;
}

bool CamClay3D::Has( const Variable<Matrix>& rThisVariable )
{
    return false;
}


int& CamClay3D::GetValue( const Variable<int>& rThisVariable, int& rValue )
{
    return rValue;
}

double& CamClay3D::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    if ( rThisVariable == PRESTRESS_FACTOR ) {
        rValue = mPrestressFactor;
        return rValue;
    }

    if(rThisVariable == YOUNG_MODULUS ) {
       rValue = mE;
       return rValue;
    }

    if ( rThisVariable == POISSON_RATIO ) {
        rValue = mNU;
        return rValue;
    }

    if (rThisVariable == DELTA_TIME) {
        rValue = sqrt(mE/mDE);
        return rValue;
    }

    return rValue;
}

Vector& CamClay3D::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    if ( rThisVariable == INSITU_STRESS || rThisVariable == PRESTRESS )
    {
        const unsigned int size = mPrestress.size();
        rValue.resize(size, false);
        rValue = mPrestressFactor*mPrestress;
        return rValue;
    }

    if ( rThisVariable == STRESSES )
    {
        const unsigned int size = mCurrentStress.size();
        rValue.resize(size, false);
        rValue  = mCurrentStress;
        return rValue;
    }

    if ( rThisVariable == PLASTIC_STRAIN_VECTOR )
    {
        rValue = ZeroVector( 6 );
        return( rValue );
    }

    if ( rThisVariable == INTERNAL_VARIABLES )
    {
        rValue = ZeroVector( 1 );
        return( rValue );
    }

    KRATOS_THROW_ERROR( std::logic_error, "Vector Variable case not considered", "" );
}

Matrix& CamClay3D::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
        KRATOS_THROW_ERROR( std::logic_error, "Matrix Variable case not considered", "" );
}

void CamClay3D::SetValue( const Variable<int>& rThisVariable, const int& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    if(rThisVariable == PARENT_ELEMENT_ID)
    {
        mParentElementId = rValue;
    }
    if(rThisVariable == INTEGRATION_POINT_INDEX)
    {
        mIntegrationPointIndex = rValue;
    }
}

void CamClay3D::SetValue( const Variable<double>& rThisVariable, const double& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
    if ( rThisVariable == PRESTRESS_FACTOR )
        mPrestressFactor = rValue;
    if ( rThisVariable == YOUNG_MODULUS )
        mE = rValue;
    if ( rThisVariable == POISSON_RATIO )
        mNU = rValue;

    if ( rThisVariable == PRECONSOLIDATION_PRESSURE_MIN )
    {
        mPck = this->getP();
        if(mPck < rValue)
            mPck = rValue;

        mKm = (1.0 + mVoidRatio) * mPck / mKappa;
        mGm = 3.0 * mKm * (1.0 - 2.0 * mNU) / 2.0 / (1.0 + mNU);

        #ifdef DEBUG_CAM_CLAY
        std::cout << "At PRECONSOLIDATION_PRESSURE_MIN, material point elem = "
                  << mParentElementId << ", int = " << mIntegrationPointIndex << ":" << std::endl;
        std::cout << " mVoidRatio: " << mVoidRatio << std::endl;
        std::cout << " this->getP(): " << this->getP() << std::endl;
        std::cout << " mPck: " << mPck << std::endl;
        std::cout << " mKappa: " << mKappa << std::endl;
        std::cout << " mNU: " << mNU << std::endl;
        std::cout << " mKm: " << mKm << std::endl;
        std::cout << " mGm: " << mGm << std::endl;
        #endif
    }

    if ( rThisVariable == PRECONSOLIDATION_PRESSURE_DEF )
    {
        mPck = rValue;

        mKm = (1.0 + mVoidRatio) * mPck / mKappa;
        mGm = 3.0 * mKm * (1.0 - 2.0 * mNU) / 2.0 / (1.0 + mNU);

        #ifdef DEBUG_CAM_CLAY
        std::cout << "At PRECONSOLIDATION_PRESSURE_DEF, material point elem = "
                  << mParentElementId << ", int = " << mIntegrationPointIndex << ":" << std::endl;
        std::cout << " mVoidRatio: " << mVoidRatio << std::endl;
        std::cout << " mPck: " << mPck << std::endl;
        std::cout << " mKappa: " << mKappa << std::endl;
        std::cout << " mNU: " << mNU << std::endl;
        std::cout << " mKm: " << mKm << std::endl;
        std::cout << " mGm: " << mGm << std::endl;
        #endif
    }
}

void CamClay3D::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
        const ProcessInfo& rCurrentProcessInfo )
{
    if ( rThisVariable == INSITU_STRESS || rThisVariable == PRESTRESS )
    {
        noalias(mPrestress) = rValue;
    }
    if ( rThisVariable == STRESSES )
    {
        noalias(mCurrentStress) = rValue;
    }
    if ( rThisVariable == CURRENT_STRAIN_VECTOR )
    {
        noalias(mCurrentStrain) = rValue;
    }
}

void CamClay3D::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                            const ProcessInfo& rCurrentProcessInfo )
{
}

void CamClay3D::InitializeMaterial( const Properties& props,
        const GeometryType& geom, const Vector& ShapeFunctionsValues )
{
    mCurrentStress = ZeroVector( 6 );
    mLastStress = ZeroVector( 6 );
    mCurrentStrain = ZeroVector( 6 );
    mPrestress = ZeroVector( 6 );
    mLastStrain = ZeroVector( 6 );
    mLastStrainIncr = ZeroVector( 6 );
    mPrestressFactor = 1.0;
    mE = props[YOUNG_MODULUS];
    mNU = props[POISSON_RATIO];
    mDE = props[DENSITY];
    mM = props[CSL_SLOPE];
    mLambda = props[VIRGIN_COMPRESSION_INDEX];
    mKappa = props[SWELL_INDEX];
    mVoidRatio = props[VOID_RATIO];
    mTheta = (1 + mVoidRatio) / (mLambda - mKappa);
    mDGamma = 0.0;
    mIsYielded = false;
}

void CamClay3D::ResetMaterial( const Properties& props,
        const GeometryType& geom, const Vector& ShapeFunctionsValues )
{
    noalias(mCurrentStress) = ZeroVector( 6 );
    noalias(mLastStress) = ZeroVector( 6 );
    noalias(mCurrentStrain) = ZeroVector( 6 );
    noalias(mPrestress) = ZeroVector( 6 );
    noalias(mLastStrain) = ZeroVector( 6 );
    noalias(mLastStrainIncr) = ZeroVector( 6 );
//    mPrestressFactor = 1.0;
    mDGamma = 0.0;
    mIsYielded = false;
}

void CamClay3D::InitializeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mPc = mPck;
    mDGamma = 0.0;

    noalias(mLastStrainIncr) = mLastStrain - mCurrentStrain;

    #ifdef DEBUG_CAM_CLAY
    if(mParentElementId == 1 && mIntegrationPointIndex == 0)
    {
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "At InitializeSolutionStep, material point elem = "
                  << mParentElementId << ", int = " << mIntegrationPointIndex << std::endl;
        std::cout << " mCurrentStrain: " << mCurrentStrain << std::endl;
        std::cout << " mLastStrain: " << mLastStrain << std::endl;
        std::cout << " mLastStrainIncr: " << mLastStrainIncr << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
    }
    #endif
}

void CamClay3D::InitializeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    #ifdef DEBUG_CAM_CLAY
    if(mParentElementId == 1 && mIntegrationPointIndex == 0)
    {
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "At InitializeNonLinearIteration, material point elem = "
                  << mParentElementId << ", int = " << mIntegrationPointIndex << std::endl;
//        std::cout << " mCurrentStrain: " << mCurrentStrain << std::endl;
//        std::cout << " mLastStrain: " << mLastStrain << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
    }
    #endif
}

void CamClay3D::CalculateMaterialResponse( const Vector& StrainVector,
        const Matrix& DeformationGradient,
        Vector& StressVector,
        Matrix& AlgorithmicTangent,
        const ProcessInfo& CurrentProcessInfo,
        const Properties& props,
        const GeometryType& geom,
        const Vector& ShapeFunctionsValues,
        bool CalculateStresses,
        int CalculateTangent,
        bool SaveInternalVariables )
{
    // compute the strain increment
    Vector strainIncr = mLastStrainIncr;

    #ifdef DEBUG_CAM_CLAY
    if(mParentElementId == 1 && mIntegrationPointIndex == 0)
    {
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "At CalculateMaterialResponse, material point elem = "
                  << mParentElementId << ", int = " << mIntegrationPointIndex << std::endl;
//        std::cout << " StrainVector: " << StrainVector << std::endl;
//        std::cout << " mCurrentStrain: " << mCurrentStrain << std::endl;
//        std::cout << " mLastStrain: " << mLastStrain << std::endl;
        std::cout << " strainIncr: " << strainIncr << std::endl;
        std::cout << " mIsYielded: " << mIsYielded << std::endl;
    }
    #endif

    // compute the elasticity matrix
    Matrix Ce(6, 6);
    this->CalculateElasticTangent(Ce, mKm, mGm);

    ///////////////////////////////////////////////////////////////////////
    // compute the tangent from the last converged stress
    if(mIsYielded)
    {
        Vector stressTr = prod(Ce, strainIncr) + mLastStress;
        Vector devStressTr(6);
        devStressTr = this->getDeviatoricComp(devStressTr, stressTr);
        double devNorm = sqrt(2.0 * this->getJ2(devStressTr));

        double devS, invDevNorm;
        if (devNorm < TOL)
        {
            devS = 1.0;
            invDevNorm = 0.0;
        }
        else
        {
            devS = sqrt(2.0 / 3) * mQk / devNorm;
            invDevNorm = 1.0 / devNorm;
        }

        Vector unitDev(6);
        noalias(unitDev) = devStressTr * invDevNorm;

        Vector fact;
        this->getFactors(fact, devS);
        #ifdef DEBUG_CAM_CLAY
        if(mParentElementId == 1 && mIntegrationPointIndex == 0)
        {
            std::cout << " stressTr: " << stressTr << std::endl;
            std::cout << " devStressTr: " << devStressTr << std::endl;
            std::cout << " devNorm: " << devNorm << std::endl;
            std::cout << " invDevNorm: " << invDevNorm << std::endl;
            std::cout << " unitDev: " << unitDev << std::endl;
            std::cout << " devS: " << devS << std::endl;
            std::cout << " fact: " << fact << std::endl;
        }
        #endif
        for(int i = 0; i < 6; ++i)
            for(int j = 0; j < 6; ++j)
                AlgorithmicTangent(i, j) = fact[0] * unit4thSym3D[i][j] +
                                           fact[1] * unit2nd3D[i] * unit2nd3D[j] +
                                           fact[2] * unit2nd3D[i] * unitDev[j] +
                                           fact[3] * unitDev[i] * unit2nd3D[j] +
                                           fact[4] * unitDev[i] * unitDev[j];
    }
    else
    {
        noalias(AlgorithmicTangent) = Ce;
    }

    // export the stress
    noalias(StressVector) = -1.0 * mCurrentStress;

    #ifdef DEBUG_CAM_CLAY
    if(mParentElementId == 1 && mIntegrationPointIndex == 0)
    {
        std::cout << " mKm: " << mKm << std::endl;
        std::cout << " mGm: " << mGm << std::endl;
//        std::cout << " Ce: " << Ce << std::endl;
        std::cout << " mCurrentStress: " << mCurrentStress << std::endl;
        std::cout << " AlgorithmicTangent: " << AlgorithmicTangent << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
    }
    #endif
}

void CamClay3D::FinalizeNonLinearIteration( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    ///////////////////////////////////////////////////////////////////////
    // stres integration

    // compute the strain increment
    Vector strainIncr = mLastStrain - mCurrentStrain;
    noalias(mLastStrainIncr) = strainIncr;

    #ifdef DEBUG_CAM_CLAY
    if(mParentElementId == 1 && mIntegrationPointIndex == 0)
    {
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "At FinalizeNonLinearIteration, material point elem = "
                  << mParentElementId << ", int = " << mIntegrationPointIndex << std::endl;
//        std::cout << " mCurrentStrain: " << mCurrentStrain << std::endl;
//        std::cout << " mLastStrain: " << mLastStrain << std::endl;
        std::cout << " strainIncr: " << strainIncr << std::endl;
    }
    #endif

    // compute the elasticity matrix
    Matrix Ce(6, 6);
    this->CalculateElasticTangent(Ce, mKm, mGm);

    Vector stressTr = prod(Ce, strainIncr) + mLastStress + mPrestressFactor * mPrestress;
    double pTr = (stressTr[0] + stressTr[1] + stressTr[2]) / 3;
    Vector devStressTr(6);
    devStressTr = this->getDeviatoricComp(devStressTr, stressTr);
    double J2 = this->getJ2( devStressTr );
    double qTr = sqrt(3.0 * J2);

    mIsYielded = qTr * qTr / (mM * mM) + pTr * (pTr - mPc) > TOL;

    if(mIsYielded)
    {
        this->returnMapping(pTr, qTr);

        double devNorm = sqrt(2.0 * J2);

        double factor;
        if( devNorm < TOL )
            factor = 0.0;
        else
            factor = mQk * sqrt(2.0 / 3) / devNorm;

        mCurrentStress[0] = factor * devStressTr[0] + mPk;
        mCurrentStress[1] = factor * devStressTr[1] + mPk;
        mCurrentStress[2] = factor * devStressTr[2] + mPk;
        mCurrentStress[3] = factor * devStressTr[3];
        mCurrentStress[4] = factor * devStressTr[4];
        mCurrentStress[5] = factor * devStressTr[5];
    }
    else
    {
        mDGamma = 0.0;
        noalias(mCurrentStress) = stressTr;
    }
    noalias(mCurrentStress) -= mPrestressFactor * mPrestress;

    #ifdef DEBUG_CAM_CLAY
    if(mParentElementId == 1 && mIntegrationPointIndex == 0)
    {
        std::cout << " mCurrentStress: " << mCurrentStress << std::endl;
        std::cout << " mIsYielded: " << mIsYielded << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
    }
    #endif

    #ifdef LOG_CAM_CLAY
    if(mIsYielded)
    {
        std::cout << "Element " << mParentElementId
                  << " at integration point " << mIntegrationPointIndex
                  << " is yielded" << std::endl;
    }
    #endif
}

void CamClay3D::FinalizeSolutionStep( const Properties& props,
        const GeometryType& geom, //this is just to give the array of nodes
        const Vector& ShapeFunctionsValues ,
        const ProcessInfo& CurrentProcessInfo )
{
    mKm = (1.0 + mVoidRatio) * this->getP() / mKappa;
    mGm = 3.0 * mKm * (1.0 - 2.0 * mNU) / 2.0 / (1.0 + mNU);

    noalias(mLastStrain) = mCurrentStrain;
    noalias(mLastStress) = mCurrentStress;

    #ifdef DEBUG_CAM_CLAY
    if(mParentElementId == 1 && mIntegrationPointIndex == 0)
    {
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "At FinalizeSolutionStep, material point elem = "
                  << mParentElementId << ", int = " << mIntegrationPointIndex << std::endl;
//        std::cout << " mCurrentStrain: " << mCurrentStrain << std::endl;
//        std::cout << " mLastStrain: " << mLastStrain << std::endl;
//        std::cout << " mKm: " << mKm << std::endl;
//        std::cout << " mGm: " << mGm << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
    }
    #endif

}

//**********************************************************************
int CamClay3D::Check( const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY

    if ( !props.Has( YOUNG_MODULUS ) || !props.Has(POISSON_RATIO) )
    {
        KRATOS_THROW_ERROR( std::logic_error, "this constitutive law requires YOUNG_MODULUS and POISSON_RATIO given as KRATOS variables", "" );
    }

    double nu = props[POISSON_RATIO];

    if ( nu > 0.499 && nu < 0.501 )
    {
        KRATOS_THROW_ERROR( std::logic_error, "invalid poisson ratio in input, close to incompressibility", "" );
        return -1;
    }

    return 0;

    KRATOS_CATCH( "" );
}

double CamClay3D::getP()
{
    double result = (mCurrentStress[0] + mCurrentStress[1] + mCurrentStress[2]) / 3;
//    if (result < 1.0)
//    {
//        result = 1.0;
//        KRATOS_THROW_ERROR(std::logic_error, "getP return value < 1.0", "")
//    }
    return result;
}

double CamClay3D::solveG(double pTr)
{
    double result = mPc;
    double dGp;

    double expFact = exp( mTheta * mDGamma * ( 2 * pTr - result ) / ( 1 + 2 * mDGamma * mKm ) );
    double Gp = mPc * expFact - result;

    int it = 0, max_it = 30;
    while(fabs( Gp ) > TOL && it < max_it)
    {
        dGp = mPc * expFact * ( -mTheta * mDGamma / ( 1 + 2 * mDGamma * mKm ) ) - 1.0;
        result -= Gp / dGp;

        expFact = exp( mTheta * mDGamma * ( 2 * pTr - result ) / ( 1 + 2 * mDGamma * mKm ) );
        Gp = mPc * expFact - result;

        ++it;
    }

    if(it == max_it)
    {
        KRATOS_WATCH(expFact)
        KRATOS_WATCH(mTheta)
        KRATOS_WATCH(mDGamma)
        KRATOS_WATCH(pTr)
        KRATOS_WATCH(mKm)
        KRATOS_WATCH(mPc)
        KRATOS_WATCH(Gp)
        KRATOS_WATCH(dGp)
        KRATOS_WATCH(result)
        KRATOS_WATCH(mParentElementId)
        KRATOS_WATCH(mIntegrationPointIndex)
        KRATOS_THROW_ERROR(std::runtime_error, "solveG does not converge in 30 steps", "")
    }

    return result;
}

void CamClay3D::getFactors(Vector& rResults, double devS)
{
    double a = 1.0 + 2.0 * mKm * mDGamma + mPck * mTheta * mDGamma;
    double a1 = ( 1.0 + mPck * mTheta * mDGamma ) / a;
    double a2 = -( 2.0 * mPk - mPck ) / a;
    double a3 = 2.0 * mPck * mTheta * mDGamma / a;
    double a4 = mTheta * mPck * (2.0 * mPk - mPck) / (a * mKm);
    double denom = 1.0 + 6.0 * mGm * mDGamma / (mM * mM);
    double a5 = sqrt(3.0 / 2) / denom;
    double a6 = -3.0 * mQk / ( denom * mM * mM );

    double b = -2.0 * mGm * 2.0 * mQk * a6 / (mM * mM) - mKm * ((2.0 * a2 - a4) * mPk - a2 * mPck );
    double b1 = -mKm * ( (a3 - 2.0 * a1) * mPk + a1 * mPck ) / b;
    double b2 = 2.0 * mGm * 2.0 * mQk * a5 / (b * mM * mM);

    if(rResults.size() != 5)
        rResults.resize(5);

    rResults[0] = 2.0 * mGm * devS;
    rResults[1] = mKm * (a1 + a2 * b1) - 2.0 * mGm * devS / 3;
    rResults[2] = mKm * a2 * b2;
    rResults[3] = 2.0 * mGm * sqrt(2.0 / 3) * (a6 * b1);
    rResults[4] = 2.0 * mGm * (sqrt(2.0 / 3) * (a5 + a6 * b2) - devS);
}

void CamClay3D::calculateF(std::vector<double>& rResults, double pTr, double qTr)
{
    mQk = qTr / ( 1 + 6 * mGm * mDGamma / (mM * mM) );
    mPck = this->solveG(pTr);
    mPk = (pTr + mDGamma * mKm * mPck) / (1 + 2 * mDGamma * mKm);

    double dFdPk  = 2.0 * mPk - mPck;
    double dFdQk  = 2.0 * mQk / (mM * mM);
    double dFdPck = -mPk;

    double dpdDGa  = -mKm * (2.0 * mPk - mPck) / ( 1.0 + (2.0 * mKm + mTheta * mPck) * mDGamma );
    double dqdDGa  = -mQk / (mDGamma + mM * mM / (6.0 * mGm));
    double dpcdDGa = mTheta * mPck * (2.0 * mPk - mPck) / ( 1.0 + (2.0 * mKm + mTheta * mPck) * mDGamma );

    if(rResults.size() != 2)
        rResults.resize(2);

    rResults[0] = mQk * mQk / (mM * mM) + mPk * (mPk - mPck);
    rResults[1] = dFdPk * dpdDGa + dFdQk * dqdDGa + dFdPck * dpcdDGa;
}

void CamClay3D::returnMapping(double pTr, double qTr)
{
    mDGamma = 0.0;
    std::vector<double> F(2);

    this->calculateF(F, pTr, qTr);
    int it = 0, max_it = 30;
    while( fabs( F[0] / pow(fabs(pTr) + fabs(qTr), 2) ) > TOL && it < max_it )
    {
        mDGamma -= F[0] / F[1];
        this->calculateF(F, pTr, qTr);
        ++it;
    }
    if(it == max_it)
    {
        KRATOS_WATCH(F[0])
        KRATOS_WATCH(F[1])
        KRATOS_WATCH(pTr)
        KRATOS_WATCH(qTr)
        KRATOS_THROW_ERROR(std::runtime_error, "returnMapping does not converge in 30 steps", "")
    }
}

void CamClay3D::CalculateElasticTangent(Matrix& C, double K, double G)
{
    double C1 = 2.0 * G;
    double C2 = K - 2.0 * G / 3;

    if(C.size1() != 6 || C.size2() != 6)
        C.resize(6, 6);

    for(int i = 0; i < 6; ++i)
        for(int j = i; j < 6; ++j)
            C(i, j) = C1 * unit4thSym3D[i][j] + C2 * unit2nd3D[i] * unit2nd3D[j];//p.622 De = 2G Is + (K - 2/3 G) I x I

    for(int j = 0; j < 5; ++j)
        for(int i = j+1; i < 6; ++i)
            C(i, j) = C(j, i);
}

Vector& CamClay3D::getDeviatoricComp(Vector& devStress, const Vector& stress)
{
    double p = ( stress[0] + stress[1] + stress[2] ) / 3;
    devStress[0] = stress[0] - p;
    devStress[1] = stress[1] - p;
    devStress[2] = stress[2] - p;
    devStress[3] = stress[3];
    devStress[4] = stress[4];
    devStress[5] = stress[5];

    return devStress;
}

double CamClay3D::getJ2(const Vector& s)
{
    return 0.5 * (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]) + s[3] * s[3] + s[4] * s[4] + s[5] * s[5];
}

} // Namespace Kratos

#undef DEBUG_CAM_CLAY
#undef LOG_CAM_CLAY
