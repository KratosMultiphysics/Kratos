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
#include "constitutive_laws/casm.h"

#include "includes/constitutive_law.h"

#include "utilities/math_utils.h"
#include "includes/variables.h"
#include "includes/process_info.h"
#include "structural_application.h"
#include "custom_utilities/sd_math_utils.h"
#include "includes/properties.h"
#include "casm.h"
#include <math.h>
namespace Kratos
{

/**
 *	TO BE TESTED!!!
 */
Casm::Casm()
        : ConstitutiveLaw()
{
}
/**
 *	TO BE TESTED!!!
 */
Casm::~Casm()
{
}

bool Casm::Has( const Variable<double>& rThisVariable )
{
    if ( rThisVariable == PLASTICITY_INDICATOR )
        return true;
    return false;
}

bool Casm::Has( const Variable<Vector>& rThisVariable )
{
    if ( rThisVariable == INSITU_STRESS )
        return true;
    return false;
}

bool Casm::Has( const Variable<Matrix>& rThisVariable )
{
    return false;
}

double& Casm::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    if ( rThisVariable == PLASTICITY_INDICATOR )
    {
//             KRATOS_WATCH(mIplastic);
//             KRATOS_WATCH(mCurrentHistoryVariables[1]);
//             KRATOS_WATCH(this->mIplastic);
        rValue = (double)mIplastic;
    }
    else if ( rThisVariable == START_TIME )
        rValue = 1.0;
    else if ( rThisVariable == DP_EPSILON )
        rValue = mCurrentHistoryVariables[0];
    return( rValue );
}

Vector& Casm::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    if ( rThisVariable == INSITU_STRESS )
    {
        rValue = ZeroVector(6);
        rValue[0]= mInSituStress[0];
        rValue[1]= mInSituStress[1];
        rValue[2]= mInSituStress[2];
        rValue[3]= mInSituStress[3];
        rValue[4]= mInSituStress[5];
        rValue[5]= mInSituStress[4];
    }
    if ( rThisVariable == INTERNAL_VARIABLES )
    {
//            double J= sqrt(1.0/6.0*((mStress[0]-mStress[1])*(mStress[0]-mStress[1])+(mStress[1]-mStress[2])*(mStress[1]-mStress[2])+                        (mStress[2]-mStress[0])*(mStress[2]-mStress[0]))+mStress[3]*mStress[3]+mStress[4]*mStress[4]+mStress[5]*mStress[5]);
//
//            rValue = ZeroVector(8);
//            rValue[0] = -1.0/3.0*(mStress[0]+mStress[1]+mStress[2]);
//            rValue[1] = J;
//            rValue[2] = 0.0;//OCR
//            rValue[3] = 0.0;//Stress Ratio
//            rValue[4] = 0.0;
//            rValue[5] = mIplastic;    //plastic indicator
//            //Iplastic=0 -> no plasticity
//            //Iplastic=1 -> plasticity on the cap
//            //Iplastic=2 -> plasticity on the cap and on the shear yield surface
//            //Iplastic=3 -> plasticity on the shear yield surface;
//            rValue[6] = 0.0;
//            rValue[7] = 0;
    }
    if ( rThisVariable == MATERIAL_PARAMETERS )
    {
        rValue = mMaterialParameters;
    }

    return rValue;
}

Matrix& Casm::GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue )
{
    if (rThisVariable == PK2_STRESS_TENSOR)
    {
        rValue.resize(3,3);
        rValue(0,0)= mStress(0);
        rValue(0,1)= mStress(3);
        rValue(0,2)= mStress(5);

        rValue(1,0)= mStress(3);
        rValue(1,1)= mStress(1);
        rValue(1,2)= mStress(4);

        rValue(2,0)= mStress(5);
        rValue(2,1)= mStress(4);
        rValue(2,2)= mStress(2);
    }
    if (rThisVariable == GREEN_LAGRANGE_STRAIN_TENSOR)
    {
        rValue.resize(3,3);
        rValue(0,0)= mCurrentStrain[0];
        rValue(0,1)= mCurrentStrain[3];
        rValue(0,2)= mCurrentStrain[5];

        rValue(1,0)= mCurrentStrain[3];
        rValue(1,1)= mCurrentStrain[1];
        rValue(1,2)= mCurrentStrain[4];

        rValue(2,0)= mCurrentStrain[5];
        rValue(2,1)= mCurrentStrain[4];
        rValue(2,2)= mCurrentStrain[2];
    }
    if (rThisVariable == GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR)
    {
        rValue.resize(3,3);
        rValue(0,0)=mCurrentPlasticStrains[0];
        rValue(0,1)=mCurrentPlasticStrains[3];
        rValue(0,2)=mCurrentPlasticStrains[5];

        rValue(1,0)=mCurrentPlasticStrains[3];
        rValue(1,1)=mCurrentPlasticStrains[1];
        rValue(1,2)=mCurrentPlasticStrains[4];

        rValue(2,0)=mCurrentPlasticStrains[5];
        rValue(2,1)=mCurrentPlasticStrains[4];
        rValue(2,2)=mCurrentPlasticStrains[2];
    }
    return rValue;
}

void Casm::SetValue( const Variable<double>& rThisVariable,
                     const double& rValue,
                     const ProcessInfo& rCurrentProcessInfo )
{
    if (rThisVariable == INSITU_STRESS_SCALE)
    {
        noalias(mInSituStress) = rValue*mInSituStressOriginal;
    }

    if (rThisVariable == SUCTION)
    {
        mCurrentSuction= rValue/1000.0;//from Pa -> kPa
    }
}

void Casm::SetValue( const Variable<array_1d<double, 3> >& rThisVariable,
                     const array_1d<double, 3>& rValue,
                     const ProcessInfo& rCurrentProcessInfo
                   )
{
}

void Casm::SetValue( const Variable<Vector>& rThisVariable,
                     const Vector& rValue, const ProcessInfo& rCurrentProcessInfo
                   )
{
    if ( rThisVariable == INSITU_STRESS )
    {
        noalias(mInSituStress) = rValue;
        noalias(mInSituStressOriginal) = rValue;
    }

    if ( rThisVariable == MATERIAL_PARAMETERS )
    {
        mMaterialParameters = rValue;
    }
}

void Casm::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                     const ProcessInfo& rCurrentProcessInfo )
{
}

/**
 *	TO BE TESTED!!!
 */

/**
 * Material Parameters for default model 1:
 *  0.- Poisson's ratio ------------------- (    poiss)
 *  1.- Slope of the unload/reload curve--- (   xkappa)
 *  2.- Slope of the normal comp curve----- (  xlambda)
 *  3.- Spacing ratio --------------------- (        r)
 *  4.- Shape parameter of the yiel funct-- (       xn)
 *  5.- Slope of the CSL------------------- (       xM)
 *  6.- Friction angle at CS (degrees) -----(      phi)
 *  7.- Initial void ratio ---------------- (       e0)
 *  8.- Initial Preconsolidation pressure-- (       Po)
 *  9.- Undrained shear strenght ---------- (       Cu)
 */

void Casm::InitializeMaterial( const Properties& props,
                               const GeometryType& geom,
                               const Vector& ShapeFunctionsValues )
{
    mStress.resize(6);
    mOldStress.resize(6);
    mOldStress = ZeroVector(6);
    mCtangent = ZeroMatrix(6,6);
    mOldStrain.resize(6);
    mOldStrain = ZeroVector(6);
    mCurrentStrain.resize(6);
    mCurrentStrain = ZeroVector(6);
    mStrain.resize(6);
    mStrain = ZeroVector(6);
    mCurrentStrainInc.resize(6);
    mCurrentStrainInc = ZeroVector(6);
    mCurrentStress.resize(6);
    mCurrentStress = ZeroVector(6);

    mMaterialParameters = props[MATERIAL_PARAMETERS];
//         ResetMaterial(props, geom, ShapeFunctionsValues);
    for (unsigned int i=0; i<6; i++)
    {
        mCurrentHistoryVariables[i]= 0.0;
        mOldHistoryVariables[i]= 0.0;
    }

//         TESTING: flexible material parameters for soil model iMod=1
    for ( int i=0; i<6; i++ )
    {
        mOldStress[i] = 0.0;            // previous effective stress components (size=6)
        mOldStrain[i] = 0.0;
        mCurrentStress[i] = 0.0;
        mCurrentStrainInc[i] = 0.0;
        mCurrentStrain[i] = 0.0;
    }
    mInSituStressOriginal.resize(6);
    mInSituStress.resize(6);
    mInSituStress = ZeroVector(6);
    mInSituStressOriginal = ZeroVector(6);
    mConsistentTangent.resize(6,6);
    mConsistentTangent = ZeroMatrix(6,6);
//         mInSituStrain = ZeroVector(6);

    mCurrentTime = 0.0;

    mCurrentSuction= 0.0;
    mOldSuction= 0.0;
//         KRATOS_WATCH("RESETTING MATERIAL AND mIplastic")
    mIplastic = 0;

    mModelData[0] = mMaterialParameters[0];               //Young modulus(E)
    mModelData[1] = mMaterialParameters[1];               //poisson ratio (mue)
    mModelData[2] = mMaterialParameters[2];               //slope of CSL (M)
    mModelData[3] = mMaterialParameters[3];               //Friction Angle (phi_cs [degrees])
    //Yield function shape
    mModelData[4] = mMaterialParameters[4];                //slope of normal compression curve (lambda)
    mModelData[5] = mMaterialParameters[5];                //slope of unload/reload curve (kappa)
    //Strength parameters
    mModelData[6] = mMaterialParameters[6];              //initial void ratio (e_0)
    mModelData[7] = mMaterialParameters[7];              //shape parameter of yield function (n)
    //Initial state of the soil
    mModelData[8] = mMaterialParameters[8];                //spacing ratio (r)
    mModelData[9] = mMaterialParameters[9];           //preconsolidation pressure (P0 [kPa])

    for (unsigned i = 0; i<11 ; i++)
        std::cout << " Look Here =============== " << mModelData[i] << std::endl;

    isYielded = false;

    phiCs = (mModelData[3])* PI / 180;

    theta = (1 + mModelData[6] / (mModelData[4] - mModelData[5]));

    lnr = log(mModelData[8]);

    CalculateElasticMatrix(mCtangent, mModelData[0], mModelData[1]);

    CalculateUnit4thSym3D();

    CalculateUnit2nd3D();



}

void Casm::ResetMaterial(const Properties& props, const GeometryType& geom, const Vector& ShapeFunctionsValues )
{

}

void Casm::InitializeSolutionStep( const Properties& props,
                                   const GeometryType& geom, //this is just to give the array of nodes
                                   const Vector& ShapeFunctionsValues ,
                                   const ProcessInfo& CurrentProcessInfo)
{
}

void Casm::FinalizeSolutionStep( const Properties& props,
                                 const GeometryType& geom, //this is just to give the array of nodes
                                 const Vector& ShapeFunctionsValues ,
                                 const ProcessInfo& CurrentProcessInfo)
{
//    //store old strains and history variables
//     KRATOS_WATCH(CurrentProcessInfo[CALCULATE_INSITU_STRESS]);
//     if ( CurrentProcessInfo[CALCULATE_INSITU_STRESS] )
//     {
//         for (unsigned int i =0; i<6; i++)
//         {
//             mOldStrain[i] = 0.0;
//             mOldStress[i] = 0.0;
//             mInSituStress(i) =  (-1.0)*mCurrentStress[i];
//         }
//     }
//     else
//     {
//         for ( int i=0; i<6; i++ )
//         {
//             mOldStrain[i] = mCurrentStrain[i];
//             mOldStress[i] = mCurrentStress[i];
//         }
//     }
//

    mOldStress = mCurrentStress;
    mOldStrain = mCurrentStrain;

}

void Casm::CalculateMaterialResponse( const Vector& StrainVector,
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
    if ( CalculateStresses )
        CalculateStress( StrainVector, StressVector );
    if ( CalculateTangent != 0 )
        CalculateConstitutiveMatrix(StrainVector, AlgorithmicTangent);
}

/**
 *	TO BE TESTED!!!
 */
/**
 * Performing all necessary calculations inside an external library
 * Arguments:
 * I/O  Type
 * iMod        I   I    : model number (1-4)
 * IsUndr      I   I    : =1 for undrained, 0 otherwise
 * iIncr       I   I    : Global loading increment
 * Nvh        I/O  I    : Number of total history variables
 * Npar       I/O  I    : Number of total parameters
 * Par         I   R()  : List with model parameters
 * Sig0        I   R()  : Stresses at start of step
 * Swp0        I   R    : Excess pore pressure start of step
 * Suction0    I   R    : Previous suction value
 * Suction     I   R    : Suction value to be applied
 * His_var0    I  R()   : State variable at start of step
 * dEps        I   R()  : Strain increment
 * D           O   R(,) : Material stiffness matrix
 * BulkW       O   R    : Bulk modulus for water (undrained only)
 * Sig         O   R()  : Resulting stresses
 * Swp         O   R    : Resulting excess pore pressure
 * His_var     O   R()  : Resulting state variables
 * iplastic    O   I    : Plasticity indicator
 * Ierror_code O   I    : to force stopping of calculation
 */

/**
 *	TO BE TESTED!!!
 */
void Casm::CalculateStress(const Vector& StrainVector, Vector& StressVector)
{

    //////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////


    if ( StressVector.size() != 6 )
    {
        StressVector.resize(6);
    }

    double hydrostatic_pressure = (mOldStress(0) + mOldStress(1) + mOldStress(2)) / 3.0;

    if (hydrostatic_pressure == 0)
    {
        hydrostatic_pressure = mModelData[9];
        pck = mModelData[9];
    }
    else if (hydrostatic_pressure < 1.0)
        hydrostatic_pressure = 1.0;

    KRATOS_WATCH(hydrostatic_pressure);


    CalculateModules(hydrostatic_pressure, Km, Gm);
    pressureDependentElasticTangent(mCtangent);
    KRATOS_WATCH(mModelData[9]);
    KRATOS_WATCH(Km);
    KRATOS_WATCH(Gm);


    KRATOS_WATCH(mStrain);


    KRATOS_WATCH(StrainVector);

    mCurrentStrain = StrainVector;
    mCurrentStrainInc = StrainVector - mStrain;

    KRATOS_WATCH(mCurrentStrainInc);

    KRATOS_WATCH(mStress);

    mStrain = StrainVector;
    noalias(StressVector) = prod(mCtangent,mCurrentStrainInc);
    KRATOS_WATCH(StressVector);

    for (unsigned int i=0; i<6 ; i++)
        StressVector(i)+=mStress(i);

    KRATOS_WATCH(StressVector);

    mCurrentStress = -StressVector;

    KRATOS_WATCH(mCurrentStress);
    KRATOS_WATCH(mOldStress);



    if (mOldStress(0)==0 && mOldStress(1)==0 && mOldStress(2)==0 && mOldStress(3)==0 && mOldStress(4)==0 && mOldStress(5)==0 )
        {}
    else
    {

        // noalias(mStress)= StressVector;

        double pTr = (mCurrentStress(0) + mCurrentStress(1) + mCurrentStress(2)) / 3.0;
        double firstInv= mCurrentStress(0)+mCurrentStress(1)+mCurrentStress(2);  


        mMth = 0;

        Vector deviatoric_Stress_Vector(6);
        deviatoric_Stress_Vector(0)= mCurrentStress(0)-pTr;
        deviatoric_Stress_Vector(1)= mCurrentStress(1)-pTr;
        deviatoric_Stress_Vector(2)= mCurrentStress(2)-pTr;
        deviatoric_Stress_Vector(3)= mCurrentStress(3);
        deviatoric_Stress_Vector(4)= mCurrentStress(4);
        deviatoric_Stress_Vector(5)= mCurrentStress(5);
        KRATOS_WATCH(deviatoric_Stress_Vector);

        getM(deviatoric_Stress_Vector, mMth);
        KRATOS_WATCH(mMth);

        double j_2 = 0.5 * (deviatoric_Stress_Vector(0) * deviatoric_Stress_Vector(0)
                            + deviatoric_Stress_Vector(1) * deviatoric_Stress_Vector(1)
                            + deviatoric_Stress_Vector(2) * deviatoric_Stress_Vector(2))
                     + deviatoric_Stress_Vector(3) * deviatoric_Stress_Vector(3)
                     + deviatoric_Stress_Vector(4) * deviatoric_Stress_Vector(4)
                     + deviatoric_Stress_Vector(5) * deviatoric_Stress_Vector(5);
        double qTr = sqrt(3.0 * j_2);

//     std::cout << " Look Here Pc=============== " << mModelData[9] << std::endl;
        std::cout << " Look Here pTr=============== " << pTr << std::endl;

        double dummy_1 = pow(qTr / (mMth * pTr), mModelData[7]);
//     KRATOS_WATCH(dummy_1);
        double dummy_2 = log(pTr / mModelData[9]) / lnr ;
//     KRATOS_WATCH(dummy_2);
        double yield_function = pow(qTr / (mMth * pTr), mModelData[7]) + log(pTr / mModelData[9]) / lnr ;
        KRATOS_WATCH(yield_function);

//         if (yield_function > 1e-5)
//             isYielded = true;
//         else
//             isYielded = false;
        isYielded = pow(qTr / (mMth * pTr), mModelData[7]) + log(pTr / mModelData[9]) / lnr > 1e-5;
        std::cout << " Look Here you are inside 1=============== " << std::endl;

        if (isYielded) {
            std::cout << " Look Here you are inside 2=============== " << std::endl;
            returnMapping(pTr, qTr);

            double factor;
            if (qTr < 1e-5)
                factor = 0.0;
            else
                factor = qk / qTr;

            mCurrentStress(0) = factor * deviatoric_Stress_Vector(0) + pk;
            mCurrentStress(1) = factor * deviatoric_Stress_Vector(1) + pk;
            mCurrentStress(2) = factor * deviatoric_Stress_Vector(2) + pk;
            mCurrentStress(3) = factor * deviatoric_Stress_Vector(3);
            mCurrentStress(4) = factor * deviatoric_Stress_Vector(4);
            mCurrentStress(5) = factor * deviatoric_Stress_Vector(5);



            double factor_1;
            if (qTr < 1e-5)
                factor_1 = 0.0;
            else
                factor_1 = 1.5 / qTr;

            Vector fact(5);
            Vector nDev(6);
            for (int i = 0; i < 6; i++)
                nDev(i) = deviatoric_Stress_Vector(i) * factor_1;

            getFactors(qTr, fact);

            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                    mConsistentTangent(i,j) = fact(0) * unit4thSym3D[i][j] +
                                              fact(1) * unit2nd3D[i] * unit2nd3D[j] +
                                              fact(2) * nDev[i] * nDev[j] +
                                              fact(3) * nDev[i] * unit2nd3D[j] +
                                              fact(4) * unit2nd3D[i] * nDev[j];


            StressVector = -mCurrentStress;

        }

    }

    mStress = StressVector;

}

void Casm::Calculate( const Variable<Vector >& rVariable,
                      Vector& rResult, const ProcessInfo& rCurrentProcessInfo )
{
}

/**
 *	TO BE REVIEWED!!!
 */
void Casm::CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult)
{

    if (isYielded)
        rResult = mConsistentTangent;
    else
        rResult = mCtangent;

}

Matrix Casm::InverseC(Matrix& InvC)
{
    if (InvC.size1()!=6 || InvC.size2()!=6)
        InvC.resize(6,6);

    noalias(InvC)= ZeroMatrix(6,6);

    double lambda= mMaterialParameters[1]*mMaterialParameters[0]/((1+mMaterialParameters[1])*(1-2*mMaterialParameters[1]));
    double mu= mMaterialParameters[0]/(2*(1+mMaterialParameters[1]));

    double a= (4*mu*mu+4*mu*lambda)/(8*mu*mu*mu+12*mu*mu*lambda);
    double b= -(2*mu*lambda)/(8*mu*mu*mu+12*mu*mu*lambda);

    InvC(0,0)= a;
    InvC(0,1)= b;
    InvC(0,2)= b;
    InvC(1,0)= b;
    InvC(1,1)= a;
    InvC(1,2)= b;
    InvC(2,0)= b;
    InvC(2,1)= b;
    InvC(2,2)= a;
    InvC(3,3)= 1/(2*mu);
    InvC(4,4)= 1/(2*mu);
    InvC(5,5)= 1/(2*mu);
    return InvC;
}


void Casm::CalculateStressAndTangentMatrix( Matrix& StressTensor,
        const Matrix& StrainTensor,
        MaterialTensorType& algorithmicTangent )
{
}

void Casm::CalculateStressAndTangentMatrix( Vector& StressVector,
        const Vector& StrainVector,
        Matrix& algorithmicTangent )
{
}
void Casm::returnMapping(double pTr, double qTr) {
    Vector dG(5);
    Vector R(2);
    Matrix dR(2,2);

    pk = pTr;
    qk = qTr;

    KRATOS_WATCH(pk);
    KRATOS_WATCH(qk);

    dEpspk = dEpsqk = 0.0;

    formPlasticPotentialDerivatives(dG);

    KRATOS_WATCH(dG);
    formResidualVector(R, dG);

    KRATOS_WATCH(pck);
    KRATOS_WATCH(R);
    double normR = (R(0)*R(0)+R(1)*R(1));
    normR = sqrt(normR);
    while (normR > 1E-5) {

        formResidualDerivativeMatrix(dR, dG);
        KRATOS_WATCH(dR);

        solveByGauss(dR, R);

        KRATOS_WATCH(R);
        dEpspk -= R(0);
        dEpsqk -= R(1);

        pk = pTr - Km * dEpspk;
        qk = qTr - 3.0 * Gm * dEpsqk;
        pck = mModelData[9] * exp(theta * dEpspk);

//        if (pk < 1E-13) {
//            System.err.println("retmap:: forming residual pk < 0 ");
//            System.exit(0);
//        }

        formPlasticPotentialDerivatives(dG);
        formResidualVector(R, dG);
        normR = (R(0)*R(0)+R(1)*R(1));
    }

}

void Casm::CalculateElasticMatrix(Matrix& C, const double E, const double NU)
{
    //setting up material matrix
    double c1 = E / ((1.00+NU)*(1-2*NU));
    double c2 = c1 * (1-NU);
    double c3 = c1 * NU;
    double c4 = c1 * 0.5 * (1 - 2*NU);
    //filling material matrix
    C(0,0) = c2;
    C(0,1) = c3;
    C(0,2) = c3;
    C(0,3) = 0.0;
    C(0,4) = 0.0;
    C(0,5) = 0.0;
    C(1,0) = c3;
    C(1,1) = c2;
    C(1,2) = c3;
    C(1,3) = 0.0;
    C(1,4) = 0.0;
    C(1,5) = 0.0;
    C(2,0) = c3;
    C(2,1) = c3;
    C(2,2) = c2;
    C(2,3) = 0.0;
    C(2,4) = 0.0;
    C(2,5) = 0.0;
    C(3,0) = 0.0;
    C(3,1) = 0.0;
    C(3,2) = 0.0;
    C(3,3) = c4;
    C(3,4) = 0.0;
    C(3,5) = 0.0;
    C(4,0) = 0.0;
    C(4,1) = 0.0;
    C(4,2) = 0.0;
    C(4,3) = 0.0;
    C(4,4) = c4;
    C(4,5) = 0.0;
    C(5,0) = 0.0;
    C(5,1) = 0.0;
    C(5,2) = 0.0;
    C(5,3) = 0.0;
    C(5,4) = 0.0;
    C(5,5) = c4;

}

void Casm::pressureDependentElasticTangent(Matrix& mCtangent) {

    double C1 = 2.0 * Gm;
    double C2 = Km - 2.0 * Gm / 3.0;

//         KRATOS_WATCH(C1);
//         KRATOS_WATCH(C2);

    for (unsigned int i = 0; i < 6; i++)
        for (unsigned int j = i; j < 6; j++)
            mCtangent(i,j) = C1 * unit4thSym3D[i][j] + C2 * unit2nd3D[i] * unit2nd3D[j];

    for (unsigned int j = 0; j < 5; j++)
        for (unsigned int i = j+1; i < 6; i++)
            mCtangent(i,j) = mCtangent(j,i);

}

void Casm::CalculateUnit4thSym3D()
{
    for (unsigned int i = 0; i < 6; i++)
        for (unsigned int j = 0; j < 6; j++)
            unit4thSym3D[i][j] = 0;

    unit4thSym3D[0][0] = 1.0;
    unit4thSym3D[1][1] = 1.0;
    unit4thSym3D[2][2] = 1.0;
    unit4thSym3D[3][3] = 0.5;
    unit4thSym3D[4][4] = 0.5;
    unit4thSym3D[5][5] = 0.5;

}

void Casm::CalculateUnit2nd3D()
{
    unit2nd3D[0] = 1.0;
    unit2nd3D[1] = 1.0;
    unit2nd3D[2] = 1.0;
    unit2nd3D[3] = 0.0;
    unit2nd3D[4] = 0.0;
    unit2nd3D[5] = 0.0;


}

void Casm::formPlasticPotentialDerivatives(Vector& dG) {

    double N1 = 3.0 + 2.0 * mMth;
    double D1 = 3.0 * pk + 2.0 * qk;

    double N2 = 3.0 - mMth;
    double D2 = 3.0 * pk - qk;

    dG(0) = 3.0 * (N1 / D1 - N2 / D2);
    dG(1) = 2.0 * N1 / D1 + N2 / D2;

    double DD2 = D2 * D2;
    double DD1 = D1 * D1;

    dG(2) = 9.0 * N2 / DD2 - 9.0 * N1 / DD1;
    dG(3) = N2 / DD2 - 4.0 * N1 / DD1;
    dG(4) = -3.0 * N2 / DD2 - 6.0 * N1 / DD1;
}

void Casm::formYieldFunctionDerivatives(Vector& dF) {

    dF(0) = pow(pk, mModelData[7] - 1) * (mModelData[7] * log(pk / pck) + 1.0) / lnr;
    dF(1) = mModelData[7] * pow(qk / mMth, mModelData[7] - 1) / mMth;
    dF(2) = - pow(pk, mModelData[7]) / (pck * lnr);

}

void Casm::formResidualDerivativeMatrix(Matrix& dR, Vector& dG)
{
    // dG = (dG/dp, dG/dq, d2G/dp2, d2G/dq2, d2G/dpdq)

    dR(0,0) = dG(1) - dEpspk * dG(4) * Km + dEpsqk * dG(2) * Km;
    dR(0,1) = -dG(0) + (dEpspk * dG(3) - dEpsqk * dG(4)) * (-3.0 * Gm);

    Vector dF(3);
    formYieldFunctionDerivatives(dF);

    dR(1,0) = dF(0) * (-Km) + dF(2) * pck * theta;
    dR(1,1) = dF(1) * (-3.0 * Gm);


}

void Casm::formResidualVector(Vector& R, Vector& dG)
{
    R(0) = dEpspk * dG(1) - dEpsqk * dG(0);

    double tmp1 = pow(qk / mMth, mModelData[7]);
    double tmp2 = pow(pk, mModelData[7]);
    double tmp3 = log(pk / pck);

    KRATOS_WATCH(tmp1);
    KRATOS_WATCH(tmp2);
    KRATOS_WATCH(tmp3);

    R(1) = tmp1 + tmp2 * tmp3 / lnr;
}

void Casm::solveByGauss(Matrix& A, Vector& B) {
    double sum;

    unsigned int N = B.size();

    //triangular decomposition
    decompose(A, N);

    //forward
    for (unsigned int i = 0; i < N; i++) {
        sum = 0.0;
        for (unsigned int j = 0; j < i; j++)
            sum += A(i,j) * B(j);

        B(i) -= sum;
    }

    //backward
    for (int i = N - 1; i > -1; i--) {
        sum = 0.0;
        for (int j = N - 1; j > i; j--)
            sum += A(i,j) * B(j);

        B(i) = (B(i) - sum) / A(i,i);
    }
}//solveGauss

void Casm::decompose(Matrix& A, int N) {
    for (int k = 0; k < N; k++) {

        for (int i = k + 1; i < N; i++)
            for (int j = k + 1; j < N; j++)
                A(i,j) = A(i,j) - A(i,k) * A(k,j) / A(k,k);

        for (int i = k + 1; i < N; i++)
            A(i,k) = A(i,k) / A(k,k);
    }
}

void Casm::getM( const Vector& devStr, double& rMth) {

    double J2, J3;


    J2 = 0.5 * (devStr[0] * devStr[0] + devStr[1] * devStr[1] + devStr[2] * devStr[2])	+
         devStr[3] * devStr[3] + devStr[4] * devStr[4] + devStr[5] * devStr[5];
    J3 = devStr[0]*devStr[1]*devStr[2] + 2.0*devStr[3]*devStr[4]*devStr[5] -
         devStr[0]*devStr[4]*devStr[4] - devStr[1]*devStr[5]*devStr[5] - devStr[2]*devStr[3]*devStr[3];


    double sin3ThL = J3 / sqrt(J2*J2*J2);
    sin3ThL *=  - 3.0 * sqrt(3.0) / 2.0;
    KRATOS_WATCH(sin3ThL);
    double alpha = (3.0 - sin(phiCs)) / (3.0 + sin(phiCs));
    double alpha4 = pow(alpha, 4);
    KRATOS_WATCH(alpha4);
    double nom = 2.0 * alpha4;
    KRATOS_WATCH(nom);
    double denom = 1.0 + alpha4 + (1.0 - alpha4) * sin3ThL;

    std::cout << " Look Here mModelData[2]=============== " <<  mModelData[2] << std::endl;

    rMth = mModelData[2] * pow((nom / denom), 0.25);
    KRATOS_WATCH(rMth);
}

void Casm::getFactors(double qTr, Vector& factor) {

    Vector dF(3);
    formYieldFunctionDerivatives(dF);

    Vector dG(5);
    formPlasticPotentialDerivatives(dG);

    double A11 = dF(0) * Km - dF(2) * pck * theta;
    double A12 = dF(1) * 3.0 * Gm;
    double A21 = dG(1) + Km * (-dEpspk * dG(4) + dEpsqk * dG(2));
    double A22 = -dG(0) + 3.0 * Gm * (-dEpspk * dG(3) + dEpsqk * dG(4));

    double B11 = dF(0) * Km;
    double B12 = dF(1) * 2.0 * Gm;
    double B21 = Km * (-dEpspk * dG(4) + dEpsqk * dG(2));
    double B22 = 2.0 * Gm * (-dEpspk * dG(3) + dEpsqk * dG(4));

    double Det = A11 * A22 - A21 * A12;

    double mpI = (A22 * B11 - A12 * B21) / Det;
    double mpN = (A22 * B12 - A12 * B22) / Det;
    double mqI = (A11 * B21 - A21 * B11) / Det;
    double mqN = (A11 * B22 - A21 * B12) / Det;


    factor(0) = 2.0 * Gm * qk / qTr;
    factor(1) = Km * (1.0 - mpI) - 2.0 * Gm * qk / qTr / 3.05;
    factor(2) = 4.0 * Gm / 3.0 * (1.0 - qk / qTr - 3.0 * mqN / 2.0);
    factor(3) = - 2.0 * Gm * mqI;
    factor(4) = - Km * mpN;

}

void Casm::CalculateModules(const double& p, double& rKm, double& rGm) {

    //double modulus[2];
    Km = (1.0 + mModelData[6]) * p / mModelData[5];
    Gm = 3.0 * Km * (1.0 - 2.0 * mModelData[1]) / 2.0 / (1.0 + mModelData[1]);
}


//**********************************************************************
//**********************************************************************
void Casm::CalculateCauchyStresses(
    Vector& rCauchy_StressVector,
    const Matrix& rF,
    const Vector& rPK2_StressVector,
    const Vector& rGreenLagrangeStrainVector)
{
}
} // Namespace Kratos
