//  KRATOS  ____       _               _
//         / ___|  ___(_)___ _ __ ___ (_) ___
//         \___ \ / _ \ / __| '_ ` _ \| |/ __|
//          ___) |  __/ \__ \ | | | | | | (__
//         |____/ \___|_|___/_| |_| |_|_|\___|
//
//  License:     BSD License
//  license:     structural_mechanics_application/license.txt
//
//  Main authors: Mahmoud Zidan
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/uniaxial_menegotto_pinto.h"
#include "seismic_application_variables.h"

namespace Kratos {


ConstitutiveLaw::Pointer UniaxialMenegottoPintoMaterialLaw::Clone() const
{
    return Kratos::make_shared<UniaxialMenegottoPintoMaterialLaw>(*this);
}

void UniaxialMenegottoPintoMaterialLaw::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mStrainSize = 1;
    rFeatures.mSpaceDimension = 3;
}

void UniaxialMenegottoPintoMaterialLaw::InitializeMaterialResponsePK2(Parameters& rValues)
{
    KRATOS_TRY
    mTangentModulus = rValues.GetMaterialProperties()[STEEL_YOUNGS_MODULUS];
    KRATOS_CATCH("")
}

void UniaxialMenegottoPintoMaterialLaw::CalculateMaterialResponsePK2(Parameters& rValues)
{
    KRATOS_TRY

    Flags& r_options = rValues.GetOptions();
    if (r_options.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        Vector& r_strain_vector = rValues.GetStrainVector();

        if (r_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
            Vector& r_stress_vector = rValues.GetStressVector();
            CalculateStressResponsePK2(rValues, r_strain_vector, r_stress_vector);
        }

        if (r_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
            CalculateConstitutiveMatrix(rValues, r_constitutive_matrix);
        }
    }
    else {
        KRATOS_ERROR << "Only ELEMENT_PROVIDED_STRAIN is implemented" << std::endl;
    }

    KRATOS_CATCH("")
}

void UniaxialMenegottoPintoMaterialLaw::FinalizeMaterialResponsePK2(Parameters& rValues)
{
    KRATOS_TRY
    mConvergedLoadingIndex  =                mLoadingIndex;
    mConvergedStrain0       =                     mStrain0;
    mConvergedStress0       =                     mStress0;
    mConvergedStrainR       =                     mStrainR;
    mConvergedStressR       =                     mStressR;
    mConvergedStrainPlastic =               mStrainPlastic;
    mConvergedStrainMax     =                   mStrainMax;
    mConvergedStrainMin     =                   mStrainMin;
    mConvergedStrain        = rValues.GetStrainVector()[0];
    mConvergedStress        = rValues.GetStressVector()[0];
    KRATOS_CATCH("")
}

void UniaxialMenegottoPintoMaterialLaw::CalculateStressResponsePK2(
    Parameters& rValues, Vector& rStrainVector, Vector& rStressVector)
{
    KRATOS_TRY
    double fy = rValues.GetMaterialProperties()[STEEL_YIELD_STRENGTH];
    double E  = rValues.GetMaterialProperties()[STEEL_YOUNGS_MODULUS];
    double b  = rValues.GetMaterialProperties()[STEEL_HARDENING_RATIO];
    double R0 = rValues.GetMaterialProperties()[STEEL_TRANSITION_VARIABLE];
    double a1 = rValues.GetMaterialProperties()[STEEL_A1_COEFFICIENT];
    double a2 = rValues.GetMaterialProperties()[STEEL_A2_COEFFICIENT];
    double E_inf = b * E;
    double epsy  = fy / E;

    mStrainMax      =      mConvergedStrainMax;
    mStrainMin      =      mConvergedStrainMin;
    mStrainPlastic  =  mConvergedStrainPlastic;
    mStrain0        =        mConvergedStrain0;
    mStress0        =        mConvergedStress0;
    mStrainR        =        mConvergedStrainR;
    mStressR        =        mConvergedStressR;
    mLoadingIndex   =   mConvergedLoadingIndex;

    double deps = rStrainVector[0] - mConvergedStrain;

    if (mLoadingIndex == 0 || mLoadingIndex == 3)
    {
        if (std::abs(deps) < std::numeric_limits<double>::epsilon())
        {
            return;
        }
        else
        {
            mStrainMax = epsy;
            mStrainMin = -epsy;
            if (deps < 0.0) {
                mLoadingIndex = 2;
                mStrain0 = mStrainMin;
                mStress0 = -fy;
                mStrainPlastic = mStrainMin;
            } else {
                mLoadingIndex = 1;
                mStrain0 = mStrainMax;
                mStress0 = fy;
                mStrainPlastic = mStrainMax;
            }
        }
    }

    if (mLoadingIndex == 2 && deps > 0.0)
    {
        // load reversal
        mLoadingIndex = 1;
        mStrainR = mConvergedStrain;
        mStressR = mConvergedStress;
        if (mConvergedStrain < mStrainMin) {
            mStrainMin = mConvergedStrain;
        }
        mStrain0 = (fy - E_inf*epsy - mStressR + E*mStrainR) / (E - E_inf);
        mStress0 = fy + E_inf*(mStrain0-epsy);
        mStrainPlastic = mStrainMax;
    }
    else if (mLoadingIndex == 1 && deps < 0.0)
    {
        // load reversal
        mLoadingIndex = 2;
        mStrainR = mConvergedStrain;
        mStressR = mConvergedStress;
        if (mConvergedStrain > mStrainMax) {
            mStrainMax = mConvergedStrain;
        }
        mStrain0 = (-fy + E_inf*epsy - mStressR + E*mStrainR) / (E - E_inf);
        mStress0 = -fy + E_inf*(mStrain0+epsy);
        mStrainPlastic = mStrainMin;
    }

    double xi = std::abs( (mStrainPlastic - mStrain0) / epsy );
    double R = R0 - a1*xi / (a2 + xi);
    double eps_star = (rStrainVector[0] - mStrainR) / (mStrain0 - mStrainR);
    double dummy1 = 1.0 + std::pow( std::abs(eps_star), R );
    double dummy2 = std::pow( dummy1, 1.0/R );
    double sg_star = b * eps_star + (1.0 - b) * eps_star / dummy2;
    rStressVector[0] = sg_star * (mStress0 - mStressR) + mStressR;
    mTangentModulus = b + (1.0 - b) / (dummy1 * dummy2);
    mTangentModulus *= (mStress0-mStressR) / (mStrain0-mStrainR);
    KRATOS_CATCH("")
}

void UniaxialMenegottoPintoMaterialLaw::CalculateConstitutiveMatrix(
    Parameters& rValues, Matrix& rConstitutiveMatrix)
{
    KRATOS_TRY
    rConstitutiveMatrix(0,0) = mTangentModulus;
    KRATOS_CATCH("")
}

std::string UniaxialMenegottoPintoMaterialLaw::Info() const {
    std::stringstream buffer;
    buffer << "UniaxialMenegottoPintoMaterialLaw";
    return buffer.str();
}

void UniaxialMenegottoPintoMaterialLaw::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "UniaxialMenegottoPintoMaterialLaw";
}

void UniaxialMenegottoPintoMaterialLaw::PrintData(std::ostream& rOStream) const
{
}

void UniaxialMenegottoPintoMaterialLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.save("mTangentModulus", mTangentModulus);
    rSerializer.save("mLoadingIndex", mLoadingIndex);
    rSerializer.save("mStrain0", mStrain0);
    rSerializer.save("mStress0", mStress0);
    rSerializer.save("mStrainR", mStrainR);
    rSerializer.save("mStressR", mStressR);
    rSerializer.save("mStrainPlastic", mStrainPlastic);
    rSerializer.save("mStrainMax", mStrainMax);
    rSerializer.save("mStrainMin", mStrainMin);
    rSerializer.save("mConvergedLoadingIndex", mConvergedLoadingIndex);
    rSerializer.save("mConvergedStrain0", mConvergedStrain0);
    rSerializer.save("mConvergedStress0", mConvergedStress0);
    rSerializer.save("mConvergedStrainR", mConvergedStrainR);
    rSerializer.save("mConvergedStressR", mConvergedStressR);
    rSerializer.save("mConvergedStrainPlastic", mConvergedStrainPlastic);
    rSerializer.save("mConvergedStrainMax", mConvergedStrainMax);
    rSerializer.save("mConvergedStrainMin", mConvergedStrainMin);
    rSerializer.save("mConvergedStrain", mConvergedStrain);
    rSerializer.save("mConvergedStress", mConvergedStress);
}
void UniaxialMenegottoPintoMaterialLaw::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.load("mTangentModulus", mTangentModulus);
    rSerializer.load("mLoadingIndex", mLoadingIndex);
    rSerializer.load("mStrain0", mStrain0);
    rSerializer.load("mStress0", mStress0);
    rSerializer.load("mStrainR", mStrainR);
    rSerializer.load("mStressR", mStressR);
    rSerializer.load("mStrainPlastic", mStrainPlastic);
    rSerializer.load("mStrainMax", mStrainMax);
    rSerializer.load("mStrainMin", mStrainMin);
    rSerializer.load("mConvergedLoadingIndex", mConvergedLoadingIndex);
    rSerializer.load("mConvergedStrain0", mConvergedStrain0);
    rSerializer.load("mConvergedStress0", mConvergedStress0);
    rSerializer.load("mConvergedStrainR", mConvergedStrainR);
    rSerializer.load("mConvergedStressR", mConvergedStressR);
    rSerializer.load("mConvergedStrainPlastic", mConvergedStrainPlastic);
    rSerializer.load("mConvergedStrainMax", mConvergedStrainMax);
    rSerializer.load("mConvergedStrainMin", mConvergedStrainMin);
    rSerializer.load("mConvergedStrain", mConvergedStrain);
    rSerializer.load("mConvergedStress", mConvergedStress);
}

} // namespace Kratos.
