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
#include "custom_constitutive/uniaxial_kent_park.h"
#include "seismic_application_variables.h"

namespace Kratos {


ConstitutiveLaw::Pointer UniaxialKentParkMaterialLaw::Clone() const
{
    return Kratos::make_shared<UniaxialKentParkMaterialLaw>(*this);
}

void UniaxialKentParkMaterialLaw::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mStrainSize = 1;
    rFeatures.mSpaceDimension = 3;
}

void UniaxialKentParkMaterialLaw::InitializeMaterialResponsePK2(Parameters& rValues)
{
    KRATOS_TRY
    mTangentModulus = 2.0 * rValues.GetMaterialProperties()[CONCRETE_YIELD_STRENGTH] / rValues.GetMaterialProperties()[CONCRETE_YIELD_STRAIN];
    KRATOS_CATCH("")
}

void UniaxialKentParkMaterialLaw::CalculateMaterialResponsePK2(Parameters& rValues)
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
            CalculateConstitutiveMatrixPK2(rValues, r_constitutive_matrix);
        }
    }
    else {
        KRATOS_ERROR << "Only ELEMENT_PROVIDED_STRAIN is implemented" << std::endl;
    }

    KRATOS_CATCH("")
}

void UniaxialKentParkMaterialLaw::FinalizeMaterialResponsePK2(Parameters& rValues)
{
    KRATOS_TRY
    mConvergedStress  = rValues.GetStressVector()[0];
    mConvergedStrain  = rValues.GetStrainVector()[0];
    mConvergedStrainR = mStrainR;
    mConvergedStrainP = mStrainP;
    KRATOS_CATCH("")
}

void UniaxialKentParkMaterialLaw::CalculateStressResponsePK2(
    Parameters& rValues, Vector& rStrainVector, Vector& rStressVector)
{
    KRATOS_TRY
    mStrainR = mConvergedStrainR;
    mStrainP = mConvergedStrainP;
    rStressVector[0] = mConvergedStress;
    double strain = rStrainVector[0];

    double deps = strain - mConvergedStrain;

    // quick return if the strain has not changed
    if (std::abs(deps) < std::numeric_limits<double>::epsilon()){
        return;
    }

    // quick return if the strain is positive
    if (strain > 0.0){
        mTangentModulus = 0.0;
        rStressVector[0] = 0.0;
        return;
    }

    double sgr_tmp = mConvergedStress + mUnloadSlope*(strain - mConvergedStrain);

    // further into compression
    if (strain < mConvergedStrain) {
        Reload(rValues, rStrainVector, rStressVector);
        if (sgr_tmp > rStressVector[0]) {
            rStressVector[0] = sgr_tmp;
            mTangentModulus = mUnloadSlope;
        }
    }
    // towards tension unloading path
    else if (strain < mStrainP) {
        rStressVector[0] = sgr_tmp;
        mTangentModulus = mUnloadSlope;
    }
    // further into tension (open crack)
    else {
        rStressVector[0] = 0.0;
        mTangentModulus = 0.0;
    }
    KRATOS_CATCH("")
}

void UniaxialKentParkMaterialLaw::Reload(
    Parameters& rValues, Vector& rStrainVector, Vector& rStressVector)
{
    KRATOS_TRY
    // loading
    if (rStrainVector[0] < mStrainR) {
        mStrainR = rStrainVector[0];
        Envelope(rValues, rStrainVector, rStressVector);
        Unload(rValues, rStrainVector, rStressVector);
    // towards compression on the reloading path
    } else if (rStrainVector[0] < mStrainP) {
        mTangentModulus = mUnloadSlope;
        rStressVector[0] = mTangentModulus * (rStrainVector[0]-mStrainP);
    // towards compression but open crack
    } else {
        rStressVector[0] = 0.0;
        mTangentModulus = 0.0;
    }
    KRATOS_CATCH("")
}

void UniaxialKentParkMaterialLaw::Envelope(
    Parameters& rValues, Vector& rStrainVector, Vector& rStressVector)
{
    KRATOS_TRY
    double strain0  = -1.0 * rValues.GetMaterialProperties()[CONCRETE_YIELD_STRAIN];
    double fc       = -1.0 * rValues.GetMaterialProperties()[CONCRETE_YIELD_STRENGTH];
    double strain_u = -1.0 * rValues.GetMaterialProperties()[CONCRETE_CRUSHING_STRAIN];
    // pre yielding
    if (rStrainVector[0] > strain0) {
        double eta = rStrainVector[0] / strain0;
        rStressVector[0] = fc * (2.0*eta - eta*eta);
        double E0 = 2.0 * fc / strain0;
        mTangentModulus = E0 * (1.0 - eta);
    // softening
    } else if (rStrainVector[0] >= strain_u) {
        mTangentModulus = 0.8*fc / (strain0-strain_u);
        rStressVector[0] = fc + mTangentModulus*(rStrainVector[0]-strain0);
    // past ultimate strain
    } else {
        rStressVector[0] = 0.2*fc;
        mTangentModulus = 0.0;
    }
    KRATOS_CATCH("")
}

void UniaxialKentParkMaterialLaw::Unload(
    Parameters& rValues, Vector& rStrainVector, Vector& rStressVector)
{
    KRATOS_TRY
    double strain0  = -1.0 * rValues.GetMaterialProperties()[CONCRETE_YIELD_STRAIN];
    double fc       = -1.0 * rValues.GetMaterialProperties()[CONCRETE_YIELD_STRENGTH];
    double strain_u = -1.0 * rValues.GetMaterialProperties()[CONCRETE_CRUSHING_STRAIN];
    double eps_tmp = mStrainR;
    if (eps_tmp < strain_u) {
        eps_tmp = strain_u;
    }
    double eta = eps_tmp / strain0;
    double ratio = 0.707 * (eta-2.0) + 0.834;
    if (eta < 2.0) {
        ratio = 0.145*eta*eta + 0.13*eta;
    }
    mStrainP = ratio * strain0;

    double E0 = 2.0 * fc / strain0;
    double tmp1 = mStrainR - mStrainP;
    double tmp2 = rStressVector[0] / E0;
    if (tmp1 > -std::numeric_limits<double>::epsilon()) {
        mUnloadSlope = E0;
    } else if (tmp1 <= tmp2) {
        mStrainP = mStrainR - tmp1;
        mUnloadSlope = rStressVector[0] / tmp1;
    } else {
        mStrainP = mStrainR - tmp2;
        mUnloadSlope = E0;
    }
    KRATOS_CATCH("")
}

void UniaxialKentParkMaterialLaw::CalculateConstitutiveMatrixPK2(
    Parameters& rValues, Matrix& rConstitutiveMatrix)
{
    KRATOS_TRY
    rConstitutiveMatrix(0,0) = mTangentModulus;
    KRATOS_CATCH("")
}

std::string UniaxialKentParkMaterialLaw::Info() const {
    std::stringstream buffer;
    buffer << "UniaxialKentParkMaterialLaw";
    return buffer.str();
}

void UniaxialKentParkMaterialLaw::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "UniaxialKentParkMaterialLaw";
}

void UniaxialKentParkMaterialLaw::PrintData(std::ostream& rOStream) const
{
}

void UniaxialKentParkMaterialLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.save("mStrainR", mStrainR);
    rSerializer.save("mStrainP", mStrainP);
    rSerializer.save("mUnloadSlope", mUnloadSlope);
    rSerializer.save("mTangentModulus", mTangentModulus);
    rSerializer.save("mConvergedStress", mConvergedStress);
    rSerializer.save("mConvergedStrain", mConvergedStrain);
    rSerializer.save("mConvergedStrainR", mConvergedStrainR);
    rSerializer.save("mConvergedStrainP", mConvergedStrainP);
}
void UniaxialKentParkMaterialLaw::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.load("mStrainR", mStrainR);
    rSerializer.load("mStrainP", mStrainP);
    rSerializer.load("mUnloadSlope", mUnloadSlope);
    rSerializer.load("mTangentModulus", mTangentModulus);
    rSerializer.load("mConvergedStress", mConvergedStress);
    rSerializer.load("mConvergedStrain", mConvergedStrain);
    rSerializer.load("mConvergedStrainR", mConvergedStrainR);
    rSerializer.load("mConvergedStrainP", mConvergedStrainP);
}

} // namespace Kratos.
