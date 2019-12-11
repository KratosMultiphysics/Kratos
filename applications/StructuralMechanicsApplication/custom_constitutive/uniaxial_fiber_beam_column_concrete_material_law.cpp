// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//  license:     structural_mechanics_application/license.txt
//
//  Main authors: Mahmoud Zidan
//

// System includes
#include<iostream>

// External includes

// Project includes
#include "custom_constitutive/uniaxial_fiber_beam_column_concrete_material_law.hpp"

namespace Kratos {

UniaxialFiberBeamColumnConcreteMaterialLaw::UniaxialFiberBeamColumnConcreteMaterialLaw()
    : UniaxialFiberBeamColumnMaterialLaw()
{
}

UniaxialFiberBeamColumnConcreteMaterialLaw::UniaxialFiberBeamColumnConcreteMaterialLaw(
    PropertiesType::Pointer pProperties )
    : UniaxialFiberBeamColumnMaterialLaw(pProperties)
{
    KRATOS_TRY
    // make sure the given parameters are negative
    mStrain0 = -std::abs(GetProperties()[CONCRETE_YIELD_STRAIN]);
    mFc = -std::abs(GetProperties()[CONCRETE_YIELD_STRENGTH]);
    mStrainUltimate = -std::abs(GetProperties()[CONCRETE_CRUSHING_STRAIN]);
    double E0 = 2.0 * mFc / mStrain0;
    // initialize material history variables
    mUnloadSlope    = E0;
    mTangentModulus = E0;
    KRATOS_CATCH("")
}


void UniaxialFiberBeamColumnConcreteMaterialLaw::Initialize()
{
    KRATOS_TRY
    KRATOS_CATCH("")
}

void UniaxialFiberBeamColumnConcreteMaterialLaw::CalculateMaterialResponse()
{
    KRATOS_TRY

    mStrainR = mConvergedStrainR;
    mStrainP = mConvergedStrainP;
    mStress  = mConvergedStress;

    double deps = mStrain - mConvergedStrain;

    // quick return if the strain has not changed
    if (std::abs(deps) < std::numeric_limits<double>::epsilon()){
        return;
    }

    // quick return if the strain is positive
    if (mStrain > 0.0){
        mTangentModulus = 0.0;
        mStress = 0.0;
        return;
    }

    double sgr_tmp = mConvergedStress + mUnloadSlope*(mStrain - mConvergedStrain);

    // further into compression
    if (mStrain < mConvergedStrain) {
        Reload();
        if (sgr_tmp > mStress) {
            mStress = sgr_tmp;
            mTangentModulus = mUnloadSlope;
        }
    }
    // towards tension unloading path
    else if (mStrain < mStrainP) {
        mStress = sgr_tmp;
        mTangentModulus = mUnloadSlope;
    }
    // further into tension (open crack)
    else {
        mStress = 0.0;
        mTangentModulus = 0.0;
    }

    KRATOS_CATCH("")
}

void UniaxialFiberBeamColumnConcreteMaterialLaw::Reload()
{
    KRATOS_TRY
    // loading
    if (mStrain < mStrainR) {
        mStrainR = mStrain;
        Envelope();
        Unload();
    // towards compression on the reloading path
    } else if (mStrain < mStrainP) {
        mTangentModulus = mUnloadSlope;
        mStress = mTangentModulus * (mStrain-mStrainP);
    // towards compression but open crack
    } else {
        mStress = 0.0;
        mTangentModulus = 0.0;
    }
    KRATOS_CATCH("")
}

void UniaxialFiberBeamColumnConcreteMaterialLaw::Envelope()
{
    KRATOS_TRY
    // pre yielding
    if (mStrain > mStrain0) {
        double eta = mStrain / mStrain0;
        mStress = mFc * (2.0*eta - eta*eta);
        double E0 = 2.0 * mFc / mStrain0;
        mTangentModulus = E0 * (1.0 - eta);
    // softening
    } else if (mStrain >= mStrainUltimate) {
        mTangentModulus = 0.8*mFc / (mStrain0-mStrainUltimate);
        mStress = mFc + mTangentModulus*(mStrain-mStrain0);
    // past ultimate strain
    } else {
        mStress = 0.2*mFc;
        mTangentModulus = 0.0;
    }
    KRATOS_CATCH("")
}

void UniaxialFiberBeamColumnConcreteMaterialLaw::Unload()
{
    KRATOS_TRY
    double eps_tmp = mStrainR;
    if (eps_tmp < mStrainUltimate) {
        eps_tmp = mStrainUltimate;
    }
    double eta = eps_tmp / mStrain0;
    double ratio = 0.707 * (eta-2.0) + 0.834;
    if (eta < 2.0) {
        ratio = 0.145*eta*eta + 0.13*eta;
    }
    mStrainP = ratio * mStrain0;

    double E0 = 2.0 * mFc / mStrain0;
    double tmp1 = mStrainR - mStrainP;
    double tmp2 = mStress / E0;
    if (tmp1 > -std::numeric_limits<double>::epsilon()) {
        mUnloadSlope = E0;
    } else if (tmp1 <= tmp2) {
        mStrainP = mStrainR - tmp1;
        mUnloadSlope = mStress / tmp1;
    } else {
        mStrainP = mStrainR - tmp2;
        mUnloadSlope = E0;
    }
    KRATOS_CATCH("")
}

void UniaxialFiberBeamColumnConcreteMaterialLaw::FinalizeMaterialResponse()
{
    KRATOS_TRY
    mConvergedStress  = mStress;
    mConvergedStrain  = mStrain;
    mConvergedStrainR = mStrainR;
    mConvergedStrainP = mStrainP;
    KRATOS_CATCH("")
}


std::string UniaxialFiberBeamColumnConcreteMaterialLaw::Info() const {
    std::stringstream buffer;
    buffer << "UniaxialFiberBeamColumnConcreteMaterialLaw";
    return buffer.str();
}

void UniaxialFiberBeamColumnConcreteMaterialLaw::PrintInfo(std::ostream& rOStream) const {
    rOStream << "UniaxialFiberBeamColumnConcreteMaterialLaw";
}

void UniaxialFiberBeamColumnConcreteMaterialLaw::PrintData(std::ostream& rOStream) const {
    if (pGetProperties() == nullptr) {
        rOStream << "This material parameters are not given." << std::endl;
    } else {
    rOStream << "    YieldStrength     : " << GetProperties()[CONCRETE_YIELD_STRENGTH] << std::endl;
    rOStream << "    YieldStrain       : " << GetProperties()[CONCRETE_YIELD_STRAIN] << std::endl;
    rOStream << "    CrushingStrain    : " << GetProperties()[CONCRETE_CRUSHING_STRAIN] << std::endl;
    }
}

void UniaxialFiberBeamColumnConcreteMaterialLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, UniaxialFiberBeamColumnMaterialLaw);
}
void UniaxialFiberBeamColumnConcreteMaterialLaw::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, UniaxialFiberBeamColumnMaterialLaw);
}

} // namespace Kratos.
