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
    mInitialTangentModulus = 2.0 * mFc / mStrain0;
    // initialize material history variables
    mUnloadSlope             = mInitialTangentModulus;
    mTangentModulus          = mInitialTangentModulus;
    mConvergedUnloadSlope    = mInitialTangentModulus;
    mConvergedTangentModulus = mInitialTangentModulus;
    KRATOS_CATCH("")
}

void UniaxialFiberBeamColumnConcreteMaterialLaw::Confine()
{
    KRATOS_TRY
    // mFc *= GetProperties()[CONCRETE_CONFINEMENT_FACTOR];
    mStrainUltimate *= GetProperties()[CONCRETE_CONFINEMENT_FACTOR];
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

    mStrainMin      = mConvergedStrainMin;
    mStrainEnd      = mConvergedStrainEnd;
    mUnloadSlope    = mConvergedUnloadSlope;
    mStress         = mConvergedStress;
    mTangentModulus = mConvergedTangentModulus;

    double deps = mStrain - mConvergedStrain;
    if (std::abs(deps) < std::numeric_limits<double>::epsilon()){
        return;
    }

    if (mStrain >= 0.0){
        mTangentModulus = 0.0;
        mStress = 0.0;
        return;
    }
    /************************************************/
    // else if (mStrain >= mStrain0) {
    //     mTangentModulus = mFc / mStrain0;
    //     mStress = mTangentModulus * mStrain;
    // }
    // else if (mStrain >= mStrainUltimate) {
    //     // mTangentModulus = -0.8*mFc / (mStrainUltimate-mStrain0);
    //     mTangentModulus = 0.2*mFc / mStrainUltimate;
    //     mStress = mTangentModulus * mStrain;
    // }
    // else {
    //     mTangentModulus = 0.0;
    //     mStress = 0.2*mFc;
    // }
    // return;
    /************************************************/

    double sgr_tmp = mConvergedStress + mUnloadSlope*mStrain - mUnloadSlope*mConvergedStrain;

    // KRATOS_WATCH(sgr_tmp)

    // further into compression
    if (mStrain < mConvergedStrain) {
        Reload();
        if (sgr_tmp > mStress) {
            mStress = sgr_tmp;
            mTangentModulus = mUnloadSlope;
        }
    }

    // towards tension
    else if (sgr_tmp <= 0.0) {
        // if (mConvergedStrain == mStrainMin) reversal = True;
        mStress = sgr_tmp;
        mTangentModulus = mUnloadSlope;
    }

    // further into tension
    else {
        mStress = 0.0;
        mTangentModulus = 0.0;
    }

    KRATOS_CATCH("")
}

void UniaxialFiberBeamColumnConcreteMaterialLaw::Reload()
{
    KRATOS_TRY
    if (mStrain < mStrainMin) {
        mStrainMin = mStrain;
        Envelope();
        Unload();
    } else if (mStrain < mStrainEnd) {
        mTangentModulus = mUnloadSlope;
        mStress = mTangentModulus * (mStrain-mStrainEnd);
    } else {
        mStress = 0.0;
        mTangentModulus = 0.0;
    }
    KRATOS_CATCH("")
}

void UniaxialFiberBeamColumnConcreteMaterialLaw::Envelope()
{
    KRATOS_TRY
    if (mStrain > mStrain0) {
        double eta = mStrain / mStrain0;
        mStress = mFc * (2.0*eta - eta*eta);
        mTangentModulus = mInitialTangentModulus * (1.0 - eta);
    } else if (mStrain >= mStrainUltimate) {
        mTangentModulus = (mFc - 0.2*mFc) / (mStrain0-mStrainUltimate);
        mStress = mFc + mTangentModulus*(mStrain-mStrain0);
    } else {
        mStress = 0.2*mFc;
        mTangentModulus = 0.0;
    }
    KRATOS_CATCH("")
}

void UniaxialFiberBeamColumnConcreteMaterialLaw::Unload()
{
    KRATOS_TRY
    double eps_tmp = mStrainMin;
    // KRATOS_WATCH(eps_tmp)
    if (eps_tmp < mStrainUltimate) {
        eps_tmp = mStrainUltimate;
    }
    // KRATOS_WATCH(eps_tmp)
    double eta = eps_tmp / mStrain0;
    double ratio = 0.707 * (eta-2.0) + 0.834;
    if (eta < 2.0) {
        ratio = 0.145*eta*eta + 0.13*eta;
    }
    mStrainEnd = ratio * mStrain0;

    mUnloadSlope = mStress / (mStrainMin - mStrainEnd);

    double tmp1 = mStrainMin - mStrainEnd;
    double tmp2 = mStress / mInitialTangentModulus;
    if (tmp1 > -std::numeric_limits<double>::epsilon()) {
        mUnloadSlope = mInitialTangentModulus;
    } else if (tmp1 <= tmp2) {
        mStrainEnd = mStrainMin - tmp1;
        mUnloadSlope = mStress / tmp1;
    } else {
        mStrainEnd = mStrainMin - tmp2;
        mUnloadSlope = mInitialTangentModulus;
    }
    KRATOS_CATCH("")
}

void UniaxialFiberBeamColumnConcreteMaterialLaw::FinalizeMaterialResponse()
{
    KRATOS_TRY
    mConvergedStress         = mStress;
    mConvergedStrain         = mStrain;
    mConvergedStrainMin      = mStrainMin;
    mConvergedStrainEnd      = mStrainEnd;
    mConvergedUnloadSlope    = mUnloadSlope;
    mConvergedTangentModulus = mTangentModulus;
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
