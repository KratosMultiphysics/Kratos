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
#include "custom_constitutive/uniaxial_fiber_beam_column_concrete_material_law2.hpp"

namespace Kratos {

UniaxialFiberBeamColumnConcreteMaterialLaw2::UniaxialFiberBeamColumnConcreteMaterialLaw2()
    : UniaxialFiberBeamColumnMaterialLaw()
{
}

UniaxialFiberBeamColumnConcreteMaterialLaw2::UniaxialFiberBeamColumnConcreteMaterialLaw2(
    PropertiesType::Pointer pProperties )
    : UniaxialFiberBeamColumnMaterialLaw(pProperties)
{
    KRATOS_TRY
    mFc = -std::abs(GetProperties()[CONCRETE_YIELD_STRENGTH]);
    mStrainU = -std::abs(GetProperties()[CONCRETE_CRUSHING_STRAIN]);
    mStrain0 = -std::abs(GetProperties()[CONCRETE_YIELD_STRAIN]);
    mTangentModulus = 2 * mFc / mStrain0;
    KRATOS_CATCH("")
}

void UniaxialFiberBeamColumnConcreteMaterialLaw2::Confine()
{
    KRATOS_TRY
    mStrainU *= GetProperties()[CONCRETE_CONFINEMENT_FACTOR];
    KRATOS_CATCH("")
}

void UniaxialFiberBeamColumnConcreteMaterialLaw2::Initialize()
{
    KRATOS_TRY
    KRATOS_CATCH("")
}

void UniaxialFiberBeamColumnConcreteMaterialLaw2::CalculateMaterialResponse()
{
    KRATOS_TRY

    // positive strain
    if (mStrain >= 0.0) {
        mStress = 0.0;
        mTangentModulus = 0.0;
        return;
    }

    double deps = mStrain - mConvStrain;

    // set current loading index
    if (std::abs(deps) < std::numeric_limits<double>::epsilon()) {
        mLoadingIndex = 3;
        // return;
    } else {
        if (deps < 0.0) {
            mLoadingIndex = 2;
        } else {
            mLoadingIndex = 1;
        }
    }

    // check reversal
    bool reversal = false;
    if (std::abs(mStrain) > std::numeric_limits<double>::epsilon()) {
        if (mConvLoadingIndex == 2 || mConvLoadingIndex == 3) {
            if (mLoadingIndex == 1) {
                reversal = true;
            }
        }
        if (mConvLoadingIndex == 1 || mConvLoadingIndex == 3) {
            if (mLoadingIndex == 2) {
                reversal = true;
            }
        }
    }

    // reverse if true and on loading path
    if (reversal)
    {
        if (mLoadingIndex != 2)
        {
            mStrainR = mConvStrain;
            mStressR = mConvStress;

            double crit = mStrainR / mStrain0;
            if (crit >= 2){
                mStrainP = mStrain0 * (0.707 * (crit - 2) + 0.834);
            } else if (crit < 2){
                mStrainP = mStrain0 * (0.145*crit*crit + 0.13*crit);
            }
        }
    }

    /*** calculate stress and tangent modulus ***/
    // loading path
    if (mStrain <= mStrainR) {

        // before yielding part
        if (mStrain >= mStrain0) {
            // double ratio = mStrain / mStrain0;
            // mStress = mK * mFc * (2 * ratio - ratio*ratio);
            // mTangentModulus = mK * mFc * (2/mStrain0 - 2*ratio*ratio);
            mTangentModulus = mFc / mStrain0;
            mStress = mTangentModulus * mStrain;
        }

        // // softening part
        // else {
        //     mStress = mK * mFc * (1 + mZ*(mStrain - mStrain0));
        //     // crushing part
        //     if (mStress <= 0.2 * mK * mFc) {
        //         mStress = 0.2 * mK * mFc;
        //         mTangentModulus = 0.0;
        //     } else {
        //         mTangentModulus = mK * mFc * mZ;
        //     }
        // }
        else if (mStrain >= mStrainU) {
            mTangentModulus = 0.0;
            mStress = mFc;
        }
        else {
            mTangentModulus = 0.0;
            mStress = 0.2*mFc;
        }
    }

    // unloading path
    else {
        // crushed concrete
        if (mStrain >= mStrainP) {
            mStress = 0.0;
            mTangentModulus = 0.0;
            return;
        }
        mStress = (mStressR*mStrain - mStrainP*mStressR) / (mStrainR - mStrainP);
        mTangentModulus = mStressR / (mStressR - mStrainP);
    }

    // mStress *= -1.0;
    // mTangentModulus *= -1.0;
    KRATOS_CATCH("")
}

void UniaxialFiberBeamColumnConcreteMaterialLaw2::FinalizeMaterialResponse()
{
    KRATOS_TRY
    mConvStress = mStress;
    mConvStrain = mStrain;
    mConvStressR = mStressR;
    mConvStrainR = mStrainR;
    mConvStrainP = mStrainP;
    mConvLoadingIndex = mLoadingIndex;
    KRATOS_CATCH("")
}


std::string UniaxialFiberBeamColumnConcreteMaterialLaw2::Info() const {
    std::stringstream buffer;
    buffer << "UniaxialFiberBeamColumnConcreteMaterialLaw2";
    return buffer.str();
}

void UniaxialFiberBeamColumnConcreteMaterialLaw2::PrintInfo(std::ostream& rOStream) const {
    rOStream << "UniaxialFiberBeamColumnConcreteMaterialLaw2";
}

void UniaxialFiberBeamColumnConcreteMaterialLaw2::PrintData(std::ostream& rOStream) const {
    if (pGetProperties() == nullptr) {
        rOStream << "This material parameters are not given." << std::endl;
    } else {
    rOStream << "    YieldStrength     : " << GetProperties()[CONCRETE_YIELD_STRENGTH] << std::endl;
    rOStream << "    YieldStrain       : " << GetProperties()[CONCRETE_YIELD_STRAIN] << std::endl;
    }
}

void UniaxialFiberBeamColumnConcreteMaterialLaw2::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, UniaxialFiberBeamColumnMaterialLaw);
}
void UniaxialFiberBeamColumnConcreteMaterialLaw2::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, UniaxialFiberBeamColumnMaterialLaw);
}

} // namespace Kratos.
