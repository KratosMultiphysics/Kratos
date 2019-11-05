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
#include <iostream>

// External includes

// Project includes
#include "custom_constitutive/uniaxial_fiber_beam_column_steel_material_law.hpp"

namespace Kratos {


UniaxialFiberBeamColumnSteelMaterialLaw::UniaxialFiberBeamColumnSteelMaterialLaw()
    : UniaxialFiberBeamColumnMaterialLaw() {}

UniaxialFiberBeamColumnSteelMaterialLaw::UniaxialFiberBeamColumnSteelMaterialLaw(
    PropertiesType::Pointer pProperties )
    : UniaxialFiberBeamColumnMaterialLaw(pProperties)
{
    mTangentModulus = GetProperties()[STEEL_YOUNGS_MODULUS];
    mStrainMax = GetProperties()[STEEL_YIELD_STRENGTH] / GetProperties()[STEEL_YOUNGS_MODULUS];
    mStrainMin = -1.0 * mStrainMax;
    mConvergedStrainMax = mStrainMax;
    mConvergedStrainMin = mStrainMin;
}

void UniaxialFiberBeamColumnSteelMaterialLaw::Initialize()
{
    KRATOS_TRY
    KRATOS_CATCH("")
}

void UniaxialFiberBeamColumnSteelMaterialLaw::CalculateMaterialResponse()
{
    KRATOS_TRY

    double deps = mStrain - mConvergedStrain;

    mStrainMax     = mConvergedStrainMax;
    mStrainMin     = mConvergedStrainMin;
    mStrainPlastic = mConvergedStrainPlastic;
    mStrain0       = mConvergedStrain0;
    mStress0       = mConvergedStress0;
    mStrainR       = mConvergedStrainR;
    mStressR       = mConvergedStressR;
    mLoadingIndex  = mConvergedLoadingIndex;

    double fy = GetProperties()[STEEL_YIELD_STRENGTH];
    double E  = GetProperties()[STEEL_YOUNGS_MODULUS];
    double b  = GetProperties()[STEEL_HARDENING_RATIO];
    double R0 = GetProperties()[STEEL_TRANSITION_VARIABLE];
    double a1 = GetProperties()[STEEL_A1_COEFFICIENT];
    double a2 = GetProperties()[STEEL_A2_COEFFICIENT];
    double E_inf = b * E;
    double epsy  = fy / E;

    if ((mLoadingIndex == 0) || (mLoadingIndex == 3))
    {
        if (std::abs(deps) < std::numeric_limits<double>::epsilon())
        {
            mTangentModulus = E;
            mStress = 0.0;
            mLoadingIndex = 3;
            return;
        }
        else
        {
            mStrainMax = epsy;
            mStrainMin = -epsy;
            if (deps < 0) {
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

    if ((mLoadingIndex == 2) && deps > 0)
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
    else if ((mLoadingIndex == 1) && (deps < 0))
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
        mStrainPlastic = mStrainMax;
    }

    double plastic_excursion = std::abs( (mStrainPlastic - mStrain0) / epsy );
    double R = R0 - a1 * ( plastic_excursion / (a2 + plastic_excursion) );
    double eps_star = (mStrain - mStrainR) / (mStrain0 - mStrainR);
    double dummy1 = 1.0 + std::pow( std::abs(eps_star), R );
    double dummy2 = std::pow( dummy1, 1.0/R );
    double sg_star = b * eps_star + (1.0 - b) * eps_star / dummy2;
    mStress = sg_star * (mStress0 - mStressR) + mStressR;
    mTangentModulus = b + (1.0 - b) / (dummy1 * dummy2);
    mTangentModulus *= (mStress0-mStressR) / (mStrain0-mStrainR);

    KRATOS_CATCH("")
}

void UniaxialFiberBeamColumnSteelMaterialLaw::FinalizeMaterialResponse()
{
    KRATOS_TRY
    mConvergedLoadingIndex  = mLoadingIndex;
    mConvergedStrain0       = mStrain0;
    mConvergedStress0       = mStress0;
    mConvergedStrainR       = mStrainR;
    mConvergedStressR       = mStressR;
    mConvergedStrainPlastic = mStrainPlastic;
    mConvergedStrainMax     = mStrainMax;
    mConvergedStrainMin     = mStrainMin;
    mConvergedStrain        = mStrain;
    mConvergedStress        = mStress;
    KRATOS_CATCH("")
}


std::string UniaxialFiberBeamColumnSteelMaterialLaw::Info() const {
    std::stringstream buffer;
    buffer << "UniaxialFiberBeamColumnSteelMaterialLaw";
    return buffer.str();
}

void UniaxialFiberBeamColumnSteelMaterialLaw::PrintInfo(std::ostream& rOStream) const {
    rOStream << "UniaxialFiberBeamColumnSteelMaterialLaw";
}

void UniaxialFiberBeamColumnSteelMaterialLaw::PrintData(std::ostream& rOStream) const {
    if (pGetProperties() == nullptr) {
        rOStream << "This material parameters are not given." << std::endl;
    } else {
    rOStream << "    YoungsModulus   : " << GetProperties()[STEEL_YOUNGS_MODULUS] << std::endl;
    rOStream << "    HardeningRation : " << GetProperties()[STEEL_HARDENING_RATIO] << std::endl;
    rOStream << "    Transition      : " << GetProperties()[STEEL_TRANSITION_VARIABLE] << std::endl;
    rOStream << "    YieldStrength   : " << GetProperties()[STEEL_YIELD_STRENGTH] << std::endl;
    rOStream << "    A1              : " << GetProperties()[STEEL_A1_COEFFICIENT] << std::endl;
    rOStream << "    A2              : " << GetProperties()[STEEL_A2_COEFFICIENT] << std::endl;
    }
}

void UniaxialFiberBeamColumnSteelMaterialLaw::save(Serializer& rSerializer) const
{
    // rSerializer.save("mId", mId);
}
void UniaxialFiberBeamColumnSteelMaterialLaw::load(Serializer& rSerializer)
{
    // rSerializer.load("mId", mId);
}

} // namespace Kratos.
