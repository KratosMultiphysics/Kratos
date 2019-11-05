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
    : UniaxialFiberBeamColumnMaterialLaw() {}

UniaxialFiberBeamColumnConcreteMaterialLaw::UniaxialFiberBeamColumnConcreteMaterialLaw(
    PropertiesType::Pointer pProperties )
    : UniaxialFiberBeamColumnMaterialLaw(pProperties)
{
    mTangentModulus = GetProperties()[CONCRETE_YIELD_STRENGTH];
}

void UniaxialFiberBeamColumnConcreteMaterialLaw::Initialize()
{
    KRATOS_TRY
    KRATOS_CATCH("")
}

void UniaxialFiberBeamColumnConcreteMaterialLaw::CalculateMaterialResponse()
{
    KRATOS_TRY
    mStress = mTangentModulus * mStrain;
    KRATOS_CATCH("")
}

void UniaxialFiberBeamColumnConcreteMaterialLaw::FinalizeMaterialResponse()
{
    KRATOS_TRY
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
    rOStream << "    ConfinementFactor : " << GetProperties()[CONCRETE_CONFINEMENT_FACTOR] << std::endl;
    rOStream << "    SofteningSlope    : " << GetProperties()[CONCRETE_SOFTENING_SLOPE] << std::endl;
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
