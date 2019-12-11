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
#include "includes/properties.h"
#include "custom_constitutive/uniaxial_fiber_elastic_material_law.hpp"
#include "structural_mechanics_application_variables.h"

namespace Kratos {


UniaxialFiberElasticMaterialLaw::UniaxialFiberElasticMaterialLaw()
    : ConstitutiveLaw() {}

UniaxialFiberElasticMaterialLaw::UniaxialFiberElasticMaterialLaw(
    const UniaxialFiberElasticMaterialLaw& rOther )
    : ConstitutiveLaw(rOther)
{
}

ConstitutiveLaw::Pointer UniaxialFiberElasticMaterialLaw::Clone() const
{
    return Kratos::make_shared<UniaxialFiberElasticMaterialLaw>(*this);
}

UniaxialFiberElasticMaterialLaw::~UniaxialFiberElasticMaterialLaw()
{
}

void UniaxialFiberElasticMaterialLaw::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mStrainSize = 1;
    rFeatures.mSpaceDimension = 3;
}

void UniaxialFiberElasticMaterialLaw::CalculateMaterialResponsePK2(Parameters& rValues)
{
    KRATOS_TRY
    Vector strain = ZeroVector(1);
    Vector stress = ZeroVector(1);
    strain = rValues.GetStrainVector();
    stress = rValues.GetMaterialProperties()[YOUNG_MODULUS] * strain;
    rValues.SetStressVector(stress);
    KRATOS_CATCH("")
}

void UniaxialFiberElasticMaterialLaw::FinalizeMaterialResponsePK2(Parameters& rValues)
{
    KRATOS_TRY
    KRATOS_CATCH("")
}

// double& UniaxialFiberElasticMaterialLaw::GetValue(const Variable<double>& rThisVariable, double& rValue)
// {
//     if (rThisVariable == PK2_STRESS_VECTOR) {}
// }


std::string UniaxialFiberElasticMaterialLaw::Info() const {
    std::stringstream buffer;
    buffer << "UniaxialFiberElasticMaterialLaw";
    return buffer.str();
}

void UniaxialFiberElasticMaterialLaw::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "UniaxialFiberElasticMaterialLaw";
}

void UniaxialFiberElasticMaterialLaw::PrintData(std::ostream& rOStream) const
{
}

void UniaxialFiberElasticMaterialLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
}
void UniaxialFiberElasticMaterialLaw::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
}

} // namespace Kratos.
