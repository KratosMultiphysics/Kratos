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
#include "custom_constitutive/uniaxial_elastic_material_law.hpp"
#include "structural_mechanics_application_variables.h"

namespace Kratos {


UniaxialElasticMaterialLaw::UniaxialElasticMaterialLaw()
    : ConstitutiveLaw() {}

UniaxialElasticMaterialLaw::UniaxialElasticMaterialLaw(
    const UniaxialElasticMaterialLaw& rOther )
    : ConstitutiveLaw(rOther)
{
}

ConstitutiveLaw::Pointer UniaxialElasticMaterialLaw::Clone() const
{
    return Kratos::make_shared<UniaxialElasticMaterialLaw>(*this);
}

UniaxialElasticMaterialLaw::~UniaxialElasticMaterialLaw()
{
}

void UniaxialElasticMaterialLaw::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mStrainSize = 1;
    rFeatures.mSpaceDimension = 3;
}

void UniaxialElasticMaterialLaw::CalculateMaterialResponsePK2(Parameters& rValues)
{
    KRATOS_TRY
    Flags& r_options = rValues.GetOptions();
    if (r_options.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        Vector& r_strain_vector = rValues.GetStrainVector();

        if (r_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
            Vector& r_stress_vector = rValues.GetStressVector();
            r_stress_vector = rValues.GetMaterialProperties()[YOUNG_MODULUS] * r_strain_vector;
        }

        if (r_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
            r_constitutive_matrix(0,0) = rValues.GetMaterialProperties()[YOUNG_MODULUS];
        }
    }
    else {
        KRATOS_ERROR << "Only element provided strain is implemented" << std::endl;
    }
    KRATOS_CATCH("")
}

void UniaxialElasticMaterialLaw::FinalizeMaterialResponsePK2(Parameters& rValues)
{
    KRATOS_TRY
    KRATOS_CATCH("")
}

// double& UniaxialElasticMaterialLaw::GetValue(const Variable<double>& rThisVariable, double& rValue)
// {
//     if (rThisVariable == PK2_STRESS_VECTOR) {}
// }


std::string UniaxialElasticMaterialLaw::Info() const {
    std::stringstream buffer;
    buffer << "UniaxialElasticMaterialLaw";
    return buffer.str();
}

void UniaxialElasticMaterialLaw::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "UniaxialElasticMaterialLaw";
}

void UniaxialElasticMaterialLaw::PrintData(std::ostream& rOStream) const
{
}

void UniaxialElasticMaterialLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
}
void UniaxialElasticMaterialLaw::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
}

} // namespace Kratos.
