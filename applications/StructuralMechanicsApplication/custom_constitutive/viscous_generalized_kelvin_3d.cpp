// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo&  Lucia Barbu
//  Collaborator:    Vicente Mataix Ferrandiz
//

// System includes

// Project includes
#include "custom_constitutive/viscous_generalized_kelvin_3d.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
void ViscousGeneralizedKelvin3D::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void ViscousGeneralizedKelvin3D::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void ViscousGeneralizedKelvin3D::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void ViscousGeneralizedKelvin3D::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // Integrate kelvin Generalized
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    Vector& integrated_stress_vector = rValues.GetStressVector(); // To be updated
    const Vector& strain_vector = rValues.GetStrainVector();
    Matrix& tangent_tensor = rValues.GetConstitutiveMatrix(); // todo modify after integration
    const ProcessInfo& r_process_info = rValues.GetProcessInfo();
    const double time_step = r_process_info[DELTA_TIME];
    const double delay_time = r_material_properties[DELAY_TIME];
    const Flags& constitutive_law_options = rValues.GetOptions();

    // Elastic Matrix
    Matrix C, inverse_C;
    double detC = 0.0;
    this->CalculateElasticMatrix(C, r_material_properties);
    MathUtils<double>::InvertMatrix(C, inverse_C, detC);

    Vector inelastic_strain = this->GetPreviousInelasticStrainVector();
    const Vector& previous_stress = this->GetPreviousStressVector();
    const IndexType number_sub_increments = 10;
    const double dt = time_step / number_sub_increments;

    Vector aux_stress_vector;
    aux_stress_vector = previous_stress;
    Vector aux = ZeroVector(6);

    Vector elastic_strain;
    for (IndexType i = 0; i < number_sub_increments; i++) {
        aux = (std::exp(-dt / delay_time) * prod(inverse_C, aux_stress_vector)) / delay_time;
        inelastic_strain = std::exp(-dt / delay_time) * inelastic_strain + aux;
        elastic_strain = strain_vector - inelastic_strain;
        noalias(aux_stress_vector) = prod(C, elastic_strain);
    }
    noalias(integrated_stress_vector) = aux_stress_vector;
    if (constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        noalias(tangent_tensor) = C;
        this->SetNonConvPreviousStressVector(integrated_stress_vector);
        this->SetNonConvPreviousInelasticStrainVector(inelastic_strain);
    }
} // End CalculateMaterialResponseCauchy

/***********************************************************************************/
/***********************************************************************************/

void ViscousGeneralizedKelvin3D::FinalizeSolutionStep(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Update the required vectors
    this->SetPreviousInelasticStrainVector(this->GetNonConvPreviousInelasticStrainVector());
    this->SetPreviousStressVector(this->GetNonConvPreviousStressVector());
}

/***********************************************************************************/
/***********************************************************************************/

void ViscousGeneralizedKelvin3D::CalculateElasticMatrix(
    Matrix& rElasticityTensor,
    const Properties& rMaterialProperties)
{
    const double E = rMaterialProperties[YOUNG_MODULUS];
    const double poisson_ratio = rMaterialProperties[POISSON_RATIO];
    const double lambda = E * poisson_ratio / ((1. + poisson_ratio) * (1.0 - 2.0 * poisson_ratio));
    const double mu = E / (2.0 + 2.0 * poisson_ratio);

    if (rElasticityTensor.size1() != 6 || rElasticityTensor.size2() != 6)
        rElasticityTensor.resize(6, 6, false);
    rElasticityTensor.clear();

    rElasticityTensor(0, 0) = lambda + 2.0 * mu;
    rElasticityTensor(0, 1) = lambda;
    rElasticityTensor(0, 2) = lambda;
    rElasticityTensor(1, 0) = lambda;
    rElasticityTensor(1, 1) = lambda + 2.0 * mu;
    rElasticityTensor(1, 2) = lambda;
    rElasticityTensor(2, 0) = lambda;
    rElasticityTensor(2, 1) = lambda;
    rElasticityTensor(2, 2) = lambda + 2.0 * mu;
    rElasticityTensor(3, 3) = mu;
    rElasticityTensor(4, 4) = mu;
    rElasticityTensor(5, 5) = mu;
}

/***********************************************************************************/
/***********************************************************************************/

void ViscousGeneralizedKelvin3D::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    BaseType::FinalizeMaterialResponsePK1(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void ViscousGeneralizedKelvin3D::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    BaseType::FinalizeMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void ViscousGeneralizedKelvin3D::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    BaseType::FinalizeMaterialResponseKirchhoff(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void ViscousGeneralizedKelvin3D::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    BaseType::FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& ViscousGeneralizedKelvin3D::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    if (rThisVariable == INTEGRATED_STRESS_TENSOR) {
        rValue = MathUtils<double>::StressVectorToTensor(this->GetPreviousStressVector());
    }

    return rValue;
}

} // namespace Kratos
