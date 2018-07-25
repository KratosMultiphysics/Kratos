// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo & Lucia Barbu
//  Collaborator:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "utilities/math_utils.h"
#include "structural_mechanics_application_variables.h"
#include "custom_constitutive/viscous_generalized_maxwell_3d.h"

namespace Kratos
{
void ViscousGeneralizedMaxwell3D::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void ViscousGeneralizedMaxwell3D::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void ViscousGeneralizedMaxwell3D::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void ViscousGeneralizedMaxwell3D::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // Integrate Stress Damage
    const Properties &material_props = rValues.GetMaterialProperties();
    Vector& integrated_stress_vector = rValues.GetStressVector(); // To be updated
    const Vector& strain_vector = rValues.GetStrainVector();
    Matrix& tangent_tensor = rValues.GetConstitutiveMatrix(); // todo modify after integration
    const ProcessInfo &process_info = rValues.GetProcessInfo();
    const double time_step = process_info[DELTA_TIME];
    const Flags &ConstitutiveLawOptions = rValues.GetOptions();

    const double viscous_parameter = material_props[VISCOUS_PARAMETER]; //  C1/Cinf
    const double delay_time = material_props[DELAY_TIME];

    // Elastic Matrix
    Matrix C;
    this->CalculateElasticMatrix(C, material_props);

    const Vector& previous_strain = this->GetPreviousStrainVector();
    const Vector& previous_stress = this->GetPreviousStressVector();
    const Vector& strain_increment = strain_vector - previous_strain;

    const double coef = viscous_parameter * time_step / ((1.0 + viscous_parameter) * 2.0 * delay_time);
    const Vector& auxiliar_strain = -(strain_vector - strain_increment) * std::exp(-time_step / delay_time) * (1.0 + coef) + strain_vector * (1.0 - coef);

    noalias(integrated_stress_vector) = previous_stress * std::exp(-time_step / delay_time) + prod(C, auxiliar_strain);

    if (ConstitutiveLawOptions.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        noalias(tangent_tensor) = C;
        this->SetNonConvPreviousStressVector(integrated_stress_vector);
        this->SetNonConvPreviousStrainVector(strain_vector);
    }
} // End CalculateMaterialResponseCauchy

/***********************************************************************************/
/***********************************************************************************/

void ViscousGeneralizedMaxwell3D::FinalizeSolutionStep(
    const Properties &rMaterialProperties,
    const GeometryType &rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo &rCurrentProcessInfo)
{
    // Update the required vectors
    this->SetPreviousStrainVector(this->GetNonConvPreviousStrainVector());
    this->SetPreviousStressVector(this->GetNonConvPreviousStressVector());
}

/***********************************************************************************/
/***********************************************************************************/

void ViscousGeneralizedMaxwell3D::CalculateElasticMatrix(
    Matrix& rElasticityTensor,
    const Properties &rMaterialProperties)
{
    const double E = rMaterialProperties[YOUNG_MODULUS];
    const double poisson_ratio = rMaterialProperties[POISSON_RATIO];
    const double lambda =
        E * poisson_ratio / ((1. + poisson_ratio) * (1.0 - 2.0 * poisson_ratio));
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

void ViscousGeneralizedMaxwell3D::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    BaseType::FinalizeMaterialResponsePK1(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void ViscousGeneralizedMaxwell3D::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    BaseType::FinalizeMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void ViscousGeneralizedMaxwell3D::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    BaseType::FinalizeMaterialResponseKirchhoff(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void ViscousGeneralizedMaxwell3D::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    BaseType::FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& ViscousGeneralizedMaxwell3D::CalculateValue(
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
