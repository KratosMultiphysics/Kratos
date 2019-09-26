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
#include "includes/checks.h"
#include "utilities/math_utils.h"
#include "structural_mechanics_application_variables.h"
#include "custom_advanced_constitutive/viscous_generalized_maxwell.h"

namespace Kratos
{
/******************************CONSTRUCTOR******************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
ViscousGeneralizedMaxwell<TElasticBehaviourLaw>::ViscousGeneralizedMaxwell()
    : BaseType()
{
}

/******************************COPY CONSTRUCTOR*************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
ViscousGeneralizedMaxwell<TElasticBehaviourLaw>::ViscousGeneralizedMaxwell(const ViscousGeneralizedMaxwell& rOther)
    : BaseType(rOther),
      mPrevStressVector(rOther.mPrevStressVector),
      mPrevStrainVector(rOther.mPrevStrainVector)
{
}

/********************************CLONE**********************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
ConstitutiveLaw::Pointer ViscousGeneralizedMaxwell<TElasticBehaviourLaw>::Clone() const
{
    return Kratos::make_shared<ViscousGeneralizedMaxwell>(*this);
}

/*******************************DESTRUCTOR******************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
ViscousGeneralizedMaxwell<TElasticBehaviourLaw>::~ViscousGeneralizedMaxwell()
{
};

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
void ViscousGeneralizedMaxwell<TElasticBehaviourLaw>::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    this->ComputeViscoElasticity(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
void ViscousGeneralizedMaxwell<TElasticBehaviourLaw>::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    this->ComputeViscoElasticity(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
void ViscousGeneralizedMaxwell<TElasticBehaviourLaw>::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    this->ComputeViscoElasticity(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
void ViscousGeneralizedMaxwell<TElasticBehaviourLaw>::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // TODO: This requires to be properly defined to be able to set the StressMeasure
    this->ComputeViscoElasticity(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
void ViscousGeneralizedMaxwell<TElasticBehaviourLaw>::FinalizeSolutionStep(
    const Properties &rMaterialProperties,
    const GeometryType &rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo &rCurrentProcessInfo)
{
    // Deprecated
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
void ViscousGeneralizedMaxwell<TElasticBehaviourLaw>::ComputeViscoElasticity(ConstitutiveLaw::Parameters& rValues)
{
    const Properties& material_props = rValues.GetMaterialProperties();
    Vector& integrated_stress_vector = rValues.GetStressVector(); // To be updated
    const ProcessInfo& process_info = rValues.GetProcessInfo();
    const double time_step = process_info[DELTA_TIME];
    const Flags& r_flags = rValues.GetOptions();

    // We compute the strain in case not provided
    if (r_flags.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        Vector& r_strain_vector = rValues.GetStrainVector();
        this->CalculateValue(rValues, STRAIN, r_strain_vector);
    }

    // We compute the stress
    if (r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS)) {
        const double viscous_parameter = material_props[VISCOUS_PARAMETER]; //  C1/Cinf
        const double delay_time = material_props[DELAY_TIME];

        // Elastic Matrix
        Matrix constitutive_matrix;
        this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, constitutive_matrix);

        const Vector& r_strain_vector = rValues.GetStrainVector();
        const Vector& r_previous_strain = this->GetPreviousStrainVector();
        const Vector& r_previous_stress = this->GetPreviousStressVector();
        const Vector& r_strain_increment = r_strain_vector - r_previous_strain;

        const double coef = viscous_parameter * time_step / ((1.0 + viscous_parameter) * 2.0 * delay_time);
        const Vector& r_auxiliar_strain = -(r_strain_vector - r_strain_increment) * std::exp(-time_step / delay_time) * (1.0 + coef) + r_strain_vector * (1.0 - coef);

        noalias(integrated_stress_vector) = r_previous_stress * std::exp(-time_step / delay_time) + prod(constitutive_matrix, r_auxiliar_strain);

        if (r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            rValues.SetConstitutiveMatrix(constitutive_matrix);
        }
    } else {
        // Elastic Matrix
        if( r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
            Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
            this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
void ViscousGeneralizedMaxwell<TElasticBehaviourLaw>::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
void ViscousGeneralizedMaxwell<TElasticBehaviourLaw>::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
void ViscousGeneralizedMaxwell<TElasticBehaviourLaw>::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
void ViscousGeneralizedMaxwell<TElasticBehaviourLaw>::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    const Properties& material_props = rValues.GetMaterialProperties();
    Vector& integrated_stress_vector = rValues.GetStressVector(); // To be updated
    const ProcessInfo& process_info = rValues.GetProcessInfo();
    const double time_step = process_info[DELTA_TIME];
    const Flags& r_flags = rValues.GetOptions();

    // We compute the strain in case not provided
    if (r_flags.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        Vector& r_strain_vector = rValues.GetStrainVector();
        this->CalculateValue(rValues, STRAIN, r_strain_vector);
    }


    const double viscous_parameter = material_props[VISCOUS_PARAMETER]; //  C1/Cinf
    const double delay_time = material_props[DELAY_TIME];

    // Elastic Matrix
    Matrix constitutive_matrix;
    this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, constitutive_matrix);

    const Vector& r_strain_vector = rValues.GetStrainVector();
    const Vector& r_previous_strain = this->GetPreviousStrainVector();
    const Vector& r_previous_stress = this->GetPreviousStressVector();
    const Vector& r_strain_increment = r_strain_vector - r_previous_strain;

    const double coef = viscous_parameter * time_step / ((1.0 + viscous_parameter) * 2.0 * delay_time);
    const Vector& r_auxiliar_strain = -(r_strain_vector - r_strain_increment) * std::exp(-time_step / delay_time) * (1.0 + coef) + r_strain_vector * (1.0 - coef);

    noalias(integrated_stress_vector) = r_previous_stress * std::exp(-time_step / delay_time) + prod(constitutive_matrix, r_auxiliar_strain);

    mPrevStressVector = integrated_stress_vector;
    mPrevStrainVector = r_strain_vector;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
Vector& ViscousGeneralizedMaxwell<TElasticBehaviourLaw>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    return BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
Matrix& ViscousGeneralizedMaxwell<TElasticBehaviourLaw>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    if (rThisVariable == INTEGRATED_STRESS_TENSOR) {
        rValue = MathUtils<double>::StressVectorToTensor(this->GetPreviousStressVector());
    } else {
        rValue = BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
int ViscousGeneralizedMaxwell<TElasticBehaviourLaw>::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const int check_base = BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);

    KRATOS_CHECK_VARIABLE_KEY(VISCOUS_PARAMETER);
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(VISCOUS_PARAMETER)) << "VISCOUS_PARAMETER is not a defined value" << std::endl;

    KRATOS_ERROR_IF_NOT(VoigtSize == this->GetStrainSize()) << "You are combining not compatible constitutive laws" << std::endl;

    return check_base;
}

/***********************************************************************************/
/***********************************************************************************/

template class ViscousGeneralizedMaxwell<ElasticIsotropic3D>;

} // namespace Kratos
