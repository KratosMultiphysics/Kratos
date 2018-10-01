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
#include "custom_constitutive/viscous_generalized_kelvin.h"
#include "custom_constitutive/elastic_isotropic_3d.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
/******************************CONSTRUCTOR******************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
ViscousGeneralizedKelvin<TElasticBehaviourLaw>::ViscousGeneralizedKelvin()
    : BaseType()
{
}

/******************************COPY CONSTRUCTOR*************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
ViscousGeneralizedKelvin<TElasticBehaviourLaw>::ViscousGeneralizedKelvin(const ViscousGeneralizedKelvin& rOther)
    : BaseType(rOther),
      mPrevStressVector(rOther.mPrevStressVector),
      mPrevInelasticStrainVector(rOther.mPrevInelasticStrainVector),
      mNonConvPrevStressVector(rOther.mNonConvPrevStressVector),
      mNonConvPrevInelasticStrainVector(rOther.mNonConvPrevInelasticStrainVector)
{
}

/********************************CLONE**********************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
ConstitutiveLaw::Pointer ViscousGeneralizedKelvin<TElasticBehaviourLaw>::Clone() const
{
    return Kratos::make_shared<ViscousGeneralizedKelvin>(*this);
}

/*******************************DESTRUCTOR******************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
ViscousGeneralizedKelvin<TElasticBehaviourLaw>::~ViscousGeneralizedKelvin()
{
};

/***********************************************************************************/
/***********************************************************************************/
    
template <class TElasticBehaviourLaw>
void ViscousGeneralizedKelvin<TElasticBehaviourLaw>::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
//     BaseType::CalculateMaterialResponsePK1(rValues);
    // TODO: This requires to be properly defined to be able to set the StressMeasure
    this->ComputeViscoElasticity(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
void ViscousGeneralizedKelvin<TElasticBehaviourLaw>::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
//     BaseType::CalculateMaterialResponsePK2(rValues);
    // TODO: This requires to be properly defined to be able to set the StressMeasure
    this->ComputeViscoElasticity(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
void ViscousGeneralizedKelvin<TElasticBehaviourLaw>::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
//     BaseType::CalculateMaterialResponseKirchhoff(rValues);
    // TODO: This requires to be properly defined to be able to set the StressMeasure
    this->ComputeViscoElasticity(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
void ViscousGeneralizedKelvin<TElasticBehaviourLaw>::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // TODO: This requires to be properly defined to be able to set the StressMeasure
    this->ComputeViscoElasticity(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
void ViscousGeneralizedKelvin<TElasticBehaviourLaw>::FinalizeSolutionStep(
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

template <class TElasticBehaviourLaw>
void ViscousGeneralizedKelvin<TElasticBehaviourLaw>::ComputeViscoElasticity(ConstitutiveLaw::Parameters& rValues)
{
    // Integrate kelvin Generalized
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const ProcessInfo& r_process_info = rValues.GetProcessInfo();
    const double time_step = r_process_info[DELTA_TIME];
    const double delay_time = r_material_properties[DELAY_TIME];
    const Flags& r_flags = rValues.GetOptions();
    
    // We compute the strain in case not provided
    if( r_flags.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN )) {
        Vector& r_strain_vector = rValues.GetStrainVector();
        this->CalculateValue(rValues, STRAIN, r_strain_vector);
    }
    
    // We compute the stress
    if( r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS ) ) {
        // Elastic Matrix
        Matrix constitutive_matrix, inverse_constitutive_matrix;
        double detC = 0.0;
        this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, constitutive_matrix);
        MathUtils<double>::InvertMatrix(constitutive_matrix, inverse_constitutive_matrix, detC);

        // Ineslastic terms
        Vector inelastic_strain = this->GetPreviousInelasticStrainVector();
        const Vector& previous_stress = this->GetPreviousStressVector();
        const IndexType number_sub_increments = 10;
        const double dt = time_step / number_sub_increments;
        
        // Auxiliar terms
        Vector aux_stress_vector(previous_stress);
        Vector aux(VoigtSize);
        Vector elastic_strain(VoigtSize);
    
        // Apply increments
        const Vector& r_strain_vector = rValues.GetStrainVector();
        for (IndexType i = 0; i < number_sub_increments; ++i) {
            noalias(aux) = (std::exp(-dt / delay_time) * prod(inverse_constitutive_matrix, aux_stress_vector)) / delay_time;
            noalias(inelastic_strain) = std::exp(-dt / delay_time) * inelastic_strain + aux;
            noalias(elastic_strain) = r_strain_vector - inelastic_strain;
            noalias(aux_stress_vector) = prod(constitutive_matrix, elastic_strain);
        }
        
        Vector& integrated_stress_vector = rValues.GetStressVector(); // To be updated
        noalias(integrated_stress_vector) = aux_stress_vector;
        if (r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            rValues.SetConstitutiveMatrix(constitutive_matrix);
            this->SetNonConvPreviousStressVector(integrated_stress_vector);
            this->SetNonConvPreviousInelasticStrainVector(inelastic_strain);
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
void ViscousGeneralizedKelvin<TElasticBehaviourLaw>::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    BaseType::FinalizeMaterialResponsePK1(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
void ViscousGeneralizedKelvin<TElasticBehaviourLaw>::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    BaseType::FinalizeMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
void ViscousGeneralizedKelvin<TElasticBehaviourLaw>::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    BaseType::FinalizeMaterialResponseKirchhoff(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
void ViscousGeneralizedKelvin<TElasticBehaviourLaw>::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    BaseType::FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TElasticBehaviourLaw>
Vector& ViscousGeneralizedKelvin<TElasticBehaviourLaw>::CalculateValue(
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
Matrix& ViscousGeneralizedKelvin<TElasticBehaviourLaw>::CalculateValue(
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
int ViscousGeneralizedKelvin<TElasticBehaviourLaw>::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const int check_base = BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);

    KRATOS_ERROR_IF_NOT(VoigtSize == this->GetStrainSize()) << "You are combining not compatible constitutive laws" << std::endl;

    return check_base;
}

/***********************************************************************************/
/***********************************************************************************/

template class ViscousGeneralizedKelvin<ElasticIsotropic3D>;

} // namespace Kratos
