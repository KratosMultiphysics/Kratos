// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:        BSD License
//                  license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//
//
// System includes
#include <iostream>
#include <set>

// External includes

// Project includes
#include "includes/checks.h"
#include "jointed_cohesive_frictional_constitutive_law.h"
#include "constitutive_laws_application_variables.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"

namespace Kratos
{

/********************************CLONE**********************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
ConstitutiveLaw::Pointer JointedCohesiveFrictionalConstitutiveLaw<TYieldSurfaceType>::Clone() const
{
    return Kratos::make_shared<JointedCohesiveFrictionalConstitutiveLaw>(*this);
}


/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void JointedCohesiveFrictionalConstitutiveLaw<TYieldSurfaceType>::CalculateMaterialResponseCauchy(
    ConstitutiveLaw::Parameters& rValues
    )
{







}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void JointedCohesiveFrictionalConstitutiveLaw<TYieldSurfaceType>::FinalizeMaterialResponseCauchy(
    ConstitutiveLaw::Parameters& rValues
    )
{






}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
bool JointedCohesiveFrictionalConstitutiveLaw<TYieldSurfaceType>::Has(
    const Variable<double>& rThisVariable
    )
{
    bool has = false;
    // if (rThisVariable == PLASTIC_DISSIPATION) {
    //     has = true;
    // } else if (rThisVariable == THRESHOLD) {
    //     has = true;
    // } else if (rThisVariable == DAMAGE) {
    //     has = true;
    // } else if (rThisVariable == DISSIPATION) {
    //     has = true;
    // }
    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
bool JointedCohesiveFrictionalConstitutiveLaw<TYieldSurfaceType>::Has(
    const Variable<Vector>& rThisVariable
    )
{
    bool has = false;
    // if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
    //     has = true;
    // }
    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
double& JointedCohesiveFrictionalConstitutiveLaw<TYieldSurfaceType>::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    rValue = 0.0;
    // if (rThisVariable == PLASTIC_DISSIPATION) {
    //     rValue = mPlasticDissipation;
    // } else if (rThisVariable == THRESHOLD) {
    //     rValue = mThreshold;
    // }  else if (rThisVariable == DAMAGE) {
    //     rValue = mDamageDissipation;
    // } else if (rThisVariable == DISSIPATION) {
    //     rValue = mPlasticDissipation + mDamageDissipation;
    // }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
Vector& JointedCohesiveFrictionalConstitutiveLaw<TYieldSurfaceType>::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    rValue.resize(VoigtSize, false);
    rValue.clear();
    // if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
    //     noalias(rValue) = mPlasticStrain;
    // }
    return rValue;
}


/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void JointedCohesiveFrictionalConstitutiveLaw<TYieldSurfaceType>::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // if (rThisVariable == PLASTIC_DISSIPATION) {
    //     mPlasticDissipation = rValue;
    // } else if (rThisVariable == THRESHOLD) {
    //     mThreshold = rValue;
    // }  else if (rThisVariable == DAMAGE) {
    //     mDamageDissipation = rValue;
    // }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void JointedCohesiveFrictionalConstitutiveLaw<TYieldSurfaceType>::SetValue(
    const Variable<Vector>& rThisVariable,
    const Vector& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
    //     noalias(mPlasticStrain) = rValue;
    // }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
double& JointedCohesiveFrictionalConstitutiveLaw<TYieldSurfaceType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    // if (rThisVariable == UNIAXIAL_STRESS) {
    //     // Get Values to compute the constitutive law:
    //     Flags& r_flags = rParameterValues.GetOptions();

    //     // Previous flags saved
    //     const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
    //     const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS );

    //     r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false );
    //     r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, true );

    //     // Calculate the stress vector
    //     CalculateMaterialResponseCauchy(rParameterValues);
    //     const Vector& r_stress_vector = rParameterValues.GetStressVector();
    //     const Vector& r_strain_vector = rParameterValues.GetStrainVector();

    //     BoundedVectorType aux_stress_vector = r_stress_vector;
    //     TYieldSurfaceType::CalculateEquivalentStress( aux_stress_vector, r_strain_vector, rValue, rParameterValues);

    //     // Previous flags restored
    //     r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
    //     r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );
    // } else {
    //     BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
    // }
    // return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void JointedCohesiveFrictionalConstitutiveLaw<TYieldSurfaceType>::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{






}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void JointedCohesiveFrictionalConstitutiveLaw<TYieldSurfaceType>::CalculateMaterialResponsePK1(
    ConstitutiveLaw::Parameters& rValues
    )
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void  JointedCohesiveFrictionalConstitutiveLaw<TYieldSurfaceType>::CalculateMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rValues
    )
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void JointedCohesiveFrictionalConstitutiveLaw<TYieldSurfaceType>::CalculateMaterialResponseKirchhoff(
    ConstitutiveLaw::Parameters& rValues
    )
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void JointedCohesiveFrictionalConstitutiveLaw<TYieldSurfaceType>::InitializeMaterialResponsePK1(
    ConstitutiveLaw::Parameters& rValues
    )
{
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void JointedCohesiveFrictionalConstitutiveLaw<TYieldSurfaceType>::InitializeMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rValues
    )
{
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void JointedCohesiveFrictionalConstitutiveLaw<TYieldSurfaceType>::InitializeMaterialResponseKirchhoff(
    ConstitutiveLaw::Parameters& rValues
    )
{
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void JointedCohesiveFrictionalConstitutiveLaw<TYieldSurfaceType>::InitializeMaterialResponseCauchy(
    ConstitutiveLaw::Parameters& rValues
    )
{

}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void JointedCohesiveFrictionalConstitutiveLaw<TYieldSurfaceType>::FinalizeMaterialResponsePK1(
    ConstitutiveLaw::Parameters& rValues
    )
{
    FinalizeMaterialResponseCauchy(rValues);
}

template<class TYieldSurfaceType>
/***********************************************************************************/
/***********************************************************************************/

void JointedCohesiveFrictionalConstitutiveLaw<TYieldSurfaceType>::FinalizeMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rValues
    )
{
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void JointedCohesiveFrictionalConstitutiveLaw<TYieldSurfaceType>::FinalizeMaterialResponseKirchhoff(
    ConstitutiveLaw::Parameters& rValues
    )
{
    FinalizeMaterialResponseCauchy(rValues);
}


/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
int JointedCohesiveFrictionalConstitutiveLaw<TYieldSurfaceType>::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    // The auxiliary output
    int aux_out = 0;
    // KRATOS_ERROR_IF(!rMaterialProperties.Has(FRACTURE_ENERGY))           << "FRACTURE_ENERGY not provided in the material properties" << std::endl;
    // KRATOS_ERROR_IF(!rMaterialProperties.Has(HARDENING_CURVE))           << "HARDENING_CURVE not provided in the material properties" << std::endl;
    // KRATOS_ERROR_IF(!rMaterialProperties.Has(PLASTIC_DAMAGE_PROPORTION)) << "PLASTIC_DAMAGE_PROPORTION not provided in the material properties" << std::endl;

    // // Checking curves
    // const int curve_type = rMaterialProperties[HARDENING_CURVE];
    // if (static_cast<typename GenericConstitutiveLawIntegratorPlasticity<TYieldSurfaceType>::HardeningCurveType>(curve_type) == GenericConstitutiveLawIntegratorPlasticity<TYieldSurfaceType>::HardeningCurveType::CurveDefinedByPoints) {
    //     KRATOS_ERROR_IF(!rMaterialProperties.Has(EQUIVALENT_STRESS_VECTOR_PLASTICITY_POINT_CURVE))  << "EQUIVALENT_STRESS_VECTOR_PLASTICITY_POINT_CURVE not provided in the material properties" << std::endl;
    //     KRATOS_ERROR_IF(!rMaterialProperties.Has(TOTAL_STRAIN_VECTOR_PLASTICITY_POINT_CURVE))       << "TOTAL_STRAIN_VECTOR_PLASTICITY_POINT_CURVE not provided in the material properties" << std::endl;
    // }
    return aux_out;
}


/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void JointedCohesiveFrictionalConstitutiveLaw<TYieldSurfaceType>::CalculateTangentTensor(
    ConstitutiveLaw::Parameters& rValues
    )
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();

    const bool consider_perturbation_threshold = r_material_properties.Has(CONSIDER_PERTURBATION_THRESHOLD) ?
        r_material_properties[CONSIDER_PERTURBATION_THRESHOLD] : true;
    const TangentOperatorEstimation tangent_operator_estimation = r_material_properties.Has(TANGENT_OPERATOR_ESTIMATION) ?
        static_cast<TangentOperatorEstimation>(r_material_properties[TANGENT_OPERATOR_ESTIMATION]) : TangentOperatorEstimation::SecondOrderPerturbation;

    if (tangent_operator_estimation == TangentOperatorEstimation::Analytic) {
        // CalculateAnalyticalTangentTensor(rValues, rPlasticDamageParameters);
        noalias(rValues.GetConstitutiveMatrix()) = rPlasticDamageParameters.TangentTensor;
    } else if (tangent_operator_estimation == TangentOperatorEstimation::FirstOrderPerturbation) {
        // Calculates the Tangent Constitutive Tensor by perturbation (first order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 1);
    } else if (tangent_operator_estimation == TangentOperatorEstimation::SecondOrderPerturbation) {
        // Calculates the Tangent Constitutive Tensor by perturbation (second order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 2);
    } else if (tangent_operator_estimation == TangentOperatorEstimation::SecondOrderPerturbationV2) {
        // Calculates the Tangent Constitutive Tensor by perturbation (second order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 4);
    } else if (tangent_operator_estimation == TangentOperatorEstimation::Secant) {
        
    }
}

} // Namespace Kratos
