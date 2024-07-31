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

ConstitutiveLaw::Pointer JointedCohesiveFrictionalConstitutiveLaw::Clone() const
{
    return Kratos::make_shared<JointedCohesiveFrictionalConstitutiveLaw>(*this);
}

/***********************************************************************************/
/***********************************************************************************/

void JointedCohesiveFrictionalConstitutiveLaw::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{


// todo



}



/***********************************************************************************/
/***********************************************************************************/

void JointedCohesiveFrictionalConstitutiveLaw::CalculateMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rValues
    )
{
    // todo
}

/***********************************************************************************/
/***********************************************************************************/

void JointedCohesiveFrictionalConstitutiveLaw::FinalizeMaterialResponseCauchy(
    ConstitutiveLaw::Parameters& rValues
    )
{



// todo


}

/***********************************************************************************/
/***********************************************************************************/

bool JointedCohesiveFrictionalConstitutiveLaw::Has(
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

bool JointedCohesiveFrictionalConstitutiveLaw::Has(
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

double& JointedCohesiveFrictionalConstitutiveLaw::GetValue(
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

Vector& JointedCohesiveFrictionalConstitutiveLaw::GetValue(
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

void JointedCohesiveFrictionalConstitutiveLaw::SetValue(
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

void JointedCohesiveFrictionalConstitutiveLaw::SetValue(
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

double& JointedCohesiveFrictionalConstitutiveLaw::CalculateValue(
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
    return rValue;
}


/***********************************************************************************/
/***********************************************************************************/

int JointedCohesiveFrictionalConstitutiveLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    // The auxiliary output
    int aux_out = 0;
    return aux_out;
}

/***********************************************************************************/
/***********************************************************************************/

void JointedCohesiveFrictionalConstitutiveLaw::CalculateTangentTensor(
    ConstitutiveLaw::Parameters& rValues
    )
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();

    const bool consider_perturbation_threshold = r_material_properties.Has(CONSIDER_PERTURBATION_THRESHOLD) ?
        r_material_properties[CONSIDER_PERTURBATION_THRESHOLD] : true;
    const TangentOperatorEstimation tangent_operator_estimation = r_material_properties.Has(TANGENT_OPERATOR_ESTIMATION) ?
        static_cast<TangentOperatorEstimation>(r_material_properties[TANGENT_OPERATOR_ESTIMATION]) : TangentOperatorEstimation::SecondOrderPerturbation;

    if (tangent_operator_estimation == TangentOperatorEstimation::Analytic) {

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

/***********************************************************************************/
/***********************************************************************************/


} // Namespace Kratos
