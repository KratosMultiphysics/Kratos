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
    KRATOS_TRY;

    Flags &r_constitutive_law_options = rValues.GetOptions();
    const auto& r_strain_vector = rValues.GetStrainVector();
    const auto strain_increment = r_strain_vector - mOldStrainVector;

    if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        // We retrieve material properties
        const double E  = rValues.GetMaterialProperties()[YOUNG_MODULUS];
        const double nu = rValues.GetMaterialProperties()[POISSON_RATIO];
        const Vector &r_joint_cl_props = rValues.GetMaterialProperties()[CURVE_FITTING_PARAMETERS];
        const double fc  = r_joint_cl_props[0];
        const double ft  = r_joint_cl_props[1];
        const double Kn  = r_joint_cl_props[2];
        const double Ks  = r_joint_cl_props[3];
        const double muy = r_joint_cl_props[4];
        const double delta_0 = r_joint_cl_props[5];
        const double H = r_joint_cl_props[6];
        const double m = r_joint_cl_props[7];
        const double muy0   = r_joint_cl_props[8];
        const double alpha0 = r_joint_cl_props[9];
        const double beta   = r_joint_cl_props[10];
        const double gamma  = r_joint_cl_props[11];
        Matrix Kc(Dimension, Dimension);
        Kc.clear();
        Kc(0, 0) = Kn;
        Kc(1, 1) = Ks;

        this->CalculateElasticMatrix(rValues.GetConstitutiveMatrix(), rValues);
        Matrix &a0 = rValues.GetConstitutiveMatrix();

        // Compute the trial stress
        Vector stress_increment_trial(VoigtSize);
        noalias(stress_increment_trial) = prod(a0, strain_increment);

        // Rotation and normal matrices
        Matrix R(Dimension, Dimension);
        Matrix n(VoigtSize, Dimension);

        // Estimate joint orientation with respect to some trial stress
        CalculateJointOrientation(R, n, mOldStressVector + stress_increment_trial, mDamage, muy0, muy, ft, m, fc);

        const Matrix C = prod(trans(n), Matrix(prod(a0, n))) / H + prod(trans(R), Matrix(prod(Kc, R)));
        Matrix inv_C(Dimension, Dimension);
        double det_C;
        MathUtils<double>::InvertMatrix(C, inv_C, det_C);
        const Matrix aux = prod(inv_C, trans(n));
        const Vector delta_u_trial = prod(aux, stress_increment_trial);
        const Vector delta_uc_trial = prod(R, delta_u_trial);
        const Vector delta_stress_trial = stress_increment_trial - prod(a0, Vector(prod(n, delta_u_trial))) / H;
        const Vector tc_trial = mLocalTraction + prod(R, Vector(prod(trans(n), delta_stress_trial)));

        const double yield = YieldSurfaceValue(tc_trial[1], mDamage, muy0, muy, tc_trial[0], ft, m, fc);

        if (yield < tolerance) { // Elastic condition
            noalias(rValues.GetStressVector()) = stress_trial;
        } else {

        }








    }


    if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {



    }

    KRATOS_CATCH("");
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
