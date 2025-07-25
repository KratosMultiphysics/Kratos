// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//  Collaborator:
//

// System includes

// Project includes
#include "utilities/math_utils.h"
#include "constitutive_laws_application_variables.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"
#include "generic_small_strain_kinematic_plasticity.h"
#include "custom_constitutive/auxiliary_files/cl_integrators/generic_cl_integrator_kinematic_plasticity.h"

// Yield surfaces
#include "custom_constitutive/auxiliary_files/yield_surfaces/generic_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/von_mises_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/plane_stress_von_mises_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/modified_mohr_coulomb_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/mohr_coulomb_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/rankine_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/simo_ju_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/drucker_prager_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/tresca_yield_surface.h"

// Plastic potentials
#include "custom_constitutive/auxiliary_files/plastic_potentials/generic_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/von_mises_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/plane_stress_von_mises_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/tresca_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/modified_mohr_coulomb_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/mohr_coulomb_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/drucker_prager_plastic_potential.h"

namespace Kratos
{
template <class TConstLawIntegratorType>
void GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::CalculateMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rValues
    )
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::CalculateMaterialResponseKirchhoff(
    ConstitutiveLaw::Parameters& rValues
    )
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::CalculateMaterialResponsePK1(
    ConstitutiveLaw::Parameters& rValues
    )
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::CalculateMaterialResponseCauchy(
    ConstitutiveLaw::Parameters& rValues
    )
{
    // Auxiliary values
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();

    // We get the constitutive tensor
    Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
    this->CalculateElasticMatrix(r_constitutive_matrix, rValues);

    // We check the current step and NL iteration
    const ProcessInfo& r_current_process_info = rValues.GetProcessInfo();
    const bool first_computation = (r_current_process_info[NL_ITERATION_NUMBER] == 1 && r_current_process_info[STEP] == 1) ? true : false;

    //NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    if (first_computation) { // First computation always pure elastic for elemements not providing the strain
        if (r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN ) ) {
            BaseType::CalculateCauchyGreenStrain( rValues, r_strain_vector);
        }
        if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS) ||
            r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            Vector& r_stress_vector = rValues.GetStressVector();
            noalias(r_stress_vector) = prod(r_constitutive_matrix, r_strain_vector);
        }
    } else { // We check for plasticity
        // Integrate Stress plasticity
        Vector& r_integrated_stress_vector = rValues.GetStressVector();
        const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());

        //NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
        if (r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
            this->CalculateValue(rValues, STRAIN, r_strain_vector);
        }

        // We compute the stress
        if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {

            // We get some variables
            double threshold = this->GetThreshold();
            double plastic_dissipation = this->GetPlasticDissipation();
            Vector plastic_strain      = this->GetPlasticStrain();
            Vector back_stress_vector  = this->GetBackStressVector();
            const Vector previous_stress_vector = this->GetPreviousStressVector();

            array_1d<double, VoigtSize> predictive_stress_vector, kin_hard_stress_vector;
            if (r_constitutive_law_options.Is(ConstitutiveLaw::U_P_LAW)) {
                predictive_stress_vector = rValues.GetStressVector();
            } else {
                // S0 = r_constitutive_matrix:(E-Ep)
                predictive_stress_vector = prod(r_constitutive_matrix, r_strain_vector - plastic_strain);
            }

            // Initialize Plastic Parameters
            double uniaxial_stress = 0.0, plastic_denominator = 0.0;
            array_1d<double, VoigtSize> f_flux; // DF/DS
            array_1d<double, VoigtSize> g_flux; // DG/DS
            array_1d<double, VoigtSize> plastic_strain_increment;
            f_flux.clear();
            g_flux.clear();
            plastic_strain_increment.clear();

            // Kinematic back stress substracted
            noalias(kin_hard_stress_vector) = predictive_stress_vector - back_stress_vector;

            const double F = TConstLawIntegratorType::CalculatePlasticParameters(
                kin_hard_stress_vector, r_strain_vector, uniaxial_stress,
                threshold, plastic_denominator, f_flux, g_flux,
                plastic_dissipation, plastic_strain_increment,
                r_constitutive_matrix, rValues, characteristic_length,
                plastic_strain, back_stress_vector);

            if (F <= std::abs(1.0e-4 * threshold)) { // Elastic case
                noalias(r_integrated_stress_vector) = predictive_stress_vector;
            } else { // Plastic case
                // While loop backward euler
                /* Inside "IntegrateStressVector" the predictive_stress_vector is updated to verify the yield criterion */
                TConstLawIntegratorType::IntegrateStressVector(
                    predictive_stress_vector, r_strain_vector, uniaxial_stress,
                    threshold, plastic_denominator, f_flux, g_flux,
                    plastic_dissipation, plastic_strain_increment,
                    r_constitutive_matrix, plastic_strain, rValues,
                    characteristic_length, back_stress_vector,
                    previous_stress_vector);

                noalias(r_integrated_stress_vector) = predictive_stress_vector;

                if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                    this->CalculateTangentTensor(rValues, plastic_strain); // this modifies the ConstitutiveMatrix
                }
            }
        }
    }
} // End CalculateMaterialResponseCauchy

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::CalculateTangentTensor(
    ConstitutiveLaw::Parameters& rValues,
    const Vector& rPlasticStrain
    )
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();

    const bool consider_perturbation_threshold = r_material_properties.Has(CONSIDER_PERTURBATION_THRESHOLD) ? r_material_properties[CONSIDER_PERTURBATION_THRESHOLD] : true;
    const TangentOperatorEstimation tangent_operator_estimation = r_material_properties.Has(TANGENT_OPERATOR_ESTIMATION) ? static_cast<TangentOperatorEstimation>(r_material_properties[TANGENT_OPERATOR_ESTIMATION]) : TangentOperatorEstimation::SecondOrderPerturbation;

    if (tangent_operator_estimation == TangentOperatorEstimation::Analytic) {
        // Already stored in rValues.GetConstitutiveMatrix()...
    } else if (tangent_operator_estimation == TangentOperatorEstimation::FirstOrderPerturbation) {
        // Calculates the Tangent Constitutive Tensor by perturbation (first order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 1);
    } else if (tangent_operator_estimation == TangentOperatorEstimation::SecondOrderPerturbation) {
        // Calculates the Tangent Constitutive Tensor by perturbation (second order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 2);
    } else if (tangent_operator_estimation == TangentOperatorEstimation::Secant) {
        const Vector num = prod(rValues.GetConstitutiveMatrix(), rPlasticStrain);
        const double denom = inner_prod(rValues.GetStrainVector(), num);
        noalias(rValues.GetConstitutiveMatrix()) -= outer_prod(num, num) / denom;
    } else if (tangent_operator_estimation == TangentOperatorEstimation::SecondOrderPerturbationV2) {
        // Calculates the Tangent Constitutive Tensor by perturbation (second order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 4);
    } else if (tangent_operator_estimation == TangentOperatorEstimation::InitialStiffness) {
        this->CalculateElasticMatrix(rValues.GetConstitutiveMatrix(), rValues);
    } else if (tangent_operator_estimation == TangentOperatorEstimation::OrthogonalSecant) {
        TangentOperatorCalculatorUtility::CalculateOrthogonalSecantTensor(rValues);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    // We construct the CL parameters
    ProcessInfo dummy_process_info;
    ConstitutiveLaw::Parameters aux_param(rElementGeometry, rMaterialProperties, dummy_process_info);

    // We call the integrator
    double initial_threshold;
    TConstLawIntegratorType::GetInitialUniaxialThreshold(aux_param, initial_threshold);
    this->SetThreshold(initial_threshold);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::FinalizeMaterialResponsePK1(
    ConstitutiveLaw::Parameters& rValues
    )
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::FinalizeMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rValues
    )
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::FinalizeMaterialResponseKirchhoff(
    ConstitutiveLaw::Parameters& rValues
    )
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::FinalizeMaterialResponseCauchy(
    ConstitutiveLaw::Parameters& rValues
    )
{
    const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();
    Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
    this->CalculateElasticMatrix(r_constitutive_matrix, rValues);

    //NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    if ( r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        this->CalculateValue(rValues, STRAIN, r_strain_vector);
    }

    // We get some variables
    double threshold = this->GetThreshold();
    double plastic_dissipation = this->GetPlasticDissipation();
    Vector plastic_strain      = this->GetPlasticStrain();
    Vector back_stress_vector  = this->GetBackStressVector();
    const Vector previous_stress_vector = this->GetPreviousStressVector();

    array_1d<double, VoigtSize> predictive_stress_vector, kin_hard_stress_vector;
    if (r_constitutive_law_options.Is(ConstitutiveLaw::U_P_LAW)) {
        predictive_stress_vector = rValues.GetStressVector();
    } else {
        // S0 = r_constitutive_matrix:(E-Ep)
        predictive_stress_vector = prod(r_constitutive_matrix, r_strain_vector - plastic_strain);
    }

    // Initialize Plastic Parameters
    double uniaxial_stress = 0.0, plastic_denominator = 0.0;
    array_1d<double, VoigtSize> f_flux; // DF/DS
    array_1d<double, VoigtSize> g_flux; // DG/DS
    array_1d<double, VoigtSize> plastic_strain_increment;
    f_flux.clear();
    g_flux.clear();
    plastic_strain_increment.clear();

    // Kinematic back stress substracted
    noalias(kin_hard_stress_vector) = predictive_stress_vector - back_stress_vector;

    const double threshold_indicator = TConstLawIntegratorType::CalculatePlasticParameters(
        kin_hard_stress_vector, r_strain_vector, uniaxial_stress,
        threshold, plastic_denominator, f_flux, g_flux,
        plastic_dissipation, plastic_strain_increment,
        r_constitutive_matrix, rValues, characteristic_length,
        plastic_strain, back_stress_vector);

    if (threshold_indicator > std::abs(1.0e-4 * threshold)) {
        // while loop backward euler
        /* Inside "IntegrateStressVector" the predictive_stress_vector is updated to verify the yield criterion */
        TConstLawIntegratorType::IntegrateStressVector(
            predictive_stress_vector, r_strain_vector, uniaxial_stress,
            threshold, plastic_denominator, f_flux, g_flux,
            plastic_dissipation, plastic_strain_increment,
            r_constitutive_matrix, plastic_strain, rValues,
            characteristic_length, back_stress_vector,
            previous_stress_vector);
    }

    mPlasticDissipation = plastic_dissipation;
    mThreshold = threshold;

    noalias(mPlasticStrain) = plastic_strain;
    noalias(mPreviousStressVector) = predictive_stress_vector;
    noalias(mBackStressVector) = back_stress_vector;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
bool GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::Has(const Variable<double>& rThisVariable)
{
    if (rThisVariable == PLASTIC_DISSIPATION) {
        return true;
    }
    return BaseType::Has(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
bool GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::Has(const Variable<Vector>& rThisVariable)
{
    if (rThisVariable == INTERNAL_VARIABLES) {
        return true;
    }
    if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
        return true;
    }
    return BaseType::Has(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
bool GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::Has(const Variable<Matrix>& rThisVariable)
{
    if (rThisVariable == PLASTIC_STRAIN_TENSOR) {
        return true;
    } else if (rThisVariable == INTEGRATED_STRESS_TENSOR) {
        return true;
    } else {
        return BaseType::Has(rThisVariable);
    }
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rThisVariable == PLASTIC_DISSIPATION) {
        mPlasticDissipation = rValue;
    } else {
        BaseType::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::SetValue(
    const Variable<Vector>& rThisVariable,
    const Vector& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (rThisVariable == INTERNAL_VARIABLES) {
        mPlasticDissipation = rValue[0];
        for (std::size_t i=0; i < VoigtSize; ++i)
            mPlasticStrain[i] = rValue[i + 1];
        return;
    }
    if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
        mPlasticStrain = rValue;
        return;
    }
    BaseType::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
double& GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if (rThisVariable == PLASTIC_DISSIPATION) {
        rValue = mPlasticDissipation;
    } else {
        BaseType::GetValue(rThisVariable, rValue);
    }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
Vector& GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    if (rThisVariable == INTERNAL_VARIABLES) {
        rValue.resize(1 + VoigtSize);
        rValue[0] = mPlasticDissipation;
        for (std::size_t i=0; i < VoigtSize; ++i)
            rValue[i + 1] = mPlasticStrain[i];
        return rValue;
    }
    if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
        rValue.resize(VoigtSize, false);
        noalias(rValue) = mPlasticStrain;
        return rValue;
    }
    return BaseType::GetValue(rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
Matrix& GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::GetValue(
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    rValue.resize(Dimension, Dimension, false);
    if (rThisVariable == PLASTIC_STRAIN_TENSOR) {
        noalias(rValue) = MathUtils<double>::StrainVectorToTensor(mPlasticStrain);
    } else if (rThisVariable == INTEGRATED_STRESS_TENSOR) {
        noalias(rValue) = MathUtils<double>::StressVectorToTensor(mPreviousStressVector);
    } else {
        return BaseType::GetValue(rThisVariable, rValue);
    }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
double& GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if (rThisVariable == UNIAXIAL_STRESS) {
        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS );

        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        // Calculate the stress vector
        CalculateMaterialResponseCauchy(rParameterValues);
        const Vector& r_stress_vector = rParameterValues.GetStressVector();
        const Vector& r_strain_vector = rParameterValues.GetStrainVector();

        BoundedArrayType aux_stress_vector = r_stress_vector;
        TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress( aux_stress_vector, r_strain_vector, rValue, rParameterValues);

        // Previous flags restored
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
    } else if (rThisVariable == EQUIVALENT_PLASTIC_STRAIN) {
        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
        const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS );

        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        // Calculate the stress vector
        CalculateMaterialResponseCauchy(rParameterValues);
        const Vector& r_stress_vector = rParameterValues.GetStressVector();

        // Previous flags restored
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);

        // Compute the equivalent plastic strain
        double uniaxial_stress;
        this->CalculateValue(rParameterValues, UNIAXIAL_STRESS, uniaxial_stress);
        TConstLawIntegratorType::CalculateEquivalentPlasticStrain(r_stress_vector, uniaxial_stress, mPlasticStrain, 0.0, rParameterValues, rValue);
        return rValue;
    } else {
        return this->GetValue(rThisVariable, rValue);
    }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
Vector& GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    if (rThisVariable == BACK_STRESS_VECTOR) {
        rValue = mBackStressVector;
    } else {
        BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
Matrix& GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    if (rThisVariable == BACK_STRESS_TENSOR) {
        Vector aux_value(VoigtSize);
        this->CalculateValue(rParameterValues, BACK_STRESS_VECTOR, aux_value);
        rValue = MathUtils<double>::StressVectorToTensor(aux_value);
    } else {
        return BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
int GenericSmallStrainKinematicPlasticity<TConstLawIntegratorType>::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    const int check_base = BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);
    const int check_integrator = TConstLawIntegratorType::Check(rMaterialProperties);
    KRATOS_ERROR_IF_NOT(VoigtSize == this->GetStrainSize()) << "You are combining not compatible constitutive laws" << std::endl;
    if ((check_base + check_integrator) > 0) return 1;
    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>>;

template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<PlaneStressVonMisesYieldSurface<PlaneStressVonMisesPlasticPotential<3>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<3>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<3>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<3>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<3>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<3>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<3>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<TrescaPlasticPotential<3>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<3>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<3>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<3>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<3>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<3>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<3>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<3>>>>;
template class GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<3>>>>;

} // namespace Kratos
