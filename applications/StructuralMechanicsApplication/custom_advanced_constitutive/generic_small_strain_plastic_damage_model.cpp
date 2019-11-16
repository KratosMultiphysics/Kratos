// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo & Sergio Jimenez
//
//

// System includes

// Project includes
#include "utilities/math_utils.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"
#include "custom_advanced_constitutive/generic_small_strain_plastic_damage_model.h"
#include "custom_advanced_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_plasticity.h"
#include "custom_advanced_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_damage.h"

// Yield surfaces
#include "custom_advanced_constitutive/yield_surfaces/generic_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/von_mises_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/modified_mohr_coulomb_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/mohr_coulomb_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/rankine_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/simo_ju_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/drucker_prager_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/tresca_yield_surface.h"

// Plastic potentials
#include "custom_advanced_constitutive/plastic_potentials/generic_plastic_potential.h"
#include "custom_advanced_constitutive/plastic_potentials/von_mises_plastic_potential.h"
#include "custom_advanced_constitutive/plastic_potentials/tresca_plastic_potential.h"
#include "custom_advanced_constitutive/plastic_potentials/modified_mohr_coulomb_plastic_potential.h"
#include "custom_advanced_constitutive/plastic_potentials/mohr_coulomb_plastic_potential.h"
#include "custom_advanced_constitutive/plastic_potentials/drucker_prager_plastic_potential.h"

namespace Kratos
{
template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // Integrate Stress Damage
    Vector& r_integrated_stress_vector = rValues.GetStressVector();
    const Flags& r_constitutive_law_options = rValues.GetOptions();
    Matrix& r_tangent_tensor = rValues.GetConstitutiveMatrix();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();

    //NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    if (r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        this->CalculateValue(rValues, STRAIN, r_strain_vector);
    }

    // Elastic Matrix
    if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);
    }
    // We compute the stress
    if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        // Elastic Matrix
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);


        if (r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
            BaseType::CalculateCauchyGreenStrain( rValues, r_strain_vector);
        }

        // Converged values
        PlasticDamageParameters cl_parameters = PlasticDamageParameters();
        cl_parameters.StrainVector = r_strain_vector;
        cl_parameters.PlasticityThreshold = mThresholdPlasticity;
        cl_parameters.DamageThreshold = mThresholdDamage;
        cl_parameters.Damage = mDamage;
        cl_parameters.PlasticDissipation = mPlasticDissipation;
        cl_parameters.PlasticStrain = mPlasticStrain;
        cl_parameters.DamageDissipation = mDamageDissipation;
        cl_parameters.PlasticConsistencyIncrement = 0.0;
        cl_parameters.CharacteristicLength = ConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLength(rValues.GetElementGeometry());

        // Stress Predictor S = (1-d)C:(E-Ep)
        array_1d<double, VoigtSize> effective_predictive_stress_vector = prod(r_constitutive_matrix, cl_parameters.StrainVector - cl_parameters.PlasticStrain);
		cl_parameters.StressVector = (1.0 - cl_parameters.Damage) * effective_predictive_stress_vector;

        cl_parameters.UndamagedFreeEnergy = 0.5 * inner_prod(cl_parameters.StrainVector - cl_parameters.PlasticStrain, effective_predictive_stress_vector);

        // Compute the plastic parameters
        cl_parameters.PlasticityIndicator = this->CalculatePlasticParameters(cl_parameters, r_constitutive_matrix, rValues);

        // Compute Damage Parameters
        cl_parameters.DamageIndicator = this->CalculateDamageParameters(cl_parameters, r_constitutive_matrix, rValues);

        // Verification threshold for the plastic-damage process
        if (cl_parameters.PlasticityIndicator >= std::abs(1.0e-4 * cl_parameters.PlasticityThreshold) &&
            cl_parameters.DamageIndicator >= std::abs(1.0e-4 * cl_parameters.DamageThreshold)) {
            bool is_converged = false;
            int number_iteration = 0;
            const int max_iter = 100;

            // Integration loop
            while (!is_converged && number_iteration <= max_iter) {
                // Plastic case
				if (cl_parameters.DamageIndicator <= std::abs(1.0e-4 * cl_parameters.DamageThreshold)) {
                    if (cl_parameters.DamageIncrement > tolerance) {
                        this->CalculateIncrementsPlasticDamageCase(cl_parameters, r_constitutive_matrix);
                    } else {
                        cl_parameters.DamageIncrement = 0.0;
                        cl_parameters.PlasticConsistencyIncrement = cl_parameters.PlasticityIndicator * cl_parameters.PlasticDenominator;
                    }
                // Damage case
                } else if (cl_parameters.PlasticityIndicator <= std::abs(1.0e-4 * cl_parameters.PlasticityThreshold)) {
                    if (cl_parameters.PlasticConsistencyIncrement > tolerance) {
                        this->CalculateIncrementsPlasticDamageCase(cl_parameters, r_constitutive_matrix);
                    } else {
                        const double denom = cl_parameters.HardeningParameterDamage + inner_prod(cl_parameters.DamageYieldFLux, effective_predictive_stress_vector);
                        cl_parameters.DamageIncrement = cl_parameters.DamageIndicator / denom;
                        cl_parameters.PlasticConsistencyIncrement = 0.0;
                    }
                } else { // Plastic Damage Case
                    if (std::abs(cl_parameters.HardeningParameterDamage) < tolerance) {
                        cl_parameters.DamageIncrement = 0.0;
                        cl_parameters.PlasticConsistencyIncrement = cl_parameters.PlasticityIndicator * cl_parameters.PlasticDenominator;
                    } else {
                        this->CalculateIncrementsPlasticDamageCase(cl_parameters, r_constitutive_matrix);
                    }
				} // Increments computed

                if (cl_parameters.DamageIncrement > tolerance) cl_parameters.Damage += cl_parameters.DamageIncrement;
                this->CheckInternalVariable(cl_parameters.Damage);
                if (cl_parameters.PlasticConsistencyIncrement > tolerance) noalias(cl_parameters.PlasticStrainIncrement) = cl_parameters.PlasticConsistencyIncrement * cl_parameters.PlasticityGFLux;
                noalias(cl_parameters.PlasticStrain) += cl_parameters.PlasticStrainIncrement;

                effective_predictive_stress_vector -= prod(r_constitutive_matrix, cl_parameters.PlasticStrainIncrement);
                cl_parameters.StressVector = (1.0 - cl_parameters.Damage) * effective_predictive_stress_vector;

                cl_parameters.UndamagedFreeEnergy = 0.5 * inner_prod(cl_parameters.StrainVector - cl_parameters.PlasticStrain, effective_predictive_stress_vector);

                // Compute the plastic parameters
                cl_parameters.PlasticityIndicator = this->CalculatePlasticParameters(cl_parameters, r_constitutive_matrix, rValues);

                // Compute Damage Parameters
                cl_parameters.DamageIndicator = this->CalculateDamageParameters(cl_parameters, r_constitutive_matrix, rValues);

                if (cl_parameters.PlasticityIndicator < std::abs(1.0e-4 * cl_parameters.PlasticityThreshold) &&
                    cl_parameters.DamageIndicator < std::abs(1.0e-4 * cl_parameters.DamageThreshold)) {
                    is_converged = true;
                } else {
                    number_iteration++;
                }
            }
            KRATOS_WARNING_IF("Backward Euler Plastic Damage", number_iteration >= max_iter) << "Max iterations reached in the return mapping of the Plastic Damage model" << std::endl;
            // Updated Values
            noalias(r_integrated_stress_vector) = cl_parameters.StressVector;
            if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                this->CalculateTangentTensor(rValues);
            }
        } else {
			noalias(r_integrated_stress_vector) = cl_parameters.StressVector;
            if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                noalias(r_tangent_tensor) = (1.0 - cl_parameters.Damage) * r_constitutive_matrix;
            }
        }
    }
} // End CalculateMaterialResponseCauchy

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::CalculateTangentTensor(ConstitutiveLaw::Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();

    const bool consider_perturbation_threshold = r_material_properties.Has(CONSIDER_PERTURBATION_THRESHOLD) ? r_material_properties[CONSIDER_PERTURBATION_THRESHOLD] : true;
    const TangentOperatorEstimation tangent_operator_estimation = r_material_properties.Has(TANGENT_OPERATOR_ESTIMATION) ? static_cast<TangentOperatorEstimation>(r_material_properties[TANGENT_OPERATOR_ESTIMATION]) : TangentOperatorEstimation::SecondOrderPerturbation;

    if (tangent_operator_estimation == TangentOperatorEstimation::Analytic) {
        KRATOS_ERROR << "Analytic solution not available" << std::endl;
    } else if (tangent_operator_estimation == TangentOperatorEstimation::FirstOrderPerturbation) {
        // Calculates the Tangent Constitutive Tensor by perturbation (first order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 1);
    } else if (tangent_operator_estimation == TangentOperatorEstimation::SecondOrderPerturbation) {
        // Calculates the Tangent Constitutive Tensor by perturbation (second order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 2);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    // We construct the CL parameters
    ProcessInfo dummy_process_info;
    ConstitutiveLaw::Parameters aux_param(rElementGeometry, rMaterialProperties, dummy_process_info);

    // // We call the integrator
    double initial_threshold_plast, initial_threshold_damage;
    TPlasticityIntegratorType::GetInitialUniaxialThreshold(aux_param, initial_threshold_plast);
    TDamageIntegratorType::GetInitialUniaxialThreshold(aux_param, initial_threshold_damage);
    mThresholdPlasticity = initial_threshold_plast;
    mThresholdDamage = initial_threshold_damage;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::InitializeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::InitializeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::InitializeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // Integrate Stress Damage
    Vector& r_integrated_stress_vector = rValues.GetStressVector();
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();

    //NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    if (r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        this->CalculateValue(rValues, STRAIN, r_strain_vector);
    }

    // Elastic Matrix
    if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);
    }
    // We compute the stress
    if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        // Elastic Matrix
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);


        if (r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
            BaseType::CalculateCauchyGreenStrain( rValues, r_strain_vector);
        }

        // Converged values
        PlasticDamageParameters cl_parameters = PlasticDamageParameters();
        cl_parameters.StrainVector = r_strain_vector;
        cl_parameters.PlasticityThreshold = mThresholdPlasticity;
        cl_parameters.DamageThreshold = mThresholdDamage;
        cl_parameters.Damage = mDamage;
        cl_parameters.PlasticDissipation = mPlasticDissipation;
        cl_parameters.PlasticStrain = mPlasticStrain;
        cl_parameters.DamageDissipation = mDamageDissipation;
        cl_parameters.PlasticConsistencyIncrement = 0.0;
        cl_parameters.CharacteristicLength = ConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLength(rValues.GetElementGeometry());

        // Stress Predictor S = (1-d)C:(E-Ep)
        array_1d<double, VoigtSize> effective_predictive_stress_vector = prod(r_constitutive_matrix, cl_parameters.StrainVector - cl_parameters.PlasticStrain);
		cl_parameters.StressVector = (1.0 - cl_parameters.Damage) * effective_predictive_stress_vector;

        cl_parameters.UndamagedFreeEnergy = 0.5 * inner_prod(cl_parameters.StrainVector - cl_parameters.PlasticStrain, effective_predictive_stress_vector);

        // Compute the plastic parameters
        cl_parameters.PlasticityIndicator = this->CalculatePlasticParameters(cl_parameters, r_constitutive_matrix, rValues);

        // Compute Damage Parameters
        cl_parameters.DamageIndicator = this->CalculateDamageParameters(cl_parameters, r_constitutive_matrix, rValues);

        // Verification threshold for the plastic-damage process
        if (cl_parameters.PlasticityIndicator >= std::abs(1.0e-4 * cl_parameters.PlasticityThreshold) &&
            cl_parameters.DamageIndicator >= std::abs(1.0e-4 * cl_parameters.DamageThreshold)) {
            bool is_converged = false;
            int number_iteration = 0;
            const int max_iter = 100;

            // Integration loop
            while (!is_converged && number_iteration <= max_iter) {
                // Plastic case
				if (cl_parameters.DamageIndicator <= std::abs(1.0e-4 * cl_parameters.DamageThreshold)) {
                    if (cl_parameters.DamageIncrement > tolerance) {
                        this->CalculateIncrementsPlasticDamageCase(cl_parameters, r_constitutive_matrix);
                    } else {
                        cl_parameters.DamageIncrement = 0.0;
                        cl_parameters.PlasticConsistencyIncrement = cl_parameters.PlasticityIndicator * cl_parameters.PlasticDenominator;
                    }
                // Damage case
                } else if (cl_parameters.PlasticityIndicator <= std::abs(1.0e-4 * cl_parameters.PlasticityThreshold)) {
                    if (cl_parameters.PlasticConsistencyIncrement > tolerance) {
                        this->CalculateIncrementsPlasticDamageCase(cl_parameters, r_constitutive_matrix);
                    } else {
                        const double denom = cl_parameters.HardeningParameterDamage + inner_prod(cl_parameters.DamageYieldFLux, effective_predictive_stress_vector);
                        cl_parameters.DamageIncrement = cl_parameters.DamageIndicator / denom;
                        cl_parameters.PlasticConsistencyIncrement = 0.0;
                    }
                } else { // Plastic Damage Case
                    if (std::abs(cl_parameters.HardeningParameterDamage) < tolerance) {
                        cl_parameters.DamageIncrement = 0.0;
                        cl_parameters.PlasticConsistencyIncrement = cl_parameters.PlasticityIndicator * cl_parameters.PlasticDenominator;
                    } else {
                        this->CalculateIncrementsPlasticDamageCase(cl_parameters, r_constitutive_matrix);
                    }
				} // Increments computed

                if (cl_parameters.DamageIncrement > tolerance) cl_parameters.Damage += cl_parameters.DamageIncrement;
                this->CheckInternalVariable(cl_parameters.Damage);
                if (cl_parameters.PlasticConsistencyIncrement > tolerance) noalias(cl_parameters.PlasticStrainIncrement) = cl_parameters.PlasticConsistencyIncrement * cl_parameters.PlasticityGFLux;
                noalias(cl_parameters.PlasticStrain) += cl_parameters.PlasticStrainIncrement;

                effective_predictive_stress_vector -= prod(r_constitutive_matrix, cl_parameters.PlasticStrainIncrement);
                cl_parameters.StressVector = (1.0 - cl_parameters.Damage) * effective_predictive_stress_vector;

                cl_parameters.UndamagedFreeEnergy = 0.5 * inner_prod(cl_parameters.StrainVector - cl_parameters.PlasticStrain, effective_predictive_stress_vector);

                // Compute the plastic parameters
                cl_parameters.PlasticityIndicator = this->CalculatePlasticParameters(cl_parameters, r_constitutive_matrix, rValues);

                // Compute Damage Parameters
                cl_parameters.DamageIndicator = this->CalculateDamageParameters(cl_parameters, r_constitutive_matrix, rValues);

                if (cl_parameters.PlasticityIndicator < std::abs(1.0e-4 * cl_parameters.PlasticityThreshold) &&
                    cl_parameters.DamageIndicator < std::abs(1.0e-4 * cl_parameters.DamageThreshold)) {
                    is_converged = true;
                } else {
                    number_iteration++;
                }
            }
            KRATOS_WARNING_IF("Backward Euler Plastic Damage", number_iteration >= max_iter) << "Max iterations reached in the return mapping of the Plastic Damage model" << std::endl;
            // Updated Values
            noalias(r_integrated_stress_vector) = cl_parameters.StressVector;
        } else {
			noalias(r_integrated_stress_vector) = cl_parameters.StressVector;
        }
        // Update internal variables
        mPlasticDissipation = cl_parameters.PlasticDissipation;
        mThresholdPlasticity = cl_parameters.PlasticityThreshold;
        mPlasticStrain = cl_parameters.PlasticStrain;
        mThresholdDamage = cl_parameters.DamageThreshold;
        mDamage = cl_parameters.Damage;
        mDamageDissipation = cl_parameters.DamageDissipation;
        TPlasticityIntegratorType::YieldSurfaceType::CalculateEquivalentStress(cl_parameters.StressVector, cl_parameters.StrainVector, mUniaxialStress, rValues);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
bool GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::Has(const Variable<double>& rThisVariable)
{
    if (rThisVariable == PLASTIC_DISSIPATION) {
        return true;
    } else if (rThisVariable == DAMAGE) {
        return true;
    } else if (rThisVariable == UNIAXIAL_STRESS) {
        return true;
    } else {
        return BaseType::Has(rThisVariable);
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
bool GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::Has(const Variable<Vector>& rThisVariable)
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

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
bool GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::Has(const Variable<Matrix>& rThisVariable)
{
    return BaseType::Has(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rThisVariable == PLASTIC_DISSIPATION) {
        mPlasticDissipation = rValue;
    } else if (rThisVariable == DAMAGE) {
        mDamage = rValue;
    } else if (rThisVariable == UNIAXIAL_STRESS) {
        mUniaxialStress = rValue;
    } else {
        BaseType::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::SetValue(
    const Variable<Vector>& rThisVariable,
    const Vector& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (rThisVariable == INTERNAL_VARIABLES) {
        mPlasticDissipation = rValue[0];
        mDamage = rValue[1];
        mUniaxialStress = rValue[2];
        for (std::size_t i=0; i < VoigtSize; ++i)
            mPlasticStrain[i] = rValue[i + 3];
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

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
double& GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if (rThisVariable == PLASTIC_DISSIPATION) {
        rValue = mPlasticDissipation;
    } else if (rThisVariable == DAMAGE) {
        rValue = mDamage;
    } else if (rThisVariable == UNIAXIAL_STRESS) {
        rValue = mUniaxialStress;
    } else {
        BaseType::GetValue(rThisVariable, rValue);
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
Vector& GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    if (rThisVariable == INTERNAL_VARIABLES) {
        rValue.resize(3 + VoigtSize);
        rValue[0] = mPlasticDissipation;
        rValue[1] = mDamage;
        rValue[2] = mUniaxialStress;
        for (std::size_t i=0; i < VoigtSize; ++i)
            rValue[i + 3] = mPlasticStrain[i];
        return rValue;
    }
    if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
        rValue = mPlasticStrain;
    }
    return BaseType::GetValue(rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
Matrix& GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::GetValue(
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    if (rThisVariable == PLASTIC_STRAIN_TENSOR) {
        rValue = MathUtils<double>::StrainVectorToTensor(mPlasticStrain);
    } else {
        return BaseType::GetValue(rThisVariable, rValue);
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
double& GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    return this->GetValue(rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
Vector& GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    return BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
Matrix& GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    if (this->Has(rThisVariable)) {
        return this->GetValue(rThisVariable, rValue);
    } else {
        return BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
    }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
int GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const int check_base = BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);
    const int check_integrator_plasticity = TPlasticityIntegratorType::Check(rMaterialProperties);
	const int check_integrator_damage = TDamageIntegratorType::Check(rMaterialProperties);
    KRATOS_ERROR_IF_NOT(VoigtSize == this->GetStrainSize()) << "You are combining not compatible constitutive laws" << std::endl;
    if ((check_base + check_integrator_plasticity + check_integrator_damage) > 0) return 1;
    return 0;
}

/***********************************************************************************/
/***********************************************************************************/
template <class TPlasticityIntegratorType, class TDamageIntegratorType>
double GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::
CalculateDamageParameters(
    PlasticDamageParameters& rParameters,
    const Matrix& rElasticMatrix,
    ConstitutiveLaw::Parameters& rValues
)
{
    array_1d<double, VoigtSize> deviator = ZeroVector(6);
    array_1d<double, VoigtSize> h_capa = ZeroVector(6);
    double J2, tensile_indicator_factor, compression_indicator_factor, suma = 0.0;
    array_1d<double, Dimension> principal_stresses;

    TDamageIntegratorType::YieldSurfaceType::CalculateEquivalentStress(rParameters.StressVector, rParameters.StrainVector, rParameters.UniaxialStressDamage, rValues);
    const double I1 = rParameters.StressVector[0] + rParameters.StressVector[1] + rParameters.StressVector[2];
    ConstitutiveLawUtilities<VoigtSize>::CalculateJ2Invariant(rParameters.StressVector, I1, deviator, J2);
    TDamageIntegratorType::YieldSurfaceType::CalculateYieldSurfaceDerivative(rParameters.StressVector, deviator, J2, rParameters.DamageYieldFLux, rValues);
    this->CalculateIndicatorsFactors(rParameters.StressVector, tensile_indicator_factor, compression_indicator_factor, suma, principal_stresses);

    auto& r_matProps = rValues.GetMaterialProperties();
    const bool has_symmetric_yield_stress = r_matProps.Has(YIELD_STRESS);
    const double yield_compression = has_symmetric_yield_stress ? r_matProps[YIELD_STRESS] : r_matProps[YIELD_STRESS_COMPRESSION];
    const double yield_tension = has_symmetric_yield_stress ? r_matProps[YIELD_STRESS] : r_matProps[YIELD_STRESS_TENSION];
    const double yield_ratio = yield_compression / yield_tension;

    const double fracture_energy_tension = r_matProps[FRACTURE_ENERGY_DAMAGE_PROCESS] / rParameters.CharacteristicLength;
    const double fracture_energy_compression = fracture_energy_tension * std::pow(yield_ratio, 2.0);

    double const0 = 0.0, const1 = 0.0;
    if (std::abs(suma) > tolerance) {
        const0 = tensile_indicator_factor * (rParameters.UniaxialStressDamage / yield_ratio) / (fracture_energy_tension * suma);
        const1 = compression_indicator_factor * rParameters.UniaxialStressDamage / (fracture_energy_compression*suma);
    }
    double constant = const0 + const1;
    const double normalized_free_energy = constant * rParameters.UndamagedFreeEnergy;
    rParameters.DamageDissipationIncrement = rParameters.DamageIncrement * normalized_free_energy;

    this->CheckInternalVariable(rParameters.DamageDissipationIncrement);
    rParameters.DamageDissipation += rParameters.DamageDissipationIncrement;
    this->CheckInternalVariable(rParameters.DamageDissipation);

    Vector slopes(2), thresholds(2);
    // Tension
    thresholds[0] = yield_tension * (1.0 - rParameters.DamageDissipation);
    slopes[0] = -yield_tension;

    // Compression
    thresholds[1] = yield_compression * (1.0 - rParameters.DamageDissipation);
    slopes[1] = -yield_compression;

    rParameters.DamageThreshold = (tensile_indicator_factor * thresholds[0]) + (compression_indicator_factor * thresholds[1]);
    const double hsigr = rParameters.DamageThreshold * (tensile_indicator_factor * slopes[0] / thresholds[0] + compression_indicator_factor * slopes[1] / thresholds[1]);
    rParameters.HardeningParameterDamage = normalized_free_energy * hsigr;

    return rParameters.UniaxialStressDamage - rParameters.DamageThreshold;
}


/***********************************************************************************/
/***********************************************************************************/
template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::
CalculateIndicatorsFactors(
    const array_1d<double, 6>& rPredictiveStressVector,
    double& rTensileIndicatorFactor,
    double& rCompressionIndicatorFactor,
    double& rSumPrincipalStresses,
    array_1d<double, 3>& rPrincipalStresses
)
{
    // We do an initial check
    if (norm_2(rPredictiveStressVector) < 1.0e-8) {
        rTensileIndicatorFactor = 0.0;
        rCompressionIndicatorFactor = 0.0;
        return;
    }

    // We proceed as usual
    rPrincipalStresses = ZeroVector(Dimension);
    ConstitutiveLawUtilities<VoigtSize>::CalculatePrincipalStresses(rPrincipalStresses, rPredictiveStressVector);

    double suma = 0.0, sumb = 0.0, sumc = 0.0;
    double aux_sa;

    for (IndexType i = 0; i < Dimension; ++i) {
        aux_sa = std::abs(rPrincipalStresses[i]);
        suma += aux_sa;
        sumb += 0.5 * (rPrincipalStresses[i] + aux_sa);
        sumc += 0.5 * (-rPrincipalStresses[i] + aux_sa);
    }
    rSumPrincipalStresses = suma;

    if (std::abs(suma) > tolerance) {
        rTensileIndicatorFactor = sumb / suma;
        rCompressionIndicatorFactor = sumc / suma;
    } else {
        rTensileIndicatorFactor = sumb;
        rCompressionIndicatorFactor = sumc;
    }

    // Final check
    if ((std::abs(rTensileIndicatorFactor) + std::abs(rCompressionIndicatorFactor)) < tolerance) {
        rTensileIndicatorFactor = 0.0;
        rCompressionIndicatorFactor = 0.0;
        return;
    }
}

/***********************************************************************************/
/***********************************************************************************/
template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::
CheckInternalVariable(
    double& rInternalVariable
)
{
    if (rInternalVariable >= 1.0) rInternalVariable = 0.99999;
    else if (rInternalVariable < tolerance) rInternalVariable = 0.0;
}

/***********************************************************************************/
/***********************************************************************************/
template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::
CalculateIncrementsPlasticDamageCase(
    PlasticDamageParameters& rParameters,
	const Matrix& rElasticMatrix
)
{
    const Vector effective_stress_vector = prod(rElasticMatrix, rParameters.StrainVector - rParameters.PlasticStrain);
    const Vector stress_vector = (1.0 - rParameters.Damage) * effective_stress_vector;
    const double inner_fluxdamage_eff_stress = inner_prod(rParameters.DamageYieldFLux, effective_stress_vector);
	const double inner_fluxplast_eff_stress = inner_prod(rParameters.PlasticityFFLux, effective_stress_vector);
    double fact1 = 0.0, inner_normstress_gflux = 0.0;

    const Vector normalized_stress = stress_vector / rParameters.UniaxialStressPlasticity;

    for (IndexType i = 0; i < VoigtSize; ++i) {
        double c = 0.0;
        for (IndexType j = 0; j < 6; ++j) {
            c += rParameters.DamageYieldFLux[j] * rElasticMatrix(i, j);
        }
        inner_normstress_gflux += normalized_stress[i] * rParameters.PlasticityGFLux[i];
        fact1 += c * rParameters.PlasticityGFLux[i];
    }
    fact1 *= (1.0 - rParameters.Damage);

    const double facta = inner_fluxdamage_eff_stress;
    const double factb = 1.0 / rParameters.PlasticDenominator;
    const double factc = inner_fluxplast_eff_stress + rParameters.HardeningParameterDamage;
    const double factd = fact1;
    const double denominator = facta*factd - factb*factc;

    if (std::abs(denominator) > tolerance) {
        rParameters.DamageIncrement = (factd * rParameters.PlasticityIndicator - factb * rParameters.DamageIndicator) / denominator;
        rParameters.PlasticConsistencyIncrement = (facta * rParameters.DamageIndicator - factc * rParameters.PlasticityIndicator) / denominator;
    } else {
        rParameters.DamageIncrement = rParameters.PlasticityIndicator / (facta + factd * rParameters.DamageDissipationIncrement / inner_normstress_gflux);
        rParameters.PlasticConsistencyIncrement = rParameters.PlasticityIndicator / (factd + facta * inner_normstress_gflux / rParameters.DamageDissipationIncrement);
    }
}

/***********************************************************************************/
/***********************************************************************************/
template <class TPlasticityIntegratorType, class TDamageIntegratorType>
double GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::
CalculatePlasticParameters(
    PlasticDamageParameters& rParameters,
    const Matrix& rConstitutiveMatrix,
    ConstitutiveLaw::Parameters& rValues
)
{
        array_1d<double, VoigtSize> deviator = ZeroVector(6);
        array_1d<double, VoigtSize> h_capa = ZeroVector(6);
        double J2, tensile_indicator_factor, compression_indicator_factor, slope, hardening_parameter, equivalent_plastic_strain;

        TPlasticityIntegratorType::YieldSurfaceType::CalculateEquivalentStress(rParameters.StressVector, rParameters.StrainVector, rParameters.UniaxialStressPlasticity, rValues);
        const double I1 = rParameters.StressVector[0] + rParameters.StressVector[1] + rParameters.StressVector[2];
        ConstitutiveLawUtilities<VoigtSize>::CalculateJ2Invariant(rParameters.StressVector, I1, deviator, J2);
        TPlasticityIntegratorType::CalculateFFluxVector(rParameters.StressVector, deviator, J2, rParameters.PlasticityFFLux, rValues);
        TPlasticityIntegratorType::CalculateGFluxVector(rParameters.StressVector, deviator, J2, rParameters.PlasticityGFLux, rValues);
        TPlasticityIntegratorType::CalculateIndicatorsFactors(rParameters.StressVector, tensile_indicator_factor, compression_indicator_factor);
        TPlasticityIntegratorType::CalculatePlasticDissipation(rParameters.StressVector, tensile_indicator_factor, compression_indicator_factor, rParameters.PlasticStrainIncrement, rParameters.PlasticDissipation, h_capa, rValues, rParameters.CharacteristicLength);
        TPlasticityIntegratorType::CalculateEquivalentPlasticStrain(rParameters.StressVector, rParameters.UniaxialStressPlasticity, rParameters.PlasticStrain, tensile_indicator_factor, rValues, equivalent_plastic_strain);
        TPlasticityIntegratorType::CalculateEquivalentStressThreshold(rParameters.PlasticDissipation, tensile_indicator_factor, compression_indicator_factor, rParameters.PlasticityThreshold, slope, rValues, equivalent_plastic_strain);
        TPlasticityIntegratorType::CalculateHardeningParameter(rParameters.PlasticityFFLux, slope, h_capa, hardening_parameter);

        // Has to be slightly modified
        this->CalculatePlasticDenominator(rParameters.PlasticityFFLux, rParameters.PlasticityGFLux, rConstitutiveMatrix, hardening_parameter, rParameters.Damage, rParameters.PlasticDenominator);

        return rParameters.UniaxialStressPlasticity - rParameters.PlasticityThreshold;
}
/***********************************************************************************/
/***********************************************************************************/
template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::
CalculatePlasticDenominator(
    const array_1d<double, VoigtSize>& rFFlux,
    const array_1d<double, VoigtSize>& rGFlux,
    const Matrix& rConstitutiveMatrix,
    double& rHardeningParameter,
    const double Damage,
    double& rPlasticDenominator
)
{
    const array_1d<double, VoigtSize> delta_vector = prod(rGFlux, rConstitutiveMatrix);
    double A1 = 0.0;

    for (IndexType i = 0; i < VoigtSize; ++i) {
        A1 += rFFlux[i] * delta_vector[i];
    }
    A1 *= (1.0 - Damage);
    const double A2 = 0.0; // Only for isotropic hard
    const double A3 = rHardeningParameter;
    if (std::abs(A1 + A2 + A3) > tolerance)
        rPlasticDenominator = 1.0 / (A1 + A2 + A3);
    else {
        rPlasticDenominator = 1.0e-3 * std::numeric_limits<double>::max(); // TODO: Discuss this!!!
    }
}
/***********************************************************************************/
/***********************************************************************************/
template class GenericSmallStrainPlasticDamageModel<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainPlasticDamageModel<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>;


} // namespace Kratos
