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
#include "custom_constitutive/generic_small_strain_plastic_damage_model.h"
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_plasticity.h"
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_damage.h"

// Yield surfaces
#include "custom_constitutive/yield_surfaces/generic_yield_surface.h"
#include "custom_constitutive/yield_surfaces/von_mises_yield_surface.h"
#include "custom_constitutive/yield_surfaces/modified_mohr_coulomb_yield_surface.h"
#include "custom_constitutive/yield_surfaces/mohr_coulomb_yield_surface.h"
#include "custom_constitutive/yield_surfaces/rankine_yield_surface.h"
#include "custom_constitutive/yield_surfaces/simo_ju_yield_surface.h"
#include "custom_constitutive/yield_surfaces/drucker_prager_yield_surface.h"
#include "custom_constitutive/yield_surfaces/tresca_yield_surface.h"

// Plastic potentials
#include "custom_constitutive/plastic_potentials/generic_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/von_mises_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/tresca_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/modified_mohr_coulomb_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/mohr_coulomb_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/drucker_prager_plastic_potential.h"

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
    // const double characteristic_length = ConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLength(rValues.GetElementGeometry());
    const double characteristic_length = 0.01; // todo
    Matrix& r_tangent_tensor = rValues.GetConstitutiveMatrix(); // todo modify after integration
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
    if(r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        // Elastic Matrix
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);
		

        if (r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
            BaseType::CalculateCauchyGreenStrain( rValues, r_strain_vector);
        }

        // Converged values
        double threshold_plasticity = mThresholdPlasticity;
        double threshold_damage = mThresholdDamage;
        double damage = mDamage;
        double plastic_dissipation = mPlasticDissipation;
        Vector plastic_strain = mPlasticStrain;
        double damage_dissipation = mDamageDissipation;
        double damage_increment = 0.0;  // dDamage
        double plastic_consistency_increment = 0.0; // dlambda
        double hardd = 0.0, damage_dissipation_increment = 0.0;

        // Stress Predictor S = (1-d)C:(E-Ep)
        array_1d<double, VoigtSize> effective_predictive_stress_vector = prod(r_constitutive_matrix, r_strain_vector - plastic_strain);
		array_1d<double, VoigtSize> predictive_stress_vector = (1.0 - damage) * effective_predictive_stress_vector;

        // Initialize Plastic Parameters
        double uniaxial_stress_plasticity = 0.0, plastic_denominator = 0.0, uniaxial_stress_damage = 0.0;
        BoundedArrayType damage_yield_flux = ZeroVector(VoigtSize); // DF/DS
        BoundedArrayType f_flux = ZeroVector(VoigtSize); // DF/DS
        BoundedArrayType g_flux = ZeroVector(VoigtSize); // DG/DS
        BoundedArrayType plastic_strain_increment = ZeroVector(VoigtSize);
        BoundedArrayType deepp = ZeroVector(VoigtSize);

        double undamaged_free_energy = 0.5 * inner_prod(r_strain_vector - plastic_strain, effective_predictive_stress_vector);

        // Compute the plastic parameters
        double plasticity_indicator = TPlasticityIntegratorType::CalculatePlasticParameters(
                predictive_stress_vector, r_strain_vector, uniaxial_stress_plasticity,
                threshold_plasticity, plastic_denominator, f_flux, g_flux,
                plastic_dissipation, plastic_strain_increment,
                r_constitutive_matrix, rValues, characteristic_length,
                plastic_strain);

        // Compute Damage Parameters
        double damage_indicator = this->CalculateDamageParameters(
                predictive_stress_vector, r_strain_vector,
                uniaxial_stress_damage, threshold_damage, 
                damage_dissipation, r_constitutive_matrix,
                rValues, characteristic_length, damage_yield_flux,
                plastic_strain, damage, damage_increment, 
                undamaged_free_energy, hardd, damage_dissipation_increment);

        // Verification threshold for the plastic-damage process
        if (plasticity_indicator >= std::abs(1.0e-4 * threshold_plasticity) && damage_indicator >= std::abs(1.0e-4 * threshold_damage)) {
            bool is_converged = false;
            int number_iteration = 0;
            const int max_iter = 100;
			
            // Integration loop
            while (!is_converged && number_iteration <= max_iter) {
                // Plastic case
				if (damage_indicator <= std::abs(1.0e-4 * threshold_damage)) {
                    if (damage_increment > tolerance) {
                        this->CalculateIncrementsPlasticDamageCase(
                                damage_yield_flux, r_strain_vector, damage,
                                f_flux, g_flux, r_constitutive_matrix, 
                                damage_indicator, plasticity_indicator, 
                                plastic_strain, damage_increment, 
                                plastic_consistency_increment, plastic_denominator, uniaxial_stress_plasticity, hardd, damage_dissipation_increment, 1);  
                    } else {
                        damage_increment = 0.0;
                        plastic_consistency_increment = plasticity_indicator * plastic_denominator;
                    }
                // Damage case
                } else if (plasticity_indicator <= std::abs(1.0e-4 * threshold_plasticity)) {
                    if (plastic_consistency_increment > tolerance) {
                        this->CalculateIncrementsPlasticDamageCase(
                                damage_yield_flux, r_strain_vector, damage,
                                f_flux, g_flux, r_constitutive_matrix, 
                                damage_indicator, plasticity_indicator, 
                                plastic_strain, damage_increment, 
                                plastic_consistency_increment, plastic_denominator, uniaxial_stress_plasticity, hardd, damage_dissipation_increment, 1);  
                    } else {
                        const double denom = hardd + inner_prod(damage_yield_flux, effective_predictive_stress_vector);
                        damage_increment = damage_indicator / denom;
                        plastic_consistency_increment = 0.0;
                    }
                } else { // Plastic Damage Case
                    if (std::abs(hardd) < tolerance) {
                        damage_increment = 0.0;
                        plastic_consistency_increment = plasticity_indicator * plastic_denominator;  
                    } else {
                        this->CalculateIncrementsPlasticDamageCase(
                                damage_yield_flux, r_strain_vector, damage,
                                f_flux, g_flux, r_constitutive_matrix, 
                                damage_indicator, plasticity_indicator, 
                                plastic_strain, damage_increment, 
                                plastic_consistency_increment, plastic_denominator, uniaxial_stress_plasticity, hardd, damage_dissipation_increment, 1);
                    }
				} // Increments computed

                if (damage_increment > tolerance) damage += damage_increment;
                this->CheckInternalVariable(damage);
                if (plastic_consistency_increment > tolerance) noalias(plastic_strain_increment) = plastic_consistency_increment * g_flux;
                noalias(deepp) += plastic_strain_increment;
                noalias(plastic_strain) += deepp;

                effective_predictive_stress_vector -= prod(r_constitutive_matrix, plastic_strain_increment);
                predictive_stress_vector = (1.0 - damage) * effective_predictive_stress_vector;

                undamaged_free_energy = 0.5 * inner_prod(r_strain_vector - plastic_strain, effective_predictive_stress_vector);

                // Compute the plastic parameters
                plasticity_indicator = TPlasticityIntegratorType::CalculatePlasticParameters(
                        predictive_stress_vector, r_strain_vector, uniaxial_stress_plasticity,
                        threshold_plasticity, plastic_denominator, f_flux, g_flux,
                        plastic_dissipation, plastic_strain_increment,
                        r_constitutive_matrix, rValues, characteristic_length,
                        plastic_strain);

                // Compute Damage Parameters
                damage_indicator = this->CalculateDamageParameters(
                        predictive_stress_vector, r_strain_vector,
                        uniaxial_stress_damage, threshold_damage, 
                        damage_dissipation, r_constitutive_matrix,
                        rValues, characteristic_length, damage_yield_flux,
                        plastic_strain, damage, damage_increment, 
                        undamaged_free_energy, hardd, damage_dissipation_increment);

                if (plasticity_indicator < std::abs(1.0e-4 * threshold_plasticity) && damage_indicator < std::abs(1.0e-4 * threshold_damage)) {
                    is_converged = true;
                } else {
                    number_iteration++;
                }
            }
            KRATOS_WARNING_IF("Backward Euler Plastic Damage", number_iteration >= max_iter) << "Max iterations reached in the return mapping of the Plastic Damage model" << std::endl; 
            // Updated Values
            noalias(r_integrated_stress_vector) = predictive_stress_vector;
            if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                this->CalculateTangentTensor(rValues);
            } 
        } else {
			// noalias(r_integrated_stress_vector) = predictive_stress_vector;
			noalias(r_integrated_stress_vector) = (1.0 - damage) * effective_predictive_stress_vector;
            if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                noalias(r_tangent_tensor) = (1.0 - damage) * r_constitutive_matrix;
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
    // const double characteristic_length = ConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLength(rValues.GetElementGeometry());
    const double characteristic_length = 0.01; // todo
    Matrix& r_tangent_tensor = rValues.GetConstitutiveMatrix(); // todo modify after integration
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
    if(r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        // Elastic Matrix
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);
		

        if (r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
            BaseType::CalculateCauchyGreenStrain( rValues, r_strain_vector);
        }

        // Converged values
        double threshold_plasticity = mThresholdPlasticity;
        double threshold_damage = mThresholdDamage;
        double damage = mDamage;
        double plastic_dissipation = mPlasticDissipation;
        Vector plastic_strain = mPlasticStrain;
        double damage_dissipation = mDamageDissipation;
        double damage_increment = 0.0;  // dDamage
        double plastic_consistency_increment = 0.0; // dlambda
        double hardd = 0.0, damage_dissipation_increment = 0.0;

        // Stress Predictor S = (1-d)C:(E-Ep)
        array_1d<double, VoigtSize> effective_predictive_stress_vector = prod(r_constitutive_matrix, r_strain_vector - plastic_strain);
		array_1d<double, VoigtSize> predictive_stress_vector = (1.0 - damage) * effective_predictive_stress_vector;

        // Initialize Plastic Parameters
        double uniaxial_stress_plasticity = 0.0, plastic_denominator = 0.0, uniaxial_stress_damage = 0.0;
        BoundedArrayType damage_yield_flux = ZeroVector(VoigtSize); // DF/DS
        BoundedArrayType f_flux = ZeroVector(VoigtSize); // DF/DS
        BoundedArrayType g_flux = ZeroVector(VoigtSize); // DG/DS
        BoundedArrayType plastic_strain_increment = ZeroVector(VoigtSize);
        BoundedArrayType deepp = ZeroVector(VoigtSize);

        double undamaged_free_energy = 0.5 * inner_prod(r_strain_vector - plastic_strain, effective_predictive_stress_vector);

        // Compute the plastic parameters
        double plasticity_indicator = TPlasticityIntegratorType::CalculatePlasticParameters(
                predictive_stress_vector, r_strain_vector, uniaxial_stress_plasticity,
                threshold_plasticity, plastic_denominator, f_flux, g_flux,
                plastic_dissipation, plastic_strain_increment,
                r_constitutive_matrix, rValues, characteristic_length,
                plastic_strain);

        // Compute Damage Parameters
        double damage_indicator = this->CalculateDamageParameters(
                predictive_stress_vector, r_strain_vector,
                uniaxial_stress_damage, threshold_damage, 
                damage_dissipation, r_constitutive_matrix,
                rValues, characteristic_length, damage_yield_flux,
                plastic_strain, damage, damage_increment, 
                undamaged_free_energy, hardd, damage_dissipation_increment);

        // Verification threshold for the plastic-damage process
        if (plasticity_indicator >= std::abs(1.0e-4 * threshold_plasticity) && damage_indicator >= std::abs(1.0e-4 * threshold_damage)) {
            bool is_converged = false;
            int number_iteration = 0;
            const int max_iter = 100;
			
            // Integration loop
            while (!is_converged && number_iteration <= max_iter) {
                // Plastic case
				if (damage_indicator <= std::abs(1.0e-4 * threshold_damage)) {
                    if (damage_increment > tolerance) {
                        this->CalculateIncrementsPlasticDamageCase(
                                damage_yield_flux, r_strain_vector, damage,
                                f_flux, g_flux, r_constitutive_matrix, 
                                damage_indicator, plasticity_indicator, 
                                plastic_strain, damage_increment, 
                                plastic_consistency_increment, plastic_denominator, uniaxial_stress_plasticity, hardd, damage_dissipation_increment, 1);  
                    } else {
                        damage_increment = 0.0;
                        plastic_consistency_increment = plasticity_indicator * plastic_denominator;
                    }
                // Damage case
                } else if (plasticity_indicator <= std::abs(1.0e-4 * threshold_plasticity)) {
                    if (plastic_consistency_increment > tolerance) {
                        this->CalculateIncrementsPlasticDamageCase(
                                damage_yield_flux, r_strain_vector, damage,
                                f_flux, g_flux, r_constitutive_matrix, 
                                damage_indicator, plasticity_indicator, 
                                plastic_strain, damage_increment, 
                                plastic_consistency_increment, plastic_denominator, uniaxial_stress_plasticity, hardd, damage_dissipation_increment, 1);  
                    } else {
                        const double denom = hardd + inner_prod(damage_yield_flux, effective_predictive_stress_vector);
                        damage_increment = damage_indicator / denom;
                        plastic_consistency_increment = 0.0;
                    }
                } else { // Plastic Damage Case
                    if (std::abs(hardd) < tolerance) {
                        damage_increment = 0.0;
                        plastic_consistency_increment = plasticity_indicator * plastic_denominator;  
                    } else {
                        this->CalculateIncrementsPlasticDamageCase(
                                damage_yield_flux, r_strain_vector, damage,
                                f_flux, g_flux, r_constitutive_matrix, 
                                damage_indicator, plasticity_indicator, 
                                plastic_strain, damage_increment, 
                                plastic_consistency_increment, plastic_denominator, uniaxial_stress_plasticity, hardd, damage_dissipation_increment, 1);
                    }
				} // Increments computed

                if (damage_increment > tolerance) damage += damage_increment;
                this->CheckInternalVariable(damage);
                if (plastic_consistency_increment > tolerance) noalias(plastic_strain_increment) = plastic_consistency_increment * g_flux;
                noalias(deepp) += plastic_strain_increment;
                noalias(plastic_strain) += deepp;

                effective_predictive_stress_vector -= prod(r_constitutive_matrix, plastic_strain_increment);
                predictive_stress_vector = (1.0 - damage) * effective_predictive_stress_vector;

                undamaged_free_energy = 0.5 * inner_prod(r_strain_vector - plastic_strain, effective_predictive_stress_vector);

                // Compute the plastic parameters
                plasticity_indicator = TPlasticityIntegratorType::CalculatePlasticParameters(
                        predictive_stress_vector, r_strain_vector, uniaxial_stress_plasticity,
                        threshold_plasticity, plastic_denominator, f_flux, g_flux,
                        plastic_dissipation, plastic_strain_increment,
                        r_constitutive_matrix, rValues, characteristic_length,
                        plastic_strain);

                // Compute Damage Parameters
                damage_indicator = this->CalculateDamageParameters(
                        predictive_stress_vector, r_strain_vector,
                        uniaxial_stress_damage, threshold_damage, 
                        damage_dissipation, r_constitutive_matrix,
                        rValues, characteristic_length, damage_yield_flux,
                        plastic_strain, damage, damage_increment, 
                        undamaged_free_energy, hardd, damage_dissipation_increment);

                if (plasticity_indicator < std::abs(1.0e-4 * threshold_plasticity) && damage_indicator < std::abs(1.0e-4 * threshold_damage)) {
                    is_converged = true;
                } else {
                    number_iteration++;
                }
            }
            KRATOS_WARNING_IF("Backward Euler Plastic Damage", number_iteration >= max_iter) << "Max iterations reached in the return mapping of the Plastic Damage model" << std::endl; 
            // Updated Values
            noalias(r_integrated_stress_vector) = predictive_stress_vector;
            if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                this->CalculateTangentTensor(rValues);
            } 
        } else {
			// noalias(r_integrated_stress_vector) = predictive_stress_vector;
			noalias(r_integrated_stress_vector) = (1.0 - damage) * effective_predictive_stress_vector;
            if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                noalias(r_tangent_tensor) = (1.0 - damage) * r_constitutive_matrix;
            }
        }
        // Update internal variables
        mPlasticDissipation = plastic_dissipation;
        mThresholdPlasticity = threshold_plasticity;
        mPlasticStrain = plastic_strain;
        mThresholdDamage = threshold_damage;
        mDamage = damage;
        mDamageDissipation = damage_dissipation;
        TPlasticityIntegratorType::YieldSurfaceType::CalculateEquivalentStress(predictive_stress_vector, r_strain_vector, mUniaxialStress, rValues);
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
    if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
        return true;
    } else {
        return BaseType::Has(rThisVariable);
    }

    return false;
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
    if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
        mPlasticStrain = rValue;
    } else {
        BaseType::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
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
    if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
        rValue = mPlasticStrain;
    } else {
        return BaseType::GetValue(rThisVariable, rValue);
    }

    return rValue;
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
	array_1d<double, 6>& rPredictiveStressVector,
	Vector& rStrainVector,
	double& rUniaxialStress,
	double& rDamageThreshold,
	double& rDamageDissipation,
	const Matrix& rConstitutiveMatrix,
	ConstitutiveLaw::Parameters& rValues,
	const double CharacteristicLength,
	array_1d<double, 6>& rDamageFlux,
	const Vector& rPlasticStrain,
	const double Damage,
	const double DamageIncrement,
	const double UndamagedFreeEnergy,
	double& rHardd,
    double& rDamageDissipationIncrement
)
{
    array_1d<double, VoigtSize> deviator = ZeroVector(6);
    array_1d<double, VoigtSize> h_capa = ZeroVector(6);
    double J2, tensile_indicator_factor, compression_indicator_factor, suma = 0.0;
    array_1d<double, 3> principal_stresses;

    TDamageIntegratorType::YieldSurfaceType::CalculateEquivalentStress(rPredictiveStressVector, rStrainVector, rUniaxialStress, rValues);
    const double I1 = rPredictiveStressVector[0] + rPredictiveStressVector[1] + rPredictiveStressVector[2];
    ConstitutiveLawUtilities<VoigtSize>::CalculateJ2Invariant(rPredictiveStressVector, I1, deviator, J2);
    TDamageIntegratorType::YieldSurfaceType::CalculateYieldSurfaceDerivative(rPredictiveStressVector, deviator, J2, rDamageFlux, rValues);
    this->CalculateIndicatorsFactors(rPredictiveStressVector, tensile_indicator_factor, compression_indicator_factor, suma, principal_stresses);

    auto& r_matProps = rValues.GetMaterialProperties();
    const bool has_symmetric_yield_stress = r_matProps.Has(YIELD_STRESS);
    const double yield_compression = has_symmetric_yield_stress ? r_matProps[YIELD_STRESS] : r_matProps[YIELD_STRESS_COMPRESSION];
    const double yield_tension = has_symmetric_yield_stress ? r_matProps[YIELD_STRESS] : r_matProps[YIELD_STRESS_TENSION];
    const double yield_ratio = yield_compression / yield_tension;

    const double fracture_energy_tension = r_matProps[FRACTURE_ENERGY_DAMAGE_PROCESS] / CharacteristicLength;
    const double fracture_energy_compression = fracture_energy_tension * std::pow(yield_ratio, 2.0);

    double const0 = 0.0, const1 = 0.0;
    if (std::abs(suma) > tolerance) {
        const0 = tensile_indicator_factor*(rUniaxialStress / yield_ratio) / (fracture_energy_tension * suma);
        const1 = compression_indicator_factor*rUniaxialStress / (fracture_energy_compression*suma);
    }
    double constant = const0 + const1;
    const double normalized_free_energy = constant*UndamagedFreeEnergy;
    rDamageDissipationIncrement = DamageIncrement * normalized_free_energy;

    this->CheckInternalVariable(rDamageDissipationIncrement);
    rDamageDissipation += rDamageDissipationIncrement;
    this->CheckInternalVariable(rDamageDissipation);

    Vector slopes(2), thresholds(2);
    // Tension
    thresholds[0] = yield_tension * (1.0 - rDamageDissipation);
    slopes[0] = -yield_tension;

    // Compression
    thresholds[1] = yield_compression * (1.0 - rDamageDissipation);
    slopes[1] = -yield_compression;

    rDamageThreshold = (tensile_indicator_factor * thresholds[0]) + (compression_indicator_factor * thresholds[1]);
    const double hsigr = rDamageThreshold * (tensile_indicator_factor * slopes[0] / thresholds[0] + compression_indicator_factor * slopes[1] / thresholds[1]);  
    rHardd = normalized_free_energy * hsigr;

    return rUniaxialStress - rDamageThreshold;
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
    if (rInternalVariable > 1.0) rInternalVariable = 0.99999;
    else if (rInternalVariable < tolerance) rInternalVariable = 0.0;
}

/***********************************************************************************/
/***********************************************************************************/
template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::
CalculateIncrementsPlasticDamageCase(
	const Vector& rFluxDamageYield,
	const Vector& rStrainVector,
	const double Damage,
	const Vector& rPlasticityFlux,
	const Vector& rPlasticityGFlux,
	const Matrix& rElasticMatrix,
	const double DamageIndicator,
	const double PlasticityIndicator,
	const Vector& rPlasticStrain,
	double& rDamageIncrement,
	double& rPlasticConsistencyIncrement,
	const double PlasticDenominator,
    const double UniaxialStressPlast,
    const double Hardd,
    const double& rDamageDissipationIncrement,
	const int PlasticDamageCase
)
{
    const Vector S0 = prod(rElasticMatrix, rStrainVector - rPlasticStrain);
    const Vector S = (1.0 - Damage) * S0;
    const double vs = inner_prod(rFluxDamageYield, S0);
	const double as = inner_prod(rPlasticityFlux, S0);
    double fact1 = 0.0, hkg = 0.0;

    const Vector HCAPA = S / UniaxialStressPlast;

    for (IndexType i = 0; i < 6; ++i) {
        double c = 0.0;
        for (IndexType j = 0; j < 6; ++j) {
            c += rFluxDamageYield[j] * rElasticMatrix(i, j);
        }
        hkg += HCAPA[i] * rPlasticityGFlux[i];
        fact1 += c * rPlasticityGFlux[i];
    }
    fact1 *= (1.0 - Damage);

    const double facta = vs;
    const double factb = 1.0 / (PlasticDenominator*(1.0- Damage));
    const double factc = as + Hardd;
    const double factd = fact1;
    const double denom = facta*factd - factb*factc;

    if (std::abs(denom) > tolerance) {
        rDamageIncrement = (factd*PlasticityIndicator - factb*DamageIndicator) / denom;
        rPlasticConsistencyIncrement = (facta*DamageIndicator - factc*PlasticityIndicator) / denom;
    } else {
        rDamageIncrement = PlasticityIndicator / (facta + factd*rDamageDissipationIncrement/hkg);
        rPlasticConsistencyIncrement = PlasticityIndicator / (factd + facta*hkg/rDamageDissipationIncrement);
    }
}

/***********************************************************************************/
/***********************************************************************************/
template <class TPlasticityIntegratorType, class TDamageIntegratorType>
double GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::
McaullyBrackets(
    const double Number
)
{
    return 0.5*(Number + std::abs(Number));
}

/***********************************************************************************/
/***********************************************************************************/
template class GenericSmallStrainPlasticDamageModel<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainPlasticDamageModel<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>;


} // namespace Kratos
