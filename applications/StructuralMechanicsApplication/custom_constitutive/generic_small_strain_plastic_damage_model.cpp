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
    const double characteristic_length = ConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLength(rValues.GetElementGeometry());
    Matrix& r_tangent_tensor = rValues.GetConstitutiveMatrix(); // todo modify after integration
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();

    //NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    if( r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN )) {
        this->CalculateValue(rValues, STRAIN, r_strain_vector);
    }

    // Elastic Matrix
    if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
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
        double hard_damage = 0.0;
        double hcapd = 0.0;
        double denominator;

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
                hard_damage, hcapd, undamaged_free_energy);

        // Verification threshold for the plastic-damage process
        if (plasticity_indicator >= std::abs(1.0e-4 * threshold_plasticity) || damage_indicator >= std::abs(1.0e-4 * threshold_damage)) {
            bool is_converged = false;
            int number_iteration = 0;
            const int max_iter = 100;

 
            // Integration loop
            while (!is_converged && number_iteration <= max_iter) {
                number_iteration++;

            // *****************************************************
                // Plastic Damage Case
                if (plasticity_indicator >= std::abs(1.0e-4 * threshold_plasticity) && damage_indicator >= std::abs(1.0e-4 * threshold_damage)) {
                    
                }









            // *****************************************************
            }
            // Updated Values
            noalias(r_integrated_stress_vector) = predictive_stress_vector;
            if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                this->CalculateTangentTensor(rValues);
            } 
        } else {
			noalias(r_integrated_stress_vector) = predictive_stress_vector;
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
    // // Integrate Stress Damage
    // Vector& r_integrated_stress_vector = rValues.GetStressVector();
    // const double characteristic_length = ConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLength(rValues.GetElementGeometry());
    // Matrix& r_tangent_tensor = rValues.GetConstitutiveMatrix(); // todo modify after integration
    // const Flags& r_constitutive_law_options = rValues.GetOptions();

    // // We get the strain vector
    // Vector& r_strain_vector = rValues.GetStrainVector();

    // //NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    // if( r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN )) {
    //     this->CalculateValue(rValues, STRAIN, r_strain_vector);
    // }

    // // Elastic Matrix
    // if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
    //     Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
    //     this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);
    // }

    // // We compute the stress
    // if(r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
    //     // Elastic Matrix
    //     Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
    //     this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);

    //     if (r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
    //         BaseType::CalculateCauchyGreenStrain( rValues, r_strain_vector);
    //     }

    //     // Converged values
    //     double threshold_plasticity = mThresholdPlasticity;
    //     double threshold_damage = mThresholdDamage;
    //     double damage = mDamage;
    //     double plastic_dissipation = mPlasticDissipation;
    //     Vector plastic_strain = mPlasticStrain;
    //     double damage_dissipation = mDamageDissipation;
    //     double damage_increment = 0.0;  // dDamage
    //     double plastic_consistency_increment = 0.0; // dlambda
    //     double hard_damage = 0.0;
    //     double hcapd = 0.0;
    //     double denominator;

    //     // Stress Predictor S = (1-d)C:(E-Ep)
    //     array_1d<double, VoigtSize> effective_predictive_stress_vector = prod(r_constitutive_matrix, r_strain_vector - plastic_strain);
	// 	array_1d<double, VoigtSize> predictive_stress_vector = (1.0 - damage) * effective_predictive_stress_vector;

    //     // Initialize Plastic Parameters
    //     double uniaxial_stress_plasticity = 0.0, plastic_denominator = 0.0, uniaxial_stress_damage = 0.0;
    //     BoundedArrayType damage_yield_flux = ZeroVector(VoigtSize); // DF/DS
    //     BoundedArrayType f_flux = ZeroVector(VoigtSize); // DF/DS
    //     BoundedArrayType g_flux = ZeroVector(VoigtSize); // DG/DS
    //     BoundedArrayType plastic_strain_increment = ZeroVector(VoigtSize);
    //     BoundedArrayType deepp = ZeroVector(VoigtSize);

    //     double undamaged_free_energy = 0.5 * inner_prod(r_strain_vector - plastic_strain, effective_predictive_stress_vector);

    //     // Compute the plastic parameters
    //     double plasticity_indicator = TPlasticityIntegratorType::CalculatePlasticParameters(
    //             predictive_stress_vector, r_strain_vector, uniaxial_stress_plasticity,
    //             threshold_plasticity, plastic_denominator, f_flux, g_flux,
    //             plastic_dissipation, plastic_strain_increment,
    //             r_constitutive_matrix, rValues, characteristic_length,
    //             plastic_strain);

    //     // Compute Damage Parameters
    //     double damage_indicator = this->CalculateDamageParameters(
    //             predictive_stress_vector, r_strain_vector,
    //             uniaxial_stress_damage, threshold_damage, 
    //             damage_dissipation, r_constitutive_matrix,
    //             rValues, characteristic_length, damage_yield_flux,
    //             plastic_strain, damage, damage_increment, 
    //             hard_damage, hcapd, undamaged_free_energy);

    //     // Verification threshold for the plastic-damage process
    //     if (plasticity_indicator >= std::abs(1.0e-4 * threshold_plasticity) && damage_indicator >= std::abs(1.0e-4 * threshold_damage)) {
    //         bool is_converged = false;
    //         int number_iteration = 0;
    //         const int max_iter = 100;

    //         // Integration loop
    //         while (!is_converged && number_iteration <= max_iter) {
    //             number_iteration++;

    //             // Damage case without plasticity
    //             if (plasticity_indicator < 1.0e-4 * threshold_plasticity) { 

    //                 if (plastic_consistency_increment > tolerance) {
    //                     this->CalculateIncrementsPlasticDamageCase(
    //                         damage_yield_flux, effective_predictive_stress_vector,
    //                         damage, f_flux, g_flux, r_constitutive_matrix,
    //                         uniaxial_stress_plasticity, damage_indicator,
    //                         plasticity_indicator, hard_damage, plastic_denominator,
    //                         hcapd, damage_increment, plastic_consistency_increment);
    //                 } else {
    //                     plastic_consistency_increment = 0.0;
    //                     denominator = hard_damage;

    //                     denominator += inner_prod(damage_yield_flux, effective_predictive_stress_vector);
    //                     damage_increment = damage_indicator / denominator;
    //                 }


    //             // Plasticity case without damage
    //             } else if (damage_indicator < std::abs(1.0e-4 * threshold_damage)) { 
    //                 if (damage_increment > 0.0) {
    //                     this->CalculateIncrementsPlasticDamageCase(
    //                         damage_yield_flux, effective_predictive_stress_vector,
    //                         damage, f_flux, g_flux, r_constitutive_matrix,
    //                         uniaxial_stress_plasticity, damage_indicator,
    //                         plasticity_indicator, hard_damage, plastic_denominator,
    //                         hcapd, damage_increment, plastic_consistency_increment);
    //                 } else {
    //                     damage_increment = 0.0;
    //                     plastic_consistency_increment = plasticity_indicator * plastic_denominator;
    //                 }

    //             // Plastic-Damage case
    //             } else {
    //                 if (hard_damage == 0.0) {
    //                     damage_increment = 0.0;
    //                     plastic_consistency_increment = plasticity_indicator * plastic_denominator;
    //                 } else {
    //                     this->CalculateIncrementsPlasticDamageCase(
    //                         damage_yield_flux, effective_predictive_stress_vector,
    //                         damage, f_flux, g_flux, r_constitutive_matrix,
    //                         uniaxial_stress_plasticity, damage_indicator,
    //                         plasticity_indicator, hard_damage, plastic_denominator,
    //                         hcapd, damage_increment, plastic_consistency_increment);
    //                 }
    //             }
                
    //             // Update internal variables damage
    //             if (damage_increment > tolerance) damage += damage_increment;
    //             this->CheckInternalVariable(damage); // Just check te upper-lower bounds

    //             // Update internals variables plasticity
    //             if (plastic_consistency_increment > tolerance) plastic_strain_increment = plastic_consistency_increment * g_flux;
    //             else plastic_consistency_increment = 0.0;

    //             // plastic_strain_increment = std::abs(plastic_consistency_increment) * g_flux;

    //             noalias(plastic_strain) += plastic_strain_increment;
    //             noalias(deepp) += plastic_strain_increment;
    //             array_1d<double, VoigtSize> delta_sigma = prod(r_constitutive_matrix, plastic_strain_increment);

    //             // Return mapping
    //             noalias(effective_predictive_stress_vector) -= delta_sigma;
    //             noalias(predictive_stress_vector) = (1.0 - damage) * effective_predictive_stress_vector;
    //             undamaged_free_energy = 0.5 * inner_prod(r_strain_vector - plastic_strain - deepp, effective_predictive_stress_vector);
                
    //             // Verification to check wether we are inside the yield surfaces
    //             plasticity_indicator = TPlasticityIntegratorType::CalculatePlasticParameters(
    //                             predictive_stress_vector, r_strain_vector, uniaxial_stress_plasticity,
    //                             threshold_plasticity, plastic_denominator, f_flux, g_flux,
    //                             plastic_dissipation, plastic_strain_increment,
    //                             r_constitutive_matrix, rValues, characteristic_length,
    //                             plastic_strain);
    //             damage_indicator = this->CalculateDamageParameters(
    //                     predictive_stress_vector, r_strain_vector,
    //                     uniaxial_stress_damage, threshold_damage, 
    //                     damage_dissipation, r_constitutive_matrix,
    //                     rValues, characteristic_length, damage_yield_flux,
    //                     plastic_strain, damage, damage_increment, 
    //                     hard_damage, hcapd, undamaged_free_energy);

    //             // Final check
    //             if (plasticity_indicator < std::abs(1.0e-4 * threshold_plasticity) && damage_indicator < std::abs(1.0e-4 * threshold_damage)) {
    //                 is_converged = true; // There is convergence
    //             }
    //         }
    //         KRATOS_WARNING_IF("Backward Euler Plastic-Damage Model", number_iteration > max_iter) << "Maximum iterations in Plastic-Damage return mapping" << std::endl;

    //         // Updated Values
    //         noalias(r_integrated_stress_vector) = predictive_stress_vector;
    //     } else {
	// 		noalias(r_integrated_stress_vector) = predictive_stress_vector;
    //     }
    //     // Update internal variables
    //     mPlasticDissipation = plastic_dissipation;
    //     mThresholdPlasticity = threshold_plasticity;
    //     mPlasticStrain = plastic_strain;
    //     mThresholdDamage = threshold_damage;
    //     mDamage = damage;
    //     mDamageDissipation = damage_dissipation;
    // }

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
    array_1d<double, 6>& rFflux,
    const Vector& rPlasticStrain,
    const double Damage,
    const double DamageIncrement,
    double& rHardd,
    double& rHcapd,
    const double UndamagedFreeEnergy
)
{
    array_1d<double, VoigtSize> deviator = ZeroVector(6);
    array_1d<double, VoigtSize> h_capa = ZeroVector(6);
    double J2, tensile_indicator_factor, compression_indicator_factor, slope, hardening_parameter, suma = 0.0;

    TDamageIntegratorType::YieldSurfaceType::CalculateEquivalentStress(rPredictiveStressVector, rStrainVector, rUniaxialStress, rValues);
    const double I1 = rPredictiveStressVector[0] + rPredictiveStressVector[1] + rPredictiveStressVector[2];
    ConstitutiveLawUtilities<VoigtSize>::CalculateJ2Invariant(rPredictiveStressVector, I1, deviator, J2);
    TDamageIntegratorType::YieldSurfaceType::CalculateYieldSurfaceDerivative(rPredictiveStressVector, deviator, J2, rFflux, rValues);
    this->CalculateIndicatorsFactors(rPredictiveStressVector, tensile_indicator_factor, compression_indicator_factor, suma);

    auto& r_matProps = rValues.GetMaterialProperties();
    const bool has_symmetric_yield_stress = r_matProps.Has(YIELD_STRESS);
    const double yield_compression = has_symmetric_yield_stress ? r_matProps[YIELD_STRESS] : r_matProps[YIELD_STRESS_COMPRESSION];
    const double yield_tension = has_symmetric_yield_stress ? r_matProps[YIELD_STRESS] : r_matProps[YIELD_STRESS_TENSION];
    const double yield_ratio = yield_compression / yield_tension;

    const double fracture_energy_tension = r_matProps[FRACTURE_ENERGY_DAMAGE_PROCESS] / CharacteristicLength;
    const double fracture_energy_compression = fracture_energy_tension * std::pow(yield_ratio, 2.0);

    const double constant1 = tensile_indicator_factor * std::abs(rUniaxialStress / yield_ratio) / (fracture_energy_tension * suma);
    const double constant2 = compression_indicator_factor * std::abs(rUniaxialStress) / (fracture_energy_compression * suma);
    const double constant = constant1 + constant2;

    // Free Energy Undamaged
    rHcapd = constant * UndamagedFreeEnergy;
    double damage_dissipation_increment = rHcapd * DamageIncrement;
    this->CheckInternalVariable(damage_dissipation_increment);
    rDamageDissipation += damage_dissipation_increment;
    if (rDamageDissipation > 1.0) rDamageDissipation = 0.99999;

    Vector slopes(2), thresholds(2);
    for (IndexType i = 0; i < 2; ++i) {
        double initial_threshold;
        TDamageIntegratorType::GetInitialUniaxialThreshold(rValues, initial_threshold);
        thresholds[i] = initial_threshold * (1.0 - rDamageDissipation);
        slopes[i] = -initial_threshold;
    }
    rDamageThreshold = (tensile_indicator_factor * thresholds[0]) + (compression_indicator_factor * thresholds[1]);
    const double hsigr = rDamageThreshold * (tensile_indicator_factor * slopes[0] / thresholds[0] + compression_indicator_factor * slopes[1] / thresholds[1]);  
    rHardd = rHcapd * hsigr;

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
    double& rSumPrincipalStresses
)
{
    // We do an initial check
    if (norm_2(rPredictiveStressVector) < 1.0e-8) {
        rTensileIndicatorFactor = 0.0;
        rCompressionIndicatorFactor = 0.0;
        return;
    }

    // We proceed as usual
    array_1d<double, Dimension> principal_stresses = ZeroVector(Dimension);
    ConstitutiveLawUtilities<VoigtSize>::CalculatePrincipalStresses(principal_stresses, rPredictiveStressVector);

    double suma = 0.0, sumb = 0.0, sumc = 0.0;
    double aux_sa;

    for (IndexType i = 0; i < Dimension; ++i) {
        aux_sa = std::abs(principal_stresses[i]);
        suma += aux_sa;
        sumb += 0.5 * (principal_stresses[i] + aux_sa);
        sumc += 0.5 * (-principal_stresses[i] + aux_sa);
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
    double& rPlasticConsistencyIncrement
)
{
    const Vector dS_ddam = -prod(rElasticMatrix, rStrainVector - rPlasticStrain - rPlasticConsistencyIncrement*rPlasticityGFlux);
	const Vector dS_dlambda = -(1.0 - Damage - rDamageIncrement)*prod(rElasticMatrix, rPlasticityGFlux);
    const double dFp_dlambda = inner_prod(rPlasticityFlux, dS_dlambda); // A
    const double dFp_ddam = inner_prod(rPlasticityFlux, dS_ddam); // B
    const double dFd_dlamba = inner_prod(rFluxDamageYield, dS_dlambda); // C
    const double dFd_ddam = inner_prod(rFluxDamageYield, dS_ddam); // D
    const double jacobian_determinant = dFp_dlambda * dFd_ddam - dFp_ddam * dFd_dlamba;

    // todo try mcaully brackets

    rPlasticConsistencyIncrement -=  (dFd_ddam * PlasticityIndicator - dFp_ddam * DamageIndicator) / jacobian_determinant;
    rDamageIncrement -= (dFp_dlambda * DamageIndicator - dFd_dlamba * PlasticityIndicator) / jacobian_determinant;
}

/***********************************************************************************/
/***********************************************************************************/
template class GenericSmallStrainPlasticDamageModel<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>;


} // namespace Kratos
