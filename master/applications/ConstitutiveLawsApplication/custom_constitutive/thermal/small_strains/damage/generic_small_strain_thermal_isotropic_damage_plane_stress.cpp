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
//

// System includes

// External includes

// Project includes
// #include "constitutive_laws_application_variables.h"
#include "generic_small_strain_thermal_isotropic_damage_plane_stress.h"
#include "custom_constitutive/auxiliary_files/cl_integrators/generic_cl_integrator_damage.h"

// Yield surfaces
#include "custom_constitutive/thermal/auxiliary_files/thermal_yield_surfaces/thermal_von_mises_yield_surface.h"
#include "custom_constitutive/thermal/auxiliary_files/thermal_yield_surfaces/thermal_modified_mohr_coulomb_yield_surface.h"
#include "custom_constitutive/thermal/auxiliary_files/thermal_yield_surfaces/thermal_mohr_coulomb_yield_surface.h"
#include "custom_constitutive/thermal/auxiliary_files/thermal_yield_surfaces/thermal_rankine_yield_surface.h"
#include "custom_constitutive/thermal/auxiliary_files/thermal_yield_surfaces/thermal_simo_ju_yield_surface.h"
#include "custom_constitutive/thermal/auxiliary_files/thermal_yield_surfaces/thermal_drucker_prager_yield_surface.h"
#include "custom_constitutive/thermal/auxiliary_files/thermal_yield_surfaces/thermal_tresca_yield_surface.h"

// Plastic potentials
#include "custom_constitutive/auxiliary_files/plastic_potentials/generic_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/von_mises_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/tresca_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/modified_mohr_coulomb_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/mohr_coulomb_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/drucker_prager_plastic_potential.h"

namespace Kratos
{


/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainThermalIsotropicDamagePlaneStress<TConstLawIntegratorType>::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    const auto& r_cl_options = rValues.GetOptions();

    // We get the strain vector
    auto& r_strain_vector = rValues.GetStrainVector();

    // Since small strains are assumed, any measure of strains is valid
    if (r_cl_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        BaseType::CalculateCauchyGreenStrain(rValues, r_strain_vector);
    }

    // We compute the stress
    if (r_cl_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        const double ref_temperature = this->GetReferenceTemperature();
        auto& r_integrated_stress_vector = rValues.GetStressVector();
        // Elastic Matrix
        auto& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        const double E = AdvCLutils::GetMaterialPropertyThroughAccessor(YOUNG_MODULUS, rValues);
        const double poisson_ratio = AdvCLutils::GetMaterialPropertyThroughAccessor(POISSON_RATIO, rValues);


        CLutils::CalculateElasticMatrixPlaneStress(r_constitutive_matrix, E, poisson_ratio);
        AdvCLutils::SubstractThermalStrain(r_strain_vector, ref_temperature, rValues);

        this->template AddInitialStrainVectorContribution<Vector>(r_strain_vector);

        // S0 = C:(E-E0) + S0
        BoundedArrayType predictive_stress_vector = prod(r_constitutive_matrix, r_strain_vector);
        this->template AddInitialStressVectorContribution<BoundedArrayType>(predictive_stress_vector);

        // Converged values
        double threshold = this->GetThreshold();
        double damage = this->GetDamage();

        double uniaxial_stress;
        TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(predictive_stress_vector, r_strain_vector, uniaxial_stress, rValues);

        // NOTE: think about yield surfaces defined with tensile and compressive yields....
        const double ref_yield = AdvCLutils::GetPropertyFromTemperatureTable(YIELD_STRESS, rValues, ref_temperature);
        const double current_yield = AdvCLutils::GetMaterialPropertyThroughAccessor(YIELD_STRESS, rValues);
        const double temperature_reduction_factor = current_yield / ref_yield;

        // We affect the stress by Temperature reduction factor
        uniaxial_stress /=  temperature_reduction_factor;
        const double F = uniaxial_stress - threshold;

        if (F <= threshold_tolerance) { // Elasticity
            noalias(r_integrated_stress_vector) = (1.0 - damage) * predictive_stress_vector;
            if (r_cl_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                r_constitutive_matrix *= (1.0 - damage);
            }
        } else { // Damage case
            const double characteristic_length = AdvCLutils::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());

            // This routine updates the PredictiveStress to verify the yield surf
            TConstLawIntegratorType::IntegrateStressVector(predictive_stress_vector, uniaxial_stress, damage, threshold, rValues, characteristic_length);

            // Updated Values
            noalias(r_integrated_stress_vector) = predictive_stress_vector;

            if (r_cl_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                this->CalculateTangentTensor(rValues);
            }
        }
    }

} // End CalculateMaterialResponseCauchy

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainThermalIsotropicDamagePlaneStress<TConstLawIntegratorType>::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    const auto& r_cl_options = rValues.GetOptions();

    // We get the strain vector
    auto& r_strain_vector = rValues.GetStrainVector();

    // Since small strains are assumed, any measure of strains is valid
    if (r_cl_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        BaseType::CalculateCauchyGreenStrain(rValues, r_strain_vector);
    }

    const double ref_temperature = this->GetReferenceTemperature();

    // Elastic Matrix
    auto& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
    const double E = AdvCLutils::GetMaterialPropertyThroughAccessor(YOUNG_MODULUS, rValues);
    const double poisson_ratio = AdvCLutils::GetMaterialPropertyThroughAccessor(POISSON_RATIO, rValues);

    CLutils::CalculateElasticMatrixPlaneStress(r_constitutive_matrix, E, poisson_ratio);
    AdvCLutils::SubstractThermalStrain(r_strain_vector, ref_temperature, rValues);

    this->template AddInitialStrainVectorContribution<Vector>(r_strain_vector);

    // S0 = C:(E-E0) + S0
    BoundedArrayType predictive_stress_vector = prod(r_constitutive_matrix, r_strain_vector);
    this->template AddInitialStressVectorContribution<BoundedArrayType>(predictive_stress_vector);

    double uniaxial_stress;
    TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(predictive_stress_vector, r_strain_vector, uniaxial_stress, rValues);

    // NOTE: think about yield surfaces defined with tensile and compressive yields....
    const double ref_yield = AdvCLutils::GetPropertyFromTemperatureTable(YIELD_STRESS, rValues, ref_temperature);
    const double current_yield = AdvCLutils::GetMaterialPropertyThroughAccessor(YIELD_STRESS, rValues);
    const double temperature_reduction_factor = current_yield / ref_yield;

    // We modify the references...
    double &r_damage = this->GetDamage();
    double &r_threshold = this->GetThreshold();

    // We affect the stress by Temperature reduction factor
    uniaxial_stress /=  temperature_reduction_factor;
    const double F = uniaxial_stress - r_threshold;

    if (F > threshold_tolerance) {
        const double characteristic_length = AdvCLutils::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());

        // This routine updates the PredictiveStress to verify the yield surf
        TConstLawIntegratorType::IntegrateStressVector(predictive_stress_vector, uniaxial_stress, r_damage, r_threshold, rValues, characteristic_length);
        r_threshold = uniaxial_stress;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template class GenericSmallStrainThermalIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ThermalVonMisesYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainThermalIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ThermalModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainThermalIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ThermalTrescaYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainThermalIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ThermalDruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainThermalIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ThermalRankineYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainThermalIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ThermalSimoJuYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainThermalIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ThermalMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>;

template class GenericSmallStrainThermalIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ThermalVonMisesYieldSurface<VonMisesPlasticPotential<3>>>>;
template class GenericSmallStrainThermalIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ThermalModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>;
template class GenericSmallStrainThermalIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ThermalTrescaYieldSurface<VonMisesPlasticPotential<3>>>>;
template class GenericSmallStrainThermalIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ThermalDruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>>;
template class GenericSmallStrainThermalIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ThermalRankineYieldSurface<VonMisesPlasticPotential<3>>>>;
template class GenericSmallStrainThermalIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ThermalSimoJuYieldSurface<VonMisesPlasticPotential<3>>>>;
template class GenericSmallStrainThermalIsotropicDamagePlaneStress <GenericConstitutiveLawIntegratorDamage<ThermalMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>;
} // namespace Kratos
