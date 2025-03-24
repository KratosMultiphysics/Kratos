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

// Project includes
#include "utilities/math_utils.h"
#include "constitutive_laws_application_variables.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"
#include "generic_small_strain_isotropic_viscoplasticity.h"
#include "custom_constitutive/auxiliary_files/cl_integrators/generic_cl_integrator_plasticity.h"

// Yield surfaces
// #include "custom_constitutive/auxiliary_files/yield_surfaces/generic_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/von_mises_yield_surface.h"
// #include "custom_constitutive/auxiliary_files/yield_surfaces/modified_mohr_coulomb_yield_surface.h"
// #include "custom_constitutive/auxiliary_files/yield_surfaces/mohr_coulomb_yield_surface.h"
// #include "custom_constitutive/auxiliary_files/yield_surfaces/rankine_yield_surface.h"
// #include "custom_constitutive/auxiliary_files/yield_surfaces/simo_ju_yield_surface.h"
// #include "custom_constitutive/auxiliary_files/yield_surfaces/drucker_prager_yield_surface.h"
// #include "custom_constitutive/auxiliary_files/yield_surfaces/tresca_yield_surface.h"

// Plastic potentials
// #include "custom_constitutive/auxiliary_files/plastic_potentials/generic_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/von_mises_plastic_potential.h"
// #include "custom_constitutive/auxiliary_files/plastic_potentials/tresca_plastic_potential.h"
// #include "custom_constitutive/auxiliary_files/plastic_potentials/modified_mohr_coulomb_plastic_potential.h"
// #include "custom_constitutive/auxiliary_files/plastic_potentials/mohr_coulomb_plastic_potential.h"
// #include "custom_constitutive/auxiliary_files/plastic_potentials/drucker_prager_plastic_potential.h"

namespace Kratos
{

template <class TConstLawIntegratorType>
void GenericSmallStrainIsotropicViscoPlasticity<TConstLawIntegratorType>::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();

    // We get the constitutive tensor
    Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();

    // Integrate Stress plasticity
    const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());

    if (r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        BaseType::CalculateCauchyGreenStrain(rValues, r_strain_vector);
    }
    this->template AddInitialStrainVectorContribution<Vector>(r_strain_vector);

    if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        CalculateElasticMatrix(r_constitutive_matrix, rValues);
    
        double threshold           = GetThreshold();
        double plastic_dissipation = GetPlasticDissipation();
        Vector plastic_strain      = GetPlasticStrain();
    
        BoundedArrayType predictive_stress_vector, deviatoric_stress_vector;
    
        noalias(predictive_stress_vector) = prod(r_constitutive_matrix, r_strain_vector - plastic_strain);
        this->template AddInitialStressVectorContribution<BoundedArrayType>(predictive_stress_vector);
    
        double equivalent_stress;
        ConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(predictive_stress_vector, r_strain_vector, equivalent_stress, rValues);
    
        const double F = equivalent_stress - threshold;
    
        if (F >= 0.0) {
            BoundedArrayType deviatoric_stress_vector;
            double I1, J2;
            ConstitutiveLawUtilities<VoigtSize>::CalculateI1Invariant<BoundedArrayType>(predictive_stress_vector, I1);
            ConstitutiveLawUtilities<VoigtSize>::CalculateJ2Invariant<BoundedArrayType>(predictive_stress_vector, I1, deviatoric_stress_vector, J2);
    
            const auto& r_props = rValues.GetMaterialProperties();
            const double mu = r_props[MIU];
            const double sensitivity = r_props[DP_EPSILON];
            const double plastic_multiplier = (std::pow(equivalent_stress / threshold, 1.0 / sensitivity) - 1.0) / mu;
    
            array_1d<double, VoigtSize> g_flux;
            ConstLawIntegratorType::YieldSurfaceType::CalculatePlasticPotentialDerivative(predictive_stress_vector, deviatoric_stress_vector, J2, g_flux, rValues);
    
            const array_1d<double, VoigtSize> plastic_strain_increment = plastic_multiplier * g_flux;
            noalias(rValues.GetStressVector()) = predictive_stress_vector - prod(r_constitutive_matrix, plastic_strain_increment);

            if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                BaseType::CalculateTangentTensor(rValues, plastic_strain); // this modifies the ConstitutiveMatrix
            }
        } else {
            noalias(rValues.GetStressVector()) = predictive_stress_vector;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainIsotropicViscoPlasticity<TConstLawIntegratorType>::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();

    // We get the constitutive tensor
    Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();

    // Integrate Stress plasticity
    const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());

    if (r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        BaseType::CalculateCauchyGreenStrain(rValues, r_strain_vector);
    }
    this->template AddInitialStrainVectorContribution<Vector>(r_strain_vector);

    CalculateElasticMatrix(r_constitutive_matrix, rValues);

    double& r_threshold           = GetThreshold();
    double& r_plastic_dissipation = GetPlasticDissipation();
    Vector& r_plastic_strain      = GetPlasticStrain();

    BoundedArrayType predictive_stress_vector, deviatoric_stress_vector;

    noalias(predictive_stress_vector) = prod(r_constitutive_matrix, r_strain_vector - r_plastic_strain);
    this->template AddInitialStressVectorContribution<BoundedArrayType>(predictive_stress_vector);

    double equivalent_stress;
    ConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(predictive_stress_vector, r_strain_vector, equivalent_stress, rValues);

    const double F = equivalent_stress - threshold;

    if (F >= 0.0) {
        BoundedArrayType deviatoric_stress_vector;
        double I1, J2;
        ConstitutiveLawUtilities<VoigtSize>::CalculateI1Invariant<BoundedArrayType>(predictive_stress_vector, I1);
        ConstitutiveLawUtilities<VoigtSize>::CalculateJ2Invariant<BoundedArrayType>(predictive_stress_vector, I1, deviatoric_stress_vector, J2);

        const auto& r_props = rValues.GetMaterialProperties();
        const double mu = r_props[MIU];
        const double sensitivity = r_props[DP_EPSILON];
        const double plastic_multiplier = (std::pow(equivalent_stress / threshold, 1.0 / sensitivity) - 1.0) / mu;

        array_1d<double, VoigtSize> g_flux;
        ConstLawIntegratorType::YieldSurfaceType::CalculatePlasticPotentialDerivative(predictive_stress_vector, deviatoric_stress_vector, J2, g_flux, rValues);

        const array_1d<double, VoigtSize> plastic_strain_increment = plastic_multiplier * g_flux;
        noalias(rValues.GetStressVector()) = predictive_stress_vector - prod(r_constitutive_matrix, plastic_strain_increment);

        const double g = r_props[FRACTURE_ENERGY] / characteristic_length;

        r_plastic_dissipation += inner_prod(rValues.GetStressVector(), plastic_strain_increment) / g;
        noalias(r_plastic_strain) += plastic_strain_increment;
        ConstLawIntegratorType::CalculateEquivalentStressThreshold(r_plastic_dissipation, 1, 0, r_threshold, 0, rValues, 0, characteristic_length);

    }
}

/***********************************************************************************/
/***********************************************************************************/

template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>>;

// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<3>>>>;

} // namespace Kratos
