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
#include "generic_small_strain_isotropic_cornejo_viscoplasticity.h"
#include "custom_constitutive/auxiliary_files/cl_integrators/generic_cl_integrator_plasticity.h"

// Yield surfaces
#include "custom_constitutive/auxiliary_files/yield_surfaces/generic_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/von_mises_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/modified_mohr_coulomb_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/mohr_coulomb_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/rankine_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/simo_ju_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/drucker_prager_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/tresca_yield_surface.h"

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
void GenericSmallStrainIsotropicCornejoViscoPlasticity<TConstLawIntegratorType>::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();

    // We get the constitutive tensor
    Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();

    if (r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        BaseType::CalculateCauchyGreenStrain(rValues, r_strain_vector);
    }
    this->template AddInitialStrainVectorContribution<Vector>(r_strain_vector);

    if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        CalculateElasticMatrix(r_constitutive_matrix, rValues);

        const auto &r_props = rValues.GetMaterialProperties();

        double threshold = GetThreshold();
        double plastic_dissipation = GetPlasticDissipation();
        Vector plastic_strain = GetPlasticStrain();
        const double time_regularization_factor = r_props.Has(TIME_REGULARIZATION) ? r_props[TIME_REGULARIZATION] : time_regularization;
        const double strain_rate_norm = time_regularization_factor * mStrainRateHistory[0] + (1.0 - time_regularization_factor) * mStrainRateHistory[1];

        BoundedArrayType predictive_stress_vector;

        noalias(predictive_stress_vector) = prod(r_constitutive_matrix, r_strain_vector - plastic_strain);
        this->template AddInitialStressVectorContribution<BoundedArrayType>(predictive_stress_vector);

        double equivalent_stress;
        CLIntegrator::YieldSurfaceType::CalculateEquivalentStress(predictive_stress_vector, r_strain_vector, equivalent_stress, rValues);

        const double viscous_eta = r_props[VISCOUS_PARAMETER];
        const double viscous_alpha = r_props.Has(VISCOUS_ALPHA) ? r_props[VISCOUS_ALPHA] : 1.0;
        const double viscous_beta = r_props.Has(VISCOUS_BETA) ? r_props[VISCOUS_BETA] : 1.0;
        const double time_delay = r_props[DELAY_TIME];
        const double yield_strain = r_props.Has(YIELD_STRESS) ? r_props[YIELD_STRESS] / r_props[YOUNG_MODULUS] : r_props[YIELD_STRESS_TENSION] / r_props[YOUNG_MODULUS];

        const double equivalent_strain = inner_prod(predictive_stress_vector, r_strain_vector) / equivalent_stress;
        const double m = r_props[COHESION];
        const double strain_regularization = (equivalent_strain - yield_strain > 0.0) ? 1.0 - std::exp(-m * (equivalent_strain - yield_strain) / yield_strain) : 0.0;
        const double viscous_overstress = strain_regularization * viscous_eta * (std::exp(viscous_alpha * strain_rate_norm / yield_strain) - 1.0) * std::exp(-viscous_beta * mViscousTime / time_delay);

        double F = equivalent_stress - threshold - viscous_overstress;

        if (F > 0.0) {
            // Integrate Stress plasticity
            const double characteristic_length = AdvCLUtils::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());

            IndexType iteration = 0;
            IndexType max_iter = r_props.Has(MAX_NUMBER_NL_CL_ITERATIONS) ? r_props[MAX_NUMBER_NL_CL_ITERATIONS] : 100;
            bool is_converged = false;

            // Initialize Plastic Parameters
            double plastic_denominator = 0.0;
            array_1d<double, VoigtSize> dF_dS; // DF/DS
            array_1d<double, VoigtSize> dG_dS; // DG/DS
            array_1d<double, VoigtSize> plastic_strain_increment;
            dF_dS.clear();
            dG_dS.clear();
            plastic_strain_increment.clear();
            double plastic_consistency_factor_increment; // lambda dot

            // Compute the plastic parameters
            F = CLIntegrator::CalculatePlasticParameters(
                    predictive_stress_vector, r_strain_vector, equivalent_stress,
                    threshold, plastic_denominator, dF_dS, dG_dS,
                    plastic_dissipation, plastic_strain_increment,
                    r_constitutive_matrix, rValues, characteristic_length,
                    plastic_strain) -
                viscous_overstress;

            // Backward Euler
            while (!is_converged && iteration <= max_iter) {
                plastic_consistency_factor_increment = (F * plastic_denominator > 0.0) ? F * plastic_denominator : 0.0;
                noalias(plastic_strain_increment) = plastic_consistency_factor_increment * dG_dS;
                noalias(plastic_strain) += plastic_strain_increment;
    
                noalias(predictive_stress_vector) -= prod(r_constitutive_matrix, plastic_strain_increment);

                F = CLIntegrator::CalculatePlasticParameters(
                    predictive_stress_vector, r_strain_vector, equivalent_stress,
                    threshold, plastic_denominator, dF_dS, dG_dS,
                    plastic_dissipation, plastic_strain_increment,
                    r_constitutive_matrix, rValues, characteristic_length,
                    plastic_strain) -
                viscous_overstress;

                if (F <= std::abs(return_mapping_tol * threshold)) { // Has converged
                    is_converged = true;
                } else {
                    iteration++;
                }
            }
            KRATOS_WARNING_IF("Backward Euler ViscoPlasticity: ", iteration > max_iter) << "Maximum number of iterations in plasticity loop, F = " << F << std::endl;

            if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                CalculateTangentTensor(rValues, plastic_strain);
            }
        }
        noalias(rValues.GetStressVector()) = predictive_stress_vector;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainIsotropicCornejoViscoPlasticity<TConstLawIntegratorType>::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();

    // We get the constitutive tensor
    Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();

    if (r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        BaseType::CalculateCauchyGreenStrain(rValues, r_strain_vector);
    }
    this->template AddInitialStrainVectorContribution<Vector>(r_strain_vector);

    if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        CalculateElasticMatrix(r_constitutive_matrix, rValues);

        const auto &r_props = rValues.GetMaterialProperties();

        // Retrieve int vars as reference now
        double& threshold = GetThreshold();
        double& plastic_dissipation = GetPlasticDissipation();
        Vector& plastic_strain = GetPlasticStrain();
    
        const double time_regularization_factor = r_props.Has(TIME_REGULARIZATION) ? r_props[TIME_REGULARIZATION] : time_regularization;
        const double strain_rate_norm = time_regularization_factor * mStrainRateHistory[0] + (1.0 - time_regularization_factor) * mStrainRateHistory[1];

        BoundedArrayType predictive_stress_vector;
        noalias(predictive_stress_vector) = prod(r_constitutive_matrix, r_strain_vector - plastic_strain);
        this->template AddInitialStressVectorContribution<BoundedArrayType>(predictive_stress_vector);

        double equivalent_stress;
        CLIntegrator::YieldSurfaceType::CalculateEquivalentStress(predictive_stress_vector, r_strain_vector, equivalent_stress, rValues);

        const double viscous_eta = r_props[VISCOUS_PARAMETER];
        const double viscous_alpha = r_props.Has(VISCOUS_ALPHA) ? r_props[VISCOUS_ALPHA] : 1.0;
        const double viscous_beta = r_props.Has(VISCOUS_BETA) ? r_props[VISCOUS_BETA] : 1.0;
        const double time_delay = r_props[DELAY_TIME];
        const double yield_strain = r_props.Has(YIELD_STRESS) ? r_props[YIELD_STRESS] / r_props[YOUNG_MODULUS] : r_props[YIELD_STRESS_TENSION] / r_props[YOUNG_MODULUS];

        const double equivalent_strain = inner_prod(predictive_stress_vector, r_strain_vector) / equivalent_stress;
        const double m = r_props[COHESION];
        const double strain_regularization = (equivalent_strain-yield_strain > 0.0) ? 1.0 - std::exp(-m * (equivalent_strain-yield_strain) / yield_strain) : 0.0;
        const double viscous_overstress = strain_regularization * viscous_eta * (std::exp(viscous_alpha * strain_rate_norm / yield_strain) - 1.0) * std::exp(-viscous_beta * mViscousTime / time_delay);

        double F = equivalent_stress - threshold - viscous_overstress;

        if (F > 0.0) {
            mViscousTime += rValues.GetProcessInfo()[DELTA_TIME];
            // Integrate Stress plasticity
            const double characteristic_length = AdvCLUtils::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());

            IndexType iteration = 0;
            IndexType max_iter = r_props.Has(MAX_NUMBER_NL_CL_ITERATIONS) ? r_props[MAX_NUMBER_NL_CL_ITERATIONS] : 100;
            bool is_converged = false;

            // Initialize Plastic Parameters
            double plastic_denominator = 0.0;
            array_1d<double, VoigtSize> dF_dS; // DF/DS
            array_1d<double, VoigtSize> dG_dS; // DG/DS
            array_1d<double, VoigtSize> plastic_strain_increment;
            dF_dS.clear();
            dG_dS.clear();
            plastic_strain_increment.clear();
            double plastic_consistency_factor_increment; // lambda dot

            // Compute the plastic parameters
            F = CLIntegrator::CalculatePlasticParameters(
                    predictive_stress_vector, r_strain_vector, equivalent_stress,
                    threshold, plastic_denominator, dF_dS, dG_dS,
                    plastic_dissipation, plastic_strain_increment,
                    r_constitutive_matrix, rValues, characteristic_length,
                    plastic_strain) -
                viscous_overstress;

            // Backward Euler
            while (!is_converged && iteration <= max_iter) {
                plastic_consistency_factor_increment = (F * plastic_denominator > 0.0) ? F * plastic_denominator : 0.0;
                noalias(plastic_strain_increment) = plastic_consistency_factor_increment * dG_dS;
                noalias(plastic_strain) += plastic_strain_increment;
    
                noalias(predictive_stress_vector) -= prod(r_constitutive_matrix, plastic_strain_increment);

                F = CLIntegrator::CalculatePlasticParameters(
                    predictive_stress_vector, r_strain_vector, equivalent_stress,
                    threshold, plastic_denominator, dF_dS, dG_dS,
                    plastic_dissipation, plastic_strain_increment,
                    r_constitutive_matrix, rValues, characteristic_length,
                    plastic_strain) -
                viscous_overstress;

                if (F <= std::abs(return_mapping_tol * threshold)) { // Has converged
                    is_converged = true;
                } else {
                    iteration++;
                }
            }
            KRATOS_WARNING_IF("Backward Euler ViscoPlasticity: ", iteration > max_iter) << "Maximum number of iterations in plasticity loop, F = " << F << std::endl;
        } else {
            // If we enter elasticity, the viscous time is reset.
            mViscousTime = 0.0;
        }

        // We update the strain rate history
        mStrainRateHistory[1] = mStrainRateHistory[0];
        mStrainRateHistory[0] = inner_prod(predictive_stress_vector, r_strain_vector - mPreviousStrain) / (equivalent_stress * rValues.GetProcessInfo()[DELTA_TIME]);
        noalias(mPreviousStrain) = r_strain_vector;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
int GenericSmallStrainIsotropicCornejoViscoPlasticity<TConstLawIntegratorType>::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    const int check_base = BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);

    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(VISCOUS_PARAMETER)) << "The VISCOUS_PARAMETER variable is not defined in the properties..." << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(DELAY_TIME)) << "The DELAY_TIME variable is not defined in the properties..." << std::endl;

    CLIntegrator::Check(rMaterialProperties);

    return check_base;
}

/***********************************************************************************/
/***********************************************************************************/

template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>>;

// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<3>>>>;
// template class GenericSmallStrainIsotropicCornejoViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<3>>>>;

} // namespace Kratos
