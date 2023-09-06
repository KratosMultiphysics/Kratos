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
//  Main authors:    Sergio Jimenez
//                   Alejandro Cornejo
//  Collaborator:    Lucia Barbu
//

// System includes

// External includes

// Project includes
#include "utilities/math_utils.h"
#include "constitutive_laws_application_variables.h"
#include "uncoupled_plastic_damage_law.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"
#include "custom_constitutive/constitutive_laws_integrators/high_cycle_fatigue_law_integrator.h"
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_plasticity_with_fatigue.h"
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_damage.h"
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_plasticity.h"

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
template<class TConstLawIntegratorType>
ConstitutiveLaw::Pointer UncoupledPlasticDamageLaw<TConstLawIntegratorType>::Create(Kratos::Parameters NewParameters) const
{
    return Kratos::make_shared<UncoupledPlasticDamageLaw>(*this);

}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UncoupledPlasticDamageLaw<TConstLawIntegratorType>::InitializeMaterialResponsePK1(
    ConstitutiveLaw::Parameters& rValues
    )
{
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UncoupledPlasticDamageLaw<TConstLawIntegratorType>::InitializeMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rValues
    )
{
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UncoupledPlasticDamageLaw<TConstLawIntegratorType>::InitializeMaterialResponseKirchhoff(
    ConstitutiveLaw::Parameters& rValues
    )
{
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UncoupledPlasticDamageLaw<TConstLawIntegratorType>::InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    auto& r_material_properties = rValues.GetMaterialProperties();
    const auto it_cl_begin = r_material_properties.GetSubProperties().begin();
    const auto& r_props_isotropic_damage_cl = *(it_cl_begin);
    const auto& r_props_plasticity_cl = *(it_cl_begin + 1);
    ConstitutiveLaw::Parameters values_plastic_damage  = rValues;
    ConstitutiveLaw::Parameters values_damage_component  = rValues;
    ConstitutiveLaw::Parameters values_plasticity_component = rValues;
    values_damage_component.SetMaterialProperties(r_props_isotropic_damage_cl);
    values_plasticity_component.SetMaterialProperties(r_props_plasticity_cl);
    //Checking which material has the fatigue properties
    if (r_props_isotropic_damage_cl.Has(VOLUMETRIC_PART)) {
        values_plastic_damage.SetMaterialProperties(r_props_isotropic_damage_cl);
    } else if (r_props_plasticity_cl.Has(VOLUMETRIC_PART)) {
        values_plastic_damage.SetMaterialProperties(r_props_plasticity_cl);
    } else {
        KRATOS_ERROR << "Volumetric participation not defined" << std::endl;
    }
    const SizeType volumetric_participation_size = values_plastic_damage.GetMaterialProperties()[VOLUMETRIC_PART].size();
    double volumetric_participation;
    const double damage = mDamage;
    const double plastic_dissipation = mPlasticDissipation;
    double reference_damage = mReferenceDamage;
    double reference_plastic_dissipation = mReferencePlasticDissipation;
    double reference_volumetric_participation = mIsotrpicDamageVolumetricParticipation;
    double reference_equivalent_plastic_strain, equivalent_plastic_strain;
    const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(values_plasticity_component.GetElementGeometry());
    GenericConstitutiveLawIntegratorPlasticityWithFatigue<TConstLawIntegratorType::YieldSurfaceType>::EquivalencyPlasticDissipationUniaxialPlasticStrain(reference_equivalent_plastic_strain,
                                                                                                                                                        values_plasticity_component, reference_plastic_dissipation,
                                                                                                                                                        1.0, characteristic_length);
    GenericConstitutiveLawIntegratorPlasticityWithFatigue<TConstLawIntegratorType::YieldSurfaceType>::EquivalencyPlasticDissipationUniaxialPlasticStrain(equivalent_plastic_strain,
                                                                                                                                                        values_plasticity_component, plastic_dissipation,
                                                                                                                                                        1.0, characteristic_length);

    if (volumetric_participation_size == 1) {
        volumetric_participation = values_plastic_damage.GetMaterialProperties()[VOLUMETRIC_PART][0];
    } else if (volumetric_participation_size == 3) { //Linear or exponential transition
        const int volumetric_participation_transition_type = values_plastic_damage.GetMaterialProperties()[VOLUMETRIC_PART][0];
        const double initial_volumetric_participation = values_plastic_damage.GetMaterialProperties()[VOLUMETRIC_PART][1];
        const double final_volumetric_participation = values_plastic_damage.GetMaterialProperties()[VOLUMETRIC_PART][2];
        if (volumetric_participation_transition_type == 0) { //Linear case
            volumetric_participation = initial_volumetric_participation * (1.0 - damage) + final_volumetric_participation * damage;
        } else { //Exponential case
            volumetric_participation = initial_volumetric_participation * std::exp(damage * std::log(final_volumetric_participation / initial_volumetric_participation));
        }
    } else if (volumetric_participation_size == 4) { //Potential or inverse potential transition
        const int volumetric_participation_transition_type = values_plastic_damage.GetMaterialProperties()[VOLUMETRIC_PART][0];
        const double initial_volumetric_participation = values_plastic_damage.GetMaterialProperties()[VOLUMETRIC_PART][1];
        const double final_volumetric_participation = values_plastic_damage.GetMaterialProperties()[VOLUMETRIC_PART][2];
        const double index_volumetric_participation = values_plastic_damage.GetMaterialProperties()[VOLUMETRIC_PART][3];
        if (volumetric_participation_transition_type == 0) { //Potential case
            volumetric_participation = (final_volumetric_participation - initial_volumetric_participation) * std::pow(damage, index_volumetric_participation) + initial_volumetric_participation;
        } else { //Inverse potential case
            volumetric_participation = (final_volumetric_participation - initial_volumetric_participation) * std::pow((damage - 1.0), index_volumetric_participation) + final_volumetric_participation;
        }
    } else {
        KRATOS_ERROR << "Wrong size VOLUMETRIC_PART variable" << std::endl;
    }
    //Limits to not recover aparent stiffness nor reduce aparent plastic strain
    const double min_volumetric_participation = (damage > 0.0) ? reference_damage / damage * reference_volumetric_participation : 0.0;
    const double max_volumetric_participation = (plastic_dissipation > 0.0) ? ((1.0 - reference_volumetric_participation) * reference_equivalent_plastic_strain - (1.0 - reference_volumetric_participation * reference_damage) * equivalent_plastic_strain)
                                                / (damage * (1.0 - reference_volumetric_participation) * reference_equivalent_plastic_strain - (1.0 - reference_volumetric_participation * reference_damage) * equivalent_plastic_strain) : 1.0;

    // volumetric_participation = (min_volumetric_participation + max_volumetric_participation) / 2.0;
    volumetric_participation = (volumetric_participation < min_volumetric_participation) ? min_volumetric_participation : volumetric_participation;
    volumetric_participation = (volumetric_participation > max_volumetric_participation) ? max_volumetric_participation : volumetric_participation;
    mMinVolumetricParticipation = min_volumetric_participation;
    mMaxVolumetricParticipation = max_volumetric_participation;
    mIsotrpicDamageVolumetricParticipation = volumetric_participation;
    mReferenceDamage = damage;
    mReferencePlasticDissipation = plastic_dissipation;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UncoupledPlasticDamageLaw<TConstLawIntegratorType>::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UncoupledPlasticDamageLaw<TConstLawIntegratorType>::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UncoupledPlasticDamageLaw<TConstLawIntegratorType>::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UncoupledPlasticDamageLaw<TConstLawIntegratorType>::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // Some auxiliar values
    const SizeType dimension = WorkingSpaceDimension();
    const SizeType voigt_size = GetStrainSize();

    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();

    // Previous flags saved
    const bool flag_strain = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN );
    const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
    const bool flag_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS );

    const Properties& r_material_properties  = rValues.GetMaterialProperties();

    // The deformation gradient
    if (rValues.IsSetDeterminantF()) {
        const double determinant_f = rValues.GetDeterminantF();
        KRATOS_ERROR_IF(determinant_f < 0.0) << "Deformation gradient determinant (detF) < 0.0 : " << determinant_f << std::endl;
    }
    // In case the element has not computed the Strain
    if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        Vector& r_strain_vector = rValues.GetStrainVector();

        Matrix F_deformation_gradient(dimension, dimension);
        this->CalculateValue(rValues, DEFORMATION_GRADIENT, F_deformation_gradient);
        const Matrix B_matrix = prod(F_deformation_gradient, trans(F_deformation_gradient));
        // Doing resize in case is needed
        if (r_strain_vector.size() != voigt_size)
            r_strain_vector.resize(voigt_size);

         // Identity matrix
        Matrix identity_matrix(dimension, dimension);
        for (IndexType i = 0; i < dimension; ++i) {
            for (IndexType j = 0; j < dimension; ++j) {
                if (i == j) identity_matrix(i, j) = 1.0;
                else identity_matrix(i, j) = 0.0;
            }
        }

        // Calculating the inverse of the left Cauchy tensor
        Matrix inverse_B_tensor(dimension, dimension);
        double aux_det_b = 0;
        MathUtils<double>::InvertMatrix(B_matrix, inverse_B_tensor, aux_det_b);

        // Calculate E matrix
        const Matrix E_matrix = 0.5 * (identity_matrix - inverse_B_tensor);
        // Almansi Strain Calculation
        r_strain_vector = MathUtils<double>::StrainTensorToVector(E_matrix, voigt_size);
    }

    if (r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        // Set new flags
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false); //Tangent tensor is not computed at the component level (only at the composite level)
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        // Total strain vector hich is equal to component CLs strains
        Vector& r_strain_vector = rValues.GetStrainVector();

        // This method integrates the stress according to each simple material CL
        Vector isotropic_damage_stress_vector, plasticity_stress_vector;
        this->IntegrateStressesOfIsotropicDamageAndPlasticityModels(rValues, r_strain_vector, r_strain_vector, isotropic_damage_stress_vector, plasticity_stress_vector);

        Vector& r_integrated_stress_vector = rValues.GetStressVector();
        noalias(r_integrated_stress_vector) = mIsotrpicDamageVolumetricParticipation * isotropic_damage_stress_vector
                                     + (1.0 - mIsotrpicDamageVolumetricParticipation) * plasticity_stress_vector;
        // Previous flags restored
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);

        if (r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            this->CalculateTangentTensor(rValues);
        }
    }

} // End CalculateMaterialResponseCauchy

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UncoupledPlasticDamageLaw<TConstLawIntegratorType>::IntegrateStressesOfIsotropicDamageAndPlasticityModels(
    ConstitutiveLaw::Parameters& rValues,
    Vector rIsotropicDamgeStrainVector,
    Vector rPlasticityStrainVector,
    Vector& rIsotropicDamageStressVector,
    Vector& rPlasticityStressVector
)
{
    auto& r_material_properties = rValues.GetMaterialProperties();
    const auto it_cl_begin = r_material_properties.GetSubProperties().begin();
    const auto& r_props_isotropic_damage_cl = *(it_cl_begin);
    const auto& r_props_plasticity_cl = *(it_cl_begin + 1);

    ConstitutiveLaw::Parameters values_damage_component  = rValues;
    ConstitutiveLaw::Parameters values_plasticity_component = rValues;

    values_damage_component.SetStrainVector(rIsotropicDamgeStrainVector);
    values_plasticity_component.SetStrainVector(rPlasticityStrainVector);


    // Integrate Stress of the damage component
    values_damage_component.SetMaterialProperties(r_props_isotropic_damage_cl);
    this->CalculateMaterialResponseIsotropicDamageModelModel(values_damage_component);
    rIsotropicDamageStressVector = values_damage_component.GetStressVector();

    // Integrate Stress of the plasticity component
    values_plasticity_component.SetMaterialProperties(r_props_plasticity_cl);
    this->CalculateMaterialResponsePlasticityModel(values_plasticity_component);
    rPlasticityStressVector = values_plasticity_component.GetStressVector();
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UncoupledPlasticDamageLaw<TConstLawIntegratorType>::CalculateMaterialResponseIsotropicDamageModelModel(
        ConstitutiveLaw::Parameters& values_damage_component
)
{
    Vector& integrated_stress_vector = values_damage_component.GetStressVector();
    array_1d<double, VoigtSize> auxiliar_integrated_stress_vector = integrated_stress_vector;
    Matrix& r_tangent_tensor = values_damage_component.GetConstitutiveMatrix(); // todo modify after integration
    const Flags& r_constitutive_law_options = values_damage_component.GetOptions();
    Matrix& r_constitutive_matrix = values_damage_component.GetConstitutiveMatrix();
    mpIsotropicDamageConstitutiveLaw->CalculateValue(values_damage_component, CONSTITUTIVE_MATRIX, r_constitutive_matrix);

    Vector& r_strain_vector = values_damage_component.GetStrainVector();
    this->template AddInitialStrainVectorContribution<Vector>(r_strain_vector);

    // Converged values
    double isotropic_damage_threshold = this->GetIsotropicDamageThreshold();
    double damage = this->GetDamage();

    // S0 = C:(E-E0) + S0
    array_1d<double, VoigtSize> predictive_stress_vector = prod(r_constitutive_matrix, r_strain_vector);
    this->template AddInitialStressVectorContribution<array_1d<double, VoigtSize>>(predictive_stress_vector);

    // Initialize Plastic Parameters
    double uniaxial_stress;
    TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(predictive_stress_vector, r_strain_vector, uniaxial_stress, values_damage_component);
    const double F = uniaxial_stress - isotropic_damage_threshold;
    const double sign_factor = HighCycleFatigueLawIntegrator<6>::CalculateTensionCompressionFactorComponentVM(predictive_stress_vector);

    if (F <= threshold_tolerance || sign_factor == -1.0) { // Elastic case
        noalias(auxiliar_integrated_stress_vector) = (1.0 - damage) * predictive_stress_vector;
		noalias(integrated_stress_vector) = auxiliar_integrated_stress_vector;

        if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            noalias(r_tangent_tensor) = (1.0 - damage) * r_constitutive_matrix;
        }
    } else { // Damage case
        // const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLength(values_damage_component.GetElementGeometry());
        const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(values_damage_component.GetElementGeometry());
        // This routine updates the PredictiveStress to verify the yield surf
        TConstLawIntegratorType::IntegrateStressVector(
            predictive_stress_vector,
            uniaxial_stress, damage,
            isotropic_damage_threshold, values_damage_component,
            characteristic_length);

        // Updated Values
        noalias(auxiliar_integrated_stress_vector) = predictive_stress_vector;

        if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            this->CalculateTangentTensor(values_damage_component);
        }
        noalias(integrated_stress_vector) = auxiliar_integrated_stress_vector;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UncoupledPlasticDamageLaw<TConstLawIntegratorType>::CalculateMaterialResponsePlasticityModel(
        ConstitutiveLaw::Parameters& values_plasticity_component
)
{
    // Auxiliar values
    const Flags& r_constitutive_law_options = values_plasticity_component.GetOptions();
    // We get the strain vector
    Vector& r_strain_vector = values_plasticity_component.GetStrainVector();
    // We get the constitutive tensor
    Matrix& r_constitutive_matrix = values_plasticity_component.GetConstitutiveMatrix();

    // We check the current step and NL iteration
    const ProcessInfo& r_current_process_info = values_plasticity_component.GetProcessInfo();
    const bool first_computation = (r_current_process_info[NL_ITERATION_NUMBER] == 1 && r_current_process_info[STEP] == 1) ? true : false;

    //NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    if (first_computation) { // First computation always pure elastic for elemements not providing the strain
        if (r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN ) ) {
            BaseType::CalculateCauchyGreenStrain( values_plasticity_component, r_strain_vector);
        }
        this->template AddInitialStrainVectorContribution<Vector>(r_strain_vector);
        if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS) ||
            r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            Vector& r_stress_vector = values_plasticity_component.GetStressVector();
            if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                BaseType::CalculateElasticMatrix( r_constitutive_matrix, values_plasticity_component);
                noalias(r_stress_vector) = prod( r_constitutive_matrix, r_strain_vector);
                this->template AddInitialStressVectorContribution<Vector>(r_stress_vector);
            } else {
                BaseType::CalculatePK2Stress( r_strain_vector, r_stress_vector, values_plasticity_component);
                this->template AddInitialStressVectorContribution<Vector>(r_stress_vector);
            }
        }
    } else { // We check for plasticity
        // Integrate Stress plasticity
        Vector& r_integrated_stress_vector = values_plasticity_component.GetStressVector();
        // const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLength(values_plasticity_component.GetElementGeometry());
        const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(values_plasticity_component.GetElementGeometry());

        if (r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
            BaseType::CalculateCauchyGreenStrain( values_plasticity_component, r_strain_vector);
        }

        this->template AddInitialStrainVectorContribution<Vector>(r_strain_vector);

        // We compute the stress or the constitutive matrix
        if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS) ||
            r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {

            // We get some variables
            double plasticity_threshold = this->GetPlasticityThreshold();
            double plastic_dissipation = this->GetPlasticDissipation();
            Vector plastic_strain = this->GetPlasticStrain();

            array_1d<double, VoigtSize> predictive_stress_vector;
            if (r_constitutive_law_options.Is( ConstitutiveLaw::U_P_LAW)) {
                noalias(predictive_stress_vector) = values_plasticity_component.GetStressVector();
            } else {
                // S0 = Elastic stress with strain (E-Ep) + S0
                Vector aux_stress = ZeroVector(VoigtSize);
                BaseType::CalculatePK2Stress(r_strain_vector - plastic_strain, aux_stress, values_plasticity_component);
                this->template AddInitialStressVectorContribution<Vector>(aux_stress);
                noalias(predictive_stress_vector) = aux_stress;
            }

            // Initialize Plastic Parameters
            double uniaxial_stress = 0.0, plastic_denominator = 0.0;
            array_1d<double, VoigtSize> f_flux = ZeroVector(VoigtSize); // DF/DS
            array_1d<double, VoigtSize> g_flux = ZeroVector(VoigtSize); // DG/DS
            array_1d<double, VoigtSize> plastic_strain_increment = ZeroVector(VoigtSize);

            // Elastic Matrix
            BaseType::CalculateElasticMatrix(r_constitutive_matrix, values_plasticity_component);

            // Compute the plastic parameters
            const double F = GenericConstitutiveLawIntegratorPlasticityWithFatigue<TConstLawIntegratorType::YieldSurfaceType>::CalculatePlasticParameters(
                predictive_stress_vector, r_strain_vector, uniaxial_stress,
                plasticity_threshold, plastic_denominator, f_flux, g_flux,
                plastic_dissipation, plastic_strain_increment,
                r_constitutive_matrix, values_plasticity_component, characteristic_length,
                plastic_strain, 1.0);

            const double sign_factor = HighCycleFatigueLawIntegrator<6>::CalculateTensionCompressionFactorComponentVM(predictive_stress_vector);
            const double yield_criteria = uniaxial_stress - plasticity_threshold;

            if (yield_criteria <= std::abs(1.0e-4 * plasticity_threshold) || sign_factor == -1.0) { // Elastic case
                noalias(r_integrated_stress_vector) = predictive_stress_vector;
            } else { // Plastic case
                // While loop backward euler
                /* Inside "IntegrateStressVector" the predictive_stress_vector is updated to verify the yield criterion */
                GenericConstitutiveLawIntegratorPlasticityWithFatigue<TConstLawIntegratorType::YieldSurfaceType>::IntegrateStressVector(
                    predictive_stress_vector, r_strain_vector, uniaxial_stress,
                    plasticity_threshold, plastic_denominator, f_flux, g_flux,
                    plastic_dissipation, plastic_strain_increment,
                    r_constitutive_matrix, plastic_strain, values_plasticity_component,
                    characteristic_length, 1.0);
                noalias(r_integrated_stress_vector) = predictive_stress_vector;

                if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                    this->CalculateTangentTensor(values_plasticity_component); // this modifies the ConstitutiveMatrix
                } else {
                    BaseType::CalculateElasticMatrix( r_constitutive_matrix, values_plasticity_component);
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UncoupledPlasticDamageLaw<TConstLawIntegratorType>::FinalizeSolutionStep(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Deprecated
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UncoupledPlasticDamageLaw<TConstLawIntegratorType>::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    this->FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UncoupledPlasticDamageLaw<TConstLawIntegratorType>::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    this->FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UncoupledPlasticDamageLaw<TConstLawIntegratorType>::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    this->FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UncoupledPlasticDamageLaw<TConstLawIntegratorType>::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    const Vector& r_strain_vector = rValues.GetStrainVector();

    // Recalculation to obtain the serial_strain_matrix and store the value
    const SizeType voigt_size = GetStrainSize();

    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();

    // Previous flags saved
    const bool flag_strain = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    const bool flag_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);

    const Properties& r_material_properties = rValues.GetMaterialProperties();

    if (r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        // Set new flags
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        // Total strain vector
        Vector& r_strain_vector = rValues.GetStrainVector();

        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);

        // We call the FinalizeMaterialResponse of the HCF and ULCF CL
        auto& r_material_properties = rValues.GetMaterialProperties();
        const auto it_cl_begin = r_material_properties.GetSubProperties().begin();
        const auto& r_props_isotropic_damage_cl = *(it_cl_begin);
        const auto& r_props_plasticity_cl = *(it_cl_begin + 1);

        ConstitutiveLaw::Parameters values_damage_component  = rValues;
        ConstitutiveLaw::Parameters values_plasticity_component = rValues;

        values_damage_component.SetMaterialProperties(r_props_isotropic_damage_cl);
        values_plasticity_component.SetMaterialProperties(r_props_plasticity_cl);

        values_damage_component.SetStrainVector(r_strain_vector);
        values_plasticity_component.SetStrainVector(r_strain_vector);

        Vector isotropic_damage_stress_vector, plasticity_stress_vector;

        this->FinalizeMaterialResponseIsotropicDamageModel(isotropic_damage_stress_vector, values_damage_component);
        this->FinalizeMaterialResponsePlastitictyModel(plasticity_stress_vector, values_plasticity_component);

        Vector r_integrated_stress_vector;
        r_integrated_stress_vector = mIsotrpicDamageVolumetricParticipation * isotropic_damage_stress_vector
                                     + (1.0 - mIsotrpicDamageVolumetricParticipation) * plasticity_stress_vector;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UncoupledPlasticDamageLaw<TConstLawIntegratorType>::FinalizeMaterialResponseIsotropicDamageModel(
        Vector& rIsotropicDamageStressVector,
        ConstitutiveLaw::Parameters& values_damage_component
)
{
    // Integrate Stress Damage
    const Flags& r_constitutive_law_options = values_damage_component.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = values_damage_component.GetStrainVector();

    // Elastic Matrix
    Matrix& r_constitutive_matrix = values_damage_component.GetConstitutiveMatrix();
    mpIsotropicDamageConstitutiveLaw->CalculateValue(values_damage_component, CONSTITUTIVE_MATRIX, r_constitutive_matrix);

    if (r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        BaseType::CalculateCauchyGreenStrain( values_damage_component, r_strain_vector);
    }

    this->template AddInitialStrainVectorContribution<Vector>(r_strain_vector);

    // Converged values
    double isotropic_damage_threshold = this->GetIsotropicDamageThreshold();
    double damage = this->GetDamage();

    // S0 = C:(E-E0) + S0
    array_1d<double, VoigtSize> predictive_stress_vector = prod(r_constitutive_matrix, r_strain_vector);
    this->template AddInitialStressVectorContribution<array_1d<double, VoigtSize>>(predictive_stress_vector);

    // Initialize Plastic Parameters
    double uniaxial_stress;
    TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(predictive_stress_vector, r_strain_vector, uniaxial_stress, values_damage_component);
    // rIsotropicDamagePredictiveUniaxialStress = uniaxial_stress;

    const double F = uniaxial_stress - isotropic_damage_threshold;
    const double sign_factor = HighCycleFatigueLawIntegrator<6>::CalculateTensionCompressionFactorComponentVM(predictive_stress_vector);

    if (F >= threshold_tolerance && sign_factor == 1.0) { // Plastic case
        // const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLength(values_damage_component.GetElementGeometry());
        const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(values_damage_component.GetElementGeometry());
        // This routine updates the PredictiveStress to verify the yield surf
        TConstLawIntegratorType::IntegrateStressVector(
            predictive_stress_vector,
            uniaxial_stress, damage,
            isotropic_damage_threshold, values_damage_component,
            characteristic_length);
        mDamage = damage;
        mIsotropicDamageThreshold = uniaxial_stress;

    } else {
        predictive_stress_vector *= (1.0 - mDamage);
    }
    rIsotropicDamageStressVector = predictive_stress_vector;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UncoupledPlasticDamageLaw<TConstLawIntegratorType>::FinalizeMaterialResponsePlastitictyModel(
        Vector& rPlasticityStressVector,
        ConstitutiveLaw::Parameters& values_plasticity_component
)
{
    // Auxiliar values
    const Flags& r_constitutive_law_options = values_plasticity_component.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = values_plasticity_component.GetStrainVector();

    // We get the constitutive tensor
    Matrix& r_constitutive_matrix = values_plasticity_component.GetConstitutiveMatrix();

    // Integrate Stress plasticity
    // const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLength(values_plasticity_component.GetElementGeometry());
    const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(values_plasticity_component.GetElementGeometry());

    if (r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        BaseType::CalculateCauchyGreenStrain( values_plasticity_component, r_strain_vector);
    }

    this->template AddInitialStrainVectorContribution<Vector>(r_strain_vector);

    // We compute the stress
    // Elastic Matrix
    this->CalculateElasticMatrix(r_constitutive_matrix, values_plasticity_component);

    // We get some variables
    double plasticity_threshold = this->GetPlasticityThreshold();
    double plastic_dissipation = this->GetPlasticDissipation();
    Vector plastic_strain = this->GetPlasticStrain();

    array_1d<double, VoigtSize> predictive_stress_vector;
    if (r_constitutive_law_options.Is( ConstitutiveLaw::U_P_LAW)) {
        noalias(predictive_stress_vector) = values_plasticity_component.GetStressVector();
    } else {
        // Spred = r_constitutive_matrix:(E-Ep) + S0
        noalias(predictive_stress_vector) = prod(r_constitutive_matrix, r_strain_vector - plastic_strain);
        this->template AddInitialStressVectorContribution<array_1d<double, VoigtSize>>(predictive_stress_vector);
    }

    // Initialize Plastic Parameters
    double uniaxial_stress = 0.0, plastic_denominator = 0.0;
    array_1d<double, VoigtSize> f_flux = ZeroVector(VoigtSize); // DF/DS
    array_1d<double, VoigtSize> g_flux = ZeroVector(VoigtSize); // DG/DS
    array_1d<double, VoigtSize> plastic_strain_increment = ZeroVector(VoigtSize);

    const double F = GenericConstitutiveLawIntegratorPlasticityWithFatigue<TConstLawIntegratorType::YieldSurfaceType>::CalculatePlasticParameters(
        predictive_stress_vector, r_strain_vector, uniaxial_stress,
        plasticity_threshold, plastic_denominator, f_flux, g_flux,
        plastic_dissipation, plastic_strain_increment,
        r_constitutive_matrix, values_plasticity_component, characteristic_length,
        plastic_strain, 1.0);

    const double sign_factor = HighCycleFatigueLawIntegrator<6>::CalculateTensionCompressionFactorComponentVM(predictive_stress_vector);
    const double yield_criteria = uniaxial_stress - plasticity_threshold;

    if (yield_criteria > std::abs(1.0e-4 * plasticity_threshold) && sign_factor == 1.0) { // Plastic case
        // While loop backward euler
        /* Inside "IntegrateStressVector" the predictive_stress_vector is updated to verify the yield criterion */
        GenericConstitutiveLawIntegratorPlasticityWithFatigue<TConstLawIntegratorType::YieldSurfaceType>::IntegrateStressVector(
            predictive_stress_vector, r_strain_vector, uniaxial_stress,
            plasticity_threshold, plastic_denominator, f_flux, g_flux,
            plastic_dissipation, plastic_strain_increment,
            r_constitutive_matrix, plastic_strain, values_plasticity_component,
            characteristic_length, 1.0);
        BaseType::CalculateElasticMatrix( r_constitutive_matrix, values_plasticity_component);
    }
    mPlasticDissipation = plastic_dissipation;
    mPlasticStrain = plastic_strain;
    mPlasticityThreshold = plasticity_threshold;
    rPlasticityStressVector = predictive_stress_vector;
}
/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
bool UncoupledPlasticDamageLaw<TConstLawIntegratorType>::Has(const Variable<bool>& rThisVariable)
{
    if (rThisVariable == CYCLE_INDICATOR) {
        return true;
    }
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
bool UncoupledPlasticDamageLaw<TConstLawIntegratorType>::Has(const Variable<int>& rThisVariable)
{
    if (rThisVariable == NUMBER_OF_CYCLES) {
        return true;
    } else if (rThisVariable == LOCAL_NUMBER_OF_CYCLES) {
        return true;
    }
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
bool UncoupledPlasticDamageLaw<TConstLawIntegratorType>::Has(const Variable<double>& rThisVariable)
{
    if (rThisVariable == DAMAGE) {
        return true;
    } else if (rThisVariable == PLASTIC_DISSIPATION) {
        return true;
    } else if (rThisVariable == VOLUMETRIC_PARTICIPATION) {
        return true;
    } else if (rThisVariable == VISCOUS_PARAMETER) {
        return true;
    } else if (rThisVariable == DELAY_TIME) {
        return true;
    } else {
        return BaseType::Has(rThisVariable);
    }

    return false;

}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
bool UncoupledPlasticDamageLaw<TConstLawIntegratorType>::Has(const Variable<Vector>& rThisVariable)
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

template <class TConstLawIntegratorType>
bool UncoupledPlasticDamageLaw<TConstLawIntegratorType>::Has(const Variable<Matrix>& rThisVariable)
{
    return BaseType::Has(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UncoupledPlasticDamageLaw<TConstLawIntegratorType>::SetValue(
    const Variable<bool>& rThisVariable,
    const bool& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
}

/***********************************************************************************/
/***********************************************************************************/
template <class TConstLawIntegratorType>
void UncoupledPlasticDamageLaw<TConstLawIntegratorType>::SetValue(
    const Variable<int>& rThisVariable,
    const int& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UncoupledPlasticDamageLaw<TConstLawIntegratorType>::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rThisVariable == PLASTIC_DISSIPATION) {
        mPlasticDissipation = rValue;
    } else if (rThisVariable == DAMAGE) {
        mDamage = rValue;
    } else if (rThisVariable == VOLUMETRIC_PARTICIPATION) {
        mIsotrpicDamageVolumetricParticipation = rValue;
    } else if (rThisVariable == VISCOUS_PARAMETER) {
        mMinVolumetricParticipation = rValue;
    } else if (rThisVariable == DELAY_TIME) {
        mMaxVolumetricParticipation = rValue;
    } else {
        BaseType::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UncoupledPlasticDamageLaw<TConstLawIntegratorType>::SetValue(
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

template <class TConstLawIntegratorType>
bool& UncoupledPlasticDamageLaw<TConstLawIntegratorType>::GetValue(
    const Variable<bool>& rThisVariable,
    bool& rValue
    )
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
int& UncoupledPlasticDamageLaw<TConstLawIntegratorType>::GetValue(
    const Variable<int>& rThisVariable,
    int& rValue
    )
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
double& UncoupledPlasticDamageLaw<TConstLawIntegratorType>::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if (rThisVariable == DAMAGE) {
        rValue = mDamage;
    } else if (rThisVariable == PLASTIC_DISSIPATION) {
        rValue = mPlasticDissipation;
    } else if (rThisVariable == VOLUMETRIC_PARTICIPATION) {
        rValue = mIsotrpicDamageVolumetricParticipation;
    } else if (rThisVariable == VISCOUS_PARAMETER) {
        rValue = mMinVolumetricParticipation;
    } else if (rThisVariable == DELAY_TIME) {
        rValue = mMaxVolumetricParticipation;
    } else {
        return BaseType::GetValue(rThisVariable, rValue);
    }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
Vector& UncoupledPlasticDamageLaw<TConstLawIntegratorType>::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
        rValue = mPlasticStrain;
        return rValue;
    } else {
        return rValue;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
Matrix& UncoupledPlasticDamageLaw<TConstLawIntegratorType>::GetValue(
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

template <class TConstLawIntegratorType>
double& UncoupledPlasticDamageLaw<TConstLawIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue)
{
    auto& r_material_properties = rParameterValues.GetMaterialProperties();
    const auto it_cl_begin = r_material_properties.GetSubProperties().begin();
    const auto& r_props_isotropic_damage_cl = *(it_cl_begin);
    const auto& r_props_plasticity_cl = *(it_cl_begin + 1);
    ConstitutiveLaw::Parameters values_damage_component  = rParameterValues;
    ConstitutiveLaw::Parameters values_plasticity_component = rParameterValues;
    values_damage_component.SetMaterialProperties(r_props_isotropic_damage_cl);
    values_plasticity_component.SetMaterialProperties(r_props_plasticity_cl);

    if (rThisVariable == UNIAXIAL_STRESS) {
        // Calculate the stress vector
        CalculateMaterialResponseCauchy(rParameterValues);

        Vector& r_strain_vector = rParameterValues.GetStrainVector();
        Vector& r_integrated_stress_vector = rParameterValues.GetStressVector();
        TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress( r_integrated_stress_vector, r_strain_vector, rValue, rParameterValues);

        return rValue;

    } else if (rThisVariable == UNIAXIAL_STRESS_HCF) {
        this->CalculateMaterialResponseIsotropicDamageModelModel(values_damage_component);
        Vector& r_strain_vector = values_damage_component.GetStrainVector();
        Vector& r_integrated_stress_vector = values_damage_component.GetStressVector();
        double uniaxial_stress;
        TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress( r_integrated_stress_vector, r_strain_vector, uniaxial_stress, values_damage_component);
        const double sign_factor = HighCycleFatigueLawIntegrator<6>::CalculateTensionCompressionFactorComponentVM(r_integrated_stress_vector);
        rValue = sign_factor * uniaxial_stress;
        return rValue;
    } else if (rThisVariable == UNIAXIAL_STRESS_ULCF) {
        this->CalculateMaterialResponsePlasticityModel(values_plasticity_component);
        Vector& r_strain_vector = rParameterValues.GetStrainVector();
        Vector& r_integrated_stress_vector = values_plasticity_component.GetStressVector();
        double uniaxial_stress;
        TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress( r_integrated_stress_vector, r_strain_vector, uniaxial_stress, values_plasticity_component);
        const double sign_factor = HighCycleFatigueLawIntegrator<6>::CalculateTensionCompressionFactorComponentVM(r_integrated_stress_vector);
        rValue = sign_factor * uniaxial_stress;
        return rValue;
    } else {
        return this->GetValue(rThisVariable, rValue);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
Vector& UncoupledPlasticDamageLaw<TConstLawIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue)
{
    return this->GetValue(rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
Matrix& UncoupledPlasticDamageLaw<TConstLawIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    // We do some special operations for constitutive matrices
    if (rThisVariable == CONSTITUTIVE_MATRIX ||
        rThisVariable == CONSTITUTIVE_MATRIX_PK2 ||
        rThisVariable == CONSTITUTIVE_MATRIX_KIRCHHOFF) {
        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_strain = r_flags.Is( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN );
        const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
        const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS );

        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, false);

         // We compute the constitutive matrix
        if (rThisVariable == CONSTITUTIVE_MATRIX) {
            this->CalculateMaterialResponse(rParameterValues, this->GetStressMeasure());
        } else if (rThisVariable == CONSTITUTIVE_MATRIX_PK2) {
            this->CalculateMaterialResponsePK2(rParameterValues);
        } else if (rThisVariable == CONSTITUTIVE_MATRIX_KIRCHHOFF) {
            this->CalculateMaterialResponsePK2(rParameterValues);
        }

        noalias(rValue) = rParameterValues.GetConstitutiveMatrix();

        // Previous flags restored
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
    } else if (rThisVariable == DEFORMATION_GRADIENT) { // TODO: Make in the future modifications for take into account different layers combinations
        noalias(rValue) = rParameterValues.GetDeformationGradientF();
    } else if (rThisVariable == CAUCHY_STRESS_TENSOR) {
        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS);

        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        // We compute the stress
        this->CalculateMaterialResponseCauchy(rParameterValues);
        rValue = MathUtils<double>::StressVectorToTensor(rParameterValues.GetStressVector());

        // Previous flags restored
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );
        return rValue;
    } else {
        Matrix aux_value;
        Properties material_properties  = rParameterValues.GetMaterialProperties();
        Properties& r_prop = material_properties.GetSubProperties(0);

        rValue.clear();
        rParameterValues.SetMaterialProperties(r_prop);
        mpIsotropicDamageConstitutiveLaw->CalculateValue(rParameterValues, rThisVariable, aux_value);
        noalias(rValue) += mIsotrpicDamageVolumetricParticipation * aux_value;

        r_prop = material_properties.GetSubProperties(1);
        rParameterValues.SetMaterialProperties(r_prop);
        mpPlasticityConstitutiveLaw->CalculateValue(rParameterValues, rThisVariable, aux_value);
        noalias(rValue) += (1.0 - mIsotrpicDamageVolumetricParticipation) * aux_value;

        // Reset properties
        rParameterValues.SetMaterialProperties(material_properties);
    }
    return(rValue);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UncoupledPlasticDamageLaw<TConstLawIntegratorType>::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues)
{
    const auto it_cl_begin = rMaterialProperties.GetSubProperties().begin();
    const auto r_props_isotropic_damage_cl = *(it_cl_begin);
    const auto r_props_plasticity_cl  = *(it_cl_begin + 1);

    KRATOS_ERROR_IF_NOT(r_props_isotropic_damage_cl.Has(CONSTITUTIVE_LAW)) << "No constitutive law set" << std::endl;
    KRATOS_ERROR_IF_NOT(r_props_plasticity_cl.Has(CONSTITUTIVE_LAW))  << "No constitutive law set" << std::endl;

    mpIsotropicDamageConstitutiveLaw = r_props_isotropic_damage_cl[CONSTITUTIVE_LAW]->Clone();
    mpPlasticityConstitutiveLaw  = r_props_plasticity_cl[CONSTITUTIVE_LAW]->Clone();
    mpIsotropicDamageConstitutiveLaw->InitializeMaterial(r_props_isotropic_damage_cl, rElementGeometry, rShapeFunctionsValues);
    mpPlasticityConstitutiveLaw ->InitializeMaterial(r_props_plasticity_cl, rElementGeometry, rShapeFunctionsValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void UncoupledPlasticDamageLaw<TConstLawIntegratorType>::CalculateTangentTensor(ConstitutiveLaw::Parameters& rValues)
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
template class UncoupledPlasticDamageLaw<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>;

} // namespace Kratos