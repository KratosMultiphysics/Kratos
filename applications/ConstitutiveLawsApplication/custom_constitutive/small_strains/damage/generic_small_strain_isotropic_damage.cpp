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
//  Main authors:    Alejandro Cornejo & Lucia Barbu
//  Collaborator:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "custom_utilities/tangent_operator_calculator_utility.h"
#include "custom_utilities/automatic_differentiation_tangent_utilities.h"
#include "constitutive_laws_application_variables.h"
#include "generic_small_strain_isotropic_damage.h"
#include "custom_constitutive/auxiliary_files/cl_integrators/generic_cl_integrator_damage.h"

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

#include "custom_constitutive/thermal/auxiliary_files/thermal_yield_surfaces/thermal_von_mises_yield_surface.h"
#include "custom_constitutive/thermal/auxiliary_files/thermal_yield_surfaces/thermal_tresca_yield_surface.h"
#include "custom_constitutive/thermal/auxiliary_files/thermal_yield_surfaces/thermal_rankine_yield_surface.h"
#include "custom_constitutive/thermal/auxiliary_files/thermal_yield_surfaces/thermal_simo_ju_yield_surface.h"
#include "custom_constitutive/thermal/auxiliary_files/thermal_yield_surfaces/thermal_mohr_coulomb_yield_surface.h"
#include "custom_constitutive/thermal/auxiliary_files/thermal_yield_surfaces/thermal_modified_mohr_coulomb_yield_surface.h"
#include "custom_constitutive/thermal/auxiliary_files/thermal_yield_surfaces/thermal_drucker_prager_yield_surface.h"

namespace Kratos
{
template <class TConstLawIntegratorType>
void GenericSmallStrainIsotropicDamage<TConstLawIntegratorType>::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainIsotropicDamage<TConstLawIntegratorType>::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainIsotropicDamage<TConstLawIntegratorType>::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainIsotropicDamage<TConstLawIntegratorType>::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // Integrate Stress Damage
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // We get the strain vector
    auto& r_strain_vector = rValues.GetStrainVector();

    //NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    if (r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        BaseType::CalculateCauchyGreenStrain(rValues, r_strain_vector);
    }

    // We compute the stress
    if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        auto& r_integrated_stress_vector = rValues.GetStressVector();
        // Elastic Matrix
        auto& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);

        this->template AddInitialStrainVectorContribution<Vector>(r_strain_vector);

        // Converged values
        double threshold = this->GetThreshold();
        double damage = this->GetDamage();

        // S0 = C:(E-E0) + S0
        BoundedArrayType predictive_stress_vector = prod(r_constitutive_matrix, r_strain_vector);
        this->template AddInitialStressVectorContribution<BoundedArrayType>(predictive_stress_vector);

        double uniaxial_stress;
        TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(predictive_stress_vector, r_strain_vector, uniaxial_stress, rValues);

        const double F = uniaxial_stress - threshold;

        if (F <= threshold_tolerance) { // Elastic case
            noalias(r_integrated_stress_vector) = (1.0 - damage) * predictive_stress_vector;

            if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                r_constitutive_matrix *= (1.0 - damage);
            }
        } else { // Damage case
            const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());
            // This routine updates the PredictiveStress to verify the yield surf
            TConstLawIntegratorType::IntegrateStressVector(predictive_stress_vector, uniaxial_stress, damage, threshold, rValues, characteristic_length);

            // Updated Values
            noalias(r_integrated_stress_vector) = predictive_stress_vector;

            if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                this->CalculateTangentTensor(rValues);
            }
        }
    }
} // End CalculateMaterialResponseCauchy

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainIsotropicDamage<TConstLawIntegratorType>::CalculateTangentTensor(ConstitutiveLaw::Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();

    const bool consider_perturbation_threshold = r_material_properties.Has(CONSIDER_PERTURBATION_THRESHOLD) ? r_material_properties[CONSIDER_PERTURBATION_THRESHOLD] : true;
    const TangentOperatorEstimation tangent_operator_estimation = r_material_properties.Has(TANGENT_OPERATOR_ESTIMATION) ? static_cast<TangentOperatorEstimation>(r_material_properties[TANGENT_OPERATOR_ESTIMATION]) : TangentOperatorEstimation::SecondOrderPerturbation;

    if (tangent_operator_estimation == TangentOperatorEstimation::Analytic) {
        const SizeType softening_type = r_material_properties[SOFTENING_TYPE];
        if (softening_type == static_cast<SizeType>(SofteningType::Linear)) {
            AutomaticDifferentiationTangentUtilities<typename TConstLawIntegratorType::YieldSurfaceType, 0>::CalculateTangentTensorIsotropicDamage(rValues);
        } else if (softening_type == static_cast<SizeType>(SofteningType::Exponential)) {
            AutomaticDifferentiationTangentUtilities<typename TConstLawIntegratorType::YieldSurfaceType, 1>::CalculateTangentTensorIsotropicDamage(rValues);
        } else {
            KRATOS_ERROR << "The analytical tangent operator is not implemented for this SOFTENING_TYPE" << std::endl;
        }

    } else if (tangent_operator_estimation == TangentOperatorEstimation::FirstOrderPerturbation) {
        // Calculates the Tangent Constitutive Tensor by perturbation (first order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 1);
    } else if (tangent_operator_estimation == TangentOperatorEstimation::SecondOrderPerturbation) {
        // Calculates the Tangent Constitutive Tensor by perturbation (second order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 2);
    } else if (tangent_operator_estimation == TangentOperatorEstimation::Secant) {
        rValues.GetConstitutiveMatrix() *= (1.0 - mDamage);
    } else if (tangent_operator_estimation == TangentOperatorEstimation::InitialStiffness) {
        return;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainIsotropicDamage<TConstLawIntegratorType>::InitializeMaterial(
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
void GenericSmallStrainIsotropicDamage<TConstLawIntegratorType>::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainIsotropicDamage<TConstLawIntegratorType>::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainIsotropicDamage<TConstLawIntegratorType>::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainIsotropicDamage<TConstLawIntegratorType>::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // Integrate Stress Damage
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();

    //NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    if (r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        BaseType::CalculateCauchyGreenStrain( rValues, r_strain_vector);
    }

    // We compute the stress
    if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        // Elastic Matrix
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);

        this->template AddInitialStrainVectorContribution<Vector>(r_strain_vector);

        // S0 = C:(E-E0) + S0
        BoundedArrayType predictive_stress_vector = prod(r_constitutive_matrix, r_strain_vector);
        this->template AddInitialStressVectorContribution<BoundedArrayType>(predictive_stress_vector);

        // Initialize Plastic Parameters
        double uniaxial_stress;
        TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(predictive_stress_vector, r_strain_vector, uniaxial_stress, rValues);

        const double F = uniaxial_stress - mThreshold;

        if (F >= threshold_tolerance) { // Damage case
            const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());
            // This routine updates the PredictiveStress to verify the yield surf
            TConstLawIntegratorType::IntegrateStressVector(predictive_stress_vector, uniaxial_stress, mDamage, mThreshold, rValues, characteristic_length);
            mThreshold = uniaxial_stress;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
bool GenericSmallStrainIsotropicDamage<TConstLawIntegratorType>::Has(const Variable<double>& rThisVariable)
{
    if (rThisVariable == DAMAGE) {
        return true;
    } else if (rThisVariable == THRESHOLD) {
        return true;
    } else {
        return BaseType::Has(rThisVariable);
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
bool GenericSmallStrainIsotropicDamage<TConstLawIntegratorType>::Has(const Variable<Vector>& rThisVariable)
{
    if(rThisVariable == INTERNAL_VARIABLES){
        return true;
    }
    return BaseType::Has(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
bool GenericSmallStrainIsotropicDamage<TConstLawIntegratorType>::Has(const Variable<Matrix>& rThisVariable)
{
    return BaseType::Has(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainIsotropicDamage<TConstLawIntegratorType>::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (rThisVariable == DAMAGE) {
        mDamage = rValue;
    } else if (rThisVariable == THRESHOLD) {
        mThreshold = rValue;
    } else {
        return BaseType::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainIsotropicDamage<TConstLawIntegratorType>::SetValue(
    const Variable<Vector>& rThisVariable,
    const Vector& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (rThisVariable == INTERNAL_VARIABLES) {
        mDamage = rValue[0];
        mThreshold = rValue[1];
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
double& GenericSmallStrainIsotropicDamage<TConstLawIntegratorType>::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if (rThisVariable == DAMAGE) {
        rValue = mDamage;
    } else if (rThisVariable == THRESHOLD) {
        rValue = mThreshold;
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
Vector& GenericSmallStrainIsotropicDamage<TConstLawIntegratorType>::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    if(rThisVariable == INTERNAL_VARIABLES){
        rValue.resize(2);
        rValue[0] = mDamage;
        rValue[1] = mThreshold;
    }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
double& GenericSmallStrainIsotropicDamage<TConstLawIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if (rThisVariable == UNIAXIAL_STRESS) {
        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
        const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS );

        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, true );

        // Calculate the stress vector
        CalculateMaterialResponseCauchy(rParameterValues);
        const Vector& r_stress_vector = rParameterValues.GetStressVector();
        const Vector& r_strain_vector = rParameterValues.GetStrainVector();

        BoundedArrayType aux_stress_vector = r_stress_vector;
        TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(aux_stress_vector, r_strain_vector, rValue, rParameterValues);

        // Previous flags restored
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
        return rValue;
    } else {
        return BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
Vector& GenericSmallStrainIsotropicDamage<TConstLawIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    return BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
Matrix& GenericSmallStrainIsotropicDamage<TConstLawIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{

    return BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
int GenericSmallStrainIsotropicDamage<TConstLawIntegratorType>::Check(
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

template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<RankinePlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>>;

template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalVonMisesYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalTrescaYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalDruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalRankineYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalSimoJuYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>;

template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>>;
template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>;
template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<TrescaPlasticPotential<3>>>>;
template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<3>>>>;
template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<RankinePlasticPotential<3>>>>;
template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<3>>>>;
template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<3>>>>;

template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalVonMisesYieldSurface<VonMisesPlasticPotential<3>>>>;
template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>;
template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalTrescaYieldSurface<VonMisesPlasticPotential<3>>>>;
template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalDruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>>;
template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalRankineYieldSurface<VonMisesPlasticPotential<3>>>>;
template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalSimoJuYieldSurface<VonMisesPlasticPotential<3>>>>;
template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ThermalMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>;
} // namespace Kratos
