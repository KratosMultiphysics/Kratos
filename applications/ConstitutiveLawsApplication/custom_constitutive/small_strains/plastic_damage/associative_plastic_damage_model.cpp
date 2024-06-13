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
//                   Sergio Jimenez
//
//
// System includes
#include <iostream>
#include <set>

// External includes

// Project includes
#include "includes/checks.h"
#include "associative_plastic_damage_model.h"
#include "constitutive_laws_application_variables.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"

// Plasticity Integrator includes
#include "custom_constitutive/auxiliary_files/cl_integrators/generic_cl_integrator_plasticity.h"

// Yield surfaces
#include "custom_constitutive/auxiliary_files/yield_surfaces/generic_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/von_mises_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/drucker_prager_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/modified_mohr_coulomb_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/rankine_yield_surface.h"

// Plastic potentials
#include "custom_constitutive/auxiliary_files/plastic_potentials/generic_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/von_mises_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/drucker_prager_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/modified_mohr_coulomb_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/rankine_plastic_potential.h"

namespace Kratos
{

/********************************CLONE**********************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
ConstitutiveLaw::Pointer AssociativePlasticDamageModel<TYieldSurfaceType>::Clone() const
{
    return Kratos::make_shared<AssociativePlasticDamageModel>(*this);
}


/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void AssociativePlasticDamageModel<TYieldSurfaceType>::CalculateMaterialResponseCauchy(
    ConstitutiveLaw::Parameters& rValues
    )
{
    // Auxiliary values
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();

    Vector& r_integrated_stress_vector = rValues.GetStressVector();
    const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());

    if (r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        BaseType::CalculateCauchyGreenStrain(rValues, r_strain_vector);
    }

    // We compute the stress or the constitutive matrix
    if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS) ||
        r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {

        PlasticDamageParameters plastic_damage_parameters = PlasticDamageParameters();
        const auto &r_mat_props = rValues.GetMaterialProperties();
        InitializePlasticDamageParameters(r_strain_vector, r_mat_props, characteristic_length, plastic_damage_parameters);

        CheckMinimumFractureEnergy(rValues, plastic_damage_parameters);

        const bool crack_reclosure = (r_mat_props.Has(CRACK_RECLOSING)) ? r_mat_props[CRACK_RECLOSING] : false;
        if (crack_reclosure) {
            // here we compute a predictor of the constitutive matrix
            BoundedVectorType predictor(VoigtSize);
            CalculateConstitutiveMatrix(rValues, plastic_damage_parameters);
            noalias(predictor) = prod(plastic_damage_parameters.ConstitutiveMatrix, plastic_damage_parameters.StrainVector - plastic_damage_parameters.PlasticStrain);
            double tension_parameter, compression_parameter;
            GenericConstitutiveLawIntegratorPlasticity<TYieldSurfaceType>::CalculateIndicatorsFactors(predictor, tension_parameter,compression_parameter);
            double det = 0.0;
            BoundedMatrixType temp;
            noalias(temp) = tension_parameter * plastic_damage_parameters.ComplianceMatrix + compression_parameter * plastic_damage_parameters.ComplianceMatrixCompression;
            MathUtils<double>::InvertMatrix(temp, plastic_damage_parameters.ConstitutiveMatrix, det);
        } else {
            CalculateConstitutiveMatrix(rValues, plastic_damage_parameters);
        }

        noalias(rValues.GetConstitutiveMatrix()) = plastic_damage_parameters.ConstitutiveMatrix;

        noalias(plastic_damage_parameters.StressVector) = prod(plastic_damage_parameters.ConstitutiveMatrix, r_strain_vector - plastic_damage_parameters.PlasticStrain);

        TYieldSurfaceType::CalculateEquivalentStress(plastic_damage_parameters.StressVector, plastic_damage_parameters.StrainVector, plastic_damage_parameters.UniaxialStress, rValues);

        plastic_damage_parameters.NonLinearIndicator = plastic_damage_parameters.UniaxialStress - mThreshold;

        if (plastic_damage_parameters.NonLinearIndicator <= std::abs(tolerance * mThreshold)) {
            noalias(r_integrated_stress_vector) = plastic_damage_parameters.StressVector;
        } else {
            IntegrateStressPlasticDamageMechanics(rValues, plastic_damage_parameters);
            noalias(r_integrated_stress_vector) = plastic_damage_parameters.StressVector;

            if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                CalculateTangentTensor(rValues, plastic_damage_parameters); // this modifies the ConstitutiveMatrix
            }
        }
    } // compute stress or constitutive matrix
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void AssociativePlasticDamageModel<TYieldSurfaceType>::FinalizeMaterialResponseCauchy(
    ConstitutiveLaw::Parameters& rValues
    )
{
    // Auxiliary values
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();

    const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::
        CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());

    if (r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        BaseType::CalculateCauchyGreenStrain(rValues, r_strain_vector);
    }

    PlasticDamageParameters plastic_damage_parameters = PlasticDamageParameters();
    const auto &r_mat_props = rValues.GetMaterialProperties();
    InitializePlasticDamageParameters(r_strain_vector, r_mat_props, characteristic_length, plastic_damage_parameters);

    CheckMinimumFractureEnergy(rValues, plastic_damage_parameters);

        const bool crack_reclosure = (r_mat_props.Has(CRACK_RECLOSING)) ? r_mat_props[CRACK_RECLOSING] : false;
        if (crack_reclosure) {
            // here we compute a predictor of the constitutive matrix
            BoundedVectorType predictor(VoigtSize);
            CalculateConstitutiveMatrix(rValues, plastic_damage_parameters);
            noalias(predictor) = prod(plastic_damage_parameters.ConstitutiveMatrix, plastic_damage_parameters.StrainVector - plastic_damage_parameters.PlasticStrain);
            double tension_parameter, compression_parameter;
            GenericConstitutiveLawIntegratorPlasticity<TYieldSurfaceType>::CalculateIndicatorsFactors(predictor, tension_parameter,compression_parameter);
            double det = 0.0;
            BoundedMatrixType temp;
            noalias(temp) = tension_parameter * plastic_damage_parameters.ComplianceMatrix + compression_parameter * plastic_damage_parameters.ComplianceMatrixCompression;
            MathUtils<double>::InvertMatrix(temp, plastic_damage_parameters.ConstitutiveMatrix, det);
        } else {
            CalculateConstitutiveMatrix(rValues, plastic_damage_parameters);
        }
    noalias(plastic_damage_parameters.StressVector) = prod(plastic_damage_parameters.ConstitutiveMatrix,
        r_strain_vector - plastic_damage_parameters.PlasticStrain);

    TYieldSurfaceType::CalculateEquivalentStress(plastic_damage_parameters.StressVector,
        plastic_damage_parameters.StrainVector, plastic_damage_parameters.UniaxialStress, rValues);

    plastic_damage_parameters.NonLinearIndicator = plastic_damage_parameters.UniaxialStress - mThreshold;

    if (plastic_damage_parameters.NonLinearIndicator > std::abs(tolerance * mThreshold)) {
        IntegrateStressPlasticDamageMechanics(rValues, plastic_damage_parameters);
        UpdateInternalVariables(plastic_damage_parameters);
    }

}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
double AssociativePlasticDamageModel<TYieldSurfaceType>::CalculateVolumetricFractureEnergy( // g_F
    const Properties& rMaterialProperties,
    PlasticDamageParameters &rPDParameters
    )
{
    double tension_parameter, compression_parameter;
    GenericConstitutiveLawIntegratorPlasticity<TYieldSurfaceType>::CalculateIndicatorsFactors(rPDParameters.StressVector, tension_parameter,compression_parameter);

    const bool has_symmetric_yield_stress = rMaterialProperties.Has(YIELD_STRESS);
    const double yield_compression = has_symmetric_yield_stress ? rMaterialProperties[YIELD_STRESS] : rMaterialProperties[YIELD_STRESS_COMPRESSION];
    const double yield_tension = has_symmetric_yield_stress ? rMaterialProperties[YIELD_STRESS] : rMaterialProperties[YIELD_STRESS_TENSION];
    const double n = yield_compression / yield_tension;
    const double fracture_energy_tension = rMaterialProperties[FRACTURE_ENERGY]; // Frac energy in tension
    const double fracture_energy_compression = rMaterialProperties.Has(FRACTURE_ENERGY_COMPRESSION) ? rMaterialProperties[FRACTURE_ENERGY_COMPRESSION] : rMaterialProperties[FRACTURE_ENERGY] * std::pow(n, 2);

    const double characteristic_fracture_energy_tension = fracture_energy_tension / rPDParameters.CharacteristicLength;
    const double characteristic_fracture_energy_compression = fracture_energy_compression / rPDParameters.CharacteristicLength;

    return 1.0 / (tension_parameter / characteristic_fracture_energy_tension + compression_parameter / characteristic_fracture_energy_compression);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void AssociativePlasticDamageModel<TYieldSurfaceType>::CalculateAnalyticalTangentTensor(
    ConstitutiveLaw::Parameters& rValues,
    PlasticDamageParameters &rPDParameters
    )
{
    const BoundedMatrixType& r_C = rPDParameters.ConstitutiveMatrix;
    const double xi = rPDParameters.PlasticDamageProportion;
    const BoundedVectorType& r_plastic_flow = rPDParameters.PlasticFlow;
    const BoundedVectorType& r_stress = rPDParameters.StressVector;
    const double denominator = CalculatePlasticDenominator(rValues, rPDParameters);
    const BoundedMatrixType aux_compliance_incr = outer_prod(r_plastic_flow,r_plastic_flow) / inner_prod(r_plastic_flow, r_stress);
    const BoundedVectorType left_vector = (1.0 - xi) * prod(r_C, r_plastic_flow) + xi * prod(Matrix(prod(r_C, aux_compliance_incr)), r_stress);
    const BoundedVectorType right_vector = prod(r_C, r_plastic_flow);
    noalias(rPDParameters.TangentTensor) = r_C -outer_prod(right_vector, left_vector) / denominator;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
double AssociativePlasticDamageModel<TYieldSurfaceType>::CalculatePlasticDenominator(
    ConstitutiveLaw::Parameters& rValues,
    PlasticDamageParameters &rParam)
{
    const double g = CalculateVolumetricFractureEnergy(rValues.GetMaterialProperties(), rParam);
    const BoundedVectorType& r_plastic_flow = rParam.PlasticFlow;
    const BoundedVectorType& r_stress = rParam.StressVector;
    const double slope = rParam.Slope; // d(Threshold) / d(plastic_dissipation)
    const BoundedMatrixType& r_C = rParam.ConstitutiveMatrix;
    const double xi = rParam.PlasticDamageProportion;

    // Plasticity terms
    const double A = inner_prod(r_plastic_flow, prod(r_C, r_plastic_flow))*(1.0 - xi);
    const double B = (1.0 - xi) * (1.0 / g) * slope * inner_prod(r_plastic_flow, r_stress);

    // Damage terms
    const BoundedMatrixType aux_compliance_incr = outer_prod(r_plastic_flow,r_plastic_flow) / inner_prod(r_plastic_flow, r_stress);
    const BoundedMatrixType aux_mat = prod(r_C, aux_compliance_incr);
    const double C = xi*inner_prod(r_plastic_flow, prod(aux_mat, r_stress));
    const double D = (0.5 * slope * xi / g)*inner_prod(r_stress, prod(aux_compliance_incr, r_stress));

    return A + B + C + D;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void AssociativePlasticDamageModel<TYieldSurfaceType>::CalculatePlasticDissipationIncrement(
    const Properties &rMaterialProperties,
    PlasticDamageParameters &rParam
    )
{
    const double g = CalculateVolumetricFractureEnergy(rMaterialProperties, rParam);
    rParam.PlasticDissipationIncrement = inner_prod(rParam.StressVector,
                            rParam.PlasticStrainIncrement) / g;
    rParam.PlasticDissipationIncrement = MacaullyBrackets(rParam.PlasticDissipationIncrement);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void AssociativePlasticDamageModel<TYieldSurfaceType>::CalculateDamageDissipationIncrement(
    const Properties &rMaterialProperties,
    PlasticDamageParameters &rParam
    )
{
    const double g = CalculateVolumetricFractureEnergy(rMaterialProperties, rParam);
    rParam.DamageDissipationIncrement = inner_prod(rParam.StressVector,
                            prod(rParam.ComplianceMatrixIncrement,rParam.StressVector)) * 0.5 / g;

    rParam.DamageDissipationIncrement = MacaullyBrackets(rParam.DamageDissipationIncrement);
}

/***********************************************************************************/
/***********************************************************************************/
template<class TYieldSurfaceType>
typename AssociativePlasticDamageModel<TYieldSurfaceType>::ResidualFunctionType
AssociativePlasticDamageModel<TYieldSurfaceType>::ExponentialSofteningImplicitFunction(
    )
{
    ResidualFunctionType implicit_function = [](const double Dissipation, const double Threshold, ConstitutiveLaw::Parameters &rValues, PlasticDamageParameters &rPDParameters)
    {
        const auto &r_mat_properties = rValues.GetMaterialProperties();
        const double xi = rPDParameters.PlasticDamageProportion;
        const double E = r_mat_properties[YOUNG_MODULUS];
        const double g = AssociativePlasticDamageModel<TYieldSurfaceType>::CalculateVolumetricFractureEnergy(r_mat_properties, rPDParameters);
        double K0;
        GenericConstitutiveLawIntegratorPlasticity<TYieldSurfaceType>::GetInitialUniaxialThreshold(rValues, K0);
        const double alpha = std::pow(K0, 2) / (2.0 * E * g);
        const double K_K0 = Threshold / K0;
        return K0 * (1.0 - Dissipation) - Threshold * (1.0 + alpha * ((1.0 - xi) * (K_K0 - 0.5 * std::log(K_K0) - 1.0) + 0.5 * std::log(K_K0)) - 0.5 * xi * std::log(K_K0));
    };
    return implicit_function;
}

/***********************************************************************************/
/***********************************************************************************/
template<class TYieldSurfaceType>
typename AssociativePlasticDamageModel<TYieldSurfaceType>::ResidualFunctionType
AssociativePlasticDamageModel<TYieldSurfaceType>::ExponentialSofteningImplicitFunctionDerivative(
    )
{
    ResidualFunctionType function_derivative = [](const double Dissipation, const double Threshold, ConstitutiveLaw::Parameters& rValues, PlasticDamageParameters &rPDParameters) {
        const auto &r_mat_properties = rValues.GetMaterialProperties();
        const double xi = rPDParameters.PlasticDamageProportion;
        const double E = r_mat_properties[YOUNG_MODULUS];
        const double g = AssociativePlasticDamageModel<TYieldSurfaceType>::CalculateVolumetricFractureEnergy(r_mat_properties, rPDParameters);
        double K0;
        GenericConstitutiveLawIntegratorPlasticity<TYieldSurfaceType>::GetInitialUniaxialThreshold(rValues, K0);
        const double alpha = std::pow(K0, 2) / (2.0 * E * g);
        const double K_K0 = Threshold / K0;
        return -(1.0 + alpha * ((1.0 - xi) * (K_K0 - 0.5 * std::log(K_K0) - 1.0) + 0.5 * std::log(K_K0)) - 0.5 * xi * std::log(K_K0)) - Threshold * (alpha * ((1.0 - xi) * (1.0 / K0 - 1.0 / (2.0 * Threshold)) + 1 / (2.0 * Threshold)) - 0.5 * xi / Threshold);
    };
    return function_derivative;
}

/***********************************************************************************/
/***********************************************************************************/
template<class TYieldSurfaceType>
typename AssociativePlasticDamageModel<TYieldSurfaceType>::ResidualFunctionType
AssociativePlasticDamageModel<TYieldSurfaceType>::ExponentialHardeningImplicitFunction(
    )
{
    ResidualFunctionType implicit_function = [](const double Dissipation, const double Threshold, ConstitutiveLaw::Parameters &rValues, PlasticDamageParameters &rPDParameters)
    {
        const auto& r_mat_properties = rValues.GetMaterialProperties();
        // we first need a plastic-damage indicator
        const double xi = rPDParameters.PlasticDamageProportion;
        double max_threshold = 0.0;
        double chi = 0.0;
        double K0;
        GenericConstitutiveLawIntegratorPlasticity<TYieldSurfaceType>::GetInitialUniaxialThreshold(rValues, K0);
        const double g = CalculateVolumetricFractureEnergy(rValues.GetMaterialProperties(), rPDParameters);
        const double E = r_mat_properties[YOUNG_MODULUS];
        const double factor = std::pow(K0, 2) / E;

        if (r_mat_properties.Has(MAXIMUM_STRESS)) {
            max_threshold = r_mat_properties[MAXIMUM_STRESS];
            chi = -std::sqrt(max_threshold / (max_threshold - K0));
        } else {
            chi = (factor + g + std::sqrt(factor * (5.0 / 4.0 * factor + 2.0 * g))) / (0.5 * factor - g);
            const double chi_square = std::pow(chi, 2);
            max_threshold = (chi_square * K0) / (chi_square - 1.0);
        }
        const double diss_indicator = factor / (2.0 * g) * (1.0 - std::pow(max_threshold / K0, 2) * (1.0 + K0 * xi / max_threshold) - xi) +
            (0.5 * factor - g) / (g * (3.0 * chi + 1.0) * (chi - 1.0)) * ((2.0 * chi + 1) - max_threshold*xi / K0 * (std::pow(chi, 2) - 1.0) * std::log(chi / (chi - 1.0)));

        const double sign = (Dissipation < diss_indicator) ? -1.0 : 1.0; // In hardening should be negative
        const double gamma = (0.5*factor - g) / (g * (3.0 * chi + 1) * (chi - 1.0));
        const double alpha = std::sqrt(std::pow(chi, 2) * (1.0 - Threshold / K0) + Threshold / K0);
        return -Dissipation + factor / (2.0 * g) * (1.0 - (std::pow(Threshold / K0, 2) * (1.0 + K0 / Threshold * xi - xi))) +
               gamma * ((1.0 + sign * alpha) * (2.0 * chi + 1.0 - sign * alpha) - Threshold * xi / K0 * (std::pow(chi, 2) - 1.0) * std::log((chi + sign * alpha) / (chi - 1.0)));
    };
    return implicit_function;
}

/***********************************************************************************/
/***********************************************************************************/
template<class TYieldSurfaceType>
typename AssociativePlasticDamageModel<TYieldSurfaceType>::ResidualFunctionType
AssociativePlasticDamageModel<TYieldSurfaceType>::ExponentialHardeningImplicitFunctionDerivative(
    )
{
    ResidualFunctionType function_derivative = [](const double Dissipation, const double Threshold, ConstitutiveLaw::Parameters &rValues, PlasticDamageParameters &rPDParameters)
    {
        const auto& r_mat_properties = rValues.GetMaterialProperties();
        // we first need a plastic-damage indicator
        const double xi = rPDParameters.PlasticDamageProportion;
        double max_threshold = 0.0;
        double chi = 0.0;
        double K0;
        GenericConstitutiveLawIntegratorPlasticity<TYieldSurfaceType>::GetInitialUniaxialThreshold(rValues, K0);
        const double g = CalculateVolumetricFractureEnergy(rValues.GetMaterialProperties(), rPDParameters);
        const double E = r_mat_properties[YOUNG_MODULUS];
        const double factor = std::pow(K0, 2) / E;
        double chi_square;

        if (r_mat_properties.Has(MAXIMUM_STRESS)) {
            max_threshold = r_mat_properties[MAXIMUM_STRESS];
            chi = -std::sqrt(max_threshold / (max_threshold - K0));
            chi_square = std::pow(chi, 2);
        } else {
            chi = (factor + g + std::sqrt(factor * (5.0 / 4.0 * factor + 2.0 * g))) / (0.5 * factor - g);
            chi_square = std::pow(chi, 2);
            max_threshold = (chi_square * K0) / (chi_square - 1.0);
        }
        const double diss_indicator = factor / (2.0 * g) * (1.0 - std::pow(max_threshold / K0, 2) * (1.0 + K0 * xi / max_threshold) - xi) +
            (0.5 * factor - g) / (g * (3.0 * chi + 1.0) * (chi - 1.0)) * ((2.0 * chi + 1) - max_threshold*xi / K0 * (std::pow(chi, 2) - 1.0) * std::log(chi / (chi - 1.0)));

        const double sign = (Dissipation < diss_indicator) ? -1.0 : 1.0; // In hardening should be negative

        const double ey = K0/E;
        return (ey * K0 * (xi / K0 - (2.0 * Threshold * ((K0 * xi) / Threshold - xi + 1.0)) / (K0 * K0))) / (2.0 * g) + ((g - (ey * K0) / 2.0) * ((-sign *
                (1.0 / K0 - chi * chi / K0) * (2.0 * chi + -sign * std::sqrt((1.0 - Threshold / K0) * (chi * chi) + Threshold / K0) + 1.0)) /
                (2.0 * std::sqrt(Threshold / K0 - (chi * chi) * (Threshold / K0 - 1.0))) + (-sign * (1.0 / K0 - (chi * chi) / K0) *
                (-sign * std::sqrt((1.0 - Threshold / K0) * (chi * chi) + Threshold / K0) - 1.0)) / (2.0 * std::sqrt(Threshold / K0 - (chi * chi) * (Threshold / K0 - 1.0))) +
                (xi * std::log((chi - (-sign) * std::sqrt((1.0 - Threshold / K0) * (chi * chi) + Threshold / K0)) / (chi - 1.0)) * ((chi * chi) - 1.0)) / K0 - (-sign * Threshold * xi *
                 (1.0 / K0 - (chi * chi) / K0) * ((chi * chi) - 1.0)) / (2.0 * K0 * std::sqrt(Threshold / K0 - (chi * chi) * (Threshold / K0 - 1)) * (chi - -sign *
                 std::sqrt((1.0 - Threshold / K0) * (chi * chi) + Threshold / K0))))) / (g * (3.0 * chi + 1.0) * (chi - 1.0));
    };
    return function_derivative;
}

/***********************************************************************************/
/***********************************************************************************/
template<class TYieldSurfaceType>
typename AssociativePlasticDamageModel<TYieldSurfaceType>::ResidualFunctionType
AssociativePlasticDamageModel<TYieldSurfaceType>::CurveByPointsHardeningImplicitFunction(
    )
{
    ResidualFunctionType implicit_function = [](const double Dissipation, const double Threshold, ConstitutiveLaw::Parameters &rValues, PlasticDamageParameters &rPDParameters)
    {
        const auto& r_mat_properties = rValues.GetMaterialProperties();
        // we first need a plastic-damage indicator
        const double xi = rPDParameters.PlasticDamageProportion;
        const double C0 = r_mat_properties[YOUNG_MODULUS];
        const double g = AssociativePlasticDamageModel<TYieldSurfaceType>::CalculateVolumetricFractureEnergy(r_mat_properties, rPDParameters);
        double K0;
        GenericConstitutiveLawIntegratorPlasticity<TYieldSurfaceType>::GetInitialUniaxialThreshold(rValues, K0);
        const double E0 = K0 / C0; // Yield strain
        const Vector& s_by_points = r_mat_properties[EQUIVALENT_STRESS_VECTOR_PLASTICITY_POINT_CURVE]; // By-Points region. Stress vector
        const Vector& e_by_points = r_mat_properties[TOTAL_STRAIN_VECTOR_PLASTICITY_POINT_CURVE]; // By-Points region. Strain vector
        const SizeType size_curve = s_by_points.size();

        // Compute volumetric fracture energy in the transition by-points -> exponential
        double g_by_points = 0.5 * K0 * E0;
        IndexType i = 1;
        for (i; i < size_curve; ++i) {
            g_by_points += 0.5 * (s_by_points(i - 1) + s_by_points(i)) * (e_by_points(i) - e_by_points(i - 1));
        }
        g_by_points = g_by_points - (0.5 * xi * s_by_points(i-1) * e_by_points(i-1) + 0.5 * (1.0 - xi) * s_by_points(i-1) * s_by_points(i-1) / C0);
        const double pd_dissipation_threshold = (g_by_points) / g;
        const double A = s_by_points(size_curve - 1) / (0.5 * (1.0 - xi) * s_by_points(size_curve - 1) * s_by_points(size_curve - 1) / C0 + 0.5 * xi * s_by_points(size_curve - 1)
                        * e_by_points(size_curve - 1) - (1.0 - pd_dissipation_threshold) * g);
        return (1.0 - xi) * (s_by_points(size_curve - 1) * s_by_points(size_curve - 1) - Threshold * Threshold) / (2.0 * g * C0) + (s_by_points(size_curve - 1) - Threshold) / g
                * (0.5 * xi * e_by_points(size_curve - 1) - 1.0 / A) - xi * Threshold / (2.0 * g * A) * std::log(Threshold / s_by_points(size_curve - 1))
                + pd_dissipation_threshold - Dissipation;
    };
    return implicit_function;
}

/***********************************************************************************/
/***********************************************************************************/
template<class TYieldSurfaceType>
typename AssociativePlasticDamageModel<TYieldSurfaceType>::ResidualFunctionType
AssociativePlasticDamageModel<TYieldSurfaceType>::CurveByPointsHardeningImplicitFunctionDerivative(
    )
{
    ResidualFunctionType function_derivative = [](const double Dissipation, const double Threshold, ConstitutiveLaw::Parameters &rValues, PlasticDamageParameters &rPDParameters)
    {
        const auto& r_mat_properties = rValues.GetMaterialProperties();
        // we first need a plastic-damage indicator
        const double xi = rPDParameters.PlasticDamageProportion;
        const double C0 = r_mat_properties[YOUNG_MODULUS];
        const double g = CalculateVolumetricFractureEnergy(rValues.GetMaterialProperties(), rPDParameters);
        double K0;
        GenericConstitutiveLawIntegratorPlasticity<TYieldSurfaceType>::GetInitialUniaxialThreshold(rValues, K0);
        const double E0 = K0 / C0; // Yield strain
        const Vector& s_by_points = r_mat_properties[EQUIVALENT_STRESS_VECTOR_PLASTICITY_POINT_CURVE]; // By-Points region. Stress vector
        const Vector& e_by_points = r_mat_properties[TOTAL_STRAIN_VECTOR_PLASTICITY_POINT_CURVE]; // By-Points region. Strain vector
        const SizeType size_curve = s_by_points.size();

        // Compute volumetric fracture energy in the transition by-points -> exponential
        double g_by_points = 0.5 * K0 * E0;
        IndexType i = 1;
        for (i; i < size_curve; ++i) {
            g_by_points += 0.5 * (s_by_points(i - 1) + s_by_points(i)) * (e_by_points(i) - e_by_points(i - 1));
        }
        g_by_points = g_by_points - (0.5 * xi * s_by_points(i-1) * e_by_points(i-1) + 0.5 * (1.0 - xi) * s_by_points(i-1) * s_by_points(i-1) / C0);
        const double pd_dissipation_threshold = (g_by_points) / g;

        const double A = s_by_points(size_curve - 1) / (0.5 * (1.0 - xi) * s_by_points(size_curve - 1) * s_by_points(size_curve - 1) / C0 + 0.5 * xi * s_by_points(size_curve - 1)
                        * e_by_points(size_curve - 1) - (1.0 - pd_dissipation_threshold) * g);

        return (1.0 / A - (1.0 - xi) * Threshold * E0 / K0 - 0.5 * xi * (e_by_points(size_curve - 1) + (1.0 + std::log(Threshold / s_by_points(size_curve - 1))) / A)) / g;
    };
    return function_derivative;
}

/***********************************************************************************/
/***********************************************************************************/
template<class TYieldSurfaceType>
void AssociativePlasticDamageModel<TYieldSurfaceType>::CalculateThresholdAndSlope(
    ConstitutiveLaw::Parameters& rValues,
    PlasticDamageParameters &rPDParameters
    )
{
    const auto& r_mat_properties = rValues.GetMaterialProperties();

    if (rPDParameters.PlasticDamageProportion == 0.0) { // full plastic behaviour
        double uniaxial_plastic_strain = 0.0;
        GenericConstitutiveLawIntegratorPlasticity<TYieldSurfaceType>::
            CalculateEquivalentPlasticStrain(rPDParameters.StressVector,
            rPDParameters.UniaxialStress, rPDParameters.PlasticStrain, 0.0, rValues, uniaxial_plastic_strain);
        double tension_parameter, compression_parameter;
        GenericConstitutiveLawIntegratorPlasticity<TYieldSurfaceType>::CalculateIndicatorsFactors(
            rPDParameters.StressVector, tension_parameter,compression_parameter);

        GenericConstitutiveLawIntegratorPlasticity<TYieldSurfaceType>::
            CalculateEquivalentStressThreshold(rPDParameters.TotalDissipation,
            tension_parameter, compression_parameter, rPDParameters.Threshold, rPDParameters.Slope, rValues,
            uniaxial_plastic_strain, rPDParameters.CharacteristicLength);
    } else { // plastic-damage combinations
        const int curve_type = r_mat_properties[HARDENING_CURVE];
        switch (static_cast<typename GenericConstitutiveLawIntegratorPlasticity<TYieldSurfaceType>::HardeningCurveType>(curve_type)) {
            case GenericConstitutiveLawIntegratorPlasticity<TYieldSurfaceType>::HardeningCurveType::LinearSoftening: {
                const double xi = rPDParameters.PlasticDamageProportion;
                double K0;
                GenericConstitutiveLawIntegratorPlasticity<TYieldSurfaceType>::GetInitialUniaxialThreshold(rValues, K0);
                if (xi == 1.0) { //Necessary to avoid indetermination
                    rPDParameters.Threshold = K0 * (1.0 - rPDParameters.TotalDissipation);
                    rPDParameters.Slope = -K0;
                } else {
                    rPDParameters.Threshold = K0 * (std::sqrt((2.0 - xi) * (2.0 - xi) - 4.0 * rPDParameters.TotalDissipation * (1.0 - xi)) - xi) / (2.0 * (1.0 - xi));
                    rPDParameters.Slope = -K0 / (std::sqrt((2.0 - xi) * (2.0 - xi) - 4.0 * rPDParameters.TotalDissipation * (1.0 - xi)));
                }
                break;
            }
            case GenericConstitutiveLawIntegratorPlasticity<TYieldSurfaceType>::HardeningCurveType::ExponentialSoftening: {
                ResidualFunctionType implicit_function   = ExponentialSofteningImplicitFunction();
                ResidualFunctionType function_derivative = ExponentialSofteningImplicitFunctionDerivative();
                rPDParameters.Threshold = CalculateThresholdImplicitExpression(implicit_function, function_derivative, rValues, rPDParameters);
                rPDParameters.Slope     = CalculateSlopeFiniteDifferences(implicit_function, function_derivative, rValues, rPDParameters);
                break;
            }
            case GenericConstitutiveLawIntegratorPlasticity<TYieldSurfaceType>::HardeningCurveType::InitialHardeningExponentialSoftening: {
                ResidualFunctionType implicit_function   = ExponentialHardeningImplicitFunction();
                ResidualFunctionType function_derivative = ExponentialHardeningImplicitFunctionDerivative();
                double K0;
                GenericConstitutiveLawIntegratorPlasticity<TYieldSurfaceType>::GetInitialUniaxialThreshold(rValues, K0);
                const double g = CalculateVolumetricFractureEnergy(rValues.GetMaterialProperties(), rPDParameters);
                const double E = rValues.GetMaterialProperties()[YOUNG_MODULUS];
                const double factor = std::pow(K0, 2) / E;
                const double chi = (factor + g + std::sqrt(factor * (5.0 / 4.0 * factor + 2.0 * g))) / (0.5 * factor - g);
                const double chi_square = std::pow(chi, 2);
                const double max_threshold = (chi_square * K0) / (chi_square - 1.0);
                const double limit_factor = 1.0-1.0e-15;
                rPDParameters.Threshold = CalculateThresholdImplicitExpression(implicit_function, function_derivative, rValues, rPDParameters, max_threshold*limit_factor);
                rPDParameters.Slope     = CalculateSlopeFiniteDifferences(implicit_function, function_derivative, rValues, rPDParameters, max_threshold*limit_factor);
                break;
            }
            case GenericConstitutiveLawIntegratorPlasticity<TYieldSurfaceType>::HardeningCurveType::CurveDefinedByPoints: {
                const double xi = rPDParameters.PlasticDamageProportion;
                const double g = CalculateVolumetricFractureEnergy(rValues.GetMaterialProperties(), rPDParameters);
                const double C0 = rValues.GetMaterialProperties()[YOUNG_MODULUS];
                double K0;
                GenericConstitutiveLawIntegratorPlasticity<TYieldSurfaceType>::GetInitialUniaxialThreshold(rValues, K0);
                const double E0 = K0 / C0; // Yield strain
                const Vector& s_by_points = rValues.GetMaterialProperties()[EQUIVALENT_STRESS_VECTOR_PLASTICITY_POINT_CURVE]; // By-Points region. Stress vector
                const Vector& e_by_points = rValues.GetMaterialProperties()[TOTAL_STRAIN_VECTOR_PLASTICITY_POINT_CURVE]; // By-Points region. Strain vector
                const SizeType size_curve = s_by_points.size();

                // Compute volumetric fracture energies of each region
                double g_by_points = 0.5 * K0 * E0;
                IndexType i = 1;
                for (i; i < size_curve; ++i) {
                    g_by_points += 0.5 * (s_by_points(i - 1) + s_by_points(i)) * (e_by_points(i) - e_by_points(i - 1));
                }
                g_by_points = g_by_points - (0.5 * xi * s_by_points(i-1) * e_by_points(i-1) + 0.5 * (1.0 - xi) * s_by_points(i-1) * s_by_points(i-1) / C0);
                const double g_exponential = g - g_by_points;

                KRATOS_ERROR_IF(g_exponential < 0.0) << "Fracture energy too low in CurveDefinedByPoints of plasticity..."  << std::endl;

                // Compute plastic-damage threshold (pd_dissipation) between by-points and exponential parts
                const double pd_dissipation_threshold = (g_by_points) / g;
                if (rPDParameters.TotalDissipation < pd_dissipation_threshold) {
                    IndexType i = 0;
                    double pd_dissipation_previous_point = 0.0;
                    double pd_dissipation_next_point = pd_dissipation_previous_point;
                    // double elastic_recovery = 0.0
                    while (rPDParameters.TotalDissipation >= pd_dissipation_next_point) {
                        i += 1;
                        pd_dissipation_previous_point = pd_dissipation_next_point;
                        pd_dissipation_next_point += 0.5 / g * ((s_by_points(i) + s_by_points(i-1)) * (e_by_points(i) - e_by_points(i-1))
                                                    - xi * (s_by_points(i) * e_by_points(i) - s_by_points(i-1) * e_by_points(i-1))
                                                    - (1.0 - xi) * (s_by_points(i) * s_by_points(i) - s_by_points(i-1) * s_by_points(i-1)) * E0 / K0);
                    }
                    if (s_by_points(i) == s_by_points(i-1)) {
                        rPDParameters.Threshold = s_by_points(i);
                        rPDParameters.Slope = 0.0;
                    } else {
                        if (xi == 1.0) {
                            rPDParameters.Threshold = s_by_points(i-1) + (2.0 * g * (rPDParameters.TotalDissipation - pd_dissipation_previous_point)) / (s_by_points(i-1) * (e_by_points(i) - e_by_points(i-1)) / (s_by_points(i) - s_by_points(i-1)) - e_by_points(i-1));
                            rPDParameters.Slope = 2.0 * g * (s_by_points(i) - s_by_points(i-1)) / (e_by_points(i) * s_by_points(i-1) - e_by_points(i-1) * s_by_points(i));
                        } else {
                            const double A = (1.0 - xi) * ((e_by_points(i) - e_by_points(i-1)) / (s_by_points(i) - s_by_points(i-1)) - E0 / K0);
                            const double B = xi * (s_by_points(i-1) * (e_by_points(i) - e_by_points(i-1)) / (s_by_points(i) - s_by_points(i-1)) - e_by_points(i-1));
                            // const double current_dissipation = (rPDParameters.TotalDissipation == 0.0) ? (pd_dissipation_next_point - pd_dissipation_previous_point) * 0.1 : rPDParameters.TotalDissipation;
                            const double C = s_by_points(i-1) * (xi * e_by_points(i-1) - s_by_points(i-1) * ((e_by_points(i) - e_by_points(i-1)) / (s_by_points(i) - s_by_points(i-1)) - (1.0-xi) * E0 / K0)) - 2.0 * g * (rPDParameters.TotalDissipation - pd_dissipation_previous_point);
                            if (s_by_points(i-1) <= s_by_points(i)) { // While hardening positive square root
                                rPDParameters.Threshold = (rPDParameters.TotalDissipation == 0.0) ? s_by_points(i-1) : (- B + std::sqrt(B * B - 4.0 * A * C)) / (2.0 * A);
                                rPDParameters.Slope = 2.0 * g / (std::sqrt(B * B - 4.0 * A * C));
                            } else { // While softening negative square root
                                rPDParameters.Threshold = (rPDParameters.TotalDissipation == 0.0) ? s_by_points(i-1) : (- B - std::sqrt(B * B - 4.0 * A * C)) / (2.0 * A);
                                rPDParameters.Slope = - 2.0 * g / (std::sqrt(B * B - 4.0 * A * C));
                            }
                        }
                    }
                } else { //Exponential branch included to achieve consistent results after full plasticity scenarios
                    ResidualFunctionType implicit_function   = CurveByPointsHardeningImplicitFunction();
                    ResidualFunctionType function_derivative = CurveByPointsHardeningImplicitFunctionDerivative();
                    rPDParameters.Threshold = CalculateThresholdImplicitExpression(implicit_function, function_derivative, rValues, rPDParameters);
                    rPDParameters.Slope     = CalculateSlopeFiniteDifferences(implicit_function, function_derivative, rValues, rPDParameters);
                }
                break;
            }
            default: {
                KRATOS_ERROR << "Invalid HARDENING_CURVE defined..." << std::endl;
            }
        }
    }
}
/***********************************************************************************/
/***********************************************************************************/
template<class TYieldSurfaceType>
double AssociativePlasticDamageModel<TYieldSurfaceType>::CalculateThresholdImplicitExpression(
    ResidualFunctionType& rF,
    ResidualFunctionType& rdF_dk,
    ConstitutiveLaw::Parameters& rValues,
    PlasticDamageParameters &rPDParameters,
    const double MaxThreshold
)
{
    double old_threshold = rPDParameters.Threshold;
    const double perturbation = 1.0e-4;
    if (std::abs(rdF_dk(rPDParameters.TotalDissipation, old_threshold, rValues, rPDParameters)) < machine_tolerance) {
        old_threshold += perturbation* rPDParameters.Threshold;
        if (old_threshold >= MaxThreshold)
            old_threshold -= 2.0 * perturbation * rPDParameters.Threshold;
    }
    double residual = 1.0;
    double rel_residual = 1.0;
    double new_threshold = 0.0;
    double derivative = 0.0;
    int max_iter = 2000;
    int iteration = 0;
    const double nr_tol = 1.0e-12;

    while (residual > nr_tol && iteration < max_iter && rel_residual > nr_tol) {
        derivative = rdF_dk(rPDParameters.TotalDissipation, old_threshold, rValues, rPDParameters);
        if (std::abs(derivative) > 0.0) {
            new_threshold = old_threshold - (1.0 / derivative) * rF(rPDParameters.TotalDissipation, old_threshold, rValues, rPDParameters);
            if (new_threshold >= MaxThreshold) {
                new_threshold = MaxThreshold;
                break;
            }
        } else {
            break;
        }
        rel_residual = std::abs(new_threshold - old_threshold);
        residual = rF(rPDParameters.TotalDissipation, new_threshold, rValues, rPDParameters);
        old_threshold = new_threshold;
        iteration++;
    }
    KRATOS_WARNING_IF("AssociativePlasticDamageModel", iteration == max_iter) << "Inner Newton-Raphson to find an updated threshold did not converge..." << " Tolerance achieved: " << residual << std::endl;

    return new_threshold;
}
/***********************************************************************************/
/***********************************************************************************/
template<class TYieldSurfaceType>
double AssociativePlasticDamageModel<TYieldSurfaceType>::CalculateSlopeFiniteDifferences(
    ResidualFunctionType& rF,
    ResidualFunctionType& rdF_dk,
    ConstitutiveLaw::Parameters& rValues,
    PlasticDamageParameters &rPDParameters,
    const double MaxThreshold
)
{
    const double current_threshold = rPDParameters.Threshold;
    const double perturbation = std::sqrt(tolerance);

    rPDParameters.TotalDissipation += perturbation;
    const double perturbed_threshold = CalculateThresholdImplicitExpression(rF, rdF_dk, rValues, rPDParameters, MaxThreshold);
    rPDParameters.TotalDissipation -= perturbation;

    return (perturbed_threshold - current_threshold) / perturbation;
}

/***********************************************************************************/
/***********************************************************************************/
template<class TYieldSurfaceType>
void AssociativePlasticDamageModel<TYieldSurfaceType>::CalculateFlowVector(
    ConstitutiveLaw::Parameters& rValues,
    PlasticDamageParameters &rPDParameters
    )
{
    BoundedVectorType deviator;
    double J2;
    const BoundedVectorType& r_stress = rPDParameters.StressVector;
    const double I1 = r_stress[0] + r_stress[1] + r_stress[2];
    ConstitutiveLawUtilities<VoigtSize>::CalculateJ2Invariant(r_stress, I1, deviator, J2);
    TYieldSurfaceType::CalculateYieldSurfaceDerivative(r_stress, deviator, J2, rPDParameters.PlasticFlow, rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void AssociativePlasticDamageModel<TYieldSurfaceType>::CalculatePlasticStrainIncrement(
    ConstitutiveLaw::Parameters& rValues,
    PlasticDamageParameters &rPDParameters
    )
{
    noalias(rPDParameters.PlasticStrainIncrement) = (1.0 - rPDParameters.PlasticDamageProportion) *
        rPDParameters.PlasticConsistencyIncrement * rPDParameters.PlasticFlow;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void AssociativePlasticDamageModel<TYieldSurfaceType>::CalculateComplianceMatrixIncrement(
    ConstitutiveLaw::Parameters& rValues,
    PlasticDamageParameters &rPDParameters
    )
{
    const BoundedVectorType& plastic_flow = rPDParameters.PlasticFlow;

    const double denominator = inner_prod(plastic_flow, rPDParameters.StressVector);

    if (std::abs(denominator) > machine_tolerance)
        noalias(rPDParameters.ComplianceMatrixIncrement) = (rPDParameters.PlasticDamageProportion) * rPDParameters.PlasticConsistencyIncrement * outer_prod(plastic_flow, plastic_flow) / denominator;
    else
        noalias(rPDParameters.ComplianceMatrixIncrement) = ZeroMatrix(VoigtSize, VoigtSize);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void AssociativePlasticDamageModel<TYieldSurfaceType>::CalculatePlasticConsistencyIncrement(
    ConstitutiveLaw::Parameters& rValues,
    PlasticDamageParameters &rPDParameters
    )
{
    const double denominator = CalculatePlasticDenominator(rValues, rPDParameters);

    if (std::abs(denominator) > machine_tolerance)
        rPDParameters.PlasticConsistencyIncrement = (rPDParameters.NonLinearIndicator) / denominator;
    else
        rPDParameters.PlasticConsistencyIncrement = 0.0;

    rPDParameters.PlasticConsistencyIncrement = MacaullyBrackets(rPDParameters.PlasticConsistencyIncrement);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void AssociativePlasticDamageModel<TYieldSurfaceType>::IntegrateStressPlasticDamageMechanics(
    ConstitutiveLaw::Parameters& rValues,
    PlasticDamageParameters &rPDParameters
    )
{
    KRATOS_TRY;
    BoundedMatrixType constitutive_matrix_increment;
    noalias(rPDParameters.TangentTensor) = rPDParameters.ConstitutiveMatrix;
    const auto& r_material_properties = rValues.GetMaterialProperties();
    const bool crack_reclosure = (r_material_properties.Has(CRACK_RECLOSING)) ? r_material_properties[CRACK_RECLOSING] : false;

    bool is_converged = false;
    IndexType iteration = 0, max_iter = r_material_properties.Has(MAX_NUMBER_NL_CL_ITERATIONS) ? r_material_properties.GetValue(MAX_NUMBER_NL_CL_ITERATIONS) : 1000;

    const double splits = 1;
    double tension_parameter, compression_parameter, det;
    noalias(rPDParameters.StressVector) = prod(rPDParameters.ConstitutiveMatrix, mOldStrain - rPDParameters.PlasticStrain);


    for (SizeType i = 0; i < splits; i++) {
        rPDParameters.StressVector += prod(rPDParameters.ConstitutiveMatrix, rPDParameters.StrainVector - mOldStrain) / splits;
        is_converged = false;
        iteration = 0;
        TYieldSurfaceType::CalculateEquivalentStress(rPDParameters.StressVector, rPDParameters.StrainVector, rPDParameters.UniaxialStress, rValues);
        rPDParameters.NonLinearIndicator = rPDParameters.UniaxialStress - rPDParameters.Threshold;
        if (rPDParameters.NonLinearIndicator > tolerance*rPDParameters.Threshold) {
            while (is_converged == false && iteration <= max_iter) {
                CalculateThresholdAndSlope(rValues, rPDParameters);
                CalculateFlowVector(rValues, rPDParameters);
                CalculatePlasticConsistencyIncrement(rValues, rPDParameters);

                // Update the analytical tangent tensor
                if (rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
                    CalculateAnalyticalTangentTensor(rValues, rPDParameters);

                // Compute the plastic strain increment
                CalculatePlasticStrainIncrement(rValues, rPDParameters);
                noalias(rPDParameters.PlasticStrain) += rPDParameters.PlasticStrainIncrement;

                // Compute the compliance increment -> C dot
                CalculateComplianceMatrixIncrement(rValues, rPDParameters);

                noalias(rPDParameters.StressVector) -= rPDParameters.PlasticConsistencyIncrement * prod(rPDParameters.ConstitutiveMatrix, rPDParameters.PlasticFlow);

                if (crack_reclosure) {
                    // now we check the tensile-compressive distribution
                    GenericConstitutiveLawIntegratorPlasticity<TYieldSurfaceType>::CalculateIndicatorsFactors(rPDParameters.StressVector, tension_parameter,compression_parameter);
                    noalias(rPDParameters.ComplianceMatrix) += tension_parameter*rPDParameters.ComplianceMatrixIncrement;
                    noalias(rPDParameters.ComplianceMatrixCompression) += compression_parameter*rPDParameters.ComplianceMatrixIncrement;
                    BoundedMatrixType temp;
                    noalias(temp) = tension_parameter * rPDParameters.ComplianceMatrix + compression_parameter * rPDParameters.ComplianceMatrixCompression;
                    MathUtils<double>::InvertMatrix(temp, rPDParameters.ConstitutiveMatrix, det);

                } else {
                    noalias(rPDParameters.ComplianceMatrix) += rPDParameters.ComplianceMatrixIncrement;
                    CalculateConstitutiveMatrix(rValues, rPDParameters);
                }

                // Compute the non-linear dissipation performed
                CalculatePlasticDissipationIncrement(r_material_properties, rPDParameters);
                CalculateDamageDissipationIncrement(r_material_properties, rPDParameters);
                AddNonLinearDissipation(rPDParameters);

                // updated uniaxial and threshold stress check
                TYieldSurfaceType::CalculateEquivalentStress(rPDParameters.StressVector, rPDParameters.StrainVector, rPDParameters.UniaxialStress, rValues);
                CalculateThresholdAndSlope(rValues, rPDParameters);
                rPDParameters.NonLinearIndicator = rPDParameters.UniaxialStress - rPDParameters.Threshold;

                if (rPDParameters.NonLinearIndicator <= tolerance*rPDParameters.Threshold || rPDParameters.TotalDissipation > 0.999) {
                    is_converged = true;
                } else {
                    iteration++;
                }
            }
        }
    }
    KRATOS_WARNING_IF("GenericConstitutiveLawIntegratorPlasticDamage", iteration > max_iter) << "Maximum number of iterations in plastic-damage loop reached..." << std::endl;
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void AssociativePlasticDamageModel<TYieldSurfaceType>::CalculateConstitutiveMatrix(
    ConstitutiveLaw::Parameters& rValues,
    PlasticDamageParameters &rPDParameters
    )
{
    double det = 0.0;
    MathUtils<double>::InvertMatrix(rPDParameters.ComplianceMatrix, rPDParameters.ConstitutiveMatrix, det);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void AssociativePlasticDamageModel<TYieldSurfaceType>::UpdateInternalVariables(
    PlasticDamageParameters &rPDParameters
    )
{
    mPlasticDissipation        = rPDParameters.PlasticDissipation;
    mDamageDissipation         = rPDParameters.DamageDissipation;
    mThreshold                 = rPDParameters.Threshold;
    noalias(mPlasticStrain)    = rPDParameters.PlasticStrain;
    noalias(mComplianceMatrix) = rPDParameters.ComplianceMatrix;
    noalias(mComplianceMatrixCompression) = rPDParameters.ComplianceMatrixCompression;
    noalias(mOldStrain) = rPDParameters.StrainVector;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void AssociativePlasticDamageModel<TYieldSurfaceType>::CheckMinimumFractureEnergy(
    ConstitutiveLaw::Parameters& rValues,
    PlasticDamageParameters &rPDParameters
    )
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const bool is_yield_symmetric = r_material_properties.Has(YIELD_STRESS_TENSION) ? false : true;

    const double young_modulus = r_material_properties[YOUNG_MODULUS];
    const double yield = (is_yield_symmetric == false) ? r_material_properties[YIELD_STRESS_TENSION] : r_material_properties[YIELD_STRESS];
    const double fracture_energy = r_material_properties[FRACTURE_ENERGY];

    const double hlim = 2.0 * young_modulus * fracture_energy / (std::pow(yield, 2));
    KRATOS_ERROR_IF(rPDParameters.CharacteristicLength > hlim) << "The Fracture Energy is to low: " << rPDParameters.CharacteristicLength << std::endl;

    if (!is_yield_symmetric) { // Check frac energy in compression
        const double yield_compression =  r_material_properties[YIELD_STRESS_COMPRESSION];
        const double fracture_energy_compr = r_material_properties[FRACTURE_ENERGY_COMPRESSION];
        const double hlim_compr = 2.0 * young_modulus * fracture_energy_compr / (std::pow(yield_compression, 2));
        KRATOS_ERROR_IF(rPDParameters.CharacteristicLength > hlim_compr) << "The Fracture Energy in compression is to low: " << rPDParameters.CharacteristicLength << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void AssociativePlasticDamageModel<TYieldSurfaceType>::CalculateElasticComplianceMatrix(
    BoundedMatrixType& rC,
    const Properties& rMaterialProperties
    )
{
    const double E = rMaterialProperties[YOUNG_MODULUS];
    const double NU = rMaterialProperties[POISSON_RATIO];

    noalias(rC) = ZeroMatrix(VoigtSize,VoigtSize);

    const double G = E / (2.0 * (1.0 + NU));
    const double c1 = 1.0 / E;
    const double c2 = -NU / E;
    const double c3 = 1.0 / G;

    rC(0,0) = c1; rC(0,1) = c2; rC(0,2) = c2;
    rC(1,0) = c2; rC(1,1) = c1; rC(1,2) = c2;
    rC(2,0) = c2; rC(2,1) = c2; rC(2,2) = c1;

    rC(3,3) = c3;
    rC(4,4) = c3;
    rC(5,5) = c3;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
bool AssociativePlasticDamageModel<TYieldSurfaceType>::Has(
    const Variable<double>& rThisVariable
    )
{
    bool has = false;
    if (rThisVariable == PLASTIC_DISSIPATION) {
        has = true;
    } else if (rThisVariable == THRESHOLD) {
        has = true;
    } else if (rThisVariable == DAMAGE) {
        has = true;
    } else if (rThisVariable == DISSIPATION) {
        has = true;
    }
    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
bool AssociativePlasticDamageModel<TYieldSurfaceType>::Has(
    const Variable<Vector>& rThisVariable
    )
{
    bool has = false;
    if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
        has = true;
    }
    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
double& AssociativePlasticDamageModel<TYieldSurfaceType>::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    rValue = 0.0;
    if (rThisVariable == PLASTIC_DISSIPATION) {
        rValue = mPlasticDissipation;
    } else if (rThisVariable == THRESHOLD) {
        rValue = mThreshold;
    }  else if (rThisVariable == DAMAGE) {
        rValue = mDamageDissipation;
    } else if (rThisVariable == DISSIPATION) {
        rValue = mPlasticDissipation + mDamageDissipation;
    }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
Vector& AssociativePlasticDamageModel<TYieldSurfaceType>::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    rValue.resize(VoigtSize, false);
    rValue.clear();
    if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
        noalias(rValue) = mPlasticStrain;
    }
    return rValue;
}


/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void AssociativePlasticDamageModel<TYieldSurfaceType>::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (rThisVariable == PLASTIC_DISSIPATION) {
        mPlasticDissipation = rValue;
    } else if (rThisVariable == THRESHOLD) {
        mThreshold = rValue;
    }  else if (rThisVariable == DAMAGE) {
        mDamageDissipation = rValue;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void AssociativePlasticDamageModel<TYieldSurfaceType>::SetValue(
    const Variable<Vector>& rThisVariable,
    const Vector& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
        noalias(mPlasticStrain) = rValue;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
double& AssociativePlasticDamageModel<TYieldSurfaceType>::CalculateValue(
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

        BoundedVectorType aux_stress_vector = r_stress_vector;
        TYieldSurfaceType::CalculateEquivalentStress( aux_stress_vector, r_strain_vector, rValue, rParameterValues);

        // Previous flags restored
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );
    } else {
        BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
    }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void AssociativePlasticDamageModel<TYieldSurfaceType>::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    // We construct the CL parameters
    ProcessInfo dummy_process_info;
    ConstitutiveLaw::Parameters aux_param(rElementGeometry, rMaterialProperties,dummy_process_info);

    // We call the integrator
    double initial_threshold;
    TYieldSurfaceType::GetInitialUniaxialThreshold(aux_param, initial_threshold);
    mThreshold = initial_threshold;

    BoundedMatrixType initial_compliance;
    CalculateElasticComplianceMatrix(initial_compliance, rMaterialProperties);
    noalias(mComplianceMatrix) = initial_compliance;
    noalias(mComplianceMatrixCompression) = initial_compliance;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void AssociativePlasticDamageModel<TYieldSurfaceType>::CalculateMaterialResponsePK1(
    ConstitutiveLaw::Parameters& rValues
    )
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void  AssociativePlasticDamageModel<TYieldSurfaceType>::CalculateMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rValues
    )
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void AssociativePlasticDamageModel<TYieldSurfaceType>::CalculateMaterialResponseKirchhoff(
    ConstitutiveLaw::Parameters& rValues
    )
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void AssociativePlasticDamageModel<TYieldSurfaceType>::InitializeMaterialResponsePK1(
    ConstitutiveLaw::Parameters& rValues
    )
{
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void AssociativePlasticDamageModel<TYieldSurfaceType>::InitializeMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rValues
    )
{
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void AssociativePlasticDamageModel<TYieldSurfaceType>::InitializeMaterialResponseKirchhoff(
    ConstitutiveLaw::Parameters& rValues
    )
{
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void AssociativePlasticDamageModel<TYieldSurfaceType>::InitializeMaterialResponseCauchy(
    ConstitutiveLaw::Parameters& rValues
    )
{
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void AssociativePlasticDamageModel<TYieldSurfaceType>::FinalizeMaterialResponsePK1(
    ConstitutiveLaw::Parameters& rValues
    )
{
    FinalizeMaterialResponseCauchy(rValues);
}

template<class TYieldSurfaceType>
/***********************************************************************************/
/***********************************************************************************/

void AssociativePlasticDamageModel<TYieldSurfaceType>::FinalizeMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rValues
    )
{
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void AssociativePlasticDamageModel<TYieldSurfaceType>::FinalizeMaterialResponseKirchhoff(
    ConstitutiveLaw::Parameters& rValues
    )
{
    FinalizeMaterialResponseCauchy(rValues);
}


/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
int AssociativePlasticDamageModel<TYieldSurfaceType>::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    // The auxiliary output
    int aux_out = 0;
    KRATOS_ERROR_IF(!rMaterialProperties.Has(FRACTURE_ENERGY))           << "FRACTURE_ENERGY not provided in the material properties" << std::endl;
    KRATOS_ERROR_IF(!rMaterialProperties.Has(HARDENING_CURVE))           << "HARDENING_CURVE not provided in the material properties" << std::endl;
    KRATOS_ERROR_IF(!rMaterialProperties.Has(PLASTIC_DAMAGE_PROPORTION)) << "PLASTIC_DAMAGE_PROPORTION not provided in the material properties" << std::endl;
    return aux_out;
}


/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void AssociativePlasticDamageModel<TYieldSurfaceType>::CalculateTangentTensor(
    ConstitutiveLaw::Parameters& rValues,
    PlasticDamageParameters &rPlasticDamageParameters
    )
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();

    const bool consider_perturbation_threshold = r_material_properties.Has(CONSIDER_PERTURBATION_THRESHOLD) ?
        r_material_properties[CONSIDER_PERTURBATION_THRESHOLD] : true;
    const TangentOperatorEstimation tangent_operator_estimation = r_material_properties.Has(TANGENT_OPERATOR_ESTIMATION) ?
        static_cast<TangentOperatorEstimation>(r_material_properties[TANGENT_OPERATOR_ESTIMATION]) : TangentOperatorEstimation::SecondOrderPerturbation;

    if (tangent_operator_estimation == TangentOperatorEstimation::Analytic) {
        // CalculateAnalyticalTangentTensor(rValues, rPlasticDamageParameters);
        noalias(rValues.GetConstitutiveMatrix()) = rPlasticDamageParameters.TangentTensor;
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
        double det = 0.0;
        MathUtils<double>::InvertMatrix(mComplianceMatrix, rValues.GetConstitutiveMatrix(), det);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template class AssociativePlasticDamageModel<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>;
template class AssociativePlasticDamageModel<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>;
template class AssociativePlasticDamageModel<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>;
template class AssociativePlasticDamageModel<RankineYieldSurface<RankinePlasticPotential<6>>>;

} // Namespace Kratos
