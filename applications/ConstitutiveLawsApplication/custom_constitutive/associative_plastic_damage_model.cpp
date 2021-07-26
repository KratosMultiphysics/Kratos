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
#include "custom_constitutive/associative_plastic_damage_model.h"
#include "constitutive_laws_application_variables.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"

// Plasticity Integrator includes
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_plasticity.h"

// Yield surfaces
#include "custom_constitutive/yield_surfaces/generic_yield_surface.h"
#include "custom_constitutive/yield_surfaces/von_mises_yield_surface.h"
#include "custom_constitutive/yield_surfaces/drucker_prager_yield_surface.h"
#include "custom_constitutive/yield_surfaces/modified_mohr_coulomb_yield_surface.h"

// Plastic potentials
#include "custom_constitutive/plastic_potentials/generic_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/von_mises_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/drucker_prager_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/modified_mohr_coulomb_plastic_potential.h"

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
    // Auxiliar values
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();

    Vector& r_integrated_stress_vector = rValues.GetStressVector();
    const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::
        CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());

    if (r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        BaseType::CalculateCauchyGreenStrain(rValues, r_strain_vector);
    }

    // We compute the stress or the constitutive matrix
    if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS) ||
        r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {

        PlasticDamageParameters plastic_damage_parameters = PlasticDamageParameters();
        InitializePlasticDamageParameters(r_strain_vector, rValues.GetMaterialProperties(),
            characteristic_length, plastic_damage_parameters);

        CheckMinimumFractureEnergy(rValues, plastic_damage_parameters);

        CalculateConstitutiveMatrix(rValues, plastic_damage_parameters);
        noalias(rValues.GetConstitutiveMatrix()) = plastic_damage_parameters.ConstitutiveMatrix;

        noalias(plastic_damage_parameters.StressVector) = prod(plastic_damage_parameters.ConstitutiveMatrix,
            r_strain_vector - plastic_damage_parameters.PlasticStrain);

        TYieldSurfaceType::CalculateEquivalentStress(plastic_damage_parameters.StressVector,
            plastic_damage_parameters.StrainVector, plastic_damage_parameters.UniaxialStress, rValues);

        plastic_damage_parameters.NonLinearIndicator = plastic_damage_parameters.UniaxialStress - mThreshold;

        if (plastic_damage_parameters.NonLinearIndicator <= std::abs(1.0e-8 * mThreshold)) {
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
    // Auxiliar values
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();

    const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::
        CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());

    if (r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        BaseType::CalculateCauchyGreenStrain(rValues, r_strain_vector);
    }

    PlasticDamageParameters plastic_damage_parameters = PlasticDamageParameters();
    InitializePlasticDamageParameters(r_strain_vector, rValues.GetMaterialProperties(),
        characteristic_length, plastic_damage_parameters);

    CheckMinimumFractureEnergy(rValues, plastic_damage_parameters);

    CalculateConstitutiveMatrix(rValues, plastic_damage_parameters);
    noalias(plastic_damage_parameters.StressVector) = prod(plastic_damage_parameters.ConstitutiveMatrix,
        r_strain_vector - plastic_damage_parameters.PlasticStrain);

    TYieldSurfaceType::CalculateEquivalentStress(plastic_damage_parameters.StressVector,
        plastic_damage_parameters.StrainVector, plastic_damage_parameters.UniaxialStress, rValues);

    plastic_damage_parameters.NonLinearIndicator = plastic_damage_parameters.UniaxialStress - mThreshold;

    if (plastic_damage_parameters.NonLinearIndicator > std::abs(1.0e-8 * mThreshold)) {
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
    GenericConstitutiveLawIntegratorPlasticity<TYieldSurfaceType>::CalculateIndicatorsFactors(
        rPDParameters.StressVector, tension_parameter,compression_parameter);

    const bool has_symmetric_yield_stress = rMaterialProperties.Has(YIELD_STRESS);
    const double yield_compression = has_symmetric_yield_stress ? rMaterialProperties[YIELD_STRESS]
        : rMaterialProperties[YIELD_STRESS_COMPRESSION];
    const double yield_tension = has_symmetric_yield_stress ? rMaterialProperties[YIELD_STRESS]
        : rMaterialProperties[YIELD_STRESS_TENSION];
    const double n = yield_compression / yield_tension;
    const double fracture_energy_tension = rMaterialProperties[FRACTURE_ENERGY]; // Frac energy in tension
    const double fracture_energy_compression = rMaterialProperties.Has(FRACTURE_ENERGY_COMPRESSION) ?
        rMaterialProperties[FRACTURE_ENERGY_COMPRESSION] : rMaterialProperties[FRACTURE_ENERGY] * std::pow(n, 2);

    const double characteristic_fracture_energy_tension = fracture_energy_tension /
        rPDParameters.CharacteristicLength;
    const double characteristic_fracture_energy_compression = fracture_energy_compression /
        rPDParameters.CharacteristicLength;

    return 1.0 / (tension_parameter / characteristic_fracture_energy_tension +
        compression_parameter / characteristic_fracture_energy_compression);
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
    const double chi = rPDParameters.PlasticDamageProportion;
    const BoundedVectorType& r_plastic_flow = rPDParameters.PlasticFlow;
    const BoundedVectorType& r_stress = rPDParameters.StressVector;
    const double denominator = CalculatePlasticDenominator(rValues, rPDParameters);

    const BoundedMatrixType aux_compliance_incr = outer_prod(r_plastic_flow,r_plastic_flow) /
        inner_prod(r_plastic_flow, r_stress);

    const BoundedVectorType left_vector = (1.0 - chi) * prod(r_C, r_plastic_flow) + chi *
        prod(Matrix(prod(r_C, aux_compliance_incr)), r_stress);
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
    const double chi = rParam.PlasticDamageProportion;

    // Plasticity terms
    const double A = inner_prod(r_plastic_flow, prod(r_C, r_plastic_flow))*(1.0 - chi);
    const double B = (1.0 - chi) * (1.0 / g) * slope * inner_prod(r_plastic_flow, r_stress);

    // Damage terms
    const BoundedMatrixType aux_compliance_incr = outer_prod(r_plastic_flow,r_plastic_flow) /
        inner_prod(r_plastic_flow, r_stress);
    const BoundedMatrixType aux_mat = prod(r_C, aux_compliance_incr);
    const double C = chi*inner_prod(r_plastic_flow, prod(aux_mat, r_stress));
    const double D = (0.5 * slope * chi / g)*inner_prod(r_stress, prod(aux_compliance_incr, r_stress));

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
void AssociativePlasticDamageModel<TYieldSurfaceType>::CalculateThresholdAndSlope(
    ConstitutiveLaw::Parameters& rValues,
    PlasticDamageParameters &rPDParameters
    )
{
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
    ConstitutiveLawUtilities<VoigtSize>::CalculateJ2Invariant(r_stress,
        I1, deviator, J2);
    TYieldSurfaceType::CalculateYieldSurfaceDerivative(r_stress,
        deviator, J2, rPDParameters.PlasticFlow, rValues);
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
        noalias(rPDParameters.ComplianceMatrixIncrement) = (rPDParameters.PlasticDamageProportion) *
            rPDParameters.PlasticConsistencyIncrement * outer_prod(plastic_flow,plastic_flow) /
            denominator;
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
    const auto& r_mat_properties = rValues.GetMaterialProperties();
    BoundedMatrixType constitutive_matrix_increment;
    CalculateConstitutiveMatrix(rValues, rPDParameters);
    noalias(rPDParameters.TangentTensor) = rPDParameters.ConstitutiveMatrix;

    bool is_converged = false;
    IndexType iteration = 0, max_iter = 1000;
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
        noalias(rPDParameters.ComplianceMatrix) += rPDParameters.ComplianceMatrixIncrement;

        noalias(rPDParameters.StressVector) -= rPDParameters.PlasticConsistencyIncrement *
            prod(rPDParameters.ConstitutiveMatrix, rPDParameters.PlasticFlow);

        CalculateConstitutiveMatrix(rValues, rPDParameters);

        // Compute the non-linear dissipation performed
        CalculatePlasticDissipationIncrement(r_mat_properties, rPDParameters);
        CalculateDamageDissipationIncrement(r_mat_properties, rPDParameters);
        AddNonLinearDissipation(rPDParameters);

        // updated uniaxial and threshold stress check
        TYieldSurfaceType::CalculateEquivalentStress(rPDParameters.StressVector,
            rPDParameters.StrainVector, rPDParameters.UniaxialStress, rValues);
        CalculateThresholdAndSlope(rValues, rPDParameters);
        rPDParameters.NonLinearIndicator = rPDParameters.UniaxialStress - rPDParameters.Threshold;

        if (rPDParameters.NonLinearIndicator <= 1.0e-8*rPDParameters.Threshold) {
            is_converged = true;
        } else {
            iteration++;
        }
    }
    KRATOS_ERROR_IF(iteration > max_iter) << "Maximum number of iterations in plasticity loop reached..." << std::endl;
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
    MathUtils<double>::InvertMatrix(rPDParameters.ComplianceMatrix,
        rPDParameters.ConstitutiveMatrix, det);
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
    const double yield = (is_yield_symmetric == false) ? r_material_properties[YIELD_STRESS_TENSION] :
        r_material_properties[YIELD_STRESS];
    const double fracture_energy = r_material_properties[FRACTURE_ENERGY];

    const double hlim = 2.0 * young_modulus * fracture_energy / (std::pow(yield, 2));
    KRATOS_ERROR_IF(rPDParameters.CharacteristicLength > hlim) << "The Fracture Energy is to low: " <<
        rPDParameters.CharacteristicLength << std::endl;

    if (is_yield_symmetric == false) { // Check frac energy in compression
        const double yield_compression =  r_material_properties[YIELD_STRESS_COMPRESSION];
        const double fracture_energy_compr = r_material_properties[FRACTURE_ENERGY_COMPRESSION];
        const double hlim_compr = 2.0 * young_modulus * fracture_energy_compr / (std::pow(yield_compression, 2));
        KRATOS_ERROR_IF(rPDParameters.CharacteristicLength > hlim_compr) << "The Fracture Energy in compression is to low: " <<
            rPDParameters.CharacteristicLength << std::endl;
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
    const Variable<bool>& rThisVariable
    )
{
    bool has = false;
    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
bool AssociativePlasticDamageModel<TYieldSurfaceType>::Has(
    const Variable<int>& rThisVariable
    )
{
    bool has = false;
    return has;
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
    // At least one layer should have the value
    bool has = false;
    if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
        has = true;
    }
    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
bool AssociativePlasticDamageModel<TYieldSurfaceType>::Has(
    const Variable<Matrix>& rThisVariable
    )
{
    // At least one layer should have the value
    bool has = false;
    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
bool AssociativePlasticDamageModel<TYieldSurfaceType>::Has(
    const Variable<array_1d<double, 3 > >& rThisVariable
    )
{
    // At least one layer should have the value
    bool has = false;
    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
bool AssociativePlasticDamageModel<TYieldSurfaceType>::Has(
    const Variable<array_1d<double, 6 > >& rThisVariable
    )
{
    // At least one layer should have the value
    bool has = false;
    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
bool& AssociativePlasticDamageModel<TYieldSurfaceType>::GetValue(
    const Variable<bool>& rThisVariable,
    bool& rValue
    )
{
    // At least one layer should have the value
    rValue = false;
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
int& AssociativePlasticDamageModel<TYieldSurfaceType>::GetValue(
    const Variable<int>& rThisVariable,
    int& rValue
    )
{
    // At least one layer should have the value
    rValue = 0;
    return rValue;
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
        rValue = mDamageDissipation + mPlasticDissipation;
    } else if (rThisVariable == THRESHOLD) {
        rValue = mThreshold;
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
    // We combine the values of the layers
    rValue.clear();
    if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
        rValue = mPlasticStrain;
    }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
Matrix& AssociativePlasticDamageModel<TYieldSurfaceType>::GetValue(
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    // We combine the values of the layers
    rValue.clear();
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
array_1d<double, 3 >& AssociativePlasticDamageModel<TYieldSurfaceType>::GetValue(
    const Variable<array_1d<double, 3 >>& rThisVariable,
    array_1d<double, 3 >& rValue
    )
{
    // We combine the values of the layers
    rValue = ZeroVector(3);
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
array_1d<double, 6 >& AssociativePlasticDamageModel<TYieldSurfaceType>::GetValue(
    const Variable<array_1d<double, 6 >>& rThisVariable,
    array_1d<double, 6 >& rValue
    )
{
    // We combine the values of the layers
    rValue = ZeroVector(6);
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void AssociativePlasticDamageModel<TYieldSurfaceType>::SetValue(
    const Variable<bool>& rThisVariable,
    const bool& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void AssociativePlasticDamageModel<TYieldSurfaceType>::SetValue(
    const Variable<int>& rThisVariable,
    const int& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
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
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void AssociativePlasticDamageModel<TYieldSurfaceType>::SetValue(
    const Variable<Matrix>& rThisVariable,
    const Matrix& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void AssociativePlasticDamageModel<TYieldSurfaceType>::SetValue(
    const Variable<array_1d<double, 3 >>& rThisVariable,
    const array_1d<double, 3 >& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void AssociativePlasticDamageModel<TYieldSurfaceType>::SetValue(
    const Variable<array_1d<double, 6 >>& rThisVariable,
    const array_1d<double, 6 >& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
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
Vector& AssociativePlasticDamageModel<TYieldSurfaceType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    return BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
Matrix& AssociativePlasticDamageModel<TYieldSurfaceType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    return BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
}
/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
array_1d<double, 3 >& AssociativePlasticDamageModel<TYieldSurfaceType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<array_1d<double, 3 >>& rThisVariable,
    array_1d<double, 3 >& rValue
    )
{
    // return BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
array_1d<double, 6 >& AssociativePlasticDamageModel<TYieldSurfaceType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<array_1d<double, 6 >>& rThisVariable,
    array_1d<double, 6 >& rValue
    )
{
    // BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
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
    ConstitutiveLaw::Parameters aux_param(rElementGeometry, rMaterialProperties,
                                          dummy_process_info);

    // We call the integrator
    double initial_threshold;
    TYieldSurfaceType::GetInitialUniaxialThreshold(aux_param, initial_threshold);
    mThreshold = initial_threshold;

    BoundedMatrixType initial_compliance;
    CalculateElasticComplianceMatrix(initial_compliance, rMaterialProperties);
    noalias(mComplianceMatrix) = initial_compliance;
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
    )
{
    // The auxiliar output
    int aux_out = 0;
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
    }
}

/***********************************************************************************/
/***********************************************************************************/

template class AssociativePlasticDamageModel<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>;
template class AssociativePlasticDamageModel<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>;
template class AssociativePlasticDamageModel<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>;

} // Namespace Kratos
