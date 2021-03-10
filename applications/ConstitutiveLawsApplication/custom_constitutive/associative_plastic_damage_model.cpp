// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
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
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"

// Plasticity Integrator includes
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_plasticity.h"

// Yield surfaces
#include "custom_constitutive/yield_surfaces/generic_yield_surface.h"
#include "custom_constitutive/yield_surfaces/von_mises_yield_surface.h"

// Plastic potentials
#include "custom_constitutive/plastic_potentials/generic_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/von_mises_plastic_potential.h"

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

    // We get the constitutive tensor
    Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();

    const ProcessInfo& r_current_process_info = rValues.GetProcessInfo();

    Vector& r_integrated_stress_vector = rValues.GetStressVector();
    const double characteristic_length = ConstitutiveLawUtilities<VoigtSize>::
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

        if (plastic_damage_parameters.NonLinearIndicator <= std::abs(1.0e-4 * mThreshold)) {
            noalias(r_integrated_stress_vector) = plastic_damage_parameters.StressVector;
        } else {
            IntegrateStressPlasticDamageMechanics(rValues, plastic_damage_parameters);
            noalias(r_integrated_stress_vector) = plastic_damage_parameters.StressVector;

            if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                CalculateTangentTensor(rValues); // this modifies the ConstitutiveMatrix
                // CalculateAnalyticalTangentTensor(rValues, plastic_damage_parameters);
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

    // We get the constitutive tensor
    Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();

    const ProcessInfo& r_current_process_info = rValues.GetProcessInfo();

    Vector& r_integrated_stress_vector = rValues.GetStressVector();
    const double characteristic_length = ConstitutiveLawUtilities<VoigtSize>::
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

    if (plastic_damage_parameters.NonLinearIndicator > std::abs(1.0e-4 * mThreshold)) {
        IntegrateStressPlasticDamageMechanics(rValues, plastic_damage_parameters);
        UpdateInternalVariables(plastic_damage_parameters);
    }

}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void AssociativePlasticDamageModel<TYieldSurfaceType>::CalculateAnalyticalTangentTensor(
    ConstitutiveLaw::Parameters& rValues,
    PlasticDamageParameters &rParam
    )
{
    const double fracture_energy = rValues.GetMaterialProperties()[FRACTURE_ENERGY];
    const double g = fracture_energy / rParam.CharacteristicLength;
    BoundedVectorType plastic_flow = rParam.PlasticFlow;

    const double denominator = (inner_prod(plastic_flow, prod(rParam.ConstitutiveMatrix, plastic_flow)) +
        rParam.Slope*inner_prod(rParam.StressVector / g, plastic_flow));

    noalias(rValues.GetConstitutiveMatrix()) = rParam.ConstitutiveMatrix -
        prod(Matrix(prod(rParam.ConstitutiveMatrix, outer_prod(rParam.PlasticFlow, rParam.PlasticFlow))),
        rParam.ConstitutiveMatrix) /denominator;
}

template<class TYieldSurfaceType>
void AssociativePlasticDamageModel<TYieldSurfaceType>::CalculatePlasticDissipationIncrement(
    const Properties &rMaterialProperties,
    PlasticDamageParameters &rParam
    )
{
    const double fracture_energy = rMaterialProperties[FRACTURE_ENERGY];
    rParam.PlasticDissipationIncrement = inner_prod(rParam.StressVector,
                            rParam.PlasticStrainIncrement) *
                            rParam.CharacteristicLength / fracture_energy;
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
    const double fracture_energy = rMaterialProperties[FRACTURE_ENERGY];
    // rParam.DamageDissipationIncrement = inner_prod(rParam.StressVector,
    //                         prod(rParam.ComplianceMatrixIncrement,rParam.StressVector)) *
    //                         0.5 * rParam.CharacteristicLength / fracture_energy;

        rParam.DamageDissipationIncrement = inner_prod(rParam.StressVector,
                            prod(rParam.ComplianceMatrixIncrement,rParam.StressVector)) *
                            rParam.CharacteristicLength / fracture_energy;

    rParam.DamageDissipationIncrement = MacaullyBrackets(rParam.DamageDissipationIncrement);
}

/***********************************************************************************/
/***********************************************************************************/

// TO BE GENERALIZED TO ANY HARDENING LAW
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
    GenericConstitutiveLawIntegratorPlasticity<TYieldSurfaceType>::
        CalculateEquivalentStressThreshold(rPDParameters.TotalDissipation,
        1.0, 0.0, rPDParameters.Threshold, rPDParameters.Slope, rValues, uniaxial_plastic_strain,
        rPDParameters.CharacteristicLength);
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
    const double fracture_energy = rValues.GetMaterialProperties()[FRACTURE_ENERGY];
    const double g = fracture_energy / rPDParameters.CharacteristicLength;
    BoundedVectorType plastic_flow = rPDParameters.PlasticFlow;

    double denominator = (inner_prod(plastic_flow, prod(rPDParameters.ConstitutiveMatrix, plastic_flow)) +
        rPDParameters.Slope*inner_prod(rPDParameters.StressVector / g, plastic_flow));

    if (std::abs(denominator) > machine_tolerance)
        rPDParameters.PlasticConsistencyIncrement = rPDParameters.NonLinearIndicator / denominator;
    else
        rPDParameters.PlasticConsistencyIncrement = 0.0;

    if (rPDParameters.PlasticConsistencyIncrement < machine_tolerance) {
        BoundedMatrixType constitutive_matrix;
        CalculateElasticMatrix(constitutive_matrix, rValues);
        denominator = (inner_prod(plastic_flow, prod(constitutive_matrix, plastic_flow)) +
            rPDParameters.Slope*inner_prod(rPDParameters.StressVector / g, plastic_flow));
        rPDParameters.PlasticConsistencyIncrement = rPDParameters.NonLinearIndicator / denominator;
    }

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
    const auto& r_mat_properties = rValues.GetMaterialProperties();
    BoundedMatrixType stiffness_increment;

    bool is_converged = false;
    IndexType iteration = 0, max_iter = 1000;
    if (rPDParameters.PlasticDissipation < 0.99999) {
        while (is_converged == false && iteration <= max_iter) {
            // Evaluate the updated Threshold and slope
            CalculateThresholdAndSlope(rValues, rPDParameters);
            CalculateFlowVector(rValues, rPDParameters);
            CalculatePlasticConsistencyIncrement(rValues, rPDParameters);

            // Compute the compliance increment -> C dot
            CalculateComplianceMatrixIncrement(rValues, rPDParameters);
            noalias(rPDParameters.ComplianceMatrix) += rPDParameters.ComplianceMatrixIncrement;
            // CalculateConstitutiveMatrix(rValues, rPDParameters);

            // Compute the plastic strain increment
            CalculatePlasticStrainIncrement(rValues, rPDParameters);
            noalias(rPDParameters.PlasticStrain) += rPDParameters.PlasticStrainIncrement;

            // Correct the stress
            // noalias(rPDParameters.StressVector) = prod(rPDParameters.ConstitutiveMatrix,
            //     rPDParameters.StrainVector - rPDParameters.PlasticStrain);
            noalias(rPDParameters.StressVector) -= prod(rPDParameters.ConstitutiveMatrix,rPDParameters.PlasticStrainIncrement);
            noalias(stiffness_increment) = -prod(rPDParameters.ConstitutiveMatrix, Matrix(prod(rPDParameters.ComplianceMatrixIncrement, rPDParameters.ConstitutiveMatrix)));
            noalias(rPDParameters.StressVector) += prod(stiffness_increment, rPDParameters.StrainVector - rPDParameters.PlasticStrain);

            CalculateConstitutiveMatrix(rValues, rPDParameters);

            // Compute the non-linear dissipation performed
            CalculatePlasticDissipationIncrement(r_mat_properties, rPDParameters);
            CalculateDamageDissipationIncrement(r_mat_properties, rPDParameters);
            rPDParameters.DamageDissipation  += rPDParameters.DamageDissipationIncrement;
            rPDParameters.PlasticDissipation += rPDParameters.PlasticDissipationIncrement;

            rPDParameters.TotalDissipation   = (rPDParameters.PlasticDissipation + rPDParameters.DamageDissipation);
            rPDParameters.TotalDissipation = (rPDParameters.TotalDissipation > 0.99999) ? 0.99999 : rPDParameters.TotalDissipation;

            // updated uniaxial and threshold stress check
            TYieldSurfaceType::CalculateEquivalentStress(rPDParameters.StressVector,
                rPDParameters.StrainVector, rPDParameters.UniaxialStress, rValues);
            CalculateThresholdAndSlope(rValues, rPDParameters);
            rPDParameters.NonLinearIndicator = rPDParameters.UniaxialStress - rPDParameters.Threshold;

            if (rPDParameters.NonLinearIndicator <= std::abs(1.0e-4 * rPDParameters.Threshold)) {
                is_converged = true;
            } else {
                iteration++;
            }
        }
        if (iteration > max_iter) {
            KRATOS_ERROR << "Maximum number of iterations in plasticity loop reached..." << std::endl;
        }
    } else {
        noalias(rPDParameters.StressVector) = ZeroVector(VoigtSize);
    }

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

    const double young_modulus = r_material_properties[YOUNG_MODULUS];
    const double yield = r_material_properties[YIELD_STRESS];
    const double fracture_energy = r_material_properties[FRACTURE_ENERGY];

    const double hlim = 2.0 * young_modulus * fracture_energy / (std::pow(yield, 2));
    KRATOS_ERROR_IF(rPDParameters.CharacteristicLength > hlim) << "The Fracture Energy is to low: " <<
        rPDParameters.CharacteristicLength << std::endl;

}

/***********************************************************************************/
/***********************************************************************************/

/***********************************************************************************/
/***********************************************************************************/

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
    ConstitutiveLaw::Parameters& rValues
    )
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

template class AssociativePlasticDamageModel<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>;
// template class AssociativePlasticDamageModel<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>;

} // Namespace Kratos
