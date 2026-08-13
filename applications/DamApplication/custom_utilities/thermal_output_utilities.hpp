//
//   Project Name:                  KratosDamApplication $
//   Last Modified by:    $Author:    DamApplication developers $
//   Date:                $Date:                           $
//   Revision:            $Revision:                 1.0 $
//

// Small Dam-only helpers shared by the thermal constitutive families when they
// provide the specialized thermo-mechanical outputs
// (THERMAL_STRAIN/STRESS_VECTOR/TENSOR, MECHANICAL_STRESS_VECTOR/TENSOR)
// through the parameter-aware CalculateValue path. No base hierarchy is
// introduced: these are free inline functions on top of the existing
// constitutive machinery.

#if !defined (KRATOS_THERMAL_OUTPUT_UTILITIES_H_INCLUDED)
#define  KRATOS_THERMAL_OUTPUT_UTILITIES_H_INCLUDED

// Application includes
#include "custom_constitutive/continuum_laws/custom_flow_rules/flow_rule.hpp"
#include "custom_constitutive/continuum_laws/custom_yield_criteria/yield_criterion.hpp"
#include "utilities/math_utils.h"

namespace Kratos
{

/// Helpers for the specialized thermo-mechanical outputs.
namespace ThermalOutputUtilities
{

/// Voigt/tensor dimension for the given strain/voigt size (6 -> 3D, else 2D).
inline std::size_t GetDimension(const std::size_t rVoigtSize)
{
    return rVoigtSize == 6 ? 3 : 2;
}

/// Assigns a strain tensor output from its Voigt vector (single source: the
/// vector overload of CalculateValue).
inline void AssignStrainTensor(Matrix& rValue, const Vector& rVoigtVector)
{
    const std::size_t dimension = GetDimension(rVoigtVector.size());
    if (rValue.size1() != dimension || rValue.size2() != dimension)
        rValue.resize(dimension, dimension, false);
    noalias(rValue) = MathUtils<double>::StrainVectorToTensor(rVoigtVector);
}

/// Assigns a stress tensor output from its Voigt vector (single source: the
/// vector overload of CalculateValue).
inline void AssignStressTensor(Matrix& rValue, const Vector& rVoigtVector)
{
    const std::size_t dimension = GetDimension(rVoigtVector.size());
    if (rValue.size1() != dimension || rValue.size2() != dimension)
        rValue.resize(dimension, dimension, false);
    noalias(rValue) = MathUtils<double>::StressVectorToTensor(rVoigtVector);
}

/// Assembles the three specialized outputs from the total strain, the thermal
/// strain, the constitutive matrix and the current damage factor (1 - d):
///   THERMAL_STRAIN_VECTOR    = epsilon_th
///   THERMAL_STRESS_VECTOR    = (1-d) * C * epsilon_th
///   MECHANICAL_STRESS_VECTOR = (1-d) * C * epsilon_total
/// so that the total stress satisfies
///   stress == MECHANICAL_STRESS_VECTOR - THERMAL_STRESS_VECTOR.
inline void AssembleOutputs(
    Vector& rThermalStrain,
    Vector& rThermalStress,
    Vector& rMechanicalStress,
    const Vector& rTotalStrain,
    const Vector& rThermalStrainSource,
    const Matrix& rConstitutiveMatrix,
    const double rDamageFactor)
{
    if (rThermalStrain.size() != rThermalStrainSource.size())
        rThermalStrain.resize(rThermalStrainSource.size(), false);
    noalias(rThermalStrain) = rThermalStrainSource;

    if (rThermalStress.size() != rTotalStrain.size())
        rThermalStress.resize(rTotalStrain.size(), false);
    noalias(rThermalStress) = rDamageFactor * prod(rConstitutiveMatrix, rThermalStrainSource);

    if (rMechanicalStress.size() != rTotalStrain.size())
        rMechanicalStress.resize(rTotalStrain.size(), false);
    noalias(rMechanicalStress) = rDamageFactor * prod(rConstitutiveMatrix, rTotalStrain);
}

/// Evaluates the CURRENT damage factor (1 - d) from the committed maximum
/// equivalent strain, exactly as the trial return-mapping does during the
/// current constitutive response. The evaluation is read-only: no flow-rule
/// state (committed history, LOCAL/NONLOCAL driving quantities) is modified.
inline double CalculateCurrentDamageFactor(
    FlowRule& rFlowRule,
    YieldCriterion& rYieldCriterion,
    const double rCharacteristicSize)
{
    const FlowRule::InternalVariables& r_internal_variables = rFlowRule.GetInternalVariables();

    YieldCriterion::Parameters yield_criterion_parameters;
    yield_criterion_parameters.SetCharacteristicSize(rCharacteristicSize);
    yield_criterion_parameters.SetDeltaGamma(r_internal_variables.EquivalentPlasticStrain);

    double damage;
    rYieldCriterion.CalculateStateFunction(damage, yield_criterion_parameters);

    return 1.0 - damage;
}

} // namespace ThermalOutputUtilities

} // namespace Kratos

#endif // KRATOS_THERMAL_OUTPUT_UTILITIES_H_INCLUDED defined
