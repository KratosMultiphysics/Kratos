// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    DamApplication developers
//

// Characterization tests for the central refactoring hypothesis: can the
// DamApplication thermal constitutive law (ThermalLinearElastic3DLaw) be used
// directly with the StructuralMechanicsApplication small-displacement element
// (Kratos::SmallDisplacement) while reproducing the thermo-mechanical
// equilibrium response of the legacy DamApplication thermo-mechanical element
// (Kratos::SmallDisplacementThermoMechanicElement)?
//
// Three scenarios are run on a Hexahedra3D8:
//   T0 - no thermal increment (TEMPERATURE == reference temperature), with a
//        small non-zero displacement field;
//   T1 - uniform restrained temperature increment (zero displacements,
//        uniform TEMPERATURE = reference + delta_T);
//   T2 - non-uniform nodal temperature field (zero displacements, different
//        temperatures per node, consistent reference-temperature field).
//
// For each scenario two independent but numerically identical ModelParts are
// created (one per element implementation), each with its own Properties and
// its own ThermalLinearElastic3DLaw instance, and the following quantities are
// compared: local-system LHS, local-system RHS, independently calculated LHS
// and RHS, strain vector, Cauchy stress vector and PK2 stress vector at every
// integration point.
//
// ThermalLinearElastic3DLaw now implements the same infinitesimal thermo-elastic
// response for the PK2, Kirchhoff and Cauchy stress measures:
//     epsilon_th = alpha * (T - T_ref) * [1,1,1,0,0,0]
//     stress     = C * (epsilon - epsilon_th)
// where temperature and reference temperature are interpolated from the nodal
// TEMPERATURE and NODAL_REFERENCE_TEMPERATURE values with the shape functions
// supplied through ConstitutiveLaw::Parameters. For an infinitesimal-strain law
// the three measures coincide, so the finite-deformation 1/detF conversion of
// the inherited HyperElastic3DLaw Cauchy implementation is not applied. The law
// is stateless, so RequiresInitializeMaterialResponse()/RequiresFinalizeMaterialResponse()
// return false.
//
// Specialized thermal outputs (THERMAL_STRAIN_VECTOR, THERMAL_STRESS_VECTOR,
// MECHANICAL_STRESS_VECTOR) are POSTPONED to a dedicated task and are only
// probed here as non-asserting diagnostics: the previous characterization
// showed that the legacy element's THERMAL_STRAIN_VECTOR returns the mechanical
// strain (a legacy bug) and that the candidate element returns zero (outputs
// not implemented). They are not part of the equilibrium-compatibility
// acceptance criteria of this task.

// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

// Project includes
#include "dam_fast_suite.h"
#include "containers/model.h"
#include "includes/constitutive_law.h"
#include "includes/element.h"
#include "includes/kratos_components.h"
#include "includes/model_part.h"
#include "includes/variables.h"

// Application includes
#include "dam_application_variables.h"
#include "custom_constitutive/thermal_linear_elastic_3D_law.hpp"

namespace Kratos
{
namespace Testing
{

namespace
{

/// Comparison tolerances (same philosophy as the mechanical characterization).
constexpr double comparison_absolute_tolerance = 1.0e-12;
constexpr double comparison_relative_tolerance = 1.0e-10;
constexpr double machine_precision_allowance = 1.0e-15;

/// Material and loading data shared by both model parts.
constexpr double test_young_modulus = 2.0e7;
constexpr double test_poisson_ratio = 0.2;
constexpr double test_density = 2400.0;
constexpr double test_thermal_expansion = 1.0e-5;
constexpr double test_reference_temperature = 20.0;
constexpr double test_uniform_temperature_t1 = 45.0; // reference + 25
constexpr double test_delta_temperature_t1 = 25.0;
constexpr std::size_t test_number_of_hexa_nodes = 8;

/// All the element responses compared between the two implementations.
struct ThermalElementResponse
{
    Matrix lhs;
    Vector rhs;
    Matrix independent_lhs;
    Vector independent_rhs;
    std::vector<Vector> strain_vectors;
    std::vector<Vector> cauchy_stress_vectors;
    std::vector<Vector> pk2_stress_vectors;
    std::vector<Vector> thermal_strain_vectors;
    std::vector<Vector> thermal_stress_vectors;
    std::vector<Vector> mechanical_stress_vectors;
};

/// Creates one of the two numerically identical 3D8N model parts.
/// @param rElementName Registered name of the element to be created
/// @param rNodeIdOffset Offset applied to the node ids
ModelPart& CreateThermalComparisonModelPart(
    Model& rModel,
    const std::string& rModelPartName,
    const std::string& rElementName,
    const ModelPart::IndexType rNodeIdOffset)
{
    KRATOS_TRY;

    ModelPart& r_model_part = rModel.CreateModelPart(rModelPartName, 2);

    ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
    r_process_info[DOMAIN_SIZE] = 3;
    r_process_info[SPACE_DIMENSION] = 3;
    r_process_info[IS_RESTARTED] = false;

    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(VELOCITY);
    r_model_part.AddNodalSolutionStepVariable(ACCELERATION);
    r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    r_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    r_model_part.AddNodalSolutionStepVariable(NODAL_REFERENCE_TEMPERATURE);

    // Nodes on a scaled and translated copy of the prototype hexahedron so
    // that both elements are tested with exactly the same geometry.
    const Element& r_prototype = KratosComponents<Element>::Get("SmallDisplacementSolidElement3D8N");
    const auto& r_geometry = r_prototype.GetGeometry();
    Matrix local_coordinates;
    r_geometry.PointsLocalCoordinates(local_coordinates);

    const double geometry_scale = 2.5;
    array_1d<double, 3> geometry_offset;
    geometry_offset[0] = 0.75;
    geometry_offset[1] = 1.25;
    geometry_offset[2] = 0.5;

    const std::size_t number_of_nodes = r_geometry.PointsNumber();
    for (std::size_t i = 0; i < number_of_nodes; ++i) {
        const double x = geometry_scale * local_coordinates(i, 0) + geometry_offset[0];
        const double y = geometry_scale * local_coordinates(i, 1) + geometry_offset[1];
        const double z = geometry_scale * local_coordinates(i, 2) + geometry_offset[2];
        r_model_part.CreateNewNode(rNodeIdOffset + i + 1, x, y, z);
    }

    auto p_prop = r_model_part.CreateNewProperties(1);
    (*p_prop)[YOUNG_MODULUS] = test_young_modulus;
    (*p_prop)[POISSON_RATIO] = test_poisson_ratio;
    (*p_prop)[DENSITY] = test_density;
    (*p_prop)[THERMAL_EXPANSION] = test_thermal_expansion;
    // Independent constitutive law instance for this model part (each element
    // clones its own during Initialize).
    p_prop->SetValue(CONSTITUTIVE_LAW, ThermalLinearElastic3DLaw().Clone());

    std::vector<ModelPart::IndexType> element_nodes(number_of_nodes);
    for (std::size_t i = 0; i < number_of_nodes; ++i) {
        element_nodes[i] = rNodeIdOffset + i + 1;
    }
    r_model_part.CreateNewElement(rElementName, 1, element_nodes, p_prop);

    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(DISPLACEMENT_X);
        r_node.AddDof(DISPLACEMENT_Y);
        r_node.AddDof(DISPLACEMENT_Z);
        r_node.FastGetSolutionStepValue(NODAL_REFERENCE_TEMPERATURE) = test_reference_temperature;
        r_node.FastGetSolutionStepValue(TEMPERATURE) = test_reference_temperature;
    }

    return r_model_part;

    KRATOS_CATCH("");
}

/// Prescribes a scenario. The nodal VOLUME_ACCELERATION is identical in both
/// model parts so that the external-force contribution is symmetric.
void PrescribeScenario(
    ModelPart& rModelPart,
    const std::size_t rScenarioIndex,
    const bool rNonZeroDisplacement,
    const double rTemperature)
{
    KRATOS_TRY;

    if (rNonZeroDisplacement) {
        // u = A*X0 + t with a small strain magnitude.
        Matrix displacement_gradient = ZeroMatrix(3, 3);
        displacement_gradient(0, 0) = 2.0e-3;
        displacement_gradient(0, 1) = 3.0e-4;
        displacement_gradient(0, 2) = 1.0e-4;
        displacement_gradient(1, 0) = 3.0e-4;
        displacement_gradient(1, 1) = -1.0e-3;
        displacement_gradient(1, 2) = 2.0e-4;
        displacement_gradient(2, 0) = 1.0e-4;
        displacement_gradient(2, 1) = 2.0e-4;
        displacement_gradient(2, 2) = 5.0e-4;
        array_1d<double, 3> translation;
        translation[0] = 1.0e-2;
        translation[1] = -2.0e-2;
        translation[2] = 5.0e-3;

        for (auto& r_node : rModelPart.Nodes()) {
            array_1d<double, 3>& r_displacement = r_node.FastGetSolutionStepValue(DISPLACEMENT);
            noalias(r_displacement) = prod(displacement_gradient, r_node.GetInitialPosition());
            noalias(r_displacement) += translation;
            // Updated coordinates so that the reference configuration (current
            // position - total displacement) is the initial one.
            r_node.X() = r_node.X0() + r_displacement[0];
            r_node.Y() = r_node.Y0() + r_displacement[1];
            r_node.Z() = r_node.Z0() + r_displacement[2];
        }
    } else {
        for (auto& r_node : rModelPart.Nodes()) {
            r_node.FastGetSolutionStepValue(DISPLACEMENT) = ZeroVector(3);
            r_node.X() = r_node.X0();
            r_node.Y() = r_node.Y0();
            r_node.Z() = r_node.Z0();
        }
    }

    std::size_t node_index = 0;
    for (auto& r_node : rModelPart.Nodes()) {
        // T0: TEMPERATURE == reference (rTemperature == test_reference_temperature).
        // T1: uniform TEMPERATURE == reference + delta_T.
        // T2: non-uniform nodal temperature field.
        double nodal_temperature = rTemperature;
        if (rScenarioIndex == 2) {
            nodal_temperature = test_reference_temperature + 5.0 * static_cast<double>(node_index);
        }
        r_node.FastGetSolutionStepValue(TEMPERATURE) = nodal_temperature;
        ++node_index;
    }

    // Identical uniform body force in both models.
    for (auto& r_node : rModelPart.Nodes()) {
        array_1d<double, 3>& r_volume_acceleration = r_node.FastGetSolutionStepValue(VOLUME_ACCELERATION);
        r_volume_acceleration[0] = 0.0;
        r_volume_acceleration[1] = -9.81;
        r_volume_acceleration[2] = 0.0;
    }

    KRATOS_CATCH("");
}

/// Runs the standard element lifecycle up to the point where responses can be
/// calculated. After ThermalLinearElastic3DLaw reported that no material-response
/// initialization/finalization is required, every lifecycle step must succeed.
void InitializeThermalElementLifecycle(ModelPart& rModelPart, const std::string& r_element_label)
{
    KRATOS_TRY;

    auto p_element = rModelPart.pGetElement(1);
    const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

    KRATOS_EXPECT_EQ(p_element->Check(r_process_info), 0);
    p_element->Initialize(r_process_info);
    p_element->InitializeSolutionStep(r_process_info);
    std::cout << "[CHARACTERIZATION] lifecycle " << r_element_label
              << ": InitializeSolutionStep OK" << std::endl;
    p_element->InitializeNonLinearIteration(r_process_info);
    std::cout << "[CHARACTERIZATION] lifecycle " << r_element_label
              << ": InitializeNonLinearIteration OK" << std::endl;

    KRATOS_CATCH("");
}

/// Calculates all the compared responses of the element of the model part.
void CalculateThermalElementResponse(ModelPart& rModelPart, ThermalElementResponse& rResponse)
{
    KRATOS_TRY;

    auto p_element = rModelPart.pGetElement(1);
    const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

    p_element->CalculateLocalSystem(rResponse.lhs, rResponse.rhs, r_process_info);
    p_element->CalculateLeftHandSide(rResponse.independent_lhs, r_process_info);
    p_element->CalculateRightHandSide(rResponse.independent_rhs, r_process_info);
    p_element->CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR, rResponse.strain_vectors, r_process_info);
    p_element->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, rResponse.cauchy_stress_vectors, r_process_info);
    p_element->CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, rResponse.pk2_stress_vectors, r_process_info);
    p_element->CalculateOnIntegrationPoints(THERMAL_STRAIN_VECTOR, rResponse.thermal_strain_vectors, r_process_info);
    p_element->CalculateOnIntegrationPoints(THERMAL_STRESS_VECTOR, rResponse.thermal_stress_vectors, r_process_info);
    p_element->CalculateOnIntegrationPoints(MECHANICAL_STRESS_VECTOR, rResponse.mechanical_stress_vectors, r_process_info);

    KRATOS_CATCH("");
}

/// Creates the two independent but numerically identical model parts (one per
/// element implementation), prescribes the scenario and calculates their
/// responses.
void CalculateThermalComparisonResponses(
    const std::size_t rScenarioIndex,
    ThermalElementResponse& rLegacyResponse,
    ThermalElementResponse& rCandidateResponse)
{
    KRATOS_TRY;

    const bool non_zero_displacement = (rScenarioIndex == 0);
    const double uniform_temperature =
        (rScenarioIndex == 0) ? test_reference_temperature : test_uniform_temperature_t1;

    Model model;
    ModelPart& r_legacy_model_part = CreateThermalComparisonModelPart(
        model, "LegacyDamThermoMechanic", "SmallDisplacementThermoMechanicElement3D8N", 0);
    ModelPart& r_candidate_model_part = CreateThermalComparisonModelPart(
        model, "CandidateStructuralSmallDisplacement", "SmallDisplacementElement3D8N", 100);

    PrescribeScenario(r_legacy_model_part, rScenarioIndex, non_zero_displacement, uniform_temperature);
    PrescribeScenario(r_candidate_model_part, rScenarioIndex, non_zero_displacement, uniform_temperature);

    InitializeThermalElementLifecycle(r_legacy_model_part, "legacy");
    InitializeThermalElementLifecycle(r_candidate_model_part, "candidate");

    CalculateThermalElementResponse(r_legacy_model_part, rLegacyResponse);
    CalculateThermalElementResponse(r_candidate_model_part, rCandidateResponse);

    KRATOS_CATCH("");
}

/// Tolerance associated with a reference value (same philosophy as the
/// mechanical characterization tests).
double ComponentTolerance(const double rReferenceValue, const double rReferenceScale)
{
    double tolerance = std::max(comparison_absolute_tolerance,
                                comparison_relative_tolerance * std::abs(rReferenceValue));
    if (rReferenceValue == 0.0) {
        tolerance = std::max(tolerance, machine_precision_allowance * rReferenceScale);
    }
    return tolerance;
}

double MaxAbsEntry(const Matrix& rMatrix)
{
    double max_abs_entry = 0.0;
    for (std::size_t i = 0; i < rMatrix.size1(); ++i) {
        for (std::size_t j = 0; j < rMatrix.size2(); ++j) {
            max_abs_entry = std::max(max_abs_entry, std::abs(rMatrix(i, j)));
        }
    }
    return max_abs_entry;
}

double MaxAbsEntry(const Vector& rVector)
{
    double max_abs_entry = 0.0;
    for (std::size_t i = 0; i < rVector.size(); ++i) {
        max_abs_entry = std::max(max_abs_entry, std::abs(rVector(i)));
    }
    return max_abs_entry;
}

/// Metrics accumulated over the components of one comparison. Purely
/// diagnostic; the pass/fail decision is the component-wise EXPECT check.
struct ComparisonMetrics
{
    bool pass = true;
    bool comparable = true;
    std::size_t component_count = 0;
    std::size_t failed_component_count = 0;
    double max_absolute_difference = 0.0;
    double max_relative_difference = 0.0;
    double max_tolerance_used = 0.0;
    double reference_scale = 0.0;
};

void AccumulateComponentMetrics(
    const double rComputedValue,
    const double rReferenceValue,
    const double rTolerance,
    ComparisonMetrics& rMetrics)
{
    const double absolute_difference = std::abs(rComputedValue - rReferenceValue);
    rMetrics.max_absolute_difference = std::max(rMetrics.max_absolute_difference, absolute_difference);
    rMetrics.max_tolerance_used = std::max(rMetrics.max_tolerance_used, rTolerance);
    rMetrics.reference_scale = std::max(rMetrics.reference_scale, std::abs(rReferenceValue));
    if (rReferenceValue != 0.0) {
        rMetrics.max_relative_difference =
            std::max(rMetrics.max_relative_difference, absolute_difference / std::abs(rReferenceValue));
    }
    if (absolute_difference > rTolerance) {
        ++rMetrics.failed_component_count;
        rMetrics.pass = false;
    }
    ++rMetrics.component_count;
}

void PrintComparisonMetrics(const std::string& rWhat, const ComparisonMetrics& rMetrics)
{
    if (!rMetrics.comparable) {
        std::cout << "[CHARACTERIZATION] " << rWhat << ": NOT COMPARABLE (output unavailable)"
                  << std::endl;
        return;
    }
    std::cout << "[CHARACTERIZATION] " << rWhat << ": "
              << (rMetrics.pass ? "PASS" : "FAIL")
              << " | max_abs_diff=" << rMetrics.max_absolute_difference
              << " | max_rel_diff=" << rMetrics.max_relative_difference
              << " | tolerance=" << rMetrics.max_tolerance_used
              << " | scale=" << rMetrics.reference_scale
              << " | failed_components=" << rMetrics.failed_component_count
              << "/" << rMetrics.component_count
              << std::endl;
}

void ExpectVectorComponentsNear(
    const Vector& rComputed,
    const Vector& rReference,
    const std::string& rWhat)
{
    ASSERT_EQ(rComputed.size(), rReference.size()) << "Size mismatch comparing " << rWhat;
    const double reference_scale = MaxAbsEntry(rReference);
    ComparisonMetrics metrics;
    for (std::size_t i = 0; i < rReference.size(); ++i) {
        AccumulateComponentMetrics(
            rComputed(i), rReference(i), ComponentTolerance(rReference(i), reference_scale), metrics);
    }
    PrintComparisonMetrics(rWhat, metrics);
    for (std::size_t i = 0; i < rReference.size(); ++i) {
        KRATOS_EXPECT_NEAR(rComputed(i), rReference(i), ComponentTolerance(rReference(i), reference_scale));
    }
}

void ExpectMatrixComponentsNear(
    const Matrix& rComputed,
    const Matrix& rReference,
    const std::string& rWhat)
{
    ASSERT_EQ(rComputed.size1(), rReference.size1()) << "Row size mismatch comparing " << rWhat;
    ASSERT_EQ(rComputed.size2(), rReference.size2()) << "Column size mismatch comparing " << rWhat;
    const double reference_scale = MaxAbsEntry(rReference);
    ComparisonMetrics metrics;
    for (std::size_t i = 0; i < rReference.size1(); ++i) {
        for (std::size_t j = 0; j < rReference.size2(); ++j) {
            AccumulateComponentMetrics(
                rComputed(i, j), rReference(i, j), ComponentTolerance(rReference(i, j), reference_scale), metrics);
        }
    }
    PrintComparisonMetrics(rWhat, metrics);
    for (std::size_t i = 0; i < rReference.size1(); ++i) {
        for (std::size_t j = 0; j < rReference.size2(); ++j) {
            KRATOS_EXPECT_NEAR(rComputed(i, j), rReference(i, j), ComponentTolerance(rReference(i, j), reference_scale));
        }
    }
}

void ExpectIntegrationPointVectorsNear(
    const std::vector<Vector>& rComputed,
    const std::vector<Vector>& rReference,
    const std::string& rWhat)
{
    ASSERT_EQ(rComputed.size(), rReference.size()) << "Integration point count mismatch comparing " << rWhat;
    for (std::size_t point_number = 0; point_number < rReference.size(); ++point_number) {
        ExpectVectorComponentsNear(
            rComputed[point_number], rReference[point_number],
            rWhat + " at integration point " + std::to_string(point_number));
    }
}

} // namespace

//************************************************************************************
// T0 - No thermal increment, small non-zero displacement field
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementThermalComparisonT0_LocalSystemLHS, KratosDamFastSuite)
{
    ThermalElementResponse legacy, candidate;
    CalculateThermalComparisonResponses(0, legacy, candidate);
    ExpectMatrixComponentsNear(legacy.lhs, candidate.lhs, "T0 local system LHS");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementThermalComparisonT0_LocalSystemRHS, KratosDamFastSuite)
{
    ThermalElementResponse legacy, candidate;
    CalculateThermalComparisonResponses(0, legacy, candidate);
    ExpectVectorComponentsNear(legacy.rhs, candidate.rhs, "T0 local system RHS");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementThermalComparisonT0_IndependentLHS, KratosDamFastSuite)
{
    ThermalElementResponse legacy, candidate;
    CalculateThermalComparisonResponses(0, legacy, candidate);
    ExpectMatrixComponentsNear(legacy.independent_lhs, candidate.independent_lhs, "T0 independent LHS");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementThermalComparisonT0_IndependentRHS, KratosDamFastSuite)
{
    ThermalElementResponse legacy, candidate;
    CalculateThermalComparisonResponses(0, legacy, candidate);
    ExpectVectorComponentsNear(legacy.independent_rhs, candidate.independent_rhs, "T0 independent RHS");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementThermalComparisonT0_StrainVector, KratosDamFastSuite)
{
    ThermalElementResponse legacy, candidate;
    CalculateThermalComparisonResponses(0, legacy, candidate);
    ExpectIntegrationPointVectorsNear(legacy.strain_vectors, candidate.strain_vectors, "T0 strain vector");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementThermalComparisonT0_CauchyStressVector, KratosDamFastSuite)
{
    ThermalElementResponse legacy, candidate;
    CalculateThermalComparisonResponses(0, legacy, candidate);
    ExpectIntegrationPointVectorsNear(legacy.cauchy_stress_vectors, candidate.cauchy_stress_vectors, "T0 Cauchy stress vector");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementThermalComparisonT0_PK2StressVector, KratosDamFastSuite)
{
    // The legacy thermo-mechanical element does not implement a PK2 stress
    // output, so the candidate PK2 stress is checked against the physically
    // expected thermo-elastic response: the legacy Cauchy stress (total stress
    // C*(epsilon - epsilon_th) with det(F) == 1) and, for an infinitesimal
    // strain law, the candidate's own Cauchy stress (PK2 == Kirchhoff == Cauchy).
    ThermalElementResponse legacy, candidate;
    CalculateThermalComparisonResponses(0, legacy, candidate);
    ExpectIntegrationPointVectorsNear(candidate.pk2_stress_vectors, legacy.cauchy_stress_vectors, "T0 candidate PK2 vs legacy Cauchy stress vector");
    ExpectIntegrationPointVectorsNear(candidate.pk2_stress_vectors, candidate.cauchy_stress_vectors, "T0 candidate PK2 vs candidate Cauchy stress vector");
}

//************************************************************************************
// T1 - Uniform restrained temperature increment, zero displacements
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementThermalComparisonT1_LocalSystemLHS, KratosDamFastSuite)
{
    ThermalElementResponse legacy, candidate;
    CalculateThermalComparisonResponses(1, legacy, candidate);
    ExpectMatrixComponentsNear(legacy.lhs, candidate.lhs, "T1 local system LHS");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementThermalComparisonT1_LocalSystemRHS, KratosDamFastSuite)
{
    ThermalElementResponse legacy, candidate;
    CalculateThermalComparisonResponses(1, legacy, candidate);
    ExpectVectorComponentsNear(legacy.rhs, candidate.rhs, "T1 local system RHS");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementThermalComparisonT1_IndependentLHS, KratosDamFastSuite)
{
    ThermalElementResponse legacy, candidate;
    CalculateThermalComparisonResponses(1, legacy, candidate);
    ExpectMatrixComponentsNear(legacy.independent_lhs, candidate.independent_lhs, "T1 independent LHS");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementThermalComparisonT1_IndependentRHS, KratosDamFastSuite)
{
    ThermalElementResponse legacy, candidate;
    CalculateThermalComparisonResponses(1, legacy, candidate);
    ExpectVectorComponentsNear(legacy.independent_rhs, candidate.independent_rhs, "T1 independent RHS");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementThermalComparisonT1_StrainVector, KratosDamFastSuite)
{
    ThermalElementResponse legacy, candidate;
    CalculateThermalComparisonResponses(1, legacy, candidate);
    ExpectIntegrationPointVectorsNear(legacy.strain_vectors, candidate.strain_vectors, "T1 strain vector");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementThermalComparisonT1_CauchyStressVector, KratosDamFastSuite)
{
    ThermalElementResponse legacy, candidate;
    CalculateThermalComparisonResponses(1, legacy, candidate);
    ExpectIntegrationPointVectorsNear(legacy.cauchy_stress_vectors, candidate.cauchy_stress_vectors, "T1 Cauchy stress vector");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementThermalComparisonT1_PK2StressVector, KratosDamFastSuite)
{
    ThermalElementResponse legacy, candidate;
    CalculateThermalComparisonResponses(1, legacy, candidate);
    ExpectIntegrationPointVectorsNear(candidate.pk2_stress_vectors, legacy.cauchy_stress_vectors, "T1 candidate PK2 vs legacy Cauchy stress vector");
    ExpectIntegrationPointVectorsNear(candidate.pk2_stress_vectors, candidate.cauchy_stress_vectors, "T1 candidate PK2 vs candidate Cauchy stress vector");
}

//************************************************************************************
// T2 - Non-uniform nodal temperature field, zero displacements
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementThermalComparisonT2_LocalSystemLHS, KratosDamFastSuite)
{
    ThermalElementResponse legacy, candidate;
    CalculateThermalComparisonResponses(2, legacy, candidate);
    ExpectMatrixComponentsNear(legacy.lhs, candidate.lhs, "T2 local system LHS");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementThermalComparisonT2_LocalSystemRHS, KratosDamFastSuite)
{
    ThermalElementResponse legacy, candidate;
    CalculateThermalComparisonResponses(2, legacy, candidate);
    ExpectVectorComponentsNear(legacy.rhs, candidate.rhs, "T2 local system RHS");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementThermalComparisonT2_IndependentLHS, KratosDamFastSuite)
{
    ThermalElementResponse legacy, candidate;
    CalculateThermalComparisonResponses(2, legacy, candidate);
    ExpectMatrixComponentsNear(legacy.independent_lhs, candidate.independent_lhs, "T2 independent LHS");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementThermalComparisonT2_IndependentRHS, KratosDamFastSuite)
{
    ThermalElementResponse legacy, candidate;
    CalculateThermalComparisonResponses(2, legacy, candidate);
    ExpectVectorComponentsNear(legacy.independent_rhs, candidate.independent_rhs, "T2 independent RHS");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementThermalComparisonT2_StrainVector, KratosDamFastSuite)
{
    ThermalElementResponse legacy, candidate;
    CalculateThermalComparisonResponses(2, legacy, candidate);
    ExpectIntegrationPointVectorsNear(legacy.strain_vectors, candidate.strain_vectors, "T2 strain vector");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementThermalComparisonT2_CauchyStressVector, KratosDamFastSuite)
{
    ThermalElementResponse legacy, candidate;
    CalculateThermalComparisonResponses(2, legacy, candidate);
    ExpectIntegrationPointVectorsNear(legacy.cauchy_stress_vectors, candidate.cauchy_stress_vectors, "T2 Cauchy stress vector");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementThermalComparisonT2_PK2StressVector, KratosDamFastSuite)
{
    ThermalElementResponse legacy, candidate;
    CalculateThermalComparisonResponses(2, legacy, candidate);
    ExpectIntegrationPointVectorsNear(candidate.pk2_stress_vectors, legacy.cauchy_stress_vectors, "T2 candidate PK2 vs legacy Cauchy stress vector");
    ExpectIntegrationPointVectorsNear(candidate.pk2_stress_vectors, candidate.cauchy_stress_vectors, "T2 candidate PK2 vs candidate Cauchy stress vector");
}

//************************************************************************************
// Specialized thermal outputs (POSTPONED to a dedicated task)
//************************************************************************************
// THERMAL_STRAIN_VECTOR, THERMAL_STRESS_VECTOR and MECHANICAL_STRESS_VECTOR are
// intentionally NOT part of this equilibrium-compatibility task. The previous
// characterization showed that:
//   - the legacy element's THERMAL_STRAIN_VECTOR returns the mechanical strain
//     (the legacy element sets THERMAL_RESPONSE_ONLY without VOLUMETRIC_TENSOR_ONLY,
//     so the constitutive law never fills the thermal strain - a legacy bug);
//   - THERMAL_STRESS_VECTOR and MECHANICAL_STRESS_VECTOR are returned by the
//     legacy element but are not implemented by the StructuralMechanics element
//     (its standard output mechanism returns zero for them).
// They are probed below as a NON-ASSERTING diagnostic only, so that the known
// specialized-output differences do not make this equilibrium task fail.

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementThermalSpecializedOutputs_PostponedDiagnostic, KratosDamFastSuite)
{
    for (std::size_t scenario = 0; scenario < 3; ++scenario) {
        ThermalElementResponse legacy, candidate;
        CalculateThermalComparisonResponses(scenario, legacy, candidate);
        std::cout << "[CHARACTERIZATION] specialized outputs scenario "
                  << scenario << " (postponed diagnostic):" << std::endl;
        std::cout << "[CHARACTERIZATION]   THERMAL_STRAIN_VECTOR: legacy_size="
                  << (legacy.thermal_strain_vectors.empty() ? 0 : legacy.thermal_strain_vectors[0].size())
                  << " candidate_size="
                  << (candidate.thermal_strain_vectors.empty() ? 0 : candidate.thermal_strain_vectors[0].size())
                  << std::endl;
        std::cout << "[CHARACTERIZATION]   THERMAL_STRESS_VECTOR: legacy_size="
                  << (legacy.thermal_stress_vectors.empty() ? 0 : legacy.thermal_stress_vectors[0].size())
                  << " candidate_size="
                  << (candidate.thermal_stress_vectors.empty() ? 0 : candidate.thermal_stress_vectors[0].size())
                  << std::endl;
        std::cout << "[CHARACTERIZATION]   MECHANICAL_STRESS_VECTOR: legacy_size="
                  << (legacy.mechanical_stress_vectors.empty() ? 0 : legacy.mechanical_stress_vectors[0].size())
                  << " candidate_size="
                  << (candidate.mechanical_stress_vectors.empty() ? 0 : candidate.mechanical_stress_vectors[0].size())
                  << std::endl;
    }
}

//************************************************************************************
// Material-response lifecycle
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementThermalLifecycle_CandidateFullLifecycle, KratosDamFastSuite)
{
    // The StructuralMechanics element with ThermalLinearElastic3DLaw must go
    // through the complete solution-step lifecycle (including FinalizeSolutionStep)
    // without any constitutive-law lifecycle error. RequiresInitializeMaterialResponse()
    // and RequiresFinalizeMaterialResponse() now return false, so the element
    // skips the material-response initialization/finalization calls that the
    // stateless thermal law does not implement.
    Model model;
    ModelPart& r_candidate = CreateThermalComparisonModelPart(
        model, "CandidateLifecycle", "SmallDisplacementElement3D8N", 0);
    PrescribeScenario(r_candidate, 1, false, test_uniform_temperature_t1);
    auto p_element = r_candidate.pGetElement(1);
    const ProcessInfo& r_process_info = r_candidate.GetProcessInfo();

    KRATOS_EXPECT_EQ(p_element->Check(r_process_info), 0);
    p_element->Initialize(r_process_info);
    p_element->InitializeSolutionStep(r_process_info);
    p_element->FinalizeSolutionStep(r_process_info);
    std::cout << "[CHARACTERIZATION] lifecycle candidate: full solution-step lifecycle OK"
              << std::endl;
}

} // namespace Testing
} // namespace Kratos
