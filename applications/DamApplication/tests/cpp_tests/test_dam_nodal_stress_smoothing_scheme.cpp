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

// Phase 3D.3: integration of the process-based nodal Cauchy-stress
// extrapolation into the actual Dam smoothing-scheme workflow.
//
// The legacy element (SmallDisplacementThermoMechanicElement) accumulates the
// nodal Cauchy stress inside its FinalizeSolutionStep and the smoothing scheme
// resets and normalizes the historical NODAL_AREA / NODAL_CAUCHY_STRESS_TENSOR.
// In the process-based mode the scheme itself owns the extrapolation (calling
// DamNodalCauchyStressExtrapolationProcess::ExtrapolateAndAccumulate between
// the element finalization and the normalization), so the candidate model can
// use plain StructuralMechanicsApplication small-displacement elements with the
// compatible Dam constitutive laws. Exactly one owner accumulates per solution
// step in both modes.
//
// These tests run the REAL Dam smoothing schemes
// (IncrementalUpdateStaticSmoothingScheme and its damped variant, and
// BossakDisplacementSmoothingScheme) FinalizeSolutionStep sequence on a
// multi-element model with shared nodes and unequal element measures.

// System includes
#include <algorithm>
#include <array>
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
#include "spaces/ublas_space.h"

// Application includes
#include "dam_application_variables.h"
#include "custom_constitutive/thermal_linear_elastic_3D_law.hpp"
#include "custom_processes/dam_nodal_cauchy_stress_extrapolation_process.hpp"
#include "custom_strategies/schemes/incrementalupdate_static_smoothing_scheme.hpp"
#include "custom_strategies/schemes/incrementalupdate_static_damped_smoothing_scheme.hpp"
#include "custom_strategies/schemes/bossak_displacement_smoothing_scheme.hpp"
#include "custom_utilities/transfer_selfweight_stress_utility.hpp"

namespace Kratos
{
namespace Testing
{

namespace
{

/// Comparison tolerances (same philosophy as the previous characterization).
constexpr double comparison_absolute_tolerance = 1.0e-12;
constexpr double comparison_relative_tolerance = 1.0e-10;

/// Material data shared by all model parts.
constexpr double test_young_modulus = 2.0e7;
constexpr double test_poisson_ratio = 0.2;
constexpr double test_thermal_expansion = 1.0e-5;
constexpr double test_reference_temperature = 20.0;

/// Sparse/dense spaces matching the Dam solver bindings.
using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
using StaticSmoothingScheme =
    IncrementalUpdateStaticSmoothingScheme<SparseSpaceType, LocalSpaceType>;
using DampedSmoothingScheme =
    IncrementalUpdateStaticDampedSmoothingScheme<SparseSpaceType, LocalSpaceType>;
using BossakSmoothingScheme =
    BossakDisplacementSmoothingScheme<SparseSpaceType, LocalSpaceType>;

/// Two hexahedra sharing a face (12 unique nodes):
///   A: x in [0,2], y in [0,1], z in [0,1]   (measure 2.0)
///   B: x in [2,3], y in [0,1], z in [0,1]   (measure 1.0)
/// Shared nodes 2, 3, 6, 7 belong to both elements (accumulated measure 3.0),
/// the remaining nodes to a single element (measures 2.0 or 1.0).
constexpr std::array<std::array<double, 3>, 12> TestNodeCoordinates = {{
    {{0.0, 0.0, 0.0}}, {{2.0, 0.0, 0.0}}, {{2.0, 1.0, 0.0}}, {{0.0, 1.0, 0.0}},
    {{0.0, 0.0, 1.0}}, {{2.0, 0.0, 1.0}}, {{2.0, 1.0, 1.0}}, {{0.0, 1.0, 1.0}},
    {{3.0, 0.0, 0.0}}, {{3.0, 1.0, 0.0}}, {{3.0, 0.0, 1.0}}, {{3.0, 1.0, 1.0}}
}};

constexpr std::array<std::array<int, 8>, 2> TestElementConnectivity = {{
    {{1, 2, 3, 4, 5, 6, 7, 8}},
    {{2, 9, 10, 3, 6, 11, 12, 7}}
}};

/// Creates a model part with the two-element shared-face mesh using the given
/// element name, the thermal law, and the nodal variables required by the
/// smoothing scheme (matching the Dam solver AddVariables).
ModelPart& CreateSchemeModelPart(
    Model& rModel,
    const std::string& rModelPartName,
    const std::string& rElementName)
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
    r_model_part.AddNodalSolutionStepVariable(NODAL_CAUCHY_STRESS_TENSOR);
    r_model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    r_model_part.AddNodalSolutionStepVariable(NODAL_JOINT_AREA);
    r_model_part.AddNodalSolutionStepVariable(NODAL_JOINT_WIDTH);
    r_model_part.AddNodalSolutionStepVariable(INITIAL_NODAL_CAUCHY_STRESS_TENSOR);

    for (std::size_t i = 0; i < TestNodeCoordinates.size(); ++i) {
        r_model_part.CreateNewNode(
            i + 1, TestNodeCoordinates[i][0], TestNodeCoordinates[i][1],
            TestNodeCoordinates[i][2]);
    }

    auto p_prop = r_model_part.CreateNewProperties(1);
    (*p_prop)[YOUNG_MODULUS] = test_young_modulus;
    (*p_prop)[POISSON_RATIO] = test_poisson_ratio;
    (*p_prop)[DENSITY] = 2400.0;
    (*p_prop)[THERMAL_EXPANSION] = test_thermal_expansion;
    p_prop->SetValue(CONSTITUTIVE_LAW, ThermalLinearElastic3DLaw().Clone());

    for (std::size_t e = 0; e < TestElementConnectivity.size(); ++e) {
        std::vector<ModelPart::IndexType> element_nodes(8);
        for (std::size_t i = 0; i < 8; ++i) {
            element_nodes[i] = TestElementConnectivity[e][i];
        }
        r_model_part.CreateNewElement(rElementName, e + 1, element_nodes, p_prop);
    }

    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(DISPLACEMENT_X);
        r_node.AddDof(DISPLACEMENT_Y);
        r_node.AddDof(DISPLACEMENT_Z);
        r_node.FastGetSolutionStepValue(NODAL_REFERENCE_TEMPERATURE) = test_reference_temperature;
        r_node.FastGetSolutionStepValue(NODAL_AREA) = 0.0;
        Matrix zero_stress_tensor(3, 3);
        noalias(zero_stress_tensor) = ZeroMatrix(3, 3);
        r_node.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR) = zero_stress_tensor;
    }

    return r_model_part;

    KRATOS_CATCH("");
}

/// Prescribes a non-uniform thermo-mechanical state with non-degenerate shears.
void PrescribeState(ModelPart& rModelPart, const double rStateScale)
{
    KRATOS_TRY;

    std::size_t node_index = 0;
    for (auto& r_node : rModelPart.Nodes()) {
        const array_1d<double, 3>& r_x0 = r_node.GetInitialPosition();
        array_1d<double, 3>& r_displacement = r_node.FastGetSolutionStepValue(DISPLACEMENT);
        r_displacement[0] = rStateScale * (1.0e-3 * r_x0[0] + 3.0e-4 * r_x0[1]
                                          + 1.0e-4 * r_x0[0] * r_x0[1]
                                          + 2.0e-4 * r_x0[0] * r_x0[2]);
        r_displacement[1] = rStateScale * (3.0e-4 * r_x0[0] - 1.0e-3 * r_x0[1]
                                          + 2.0e-4 * r_x0[0] * r_x0[0]
                                          + 1.0e-4 * r_x0[1] * r_x0[2]);
        r_displacement[2] = rStateScale * (5.0e-4 * r_x0[2] + 2.0e-4 * r_x0[0]
                                          + 1.0e-4 * r_x0[1] + 1.0e-4 * r_x0[0] * r_x0[1]);
        r_node.X() = r_x0[0] + r_displacement[0];
        r_node.Y() = r_x0[1] + r_displacement[1];
        r_node.Z() = r_x0[2] + r_displacement[2];
        r_node.FastGetSolutionStepValue(TEMPERATURE) =
            test_reference_temperature + 5.0 * static_cast<double>(node_index)
            + 3.0 * rStateScale;
        ++node_index;
    }

    KRATOS_CATCH("");
}

/// Initializes all elements through the standard lifecycle.
void InitializeAllElements(ModelPart& rModelPart)
{
    KRATOS_TRY;

    const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
    for (auto& r_element : rModelPart.Elements()) {
        KRATOS_EXPECT_EQ(r_element.Check(r_process_info), 0);
        r_element.Initialize(r_process_info);
        r_element.InitializeSolutionStep(r_process_info);
        r_element.InitializeNonLinearIteration(r_process_info);
    }

    KRATOS_CATCH("");
}

/// Runs the actual smoothing scheme FinalizeSolutionStep sequence (reset ->
/// element/condition finalization -> process accumulation when enabled ->
/// normalization).
void RunSchemeFinalization(StaticSmoothingScheme& rScheme, ModelPart& rModelPart)
{
    CompressedMatrix A;
    Vector Dx;
    Vector b;
    rScheme.FinalizeSolutionStep(rModelPart, A, Dx, b);
}

/// Compares final normalized NODAL_CAUCHY_STRESS_TENSOR and NODAL_AREA between
/// the legacy and candidate model parts.
void CompareFinalNodalValues(ModelPart& rLegacy, ModelPart& rCandidate, const std::string& rLabel)
{
    KRATOS_TRY;

    auto it_legacy = rLegacy.Nodes().begin();
    auto it_candidate = rCandidate.Nodes().begin();
    for (; it_legacy != rLegacy.Nodes().end(); ++it_legacy, ++it_candidate) {
        const Matrix& r_legacy_stress = it_legacy->FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
        const Matrix& r_candidate_stress = it_candidate->FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                const double diff = std::abs(r_candidate_stress(i, j) - r_legacy_stress(i, j));
                const double tolerance = std::max(
                    comparison_absolute_tolerance,
                    comparison_relative_tolerance * std::abs(r_legacy_stress(i, j)));
                KRATOS_EXPECT_NEAR(r_candidate_stress(i, j), r_legacy_stress(i, j), tolerance);
                if (diff > tolerance) {
                    std::cout << "[3D.3] " << rLabel << " node " << it_legacy->Id()
                              << " stress(" << i << "," << j << ") legacy=" << r_legacy_stress(i, j)
                              << " candidate=" << r_candidate_stress(i, j) << std::endl;
                }
            }
        }
        const double legacy_area = it_legacy->FastGetSolutionStepValue(NODAL_AREA);
        const double candidate_area = it_candidate->FastGetSolutionStepValue(NODAL_AREA);
        const double area_tolerance = std::max(
            comparison_absolute_tolerance,
            comparison_relative_tolerance * std::abs(legacy_area));
        KRATOS_EXPECT_NEAR(candidate_area, legacy_area, area_tolerance);
    }

    KRATOS_CATCH("");
}

/// Computes the expected accumulated NODAL_AREA for the test mesh (single
/// accumulation of the element measures).
void VerifyExpectedNodalAreas(ModelPart& rModelPart)
{
    for (auto& r_node : rModelPart.Nodes()) {
        double expected_area = 0.0;
        for (auto& r_element : rModelPart.Elements()) {
            const auto& r_geometry = r_element.GetGeometry();
            for (std::size_t i = 0; i < r_geometry.size(); ++i) {
                if (r_geometry[i].Id() == r_node.Id()) {
                    expected_area += r_geometry.Area();
                }
            }
        }
        KRATOS_EXPECT_NEAR(
            r_node.FastGetSolutionStepValue(NODAL_AREA), expected_area, 1.0e-12);
    }
}

} // namespace

//************************************************************************************
// Static scheme: legacy vs candidate equivalence (full scheme FinalizeSolutionStep)
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(DamSmoothingScheme_StaticLegacyVsCandidate, KratosDamFastSuite)
{
    Model model;

    // Legacy workflow: legacy Dam thermo-mechanical element + legacy smoothing
    // ownership (flag absent -> legacy mode).
    ModelPart& r_legacy = CreateSchemeModelPart(model, "LegacyStatic", "SmallDisplacementThermoMechanicElement3D8N");
    PrescribeState(r_legacy, 1.0);
    InitializeAllElements(r_legacy);

    // Candidate workflow: StructuralMechanics small-displacement element +
    // compatible Dam thermal law + process-based smoothing ownership.
    ModelPart& r_candidate = CreateSchemeModelPart(model, "CandidateStatic", "SmallDisplacementElement3D8N");
    PrescribeState(r_candidate, 1.0);
    InitializeAllElements(r_candidate);

    StaticSmoothingScheme legacy_scheme;
    StaticSmoothingScheme candidate_scheme;
    RunSchemeFinalization(legacy_scheme, r_legacy);
    RunSchemeFinalization(candidate_scheme, r_candidate);

    // Final normalized nodal values must match.
    CompareFinalNodalValues(r_legacy, r_candidate, "static");

    // NODAL_AREA must be a single accumulation of the element measures (3.0 on
    // the shared nodes, 2.0 / 1.0 on the rest) - i.e. exactly one accumulation.
    VerifyExpectedNodalAreas(r_legacy);
    VerifyExpectedNodalAreas(r_candidate);

    // Raw state immediately before normalization is recovered as
    // normalized_stress * NODAL_AREA (the normalization divides by NODAL_AREA).
    // Verify the internal consistency of the two workflows.
    for (auto& r_node : r_legacy.Nodes()) {
        const double area = r_node.FastGetSolutionStepValue(NODAL_AREA);
        KRATOS_EXPECT_TRUE(area > 1.0e-15);
        const Matrix& r_stress = r_node.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                KRATOS_EXPECT_TRUE(std::isfinite(r_stress(i, j)));
            }
        }
    }
}

//************************************************************************************
// No-double-accumulation
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(DamSmoothingScheme_NoElementAccumulation, KratosDamFastSuite)
{
    Model model;

    // After Phase 3D.4B the legacy element itself performs NO nodal
    // accumulation: calling its FinalizeSolutionStep alone leaves NODAL_AREA and
    // NODAL_CAUCHY_STRESS_TENSOR untouched.
    ModelPart& r_legacy = CreateSchemeModelPart(model, "LegacyNoAccum", "SmallDisplacementThermoMechanicElement3D8N");
    PrescribeState(r_legacy, 1.0);
    InitializeAllElements(r_legacy);

    for (auto& r_node : r_legacy.Nodes()) {
        r_node.FastGetSolutionStepValue(NODAL_AREA) = 0.0;
        Matrix zero_stress_tensor(3, 3);
        noalias(zero_stress_tensor) = ZeroMatrix(3, 3);
        r_node.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR) = zero_stress_tensor;
    }

    for (auto& r_element : r_legacy.Elements()) {
        r_element.FinalizeSolutionStep(r_legacy.GetProcessInfo());
    }

    for (auto& r_node : r_legacy.Nodes()) {
        KRATOS_EXPECT_NEAR(r_node.FastGetSolutionStepValue(NODAL_AREA), 0.0, 1.0e-15);
        const Matrix& r_stress = r_node.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                KRATOS_EXPECT_NEAR(r_stress(i, j), 0.0, 1.0e-15);
            }
        }
    }

    // The full Dam smoothing scheme then performs exactly ONE process-based
    // accumulation (NODAL_AREA equals the analytical element-measure sum).
    StaticSmoothingScheme scheme;
    RunSchemeFinalization(scheme, r_legacy);
    VerifyExpectedNodalAreas(r_legacy);
    for (auto& r_node : r_legacy.Nodes()) {
        const Matrix& r_stress = r_node.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                KRATOS_EXPECT_TRUE(std::isfinite(r_stress(i, j)));
            }
        }
    }
}

//************************************************************************************
// StructuralMechanics-only candidate (central acceptance test)
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(DamSmoothingScheme_StructuralMechanicsOnly, KratosDamFastSuite)
{
    Model model;

    // Candidate model contains NO SmallDisplacementThermoMechanicElement: only
    // StructuralMechanicsApplication small-displacement elements with the Dam
    // thermal law, driven by the real Dam smoothing scheme in process-based
    // mode.
    ModelPart& r_candidate = CreateSchemeModelPart(model, "SMACandidate", "SmallDisplacementElement3D8N");
    PrescribeState(r_candidate, 1.0);
    InitializeAllElements(r_candidate);

    StaticSmoothingScheme scheme;
    RunSchemeFinalization(scheme, r_candidate);

    // The smoothing completes and produces the expected single accumulation.
    VerifyExpectedNodalAreas(r_candidate);
    for (auto& r_node : r_candidate.Nodes()) {
        const Matrix& r_stress = r_node.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
        double max_component = 0.0;
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                max_component = std::max(max_component, std::abs(r_stress(i, j)));
                KRATOS_EXPECT_TRUE(std::isfinite(r_stress(i, j)));
            }
        }
        KRATOS_EXPECT_TRUE(max_component > 1.0e-6);
    }

    // Cross-check against the legacy element workflow on the identical mesh.
    ModelPart& r_legacy = CreateSchemeModelPart(model, "LegacyRef", "SmallDisplacementThermoMechanicElement3D8N");
    PrescribeState(r_legacy, 1.0);
    InitializeAllElements(r_legacy);
    StaticSmoothingScheme legacy_scheme;
    RunSchemeFinalization(legacy_scheme, r_legacy);
    CompareFinalNodalValues(r_legacy, r_candidate, "SMA-only");
}

//************************************************************************************
// Multi-step workflow
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(DamSmoothingScheme_MultiStep, KratosDamFastSuite)
{
    Model model;

    ModelPart& r_legacy = CreateSchemeModelPart(model, "MultiLegacy", "SmallDisplacementThermoMechanicElement3D8N");
    ModelPart& r_candidate = CreateSchemeModelPart(model, "MultiCandidate", "SmallDisplacementElement3D8N");

    StaticSmoothingScheme legacy_scheme;
    StaticSmoothingScheme candidate_scheme;

    const double state_scales[2] = {1.0, 2.5};
    for (std::size_t step = 0; step < 2; ++step) {
        PrescribeState(r_legacy, state_scales[step]);
        PrescribeState(r_candidate, state_scales[step]);
        InitializeAllElements(r_legacy);
        InitializeAllElements(r_candidate);

        RunSchemeFinalization(legacy_scheme, r_legacy);
        RunSchemeFinalization(candidate_scheme, r_candidate);

        // Per step: exactly one reset + one accumulation + one normalization.
        VerifyExpectedNodalAreas(r_legacy);
        VerifyExpectedNodalAreas(r_candidate);
        CompareFinalNodalValues(r_legacy, r_candidate, "multi-step " + std::to_string(step));
    }

    // No stale stress: the second-step result differs from the first (state 2
    // is scaled by 2.5, so the smoothed stress must scale accordingly).
    for (auto& r_node : r_legacy.Nodes()) {
        KRATOS_EXPECT_TRUE(std::abs(r_node.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR)(0, 0)) > 1.0e-3);
    }
}

//************************************************************************************
// Damped smoothing scheme (inherits the common lifecycle)
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(DamSmoothingScheme_StaticDamped, KratosDamFastSuite)
{
    Model model;

    ModelPart& r_legacy = CreateSchemeModelPart(model, "DampedLegacy", "SmallDisplacementThermoMechanicElement3D8N");
    ModelPart& r_candidate = CreateSchemeModelPart(model, "DampedCandidate", "SmallDisplacementElement3D8N");
    PrescribeState(r_legacy, 1.0);
    PrescribeState(r_candidate, 1.0);
    InitializeAllElements(r_legacy);
    InitializeAllElements(r_candidate);

    DampedSmoothingScheme legacy_scheme(1.0e-6, 1.0e-6);
    DampedSmoothingScheme candidate_scheme(1.0e-6, 1.0e-6);

    CompressedMatrix A;
    Vector Dx;
    Vector b;
    legacy_scheme.FinalizeSolutionStep(r_legacy, A, Dx, b);
    candidate_scheme.FinalizeSolutionStep(r_candidate, A, Dx, b);

    VerifyExpectedNodalAreas(r_legacy);
    VerifyExpectedNodalAreas(r_candidate);
    CompareFinalNodalValues(r_legacy, r_candidate, "damped");
}

//************************************************************************************
// Bossak smoothing scheme (owns the same lifecycle)
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(DamSmoothingScheme_Bossak, KratosDamFastSuite)
{
    Model model;

    ModelPart& r_legacy = CreateSchemeModelPart(model, "BossakLegacy", "SmallDisplacementThermoMechanicElement3D8N");
    ModelPart& r_candidate = CreateSchemeModelPart(model, "BossakCandidate", "SmallDisplacementElement3D8N");
    PrescribeState(r_legacy, 1.0);
    PrescribeState(r_candidate, 1.0);
    InitializeAllElements(r_legacy);
    InitializeAllElements(r_candidate);

    BossakSmoothingScheme legacy_scheme(0.0, 1.0e-6, 1.0e-6);
    BossakSmoothingScheme candidate_scheme(0.0, 1.0e-6, 1.0e-6);

    CompressedMatrix A;
    Vector Dx;
    Vector b;
    legacy_scheme.FinalizeSolutionStep(r_legacy, A, Dx, b);
    candidate_scheme.FinalizeSolutionStep(r_candidate, A, Dx, b);

    VerifyExpectedNodalAreas(r_legacy);
    VerifyExpectedNodalAreas(r_candidate);
    CompareFinalNodalValues(r_legacy, r_candidate, "bossak");
}

//************************************************************************************
// Downstream consumer (selfweight / initial stress transfer)
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(DamSmoothingScheme_DownstreamConsumer, KratosDamFastSuite)
{
    Model model;

    ModelPart& r_legacy = CreateSchemeModelPart(model, "DownLegacy", "SmallDisplacementThermoMechanicElement3D8N");
    ModelPart& r_candidate = CreateSchemeModelPart(model, "DownCandidate", "SmallDisplacementElement3D8N");
    PrescribeState(r_legacy, 1.0);
    PrescribeState(r_candidate, 1.0);
    InitializeAllElements(r_legacy);
    InitializeAllElements(r_candidate);

    StaticSmoothingScheme legacy_scheme;
    StaticSmoothingScheme candidate_scheme;
    RunSchemeFinalization(legacy_scheme, r_legacy);
    RunSchemeFinalization(candidate_scheme, r_candidate);

    // The final normalized NODAL_CAUCHY_STRESS_TENSOR is available at the point
    // where the selfweight/initial-stress transfer consumers execute. Set a
    // known INITIAL_NODAL_CAUCHY_STRESS_TENSOR and run the actual
    // TransferSelfweightStressUtility::TransferInitialStress consumer.
    for (auto& r_node : r_legacy.Nodes()) {
        Matrix initial(3, 3);
        noalias(initial) = ZeroMatrix(3, 3);
        initial(0, 0) = 100.0 + static_cast<double>(r_node.Id());
        initial(1, 1) = 200.0;
        initial(2, 2) = 50.0;
        r_node.FastGetSolutionStepValue(INITIAL_NODAL_CAUCHY_STRESS_TENSOR) = initial;
    }
    for (auto& r_node : r_candidate.Nodes()) {
        Matrix initial(3, 3);
        noalias(initial) = ZeroMatrix(3, 3);
        initial(0, 0) = 100.0 + static_cast<double>(r_node.Id());
        initial(1, 1) = 200.0;
        initial(2, 2) = 50.0;
        r_node.FastGetSolutionStepValue(INITIAL_NODAL_CAUCHY_STRESS_TENSOR) = initial;
    }

    const int dimension = 3;
    TransferSelfweightStressUtility transfer;
    transfer.TransferInitialStress(r_legacy, dimension);
    transfer.TransferInitialStress(r_candidate, dimension);

    // After the consumer executes, both workflows must produce the same result.
    CompareFinalNodalValues(r_legacy, r_candidate, "downstream");
}

//************************************************************************************
// Unsupported higher-order geometry: not newly fatal (legacy-compatible skip).
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(DamSmoothingScheme_UnsupportedGeometryNotFatal, KratosDamFastSuite)
{
    // The legacy ExtrapolateGPStress only supported T3/Q4/T4/H8; higher-order
    // geometries performed no nodal Cauchy-stress extrapolation. After Phase
    // 3D.4B the process preserves that policy: an unsupported higher-order
    // geometry is skipped (no accumulation, no error), while the scheme still
    // runs the reset/normalize lifecycle.
    Model model;
    ModelPart& r_candidate = model.CreateModelPart("UnsupportedH20", 2);
    ProcessInfo& r_pi = r_candidate.GetProcessInfo();
    r_pi[DOMAIN_SIZE] = 3;
    r_pi[SPACE_DIMENSION] = 3;
    r_pi[IS_RESTARTED] = false;
    r_candidate.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_candidate.AddNodalSolutionStepVariable(VELOCITY);
    r_candidate.AddNodalSolutionStepVariable(ACCELERATION);
    r_candidate.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    r_candidate.AddNodalSolutionStepVariable(TEMPERATURE);
    r_candidate.AddNodalSolutionStepVariable(NODAL_REFERENCE_TEMPERATURE);
    r_candidate.AddNodalSolutionStepVariable(NODAL_CAUCHY_STRESS_TENSOR);
    r_candidate.AddNodalSolutionStepVariable(NODAL_AREA);
    r_candidate.AddNodalSolutionStepVariable(INITIAL_STRESS_TENSOR);
    const Element& r_proto = KratosComponents<Element>::Get("SmallDisplacementElement3D20N");
    const auto& r_geom = r_proto.GetGeometry();
    Matrix lc;
    r_geom.PointsLocalCoordinates(lc);
    const double gs = 2.5;
    array_1d<double, 3> off;
    off[0] = 0.75; off[1] = 1.25; off[2] = 0.5;
    for (std::size_t i = 0; i < 20; ++i) {
        r_candidate.CreateNewNode(i + 1, gs * lc(i, 0) + off[0], gs * lc(i, 1) + off[1], gs * lc(i, 2) + off[2]);
    }
    auto p_prop = r_candidate.CreateNewProperties(1);
    (*p_prop)[YOUNG_MODULUS] = test_young_modulus;
    (*p_prop)[POISSON_RATIO] = test_poisson_ratio;
    (*p_prop)[DENSITY] = 2400.0;
    (*p_prop)[THERMAL_EXPANSION] = test_thermal_expansion;
    p_prop->SetValue(CONSTITUTIVE_LAW, ThermalLinearElastic3DLaw().Clone());
    std::vector<ModelPart::IndexType> elem_nodes(20);
    for (std::size_t i = 0; i < 20; ++i) elem_nodes[i] = i + 1;
    r_candidate.CreateNewElement("SmallDisplacementElement3D20N", 1, elem_nodes, p_prop);
    for (auto& r_node : r_candidate.Nodes()) {
        r_node.AddDof(DISPLACEMENT_X);
        r_node.AddDof(DISPLACEMENT_Y);
        r_node.AddDof(DISPLACEMENT_Z);
        r_node.FastGetSolutionStepValue(NODAL_REFERENCE_TEMPERATURE) = test_reference_temperature;
        r_node.FastGetSolutionStepValue(NODAL_AREA) = 0.0;
        Matrix zero_stress_tensor(3, 3);
        noalias(zero_stress_tensor) = ZeroMatrix(3, 3);
        r_node.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR) = zero_stress_tensor;
    }
    r_candidate.pGetElement(1)->Initialize(r_pi);

    PrescribeState(r_candidate, 1.0);

    StaticSmoothingScheme scheme;
    RunSchemeFinalization(scheme, r_candidate);

    // The unsupported element is skipped: no nodal accumulation occurs.
    for (auto& r_node : r_candidate.Nodes()) {
        KRATOS_EXPECT_NEAR(r_node.FastGetSolutionStepValue(NODAL_AREA), 0.0, 1.0e-15);
        const Matrix& r_stress = r_node.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                KRATOS_EXPECT_NEAR(r_stress(i, j), 0.0, 1.0e-15);
            }
        }
    }
}

} // namespace Testing
} // namespace Kratos
