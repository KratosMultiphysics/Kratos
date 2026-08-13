// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    DamApplication developers

// Permanent Dam smoothing-scheme integration contract. The Dam smoothing scheme
// owns the nodal Cauchy-stress extrapolation: the element performs no nodal
// accumulation (single-owner), and unsupported higher-order geometries are
// skipped without error. The historical NODAL_AREA / NODAL_CAUCHY_STRESS_TENSOR
// variables are reset and normalized by the scheme.
//
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

/// Comparison tolerances.
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
// No-double-accumulation
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(DamSmoothingScheme_NoElementAccumulation, KratosDamFastSuite)
{
    Model model;

    // The element itself performs NO nodal accumulation: calling its
    // FinalizeSolutionStep alone leaves NODAL_AREA and
    // NODAL_CAUCHY_STRESS_TENSOR untouched.
    ModelPart& r_mp = CreateSchemeModelPart(model, "NoAccum", "SmallDisplacementThermoMechanicElement3D8N");
    PrescribeState(r_mp, 1.0);
    InitializeAllElements(r_mp);

    for (auto& r_node : r_mp.Nodes()) {
        r_node.FastGetSolutionStepValue(NODAL_AREA) = 0.0;
        Matrix zero_stress_tensor(3, 3);
        noalias(zero_stress_tensor) = ZeroMatrix(3, 3);
        r_node.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR) = zero_stress_tensor;
    }

    for (auto& r_element : r_mp.Elements()) {
        r_element.FinalizeSolutionStep(r_mp.GetProcessInfo());
    }

    for (auto& r_node : r_mp.Nodes()) {
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
    RunSchemeFinalization(scheme, r_mp);
    VerifyExpectedNodalAreas(r_mp);
    for (auto& r_node : r_mp.Nodes()) {
        const Matrix& r_stress = r_node.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                KRATOS_EXPECT_TRUE(std::isfinite(r_stress(i, j)));
            }
        }
    }
}












//************************************************************************************
// Unsupported higher-order geometry is skipped without error.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(DamSmoothingScheme_UnsupportedGeometryNotFatal, KratosDamFastSuite)
{
    // Extrapolation supports only T3/Q4/T4/H8; an unsupported higher-order
    // geometry is skipped (no accumulation, no error), while the scheme still
    // runs the reset/normalize lifecycle.
    Model model;
    ModelPart& r_mp = model.CreateModelPart("UnsupportedH20", 2);
    ProcessInfo& r_pi = r_mp.GetProcessInfo();
    r_pi[DOMAIN_SIZE] = 3;
    r_pi[SPACE_DIMENSION] = 3;
    r_pi[IS_RESTARTED] = false;
    r_mp.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_mp.AddNodalSolutionStepVariable(VELOCITY);
    r_mp.AddNodalSolutionStepVariable(ACCELERATION);
    r_mp.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    r_mp.AddNodalSolutionStepVariable(TEMPERATURE);
    r_mp.AddNodalSolutionStepVariable(NODAL_REFERENCE_TEMPERATURE);
    r_mp.AddNodalSolutionStepVariable(NODAL_CAUCHY_STRESS_TENSOR);
    r_mp.AddNodalSolutionStepVariable(NODAL_AREA);
    r_mp.AddNodalSolutionStepVariable(INITIAL_STRESS_TENSOR);
    const Element& r_proto = KratosComponents<Element>::Get("SmallDisplacementElement3D20N");
    const auto& r_geom = r_proto.GetGeometry();
    Matrix lc;
    r_geom.PointsLocalCoordinates(lc);
    const double gs = 2.5;
    array_1d<double, 3> off;
    off[0] = 0.75; off[1] = 1.25; off[2] = 0.5;
    for (std::size_t i = 0; i < 20; ++i) {
        r_mp.CreateNewNode(i + 1, gs * lc(i, 0) + off[0], gs * lc(i, 1) + off[1], gs * lc(i, 2) + off[2]);
    }
    auto p_prop = r_mp.CreateNewProperties(1);
    (*p_prop)[YOUNG_MODULUS] = test_young_modulus;
    (*p_prop)[POISSON_RATIO] = test_poisson_ratio;
    (*p_prop)[DENSITY] = 2400.0;
    (*p_prop)[THERMAL_EXPANSION] = test_thermal_expansion;
    p_prop->SetValue(CONSTITUTIVE_LAW, ThermalLinearElastic3DLaw().Clone());
    std::vector<ModelPart::IndexType> elem_nodes(20);
    for (std::size_t i = 0; i < 20; ++i) elem_nodes[i] = i + 1;
    r_mp.CreateNewElement("SmallDisplacementElement3D20N", 1, elem_nodes, p_prop);
    for (auto& r_node : r_mp.Nodes()) {
        r_node.AddDof(DISPLACEMENT_X);
        r_node.AddDof(DISPLACEMENT_Y);
        r_node.AddDof(DISPLACEMENT_Z);
        r_node.FastGetSolutionStepValue(NODAL_REFERENCE_TEMPERATURE) = test_reference_temperature;
        r_node.FastGetSolutionStepValue(NODAL_AREA) = 0.0;
        Matrix zero_stress_tensor(3, 3);
        noalias(zero_stress_tensor) = ZeroMatrix(3, 3);
        r_node.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR) = zero_stress_tensor;
    }
    r_mp.pGetElement(1)->Initialize(r_pi);

    PrescribeState(r_mp, 1.0);

    StaticSmoothingScheme scheme;
    RunSchemeFinalization(scheme, r_mp);

    // The unsupported element is skipped: no nodal accumulation occurs.
    for (auto& r_node : r_mp.Nodes()) {
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
