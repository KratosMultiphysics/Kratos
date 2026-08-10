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

// Lifecycle and state-robustness characterization of the thermo-mechanical
// compatibility between:
//   StructuralMechanicsApplication::SmallDisplacement
//   + DamApplication::ThermalLinearElastic3DLaw
// against the legacy:
//   DamApplication::SmallDisplacementThermoMechanicElement
//   + DamApplication::ThermalLinearElastic3DLaw
//
// These tests verify that the PK2/Kirchhoff/Cauchy material-response methods
// and the RequiresInitializeMaterialResponse()/RequiresFinalizeMaterialResponse()
// overrides of ThermalLinearElastic3DLaw remain correct:
//   - over multiple consecutive solution steps (no hidden state accumulation);
//   - when the constitutive response is evaluated repeatedly with identical
//     inputs within one step, before and after the lifecycle callbacks;
//   - after cloning an initialized law;
//   - after serialization/deserialization (restart).
//
// A Hexahedra3D8 (3D8N) is used for all tests.

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
#include "includes/stream_serializer.h"
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

/// Comparison tolerances (same philosophy as the other characterization tests).
constexpr double comparison_absolute_tolerance = 1.0e-12;
constexpr double comparison_relative_tolerance = 1.0e-10;
constexpr double machine_precision_allowance = 1.0e-15;

/// Material data shared by both model parts.
constexpr double test_young_modulus = 2.0e7;
constexpr double test_poisson_ratio = 0.2;
constexpr double test_density = 2400.0;
constexpr double test_thermal_expansion = 1.0e-5;
constexpr double test_reference_temperature = 20.0;

/// Element responses compared after every solution step.
struct StepResponse
{
    Matrix lhs;
    Vector rhs;
    std::vector<Vector> strain_vectors;
    std::vector<Vector> cauchy_stress_vectors;
    std::vector<Vector> pk2_stress_vectors;
};

/// Creates a 3D8N model part with the given element and the thermal law.
/// The nodal solution-step data includes the variables required by the legacy
/// element FinalizeSolutionStep (NODAL_CAUCHY_STRESS_TENSOR, NODAL_AREA) and by
/// the thermal law (TEMPERATURE, NODAL_REFERENCE_TEMPERATURE).
ModelPart& CreateStateModelPart(
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
    r_process_info[DELTA_TIME] = 0.1;

    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(VELOCITY);
    r_model_part.AddNodalSolutionStepVariable(ACCELERATION);
    r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    r_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    r_model_part.AddNodalSolutionStepVariable(NODAL_REFERENCE_TEMPERATURE);
    // Required by the legacy element FinalizeSolutionStep (nodal stress
    // extrapolation).
    r_model_part.AddNodalSolutionStepVariable(NODAL_CAUCHY_STRESS_TENSOR);
    r_model_part.AddNodalSolutionStepVariable(NODAL_AREA);

    // Nodes on a scaled and translated copy of the prototype hexahedron so that
    // both elements are tested with exactly the same geometry.
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
        r_node.FastGetSolutionStepValue(NODAL_AREA) = 0.0;
        Matrix zero_stress_tensor(3, 3);
        noalias(zero_stress_tensor) = ZeroMatrix(3, 3);
        r_node.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR) = zero_stress_tensor;
    }

    return r_model_part;

    KRATOS_CATCH("");
}

/// Prescribes the state of a given solution step. Three consecutive steps are
/// defined (indices 0, 1, 2), each with a different displacement field and a
/// different temperature field:
///   step 0: reference temperature, small non-zero mechanical displacement;
///   step 1: uniform temperature increment, changed displacement field;
///   step 2: non-uniform nodal temperature field, another displacement field.
void PrescribeStepState(ModelPart& rModelPart, const std::size_t rStepIndex)
{
    KRATOS_TRY;

    Matrix displacement_gradient = ZeroMatrix(3, 3);
    array_1d<double, 3> translation;
    switch (rStepIndex) {
    case 0:
        displacement_gradient(0, 0) = 2.0e-3;
        displacement_gradient(0, 1) = 3.0e-4;
        displacement_gradient(0, 2) = 1.0e-4;
        displacement_gradient(1, 0) = 3.0e-4;
        displacement_gradient(1, 1) = -1.0e-3;
        displacement_gradient(1, 2) = 2.0e-4;
        displacement_gradient(2, 0) = 1.0e-4;
        displacement_gradient(2, 1) = 2.0e-4;
        displacement_gradient(2, 2) = 5.0e-4;
        translation[0] = 1.0e-2;
        translation[1] = -2.0e-2;
        translation[2] = 5.0e-3;
        break;
    case 1:
        displacement_gradient(0, 0) = 1.0e-3;
        displacement_gradient(0, 1) = -2.0e-4;
        displacement_gradient(0, 2) = 5.0e-5;
        displacement_gradient(1, 0) = -2.0e-4;
        displacement_gradient(1, 1) = 2.0e-3;
        displacement_gradient(1, 2) = -1.0e-4;
        displacement_gradient(2, 0) = 5.0e-5;
        displacement_gradient(2, 1) = -1.0e-4;
        displacement_gradient(2, 2) = -8.0e-4;
        translation[0] = -5.0e-3;
        translation[1] = 1.0e-2;
        translation[2] = -1.0e-2;
        break;
    default:
        displacement_gradient(0, 0) = -1.5e-3;
        displacement_gradient(0, 1) = 1.0e-4;
        displacement_gradient(0, 2) = -2.0e-4;
        displacement_gradient(1, 0) = 1.0e-4;
        displacement_gradient(1, 1) = 1.2e-3;
        displacement_gradient(1, 2) = 3.0e-4;
        displacement_gradient(2, 0) = -2.0e-4;
        displacement_gradient(2, 1) = 3.0e-4;
        displacement_gradient(2, 2) = 6.0e-4;
        translation[0] = 2.0e-2;
        translation[1] = 5.0e-3;
        translation[2] = -1.5e-2;
        break;
    }

    std::size_t node_index = 0;
    for (auto& r_node : rModelPart.Nodes()) {
        array_1d<double, 3>& r_displacement = r_node.FastGetSolutionStepValue(DISPLACEMENT);
        noalias(r_displacement) = prod(displacement_gradient, r_node.GetInitialPosition());
        noalias(r_displacement) += translation;
        // Updated coordinates so that the reference configuration (current
        // position - total displacement) is the initial one.
        r_node.X() = r_node.X0() + r_displacement[0];
        r_node.Y() = r_node.Y0() + r_displacement[1];
        r_node.Z() = r_node.Z0() + r_displacement[2];

        double nodal_temperature = test_reference_temperature;
        if (rStepIndex == 1) {
            nodal_temperature = test_reference_temperature + 25.0;
        } else if (rStepIndex == 2) {
            nodal_temperature = test_reference_temperature + 5.0 * static_cast<double>(node_index);
        }
        r_node.FastGetSolutionStepValue(TEMPERATURE) = nodal_temperature;

        array_1d<double, 3>& r_volume_acceleration = r_node.FastGetSolutionStepValue(VOLUME_ACCELERATION);
        r_volume_acceleration[0] = 0.0;
        r_volume_acceleration[1] = -9.81;
        r_volume_acceleration[2] = 0.0;

        ++node_index;
    }

    KRATOS_CATCH("");
}

/// Calculates the compared responses of the element of the model part.
void CalculateStepResponses(ModelPart& rModelPart, StepResponse& rResponse)
{
    KRATOS_TRY;

    auto p_element = rModelPart.pGetElement(1);
    const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

    p_element->CalculateLocalSystem(rResponse.lhs, rResponse.rhs, r_process_info);
    p_element->CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR, rResponse.strain_vectors, r_process_info);
    p_element->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, rResponse.cauchy_stress_vectors, r_process_info);
    p_element->CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, rResponse.pk2_stress_vectors, r_process_info);

    KRATOS_CATCH("");
}

/// Tolerance associated with a reference value (same philosophy as the other
/// characterization tests).
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

/// Metrics accumulated over the components of one comparison (diagnostic only;
/// the pass/fail decision is the component-wise EXPECT check).
struct ComparisonMetrics
{
    bool pass = true;
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

/// Compares the legacy and candidate responses for one solution step.
/// The candidate PK2 stress is checked against the physically equivalent legacy
/// thermo-elastic stress (the legacy Cauchy stress).
void CompareStepResponses(
    const std::size_t rStepIndex,
    const StepResponse& rLegacy,
    const StepResponse& rCandidate)
{
    const std::string step = "step" + std::to_string(rStepIndex + 1);
    ExpectMatrixComponentsNear(rLegacy.lhs, rCandidate.lhs, step + " local system LHS");
    ExpectVectorComponentsNear(rLegacy.rhs, rCandidate.rhs, step + " local system RHS");
    ExpectIntegrationPointVectorsNear(rLegacy.strain_vectors, rCandidate.strain_vectors, step + " strain vector");
    ExpectIntegrationPointVectorsNear(rLegacy.cauchy_stress_vectors, rCandidate.cauchy_stress_vectors, step + " Cauchy stress vector");
    ExpectIntegrationPointVectorsNear(rCandidate.pk2_stress_vectors, rLegacy.cauchy_stress_vectors, step + " candidate PK2 vs legacy Cauchy stress vector");
}

/// Owns the vectors bound to a ConstitutiveLaw::Parameters (the Parameters
/// stores pointers to them, so they must outlive the response evaluation).
struct LawParametersBundle
{
    Vector strain;
    Vector stress;
    Matrix constitutive_matrix;
    Vector shape_function_values;
    ConstitutiveLaw::Parameters values;

    LawParametersBundle(
        const Geometry<Node>& rGeometry,
        const Properties& rProperties,
        const ProcessInfo& rProcessInfo,
        const Vector& rStrain)
        : strain(rStrain),
          stress(rStrain.size()),
          constitutive_matrix(rStrain.size(), rStrain.size()),
          shape_function_values(rGeometry.PointsNumber()),
          values(rGeometry, rProperties, rProcessInfo)
    {
        noalias(stress) = ZeroVector(rStrain.size());
        noalias(constitutive_matrix) = ZeroMatrix(rStrain.size(), rStrain.size());
        noalias(shape_function_values) = row(rGeometry.ShapeFunctionsValues(), 0);

        auto& r_options = values.GetOptions();
        r_options.Set(ConstitutiveLaw::COMPUTE_STRESS);
        r_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        r_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
        values.SetStrainVector(strain);
        values.SetStressVector(stress);
        values.SetConstitutiveMatrix(constitutive_matrix);
        values.SetShapeFunctionsValues(shape_function_values);
    }
};

} // namespace

//************************************************************************************
// Constitutive-law features
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalStateRobustness_GetLawFeatures, KratosDamFastSuite)
{
    // GetLawFeatures() is inherited from LinearElastic3DLaw (Poromechanics):
    // the law advertises an infinitesimal-strain formulation, which is
    // consistent with PK2 == Kirchhoff == Cauchy returning the same
    // infinitesimal constitutive stress.
    ThermalLinearElastic3DLaw law;
    ConstitutiveLaw::Features features;
    law.GetLawFeatures(features);

    bool has_infinitesimal_measure = false;
    for (const auto& r_measure : features.mStrainMeasures) {
        if (r_measure == ConstitutiveLaw::StrainMeasure_Infinitesimal) {
            has_infinitesimal_measure = true;
        }
    }
    KRATOS_EXPECT_TRUE(has_infinitesimal_measure);
    KRATOS_EXPECT_TRUE(features.mOptions.Is(ConstitutiveLaw::INFINITESIMAL_STRAINS));
    KRATOS_EXPECT_EQ(features.mStrainSize, 6);
    KRATOS_EXPECT_EQ(features.mSpaceDimension, 3);
    std::cout << "[CHARACTERIZATION] ThermalLinearElastic3DLaw features: "
              << "strain_measure=StrainMeasure_Infinitesimal, "
              << "strain_size=" << features.mStrainSize
              << ", space_dimension=" << features.mSpaceDimension
              << std::endl;
}

//************************************************************************************
// Multi-step characterization
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalStateRobustness_MultiStep3D8N, KratosDamFastSuite)
{
    Model model;
    ModelPart& r_legacy = CreateStateModelPart(
        model, "LegacyMultiStep", "SmallDisplacementThermoMechanicElement3D8N", 0);
    ModelPart& r_candidate = CreateStateModelPart(
        model, "CandidateMultiStep", "SmallDisplacementElement3D8N", 100);

    auto p_legacy_element = r_legacy.pGetElement(1);
    auto p_candidate_element = r_candidate.pGetElement(1);
    ProcessInfo& r_legacy_pi = r_legacy.GetProcessInfo();
    ProcessInfo& r_candidate_pi = r_candidate.GetProcessInfo();

    KRATOS_EXPECT_EQ(p_legacy_element->Check(r_legacy_pi), 0);
    KRATOS_EXPECT_EQ(p_candidate_element->Check(r_candidate_pi), 0);
    p_legacy_element->Initialize(r_legacy_pi);
    p_candidate_element->Initialize(r_candidate_pi);

    // Three consecutive solution steps with different mechanical and thermal
    // states. ProcessInfo time/step data is advanced consistently.
    for (std::size_t step = 0; step < 3; ++step) {
        r_legacy_pi[STEP] = step + 1;
        r_legacy_pi[TIME] = 0.1 * static_cast<double>(step + 1);
        r_candidate_pi[STEP] = step + 1;
        r_candidate_pi[TIME] = 0.1 * static_cast<double>(step + 1);

        PrescribeStepState(r_legacy, step);
        PrescribeStepState(r_candidate, step);

        p_legacy_element->InitializeSolutionStep(r_legacy_pi);
        p_candidate_element->InitializeSolutionStep(r_candidate_pi);
        p_legacy_element->InitializeNonLinearIteration(r_legacy_pi);
        p_candidate_element->InitializeNonLinearIteration(r_candidate_pi);

        StepResponse legacy_response, candidate_response;
        CalculateStepResponses(r_legacy, legacy_response);
        CalculateStepResponses(r_candidate, candidate_response);

        p_legacy_element->FinalizeNonLinearIteration(r_legacy_pi);
        p_candidate_element->FinalizeNonLinearIteration(r_candidate_pi);
        p_legacy_element->FinalizeSolutionStep(r_legacy_pi);
        p_candidate_element->FinalizeSolutionStep(r_candidate_pi);

        CompareStepResponses(step, legacy_response, candidate_response);
    }
}

//************************************************************************************
// Repeated-response test
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalStateRobustness_RepeatedResponse, KratosDamFastSuite)
{
    // The constitutive response must not depend on the number of previous
    // identical evaluations, nor on the lifecycle callbacks. The legacy element
    // is included because its FinalizeSolutionStep invokes
    // FinalizeMaterialResponseCauchy, which updates inherited (unused) state.
    for (const char* p_element_name : {"SmallDisplacementThermoMechanicElement3D8N",
                                        "SmallDisplacementElement3D8N"}) {
        const std::string r_element_name(p_element_name);
        Model model;
        ModelPart& r_model_part = CreateStateModelPart(
            model, "Repeated" + r_element_name, r_element_name, 0);
        PrescribeStepState(r_model_part, 1);
        auto p_element = r_model_part.pGetElement(1);
        const ProcessInfo& r_pi = r_model_part.GetProcessInfo();

        KRATOS_EXPECT_EQ(p_element->Check(r_pi), 0);
        p_element->Initialize(r_pi);

        std::vector<Vector> reference;
        p_element->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, reference, r_pi);

        // Repeated evaluation with identical inputs.
        std::vector<Vector> repeated;
        p_element->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, repeated, r_pi);
        ExpectIntegrationPointVectorsNear(
            repeated, reference, r_element_name + " repeated response (no lifecycle)");

        // After InitializeSolutionStep + InitializeNonLinearIteration.
        p_element->InitializeSolutionStep(r_pi);
        p_element->InitializeNonLinearIteration(r_pi);
        std::vector<Vector> after_initialize;
        p_element->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, after_initialize, r_pi);
        ExpectIntegrationPointVectorsNear(
            after_initialize, reference, r_element_name + " repeated response (after initialize)");

        // After FinalizeNonLinearIteration + FinalizeSolutionStep.
        p_element->FinalizeNonLinearIteration(r_pi);
        p_element->FinalizeSolutionStep(r_pi);
        std::vector<Vector> after_finalize;
        p_element->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, after_finalize, r_pi);
        ExpectIntegrationPointVectorsNear(
            after_finalize, reference, r_element_name + " repeated response (after finalize)");
    }
}

//************************************************************************************
// Clone test
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalStateRobustness_CloneLaw, KratosDamFastSuite)
{
    // Clone an initialized ThermalLinearElastic3DLaw and verify that the
    // original and the clone produce identical responses for identical
    // ConstitutiveLaw::Parameters (including a non-zero temperature increment).
    Model model;
    ModelPart& r_model_part = CreateStateModelPart(
        model, "CloneLaw", "SmallDisplacementElement3D8N", 0);
    auto p_element = r_model_part.pGetElement(1);
    const auto& r_geometry = p_element->GetGeometry();
    const Properties& r_properties = p_element->GetProperties();
    const ProcessInfo& r_pi = r_model_part.GetProcessInfo();

    // Initialize a fresh law (the property's law is cloned by the element
    // during Initialize; initialize it here explicitly for the clone test).
    ConstitutiveLaw::Pointer p_law = r_properties.GetValue(CONSTITUTIVE_LAW);
    Vector shape_function_values(r_geometry.PointsNumber());
    noalias(shape_function_values) = row(r_geometry.ShapeFunctionsValues(), 0);
    p_law->InitializeMaterial(r_properties, r_geometry, shape_function_values);

    ConstitutiveLaw::Pointer p_clone = p_law->Clone();
    KRATOS_EXPECT_TRUE(p_clone != nullptr);

    // Non-zero strain and a non-zero temperature increment (the model part nodes
    // carry TEMPERATURE = reference + 25 in step 1 state).
    PrescribeStepState(r_model_part, 1);
    Vector strain(6);
    noalias(strain) = ZeroVector(6);
    strain[0] = 2.0e-3;
    strain[1] = -1.0e-3;
    strain[2] = 5.0e-4;
    strain[3] = 1.0e-4;
    strain[4] = -2.0e-4;
    strain[5] = 3.0e-4;

    LawParametersBundle original(r_geometry, r_properties, r_pi, strain);
    LawParametersBundle clone(r_geometry, r_properties, r_pi, strain);

    p_law->CalculateMaterialResponse(original.values, ConstitutiveLaw::StressMeasure_PK2);
    p_clone->CalculateMaterialResponse(clone.values, ConstitutiveLaw::StressMeasure_PK2);

    ExpectVectorComponentsNear(
        original.values.GetStressVector(), clone.values.GetStressVector(), "clone PK2 stress");
    ExpectMatrixComponentsNear(
        original.values.GetConstitutiveMatrix(), clone.values.GetConstitutiveMatrix(), "clone constitutive matrix");

    // Also compare the Cauchy measure of the clone with its PK2 stress (for an
    // infinitesimal-strain law the measures coincide).
    LawParametersBundle clone_cauchy(r_geometry, r_properties, r_pi, strain);
    p_clone->CalculateMaterialResponse(clone_cauchy.values, ConstitutiveLaw::StressMeasure_Cauchy);
    ExpectVectorComponentsNear(
        clone_cauchy.values.GetStressVector(), clone.values.GetStressVector(),
        "clone Cauchy stress equals clone PK2 stress");
}

//************************************************************************************
// Serialization / restart test
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalStateRobustness_Serialization, KratosDamFastSuite)
{
    // Serialize and deserialize the constitutive law with the standard Kratos
    // in-memory StreamSerializer (no external restart files) and verify that
    // the deserialized law reproduces the same thermo-mechanical response.
    Model model;
    ModelPart& r_model_part = CreateStateModelPart(
        model, "SerializationLaw", "SmallDisplacementElement3D8N", 0);
    auto p_element = r_model_part.pGetElement(1);
    const auto& r_geometry = p_element->GetGeometry();
    const Properties& r_properties = p_element->GetProperties();
    const ProcessInfo& r_pi = r_model_part.GetProcessInfo();

    ConstitutiveLaw::Pointer p_law = r_properties.GetValue(CONSTITUTIVE_LAW);
    Vector shape_function_values(r_geometry.PointsNumber());
    noalias(shape_function_values) = row(r_geometry.ShapeFunctionsValues(), 0);
    p_law->InitializeMaterial(r_properties, r_geometry, shape_function_values);

    // Non-zero temperature increment.
    PrescribeStepState(r_model_part, 1);
    Vector strain(6);
    noalias(strain) = ZeroVector(6);
    strain[0] = 2.0e-3;
    strain[1] = -1.0e-3;
    strain[2] = 5.0e-4;
    strain[3] = 1.0e-4;
    strain[4] = -2.0e-4;
    strain[5] = 3.0e-4;

    LawParametersBundle values_before(r_geometry, r_properties, r_pi, strain);
    p_law->CalculateMaterialResponse(values_before.values, ConstitutiveLaw::StressMeasure_PK2);

    // Serialize / deserialize.
    StreamSerializer serializer;
    serializer.save("ThermalLinearElastic3DLaw", p_law);
    ConstitutiveLaw::Pointer p_loaded;
    serializer.load("ThermalLinearElastic3DLaw", p_loaded);
    KRATOS_EXPECT_TRUE(p_loaded != nullptr);

    LawParametersBundle values_after(r_geometry, r_properties, r_pi, strain);
    p_loaded->CalculateMaterialResponse(values_after.values, ConstitutiveLaw::StressMeasure_PK2);

    // Compare stress and constitutive matrix.
    ExpectVectorComponentsNear(
        values_before.values.GetStressVector(), values_after.values.GetStressVector(),
        "serialization PK2 stress (after deserialization)");
    ExpectMatrixComponentsNear(
        values_before.values.GetConstitutiveMatrix(), values_after.values.GetConstitutiveMatrix(),
        "serialization constitutive matrix (after deserialization)");

    // Apply the next thermo-mechanical state (a different temperature field) on
    // the deserialized law and compare with the non-serialized reference.
    PrescribeStepState(r_model_part, 2);
    Vector strain_next(6);
    noalias(strain_next) = ZeroVector(6);
    strain_next[0] = -1.0e-3;
    strain_next[1] = 1.5e-3;
    strain_next[2] = 4.0e-4;
    strain_next[3] = -5.0e-4;
    strain_next[4] = 2.0e-4;
    strain_next[5] = 6.0e-4;

    LawParametersBundle values_before_next(r_geometry, r_properties, r_pi, strain_next);
    p_law->CalculateMaterialResponse(values_before_next.values, ConstitutiveLaw::StressMeasure_PK2);

    LawParametersBundle values_after_next(r_geometry, r_properties, r_pi, strain_next);
    p_loaded->CalculateMaterialResponse(values_after_next.values, ConstitutiveLaw::StressMeasure_PK2);

    ExpectVectorComponentsNear(
        values_before_next.values.GetStressVector(), values_after_next.values.GetStressVector(),
        "serialization PK2 stress (next state, after deserialization)");
    ExpectMatrixComponentsNear(
        values_before_next.values.GetConstitutiveMatrix(), values_after_next.values.GetConstitutiveMatrix(),
        "serialization constitutive matrix (next state, after deserialization)");
}

} // namespace Testing
} // namespace Kratos
