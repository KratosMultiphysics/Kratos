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

// Phase 3A characterization: the thermo-mechanical vector outputs
//   THERMAL_STRAIN_VECTOR, THERMAL_STRESS_VECTOR, MECHANICAL_STRESS_VECTOR
// are now provided directly by the thermal constitutive laws through the
// standard BaseSolidElement::CalculateOnIntegrationPoints infrastructure of the
// StructuralMechanicsApplication SmallDisplacement element.
//
// Output semantics (implemented in ThermalLinearElastic3DLaw::CalculateValue):
//   THERMAL_STRAIN_VECTOR    = epsilon_th
//   THERMAL_STRESS_VECTOR    = C * epsilon_th
//   MECHANICAL_STRESS_VECTOR = C * epsilon
// so that
//   CAUCHY_STRESS_VECTOR = MECHANICAL_STRESS_VECTOR - THERMAL_STRESS_VECTOR
//   PK2_STRESS_VECTOR    = CAUCHY_STRESS_VECTOR       (small strains)
//
// The legacy SmallDisplacementThermoMechanicElement does not return the actual
// thermal strain for THERMAL_STRAIN_VECTOR (it returns the element strain). The
// new constitutive-law output intentionally fixes that: THERMAL_STRAIN_VECTOR
// is compared against the analytical epsilon_th and never against the broken
// legacy output.

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
#include "custom_constitutive/thermal_linear_elastic_2D_plane_strain.hpp"
#include "custom_constitutive/thermal_linear_elastic_2D_plane_stress.hpp"

namespace Kratos
{
namespace Testing
{

namespace
{

/// Comparison tolerances (same philosophy as the previous characterization).
constexpr double comparison_absolute_tolerance = 1.0e-12;
constexpr double comparison_relative_tolerance = 1.0e-10;
constexpr double machine_precision_allowance = 1.0e-15;

/// Material data shared by all model parts.
constexpr double test_young_modulus = 2.0e7;
constexpr double test_poisson_ratio = 0.2;
constexpr double test_density = 2400.0;
constexpr double test_thermal_expansion = 1.0e-5;
constexpr double test_reference_temperature = 20.0;
constexpr double test_thickness = 0.15;

/// Creates the constitutive law for the given configuration.
ConstitutiveLaw::Pointer CreateThermalLaw(const std::string& rLawName)
{
    if (rLawName == "ThermalLinearElastic2DPlaneStrain") {
        return ThermalLinearElastic2DPlaneStrain().Clone();
    } else if (rLawName == "ThermalLinearElastic2DPlaneStress") {
        return ThermalLinearElastic2DPlaneStress().Clone();
    } else {
        return ThermalLinearElastic3DLaw().Clone();
    }
}

/// Creates one of the two numerically identical model parts.
ModelPart& CreateVectorModelPart(
    Model& rModel,
    const std::string& rModelPartName,
    const std::string& rElementName,
    const std::string& rLawName,
    const std::size_t rDimension,
    const ModelPart::IndexType rNodeIdOffset)
{
    KRATOS_TRY;

    ModelPart& r_model_part = rModel.CreateModelPart(rModelPartName, 2);

    ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
    r_process_info[DOMAIN_SIZE] = rDimension;
    r_process_info[SPACE_DIMENSION] = rDimension;
    r_process_info[IS_RESTARTED] = false;

    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(VELOCITY);
    r_model_part.AddNodalSolutionStepVariable(ACCELERATION);
    r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    r_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    r_model_part.AddNodalSolutionStepVariable(NODAL_REFERENCE_TEMPERATURE);
    r_model_part.AddNodalSolutionStepVariable(NODAL_CAUCHY_STRESS_TENSOR);
    r_model_part.AddNodalSolutionStepVariable(NODAL_AREA);

    const Element& r_prototype = KratosComponents<Element>::Get(rElementName);
    const auto& r_geometry = r_prototype.GetGeometry();
    Matrix local_coordinates;
    r_geometry.PointsLocalCoordinates(local_coordinates);

    const double geometry_scale = 2.5;
    array_1d<double, 3> geometry_offset;
    geometry_offset[0] = 0.75;
    geometry_offset[1] = 1.25;
    geometry_offset[2] = (rDimension == 3) ? 0.5 : 0.0;

    const std::size_t number_of_nodes = r_geometry.PointsNumber();
    for (std::size_t i = 0; i < number_of_nodes; ++i) {
        const double x = geometry_scale * local_coordinates(i, 0) + geometry_offset[0];
        const double y = geometry_scale * local_coordinates(i, 1) + geometry_offset[1];
        const double z = (rDimension == 3) ? geometry_scale * local_coordinates(i, 2) + geometry_offset[2] : 0.0;
        r_model_part.CreateNewNode(rNodeIdOffset + i + 1, x, y, z);
    }

    auto p_prop = r_model_part.CreateNewProperties(1);
    (*p_prop)[YOUNG_MODULUS] = test_young_modulus;
    (*p_prop)[POISSON_RATIO] = test_poisson_ratio;
    (*p_prop)[DENSITY] = test_density;
    (*p_prop)[THERMAL_EXPANSION] = test_thermal_expansion;
    if (rDimension == 2) {
        (*p_prop)[THICKNESS] = test_thickness;
    }
    p_prop->SetValue(CONSTITUTIVE_LAW, CreateThermalLaw(rLawName));

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

/// Prescribes a state exercising mechanical and thermal contributions
/// simultaneously: non-zero displacement strain and a non-uniform nodal
/// temperature field.
void PrescribeMixedScenario(ModelPart& rModelPart, const std::size_t rDimension)
{
    KRATOS_TRY;

    Matrix displacement_gradient = ZeroMatrix(3, 3);
    array_1d<double, 3> translation;
    if (rDimension == 2) {
        displacement_gradient(0, 0) = 2.0e-3;
        displacement_gradient(0, 1) = 3.0e-4;
        displacement_gradient(1, 0) = 5.0e-4;
        displacement_gradient(1, 1) = -1.0e-3;
        translation[0] = 1.0e-2;
        translation[1] = -2.0e-2;
    } else {
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
    }

    std::size_t node_index = 0;
    for (auto& r_node : rModelPart.Nodes()) {
        array_1d<double, 3>& r_displacement = r_node.FastGetSolutionStepValue(DISPLACEMENT);
        noalias(r_displacement) = prod(displacement_gradient, r_node.GetInitialPosition());
        noalias(r_displacement) += translation;
        r_node.X() = r_node.X0() + r_displacement[0];
        r_node.Y() = r_node.Y0() + r_displacement[1];
        r_node.Z() = r_node.Z0() + r_displacement[2];

        r_node.FastGetSolutionStepValue(TEMPERATURE) =
            test_reference_temperature + 5.0 * static_cast<double>(node_index);

        array_1d<double, 3>& r_volume_acceleration = r_node.FastGetSolutionStepValue(VOLUME_ACCELERATION);
        r_volume_acceleration[0] = 0.0;
        r_volume_acceleration[1] = -9.81;
        r_volume_acceleration[2] = 0.0;

        ++node_index;
    }

    KRATOS_CATCH("");
}

/// Initializes an element through the standard lifecycle.
void InitializeVectorOutputElement(ModelPart& rModelPart)
{
    KRATOS_TRY;

    auto p_element = rModelPart.pGetElement(1);
    const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

    KRATOS_EXPECT_EQ(p_element->Check(r_process_info), 0);
    p_element->Initialize(r_process_info);
    p_element->InitializeSolutionStep(r_process_info);
    p_element->InitializeNonLinearIteration(r_process_info);

    KRATOS_CATCH("");
}

/// Responses queried through the standard integration-point mechanism.
struct VectorOutputResponse
{
    std::vector<Vector> strain_vectors;
    std::vector<Vector> cauchy_stress_vectors;
    std::vector<Vector> pk2_stress_vectors;
    std::vector<Vector> thermal_strain_vectors;
    std::vector<Vector> thermal_stress_vectors;
    std::vector<Vector> mechanical_stress_vectors;
};

void CalculateVectorOutputs(ModelPart& rModelPart, VectorOutputResponse& rResponse)
{
    KRATOS_TRY;

    auto p_element = rModelPart.pGetElement(1);
    const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

    p_element->CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR, rResponse.strain_vectors, r_process_info);
    p_element->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, rResponse.cauchy_stress_vectors, r_process_info);
    p_element->CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, rResponse.pk2_stress_vectors, r_process_info);
    p_element->CalculateOnIntegrationPoints(THERMAL_STRAIN_VECTOR, rResponse.thermal_strain_vectors, r_process_info);
    p_element->CalculateOnIntegrationPoints(THERMAL_STRESS_VECTOR, rResponse.thermal_stress_vectors, r_process_info);
    p_element->CalculateOnIntegrationPoints(MECHANICAL_STRESS_VECTOR, rResponse.mechanical_stress_vectors, r_process_info);

    KRATOS_CATCH("");
}

/// Analytical isotropic constitutive matrix for the given law.
Matrix AnalyticalConstitutiveMatrix(const std::string& rLawName)
{
    const double E = test_young_modulus;
    const double nu = test_poisson_ratio;
    if (rLawName == "ThermalLinearElastic2DPlaneStrain") {
        Matrix C(3, 3);
        C.clear();
        const double c11 = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
        const double c12 = c11 * nu / (1.0 - nu);
        const double c33 = c11 * (1.0 - 2.0 * nu) / (2.0 * (1.0 - nu));
        C(0, 0) = c11;
        C(1, 1) = c11;
        C(2, 2) = c33;
        C(0, 1) = c12;
        C(1, 0) = c12;
        return C;
    } else if (rLawName == "ThermalLinearElastic2DPlaneStress") {
        Matrix C(3, 3);
        C.clear();
        const double c11 = E / (1.0 - nu * nu);
        const double c12 = c11 * nu;
        const double c33 = c11 * (1.0 - nu) * 0.5;
        C(0, 0) = c11;
        C(1, 1) = c11;
        C(2, 2) = c33;
        C(0, 1) = c12;
        C(1, 0) = c12;
        return C;
    } else {
        Matrix C(6, 6);
        C.clear();
        const double c11 = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
        const double c12 = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        const double c33 = E / (2.0 * (1.0 + nu));
        for (std::size_t i = 0; i < 3; ++i) {
            C(i, i) = c11;
            for (std::size_t j = 0; j < 3; ++j) {
                if (i != j) C(i, j) = c12;
            }
        }
        for (std::size_t i = 3; i < 6; ++i) {
            C(i, i) = c33;
        }
        return C;
    }
}

/// Analytical thermal strain for a given delta_T.
Vector AnalyticalThermalStrain(const std::string& rLawName, const double rDeltaTemperature)
{
    Vector epsilon_th;
    if (rLawName == "ThermalLinearElastic2DPlaneStrain") {
        epsilon_th.resize(3);
        epsilon_th[0] = (1.0 + test_poisson_ratio) * test_thermal_expansion * rDeltaTemperature;
        epsilon_th[1] = epsilon_th[0];
        epsilon_th[2] = 0.0;
    } else if (rLawName == "ThermalLinearElastic2DPlaneStress") {
        epsilon_th.resize(3);
        epsilon_th[0] = test_thermal_expansion * rDeltaTemperature;
        epsilon_th[1] = epsilon_th[0];
        epsilon_th[2] = 0.0;
    } else {
        epsilon_th.resize(6);
        epsilon_th[0] = test_thermal_expansion * rDeltaTemperature;
        epsilon_th[1] = epsilon_th[0];
        epsilon_th[2] = epsilon_th[0];
        epsilon_th[3] = 0.0;
        epsilon_th[4] = 0.0;
        epsilon_th[5] = 0.0;
    }
    return epsilon_th;
}

/// Interpolates a nodal scalar at a Gauss point with the shape functions.
double InterpolateScalarAtGP(
    const Geometry<Node>& rGeometry,
    const Variable<double>& rVariable,
    const std::size_t rPointNumber,
    const GeometryData::IntegrationMethod rMethod)
{
    const auto& r_integration_points = rGeometry.IntegrationPoints(rMethod);
    Vector N(rGeometry.PointsNumber());
    rGeometry.ShapeFunctionsValues(N, r_integration_points[rPointNumber].Coordinates());
    double value = 0.0;
    for (std::size_t i = 0; i < rGeometry.PointsNumber(); ++i) {
        value += N[i] * rGeometry[i].FastGetSolutionStepValue(rVariable);
    }
    return value;
}

/// Tolerance associated with a reference value.
double ComponentTolerance(const double rReferenceValue, const double rReferenceScale)
{
    double tolerance = std::max(comparison_absolute_tolerance,
                                comparison_relative_tolerance * std::abs(rReferenceValue));
    if (rReferenceValue == 0.0) {
        tolerance = std::max(tolerance, machine_precision_allowance * rReferenceScale);
    }
    return tolerance;
}

double MaxAbsEntry(const Vector& rVector)
{
    double max_abs_entry = 0.0;
    for (std::size_t i = 0; i < rVector.size(); ++i) {
        max_abs_entry = std::max(max_abs_entry, std::abs(rVector(i)));
    }
    return max_abs_entry;
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

/// Verifies the three vector outputs for one configuration over a mixed
/// mechanical+thermal state.
void VerifyThermalVectorOutputs(
    const std::string& rLegacyElementName,
    const std::string& rCandidateElementName,
    const std::string& rLawName,
    const std::size_t rDimension,
    const std::string& rLabel)
{
    Model model;
    ModelPart& r_legacy = CreateVectorModelPart(
        model, "Legacy" + rLabel, rLegacyElementName, rLawName, rDimension, 0);
    ModelPart& r_candidate = CreateVectorModelPart(
        model, "Candidate" + rLabel, rCandidateElementName, rLawName, rDimension, 100);

    PrescribeMixedScenario(r_legacy, rDimension);
    PrescribeMixedScenario(r_candidate, rDimension);
    InitializeVectorOutputElement(r_legacy);
    InitializeVectorOutputElement(r_candidate);

    VectorOutputResponse legacy_response, candidate_response;
    CalculateVectorOutputs(r_legacy, legacy_response);
    CalculateVectorOutputs(r_candidate, candidate_response);

    auto p_candidate_element = r_candidate.pGetElement(1);
    const auto& r_geometry = p_candidate_element->GetGeometry();
    const GeometryData::IntegrationMethod integration_method = p_candidate_element->GetIntegrationMethod();

    const Matrix C = AnalyticalConstitutiveMatrix(rLawName);
    const std::size_t expected_size = (rDimension == 2) ? 3 : 6;

    KRATOS_EXPECT_EQ(candidate_response.thermal_strain_vectors[0].size(), expected_size);
    KRATOS_EXPECT_EQ(candidate_response.thermal_stress_vectors[0].size(), expected_size);
    KRATOS_EXPECT_EQ(candidate_response.mechanical_stress_vectors[0].size(), expected_size);
    KRATOS_EXPECT_EQ(candidate_response.cauchy_stress_vectors[0].size(), expected_size);

    const std::size_t number_of_gps = candidate_response.thermal_strain_vectors.size();
    for (std::size_t gp = 0; gp < number_of_gps; ++gp) {
        const double delta_temperature =
            InterpolateScalarAtGP(r_geometry, TEMPERATURE, gp, integration_method) -
            InterpolateScalarAtGP(r_geometry, NODAL_REFERENCE_TEMPERATURE, gp, integration_method);
        const Vector epsilon_th = AnalyticalThermalStrain(rLawName, delta_temperature);
        const Vector& r_epsilon = candidate_response.strain_vectors[gp];

        // THERMAL_STRAIN_VECTOR == epsilon_th (intentional fix of the legacy bug).
        ExpectVectorComponentsNear(
            candidate_response.thermal_strain_vectors[gp], epsilon_th,
            rLabel + " thermal strain == analytical epsilon_th");

        // THERMAL_STRESS_VECTOR == C*epsilon_th.
        const Vector expected_thermal_stress = prod(C, epsilon_th);
        ExpectVectorComponentsNear(
            candidate_response.thermal_stress_vectors[gp], expected_thermal_stress,
            rLabel + " thermal stress == C*epsilon_th");

        // MECHANICAL_STRESS_VECTOR == C*epsilon.
        const Vector expected_mechanical_stress = prod(C, r_epsilon);
        ExpectVectorComponentsNear(
            candidate_response.mechanical_stress_vectors[gp], expected_mechanical_stress,
            rLabel + " mechanical stress == C*epsilon");

        // Cross-check with the legacy element outputs (correct for these two).
        ExpectVectorComponentsNear(
            candidate_response.thermal_stress_vectors[gp], legacy_response.thermal_stress_vectors[gp],
            rLabel + " thermal stress legacy cross-check");
        ExpectVectorComponentsNear(
            candidate_response.mechanical_stress_vectors[gp], legacy_response.mechanical_stress_vectors[gp],
            rLabel + " mechanical stress legacy cross-check");

        // Decomposition identity: total == mechanical - thermal.
        Vector decomposed = candidate_response.mechanical_stress_vectors[gp];
        decomposed -= candidate_response.thermal_stress_vectors[gp];
        ExpectVectorComponentsNear(
            candidate_response.cauchy_stress_vectors[gp], decomposed,
            rLabel + " cauchy == mechanical - thermal");
        ExpectVectorComponentsNear(
            candidate_response.pk2_stress_vectors[gp], decomposed,
            rLabel + " pk2 == mechanical - thermal");
    }
}

} // namespace

//************************************************************************************
// Per-configuration vector-output verification
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalVectorOutputs_Triangle2D3_PlaneStrain, KratosDamFastSuite)
{
    VerifyThermalVectorOutputs(
        "SmallDisplacementThermoMechanicElement2D3N", "SmallDisplacementElement2D3N",
        "ThermalLinearElastic2DPlaneStrain", 2, "Triangle2D3PlaneStrain");
}

KRATOS_TEST_CASE_IN_SUITE(ThermalVectorOutputs_Quadrilateral2D4_PlaneStress, KratosDamFastSuite)
{
    VerifyThermalVectorOutputs(
        "SmallDisplacementThermoMechanicElement2D4N", "SmallDisplacementElement2D4N",
        "ThermalLinearElastic2DPlaneStress", 2, "Quadrilateral2D4PlaneStress");
}

KRATOS_TEST_CASE_IN_SUITE(ThermalVectorOutputs_Tetrahedra3D4, KratosDamFastSuite)
{
    VerifyThermalVectorOutputs(
        "SmallDisplacementThermoMechanicElement3D4N", "SmallDisplacementElement3D4N",
        "ThermalLinearElastic3DLaw", 3, "Tetrahedra3D4");
}

KRATOS_TEST_CASE_IN_SUITE(ThermalVectorOutputs_Hexahedra3D8, KratosDamFastSuite)
{
    VerifyThermalVectorOutputs(
        "SmallDisplacementThermoMechanicElement3D8N", "SmallDisplacementElement3D8N",
        "ThermalLinearElastic3DLaw", 3, "Hexahedra3D8");
}

//************************************************************************************
// Intentional THERMAL_STRAIN_VECTOR bug fix regression
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalVectorOutputs_ThermalStrainBugFix, KratosDamFastSuite)
{
    // The legacy element returns the mechanical (element) strain for
    // THERMAL_STRAIN_VECTOR; the new constitutive-law output returns the actual
    // thermal strain. This is an intentional bug fix, documented here.
    Model model;
    ModelPart& r_legacy = CreateVectorModelPart(
        model, "BugFixLegacy", "SmallDisplacementThermoMechanicElement3D8N",
        "ThermalLinearElastic3DLaw", 3, 0);
    ModelPart& r_candidate = CreateVectorModelPart(
        model, "BugFixCandidate", "SmallDisplacementElement3D8N",
        "ThermalLinearElastic3DLaw", 3, 100);

    PrescribeMixedScenario(r_legacy, 3);
    PrescribeMixedScenario(r_candidate, 3);
    InitializeVectorOutputElement(r_legacy);
    InitializeVectorOutputElement(r_candidate);

    std::vector<Vector> candidate_thermal_strain, legacy_thermal_strain;
    r_candidate.pGetElement(1)->CalculateOnIntegrationPoints(
        THERMAL_STRAIN_VECTOR, candidate_thermal_strain, r_candidate.GetProcessInfo());
    r_legacy.pGetElement(1)->CalculateOnIntegrationPoints(
        THERMAL_STRAIN_VECTOR, legacy_thermal_strain, r_legacy.GetProcessInfo());

    auto p_candidate_element = r_candidate.pGetElement(1);
    const auto& r_geometry = p_candidate_element->GetGeometry();
    const GeometryData::IntegrationMethod integration_method = p_candidate_element->GetIntegrationMethod();
    const double delta_temperature =
        InterpolateScalarAtGP(r_geometry, TEMPERATURE, 0, integration_method) -
        InterpolateScalarAtGP(r_geometry, NODAL_REFERENCE_TEMPERATURE, 0, integration_method);
    const Vector expected_epsilon_th = AnalyticalThermalStrain("ThermalLinearElastic3DLaw", delta_temperature);

    // The candidate returns the actual thermal strain.
    ExpectVectorComponentsNear(
        candidate_thermal_strain[0], expected_epsilon_th,
        "thermal strain bug fix: candidate == analytical epsilon_th");
    // The legacy element returns the mechanical strain (documented legacy bug),
    // which differs from the thermal strain in this mixed scenario.
    KRATOS_EXPECT_TRUE(legacy_thermal_strain[0][0] != expected_epsilon_th[0]);
    std::cout << "[CHARACTERIZATION] legacy THERMAL_STRAIN_VECTOR[0] = "
              << legacy_thermal_strain[0][0]
              << " (mechanical strain, documented legacy bug) vs candidate epsilon_th[0] = "
              << expected_epsilon_th[0] << std::endl;
}

//************************************************************************************
// Zero-component limiting behavior
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalVectorOutputs_ZeroComponents, KratosDamFastSuite)
{
    // delta_T = 0 with non-zero displacement: thermal outputs vanish.
    {
        Model model;
        ModelPart& r_model_part = CreateVectorModelPart(
            model, "ZeroThermal", "SmallDisplacementElement3D8N", "ThermalLinearElastic3DLaw", 3, 0);
        PrescribeMixedScenario(r_model_part, 3);
        // Zero the thermal increment (TEMPERATURE = reference everywhere).
        for (auto& r_node : r_model_part.Nodes()) {
            r_node.FastGetSolutionStepValue(TEMPERATURE) = test_reference_temperature;
        }
        InitializeVectorOutputElement(r_model_part);
        VectorOutputResponse response;
        CalculateVectorOutputs(r_model_part, response);
        Vector zero(6);
        noalias(zero) = ZeroVector(6);
        for (std::size_t gp = 0; gp < response.thermal_strain_vectors.size(); ++gp) {
            ExpectVectorComponentsNear(response.thermal_strain_vectors[gp], zero, "zero delta_T thermal strain");
            ExpectVectorComponentsNear(response.thermal_stress_vectors[gp], zero, "zero delta_T thermal stress");
            KRATOS_EXPECT_TRUE(MaxAbsEntry(response.mechanical_stress_vectors[gp]) > 0.0);
        }
    }
    // displacement strain = 0 with non-zero delta_T: mechanical stress vanishes.
    {
        Model model;
        ModelPart& r_model_part = CreateVectorModelPart(
            model, "ZeroMechanical", "SmallDisplacementElement3D8N", "ThermalLinearElastic3DLaw", 3, 0);
        // Zero displacement, uniform temperature increment.
        for (auto& r_node : r_model_part.Nodes()) {
            r_node.FastGetSolutionStepValue(DISPLACEMENT) = ZeroVector(3);
            r_node.FastGetSolutionStepValue(TEMPERATURE) = test_reference_temperature + 25.0;
        }
        InitializeVectorOutputElement(r_model_part);
        VectorOutputResponse response;
        CalculateVectorOutputs(r_model_part, response);
        Vector zero(6);
        noalias(zero) = ZeroVector(6);
        for (std::size_t gp = 0; gp < response.mechanical_stress_vectors.size(); ++gp) {
            ExpectVectorComponentsNear(response.mechanical_stress_vectors[gp], zero, "zero strain mechanical stress");
            KRATOS_EXPECT_TRUE(MaxAbsEntry(response.thermal_stress_vectors[gp]) > 0.0);
        }
    }
    // Both zero: all three outputs vanish.
    {
        Model model;
        ModelPart& r_model_part = CreateVectorModelPart(
            model, "ZeroBoth", "SmallDisplacementElement3D8N", "ThermalLinearElastic3DLaw", 3, 0);
        for (auto& r_node : r_model_part.Nodes()) {
            r_node.FastGetSolutionStepValue(DISPLACEMENT) = ZeroVector(3);
            r_node.FastGetSolutionStepValue(TEMPERATURE) = test_reference_temperature;
        }
        InitializeVectorOutputElement(r_model_part);
        VectorOutputResponse response;
        CalculateVectorOutputs(r_model_part, response);
        Vector zero(6);
        noalias(zero) = ZeroVector(6);
        for (std::size_t gp = 0; gp < response.thermal_strain_vectors.size(); ++gp) {
            ExpectVectorComponentsNear(response.thermal_strain_vectors[gp], zero, "both zero thermal strain");
            ExpectVectorComponentsNear(response.thermal_stress_vectors[gp], zero, "both zero thermal stress");
            ExpectVectorComponentsNear(response.mechanical_stress_vectors[gp], zero, "both zero mechanical stress");
        }
    }
}

//************************************************************************************
// Side-effect / read-only requirement
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalVectorOutputs_SideEffect, KratosDamFastSuite)
{
    // Requesting the vector outputs must not alter the subsequent equilibrium
    // response or the constitutive state.
    for (const char* p_element_name : {"SmallDisplacementElement2D4N",
                                       "SmallDisplacementElement3D8N"}) {
        const std::string r_element_name(p_element_name);
        const std::size_t dimension = (r_element_name.find("2D") != std::string::npos) ? 2 : 3;
        const std::string r_law_name =
            (dimension == 2) ? "ThermalLinearElastic2DPlaneStress" : "ThermalLinearElastic3DLaw";

        Model model;
        ModelPart& r_model_part = CreateVectorModelPart(
            model, "SideEffect" + r_element_name, r_element_name, r_law_name, dimension, 0);
        PrescribeMixedScenario(r_model_part, dimension);
        InitializeVectorOutputElement(r_model_part);
        auto p_element = r_model_part.pGetElement(1);
        const ProcessInfo& r_pi = r_model_part.GetProcessInfo();

        Matrix lhs_before, lhs_after;
        Vector rhs_before, rhs_after;
        std::vector<Vector> cauchy_before, cauchy_after, pk2_before, pk2_after;
        p_element->CalculateLocalSystem(lhs_before, rhs_before, r_pi);
        p_element->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, cauchy_before, r_pi);
        p_element->CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, pk2_before, r_pi);

        // Repeatedly request all three vector outputs.
        for (int repeat = 0; repeat < 3; ++repeat) {
            std::vector<Vector> dummy;
            p_element->CalculateOnIntegrationPoints(THERMAL_STRAIN_VECTOR, dummy, r_pi);
            p_element->CalculateOnIntegrationPoints(THERMAL_STRESS_VECTOR, dummy, r_pi);
            p_element->CalculateOnIntegrationPoints(MECHANICAL_STRESS_VECTOR, dummy, r_pi);
        }

        p_element->CalculateLocalSystem(lhs_after, rhs_after, r_pi);
        p_element->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, cauchy_after, r_pi);
        p_element->CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, pk2_after, r_pi);

        ExpectMatrixComponentsNear(lhs_after, lhs_before, r_element_name + " side effect LHS");
        ExpectVectorComponentsNear(rhs_after, rhs_before, r_element_name + " side effect RHS");
        ExpectIntegrationPointVectorsNear(cauchy_after, cauchy_before, r_element_name + " side effect Cauchy");
        ExpectIntegrationPointVectorsNear(pk2_after, pk2_before, r_element_name + " side effect PK2");
    }
}

} // namespace Testing
} // namespace Kratos
