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

// Phase 3B characterization: the thermo-mechanical tensor outputs
//   THERMAL_STRAIN_TENSOR, THERMAL_STRESS_TENSOR, MECHANICAL_STRESS_TENSOR
// are provided by ThermalLinearElastic3DLaw through the standard
// BaseSolidElement::CalculateOnIntegrationPoints infrastructure, as the tensor
// representations of the phase-3A vector outputs:
//   THERMAL_STRAIN_TENSOR    = Tensor(THERMAL_STRAIN_VECTOR)
//   THERMAL_STRESS_TENSOR    = Tensor(THERMAL_STRESS_VECTOR)
//   MECHANICAL_STRESS_TENSOR = Tensor(MECHANICAL_STRESS_VECTOR)
// using the standard MathUtils Voigt-to-tensor conversions.
//
// Dimensions: 2x2 for the 2D laws, 3x3 for the 3D law. The strain tensor
// applies the 1/2 factor to engineering-shear components; the stress tensor
// copies the shear components without the factor.
//
// The legacy THERMAL_STRAIN_TENSOR is derived from the known-broken legacy
// THERMAL_STRAIN_VECTOR and is therefore not used as a reference; the candidate
// output is validated against the analytical thermal strain instead.

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
#include "utilities/math_utils.h"

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
ModelPart& CreateTensorModelPart(
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
void InitializeTensorOutputElement(ModelPart& rModelPart)
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

void ExpectMatrixSymmetric(const Matrix& rMatrix, const std::string& rWhat)
{
    ASSERT_EQ(rMatrix.size1(), rMatrix.size2()) << "Not square comparing " << rWhat;
    for (std::size_t i = 0; i < rMatrix.size1(); ++i) {
        for (std::size_t j = 0; j < rMatrix.size1(); ++j) {
            KRATOS_EXPECT_NEAR(rMatrix(i, j), rMatrix(j, i), comparison_absolute_tolerance);
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
// Per-configuration tensor-output verification
//************************************************************************************

namespace
{

/// Verifies the three tensor outputs for one configuration over a mixed
/// mechanical+thermal state.
void VerifyTensorOutputs(
    const std::string& rLegacyElementName,
    const std::string& rCandidateElementName,
    const std::string& rLawName,
    const std::size_t rDimension,
    const std::string& rLabel)
{
    Model model;
    ModelPart& r_legacy = CreateTensorModelPart(
        model, "Legacy" + rLabel, rLegacyElementName, rLawName, rDimension, 0);
    ModelPart& r_candidate = CreateTensorModelPart(
        model, "Candidate" + rLabel, rCandidateElementName, rLawName, rDimension, 100);

    PrescribeMixedScenario(r_legacy, rDimension);
    PrescribeMixedScenario(r_candidate, rDimension);
    InitializeTensorOutputElement(r_legacy);
    InitializeTensorOutputElement(r_candidate);

    auto p_candidate = r_candidate.pGetElement(1);
    auto p_legacy = r_legacy.pGetElement(1);
    const ProcessInfo& r_candidate_pi = r_candidate.GetProcessInfo();
    const ProcessInfo& r_legacy_pi = r_legacy.GetProcessInfo();

    // Candidate outputs (tensors + vectors).
    std::vector<Matrix> cauchy_tensor, pk2_tensor, thermal_strain_tensor, thermal_stress_tensor, mechanical_stress_tensor;
    std::vector<Vector> thermal_strain_vector, thermal_stress_vector, mechanical_stress_vector, total_strain_vector;
    p_candidate->CalculateOnIntegrationPoints(CAUCHY_STRESS_TENSOR, cauchy_tensor, r_candidate_pi);
    p_candidate->CalculateOnIntegrationPoints(PK2_STRESS_TENSOR, pk2_tensor, r_candidate_pi);
    p_candidate->CalculateOnIntegrationPoints(THERMAL_STRAIN_TENSOR, thermal_strain_tensor, r_candidate_pi);
    p_candidate->CalculateOnIntegrationPoints(THERMAL_STRESS_TENSOR, thermal_stress_tensor, r_candidate_pi);
    p_candidate->CalculateOnIntegrationPoints(MECHANICAL_STRESS_TENSOR, mechanical_stress_tensor, r_candidate_pi);
    p_candidate->CalculateOnIntegrationPoints(THERMAL_STRAIN_VECTOR, thermal_strain_vector, r_candidate_pi);
    p_candidate->CalculateOnIntegrationPoints(THERMAL_STRESS_VECTOR, thermal_stress_vector, r_candidate_pi);
    p_candidate->CalculateOnIntegrationPoints(MECHANICAL_STRESS_VECTOR, mechanical_stress_vector, r_candidate_pi);
    p_candidate->CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR, total_strain_vector, r_candidate_pi);

    // Legacy outputs (correct for the two stress tensors only).
    std::vector<Matrix> legacy_thermal_stress_tensor, legacy_mechanical_stress_tensor;
    p_legacy->CalculateOnIntegrationPoints(THERMAL_STRESS_TENSOR, legacy_thermal_stress_tensor, r_legacy_pi);
    p_legacy->CalculateOnIntegrationPoints(MECHANICAL_STRESS_TENSOR, legacy_mechanical_stress_tensor, r_legacy_pi);

    const auto& r_geometry = p_candidate->GetGeometry();
    const GeometryData::IntegrationMethod integration_method = p_candidate->GetIntegrationMethod();
    const Matrix C = AnalyticalConstitutiveMatrix(rLawName);
    const std::size_t expected_size = (rDimension == 2) ? 2 : 3;

    KRATOS_EXPECT_EQ(thermal_strain_tensor[0].size1(), expected_size);
    KRATOS_EXPECT_EQ(thermal_stress_tensor[0].size1(), expected_size);
    KRATOS_EXPECT_EQ(mechanical_stress_tensor[0].size1(), expected_size);

    const std::size_t number_of_gps = thermal_strain_tensor.size();
    for (std::size_t gp = 0; gp < number_of_gps; ++gp) {
        const double delta_temperature =
            InterpolateScalarAtGP(r_geometry, TEMPERATURE, gp, integration_method) -
            InterpolateScalarAtGP(r_geometry, NODAL_REFERENCE_TEMPERATURE, gp, integration_method);
        const Vector epsilon_th = AnalyticalThermalStrain(rLawName, delta_temperature);
        const Vector& r_epsilon = total_strain_vector[gp];

        const Vector expected_thermal_stress = prod(C, epsilon_th);
        const Vector expected_mechanical_stress = prod(C, r_epsilon);
        const Matrix expected_thermal_strain_tensor = MathUtils<double>::StrainVectorToTensor(epsilon_th);
        const Matrix expected_thermal_stress_tensor = MathUtils<double>::StressVectorToTensor(expected_thermal_stress);
        const Matrix expected_mechanical_stress_tensor = MathUtils<double>::StressVectorToTensor(expected_mechanical_stress);

        // Analytical tensor representations.
        ExpectMatrixComponentsNear(
            thermal_strain_tensor[gp], expected_thermal_strain_tensor,
            rLabel + " thermal strain tensor == analytical tensor");
        ExpectMatrixComponentsNear(
            thermal_stress_tensor[gp], expected_thermal_stress_tensor,
            rLabel + " thermal stress tensor == analytical tensor");
        ExpectMatrixComponentsNear(
            mechanical_stress_tensor[gp], expected_mechanical_stress_tensor,
            rLabel + " mechanical stress tensor == analytical tensor");

        // Tensor exactly consistent with the corresponding vector output.
        ExpectMatrixComponentsNear(
            thermal_strain_tensor[gp], MathUtils<double>::StrainVectorToTensor(thermal_strain_vector[gp]),
            rLabel + " thermal strain tensor == Tensor(vector)");
        ExpectMatrixComponentsNear(
            thermal_stress_tensor[gp], MathUtils<double>::StressVectorToTensor(thermal_stress_vector[gp]),
            rLabel + " thermal stress tensor == Tensor(vector)");
        ExpectMatrixComponentsNear(
            mechanical_stress_tensor[gp], MathUtils<double>::StressVectorToTensor(mechanical_stress_vector[gp]),
            rLabel + " mechanical stress tensor == Tensor(vector)");

        // Legacy cross-check (correct stress tensors).
        ExpectMatrixComponentsNear(
            thermal_stress_tensor[gp], legacy_thermal_stress_tensor[gp],
            rLabel + " thermal stress tensor legacy cross-check");
        ExpectMatrixComponentsNear(
            mechanical_stress_tensor[gp], legacy_mechanical_stress_tensor[gp],
            rLabel + " mechanical stress tensor legacy cross-check");

        // Decomposition: total == mechanical - thermal.
        Matrix decomposed = mechanical_stress_tensor[gp];
        decomposed -= thermal_stress_tensor[gp];
        ExpectMatrixComponentsNear(
            cauchy_tensor[gp], decomposed, rLabel + " cauchy tensor == mechanical - thermal");
        ExpectMatrixComponentsNear(
            pk2_tensor[gp], decomposed, rLabel + " pk2 tensor == mechanical - thermal");

        // Symmetry.
        ExpectMatrixSymmetric(thermal_strain_tensor[gp], rLabel + " thermal strain tensor symmetry");
        ExpectMatrixSymmetric(thermal_stress_tensor[gp], rLabel + " thermal stress tensor symmetry");
        ExpectMatrixSymmetric(mechanical_stress_tensor[gp], rLabel + " mechanical stress tensor symmetry");

        // Shear convention: strain tensor uses half the engineering shear; the
        // stress tensor copies the shear without the factor.
        if (expected_size == 3) {
            KRATOS_EXPECT_NEAR(thermal_strain_tensor[gp](0, 1), 0.5 * thermal_strain_vector[gp][3], 1.0e-12);
            KRATOS_EXPECT_NEAR(thermal_stress_tensor[gp](0, 1), thermal_stress_vector[gp][3], 1.0e-12);
        } else {
            KRATOS_EXPECT_NEAR(thermal_strain_tensor[gp](0, 1), 0.5 * thermal_strain_vector[gp][2], 1.0e-12);
            KRATOS_EXPECT_NEAR(thermal_stress_tensor[gp](0, 1), thermal_stress_vector[gp][2], 1.0e-12);
        }
    }
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(ThermalTensorOutputs_Triangle2D3_PlaneStrain, KratosDamFastSuite)
{
    VerifyTensorOutputs(
        "SmallDisplacementThermoMechanicElement2D3N", "SmallDisplacementElement2D3N",
        "ThermalLinearElastic2DPlaneStrain", 2, "Triangle2D3PlaneStrain");
}

KRATOS_TEST_CASE_IN_SUITE(ThermalTensorOutputs_Quadrilateral2D4_PlaneStress, KratosDamFastSuite)
{
    VerifyTensorOutputs(
        "SmallDisplacementThermoMechanicElement2D4N", "SmallDisplacementElement2D4N",
        "ThermalLinearElastic2DPlaneStress", 2, "Quadrilateral2D4PlaneStress");
}

KRATOS_TEST_CASE_IN_SUITE(ThermalTensorOutputs_Tetrahedra3D4, KratosDamFastSuite)
{
    VerifyTensorOutputs(
        "SmallDisplacementThermoMechanicElement3D4N", "SmallDisplacementElement3D4N",
        "ThermalLinearElastic3DLaw", 3, "Tetrahedra3D4");
}

KRATOS_TEST_CASE_IN_SUITE(ThermalTensorOutputs_Hexahedra3D8, KratosDamFastSuite)
{
    VerifyTensorOutputs(
        "SmallDisplacementThermoMechanicElement3D8N", "SmallDisplacementElement3D8N",
        "ThermalLinearElastic3DLaw", 3, "Hexahedra3D8");
}

//************************************************************************************
// Intentional THERMAL_STRAIN_TENSOR bug fix regression
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalTensorOutputs_ThermalStrainTensorBugFix, KratosDamFastSuite)
{
    // Tensor-level continuation of the phase-3A/5B.2 bug fix: since Phase 5B.2
    // BOTH the candidate and the legacy element derive THERMAL_STRAIN_TENSOR
    // from the physically correct thermal-strain vector (the legacy element no
    // longer returns the mechanical-strain tensor).
    Model model;
    ModelPart& r_legacy = CreateTensorModelPart(
        model, "TensorBugFixLegacy", "SmallDisplacementThermoMechanicElement3D8N",
        "ThermalLinearElastic3DLaw", 3, 0);
    ModelPart& r_candidate = CreateTensorModelPart(
        model, "TensorBugFixCandidate", "SmallDisplacementElement3D8N",
        "ThermalLinearElastic3DLaw", 3, 100);

    PrescribeMixedScenario(r_legacy, 3);
    PrescribeMixedScenario(r_candidate, 3);
    InitializeTensorOutputElement(r_legacy);
    InitializeTensorOutputElement(r_candidate);

    std::vector<Matrix> candidate_thermal_strain_tensor, legacy_thermal_strain_tensor;
    r_candidate.pGetElement(1)->CalculateOnIntegrationPoints(
        THERMAL_STRAIN_TENSOR, candidate_thermal_strain_tensor, r_candidate.GetProcessInfo());
    r_legacy.pGetElement(1)->CalculateOnIntegrationPoints(
        THERMAL_STRAIN_TENSOR, legacy_thermal_strain_tensor, r_legacy.GetProcessInfo());

    auto p_candidate = r_candidate.pGetElement(1);
    const auto& r_geometry = p_candidate->GetGeometry();
    const GeometryData::IntegrationMethod integration_method = p_candidate->GetIntegrationMethod();
    const double delta_temperature =
        InterpolateScalarAtGP(r_geometry, TEMPERATURE, 0, integration_method) -
        InterpolateScalarAtGP(r_geometry, NODAL_REFERENCE_TEMPERATURE, 0, integration_method);
    const Vector expected_epsilon_th = AnalyticalThermalStrain("ThermalLinearElastic3DLaw", delta_temperature);
    const Matrix expected_tensor = MathUtils<double>::StrainVectorToTensor(expected_epsilon_th);

    ExpectMatrixComponentsNear(
        candidate_thermal_strain_tensor[0], expected_tensor,
        "thermal strain tensor bug fix: candidate == analytical tensor");
    ExpectMatrixComponentsNear(
        legacy_thermal_strain_tensor[0], expected_tensor,
        "thermal strain tensor bug fix: legacy == analytical tensor");
    std::cout << "[CHARACTERIZATION] legacy THERMAL_STRAIN_TENSOR[0](0,0) = "
              << legacy_thermal_strain_tensor[0](0, 0)
              << " == analytical epsilon_th tensor (0,0) = "
              << expected_tensor(0, 0)
              << " (intentional bug fix applied to the legacy element)" << std::endl;
}

//************************************************************************************
// Zero-component limiting behavior
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalTensorOutputs_ZeroComponents, KratosDamFastSuite)
{
    // delta_T = 0 with non-zero displacement: thermal tensors vanish.
    {
        Model model;
        ModelPart& r_model_part = CreateTensorModelPart(
            model, "TensorZeroThermal", "SmallDisplacementElement3D8N", "ThermalLinearElastic3DLaw", 3, 0);
        PrescribeMixedScenario(r_model_part, 3);
        for (auto& r_node : r_model_part.Nodes()) {
            r_node.FastGetSolutionStepValue(TEMPERATURE) = test_reference_temperature;
        }
        InitializeTensorOutputElement(r_model_part);
        std::vector<Matrix> thermal_strain_tensor, thermal_stress_tensor, mechanical_stress_tensor;
        auto p_element = r_model_part.pGetElement(1);
        const ProcessInfo& r_pi = r_model_part.GetProcessInfo();
        p_element->CalculateOnIntegrationPoints(THERMAL_STRAIN_TENSOR, thermal_strain_tensor, r_pi);
        p_element->CalculateOnIntegrationPoints(THERMAL_STRESS_TENSOR, thermal_stress_tensor, r_pi);
        p_element->CalculateOnIntegrationPoints(MECHANICAL_STRESS_TENSOR, mechanical_stress_tensor, r_pi);
        const Matrix zero = ZeroMatrix(3, 3);
        for (std::size_t gp = 0; gp < thermal_strain_tensor.size(); ++gp) {
            ExpectMatrixComponentsNear(thermal_strain_tensor[gp], zero, "zero delta_T thermal strain tensor");
            ExpectMatrixComponentsNear(thermal_stress_tensor[gp], zero, "zero delta_T thermal stress tensor");
            KRATOS_EXPECT_TRUE(MaxAbsEntry(mechanical_stress_tensor[gp]) > 0.0);
        }
    }
    // displacement strain = 0 with non-zero delta_T: mechanical tensor vanishes.
    {
        Model model;
        ModelPart& r_model_part = CreateTensorModelPart(
            model, "TensorZeroMechanical", "SmallDisplacementElement3D8N", "ThermalLinearElastic3DLaw", 3, 0);
        for (auto& r_node : r_model_part.Nodes()) {
            r_node.FastGetSolutionStepValue(DISPLACEMENT) = ZeroVector(3);
            r_node.FastGetSolutionStepValue(TEMPERATURE) = test_reference_temperature + 25.0;
        }
        InitializeTensorOutputElement(r_model_part);
        std::vector<Matrix> mechanical_stress_tensor, thermal_stress_tensor;
        auto p_element = r_model_part.pGetElement(1);
        const ProcessInfo& r_pi = r_model_part.GetProcessInfo();
        p_element->CalculateOnIntegrationPoints(MECHANICAL_STRESS_TENSOR, mechanical_stress_tensor, r_pi);
        p_element->CalculateOnIntegrationPoints(THERMAL_STRESS_TENSOR, thermal_stress_tensor, r_pi);
        const Matrix zero = ZeroMatrix(3, 3);
        for (std::size_t gp = 0; gp < mechanical_stress_tensor.size(); ++gp) {
            ExpectMatrixComponentsNear(mechanical_stress_tensor[gp], zero, "zero strain mechanical stress tensor");
            KRATOS_EXPECT_TRUE(MaxAbsEntry(thermal_stress_tensor[gp]) > 0.0);
        }
    }
    // Both zero: all three tensors vanish.
    {
        Model model;
        ModelPart& r_model_part = CreateTensorModelPart(
            model, "TensorZeroBoth", "SmallDisplacementElement3D8N", "ThermalLinearElastic3DLaw", 3, 0);
        for (auto& r_node : r_model_part.Nodes()) {
            r_node.FastGetSolutionStepValue(DISPLACEMENT) = ZeroVector(3);
            r_node.FastGetSolutionStepValue(TEMPERATURE) = test_reference_temperature;
        }
        InitializeTensorOutputElement(r_model_part);
        std::vector<Matrix> thermal_strain_tensor, thermal_stress_tensor, mechanical_stress_tensor;
        auto p_element = r_model_part.pGetElement(1);
        const ProcessInfo& r_pi = r_model_part.GetProcessInfo();
        p_element->CalculateOnIntegrationPoints(THERMAL_STRAIN_TENSOR, thermal_strain_tensor, r_pi);
        p_element->CalculateOnIntegrationPoints(THERMAL_STRESS_TENSOR, thermal_stress_tensor, r_pi);
        p_element->CalculateOnIntegrationPoints(MECHANICAL_STRESS_TENSOR, mechanical_stress_tensor, r_pi);
        const Matrix zero = ZeroMatrix(3, 3);
        for (std::size_t gp = 0; gp < thermal_strain_tensor.size(); ++gp) {
            ExpectMatrixComponentsNear(thermal_strain_tensor[gp], zero, "both zero thermal strain tensor");
            ExpectMatrixComponentsNear(thermal_stress_tensor[gp], zero, "both zero thermal stress tensor");
            ExpectMatrixComponentsNear(mechanical_stress_tensor[gp], zero, "both zero mechanical stress tensor");
        }
    }
}

//************************************************************************************
// Side-effect / read-only requirement
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(ThermalTensorOutputs_SideEffect, KratosDamFastSuite)
{
    // Requesting the tensor outputs must not alter the subsequent equilibrium
    // response, the constitutive state, or the vector outputs.
    for (const char* p_element_name : {"SmallDisplacementElement2D4N",
                                       "SmallDisplacementElement3D8N"}) {
        const std::string r_element_name(p_element_name);
        const std::size_t dimension = (r_element_name.find("2D") != std::string::npos) ? 2 : 3;
        const std::string r_law_name =
            (dimension == 2) ? "ThermalLinearElastic2DPlaneStress" : "ThermalLinearElastic3DLaw";

        Model model;
        ModelPart& r_model_part = CreateTensorModelPart(
            model, "TensorSideEffect" + r_element_name, r_element_name, r_law_name, dimension, 0);
        PrescribeMixedScenario(r_model_part, dimension);
        InitializeTensorOutputElement(r_model_part);
        auto p_element = r_model_part.pGetElement(1);
        const ProcessInfo& r_pi = r_model_part.GetProcessInfo();

        Matrix lhs_before, lhs_after;
        Vector rhs_before, rhs_after;
        std::vector<Vector> cauchy_before, cauchy_after, pk2_before, pk2_after;
        std::vector<Vector> thermal_strain_before, thermal_strain_after;
        p_element->CalculateLocalSystem(lhs_before, rhs_before, r_pi);
        p_element->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, cauchy_before, r_pi);
        p_element->CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, pk2_before, r_pi);
        p_element->CalculateOnIntegrationPoints(THERMAL_STRAIN_VECTOR, thermal_strain_before, r_pi);

        for (int repeat = 0; repeat < 3; ++repeat) {
            std::vector<Matrix> dummy;
            p_element->CalculateOnIntegrationPoints(THERMAL_STRAIN_TENSOR, dummy, r_pi);
            p_element->CalculateOnIntegrationPoints(THERMAL_STRESS_TENSOR, dummy, r_pi);
            p_element->CalculateOnIntegrationPoints(MECHANICAL_STRESS_TENSOR, dummy, r_pi);
        }

        p_element->CalculateLocalSystem(lhs_after, rhs_after, r_pi);
        p_element->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, cauchy_after, r_pi);
        p_element->CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, pk2_after, r_pi);
        p_element->CalculateOnIntegrationPoints(THERMAL_STRAIN_VECTOR, thermal_strain_after, r_pi);

        ExpectMatrixComponentsNear(lhs_after, lhs_before, r_element_name + " tensor side effect LHS");
        ExpectVectorComponentsNear(rhs_after, rhs_before, r_element_name + " tensor side effect RHS");
        ExpectIntegrationPointVectorsNear(cauchy_after, cauchy_before, r_element_name + " tensor side effect Cauchy");
        ExpectIntegrationPointVectorsNear(pk2_after, pk2_before, r_element_name + " tensor side effect PK2");
        ExpectIntegrationPointVectorsNear(thermal_strain_after, thermal_strain_before, r_element_name + " tensor side effect vector output");
    }
}

} // namespace Testing
} // namespace Kratos
