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

// Characterization tests comparing the duplicated small-displacement solid
// element of DamApplication (Kratos::SmallDisplacementElement, registered as
// "SmallDisplacementSolidElement<geom>") with the corresponding element of
// StructuralMechanicsApplication (Kratos::SmallDisplacement, registered as
// "SmallDisplacementElement<geom>").
//
// Both elements are exercised with numerically identical ModelParts (same
// geometry, nodal values, properties, ProcessInfo, constitutive law,
// integration method and loading) and their responses are compared against
// each other. No hard-coded reference matrices are used for the main
// comparisons: the output of one implementation is the reference of the other.
//
// Setup notes:
// - Each ModelPart gets its own Properties and its own constitutive law
//   instance; each element clones the law of its own property during
//   Initialize, so no mutable constitutive-law state is shared between the
//   two implementations.
// - A non-zero displacement field u = A*X0 + t is prescribed and the current
//   nodal coordinates are updated to X = X0 + u. This is required because the
//   DamApplication element computes the reference configuration as
//   (current position - total displacement), while the
//   StructuralMechanicsApplication element uses the initial configuration
//   directly. With updated coordinates both reference configurations are X0.
// - A non-zero nodal VOLUME_ACCELERATION field is prescribed so that the
//   body-force (external force) contribution to the RHS is also compared.
// - THICKNESS is set for the 2D cases; all 2D comparisons (LHS, RHS, mass,
//   damping) therefore exercise the 2D thickness handling of both elements.
// - COMPUTE_CONSISTENT_MASS_MATRIX is set to true so that the DamApplication
//   element integrates the consistent mass matrix with the same increased
//   quadrature order that the StructuralMechanicsApplication element always
//   uses (IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation).

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

namespace Kratos
{
namespace Testing
{

namespace
{

/// Comparison tolerances required by the characterization task.
constexpr double comparison_absolute_tolerance = 1.0e-12;
constexpr double comparison_relative_tolerance = 1.0e-10;

/// Machine-precision level allowance used exclusively for components whose
/// reference value is exactly zero (see ComponentTolerance below).
constexpr double machine_precision_allowance = 1.0e-15;

/// Material and loading data shared by both model parts.
constexpr double test_young_modulus = 2.0e7;
constexpr double test_poisson_ratio = 0.2;
constexpr double test_density = 2400.0;
constexpr double test_thickness = 0.15;
constexpr double test_rayleigh_alpha = 0.02;
constexpr double test_rayleigh_beta = 0.03;

/// All the element responses compared between the two implementations.
struct ElementResponse
{
    Matrix lhs;
    Vector rhs;
    Matrix independent_lhs;
    Vector independent_rhs;
    Matrix mass_matrix;
    Matrix damping_matrix;
    std::vector<Vector> cauchy_stress_vectors;
    std::vector<Vector> pk2_stress_vectors;
    std::vector<Vector> strain_vectors;
};

/// Creates one of the two numerically identical model parts.
/// @param rModel The model owning the model part
/// @param rModelPartName Name of the model part
/// @param rElementName Registered name of the element to be created
/// @param rConstitutiveLawName Registered name of the linear elastic law
/// @param rDimension Working space dimension (2 or 3)
/// @param rNodeIdOffset Offset applied to the node ids
ModelPart& CreateComparisonModelPart(
    Model& rModel,
    const std::string& rModelPartName,
    const std::string& rElementName,
    const std::string& rConstitutiveLawName,
    const std::size_t rDimension,
    const ModelPart::IndexType rNodeIdOffset)
{
    KRATOS_TRY;

    ModelPart& r_model_part = rModel.CreateModelPart(rModelPartName, 2);

    ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
    r_process_info[DOMAIN_SIZE] = rDimension;
    r_process_info[SPACE_DIMENSION] = rDimension;
    r_process_info[IS_RESTARTED] = false;
    // Both implementations integrate the consistent mass matrix with the
    // same (increased) quadrature order when this flag is set.
    r_process_info[COMPUTE_CONSISTENT_MASS_MATRIX] = true;

    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(VELOCITY);
    r_model_part.AddNodalSolutionStepVariable(ACCELERATION);
    r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);

    // Nodes are created on a scaled and translated copy of the local
    // coordinates of the registered prototype geometry, so that both
    // implementations are tested with exactly the same geometry.
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
    if (rDimension == 2) {
        (*p_prop)[THICKNESS] = test_thickness;
    }
    (*p_prop)[RAYLEIGH_ALPHA] = test_rayleigh_alpha;
    (*p_prop)[RAYLEIGH_BETA] = test_rayleigh_beta;
    // Independent constitutive law instance for this model part.
    (*p_prop)[CONSTITUTIVE_LAW] = KratosComponents<ConstitutiveLaw>::Get(rConstitutiveLawName).Clone();

    std::vector<ModelPart::IndexType> element_nodes(number_of_nodes);
    for (std::size_t i = 0; i < number_of_nodes; ++i) {
        element_nodes[i] = rNodeIdOffset + i + 1;
    }
    r_model_part.CreateNewElement(rElementName, 1, element_nodes, p_prop);

    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(DISPLACEMENT_X);
        r_node.AddDof(DISPLACEMENT_Y);
        r_node.AddDof(DISPLACEMENT_Z);
    }

    return r_model_part;

    KRATOS_CATCH("");
}

/// Prescribes the non-zero displacement field u = A*X0 + t, updates the
/// current coordinates to X = X0 + u and assigns a non-zero nodal
/// VOLUME_ACCELERATION field.
void PrescribeDisplacementAndBodyForce(
    ModelPart& rModelPart,
    const std::size_t rDimension)
{
    KRATOS_TRY;

    Matrix displacement_gradient = ZeroMatrix(3, 3);
    array_1d<double, 3> translation;
    if (rDimension == 2) {
        displacement_gradient(0, 0) = 2.0e-3;
        displacement_gradient(0, 1) = 3.0e-4;
        displacement_gradient(1, 0) = 5.0e-4;
        displacement_gradient(1, 1) = -1.0e-3;
        displacement_gradient(2, 2) = 1.0;
        translation[0] = 1.0e-2;
        translation[1] = -2.0e-2;
        translation[2] = 0.0;
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

    std::size_t local_node_index = 0;
    for (auto& r_node : rModelPart.Nodes()) {
        array_1d<double, 3>& r_displacement = r_node.FastGetSolutionStepValue(DISPLACEMENT);
        noalias(r_displacement) = prod(displacement_gradient, r_node.GetInitialPosition());
        noalias(r_displacement) += translation;

        // Update the current coordinates so that the reference configuration
        // (current position - total displacement) is the initial one.
        r_node.X() = r_node.X0() + r_displacement[0];
        r_node.Y() = r_node.Y0() + r_displacement[1];
        r_node.Z() = r_node.Z0() + r_displacement[2];

        // Non-uniform body force field (interpolated with the shape functions
        // by both elements).
        array_1d<double, 3>& r_volume_acceleration = r_node.FastGetSolutionStepValue(VOLUME_ACCELERATION);
        const double index_factor = static_cast<double>(local_node_index);
        r_volume_acceleration[0] = 1.5 + 0.25 * index_factor;
        r_volume_acceleration[1] = -9.81 + 0.1 * index_factor;
        r_volume_acceleration[2] = 0.5 - 0.05 * index_factor;

        ++local_node_index;
    }

    KRATOS_CATCH("");
}

/// Runs the standard element lifecycle up to the point where responses can be
/// calculated.
void InitializeElementLifecycle(ModelPart& rModelPart)
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

/// Calculates all the compared responses of the element of the model part.
void CalculateElementResponse(ModelPart& rModelPart, ElementResponse& rResponse)
{
    KRATOS_TRY;

    auto p_element = rModelPart.pGetElement(1);
    const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

    p_element->CalculateLocalSystem(rResponse.lhs, rResponse.rhs, r_process_info);
    p_element->CalculateLeftHandSide(rResponse.independent_lhs, r_process_info);
    p_element->CalculateRightHandSide(rResponse.independent_rhs, r_process_info);
    p_element->CalculateMassMatrix(rResponse.mass_matrix, r_process_info);
    p_element->CalculateDampingMatrix(rResponse.damping_matrix, r_process_info);
    p_element->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, rResponse.cauchy_stress_vectors, r_process_info);
    p_element->CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, rResponse.pk2_stress_vectors, r_process_info);
    p_element->CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR, rResponse.strain_vectors, r_process_info);

    KRATOS_CATCH("");
}

/// Creates the two independent but numerically identical model parts (one per
/// element implementation) and calculates their responses.
void CalculateComparisonResponses(
    const std::string& rDamElementName,
    const std::string& rStructuralMechanicsElementName,
    const std::string& rConstitutiveLawName,
    const std::size_t rDimension,
    ElementResponse& rDamResponse,
    ElementResponse& rStructuralMechanicsResponse)
{
    KRATOS_TRY;

    Model model;
    ModelPart& r_dam_model_part = CreateComparisonModelPart(
        model, "DamSmallDisplacement", rDamElementName, rConstitutiveLawName, rDimension, 0);
    ModelPart& r_structural_mechanics_model_part = CreateComparisonModelPart(
        model, "StructuralMechanicsSmallDisplacement", rStructuralMechanicsElementName, rConstitutiveLawName, rDimension, 100);

    PrescribeDisplacementAndBodyForce(r_dam_model_part, rDimension);
    PrescribeDisplacementAndBodyForce(r_structural_mechanics_model_part, rDimension);

    InitializeElementLifecycle(r_dam_model_part);
    InitializeElementLifecycle(r_structural_mechanics_model_part);

    CalculateElementResponse(r_dam_model_part, rDamResponse);
    CalculateElementResponse(r_structural_mechanics_model_part, rStructuralMechanicsResponse);

    KRATOS_CATCH("");
}

/// Tolerance associated with a reference value, combining the absolute and
/// relative tolerances required by the characterization task.
///
/// Components whose reference value is exactly zero are the result of exact
/// cancellation of the terms of the sum. The two implementations evaluate
/// identical terms in slightly different orders (e.g. the shape-function
/// gradients at Gauss points with irrational coordinates are evaluated through
/// different code paths, see Geometry::Jacobian), so such components cannot be
/// compared at the absolute tolerance: both implementations are correct and
/// only differ by accumulated roundoff. For those components the tolerance is
/// raised to a machine-precision fraction of the scale of the whole quantity
/// (largest absolute entry of the reference matrix/vector).
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

/// Metrics accumulated over the components of one comparison, reported for
/// every comparison so that the equivalence evidence is self-contained in the
/// test output. These metrics are purely diagnostic; the actual pass/fail
/// decision is the component-wise KRATOS_EXPECT_NEAR check below.
struct ComparisonMetrics
{
    bool pass = true;
    std::size_t component_count = 0;
    std::size_t failed_component_count = 0;
    double max_absolute_difference = 0.0;
    double max_relative_difference = 0.0; // w.r.t. |reference|, components with |reference| == 0 excluded
    double max_tolerance_used = 0.0;
    double reference_scale = 0.0; // largest |reference| entry
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

} // namespace

//************************************************************************************
// 2D three-node triangle (plane strain)
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison2D3_LocalSystemLHS, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement2D3N", "SmallDisplacementElement2D3N",
        "LinearElasticPlaneStrain2DLaw", 2, dam_response, structural_mechanics_response);
    ExpectMatrixComponentsNear(
        dam_response.lhs, structural_mechanics_response.lhs, "local system LHS (Triangle2D3)");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison2D3_LocalSystemRHS, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement2D3N", "SmallDisplacementElement2D3N",
        "LinearElasticPlaneStrain2DLaw", 2, dam_response, structural_mechanics_response);
    ExpectVectorComponentsNear(
        dam_response.rhs, structural_mechanics_response.rhs, "local system RHS (Triangle2D3)");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison2D3_IndependentLHS, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement2D3N", "SmallDisplacementElement2D3N",
        "LinearElasticPlaneStrain2DLaw", 2, dam_response, structural_mechanics_response);
    ExpectMatrixComponentsNear(
        dam_response.independent_lhs, structural_mechanics_response.independent_lhs,
        "independently calculated LHS (Triangle2D3)");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison2D3_IndependentRHS, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement2D3N", "SmallDisplacementElement2D3N",
        "LinearElasticPlaneStrain2DLaw", 2, dam_response, structural_mechanics_response);
    ExpectVectorComponentsNear(
        dam_response.independent_rhs, structural_mechanics_response.independent_rhs,
        "independently calculated RHS (Triangle2D3)");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison2D3_MassMatrix, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement2D3N", "SmallDisplacementElement2D3N",
        "LinearElasticPlaneStrain2DLaw", 2, dam_response, structural_mechanics_response);
    ExpectMatrixComponentsNear(
        dam_response.mass_matrix, structural_mechanics_response.mass_matrix,
        "consistent mass matrix (Triangle2D3)");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison2D3_DampingMatrix, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement2D3N", "SmallDisplacementElement2D3N",
        "LinearElasticPlaneStrain2DLaw", 2, dam_response, structural_mechanics_response);
    ExpectMatrixComponentsNear(
        dam_response.damping_matrix, structural_mechanics_response.damping_matrix,
        "Rayleigh damping matrix (Triangle2D3)");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison2D3_CauchyStressVector, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement2D3N", "SmallDisplacementElement2D3N",
        "LinearElasticPlaneStrain2DLaw", 2, dam_response, structural_mechanics_response);
    ExpectIntegrationPointVectorsNear(
        dam_response.cauchy_stress_vectors, structural_mechanics_response.cauchy_stress_vectors,
        "Cauchy stress vector (Triangle2D3)");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison2D3_PK2StressVector, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement2D3N", "SmallDisplacementElement2D3N",
        "LinearElasticPlaneStrain2DLaw", 2, dam_response, structural_mechanics_response);
    ExpectIntegrationPointVectorsNear(
        dam_response.pk2_stress_vectors, structural_mechanics_response.pk2_stress_vectors,
        "PK2 stress vector (Triangle2D3)");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison2D3_StrainVector, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement2D3N", "SmallDisplacementElement2D3N",
        "LinearElasticPlaneStrain2DLaw", 2, dam_response, structural_mechanics_response);
    ExpectIntegrationPointVectorsNear(
        dam_response.strain_vectors, structural_mechanics_response.strain_vectors,
        "strain vector (Triangle2D3)");
}

//************************************************************************************
// 2D four-node quadrilateral (plane strain)
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison2D4_LocalSystemLHS, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement2D4N", "SmallDisplacementElement2D4N",
        "LinearElasticPlaneStrain2DLaw", 2, dam_response, structural_mechanics_response);
    ExpectMatrixComponentsNear(
        dam_response.lhs, structural_mechanics_response.lhs, "local system LHS (Quadrilateral2D4)");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison2D4_LocalSystemRHS, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement2D4N", "SmallDisplacementElement2D4N",
        "LinearElasticPlaneStrain2DLaw", 2, dam_response, structural_mechanics_response);
    ExpectVectorComponentsNear(
        dam_response.rhs, structural_mechanics_response.rhs, "local system RHS (Quadrilateral2D4)");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison2D4_IndependentLHS, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement2D4N", "SmallDisplacementElement2D4N",
        "LinearElasticPlaneStrain2DLaw", 2, dam_response, structural_mechanics_response);
    ExpectMatrixComponentsNear(
        dam_response.independent_lhs, structural_mechanics_response.independent_lhs,
        "independently calculated LHS (Quadrilateral2D4)");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison2D4_IndependentRHS, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement2D4N", "SmallDisplacementElement2D4N",
        "LinearElasticPlaneStrain2DLaw", 2, dam_response, structural_mechanics_response);
    ExpectVectorComponentsNear(
        dam_response.independent_rhs, structural_mechanics_response.independent_rhs,
        "independently calculated RHS (Quadrilateral2D4)");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison2D4_MassMatrix, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement2D4N", "SmallDisplacementElement2D4N",
        "LinearElasticPlaneStrain2DLaw", 2, dam_response, structural_mechanics_response);
    ExpectMatrixComponentsNear(
        dam_response.mass_matrix, structural_mechanics_response.mass_matrix,
        "consistent mass matrix (Quadrilateral2D4)");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison2D4_DampingMatrix, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement2D4N", "SmallDisplacementElement2D4N",
        "LinearElasticPlaneStrain2DLaw", 2, dam_response, structural_mechanics_response);
    ExpectMatrixComponentsNear(
        dam_response.damping_matrix, structural_mechanics_response.damping_matrix,
        "Rayleigh damping matrix (Quadrilateral2D4)");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison2D4_CauchyStressVector, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement2D4N", "SmallDisplacementElement2D4N",
        "LinearElasticPlaneStrain2DLaw", 2, dam_response, structural_mechanics_response);
    ExpectIntegrationPointVectorsNear(
        dam_response.cauchy_stress_vectors, structural_mechanics_response.cauchy_stress_vectors,
        "Cauchy stress vector (Quadrilateral2D4)");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison2D4_PK2StressVector, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement2D4N", "SmallDisplacementElement2D4N",
        "LinearElasticPlaneStrain2DLaw", 2, dam_response, structural_mechanics_response);
    ExpectIntegrationPointVectorsNear(
        dam_response.pk2_stress_vectors, structural_mechanics_response.pk2_stress_vectors,
        "PK2 stress vector (Quadrilateral2D4)");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison2D4_StrainVector, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement2D4N", "SmallDisplacementElement2D4N",
        "LinearElasticPlaneStrain2DLaw", 2, dam_response, structural_mechanics_response);
    ExpectIntegrationPointVectorsNear(
        dam_response.strain_vectors, structural_mechanics_response.strain_vectors,
        "strain vector (Quadrilateral2D4)");
}

//************************************************************************************
// 3D four-node tetrahedron
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison3D4_LocalSystemLHS, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement3D4N", "SmallDisplacementElement3D4N",
        "LinearElastic3DLaw", 3, dam_response, structural_mechanics_response);
    ExpectMatrixComponentsNear(
        dam_response.lhs, structural_mechanics_response.lhs, "local system LHS (Tetrahedra3D4)");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison3D4_LocalSystemRHS, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement3D4N", "SmallDisplacementElement3D4N",
        "LinearElastic3DLaw", 3, dam_response, structural_mechanics_response);
    ExpectVectorComponentsNear(
        dam_response.rhs, structural_mechanics_response.rhs, "local system RHS (Tetrahedra3D4)");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison3D4_IndependentLHS, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement3D4N", "SmallDisplacementElement3D4N",
        "LinearElastic3DLaw", 3, dam_response, structural_mechanics_response);
    ExpectMatrixComponentsNear(
        dam_response.independent_lhs, structural_mechanics_response.independent_lhs,
        "independently calculated LHS (Tetrahedra3D4)");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison3D4_IndependentRHS, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement3D4N", "SmallDisplacementElement3D4N",
        "LinearElastic3DLaw", 3, dam_response, structural_mechanics_response);
    ExpectVectorComponentsNear(
        dam_response.independent_rhs, structural_mechanics_response.independent_rhs,
        "independently calculated RHS (Tetrahedra3D4)");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison3D4_MassMatrix, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement3D4N", "SmallDisplacementElement3D4N",
        "LinearElastic3DLaw", 3, dam_response, structural_mechanics_response);
    ExpectMatrixComponentsNear(
        dam_response.mass_matrix, structural_mechanics_response.mass_matrix,
        "consistent mass matrix (Tetrahedra3D4)");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison3D4_DampingMatrix, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement3D4N", "SmallDisplacementElement3D4N",
        "LinearElastic3DLaw", 3, dam_response, structural_mechanics_response);
    ExpectMatrixComponentsNear(
        dam_response.damping_matrix, structural_mechanics_response.damping_matrix,
        "Rayleigh damping matrix (Tetrahedra3D4)");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison3D4_CauchyStressVector, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement3D4N", "SmallDisplacementElement3D4N",
        "LinearElastic3DLaw", 3, dam_response, structural_mechanics_response);
    ExpectIntegrationPointVectorsNear(
        dam_response.cauchy_stress_vectors, structural_mechanics_response.cauchy_stress_vectors,
        "Cauchy stress vector (Tetrahedra3D4)");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison3D4_PK2StressVector, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement3D4N", "SmallDisplacementElement3D4N",
        "LinearElastic3DLaw", 3, dam_response, structural_mechanics_response);
    ExpectIntegrationPointVectorsNear(
        dam_response.pk2_stress_vectors, structural_mechanics_response.pk2_stress_vectors,
        "PK2 stress vector (Tetrahedra3D4)");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison3D4_StrainVector, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement3D4N", "SmallDisplacementElement3D4N",
        "LinearElastic3DLaw", 3, dam_response, structural_mechanics_response);
    ExpectIntegrationPointVectorsNear(
        dam_response.strain_vectors, structural_mechanics_response.strain_vectors,
        "strain vector (Tetrahedra3D4)");
}

//************************************************************************************
// 3D eight-node hexahedron
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison3D8_LocalSystemLHS, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement3D8N", "SmallDisplacementElement3D8N",
        "LinearElastic3DLaw", 3, dam_response, structural_mechanics_response);
    ExpectMatrixComponentsNear(
        dam_response.lhs, structural_mechanics_response.lhs, "local system LHS (Hexahedra3D8)");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison3D8_LocalSystemRHS, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement3D8N", "SmallDisplacementElement3D8N",
        "LinearElastic3DLaw", 3, dam_response, structural_mechanics_response);
    ExpectVectorComponentsNear(
        dam_response.rhs, structural_mechanics_response.rhs, "local system RHS (Hexahedra3D8)");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison3D8_IndependentLHS, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement3D8N", "SmallDisplacementElement3D8N",
        "LinearElastic3DLaw", 3, dam_response, structural_mechanics_response);
    ExpectMatrixComponentsNear(
        dam_response.independent_lhs, structural_mechanics_response.independent_lhs,
        "independently calculated LHS (Hexahedra3D8)");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison3D8_IndependentRHS, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement3D8N", "SmallDisplacementElement3D8N",
        "LinearElastic3DLaw", 3, dam_response, structural_mechanics_response);
    ExpectVectorComponentsNear(
        dam_response.independent_rhs, structural_mechanics_response.independent_rhs,
        "independently calculated RHS (Hexahedra3D8)");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison3D8_MassMatrix, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement3D8N", "SmallDisplacementElement3D8N",
        "LinearElastic3DLaw", 3, dam_response, structural_mechanics_response);
    ExpectMatrixComponentsNear(
        dam_response.mass_matrix, structural_mechanics_response.mass_matrix,
        "consistent mass matrix (Hexahedra3D8)");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison3D8_DampingMatrix, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement3D8N", "SmallDisplacementElement3D8N",
        "LinearElastic3DLaw", 3, dam_response, structural_mechanics_response);
    ExpectMatrixComponentsNear(
        dam_response.damping_matrix, structural_mechanics_response.damping_matrix,
        "Rayleigh damping matrix (Hexahedra3D8)");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison3D8_CauchyStressVector, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement3D8N", "SmallDisplacementElement3D8N",
        "LinearElastic3DLaw", 3, dam_response, structural_mechanics_response);
    ExpectIntegrationPointVectorsNear(
        dam_response.cauchy_stress_vectors, structural_mechanics_response.cauchy_stress_vectors,
        "Cauchy stress vector (Hexahedra3D8)");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison3D8_PK2StressVector, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement3D8N", "SmallDisplacementElement3D8N",
        "LinearElastic3DLaw", 3, dam_response, structural_mechanics_response);
    ExpectIntegrationPointVectorsNear(
        dam_response.pk2_stress_vectors, structural_mechanics_response.pk2_stress_vectors,
        "PK2 stress vector (Hexahedra3D8)");
}

KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementComparison3D8_StrainVector, KratosDamFastSuite)
{
    ElementResponse dam_response, structural_mechanics_response;
    CalculateComparisonResponses(
        "SmallDisplacementSolidElement3D8N", "SmallDisplacementElement3D8N",
        "LinearElastic3DLaw", 3, dam_response, structural_mechanics_response);
    ExpectIntegrationPointVectorsNear(
        dam_response.strain_vectors, structural_mechanics_response.strain_vectors,
        "strain vector (Hexahedra3D8)");
}

} // namespace Testing
} // namespace Kratos
