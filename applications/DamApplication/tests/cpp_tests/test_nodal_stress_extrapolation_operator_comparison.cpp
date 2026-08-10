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

// Phase 3C (focused) characterization: element-level extrapolation operator of
// the DamApplication legacy nodal stress recovery vs the GeoMechanics
// LinearNodalExtrapolator / ExtrapolationUtilities.
//
// Legacy Dam element-level operator (Quadrilateral2D4 / Hexahedra3D8):
//   nodal_element = E_Dam * integration_point_values
// where E_Dam is the hard-coded bilinear/trilinear nodal extrapolation matrix
// from PoroElementUtilities::Calculate2D/3DExtrapolationMatrix. For the
// single-Gauss-point geometries (Triangle2D3 / Tetrahedra3D4) the legacy assigns
// the single GP value to every node (E = ones).
//
// GeoMechanics element-level operator (LinearNodalExtrapolator):
//   nodal_element = E_Geo * integration_point_values
// where E_Geo is built as a least-squares (L2-projection) matrix
//   M = sum_gp N(gp) N(gp)^T * W_gp * J_gp ;  B(:,gp) = N(gp) * W_gp * J_gp
//   E_Geo = M^{-1} * B
// and for a single integration point E_Geo = ones.
//
// Global nodal averaging:
//   Legacy Dam: sum(element_measure * nodal_element) / sum(element_measure)
//   GeoMechanics: sum(nodal_element) / number_of_connected_elements
//
// These tests replicate the GeoMechanics least-squares matrix with the formula
// read from LinearNodalExtrapolator::CalculateExtrapolationMatrixForCornerNodes
// (to avoid introducing a DamApplication dependency on GeoMechanicsApplication),
// and compare the element-level operators and the global averaging policies.

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
#include "custom_utilities/poro_element_utilities.hpp"

namespace Kratos
{
namespace Testing
{

namespace
{

/// Material data (shared).
constexpr double test_young_modulus = 2.0e7;
constexpr double test_poisson_ratio = 0.2;
constexpr double test_density = 2400.0;
constexpr double test_thermal_expansion = 1.0e-5;
constexpr double test_reference_temperature = 20.0;
constexpr double test_thickness = 0.15;

/// Creates a model part with a single element of the given type and the thermal
/// law, with a prescribed varying thermo-mechanical state.
ModelPart& CreateOperatorModelPart(
    Model& rModel,
    const std::string& rModelPartName,
    const std::string& rElementName,
    const std::string& rLawName,
    const std::size_t rDimension)
{
    KRATOS_TRY;

    ModelPart& r_model_part = rModel.CreateModelPart(rModelPartName, 2);
    ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
    r_process_info[DOMAIN_SIZE] = rDimension;
    r_process_info[SPACE_DIMENSION] = rDimension;
    r_process_info[IS_RESTARTED] = false;

    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    r_model_part.AddNodalSolutionStepVariable(NODAL_REFERENCE_TEMPERATURE);
    r_model_part.AddNodalSolutionStepVariable(CAUCHY_STRESS_TENSOR);
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
        r_model_part.CreateNewNode(i + 1, x, y, z);
    }

    auto p_prop = r_model_part.CreateNewProperties(1);
    (*p_prop)[YOUNG_MODULUS] = test_young_modulus;
    (*p_prop)[POISSON_RATIO] = test_poisson_ratio;
    (*p_prop)[DENSITY] = test_density;
    (*p_prop)[THERMAL_EXPANSION] = test_thermal_expansion;
    if (rDimension == 2) {
        (*p_prop)[THICKNESS] = test_thickness;
    }
    if (rLawName == "ThermalLinearElastic2DPlaneStrain") {
        p_prop->SetValue(CONSTITUTIVE_LAW, ThermalLinearElastic2DPlaneStrain().Clone());
    } else if (rLawName == "ThermalLinearElastic2DPlaneStress") {
        p_prop->SetValue(CONSTITUTIVE_LAW, ThermalLinearElastic2DPlaneStress().Clone());
    } else {
        p_prop->SetValue(CONSTITUTIVE_LAW, ThermalLinearElastic3DLaw().Clone());
    }

    std::vector<ModelPart::IndexType> element_nodes(number_of_nodes);
    for (std::size_t i = 0; i < number_of_nodes; ++i) {
        element_nodes[i] = i + 1;
    }
    r_model_part.CreateNewElement(rElementName, 1, element_nodes, p_prop);

    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(DISPLACEMENT_X);
        r_node.AddDof(DISPLACEMENT_Y);
        r_node.AddDof(DISPLACEMENT_Z);
        r_node.FastGetSolutionStepValue(NODAL_REFERENCE_TEMPERATURE) = test_reference_temperature;
        r_node.FastGetSolutionStepValue(NODAL_AREA) = 0.0;
        Matrix zero_tensor(3, 3);
        noalias(zero_tensor) = ZeroMatrix(3, 3);
        r_node.FastGetSolutionStepValue(CAUCHY_STRESS_TENSOR) = zero_tensor;
        r_node.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR) = zero_tensor;
    }

    // Prescribe a state with a non-uniform integration-point stress field.
    std::size_t node_index = 0;
    for (auto& r_node : r_model_part.Nodes()) {
        const array_1d<double, 3>& r_x0 = r_node.GetInitialPosition();
        array_1d<double, 3>& r_disp = r_node.FastGetSolutionStepValue(DISPLACEMENT);
        r_disp[0] = 1.0e-3 * r_x0[0] + 3.0e-4 * r_x0[1] + 1.0e-4 * r_x0[0] * r_x0[1];
        r_disp[1] = 3.0e-4 * r_x0[0] - 1.0e-3 * r_x0[1] + 2.0e-4 * r_x0[0] * r_x0[0];
        r_disp[2] = 5.0e-4 * r_x0[2];
        r_node.X() = r_x0[0] + r_disp[0];
        r_node.Y() = r_x0[1] + r_disp[1];
        r_node.Z() = r_x0[2] + r_disp[2];
        r_node.FastGetSolutionStepValue(TEMPERATURE) =
            test_reference_temperature + 5.0 * static_cast<double>(node_index);
        ++node_index;
    }

    auto p_element = r_model_part.pGetElement(1);
    KRATOS_EXPECT_EQ(p_element->Check(r_process_info), 0);
    p_element->Initialize(r_process_info);
    p_element->InitializeSolutionStep(r_process_info);
    p_element->InitializeNonLinearIteration(r_process_info);

    return r_model_part;

    KRATOS_CATCH("");
}

/// Replicates the GeoMechanics least-squares (L2 projection) extrapolation matrix
/// E_Geo = M^{-1} * B, following
/// LinearNodalExtrapolator::CalculateExtrapolationMatrixForCornerNodes:
///   M = sum_gp N(gp) N(gp)^T * J_gp * W_gp ; B(:,gp) = N(gp) * J_gp * W_gp
/// For a single integration point the operator is the identity column of ones.
Matrix CalculateLeastSquaresExtrapolationMatrix(
    const Geometry<Node>& rGeometry,
    const GeometryData::IntegrationMethod rMethod)
{
    const auto& r_integration_points = rGeometry.IntegrationPoints(rMethod);
    const std::size_t number_of_gps = r_integration_points.size();
    const std::size_t number_of_nodes = rGeometry.size();

    if (number_of_gps == 1) {
        Matrix ones(number_of_nodes, 1);
        noalias(ones) = scalar_matrix<double>(number_of_nodes, 1, 1.0);
        return ones;
    }

    Vector determinants(number_of_gps);
    rGeometry.DeterminantOfJacobian(determinants, rMethod);

    Matrix quasi_mass(number_of_nodes, number_of_nodes);
    noalias(quasi_mass) = ZeroMatrix(number_of_nodes, number_of_nodes);
    Matrix node_coefficients(number_of_nodes, number_of_gps);
    noalias(node_coefficients) = ZeroMatrix(number_of_nodes, number_of_gps);

    for (std::size_t gp = 0; gp < number_of_gps; ++gp) {
        Vector N(number_of_nodes);
        rGeometry.ShapeFunctionsValues(N, r_integration_points[gp].Coordinates());
        const double coefficient = determinants[gp] * r_integration_points[gp].Weight();
        noalias(quasi_mass) += outer_prod(N, N) * coefficient;
        column(node_coefficients, gp) = N * coefficient;
    }

    Matrix quasi_mass_inverse;
    double determinant;
    MathUtils<double>::InvertMatrix(quasi_mass, quasi_mass_inverse, determinant);
    return prod(quasi_mass_inverse, node_coefficients);
}

/// Reports the maximum absolute and relative differences between two vectors.
void ReportVectorDifference(
    const Vector& rComputed,
    const Vector& rReference,
    const std::string& rWhat)
{
    double max_abs = 0.0;
    double max_rel = 0.0;
    for (std::size_t i = 0; i < rReference.size(); ++i) {
        const double abs_diff = std::abs(rComputed(i) - rReference(i));
        max_abs = std::max(max_abs, abs_diff);
        if (std::abs(rReference(i)) > 0.0)
            max_rel = std::max(max_rel, abs_diff / std::abs(rReference(i)));
    }
    std::cout << "[CHARACTERIZATION] " << rWhat << ": max_abs_diff=" << max_abs
              << " max_rel_diff=" << max_rel << std::endl;
}

} // namespace

//************************************************************************************
// Element-level operator comparison (E_Dam vs E_Geo) on a non-uniform GP field
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(NodalExtrapolationOperator_ElementLevel, KratosDamFastSuite)
{
    // Configurations: element name, law, dimension, label.
    struct Config
    {
        const char* element_name;
        const char* law_name;
        std::size_t dimension;
        const char* label;
    };
    const Config configs[] = {
        {"SmallDisplacementElement2D3N", "ThermalLinearElastic2DPlaneStrain", 2, "Triangle2D3"},
        {"SmallDisplacementElement2D4N", "ThermalLinearElastic2DPlaneStrain", 2, "Quadrilateral2D4"},
        {"SmallDisplacementElement3D4N", "ThermalLinearElastic3DLaw", 3, "Tetrahedra3D4"},
        {"SmallDisplacementElement3D8N", "ThermalLinearElastic3DLaw", 3, "Hexahedra3D8"}};

    for (const Config& r_config : configs) {
        Model model;
        ModelPart& r_model_part = CreateOperatorModelPart(
            model, std::string("Operator") + r_config.label, r_config.element_name,
            r_config.law_name, r_config.dimension);
        auto p_element = r_model_part.pGetElement(1);
        const ProcessInfo& r_pi = r_model_part.GetProcessInfo();

        // Gauss-point Cauchy stress tensor (2x2 in 2D, 3x3 in 3D) per GP.
        std::vector<Matrix> gauss_stress;
        p_element->CalculateOnIntegrationPoints(CAUCHY_STRESS_TENSOR, gauss_stress, r_pi);

        const auto& r_geometry = p_element->GetGeometry();
        const GeometryData::IntegrationMethod method = p_element->GetIntegrationMethod();
        const auto& r_integration_points = r_geometry.IntegrationPoints(method);
        const std::size_t number_of_gps = r_integration_points.size();
        const std::size_t number_of_nodes = r_geometry.size();

        // E_Dam: hard-coded for Q4/H8; for single-GP geometries the legacy
        // assigns the single GP value to every node (E = ones).
        Matrix e_dam;
        if (r_config.dimension == 2) {
            BoundedMatrix<double, 4, 4> q4;
            PoroElementUtilities::Calculate2DExtrapolationMatrix(q4);
            e_dam = Matrix(q4);
        } else {
            BoundedMatrix<double, 8, 8> h8;
            PoroElementUtilities::Calculate3DExtrapolationMatrix(h8);
            e_dam = Matrix(h8);
        }
        if (number_of_gps == 1) {
            e_dam = scalar_matrix<double>(number_of_nodes, 1, 1.0);
        }

        // E_Geo: least-squares matrix (replicated from LinearNodalExtrapolator).
        const Matrix e_geo = CalculateLeastSquaresExtrapolationMatrix(r_geometry, method);

        // Compare the operators on a synthetic non-uniform GP field (the xx
        // component of the actual stress).
        Vector sigma_gp(number_of_gps);
        for (std::size_t gp = 0; gp < number_of_gps; ++gp)
            sigma_gp(gp) = gauss_stress[gp](0, 0);

        const Vector nodal_dam = prod(e_dam, sigma_gp);
        const Vector nodal_geo = prod(e_geo, sigma_gp);

        std::cout << "[CHARACTERIZATION] " << r_config.label
                  << " (nGP=" << number_of_gps << ") gauss_xx = " << sigma_gp << std::endl;
        ReportVectorDifference(nodal_geo, nodal_dam,
                               std::string(r_config.label) + " element-level E_Geo*sigma vs E_Dam*sigma");

        if (number_of_gps == 1) {
            // Both operators coincide for single-GP geometries (nodal = sigma).
            for (std::size_t i = 0; i < number_of_nodes; ++i)
                KRATOS_EXPECT_NEAR(nodal_geo(i), nodal_dam(i), 1.0e-12);
        } else {
            // Report whether the element-level operators coincide; no assertion is
            // imposed here (the result is the characterization finding).
            std::cout << "[CHARACTERIZATION] " << r_config.label
                      << ": element-level operator equivalence = "
                      << (norm_inf(nodal_geo - nodal_dam) < 1.0e-8 ? "coincide" : "differ") << std::endl;
        }
    }
}

//************************************************************************************
// Shared-node global averaging: legacy (area-weighted) vs Geo (count-based)
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(NodalExtrapolationOperator_SharedNodeAveraging, KratosDamFastSuite)
{
    // Two connected quads sharing nodes 2 and 3, areas 1 and 2 (see phase 3C).
    const double coords[6][2] = {
        {0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}, {3.0, 0.0}, {3.0, 1.0}};

    auto build_model = [&coords](Model& r_model, const std::string& r_element_name) -> ModelPart& {
        ModelPart& r_model_part = r_model.CreateModelPart("SharedNodes" + r_element_name, 2);
        ProcessInfo& r_pi = r_model_part.GetProcessInfo();
        r_pi[DOMAIN_SIZE] = 2;
        r_pi[SPACE_DIMENSION] = 2;
        r_pi[IS_RESTARTED] = false;
        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
        r_model_part.AddNodalSolutionStepVariable(NODAL_REFERENCE_TEMPERATURE);
        r_model_part.AddNodalSolutionStepVariable(CAUCHY_STRESS_TENSOR);
        r_model_part.AddNodalSolutionStepVariable(NODAL_CAUCHY_STRESS_TENSOR);
        r_model_part.AddNodalSolutionStepVariable(NODAL_AREA);
        for (std::size_t i = 0; i < 6; ++i)
            r_model_part.CreateNewNode(i + 1, coords[i][0], coords[i][1], 0.0);
        for (auto& r_node : r_model_part.Nodes()) {
            r_node.AddDof(DISPLACEMENT_X);
            r_node.AddDof(DISPLACEMENT_Y);
            r_node.AddDof(DISPLACEMENT_Z);
            r_node.FastGetSolutionStepValue(NODAL_REFERENCE_TEMPERATURE) = test_reference_temperature;
            r_node.FastGetSolutionStepValue(NODAL_AREA) = 0.0;
            Matrix zero_tensor(3, 3);
            noalias(zero_tensor) = ZeroMatrix(3, 3);
            r_node.FastGetSolutionStepValue(CAUCHY_STRESS_TENSOR) = zero_tensor;
            r_node.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR) = zero_tensor;
        }
        auto p_prop = r_model_part.CreateNewProperties(1);
        (*p_prop)[YOUNG_MODULUS] = test_young_modulus;
        (*p_prop)[POISSON_RATIO] = test_poisson_ratio;
        (*p_prop)[DENSITY] = test_density;
        (*p_prop)[THERMAL_EXPANSION] = test_thermal_expansion;
        (*p_prop)[THICKNESS] = test_thickness;
        p_prop->SetValue(CONSTITUTIVE_LAW, ThermalLinearElastic2DPlaneStrain().Clone());
        r_model_part.CreateNewElement(r_element_name, 1, {1, 2, 3, 4}, p_prop);
        r_model_part.CreateNewElement(r_element_name, 2, {2, 5, 6, 3}, p_prop);
        for (auto& r_node : r_model_part.Nodes()) {
            const array_1d<double, 3>& r_x0 = r_node.GetInitialPosition();
            array_1d<double, 3>& r_disp = r_node.FastGetSolutionStepValue(DISPLACEMENT);
            r_disp[0] = 1.0e-3 * r_x0[0] + 3.0e-4 * r_x0[1] + 1.0e-4 * r_x0[0] * r_x0[1];
            r_disp[1] = 3.0e-4 * r_x0[0] - 1.0e-3 * r_x0[1] + 2.0e-4 * r_x0[0] * r_x0[0];
            r_node.X() = r_x0[0] + r_disp[0];
            r_node.Y() = r_x0[1] + r_disp[1];
            r_node.FastGetSolutionStepValue(TEMPERATURE) =
                test_reference_temperature + 3.0 * r_x0[0] + 2.0 * r_x0[1];
        }
        for (auto& r_element : r_model_part.Elements()) {
            KRATOS_EXPECT_EQ(r_element.Check(r_pi), 0);
            r_element.Initialize(r_pi);
            r_element.InitializeSolutionStep(r_pi);
            r_element.InitializeNonLinearIteration(r_pi);
        }
        return r_model_part;
    };

    Model model;
    ModelPart& r_model_part = build_model(model, "SmallDisplacementElement2D4N");

    // Element-level extrapolation with both operators, plus the two global
    // averaging policies, computed analytically for the shared-node mesh.
    // Quad 1: nodes {1,2,3,4} area 1; Quad 2: nodes {2,5,6,3} area 2.
    const std::vector<std::vector<std::size_t>> element_nodes = {{1, 2, 3, 4}, {2, 5, 6, 3}};
    const double element_areas[2] = {1.0, 2.0};

    std::vector<Matrix> nodal_dam_avg(6, ZeroMatrix(2, 2));
    std::vector<Matrix> nodal_geo_avg(6, ZeroMatrix(2, 2));
    std::vector<std::size_t> connected_count(6, 0);

    for (std::size_t e = 0; e < 2; ++e) {
        auto p_element = r_model_part.pGetElement(e + 1);
        const ProcessInfo& r_pi = r_model_part.GetProcessInfo();
        std::vector<Matrix> gauss_stress;
        p_element->CalculateOnIntegrationPoints(CAUCHY_STRESS_TENSOR, gauss_stress, r_pi);

        BoundedMatrix<double, 4, 4> q4;
        PoroElementUtilities::Calculate2DExtrapolationMatrix(q4);
        const Matrix e_dam = Matrix(q4);
        const Matrix e_geo = CalculateLeastSquaresExtrapolationMatrix(
            p_element->GetGeometry(), p_element->GetIntegrationMethod());

        // Per-node element-level extrapolated stress (each GP gives a 2x2 tensor,
        // applied component-wise by the extrapolation matrix).
        const std::size_t n_gp = gauss_stress.size();
        for (std::size_t node_index = 0; node_index < 4; ++node_index) {
            const std::size_t global_node = element_nodes[e][node_index];
            ++connected_count[global_node - 1];
            for (std::size_t comp = 0; comp < 4; ++comp) { // 2x2 components
                const std::size_t i = comp / 2;
                const std::size_t j = comp % 2;
                Vector sigma_gp(n_gp);
                for (std::size_t gp = 0; gp < n_gp; ++gp)
                    sigma_gp(gp) = gauss_stress[gp](i, j);
                const double dam_nodal = inner_prod(row(e_dam, node_index), sigma_gp);
                const double geo_nodal = inner_prod(row(e_geo, node_index), sigma_gp);
                // Legacy: area-weighted accumulation.
                nodal_dam_avg[global_node - 1](i, j) += element_areas[e] * dam_nodal;
                // Geo: count-based accumulation.
                nodal_geo_avg[global_node - 1](i, j) += geo_nodal;
            }
        }
    }
    const double total_area[6] = {1.0, 3.0, 3.0, 1.0, 2.0, 2.0};
    for (std::size_t node = 0; node < 6; ++node) {
        nodal_dam_avg[node] = nodal_dam_avg[node] / total_area[node];
        nodal_geo_avg[node] = nodal_geo_avg[node] / static_cast<double>(connected_count[node]);
        std::cout << "[CHARACTERIZATION] shared node " << (node + 1)
                  << " legacy_area_avg(0,0)=" << nodal_dam_avg[node](0, 0)
                  << " geo_count_avg(0,0)=" << nodal_geo_avg[node](0, 0)
                  << " connected=" << connected_count[node]
                  << " area=" << total_area[node] << std::endl;
        ReportVectorDifference(
            row(nodal_geo_avg[node], 0), row(nodal_dam_avg[node], 0),
            "node " + std::to_string(node + 1) + " final geo(count) vs legacy(area) row0");
    }
    std::cout << "[CHARACTERIZATION] Shared-node conclusion: legacy uses element-area weights, "
              << "GeoMechanics uses count-of-connected-elements weights; the element-level operators "
              << "E_Dam and E_Geo are compared in the element-level test." << std::endl;
}

} // namespace Testing
} // namespace Kratos
