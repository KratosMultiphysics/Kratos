//     |  /           |
//     ' /   __| _` | __|  _ \   __|
//     . \  |   (   | |   (   |\__ `
//    _|\_\_|  \__,_|\__|\___/ ____/
//                    Multi-Physics
//
//   Fluid Dynamics Application
//
//   License:         BSD License
//                    Kratos default license: kratos/license.txt
//
//   Main author:     Nicolò Antonelli
//

// System includes
#include <array>

// Project includes
#include "containers/model.h"
#include "integration/line_gauss_legendre_integration_points.h"
#include "utilities/gap_quadrature_utility.h"

// Application includes
#include "custom_constitutive/newtonian_2d_law.h"
#include "custom_utilities/qsvms_gap_element_utility.h"
#include "fluid_dynamics_application_variables.h"
#include "tests/cpp_tests/fluid_dynamics_fast_suite.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(TimeIntegratedQSVMSQ1GapQuadraturePartition, FluidDynamicsApplicationFastSuite)
{
    Model model;
    auto& r_part = model.CreateModelPart("Fluid", 3);
    for (const auto* p_variable : {&VELOCITY, &MESH_VELOCITY, &ACCELERATION, &BODY_FORCE, &ADVPROJ}) {
        r_part.AddNodalSolutionStepVariable(*p_variable);
    }
    r_part.AddNodalSolutionStepVariable(PRESSURE);
    r_part.AddNodalSolutionStepVariable(DIVPROJ);
    r_part.GetProcessInfo()[DELTA_TIME] = 0.1;
    r_part.GetProcessInfo()[DYNAMIC_TAU] = 1.0;
    r_part.GetProcessInfo()[OSS_SWITCH] = 0;
    Vector bdf(3);
    bdf[0] = 15.0; bdf[1] = -20.0; bdf[2] = 5.0;
    r_part.GetProcessInfo()[BDF_COEFFICIENTS] = bdf;

    auto p_material = r_part.CreateNewProperties(1);
    p_material->SetValue(DENSITY, 1.0);
    p_material->SetValue(DYNAMIC_VISCOSITY, 0.1);
    p_material->SetValue(CONSTITUTIVE_LAW, Kratos::make_shared<Newtonian2DLaw>());
    r_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_part.CreateNewNode(3, 1.0, 1.0, 0.0);
    r_part.CreateNewNode(4, 0.0, 1.0, 0.0);
    for (auto& r_node : r_part.Nodes()) {
        r_node.AddDof(VELOCITY_X);
        r_node.AddDof(VELOCITY_Y);
        r_node.AddDof(PRESSURE);
        for (std::size_t step = 0; step < 3; ++step) {
            // Constant nonzero relative convection keeps tau constant while
            // exercising the convective QSVMS operator on custom quadrature.
            // Different history values give a nonzero BDF derivative. The
            // gap must use these same fixed-background nodal histories.
            r_node.FastGetSolutionStepValue(VELOCITY, step)[0] = 0.4 - 0.03*step;
            r_node.FastGetSolutionStepValue(MESH_VELOCITY, step)[0] = 0.2;
        }
    }
    auto p_parent = r_part.CreateNewElement("TimeIntegratedQSVMS2D4N", 1, {1, 2, 3, 4}, p_material);
    p_parent->SetValue(C_SMAGORINSKY, 0.13);
    p_parent->Initialize(r_part.GetProcessInfo());
    Matrix expected_lhs;
    Vector expected_rhs;
    p_parent->CalculateLocalSystem(expected_lhs, expected_rhs, r_part.GetProcessInfo());

    // The two integration triangles partition the parent. The custom
    // quadrature geometries carry the same parent nodes and Q1 gradients.
    const std::array<std::array<std::array<double, 2>, 3>, 2> triangles{{
        {{{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}}},
        {{{0.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}}}}};
    Matrix total_lhs = ZeroMatrix(12, 12);
    Vector total_rhs = ZeroVector(12);
    for (const auto& r_triangle : triangles) {
        Geometry<Node>::IntegrationPointsArrayType points;
        for (const auto& r_a : LineGaussLegendreIntegrationPoints3::IntegrationPoints()) {
            for (const auto& r_b : LineGaussLegendreIntegrationPoints3::IntegrationPoints()) {
                const double u = 0.5 * (1.0 + r_a.X());
                const double v = 0.5 * (1.0 + r_b.X());
                const double x = r_triangle[0][0] + u*(r_triangle[1][0]-r_triangle[0][0])
                    + (1.0-u)*v*(r_triangle[2][0]-r_triangle[0][0]);
                const double y = r_triangle[0][1] + u*(r_triangle[1][1]-r_triangle[0][1])
                    + (1.0-u)*v*(r_triangle[2][1]-r_triangle[0][1]);
                points.emplace_back(x, y, 0.0, 0.25*(1.0-u)*r_a.Weight()*r_b.Weight());
            }
        }
        auto p_geometry = GapQuadratureUtility::Create(p_parent->GetGeometry(), points);
        auto p_gap = QSVMSGapElementUtility::Create(2, *p_parent, p_geometry);
        KRATOS_EXPECT_NEAR(p_gap->GetValue(C_SMAGORINSKY), 0.13, 1e-12);
        p_gap->Initialize(r_part.GetProcessInfo());
        Matrix lhs;
        Vector rhs;
        p_gap->CalculateLocalSystem(lhs, rhs, r_part.GetProcessInfo());
        noalias(total_lhs) += lhs;
        noalias(total_rhs) += rhs;
    }
    KRATOS_EXPECT_MATRIX_NEAR(total_lhs, expected_lhs, 1e-10);
    KRATOS_EXPECT_VECTOR_NEAR(total_rhs, expected_rhs, 1e-10);
}

} // namespace Kratos::Testing
