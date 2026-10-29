//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes
#include <cmath>
#include <utility>
#include <vector>

// External includes

// Project includes
#include "containers/model.h"
#include "fluid_dynamics_application_variables.h"
#include "tests/cpp_tests/fluid_dynamics_fast_suite.h"

namespace Kratos::Testing
{
namespace
{

template<unsigned int TDim>
void CheckFractionalStepDiscontinuousEmbeddedVelocity(const int FractionalStep)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("Main", 2);
    r_model_part.AddNodalSolutionStepVariable(VELOCITY);
    r_model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    r_model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    r_model_part.AddNodalSolutionStepVariable(PRESSURE);
    r_model_part.AddNodalSolutionStepVariable(DENSITY);
    r_model_part.AddNodalSolutionStepVariable(VISCOSITY);
    r_model_part.AddNodalSolutionStepVariable(CONV_PROJ);
    r_model_part.AddNodalSolutionStepVariable(PRESS_PROJ);
    r_model_part.AddNodalSolutionStepVariable(DIVPROJ);

    auto& r_process_info = r_model_part.GetProcessInfo();
    r_process_info[FRACTIONAL_STEP] = FractionalStep;
    r_process_info[DELTA_TIME] = 0.1;
    r_process_info[DYNAMIC_TAU] = 1.0;
    r_process_info[FS_PRESSURE_GRADIENT_RELAXATION_FACTOR] = 1.0;
    Vector bdf_coefficients(2);
    bdf_coefficients[0] = 10.0;
    bdf_coefficients[1] = -10.0;
    r_process_info[BDF_COEFFICIENTS] = bdf_coefficients;

    // Put the origin last: the two cut edges then have matching indices in
    // the triangle splitting utility and the element's edge traversal.
    std::vector<ModelPart::IndexType> node_ids;
    for (unsigned int d = 0; d < TDim; ++d) {
        array_1d<double, 3> coordinates = ZeroVector(3);
        coordinates[d] = 1.0;
        r_model_part.CreateNewNode(d + 1, coordinates[0], coordinates[1], coordinates[2]);
        node_ids.push_back(d + 1);
    }
    if constexpr (TDim == 3) {
        std::swap(node_ids[0], node_ids[1]);
    }
    r_model_part.CreateNewNode(TDim + 1, 0.0, 0.0, 0.0);
    node_ids.push_back(TDim + 1);
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.FastGetSolutionStepValue(DENSITY) = 1.0;
        r_node.FastGetSolutionStepValue(VISCOSITY) = 1.0;
    }
    const auto p_properties = r_model_part.CreateNewProperties(0);
    const auto p_element = r_model_part.CreateNewElement(
        TDim == 2 ? "FractionalStepDiscontinuous2D" : "FractionalStepDiscontinuous3D",
        1, node_ids, p_properties);
    p_element->SetValue(SPLIT_ELEMENT, true);
    p_element->Set(SLIP, false);
    p_element->SetValue(C_SMAGORINSKY, 0.0);
    Vector distances(TDim + 1, 0.5);
    distances[TDim] = -0.5;
    p_element->SetValue(ELEMENTAL_DISTANCES, distances);

    array_1d<double, 3> embedded_velocity = ZeroVector(3);
    for (unsigned int d = 0; d < TDim; ++d) {
        embedded_velocity[d] = d + 1.0;
    }
    p_element->SetValue(EMBEDDED_VELOCITY, embedded_velocity);
    Matrix lhs;
    Vector rhs;
    p_element->CalculateLocalSystem(lhs, rhs, r_process_info);

    const unsigned int local_size = FractionalStep == 5 ? TDim + 1 : TDim * (TDim + 1);
    KRATOS_EXPECT_EQ(lhs.size1(), local_size);
    KRATOS_EXPECT_EQ(lhs.size2(), local_size);
    KRATOS_EXPECT_EQ(rhs.size(), local_size);
    for (unsigned int i = 0; i < local_size; ++i) {
        KRATOS_EXPECT_TRUE(std::isfinite(rhs[i]));
        for (unsigned int j = 0; j < local_size; ++j) {
            KRATOS_EXPECT_TRUE(std::isfinite(lhs(i, j)));
        }
    }

    // The cut joins the midpoints of the edges from the origin. Its normal
    // is (1,...,1)/sqrt(TDim), and each cut edge receives 1/TDim of its measure.
    const double interface_measure = TDim == 2 ? std::sqrt(2.0) / 2.0 : std::sqrt(3.0) / 8.0;
    const double edge_measure = interface_measure / TDim;
    const double mean_velocity = (TDim + 1.0) / 2.0;
    const double normal_velocity = mean_velocity * std::sqrt(static_cast<double>(TDim));
    const double penalty = 1.0 + norm_2(embedded_velocity); // Density, viscosity and minimum edge length are one.
    Vector expected_rhs = ZeroVector(local_size);
    for (unsigned int i = 0; i < TDim + 1; ++i) {
        const double nodal_measure = i == TDim ? interface_measure : edge_measure;
        if (FractionalStep == 5) {
            expected_rhs[i] = (i == TDim ? -1.0 : 1.0) * nodal_measure * normal_velocity;
        } else {
            for (unsigned int d = 0; d < TDim; ++d) {
                const double tangential_velocity = embedded_velocity[d] - mean_velocity;
                expected_rhs[i * TDim + d] = nodal_measure *
                    (2.0 * tangential_velocity + penalty * mean_velocity);
            }
        }
    }
    KRATOS_EXPECT_VECTOR_NEAR(rhs, expected_rhs, 1e-10);

    // Reversing the wall velocity preserves the penalty and LHS, and reverses
    // the boundary forcing because the fluid and all volume sources are zero.
    embedded_velocity *= -1.0;
    p_element->SetValue(EMBEDDED_VELOCITY, embedded_velocity);
    Matrix reversed_lhs;
    Vector reversed_rhs;
    p_element->CalculateLocalSystem(reversed_lhs, reversed_rhs, r_process_info);
    expected_rhs *= -1.0;
    KRATOS_EXPECT_MATRIX_NEAR(reversed_lhs, lhs, 1e-10);
    KRATOS_EXPECT_VECTOR_NEAR(reversed_rhs, expected_rhs, 1e-10);
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(FractionalStepDiscontinuous2DPressure, FluidDynamicsApplicationFastSuite)
{
    CheckFractionalStepDiscontinuousEmbeddedVelocity<2>(5);
}

KRATOS_TEST_CASE_IN_SUITE(FractionalStepDiscontinuous3DPressure, FluidDynamicsApplicationFastSuite)
{
    CheckFractionalStepDiscontinuousEmbeddedVelocity<3>(5);
}

KRATOS_TEST_CASE_IN_SUITE(FractionalStepDiscontinuous2DMomentum, FluidDynamicsApplicationFastSuite)
{
    CheckFractionalStepDiscontinuousEmbeddedVelocity<2>(1);
}

KRATOS_TEST_CASE_IN_SUITE(FractionalStepDiscontinuous3DMomentum, FluidDynamicsApplicationFastSuite)
{
    CheckFractionalStepDiscontinuousEmbeddedVelocity<3>(1);
}

} // namespace Kratos::Testing
