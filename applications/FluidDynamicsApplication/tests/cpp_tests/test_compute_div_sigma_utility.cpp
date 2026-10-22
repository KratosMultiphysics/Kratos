//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Antonelli
//

// Project includes
#include "custom_utilities/compute_div_sigma_utility.h"
#include "tests/cpp_tests/fluid_dynamics_fast_suite.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(ComputeDivSigmaUtility2D, FluidDynamicsApplicationFastSuite)
{
    Matrix sigma_values(2, 3);
    sigma_values(0, 0) = 0.0;
    sigma_values(0, 1) = 9.0;
    sigma_values(0, 2) = 7.0;
    sigma_values(1, 0) = 1.0;
    sigma_values(1, 1) = 23.0;
    sigma_values(1, 2) = 16.0;

    Matrix shape_function_values(2, 2);
    shape_function_values(0, 0) = 1.0;
    shape_function_values(0, 1) = 2.0;
    shape_function_values(1, 0) = 3.0;
    shape_function_values(1, 1) = 5.0;

    Matrix shape_function_values_dx(2, 2);
    shape_function_values_dx(0, 0) = 1.0;
    shape_function_values_dx(0, 1) = 0.5;
    shape_function_values_dx(1, 0) = -1.0;
    shape_function_values_dx(1, 1) = 2.0;

    Matrix shape_function_values_dy(2, 2);
    shape_function_values_dy(0, 0) = 0.0;
    shape_function_values_dy(0, 1) = 2.0;
    shape_function_values_dy(1, 0) = 3.0;
    shape_function_values_dy(1, 1) = -1.0;

    ComputeDivSigmaUtility div_sigma_utility;
    const Matrix divergence = div_sigma_utility.ComputeDivergence(
        sigma_values,
        shape_function_values,
        shape_function_values_dx,
        shape_function_values_dy);

    Matrix expected_divergence(2, 2);
    expected_divergence(0, 0) = 11.5;
    expected_divergence(0, 1) = 7.5;
    expected_divergence(1, 0) = -18.0;
    expected_divergence(1, 1) = 12.0;

    KRATOS_EXPECT_MATRIX_NEAR(divergence, expected_divergence, 1.0e-12);
}

KRATOS_TEST_CASE_IN_SUITE(ComputeDivSigmaUtility3D, FluidDynamicsApplicationFastSuite)
{
    Matrix sigma_values(3, 6);
    sigma_values(0, 0) = 1.0;
    sigma_values(0, 1) = 5.0;
    sigma_values(0, 2) = 4.0;
    sigma_values(0, 3) = 2.0;
    sigma_values(0, 4) = 7.0;
    sigma_values(0, 5) = 0.0;
    sigma_values(1, 0) = 2.0;
    sigma_values(1, 1) = -2.0;
    sigma_values(1, 2) = -3.0;
    sigma_values(1, 3) = 1.0;
    sigma_values(1, 4) = -1.0;
    sigma_values(1, 5) = 4.0;
    sigma_values(2, 0) = 4.0;
    sigma_values(2, 1) = 2.0;
    sigma_values(2, 2) = 3.0;
    sigma_values(2, 3) = 7.0;
    sigma_values(2, 4) = 9.0;
    sigma_values(2, 5) = 2.0;

    Matrix shape_function_values(3, 3);
    shape_function_values(0, 0) = 1.0;
    shape_function_values(0, 1) = 0.0;
    shape_function_values(0, 2) = 2.0;
    shape_function_values(1, 0) = 0.0;
    shape_function_values(1, 1) = 1.0;
    shape_function_values(1, 2) = -1.0;
    shape_function_values(2, 0) = 2.0;
    shape_function_values(2, 1) = 1.0;
    shape_function_values(2, 2) = 1.0;

    Matrix shape_function_values_dx(3, 3);
    shape_function_values_dx(0, 0) = 1.0;
    shape_function_values_dx(0, 1) = 0.0;
    shape_function_values_dx(0, 2) = 2.0;
    shape_function_values_dx(1, 0) = -1.0;
    shape_function_values_dx(1, 1) = 1.0;
    shape_function_values_dx(1, 2) = 0.0;
    shape_function_values_dx(2, 0) = 0.5;
    shape_function_values_dx(2, 1) = -0.5;
    shape_function_values_dx(2, 2) = 1.0;

    Matrix shape_function_values_dy(3, 3);
    shape_function_values_dy(0, 0) = 0.0;
    shape_function_values_dy(0, 1) = 1.0;
    shape_function_values_dy(0, 2) = 1.0;
    shape_function_values_dy(1, 0) = 2.0;
    shape_function_values_dy(1, 1) = -1.0;
    shape_function_values_dy(1, 2) = 0.0;
    shape_function_values_dy(2, 0) = 1.0;
    shape_function_values_dy(2, 1) = 1.0;
    shape_function_values_dy(2, 2) = -1.0;

    Matrix shape_function_values_dz(3, 3);
    shape_function_values_dz(0, 0) = 1.0;
    shape_function_values_dz(0, 1) = 2.0;
    shape_function_values_dz(0, 2) = 0.0;
    shape_function_values_dz(1, 0) = 0.0;
    shape_function_values_dz(1, 1) = -1.0;
    shape_function_values_dz(1, 2) = 3.0;
    shape_function_values_dz(2, 0) = 2.0;
    shape_function_values_dz(2, 1) = 0.0;
    shape_function_values_dz(2, 2) = 1.0;

    ComputeDivSigmaUtility div_sigma_utility;
    const Matrix divergence = div_sigma_utility.ComputeDivergence(
        sigma_values,
        shape_function_values,
        shape_function_values_dx,
        shape_function_values_dy,
        shape_function_values_dz);

    Matrix expected_divergence(3, 3);
    expected_divergence(0, 0) = 8.0;
    expected_divergence(0, 1) = 11.0;
    expected_divergence(0, 2) = 1.0;
    expected_divergence(1, 0) = 7.0;
    expected_divergence(1, 1) = -2.0;
    expected_divergence(1, 2) = 17.0;
    expected_divergence(2, 0) = 1.5;
    expected_divergence(2, 1) = 6.0;
    expected_divergence(2, 2) = 4.5;

    KRATOS_EXPECT_MATRIX_NEAR(divergence, expected_divergence, 1.0e-12);
}

} // namespace Kratos::Testing
