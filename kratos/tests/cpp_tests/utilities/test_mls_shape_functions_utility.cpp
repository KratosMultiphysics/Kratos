//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>


// External includes

// Project includes
#include "utilities/mls_shape_functions_utility.h"
#include "testing/testing.h"

namespace Kratos
{
namespace Testing
{

    void SetPointsCoordinatesMatrix(Matrix& rPointsCoordinates)
    {
        rPointsCoordinates = ZeroMatrix(4,3);
        rPointsCoordinates(0,0) = -1.0; rPointsCoordinates(0,1) = -1.0;
        rPointsCoordinates(1,0) = 1.0; rPointsCoordinates(1,1) = -1.0;
        rPointsCoordinates(2,0) = 1.0; rPointsCoordinates(2,1) = 1.0;
        rPointsCoordinates(3,0) = -1.0; rPointsCoordinates(3,1) = 1.0;
    }

    KRATOS_TEST_CASE_IN_SUITE(MLSShapeFunctionsUtilityCalculateKernelCoordinates, KratosCoreFastSuite)
    {
        // Set the points cloud
        Matrix pt_coords;
        SetPointsCoordinatesMatrix(pt_coords);

        // Set the point of interest
        const array_1d<double,3> ref_pt = ZeroVector(3);

        // Calculate kernel values
        const double h = 1.0;
        const double tol = 1.0e-12;
        const std::array<double, 4> expected_values = {0.0430785586037,0.0430785586037,0.0430785586037,0.0430785586037};
        for (std::size_t i_pt = 0; i_pt < 4; ++i_pt) {
            const double kernel = MLSShapeFunctionsUtility::CalculateKernel(ref_pt - row(pt_coords, i_pt), h);
            KRATOS_CHECK_NEAR(kernel, expected_values[i_pt], tol);
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(MLSShapeFunctionsUtilityCalculateKernelZeroDistance, KratosCoreFastSuite)
    {
        // Calculate kernel value
        const double h = 0.25;
        const double tol = 1.0e-12;
        const double kernel = MLSShapeFunctionsUtility::CalculateKernel(ZeroVector(3), h);
        KRATOS_CHECK_NEAR(kernel, 5.09295817894, tol);
    }

    KRATOS_TEST_CASE_IN_SUITE(MLSShapeFunctionsUtilityCalculateKernelDerivative, KratosCoreFastSuite)
    {
        // Set the points cloud
        Matrix pt_coords;
        SetPointsCoordinatesMatrix(pt_coords);

        // Set the point of interest
        const array_1d<double,3> ref_pt = ZeroVector(3);

        // Calculate kernel values
        const double h = 1.0;
        const double tol = 1.0e-12;
        array_1d<double,2> kernel_der;
        const std::array<double, 8> expected_values = {
            -0.0861571172074, -0.0861571172074,
            0.0861571172074, -0.0861571172074,
            0.0861571172074, 0.0861571172074,
            -0.0861571172074,0.0861571172074};
        for (std::size_t i_pt = 0; i_pt < 4; ++i_pt) {
            MLSShapeFunctionsUtility::CalculateKernelDerivative<2>(ref_pt - row(pt_coords, i_pt), h, kernel_der);
            KRATOS_CHECK_NEAR(kernel_der[0], expected_values[2*i_pt], tol);
            KRATOS_CHECK_NEAR(kernel_der[1], expected_values[2*i_pt+1], tol);
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(MLSShapeFunctionsUtilityCalculateShapeFunctionsAndGradients, KratosCoreFastSuite)
    {
        // Set the points cloud
        Matrix pt_coords;
        SetPointsCoordinatesMatrix(pt_coords);

        // Set the point of interest
        const array_1d<double,3> ref_pt = ZeroVector(3);

        // Calculate kernel values
        const double h = 1.0;
        Vector N_container;
        Matrix DN_DX_container;
        MLSShapeFunctionsUtility::CalculateShapeFunctionsAndGradients(
            pt_coords,
            ref_pt,
            h,
            N_container,
            DN_DX_container);

        // Check results
        const double tol = 1.0e-12;
        const std::array<double, 4> expected_N = {0.25,0.25,0.25,0.25};
        const std::array<double, 8> expected_DN_DX = {-0.25, -0.25, 0.25, -0.25, 0.25, 0.25, -0.25,0.25};
        for (std::size_t i_pt = 0; i_pt < 4; ++i_pt) {
            KRATOS_CHECK_NEAR(N_container(i_pt), expected_N[i_pt], tol);
            KRATOS_CHECK_NEAR(DN_DX_container(i_pt, 0), expected_DN_DX[2*i_pt], tol);
            KRATOS_CHECK_NEAR(DN_DX_container(i_pt, 1), expected_DN_DX[2*i_pt+1], tol);
        }
    }

// std::cout << std::setprecision(12) << kernel_der[0] << std::endl;
// std::cout << std::setprecision(12) << kernel_der[1] << std::endl;

} // namespace Testing
}  // namespace Kratos.
