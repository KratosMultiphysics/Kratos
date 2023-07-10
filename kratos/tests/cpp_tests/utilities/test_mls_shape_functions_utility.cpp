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

// External includes

// Project includes
#include "containers/model.h"
#include "geometries/quadrilateral_2d_4.h"
#include "processes/structured_mesh_generator_process.h"
#include "utilities/mls_shape_functions_utility.h"
#include "testing/testing.h"

namespace Kratos
{
namespace Testing
{

    void SetPointsCoordinatesMatrix2D(Matrix& rPointsCoordinates)
    {
        rPointsCoordinates = ZeroMatrix(4,3);
        rPointsCoordinates(0,0) = -1.0; rPointsCoordinates(0,1) = -1.0;
        rPointsCoordinates(1,0) = 1.0; rPointsCoordinates(1,1) = -1.0;
        rPointsCoordinates(2,0) = 1.0; rPointsCoordinates(2,1) = 1.0;
        rPointsCoordinates(3,0) = -1.0; rPointsCoordinates(3,1) = 1.0;
    }

    void SetPointsCoordinatesMatrix3D(Matrix& rPointsCoordinates)
    {
        rPointsCoordinates = ZeroMatrix(8,3);
        rPointsCoordinates(0,0) = -1.0; rPointsCoordinates(0,1) = -1.0; rPointsCoordinates(0,2) = -1.0;
        rPointsCoordinates(1,0) = 1.0; rPointsCoordinates(1,1) = -1.0; rPointsCoordinates(1,2) = -1.0;
        rPointsCoordinates(2,0) = 1.0; rPointsCoordinates(2,1) = 1.0; rPointsCoordinates(2,2) = -1.0;
        rPointsCoordinates(3,0) = -1.0; rPointsCoordinates(3,1) = 1.0; rPointsCoordinates(3,2) = -1.0;
        rPointsCoordinates(4,0) = -1.0; rPointsCoordinates(4,1) = -1.0; rPointsCoordinates(4,2) = 1.0;
        rPointsCoordinates(5,0) = 1.0; rPointsCoordinates(5,1) = -1.0; rPointsCoordinates(5,2) = 1.0;
        rPointsCoordinates(6,0) = 1.0; rPointsCoordinates(6,1) = 1.0; rPointsCoordinates(6,2) = 1.0;
        rPointsCoordinates(7,0) = -1.0; rPointsCoordinates(7,1) = 1.0; rPointsCoordinates(7,2) = 1.0;
    }

    void MLSShapeFunctionsUtilityConsistencyCheck(
        std::function<double(array_1d<double,3>)>& rAnalyticalSolutionFunction,
        std::function<array_1d<double,2>(array_1d<double,3>)>& rAnalyticalSolutionGradientFunction,
        std::function<void(const Matrix&, const array_1d<double,3>&, const double, Vector&, Matrix&)>& rMLSShapeFunctionsAndGradientsFunction)
    {
        // Set the background mesh model part
        Model model;
        auto& r_model_part = model.CreateModelPart("TestModelPart");
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
        r_model_part.AddNodalSolutionStepVariable(DENSITY);

        // Generate the background mesh (done with the StructuredMeshGeneratorProcess)
        auto p_point_1 = Kratos::make_intrusive<Node>(1, -0.5, -0.5,  0.0);
        auto p_point_2 = Kratos::make_intrusive<Node>(2, -0.5,  0.5,  0.0);
        auto p_point_3 = Kratos::make_intrusive<Node>(3,  0.5,  0.5,  0.0);
        auto p_point_4 = Kratos::make_intrusive<Node>(4,  0.5, -0.5,  0.0);
        Quadrilateral2D4<Node> geometry(p_point_1, p_point_2, p_point_3, p_point_4);
        Parameters mesher_parameters(R"(
        {
            "number_of_divisions": 3,
            "element_name": "Element2D3N",
            "condition_name": "LineCondition"
        })");
        StructuredMeshGeneratorProcess(geometry, r_model_part, mesher_parameters).Execute();

        // Set the linear field (x+y) to be interpolated and the coordinates matrix
        const std::size_t n_nodes = r_model_part.NumberOfNodes();
        Matrix nodes_coords = ZeroMatrix(n_nodes, 3);
        std::size_t i = 0;
        for (auto& r_node : r_model_part.Nodes()) {
            nodes_coords(i,0) = r_node.X();
            nodes_coords(i,1) = r_node.Y();
            r_node.FastGetSolutionStepValue(DENSITY) = rAnalyticalSolutionFunction(r_node.Coordinates());
            ++i;
        }

        // Calculate the error of the MLS interpolation
        double val_error = 0.0;
        double val_dx_error = 0.0;
        double val_dy_error = 0.0;
        const double h = 1.0;
        array_1d<double,3> i_gauss_N;
        array_1d<double,3> i_gauss_coords;
        Vector MLS_N_container;
        Matrix MLS_DN_DX_container;
        for (auto& r_element : r_model_part.Elements()) {
            const auto& r_geom = r_element.GetGeometry();
            const auto N_container = r_geom.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_2);
            const std::size_t n_nodes = N_container.size1();
            const std::size_t n_gauss = N_container.size2();
            for (std::size_t i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
                // Get current Gauss pt. coordinates
                i_gauss_N = row(N_container, i_gauss);
                i_gauss_coords = ZeroVector(3);
                for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
                    noalias(i_gauss_coords) += i_gauss_N[i_node] * r_geom[i_node].Coordinates();
                }

                // Calculate the MLS basis in the current Gauss pt.
                rMLSShapeFunctionsAndGradientsFunction(
                    nodes_coords,
                    i_gauss_coords,
                    h,
                    MLS_N_container,
                    MLS_DN_DX_container);

                // Calculate the MLS interpolation at current Gauss pt.
                double mls_val = 0.0;
                double mls_grad_x = 0.0;
                double mls_grad_y = 0.0;
                std::size_t i_mls = 0;
                for (const auto& r_node : r_model_part.Nodes()) {
                    const double node_val = r_node.FastGetSolutionStepValue(DENSITY);
                    mls_val += MLS_N_container[i_mls] * node_val;
                    mls_grad_x += MLS_DN_DX_container(i_mls,0) * node_val;
                    mls_grad_y += MLS_DN_DX_container(i_mls,1) * node_val;
                    ++i_mls;
                }

                // Calculate the error at current Gauss pt.
                const double sol = rAnalyticalSolutionFunction(i_gauss_coords);
                const array_1d<double,2> grad_sol = rAnalyticalSolutionGradientFunction(i_gauss_coords);
                const double weight = r_geom.Area() / n_gauss;
                val_error += weight * (sol - mls_val);
                val_dx_error += weight * (grad_sol[0] - mls_grad_x);
                val_dy_error += weight * (grad_sol[1] - mls_grad_y);
            }
        }

        const double zero_tol = 1.0e-15;
        KRATOS_CHECK_LESS_EQUAL(std::abs(val_error), zero_tol);
        KRATOS_CHECK_LESS_EQUAL(std::abs(val_dx_error), zero_tol);
        KRATOS_CHECK_LESS_EQUAL(std::abs(val_dy_error), zero_tol);
    }

    KRATOS_TEST_CASE_IN_SUITE(MLSShapeFunctionsUtilityCalculateKernelCoordinates, KratosCoreFastSuite)
    {
        // Set the points cloud
        Matrix pt_coords;
        SetPointsCoordinatesMatrix2D(pt_coords);

        // Set the point of interest
        const array_1d<double,3> ref_pt = ZeroVector(3);

        // Calculate kernel values
        const double h = 1.0;
        const double tol = 1.0e-12;
        const std::array<double, 4> expected_values = {0.367879441171,0.367879441171,0.367879441171,0.367879441171};
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
        KRATOS_CHECK_NEAR(kernel, 1.0, tol);
    }

    KRATOS_TEST_CASE_IN_SUITE(MLSShapeFunctionsUtilityCalculateKernelDerivative, KratosCoreFastSuite)
    {
        // Set the points cloud
        Matrix pt_coords;
        SetPointsCoordinatesMatrix2D(pt_coords);

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

    KRATOS_TEST_CASE_IN_SUITE(MLSShapeFunctionsUtilityCalculateShapeFunctions2D, KratosCoreFastSuite)
    {
        // Set the points cloud
        Matrix pt_coords;
        SetPointsCoordinatesMatrix2D(pt_coords);

        // Set the point of interest
        const array_1d<double,3> ref_pt = ZeroVector(3);

        // Calculate kernel values
        const double h = 1.0;
        Vector N_container;
        MLSShapeFunctionsUtility::CalculateShapeFunctions<2,1>(
            pt_coords,
            ref_pt,
            h,
            N_container);

        // Check results
        const double tol = 1.0e-12;
        const std::array<double, 4> expected_N = {0.25,0.25,0.25,0.25};
        for (std::size_t i_pt = 0; i_pt < 4; ++i_pt) {
            KRATOS_CHECK_NEAR(N_container(i_pt), expected_N[i_pt], tol);
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(MLSShapeFunctionsUtilityCalculateShapeFunctionsAndGradients2D, KratosCoreFastSuite)
    {
        // Set the points cloud
        Matrix pt_coords;
        SetPointsCoordinatesMatrix2D(pt_coords);

        // Set the point of interest
        const array_1d<double,3> ref_pt = ZeroVector(3);

        // Calculate kernel values
        const double h = 1.0;
        Vector N_container;
        Matrix DN_DX_container;
        MLSShapeFunctionsUtility::CalculateShapeFunctionsAndGradients<2,1>(
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

    KRATOS_TEST_CASE_IN_SUITE(MLSShapeFunctionsUtilityCalculateShapeFunctions3D, KratosCoreFastSuite)
    {
        // Set the points cloud
        Matrix pt_coords;
        SetPointsCoordinatesMatrix3D(pt_coords);

        // Set the point of interest
        const array_1d<double,3> ref_pt = ZeroVector(3);

        // Calculate kernel values
        const double h = 1.0;
        Vector N_container;
        MLSShapeFunctionsUtility::CalculateShapeFunctions<3,1>(
            pt_coords,
            ref_pt,
            h,
            N_container);

        // Check results
        const double tol = 1.0e-12;
        const std::array<double, 8> expected_N = {0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125};
        for (std::size_t i_pt = 0; i_pt < 8; ++i_pt) {
            KRATOS_CHECK_NEAR(N_container(i_pt), expected_N[i_pt], tol);
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(MLSShapeFunctionsUtilityCalculateShapeFunctionsAndGradients3D, KratosCoreFastSuite)
    {
        // Set the points cloud
        Matrix pt_coords;
        SetPointsCoordinatesMatrix3D(pt_coords);

        // Set the point of interest
        const array_1d<double,3> ref_pt = ZeroVector(3);

        // Calculate kernel values
        const double h = 1.0;
        Vector N_container;
        Matrix DN_DX_container;
        MLSShapeFunctionsUtility::CalculateShapeFunctionsAndGradients<3,1>(
            pt_coords,
            ref_pt,
            h,
            N_container,
            DN_DX_container);

        // Check results
        const double tol = 1.0e-12;
        const std::array<double, 8> expected_N = {0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125};
        const std::array<double, 24> expected_DN_DX = {
            -0.125, -0.125, -0.125,
            0.125, -0.125, -0.125,
            0.125, 0.125, -0.125,
            -0.125,0.125, -0.125,
            -0.125, -0.125, 0.125,
            0.125, -0.125, 0.125,
            0.125, 0.125, 0.125,
            -0.125,0.125, 0.125};
        for (std::size_t i_pt = 0; i_pt < 8; ++i_pt) {
            KRATOS_CHECK_NEAR(N_container(i_pt), expected_N[i_pt], tol);
            KRATOS_CHECK_NEAR(DN_DX_container(i_pt, 0), expected_DN_DX[3*i_pt], tol);
            KRATOS_CHECK_NEAR(DN_DX_container(i_pt, 1), expected_DN_DX[3*i_pt+1], tol);
            KRATOS_CHECK_NEAR(DN_DX_container(i_pt, 2), expected_DN_DX[3*i_pt+2], tol);
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(MLSShapeFunctionsUtility1stOrderConsistency, KratosCoreFastSuite)
    {
        // Set the analytical solution and its gradient functions
        std::function<double(array_1d<double,3>)> sol_func = [](const array_1d<double,3>& rCoords){return rCoords[0] + rCoords[1];};
        std::function<array_1d<double,2>(array_1d<double,3>)> grad_sol_func = [](const array_1d<double,3>& rCoords){
            array_1d<double,2> grad_sol;
            grad_sol[0] = 1.0;
            grad_sol[1] = 1.0;
            return grad_sol;
        };

        // Set the MLS shape functions and gradients call
        std::function<void(const Matrix&, const array_1d<double,3>&, const double, Vector&, Matrix&)> MLS_call_1st_order = [&](
            const Matrix& rPoints,
            const array_1d<double,3>& rX,
            const double h, Vector& rN,
            Matrix& rDNDX)
            {
                MLSShapeFunctionsUtility::CalculateShapeFunctionsAndGradients<2,1>(
                    rPoints,
                    rX,
                    h,
                    rN,
                    rDNDX);
            };

        // Call the auxiliary function performing the check for the fields above
        MLSShapeFunctionsUtilityConsistencyCheck(sol_func, grad_sol_func, MLS_call_1st_order);
    }

    KRATOS_TEST_CASE_IN_SUITE(MLSShapeFunctionsUtility2ndOrderConsistency, KratosCoreFastSuite)
    {
        // Set the analytical solution and its gradient functions
        std::function<double(array_1d<double,3>)> sol_func = [](const array_1d<double,3>& rCoords){return std::pow(rCoords[0],2) + std::pow(rCoords[1],2);};
        std::function<array_1d<double,2>(array_1d<double,3>)> grad_sol_func = [](const array_1d<double,3>& rCoords){
            array_1d<double,2> grad_sol;
            grad_sol[0] = 2.0*rCoords[0];
            grad_sol[1] = 2.0*rCoords[1];
            return grad_sol;
        };

        // Set the MLS shape functions and gradients call
        std::function<void(const Matrix&, const array_1d<double,3>&, const double, Vector&, Matrix&)> MLS_call_2nd_order = [&](
            const Matrix& rPoints,
            const array_1d<double,3>& rX,
            const double h, Vector& rN,
            Matrix& rDNDX)
            {
                MLSShapeFunctionsUtility::CalculateShapeFunctionsAndGradients<2,2>(
                    rPoints,
                    rX,
                    h,
                    rN,
                    rDNDX);
            };

        // Call the auxiliary function performing the check for the fields above
        MLSShapeFunctionsUtilityConsistencyCheck(sol_func, grad_sol_func, MLS_call_2nd_order);
    }

} // namespace Testing
}  // namespace Kratos.
