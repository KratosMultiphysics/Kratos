//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Suneth Warnakulasuriya
//
//

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/expect.h"
#include "utilities/integration_utilities.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/hexahedra_3d_8.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(IntegrationUtilitiesComputeArea2DGeometry, KratosCoreFastSuite)
{
    Model model;
    auto& model_part = model.CreateModelPart("test");

    // create new nodes.
    auto p_node1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto p_node2 = model_part.CreateNewNode(2, 3.0, 0.0, 0.0);
    auto p_node3 = model_part.CreateNewNode(3, 3.0, 3.0, 0.0);
    auto p_node4 = model_part.CreateNewNode(4, 0.0, 3.0, 0.0);

    // create new geometry.
    auto p_geometry = Kratos::make_shared<Quadrilateral2D4<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({p_node1, p_node2, p_node3, p_node4})});
    auto& r_geometry = *p_geometry;

    // calculates reference value.
    const double ref_value = IntegrationUtilities::ComputeArea2DGeometry(r_geometry);

    KRATOS_EXPECT_NEAR(ref_value, 9.0, 1e-9);

    // now calculate the derivatives
    for (IndexType i = 0; i < 4; ++i) {
        auto& coordinates = r_geometry[i].Coordinates();
        for (IndexType j = 0; j < 2; ++j) {
            const double analytical_derivative = IntegrationUtilities::ComputeArea2DGeometryDerivative(i, j, r_geometry);

            // finite difference derivative
            coordinates[j] += 1e-4;
            const double value = IntegrationUtilities::ComputeArea2DGeometry(r_geometry);
            coordinates[j] -= 1e-4;

            const double fd_derivative = (value - ref_value) / 1e-4;

            KRATOS_EXPECT_NEAR(analytical_derivative, fd_derivative, 1e-9);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(IntegrationUtilitiesComputeVolume3DGeometryDerivative, KratosCoreFastSuite)
{
    Model model;
    auto& model_part = model.CreateModelPart("test");

    // create new nodes.
    auto p_node1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto p_node2 = model_part.CreateNewNode(2, 3.0, 0.0, 0.0);
    auto p_node3 = model_part.CreateNewNode(3, 3.0, 3.0, 0.0);
    auto p_node4 = model_part.CreateNewNode(4, 0.0, 3.0, 0.0);
    auto p_node5 = model_part.CreateNewNode(5, 0.0, 0.0, 3.0);
    auto p_node6 = model_part.CreateNewNode(6, 3.0, 0.0, 3.0);
    auto p_node7 = model_part.CreateNewNode(7, 3.0, 3.0, 3.0);
    auto p_node8 = model_part.CreateNewNode(8, 0.0, 3.0, 3.0);

    // create new geometry.
    auto p_geometry = Kratos::make_shared<Hexahedra3D8<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({p_node1, p_node2, p_node3, p_node4, p_node5, p_node6, p_node7, p_node8})});
    auto& r_geometry = *p_geometry;

    // calculates reference value.
    const double ref_value = IntegrationUtilities::ComputeVolume3DGeometry(r_geometry);

    KRATOS_EXPECT_NEAR(ref_value, 27.0, 1e-9);

    // now calculate the derivatives
    for (IndexType i = 0; i < 8; ++i) {
        auto& coordinates = r_geometry[i].Coordinates();
        for (IndexType j = 0; j < 3; ++j) {
            const double analytical_derivative = IntegrationUtilities::ComputeVolume3DGeometryDerivative(i, j, r_geometry);

            // finite difference derivative
            coordinates[j] += 1e-4;
            const double value = IntegrationUtilities::ComputeVolume3DGeometry(r_geometry);
            coordinates[j] -= 1e-4;

            const double fd_derivative = (value - ref_value) / 1e-4;

            KRATOS_EXPECT_NEAR(analytical_derivative, fd_derivative, 1e-9);
        }
    }
}

}  // namespace Kratos::Testing.
