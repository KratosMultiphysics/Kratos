//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//
//

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "utilities/brute_force_point_locator.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(BruteForcePointLocatorTriangleElement, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& model_part = current_model.CreateModelPart("Triangle");
    model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);

    model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    Properties::Pointer p_properties(new Properties(0));
    const int elem_id_to_be_found = 43;
    model_part.CreateNewElement("Element2D3N", elem_id_to_be_found, {1, 2, 3}, p_properties);

    BruteForcePointLocator point_locator(model_part);

    Point the_point(0.1, 0.25, 0.0);

    Vector shape_function_values;

    const int found_id = point_locator.FindElement(the_point, shape_function_values);

    KRATOS_CHECK_EQUAL(found_id, elem_id_to_be_found);
    KRATOS_CHECK_EQUAL(shape_function_values.size(), 3);

    KRATOS_CHECK_NEAR(shape_function_values[0], 0.65, 1e-06);
    KRATOS_CHECK_NEAR(shape_function_values[1], 0.1,  1e-06);
    KRATOS_CHECK_NEAR(shape_function_values[2], 0.25, 1e-06);
}

KRATOS_TEST_CASE_IN_SUITE(BruteForcePointLocatorQuadrilateralElement, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& model_part = current_model.CreateModelPart("Quadrilateral");
    model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);

    model_part.CreateNewNode(1, 0.0, 0.0,  0.0);
    model_part.CreateNewNode(2, 1.2, 0.1,  0.0);
    model_part.CreateNewNode(3, 1.3, 1.02, 0.0);
    model_part.CreateNewNode(4, 0.1, 1.0,  0.0);

    Properties::Pointer p_properties(new Properties(0));
    const int elem_id_to_be_found = 181;
    model_part.CreateNewElement("Element2D4N", elem_id_to_be_found, {1, 2, 3, 4}, p_properties);

    BruteForcePointLocator point_locator(model_part);

    Point the_point(0.13, 0.52, 0.0);

    Vector shape_function_values;

    const int found_id = point_locator.FindElement(the_point, shape_function_values);

    KRATOS_CHECK_EQUAL(found_id, elem_id_to_be_found);
    KRATOS_CHECK_EQUAL(shape_function_values.size(), 4);

    KRATOS_CHECK_NEAR(shape_function_values[0], 0.452231, 1e-06);
    KRATOS_CHECK_NEAR(shape_function_values[1], 0.0316039,  1e-06);
    KRATOS_CHECK_NEAR(shape_function_values[2], 0.0337157, 1e-06);
    KRATOS_CHECK_NEAR(shape_function_values[3], 0.48245, 1e-06);
}

KRATOS_TEST_CASE_IN_SUITE(BruteForcePointLocatorTetrahedraElement, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& model_part = current_model.CreateModelPart("Tetrahedral");
    model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

    model_part.CreateNewNode(1, 0.0, 0.04, 0.0);
    model_part.CreateNewNode(2, 1.01, 0.1, 0.2);
    model_part.CreateNewNode(3, 1.3, 1.02, 0.0);
    model_part.CreateNewNode(4, 0.1, 0.03, 1.09);

    Properties::Pointer p_properties(new Properties(0));
    const int elem_id_to_be_found = 23;
    model_part.CreateNewElement("Element3D4N", elem_id_to_be_found, {1, 2, 3, 4}, p_properties);

    BruteForcePointLocator point_locator(model_part);

    Point the_point(0.25, 0.12, 0.32);

    Vector shape_function_values;

    const int found_id = point_locator.FindElement(the_point, shape_function_values);

    KRATOS_CHECK_EQUAL(found_id, elem_id_to_be_found);
    KRATOS_CHECK_EQUAL(shape_function_values.size(), 4);

    KRATOS_CHECK_NEAR(shape_function_values[0], 0.530166, 1e-06);
    KRATOS_CHECK_NEAR(shape_function_values[1], 0.121616,  1e-06);
    KRATOS_CHECK_NEAR(shape_function_values[2], 0.0769547, 1e-06);
    KRATOS_CHECK_NEAR(shape_function_values[3], 0.271263, 1e-06);
}

KRATOS_TEST_CASE_IN_SUITE(BruteForcePointLocatorHexahedraElement, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& model_part = current_model.CreateModelPart("Hexahedral");
    model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

    model_part.CreateNewNode(1, 0.0, 0.04, 0.0);
    model_part.CreateNewNode(2, 1.01, 0.1, 0.2);
    model_part.CreateNewNode(3, 1.3, 1.02, 0.0);
    model_part.CreateNewNode(4, 0.1, 1.05, 0.09);
    model_part.CreateNewNode(5, 0.02, 0.03, 1.09);
    model_part.CreateNewNode(6, 1.08, 0.03, 1.09);
    model_part.CreateNewNode(7, 1.1, 1.03, 1.22);
    model_part.CreateNewNode(8, 0.1, 0.97, 1.09);

    Properties::Pointer p_properties(new Properties(0));
    const int elem_id_to_be_found = 69;
    model_part.CreateNewElement("Element3D8N", elem_id_to_be_found, {1,2,3,4,5,6,7,8}, p_properties);

    BruteForcePointLocator point_locator(model_part);

    Point the_point(0.28, 0.32, 0.72);

    Vector shape_function_values;

    const int found_id = point_locator.FindElement(the_point, shape_function_values);

    KRATOS_CHECK_EQUAL(found_id, elem_id_to_be_found);
    KRATOS_CHECK_EQUAL(shape_function_values.size(), 8);

    KRATOS_CHECK_NEAR(shape_function_values[0], 0.197609, 1e-06);
    KRATOS_CHECK_NEAR(shape_function_values[1], 0.0590793,  1e-06);
    KRATOS_CHECK_NEAR(shape_function_values[2], 0.0242583, 1e-06);
    KRATOS_CHECK_NEAR(shape_function_values[3], 0.0811397, 1e-06);
    KRATOS_CHECK_NEAR(shape_function_values[4], 0.348142, 1e-06);
    KRATOS_CHECK_NEAR(shape_function_values[5], 0.104084, 1e-06);
    KRATOS_CHECK_NEAR(shape_function_values[6], 0.0427376, 1e-06);
    KRATOS_CHECK_NEAR(shape_function_values[7], 0.142949, 1e-06);
}

KRATOS_TEST_CASE_IN_SUITE(BruteForcePointLocatorNode, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& model_part = current_model.CreateModelPart("ForTest");

    // Fill the model part geometry data
    const int node_id_to_be_found = 517;
    const double coord_x = 11.258;
    const double coord_y = -789.2368;
    const double coord_z = 0.863;
    model_part.CreateNewNode(12, 0.0, 0.1002, 0.0);
    model_part.CreateNewNode(24, 1.0, 0.0, 47.421);
    model_part.CreateNewNode(node_id_to_be_found, coord_x, coord_y, coord_z);
    model_part.CreateNewNode(5123, coord_x, coord_y, coord_z+0.001);

    BruteForcePointLocator point_locator(model_part);

    Point the_point(coord_x, coord_y, coord_z);

    const int found_id = point_locator.FindNode(the_point);

    KRATOS_CHECK_EQUAL(found_id, node_id_to_be_found);
}

KRATOS_TEST_CASE_IN_SUITE(BruteForcePointLocatorQuadrilateralCondition, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& model_part = current_model.CreateModelPart("Quadrilateral");
    model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

    model_part.CreateNewNode(1, 0.0, 0.0,  0.0);
    model_part.CreateNewNode(2, 1.2, 0.1,  0.0);
    model_part.CreateNewNode(3, 1.3, 1.02, 0.0);
    model_part.CreateNewNode(4, 0.1, 1.0,  0.0);

    Properties::Pointer p_properties(new Properties(0));
    const int cond_id_to_be_found = 86;
    model_part.CreateNewCondition("SurfaceCondition3D4N", cond_id_to_be_found, {1, 2, 3, 4}, p_properties);

    BruteForcePointLocator point_locator(model_part);

    Point the_point(0.13, 0.52, 0.0);

    Vector shape_function_values;

    const int found_id = point_locator.FindCondition(the_point, shape_function_values);

    KRATOS_CHECK_EQUAL(found_id, cond_id_to_be_found);
    KRATOS_CHECK_EQUAL(shape_function_values.size(), 4);

    KRATOS_CHECK_NEAR(shape_function_values[0], 0.452231, 1e-06);
    KRATOS_CHECK_NEAR(shape_function_values[1], 0.0316039,  1e-06);
    KRATOS_CHECK_NEAR(shape_function_values[2], 0.0337157, 1e-06);
    KRATOS_CHECK_NEAR(shape_function_values[3], 0.48245, 1e-06);
}

} // namespace Testing
}  // namespace Kratos.
