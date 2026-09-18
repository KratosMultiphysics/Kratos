//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "spatial_containers/stl_bvh_tree.h"
#include "containers/model.h"
#include "tests/test_utilities/cpp_tests_utilities.h"

namespace Kratos::Testing
{

/** Checks StlBvhTree search in radius
*/
KRATOS_TEST_CASE_IN_SUITE(StlBvhTreeSearchInRadius, KratosCoreFastSuite)
{
    Model current_model;

    // Generate the cube skin (12 triangles, half-extents 0.6, 0.9, 0.3)
    ModelPart& r_skin_part = CppTestsUtilities::CreateCubeSkinModelPart(current_model);

    StlBvhTree bvh(r_skin_part.ElementsBegin(), r_skin_part.ElementsEnd());

    array_1d<double, 3> center_point{0.0, 0.0, 0.0};

    auto count_in_radius = [&](double radius) {
        std::size_t count = 0;
        bvh.ApplyInRadius(center_point, radius, [&](GeometricalObject&) {
            ++count;
            return true; // continue search
        });
        return count;
    };

    KRATOS_EXPECT_EQ(count_in_radius(0.29), 0);
    KRATOS_EXPECT_EQ(count_in_radius(0.30), 4);
    KRATOS_EXPECT_EQ(count_in_radius(0.40), 4);
    KRATOS_EXPECT_EQ(count_in_radius(0.60), 8);
    KRATOS_EXPECT_EQ(count_in_radius(0.70), 8);
    KRATOS_EXPECT_EQ(count_in_radius(0.90), 12);
}

/** Checks StlBvhTree ApplyInRadius early exit
*/
KRATOS_TEST_CASE_IN_SUITE(StlBvhTreeSearchInRadiusEarlyExit, KratosCoreFastSuite)
{
    Model current_model;

    ModelPart& r_skin_part = CppTestsUtilities::CreateCubeSkinModelPart(current_model);

    StlBvhTree bvh(r_skin_part.ElementsBegin(), r_skin_part.ElementsEnd());

    array_1d<double, 3> center{0.0, 0.0, 0.0};

    // Stop after the first hit — total visits must be exactly 1
    std::size_t visits = 0;
    bvh.ApplyInRadius(center, 0.9, [&](GeometricalObject&) {
        ++visits;
        return false; // stop immediately
    });

    KRATOS_EXPECT_EQ(visits, 1);
}

/** Checks StlBvhTree with an empty model part
*/
KRATOS_TEST_CASE_IN_SUITE(StlBvhTreeEmptyModelPart, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_empty = current_model.CreateModelPart("Empty");

    StlBvhTree bvh(r_empty.ElementsBegin(), r_empty.ElementsEnd());

    array_1d<double, 3> center{0.0, 0.0, 0.0};
    std::size_t visits = 0;
    bvh.ApplyInRadius(center, 1.0, [&](GeometricalObject&) { ++visits; return true; });

    KRATOS_EXPECT_EQ(visits, 0);
}

/** Checks StlBvhTree SearchNearest returns the closest triangle and correct distance
*/
KRATOS_TEST_CASE_IN_SUITE(StlBvhTreeSearchNearest, KratosCoreFastSuite)
{
    constexpr double tolerance = 1e-5; // float-precision BVH, so tolerance is relaxed vs double bins

    Model current_model;

    // Cube half-extents: x=0.6, y=0.9, z=0.3  ->  nearest face to the origin is at z=0.3
    const double cube_z = 0.3;
    ModelPart& r_skin_part = CppTestsUtilities::CreateCubeSkinModelPart(current_model, 0.6, 0.9, cube_z);

    StlBvhTree bvh(r_skin_part.ElementsBegin(), r_skin_part.ElementsEnd());

    const double epsilon = 1.0e-6;
    array_1d<double, 3> near_point{epsilon, epsilon, epsilon};

    auto result = bvh.SearchNearest(near_point);

    KRATOS_EXPECT_NE(result.p_object, nullptr);
    KRATOS_EXPECT_NEAR(result.Distance, cube_z - epsilon, tolerance);
}

/** Checks StlBvhTree SearchNearest on an empty tree returns no object
*/
KRATOS_TEST_CASE_IN_SUITE(StlBvhTreeSearchNearestEmpty, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_empty = current_model.CreateModelPart("Empty");

    StlBvhTree bvh(r_empty.ElementsBegin(), r_empty.ElementsEnd());

    array_1d<double, 3> center{0.0, 0.0, 0.0};
    auto result = bvh.SearchNearest(center);

    KRATOS_EXPECT_EQ(result.p_object, nullptr);
    KRATOS_EXPECT_EQ(result.Distance, std::numeric_limits<double>::max());
}

} // namespace Kratos::Testing
