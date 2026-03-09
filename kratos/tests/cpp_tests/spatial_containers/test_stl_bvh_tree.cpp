//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Altair Engineering
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

} // namespace Kratos::Testing
