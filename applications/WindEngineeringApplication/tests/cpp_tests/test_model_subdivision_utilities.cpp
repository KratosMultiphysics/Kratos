//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

// Project Includes
#include "containers/model.h"
#include "testing/testing.h"

// Utility includes
#include "custom_utilities/model_subdivision_utilities.h"


namespace Kratos::Testing {


KRATOS_TEST_CASE_IN_SUITE(ModelSubdivisionUtilitiesSortNodesBySlabs, KratosWindEngineeringFastSuite)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("main");

    // Fill a cartesian domain with nodes
    // Domain corners: {-2.0, -2.0, -2.0} {2.0, 2.0, 2.0}
    const std::size_t resolution = 6;
    const double distance = 4.0 / (resolution - 1);
    std::size_t node_counter = 0;

    for (std::size_t i=0; i<resolution; ++i) {
        for (std::size_t j=0; j<resolution; ++j) {
            for (std::size_t k=0; k<resolution; ++k) {
                r_model_part.CreateNewNode(
                    ++node_counter,
                    -2.0 + i * distance,
                    -2.0 + j * distance,
                    -2.0 + k * distance
                );
            }
        }
    }

    // Sort nodes by z:
    //  z in [-1, -1/3[     -> 0
    //  z in [-1/3, 1/3[    -> 1 (empty)
    //  z in [1/3, 1]       -> 2
    auto sub_model_parts = Wind::ModelSubdivisionUtilities::SortNodesBySlabs(
        r_model_part,
        array_1d<double,3>({0.0, 0.0, -1.0}),
        array_1d<double,3>({0.0, 0.0, 1.0}),
        3,
        false
    );

    KRATOS_EXPECT_EQ(sub_model_parts.size(), 3);

    KRATOS_EXPECT_NE(sub_model_parts[0], nullptr);
    KRATOS_EXPECT_EQ(sub_model_parts[0]->Nodes().size(), resolution * resolution);
    for (const auto& r_node : sub_model_parts[0]->Nodes()) {
        KRATOS_EXPECT_TRUE(-1.0 <= r_node.Z());
        KRATOS_EXPECT_TRUE(r_node.Z() < -1.0/3.0);
    }

    KRATOS_EXPECT_NE(sub_model_parts[1], nullptr);
    KRATOS_EXPECT_EQ(sub_model_parts[1]->Nodes().size(), 0);

    KRATOS_EXPECT_NE(sub_model_parts[2], nullptr);
    KRATOS_EXPECT_EQ(sub_model_parts[2]->Nodes().size(), resolution * resolution);
    for (const auto& r_node : sub_model_parts[2]->Nodes()) {
        KRATOS_EXPECT_TRUE(1.0/3.0 <= r_node.Z());
        KRATOS_EXPECT_TRUE(r_node.Z() <= 1.0);
    }
}


// Make protected content accessible
struct MockUtility : public Wind::ModelSubdivisionUtilities
{
    using Wind::ModelSubdivisionUtilities::Slab;
    using Wind::ModelSubdivisionUtilities::SlabStack;
};


KRATOS_TEST_CASE_IN_SUITE(ModelSubdivisionUtilitiesSlab, KratosCoreFastSuite)
{
    const array_1d<double,3> bottom({-1.0, -1.0, 0.0});
    const array_1d<double,3> top({1.0, 1.0, 0.0});

    // Check invalid construction
    auto make_degenerate_slab = [&](){MockUtility::Slab {bottom, bottom, false};};
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(make_degenerate_slab(), "Degenerate slab!");

    // Check point location
    const MockUtility::Slab closed_slab {bottom, top, false};
    const MockUtility::Slab open_slab {bottom, top, true};

    const array_1d<double,3> test_below({-1.0, -2.0, 100.0});
    const array_1d<double,3> test_bottom_boundary = bottom;
    const array_1d<double,3> test_inside({0.0, 0.0, 100.0});
    const array_1d<double,3> test_top_boundary = top;
    const array_1d<double,3> test_above({1.0, 2.0, 100.0});

    // Closed slab
    KRATOS_EXPECT_TRUE(closed_slab.IsBelow(test_below));
    KRATOS_EXPECT_TRUE(!closed_slab.IsInside(test_below));
    KRATOS_EXPECT_TRUE(!closed_slab.IsAbove(test_below));

    KRATOS_EXPECT_TRUE(!closed_slab.IsBelow(test_bottom_boundary));
    KRATOS_EXPECT_TRUE(closed_slab.IsInside(test_bottom_boundary));
    KRATOS_EXPECT_TRUE(!closed_slab.IsAbove(test_bottom_boundary));

    KRATOS_EXPECT_TRUE(!closed_slab.IsBelow(test_inside));
    KRATOS_EXPECT_TRUE(closed_slab.IsInside(test_inside));
    KRATOS_EXPECT_TRUE(!closed_slab.IsAbove(test_inside));

    KRATOS_EXPECT_TRUE(!closed_slab.IsBelow(test_top_boundary));
    KRATOS_EXPECT_TRUE(closed_slab.IsInside(test_top_boundary));
    KRATOS_EXPECT_TRUE(!closed_slab.IsAbove(test_top_boundary));

    KRATOS_EXPECT_TRUE(!closed_slab.IsBelow(test_above));
    KRATOS_EXPECT_TRUE(!closed_slab.IsInside(test_above));
    KRATOS_EXPECT_TRUE(closed_slab.IsAbove(test_above));

    // Open slab
    KRATOS_EXPECT_TRUE(open_slab.IsBelow(test_below));
    KRATOS_EXPECT_TRUE(!open_slab.IsInside(test_below));
    KRATOS_EXPECT_TRUE(!open_slab.IsAbove(test_below));

    KRATOS_EXPECT_TRUE(open_slab.IsBelow(test_bottom_boundary));
    KRATOS_EXPECT_TRUE(!open_slab.IsInside(test_bottom_boundary));
    KRATOS_EXPECT_TRUE(!open_slab.IsAbove(test_bottom_boundary));

    KRATOS_EXPECT_TRUE(!open_slab.IsBelow(test_inside));
    KRATOS_EXPECT_TRUE(open_slab.IsInside(test_inside));
    KRATOS_EXPECT_TRUE(!open_slab.IsAbove(test_inside));

    KRATOS_EXPECT_TRUE(!open_slab.IsBelow(test_top_boundary));
    KRATOS_EXPECT_TRUE(!open_slab.IsInside(test_top_boundary));
    KRATOS_EXPECT_TRUE(open_slab.IsAbove(test_top_boundary));

    KRATOS_EXPECT_TRUE(!open_slab.IsBelow(test_above));
    KRATOS_EXPECT_TRUE(!open_slab.IsInside(test_above));
    KRATOS_EXPECT_TRUE(open_slab.IsAbove(test_above));
}


KRATOS_TEST_CASE_IN_SUITE(ModelSubdivisionUtilitiesSlabStack, KratosCoreFastSuite)
{
    const array_1d<double,3> bottom {-1.0, -1.0, 0.0};
    const array_1d<double,3> top {1.0, 1.0, 0.0};

    // Check invalid construction
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(MockUtility::SlabStack(bottom, top, 0, false), "Number of slabs in a stack must be at least 1");

    // Check point location
    const MockUtility::SlabStack slab_stack(bottom, top, 4, false);

    KRATOS_EXPECT_EQ(slab_stack.size(), 4);

    const array_1d<double,3> test_below({-1.0, -3.0, 100.0});
    const array_1d<double,3> test_in_0({-7.0/4.0, 1.0/4.0, 100.0});
    const array_1d<double,3> test_in_1({-1.0/4.0, -1.0/4.0, 100.0});
    const array_1d<double,3> test_in_2({1.0/4.0, 1.0/4.0, 100.0});
    const array_1d<double,3> test_in_3({7.0/4.0, -1.0/4.0, 100.0});
    const array_1d<double,3> test_above({1.0, 3.0, 100.0});

    //KRATOS_EXPECT_EXCEPTION_IS_THROWN(slab_stack.SlabIndex(test_below), "");
    KRATOS_EXPECT_EQ(slab_stack.SlabIndex(test_in_0), 0);
    KRATOS_EXPECT_EQ(slab_stack.SlabIndex(test_in_1), 1);
    KRATOS_EXPECT_EQ(slab_stack.SlabIndex(test_in_2), 2);
    KRATOS_EXPECT_EQ(slab_stack.SlabIndex(test_in_3), 3);
    //KRATOS_EXPECT_EXCEPTION_IS_THROWN(slab_stack.SlabIndex(test_above), "");
}


} // namespace Kratos::Testing
