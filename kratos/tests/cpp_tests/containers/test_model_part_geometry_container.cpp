//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Teschemacher
//                   Pooyan Dadvand
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/line_3d_2.h"

#include "containers/model.h"

namespace Kratos {
namespace Testing {

    Line3D2<Node>::Pointer GenerateLineModelPartGeometryContainer() {
        PointerVector<Node> points;

        points.push_back(Node::Pointer(new Node(1, 0, 5, 0)));
        points.push_back(Node::Pointer(new Node(2, 5, 5, 0)));

        return Kratos::make_shared<Line3D2<Node>>(
            points
            );
    }

    ///// Test Geometry Container
    KRATOS_TEST_CASE_IN_SUITE(TestModelPartGeometryContainer, KratosCoreGeometryContainerFastSuite) {
        Model model;

        auto& model_part = model.CreateModelPart("Main");

        auto& model_part_lines = model_part.CreateSubModelPart("Lines");

        auto& model_part_no_lines = model_part.CreateSubModelPart("NoLines");

        auto& model_part_sub_1 = model_part.CreateSubModelPart("sub_1");
        auto& model_part_sub_2 = model_part.CreateSubModelPart("sub_2");
        auto& model_part_sub_sub_1 = model_part_sub_1.CreateSubModelPart("subsub_1");
        auto& model_part_sub_sub_2 = model_part_sub_2.CreateSubModelPart("subsub_2");

        auto p_line_1 = GenerateLineModelPartGeometryContainer();
        p_line_1->SetId(1);

        model_part_lines.AddGeometry(p_line_1);
        KRATOS_EXPECT_EQ(model_part_lines.NumberOfGeometries(), 1);
        model_part_lines.AddGeometry(p_line_1); // adding same geomerty does not fail
        KRATOS_EXPECT_EQ(model_part_lines.NumberOfGeometries(), 1);

        auto p_line_2 = GenerateLineModelPartGeometryContainer();
        p_line_2->SetId(1);

        // check correct error if multiple geometries with same id are added
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(
            model_part_lines.AddGeometry(p_line_2),
            "Error: Attempting to add Geometry with Id: 1, unfortunately a (different) geometry with the same Id already exists!");

        p_line_2->SetId(2);
        model_part_lines.AddGeometry(p_line_2);

        // check correct number of geometries
        KRATOS_EXPECT_EQ(model_part_lines.NumberOfGeometries(), 2);
        KRATOS_EXPECT_EQ(model_part.NumberOfGeometries(), 2);
        KRATOS_EXPECT_EQ(model_part_no_lines.NumberOfGeometries(), 0);

        // check adding to SubModelPart and SubSubModelPart
        // first to SubSub
        model_part_sub_sub_1.AddGeometry(p_line_1);
        KRATOS_EXPECT_EQ(model_part_sub_1.NumberOfGeometries(), 1);
        KRATOS_EXPECT_EQ(model_part_sub_sub_1.NumberOfGeometries(), 1);

        // first to Sub then to SubSub
        model_part_sub_2.AddGeometry(p_line_1);
        KRATOS_EXPECT_EQ(model_part_sub_2.NumberOfGeometries(), 1);
        KRATOS_EXPECT_EQ(model_part_sub_sub_2.NumberOfGeometries(), 0);
        model_part_sub_sub_2.AddGeometry(p_line_1);
        KRATOS_EXPECT_EQ(model_part_sub_2.NumberOfGeometries(), 1);
        KRATOS_EXPECT_EQ(model_part_sub_sub_2.NumberOfGeometries(), 1);

        // check adding with string
        auto p_line_3 = GenerateLineModelPartGeometryContainer();
        p_line_3->SetId("GeometryLine1");
        model_part_lines.AddGeometry(p_line_3);

        KRATOS_EXPECT_EQ(model_part_lines.NumberOfGeometries(), 3);

        // check if correct element is returned
        KRATOS_EXPECT_EQ(model_part_lines.GetGeometry(1).Id(), 1);
        KRATOS_EXPECT_EQ(model_part_lines.pGetGeometry(1)->Id(), 1);

        // check if correct element is returned from root
        KRATOS_EXPECT_EQ(model_part.GetGeometry(1).Id(), 1);
        KRATOS_EXPECT_EQ(model_part.pGetGeometry(1)->Id(), 1);

        // check remove functions
        model_part_lines.RemoveGeometry("GeometryLine1");
        model_part_lines.RemoveGeometryFromAllLevels(1);
        KRATOS_EXPECT_EQ(model_part_lines.NumberOfGeometries(), 1);
        KRATOS_EXPECT_EQ(model_part.NumberOfGeometries(), 2);
        KRATOS_EXPECT_EQ(model_part.NumberOfGeometries(), 2);

        // check if correct geometries are removed
        KRATOS_EXPECT_FALSE(model_part_lines.HasGeometry("GeometryLine1"));
        KRATOS_EXPECT_TRUE(model_part.HasGeometry("GeometryLine1"));

        // check if correct geometries are removed
        KRATOS_EXPECT_FALSE(model_part_lines.HasGeometry(1));
        KRATOS_EXPECT_FALSE(model_part.HasGeometry(1));
    }
} // namespace Testing.
} // namespace Kratos.
