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
#include "geometries/triangle_2d_3.h"

#include "containers/model.h"

namespace Kratos::Testing {

namespace
{

    Line3D2<Node>::Pointer GenerateLineModelPartGeometryContainer() {
        PointerVector<Node> points;

        points.push_back(Node::Pointer(new Node(1, 0, 5, 0)));
        points.push_back(Node::Pointer(new Node(2, 5, 5, 0)));

        return Kratos::make_shared<Line3D2<Node>>(
            points
            );
    }

void SetUpTestModelPart(Model& rModel)
{
    auto& r_model_part = rModel.CreateModelPart("TestModelPart");
    r_model_part.CreateSubModelPart("SubModelPart1");
    auto& r_sub_model_part_2 = r_model_part.CreateSubModelPart("SubModelPart2");
    r_sub_model_part_2.CreateSubModelPart("SubSubModelPart");

    auto p_node_1 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto p_node_2 = r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    auto p_node_3 = r_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    auto p_node_4 = r_model_part.CreateNewNode(4, 1.0, 1.0, 0.0);

    rModel.GetModelPart("TestModelPart.SubModelPart1").CreateNewGeometry("Line2D2", 1, {{1, 2}});
    rModel.GetModelPart("TestModelPart.SubModelPart2").CreateNewGeometry("Triangle2D3", 2, {{1, 2, 4}});
    rModel.GetModelPart("TestModelPart.SubModelPart2.SubSubModelPart").CreateNewGeometry("Triangle2D3", 3, {{1, 4, 3}});
}

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

        // check that we do nothing if the same geometry (same id and connectivities) is to be added multiple times
        model_part_lines.AddGeometry(p_line_2);
        KRATOS_EXPECT_EQ(model_part_lines.NumberOfGeometries(), 1);

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

    KRATOS_TEST_CASE_IN_SUITE(TestModelPartGeometryContainerRepeatedGeometries1, KratosCoreGeometryContainerFastSuite)
    {
        // Set up test model part
        Model model;
        SetUpTestModelPart(model);

        // Add a same type and connectivity geometry with different Id (valid)
        model.GetModelPart("TestModelPart.SubModelPart2").CreateNewGeometry("Triangle2D3", 4, {1,4,3});

        // Check results
        KRATOS_EXPECT_EQ(model.GetModelPart("TestModelPart").NumberOfGeometries(), 4);
        KRATOS_EXPECT_EQ(model.GetModelPart("TestModelPart.SubModelPart1").NumberOfGeometries(), 1);
        KRATOS_EXPECT_EQ(model.GetModelPart("TestModelPart.SubModelPart2").NumberOfGeometries(), 3);
        KRATOS_EXPECT_EQ(model.GetModelPart("TestModelPart.SubModelPart2.SubSubModelPart").NumberOfGeometries(), 1);
    }

    KRATOS_TEST_CASE_IN_SUITE(TestModelPartGeometryContainerRepeatedGeometries2, KratosCoreGeometryContainerFastSuite)
    {
        // Set up test model part
        Model model;
        SetUpTestModelPart(model);

        // Check that adding the same type and connectivity geometry with same Id returns the existing one (valid)
        auto& r_model_part = model.GetModelPart("TestModelPart");
        auto p_geom_3_old = r_model_part.pGetGeometry(3);
        auto p_geom_3_new = r_model_part.CreateNewGeometry("Triangle2D3", 3, {1,4,3});

        // Check that adding the same type and connectivity geometry with same Id returns the existing one (valid)
        auto p_node_1 = r_model_part.pGetNode(1);
        auto p_node_4 = r_model_part.pGetNode(4);
        auto p_node_3 = r_model_part.pGetNode(3);
        auto p_aux_geom = Kratos::make_shared<Triangle2D3<Node>>(p_node_1, p_node_4, p_node_3);
        auto p_geom_3_new_2 = r_model_part.CreateNewGeometry("Triangle2D3", 3, p_aux_geom);

        // Check results
        KRATOS_EXPECT_EQ(p_geom_3_new, p_geom_3_old);
        KRATOS_EXPECT_EQ(p_geom_3_new_2, p_geom_3_old);
        KRATOS_EXPECT_EQ(model.GetModelPart("TestModelPart").NumberOfGeometries(), 3);
        KRATOS_EXPECT_EQ(model.GetModelPart("TestModelPart.SubModelPart1").NumberOfGeometries(), 1);
        KRATOS_EXPECT_EQ(model.GetModelPart("TestModelPart.SubModelPart2").NumberOfGeometries(), 2);
        KRATOS_EXPECT_EQ(model.GetModelPart("TestModelPart.SubModelPart2.SubSubModelPart").NumberOfGeometries(), 1);
    }

    KRATOS_TEST_CASE_IN_SUITE(TestModelPartGeometryContainerRepeatedGeometries3, KratosCoreGeometryContainerFastSuite)
    {
        // Set up test model part
        Model model;
        SetUpTestModelPart(model);

        // Check that adding a different geometry with the same Id throws an error
        auto& r_model_part = model.GetModelPart("TestModelPart");
        auto p_geom_3_old = r_model_part.pGetGeometry(3);
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(r_model_part.CreateNewGeometry("Quadrilateral2D4", 3, {1,2,4,3}), "Attempting to add geometry with Id: 3. A different geometry with the same Id already exists.")

        // Check that adding a different geometry with the same Id throws an error
        auto p_node_1 = r_model_part.pGetNode(1);
        auto p_node_2 = r_model_part.pGetNode(2);
        auto p_node_3 = r_model_part.pGetNode(3);
        auto p_node_4 = r_model_part.pGetNode(4);
        auto p_aux_geom = Kratos::make_shared<Quadrilateral2D4<Node>>(p_node_1, p_node_2, p_node_4, p_node_3);
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(r_model_part.CreateNewGeometry("Quadrilateral2D4", 3, p_aux_geom), "Attempting to add geometry with Id: 3. A different geometry with the same Id already exists.")

        // Check results
        KRATOS_EXPECT_EQ(model.GetModelPart("TestModelPart").NumberOfGeometries(), 3);
        KRATOS_EXPECT_EQ(model.GetModelPart("TestModelPart.SubModelPart1").NumberOfGeometries(), 1);
        KRATOS_EXPECT_EQ(model.GetModelPart("TestModelPart.SubModelPart2").NumberOfGeometries(), 2);
        KRATOS_EXPECT_EQ(model.GetModelPart("TestModelPart.SubModelPart2.SubSubModelPart").NumberOfGeometries(), 1);
    }

    KRATOS_TEST_CASE_IN_SUITE(TestModelPartGeometryContainerRepeatedGeometries4, KratosCoreGeometryContainerFastSuite)
    {
        // Set up test model part
        Model model;
        SetUpTestModelPart(model);

        // Check that adding the same geometry with the same Id but different connectivities throws an error
        auto& r_model_part = model.GetModelPart("TestModelPart");
        auto p_geom_3_old = r_model_part.pGetGeometry(3);
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(r_model_part.CreateNewGeometry("Triangle2D3", 3, {1,3,4}), "Attempting to add a new geometry with Id: 3. A same type geometry with same Id but different connectivities already exists")

        // Check that adding the same geometry with the same Id but different connectivities throws an error
        auto p_node_1 = r_model_part.pGetNode(1);
        auto p_node_3 = r_model_part.pGetNode(3);
        auto p_node_4 = r_model_part.pGetNode(4);
        auto p_aux_geom = Kratos::make_shared<Triangle2D3<Node>>(p_node_1, p_node_3, p_node_4);
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(r_model_part.CreateNewGeometry("Triangle2D3", 3, p_aux_geom), "Attempting to add a new geometry with Id: 3. A same type geometry with same Id but different connectivities already exists")

        // Check results
        KRATOS_EXPECT_EQ(model.GetModelPart("TestModelPart").NumberOfGeometries(), 3);
        KRATOS_EXPECT_EQ(model.GetModelPart("TestModelPart.SubModelPart1").NumberOfGeometries(), 1);
        KRATOS_EXPECT_EQ(model.GetModelPart("TestModelPart.SubModelPart2").NumberOfGeometries(), 2);
        KRATOS_EXPECT_EQ(model.GetModelPart("TestModelPart.SubModelPart2.SubSubModelPart").NumberOfGeometries(), 1);
    }

    KRATOS_TEST_CASE_IN_SUITE(TestModelPartGeometryContainerRepeatedGeometries5, KratosCoreGeometryContainerFastSuite)
    {
        // Set up test model part
        Model model;
        SetUpTestModelPart(model);

        // Check that creating the same geometry twice returns the already existing one
        auto& r_model_part = model.GetModelPart("TestModelPart");
        auto p_geom_3_string = r_model_part.CreateNewGeometry("Triangle2D3", "3", {1,3,4});
        KRATOS_EXPECT_NE(p_geom_3_string, r_model_part.pGetGeometry(3));

        // Check that adding the same geometry with the same string identifier but different connectivities throws an error
        auto p_node_1 = r_model_part.pGetNode(1);
        auto p_node_3 = r_model_part.pGetNode(3);
        auto p_node_4 = r_model_part.pGetNode(4);
        auto p_aux_geom = Kratos::make_shared<Triangle2D3<Node>>(p_node_1, p_node_4, p_node_3);
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(r_model_part.CreateNewGeometry("Triangle2D3", "3", p_aux_geom), "Attempting to add a new geometry with Id: 3. A same type geometry with same Id but different connectivities already exists")

        // Check results
        KRATOS_EXPECT_EQ(model.GetModelPart("TestModelPart").NumberOfGeometries(), 4);
        KRATOS_EXPECT_EQ(model.GetModelPart("TestModelPart.SubModelPart1").NumberOfGeometries(), 1);
        KRATOS_EXPECT_EQ(model.GetModelPart("TestModelPart.SubModelPart2").NumberOfGeometries(), 2);
        KRATOS_EXPECT_EQ(model.GetModelPart("TestModelPart.SubModelPart2.SubSubModelPart").NumberOfGeometries(), 1);
    }

} // namespace Kratos::Testing.
