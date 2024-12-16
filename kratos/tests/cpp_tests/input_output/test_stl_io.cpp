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
#include <fstream>
#include <sstream>

// External includes

// Project includes
#include "includes/model_part.h"
#include "testing/testing.h"
#include "containers/model.h"
#include "input_output/stl_io.h"

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(ReadTriangleFromSTL, KratosCoreFastSuite)
{
    Kratos::shared_ptr<std::stringstream> p_input = Kratos::make_shared<std::stringstream>(R"input(
    solid 1_triangle
        facet normal  1.000000 0.000000 0.000000 
            outer loop 
            vertex 0.1 -2.56114e-08 0.1
            vertex 0.1 -0.499156 -0.0352136
            vertex 0.1 -0.473406 -0.0446259
            endloop 
        endfacet 
    endsolid 1_triangle
    )input");  

    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");

    StlIO stl_io(p_input);
    stl_io.ReadModelPart(r_model_part);

    KRATOS_EXPECT_TRUE(r_model_part.HasSubModelPart("1_triangle"));
    KRATOS_EXPECT_EQ(r_model_part.GetSubModelPart("1_triangle").NumberOfNodes(), 3);
    KRATOS_EXPECT_EQ(r_model_part.GetSubModelPart("1_triangle").NumberOfGeometries(), 1);
}

KRATOS_TEST_CASE_IN_SUITE(ReadTriangleFromSTLAsElement, KratosCoreFastSuite)
{
    Kratos::shared_ptr<std::stringstream> p_input = Kratos::make_shared<std::stringstream>(R"input(
    solid 1_triangle
        facet normal  1.000000 0.000000 0.000000 
            outer loop 
            vertex 0.1 -2.56114e-08 0.1
            vertex 0.1 -0.499156 -0.0352136
            vertex 0.1 -0.473406 -0.0446259
            endloop 
        endfacet 
    endsolid 1_triangle
    )input");  

    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");

    Parameters settings(R"({
        "new_entity_type" : "element"
    })");
    StlIO stl_io(p_input,settings);
    stl_io.ReadModelPart(r_model_part);

    KRATOS_EXPECT_TRUE(r_model_part.HasSubModelPart("1_triangle"));
    KRATOS_EXPECT_EQ(r_model_part.GetSubModelPart("1_triangle").NumberOfNodes(), 3);
    KRATOS_EXPECT_EQ(r_model_part.GetSubModelPart("1_triangle").NumberOfGeometries(), 0);
    KRATOS_EXPECT_EQ(r_model_part.GetSubModelPart("1_triangle").NumberOfElements(), 1);
}

KRATOS_TEST_CASE_IN_SUITE(ReadSolidNameWithComments, KratosCoreFastSuite)
{

    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");
    Parameters settings(R"({
        "new_entity_type" : "element"
    })");

    Kratos::shared_ptr<std::stringstream> p_input = Kratos::make_shared<std::stringstream>(R"input(
    solid triangle_1 ;Do not read after this
        facet normal  1.000000 0.000000 0.000000
            outer loop
            vertex 0.1 -2.56114e-08 0.1
            vertex 0.1 -0.499156 -0.0352136
            vertex 0.1 -0.473406 -0.0446259
            endloop
        endfacet
    endsolid triangle_1
    )input");

    StlIO stl_io(p_input,settings);
    stl_io.ReadModelPart(r_model_part);
    KRATOS_EXPECT_TRUE(r_model_part.HasSubModelPart("triangle_1"));

    Kratos::shared_ptr<std::stringstream> p_input2 = Kratos::make_shared<std::stringstream>(R"input(
    solid triangle_2    COMMENT: This should be skipped
        facet normal  1.000000 0.000000 0.000000
            outer loop
            vertex 0.1 -2.56114e-08 0.1
            vertex 0.1 -0.499156 -0.0352136
            vertex 0.1 -0.473406 -0.0446259
            endloop
        endfacet
    endsolid triangle_2
    )input");

    StlIO stl_io2(p_input2,settings);
    stl_io2.ReadModelPart(r_model_part);
    KRATOS_EXPECT_TRUE(r_model_part.HasSubModelPart("triangle_2"));

    Kratos::shared_ptr<std::stringstream> p_input3 = Kratos::make_shared<std::stringstream>(R"input(
    solid triangle _ 3    COMMENT: This should be skipped
        facet normal  1.000000 0.000000 0.000000
            outer loop
            vertex 0.1 -2.56114e-08 0.1
            vertex 0.1 -0.499156 -0.0352136
            vertex 0.1 -0.473406 -0.0446259
            endloop
        endfacet
    endsolid triangle _ 3
    )input");

    StlIO stl_io3(p_input3,settings);
    stl_io3.ReadModelPart(r_model_part);
    KRATOS_EXPECT_TRUE(r_model_part.HasSubModelPart("triangle_3"));

    Kratos::shared_ptr<std::stringstream> p_input4 = Kratos::make_shared<std::stringstream>(R"input(
    solid core_851 COMMENT: Exported from Inspire HWVERSION_2023.0.0.21_Jun  6 2023_09:37:24 Build 4500, Altair Inc. Unit : mm
    facet normal  0.999872 0.000003 -0.016028
        outer loop
           vertex 25.4997 58.7925 4.68075e-11
           vertex 25.4961 58.7925 -0.226466
           vertex 25.4959 135 -0.226461
        endloop
    endfacet
    endsolid core_851
    )input");

    StlIO stl_io4(p_input4,settings);
    stl_io4.ReadModelPart(r_model_part);
    KRATOS_EXPECT_TRUE(r_model_part.HasSubModelPart("core_851"));
}

KRATOS_TEST_CASE_IN_SUITE(ReadMultipleTrianglesFromSTL, KratosCoreFastSuite)
{
    Kratos::shared_ptr<std::stringstream> p_input = Kratos::make_shared<std::stringstream>(R"input(
    solid 3_triangles
        facet normal  1.000000 0.000000 0.000000 
            outer loop 
            vertex 0.1 -2.56114e-08 0.1
            vertex 0.1 -0.499156 -0.0352136
            vertex 0.1 -0.473406 -0.0446259
            endloop 
        endfacet 
        facet normal  1.000000 -0.000000 0.000000 
            outer loop 
            vertex 0.1 -0.473406 -0.0446259
            vertex 0.1 -0.447464 -0.0534931
            vertex 0.1 -2.56114e-08 0.1
            endloop 
        endfacet 
        facet normal  1.000000 0.000000 0.000000 
            outer loop 
            vertex 0.1 -0.6 0.1
            vertex 0.1 -0.524702 -0.0252604
            vertex 0.1 -0.499156 -0.0352136
            endloop 
        endfacet 
    endsolid 3_triangles
    )input");  

    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");

    StlIO stl_io(p_input);
    stl_io.ReadModelPart(r_model_part);

    KRATOS_EXPECT_TRUE(r_model_part.HasSubModelPart("3_triangles"));
    KRATOS_EXPECT_EQ(r_model_part.GetSubModelPart("3_triangles").NumberOfNodes(), 9);
    KRATOS_EXPECT_EQ(r_model_part.GetSubModelPart("3_triangles").NumberOfGeometries(), 3);
}

KRATOS_TEST_CASE_IN_SUITE(WriteTriangleToSTL, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");

    // create a few nodes
    r_model_part.CreateNewNode(1, 0.1, -2.56114e-08, 0.1);
    r_model_part.CreateNewNode(2, 0.1, -0.499156, -0.0352136);
    r_model_part.CreateNewNode(3, 0.1, -0.473406, -0.0446259);
    r_model_part.CreateNewNode(4, 0.1, -2.56114e-08, 0.1);

    // create a triangle element
    Properties::Pointer p_properties(new Properties(0));
    r_model_part.CreateNewElement("Element3D3N", 101, {1, 2, 3}, p_properties);
    //this element should be ignored as its area is zero.
    r_model_part.CreateNewElement("Element3D3N", 102, {1, 2, 4}, p_properties);

    Parameters settings(R"({
        "open_mode" : "write"
    })");

    // write stl
    std::filesystem::path filename = "test_stl_write.stl";
    {
        StlIO stl_write(filename, settings);
        stl_write.WriteModelPart(r_model_part);
        // force to write to file
    }
    // read the stl back...
    ModelPart & r_output_model_part = current_model.CreateModelPart("OutputModelPart");
    {
        StlIO stl_read (filename);
        stl_read.ReadModelPart(r_output_model_part);
    }

    // ... and check
    KRATOS_EXPECT_TRUE(r_output_model_part.HasSubModelPart("Main"));
    KRATOS_EXPECT_EQ(r_output_model_part.NumberOfProperties(), 1);
    KRATOS_EXPECT_EQ(r_output_model_part.NumberOfSubModelParts() ,1);
    KRATOS_EXPECT_EQ(r_output_model_part.NumberOfNodes(), 3);
    KRATOS_EXPECT_EQ(r_output_model_part.NumberOfGeometries(), 1);
    KRATOS_EXPECT_EQ(r_output_model_part.NumberOfConditions(), 0);

    KRATOS_EXPECT_EQ(r_output_model_part.GetSubModelPart("Main").NumberOfNodes(), 3);
    KRATOS_EXPECT_EQ(r_output_model_part.GetSubModelPart("Main").NumberOfGeometries(), 1);
    KRATOS_EXPECT_EQ(r_output_model_part.GetSubModelPart("Main").NumberOfElements(), 0);
    KRATOS_EXPECT_EQ(r_output_model_part.GetSubModelPart("Main").NumberOfConditions(), 0);

    // remove the generated files
    if (std::filesystem::remove(filename) != true) {
        KRATOS_ERROR << "Error deleting test output file: " << filename << "\n";
    }
}

}  // namespace Kratos::Testing.
