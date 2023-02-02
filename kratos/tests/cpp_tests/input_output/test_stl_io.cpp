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
#include <sstream>

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "input_output/stl_io.h"

namespace Kratos::Testing {

TEST(ReadTriangleFromSTL, KratosCoreFastSuite)
{
    std::stringstream* p_input = new std::stringstream(R"input(
    solid 1 triangle
        facet normal  1.000000 0.000000 0.000000 
            outer loop 
            vertex 0.1 -2.56114e-08 0.1
            vertex 0.1 -0.499156 -0.0352136
            vertex 0.1 -0.473406 -0.0446259
            endloop 
        endfacet 
    endsolid 1 triangle
    )input");  

    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");

    StlIO stl_io(p_input);
    stl_io.ReadModelPart(r_model_part);

    std::cout << r_model_part << std::endl;

    KRATOS_EXPECT_TRUE(r_model_part.HasSubModelPart("1 triangle"));
    KRATOS_EXPECT_EQ(r_model_part.GetSubModelPart("1 triangle").NumberOfNodes(), 3);
    KRATOS_EXPECT_EQ(r_model_part.GetSubModelPart("1 triangle").NumberOfElements(), 1);
}

TEST(ReadMultipleTrianglesFromSTL, KratosCoreFastSuite)
{
    std::stringstream* p_input = new std::stringstream(R"input(
    solid 3 triangles
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
    endsolid 3 triangles
    )input");  

    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");

    StlIO stl_io(p_input);
    stl_io.ReadModelPart(r_model_part);

    std::cout << r_model_part << std::endl;

    KRATOS_EXPECT_TRUE(r_model_part.HasSubModelPart("3 triangles"));
    KRATOS_EXPECT_EQ(r_model_part.GetSubModelPart("3 triangles").NumberOfNodes(), 9);
    KRATOS_EXPECT_EQ(r_model_part.GetSubModelPart("3 triangles").NumberOfElements(), 3);
}

}  // namespace Kratos::Testing.
