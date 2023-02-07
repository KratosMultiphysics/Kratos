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

KRATOS_TEST_CASE_IN_SUITE(ReadTriangleFromSTL, KratosCoreFastSuite)
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

    KRATOS_CHECK(r_model_part.HasSubModelPart("1 triangle"));
    KRATOS_CHECK_EQUAL(r_model_part.GetSubModelPart("1 triangle").NumberOfNodes(), 3);
    KRATOS_CHECK_EQUAL(r_model_part.GetSubModelPart("1 triangle").NumberOfElements(), 1);
}

KRATOS_TEST_CASE_IN_SUITE(ReadMultipleTrianglesFromSTL, KratosCoreFastSuite)
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

    KRATOS_CHECK(r_model_part.HasSubModelPart("3 triangles"));
    KRATOS_CHECK_EQUAL(r_model_part.GetSubModelPart("3 triangles").NumberOfNodes(), 9);
    KRATOS_CHECK_EQUAL(r_model_part.GetSubModelPart("3 triangles").NumberOfElements(), 3);
}

KRATOS_TEST_CASE_IN_SUITE(WriteTriangleToSTL, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");

    // create a few nodes
    r_model_part.CreateNewNode(1, 0.1, -2.56114e-08, 0.1);
    r_model_part.CreateNewNode(2, 0.1, -0.499156, -0.0352136);
    r_model_part.CreateNewNode(3, 0.1, -0.473406, -0.0446259);

    // create a triangle element
    Properties::Pointer p_properties(new Properties(0));
    r_model_part.CreateNewElement("Element3D3N", 101, {1, 2, 3}, p_properties);

    // 
    StlIO stl_io ("test_stl_write.stl");
    stl_io.WriteModelPart(r_model_part);

    // intentional failure
    KRATOS_CHECK_EQUAL(true, false);
}

}  // namespace Kratos::Testing.
