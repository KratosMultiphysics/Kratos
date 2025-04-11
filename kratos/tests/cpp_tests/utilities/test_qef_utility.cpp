//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Ariadna Cortes
//
//
// System includes

// External includes

// Project includes
#include "includes/expect.h"
#include "testing/testing.h"
#include "containers/model.h"
#include "utilities/qef_utility.h"

namespace Kratos::Testing {
namespace {
    using NodeType = Node;
    using NodePtrType = Node::Pointer;
    using GeometryType = Geometry<NodeType>;
    using GeometryPtrType = GeometryType::Pointer;
    using GeometryArrayType = GeometryType::GeometriesArrayType;
    using PointsArrayType = GeometryType::PointsArrayType;

    /*this functions are creating a model for each element used in the test in order to avoid trouble with repeated Ids.
    This is not optimal but since it is just a test (in real cases, funcions will be used with well-defined models)...
    */
    GeometryPtrType QuadraticErrorFunctionGenerateHexahedra3D8(const std::vector<double>& rDistances)
    {
        Model my_model;
        ModelPart& r_voxel = my_model.CreateModelPart("Voxel");
        r_voxel.AddNodalSolutionStepVariable(DISTANCE);

        r_voxel.CreateNewNode(1, -1, -1, -1);
        r_voxel.CreateNewNode(2,  1, -1, -1);
        r_voxel.CreateNewNode(3, 1,  1, -1);
        r_voxel.CreateNewNode(4, -1,  1, -1);
        r_voxel.CreateNewNode(5, -1, -1,  1);
        r_voxel.CreateNewNode(6, 1, -1,  1);
        r_voxel.CreateNewNode(7, 1,  1,  1);
        r_voxel.CreateNewNode(8, -1,  1,  1);
        Properties::Pointer p_properties_0(new Properties(0));
        Element::Pointer pElement = r_voxel.CreateNewElement("Element3D8N", 1, {1, 2, 3, 4, 5, 6, 7, 8}, p_properties_0);  
        GeometryPtrType p_voxel = pElement->pGetGeometry();

        //Add the distances
        PointsArrayType nodes = p_voxel->Points();

        for (unsigned int i = 0; i < 8; i++) {
            nodes[i].FastGetSolutionStepValue(DISTANCE) = rDistances[i];
        }
        return p_voxel;
    }

    GeometryPtrType QuadraticErrorFunctionGenerateUncenteredHexahedra3D8(const std::vector<double>& rDistances) 
    {
        Model my_model;
        ModelPart& r_voxel = my_model.CreateModelPart("Voxel");
        r_voxel.AddNodalSolutionStepVariable(DISTANCE);

        r_voxel.CreateNewNode(1, 0, 0, 0);
        r_voxel.CreateNewNode(2,  2, 0, 0);
        r_voxel.CreateNewNode(3, 2,  2, 0);
        r_voxel.CreateNewNode(4, 0,  2, 0);
        r_voxel.CreateNewNode(5, 0, 0, 2);
        r_voxel.CreateNewNode(6, 2, 0,  2);
        r_voxel.CreateNewNode(7, 2,  2,  2);
        r_voxel.CreateNewNode(8, 0,  2,  2);
        Properties::Pointer p_properties_0(new Properties(0));
        Element::Pointer pElement = r_voxel.CreateNewElement("Element3D8N", 1, {1, 2, 3, 4, 5, 6, 7, 8}, p_properties_0);  
        GeometryPtrType p_voxel = pElement->pGetGeometry();

        // Add the distances
        PointsArrayType nodes = p_voxel->Points();

        for (unsigned int i = 0; i < 8; i++) {
            nodes[i].FastGetSolutionStepValue(DISTANCE) = rDistances[i];
        }
        return p_voxel;
    }

    //rNodes is a 3*3 matrix representin de (x,y,z) coordinates of each of the 3 triangle nodes
    GeometryPtrType QuadraticErrorFunctionGenerateTriangle3D3(std::vector<std::vector<double>>& rNodes)
    {
        Model my_model;
        ModelPart& r_triangles = my_model.CreateModelPart("Triangle");
        r_triangles.CreateNewNode(1, rNodes[0][0], rNodes[0][1], rNodes[0][2]);
        r_triangles.CreateNewNode(2, rNodes[1][0], rNodes[1][1], rNodes[1][2]);
        r_triangles.CreateNewNode(3, rNodes[2][0], rNodes[2][1], rNodes[2][2]);
        Properties::Pointer p_properties_1(new Properties(0));
        Element::Pointer p_triangle = r_triangles.CreateNewElement("Element3D3N", 1, {1, 2, 3}, p_properties_1);
        return p_triangle->pGetGeometry();
    }
} //unnamed namespace

KRATOS_TEST_CASE_IN_SUITE(QuadraticErrorFunction0dof, KratosCoreFastSuite)
{
    // A voxel crossed by a straight plane with only 1 nodes inside the volume
    std::vector<double> distances{1, -1, -1, -1, -1, -1, -1, -1};
    GeometryPtrType p_voxel = QuadraticErrorFunctionGenerateHexahedra3D8(distances);

    // Generate the intersecting triangles
    std::vector<std::vector<double>> triangle1{{-1,0.5,-0.95},{-0.95,0.5,-1.05},{-1.05,0.5,-1.05}};
    std::vector<std::vector<double>> triangle2{{-1,-0.95,0},{-0.95,-1.05,0},{-1.05,-1.05,0}};
    std::vector<std::vector<double>> triangle3{{0.75,-0.95,-1},{0.75,-1.05,-0.95},{0.75,-1.05,-1.05}};

    GeometryPtrType p_triangle1 = QuadraticErrorFunctionGenerateTriangle3D3(triangle1);
    GeometryPtrType p_triangle2 = QuadraticErrorFunctionGenerateTriangle3D3(triangle2);
    GeometryPtrType p_triangle3 = QuadraticErrorFunctionGenerateTriangle3D3(triangle3);

    GeometryArrayType array1;
    array1.push_back(p_triangle1);
    array1.push_back(p_triangle2);
    array1.push_back(p_triangle3);

    const array_1d<double,3> point = QuadraticErrorFunction::QuadraticErrorFunctionPoint(*p_voxel,array1);

    KRATOS_EXPECT_NEAR(point[0], 0.75, 1e-8);
    KRATOS_EXPECT_NEAR(point[1], 0.5, 1e-8);
    KRATOS_EXPECT_NEAR(point[2], 0.0, 1e-8);
}

KRATOS_TEST_CASE_IN_SUITE(QuadraticErrorFunction1dof, KratosCoreFastSuite)
{
    //A voxel crossed by a straight plane with only 2 nodes inside the volume
    std::vector<double> distances{1, 1, -1, -1, -1, -1, -1, -1};
    GeometryPtrType p_voxel = QuadraticErrorFunctionGenerateHexahedra3D8(distances);

    //Generate the intersecting triangles
    std::vector<std::vector<double>> triangle1{{-1,0.5,-0.95},{-0.95,0.5,-1.05},{-1.05,0.5,-1.05}};
    std::vector<std::vector<double>> triangle2{{-1,-0.95,0},{-0.95,-1.05,0},{-1.05,-1.05,0}};
    std::vector<std::vector<double>> triangle3{{1,-0.95,0},{1.05,-1.05,0},{0.95,-1.05,0}};
    std::vector<std::vector<double>> triangle4{{1,0.5,-0.95},{1.05,0.5,-1.05},{0.95,0.5,-1.05}};

    GeometryPtrType p_triangle1 = QuadraticErrorFunctionGenerateTriangle3D3(triangle1);
    GeometryPtrType p_triangle2 = QuadraticErrorFunctionGenerateTriangle3D3(triangle2);
    GeometryPtrType p_triangle3 = QuadraticErrorFunctionGenerateTriangle3D3(triangle3);
    GeometryPtrType p_triangle4 = QuadraticErrorFunctionGenerateTriangle3D3(triangle4);

    GeometryArrayType array1;
    array1.push_back(p_triangle1);
    array1.push_back(p_triangle2);
    array1.push_back(p_triangle3);
    array1.push_back(p_triangle4);

    const array_1d<double,3> point = QuadraticErrorFunction::QuadraticErrorFunctionPoint(*p_voxel,array1);

    KRATOS_EXPECT_NEAR(point[0], 0.0, 1e-8);
    KRATOS_EXPECT_NEAR(point[1], 0.5, 1e-8);
    KRATOS_EXPECT_NEAR(point[2], 0.0, 1e-8);
}

KRATOS_TEST_CASE_IN_SUITE(QuadraticErrorFunction0dofExtremeCase, KratosCoreFastSuite)
{
    // A voxel crossed by a straight plane with only 1 nodes inside the volume
    std::vector<double> distances{1, -1, -1, -1, -1, -1, -1, -1};
    GeometryPtrType p_voxel = QuadraticErrorFunctionGenerateHexahedra3D8(distances);

    // Generate the intersecting triangles
    std::vector<std::vector<double>> triangle1{{-1,0.5,-0.95},{-0.95,0.5,-1.05},{-1.05,0.5,-1.05}};
    std::vector<std::vector<double>> triangle2{{-1,-0.95,0},{-0.95,-1.05,0},{-1.05,-1.05,0}};
    std::vector<std::vector<double>> triangle3{{0.9999,-0.95,-1},{0.9999,-1.05,-0.95},{0.9999,-1.05,-1.05}}; 

    GeometryPtrType p_triangle1 = QuadraticErrorFunctionGenerateTriangle3D3(triangle1);
    GeometryPtrType p_triangle2 = QuadraticErrorFunctionGenerateTriangle3D3(triangle2);
    GeometryPtrType p_triangle3 = QuadraticErrorFunctionGenerateTriangle3D3(triangle3);

    GeometryArrayType array1;
    array1.push_back(p_triangle1); 
    array1.push_back(p_triangle2);
    array1.push_back(p_triangle3);

    const array_1d<double,3> point = QuadraticErrorFunction::QuadraticErrorFunctionPoint(*p_voxel,array1);

    KRATOS_EXPECT_NEAR(point[0], 0.9999, 1e-8);
    KRATOS_EXPECT_NEAR(point[1], 0.5, 1e-8);
    KRATOS_EXPECT_NEAR(point[2], 0.0, 1e-8);
}

KRATOS_TEST_CASE_IN_SUITE(QuadraticErrorFunction1dofExtremeCase, KratosCoreFastSuite)
{
    // A voxel crossed by a straight plane with only 2 nodes inside the volume (not good case approximation)
    std::vector<double> distances{1, 1, -1, -1, -1, -1, -1, -1};
    GeometryPtrType p_voxel = QuadraticErrorFunctionGenerateHexahedra3D8(distances);

    // Generate the intersecting triangles
    std::vector<std::vector<double>> triangle1{{-1,0.5,-0.95},{-0.95,0.5,-1.05},{-1.05,0.5,-1.05}};
    std::vector<std::vector<double>> triangle2{{-1,-0.95,0.9999},{-0.95,-1.05,0.9999},{-1.05,-1.05,0.9999}};
    std::vector<std::vector<double>> triangle3{{1,-0.95,0.9999},{1.05,-1.05,0.9999},{0.95,-1.05,0.9999}};
    std::vector<std::vector<double>> triangle4{{1,0.5,-0.95},{1.05,0.5,-1.05},{0.95,0.5,-1.05}};

    GeometryPtrType p_triangle1 = QuadraticErrorFunctionGenerateTriangle3D3(triangle1);
    GeometryPtrType p_triangle2 = QuadraticErrorFunctionGenerateTriangle3D3(triangle2);
    GeometryPtrType p_triangle3 = QuadraticErrorFunctionGenerateTriangle3D3(triangle3);
    GeometryPtrType p_triangle4 = QuadraticErrorFunctionGenerateTriangle3D3(triangle4);

    GeometryArrayType array1;
    array1.push_back(p_triangle1);
    array1.push_back(p_triangle2);
    array1.push_back(p_triangle3);
    array1.push_back(p_triangle4);

    const array_1d<double,3> point = QuadraticErrorFunction::QuadraticErrorFunctionPoint(*p_voxel,array1);

    KRATOS_EXPECT_NEAR(point[0], 0.0, 1e-8);
    KRATOS_EXPECT_NEAR(point[1], 0.5, 1e-8);
    KRATOS_EXPECT_NEAR(point[2], 0.9999, 1e-8);
}

KRATOS_TEST_CASE_IN_SUITE(QuadraticErrorFunction2dof, KratosCoreFastSuite)
{
    // A voxel crossed by a straight plane with 4 nodes inside the volume
    std::vector<double> distances{1, 1, -1, -1, 1, 1, -1, -1};
    GeometryPtrType p_voxel = QuadraticErrorFunctionGenerateHexahedra3D8(distances);

    // Generate the intersecting triangles
    std::vector<std::vector<double>> triangle1{{1,0.5,1.05},{1.05,0.5,0.95},{0.95,0.5,0.95}};
    std::vector<std::vector<double>> triangle2{{1,0.5,-0.95},{1.05,0.5,-1.05},{0.95,0.5,-1.05}};
    std::vector<std::vector<double>> triangle3{{-1,0.5,-0.95},{-0.95,0.5,-1.05},{-1.05,0.5,-1.05}};
    std::vector<std::vector<double>> triangle4{{-1,0.5,1.05},{-0.95,0.5,0.95},{-1.05,0.5,0.95}};

    GeometryPtrType p_triangle1 = QuadraticErrorFunctionGenerateTriangle3D3(triangle1);
    GeometryPtrType p_triangle2 = QuadraticErrorFunctionGenerateTriangle3D3(triangle2);
    GeometryPtrType p_triangle3 = QuadraticErrorFunctionGenerateTriangle3D3(triangle3);
    GeometryPtrType p_triangle4 = QuadraticErrorFunctionGenerateTriangle3D3(triangle4);

    GeometryArrayType array1;
    array1.push_back(p_triangle1);
    array1.push_back(p_triangle2);
    array1.push_back(p_triangle3);
    array1.push_back(p_triangle4);

    const array_1d<double,3> point = QuadraticErrorFunction::QuadraticErrorFunctionPoint(*p_voxel,array1);

    KRATOS_EXPECT_NEAR(point[0], 0.0, 1e-8);
    KRATOS_EXPECT_NEAR(point[1], 0.5, 1e-8);
    KRATOS_EXPECT_NEAR(point[2], 0.0, 1e-8);
}

KRATOS_TEST_CASE_IN_SUITE(QuadraticErrorFunctiondoublePlain, KratosCoreFastSuite)
{
    // A voxel crossed by two straight planes (enclosing half of the volume) with no nodes inside the volume
    std::vector<double> distances{-1, -1, -1, -1, -1, -1, -1, -1};
    GeometryPtrType p_voxel = QuadraticErrorFunctionGenerateHexahedra3D8(distances);

    // Generate the intersecting triangles
    std::vector<std::vector<double>> triangle1{{1,0.5,1.05},{1.05,0.5,0.95},{0.95,0.5,0.95}};
    std::vector<std::vector<double>> triangle2{{1,0.5,-0.95},{1.05,0.5,-1.05},{0.95,0.5,-1.05}};
    std::vector<std::vector<double>> triangle3{{-1,0.5,-0.95},{-0.95,0.5,-1.05},{-1.05,0.5,-1.05}};
    std::vector<std::vector<double>> triangle4{{-1,0.5,1.05},{-0.95,0.5,0.95},{-1.05,0.5,0.95}};
    std::vector<std::vector<double>> triangle5{{1,-0.5,1.05},{1.05,-0.5,0.95},{0.95,-0.5,0.95}};
    std::vector<std::vector<double>> triangle6{{1,-0.5,-0.95},{1.05,-0.5,-1.05},{0.95,-0.5,-1.05}};
    std::vector<std::vector<double>> triangle7{{-1,-0.5,-0.95},{-0.95,-0.5,-1.05},{-1.05,-0.5,-1.05}};
    std::vector<std::vector<double>> triangle8{{-1,-0.5,1.05},{-0.95,-0.5,0.95},{-1.05,-0.5,0.95}};

    GeometryPtrType p_triangle1 = QuadraticErrorFunctionGenerateTriangle3D3(triangle1);
    GeometryPtrType p_triangle2 = QuadraticErrorFunctionGenerateTriangle3D3(triangle2);
    GeometryPtrType p_triangle3 = QuadraticErrorFunctionGenerateTriangle3D3(triangle3);
    GeometryPtrType p_triangle4 = QuadraticErrorFunctionGenerateTriangle3D3(triangle4);
    GeometryPtrType p_triangle5 = QuadraticErrorFunctionGenerateTriangle3D3(triangle5);
    GeometryPtrType p_triangle6 = QuadraticErrorFunctionGenerateTriangle3D3(triangle6);
    GeometryPtrType p_triangle7 = QuadraticErrorFunctionGenerateTriangle3D3(triangle7);
    GeometryPtrType p_triangle8 = QuadraticErrorFunctionGenerateTriangle3D3(triangle8);

    GeometryArrayType array1;
    array1.push_back(p_triangle1);
    array1.push_back(p_triangle2);
    array1.push_back(p_triangle3);
    array1.push_back(p_triangle4);
    array1.push_back(p_triangle5);
    array1.push_back(p_triangle6);
    array1.push_back(p_triangle7);
    array1.push_back(p_triangle8);

    // Call the point utility
    const array_1d<double,3> point = QuadraticErrorFunction::QuadraticErrorFunctionPoint(*p_voxel,array1);

    KRATOS_EXPECT_NEAR(point[0], 0.0, 1e-8);
    KRATOS_EXPECT_NEAR(point[1], 0.0, 1e-8);
    KRATOS_EXPECT_NEAR(point[2], 0.0, 1e-8);
}

KRATOS_TEST_CASE_IN_SUITE(QuadraticErrorFunction0dofmovedCenter, KratosCoreFastSuite)
{
    // A voxel with one node inside the volume but centered at (1,1,1)
    std::vector<double> distances{1, -1, -1, -1, -1, -1, -1, -1};
    GeometryPtrType p_voxel = QuadraticErrorFunctionGenerateUncenteredHexahedra3D8(distances);

    std::vector<std::vector<double>> triangle1{{0,1.5,0.05},{0.05,1.5,-0.05},{-0.05,1.5,-0.05}};
    std::vector<std::vector<double>> triangle2{{0,0.05,1},{0.05,-0.05,1},{-0.05,-0.05,1}};
    std::vector<std::vector<double>> triangle3{{1.75,0.05,0},{1.75,-0.05,0.05},{1.75,-0.05,-0.05}};

    GeometryPtrType p_triangle1 = QuadraticErrorFunctionGenerateTriangle3D3(triangle1);
    GeometryPtrType p_triangle2 = QuadraticErrorFunctionGenerateTriangle3D3(triangle2);
    GeometryPtrType p_triangle3 = QuadraticErrorFunctionGenerateTriangle3D3(triangle3);

    GeometryArrayType array1;
    array1.push_back(p_triangle1);
    array1.push_back(p_triangle2);
    array1.push_back(p_triangle3);

    const array_1d<double,3> point = QuadraticErrorFunction::QuadraticErrorFunctionPoint(*p_voxel,array1);

    KRATOS_EXPECT_NEAR(point[0], 1.75, 1e-8);
    KRATOS_EXPECT_NEAR(point[1], 1.5, 1e-8);
    KRATOS_EXPECT_NEAR(point[2], 1.0, 1e-8);
}

KRATOS_TEST_CASE_IN_SUITE(QuadraticErrorFunction2dofMovedCenter, KratosCoreFastSuite)
{
    std::vector<double> distances{1, 1, -1, -1, 1, 1, -1, -1};
    GeometryPtrType p_voxel = QuadraticErrorFunctionGenerateUncenteredHexahedra3D8(distances);

    std::vector<std::vector<double>> triangle1{{2,1.5,2.05},{2.05,1.5,1.95},{1.95,1.5,1.95}};
    std::vector<std::vector<double>> triangle2{{2,1.5,0.05},{2.05,1.5,-0.05},{1.95,1.5,-0.05}};
    std::vector<std::vector<double>> triangle3{{0,1.5,0.05},{0.05,1.5,-0.05},{-0.05,1.5,-0.05}};
    std::vector<std::vector<double>> triangle4{{0,1.5,2.05},{0.05,1.5,1.95},{-0.05,1.5,1.95}};

    GeometryPtrType p_triangle1 = QuadraticErrorFunctionGenerateTriangle3D3(triangle1);
    GeometryPtrType p_triangle2 = QuadraticErrorFunctionGenerateTriangle3D3(triangle2);
    GeometryPtrType p_triangle3 = QuadraticErrorFunctionGenerateTriangle3D3(triangle3);
    GeometryPtrType p_triangle4 = QuadraticErrorFunctionGenerateTriangle3D3(triangle4);

    GeometryArrayType array1;
    array1.push_back(p_triangle1);
    array1.push_back(p_triangle2);
    array1.push_back(p_triangle3);
    array1.push_back(p_triangle4);

    const array_1d<double,3> point = QuadraticErrorFunction::QuadraticErrorFunctionPoint(*p_voxel,array1);

    KRATOS_EXPECT_NEAR(point[0], 1.0, 1e-8);
    KRATOS_EXPECT_NEAR(point[1], 1.5, 1e-8);
    KRATOS_EXPECT_NEAR(point[2], 1.0, 1e-8);
}

}  // namespace Kratos::Testing.