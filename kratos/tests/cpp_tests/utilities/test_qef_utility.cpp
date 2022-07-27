//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/license.txt
//
//  Main authors:    Ariadna Cortés
//
//

// Project includes
#include "includes/checks.h"
#include "testing/testing.h"
#include "utilities/qef_utility.h"

namespace Kratos {
namespace Testing {

    typedef Node<3> NodeType;
    typedef Node<3>::Pointer NodePtrType;
    typedef Geometry<NodeType> GeometryType;
    typedef GeometryType::Pointer GeometryPtrType;
    typedef GeometryType::GeometriesArrayType GeometryArrayType;
    typedef GeometryType::PointsArrayType PointsArrayType;

    /*this functions are creating a model for each element used in the test in order to avoid trouble with repeated Ids.
    This is not optimal but since it is just a test (in real cases, funcions will be used with well-defined models)...
    */
    GeometryPtrType GenerateHexahedra3D8(const std::vector<double>& rDistances) 
    { 
        Model my_model;
        ModelPart &voxel = my_model.CreateModelPart("Voxel");  
        voxel.AddNodalSolutionStepVariable(DISTANCE);  

        voxel.CreateNewNode(1, -1, -1, -1);
        voxel.CreateNewNode(2,  1, -1, -1);
        voxel.CreateNewNode(3, 1,  1, -1);
        voxel.CreateNewNode(4, -1,  1, -1);
        voxel.CreateNewNode(5, -1, -1,  1);
        voxel.CreateNewNode(6, 1, -1,  1);
        voxel.CreateNewNode(7, 1,  1,  1);
        voxel.CreateNewNode(8, -1,  1,  1); 
        Properties::Pointer p_properties_0(new Properties(0));
        Element::Pointer pElement = voxel.CreateNewElement("Element3D8N", 1, {1, 2, 3, 4, 5, 6, 7, 8}, p_properties_0);  
        GeometryPtrType pVoxel = pElement->pGetGeometry();

        //Add the distances
        PointsArrayType nodes = pVoxel->Points();   
        
        for (int i = 0; i < 8; i++) {
            nodes[i].FastGetSolutionStepValue(DISTANCE) = rDistances[i];
        }
        return pVoxel;  
    }
    GeometryPtrType GenerateUncenteredHexahedra3D8(const std::vector<double>& rDistances) 
    { 
        Model my_model;
        ModelPart &voxel = my_model.CreateModelPart("Voxel");  
        voxel.AddNodalSolutionStepVariable(DISTANCE);  

        voxel.CreateNewNode(1, 0, 0, 0);
        voxel.CreateNewNode(2,  2, 0, 0);
        voxel.CreateNewNode(3, 2,  2, 0);
        voxel.CreateNewNode(4, 0,  2, 0);
        voxel.CreateNewNode(5, 0, 0, 2);
        voxel.CreateNewNode(6, 2, 0,  2);
        voxel.CreateNewNode(7, 2,  2,  2);
        voxel.CreateNewNode(8, 0,  2,  2); 
        Properties::Pointer p_properties_0(new Properties(0));
        Element::Pointer pElement = voxel.CreateNewElement("Element3D8N", 1, {1, 2, 3, 4, 5, 6, 7, 8}, p_properties_0);  
        GeometryPtrType pVoxel = pElement->pGetGeometry();

        //Add the distances
        PointsArrayType nodes = pVoxel->Points();   
        
        for (int i = 0; i < 8; i++) {
            nodes[i].FastGetSolutionStepValue(DISTANCE) = rDistances[i];
        }
        return pVoxel;  
    }

    //rNodes is a 3*3 matrix representin de (x,y,z) coordinates of each of the 3 triangle nodes
    GeometryPtrType GenerateTriangle3D3(std::vector<std::vector<double>>& rNodes)
    {
        Model my_model;
        ModelPart& triangles = my_model.CreateModelPart("Triangle"); 
        triangles.CreateNewNode(1, rNodes[0][0], rNodes[0][1], rNodes[0][2]);
        triangles.CreateNewNode(2, rNodes[1][0], rNodes[1][1], rNodes[1][2]);
        triangles.CreateNewNode(3, rNodes[2][0], rNodes[2][1], rNodes[2][2]);
        Properties::Pointer p_properties_1(new Properties(0)); 
        Element::Pointer pTriangle = triangles.CreateNewElement("Element3D3N", 1, {1, 2, 3}, p_properties_1);
        return pTriangle->pGetGeometry();
    } 

    KRATOS_TEST_CASE_IN_SUITE(CalculateCenter, KratosCoreFastSuite) 
    {

        //Generate the HEXAHEDRA3D8
        std::vector<double> distances{1, -1, -1, -1, -1, -1, -1, -1,};   
        GeometryPtrType pVoxel = GenerateHexahedra3D8(distances);

        array_1d<double,3> center = QEF::CalculateCenter(*pVoxel);

        //Expected output of the function
        KRATOS_CHECK_EQUAL(center[0], 0.0);
        KRATOS_CHECK_EQUAL(center[1], 0.0);
        KRATOS_CHECK_EQUAL(center[2], 0.0);
    }

    KRATOS_TEST_CASE_IN_SUITE(CalculateNormal, KratosCoreFastSuite) 
    {

        //Generate a triangles
        std::vector<std::vector<double>> triangle1{{0,0,0},{1,0,0},{0,1,0}};    
        GeometryPtrType pTriangle = GenerateTriangle3D3(triangle1);

        //Call the normal utility
        array_1d<double,3> normal = QEF::CalculateNormal(*pTriangle);

        //Expected output of the function
        KRATOS_CHECK_EQUAL(normal[0], 0.0);
        KRATOS_CHECK_EQUAL(normal[1], 0.0);
        KRATOS_CHECK_EQUAL(normal[2], 1.0);
    }

    KRATOS_TEST_CASE_IN_SUITE(QEF0dof, KratosCoreFastSuite) {
        //A voxel crossed by a straight plane with only 2 nodes inside the volume (not good case approximation)
        std::vector<double> distances{1, -1, -1, -1, -1, -1, -1, -1};   
        GeometryPtrType pVoxel = GenerateHexahedra3D8(distances);

        //Generate the intersecting triangles
        std::vector<std::vector<double>> triangle1{{-1,0.5,-0.95},{-0.95,0.5,-1.05},{-1.05,0.5,-1.05}}; 
        std::vector<std::vector<double>> triangle2{{-1,-0.95,0},{-0.95,-1.05,0},{-1.05,-1.05,0}}; 
        std::vector<std::vector<double>> triangle3{{0.75,-0.95,-1},{0.75,-1.05,-0.95},{0.75,-1.05,-1.05}}; 

        
        
        GeometryPtrType pTriangle1 = GenerateTriangle3D3(triangle1);
        GeometryPtrType pTriangle2 = GenerateTriangle3D3(triangle2);
        GeometryPtrType pTriangle3 = GenerateTriangle3D3(triangle3);

        GeometryArrayType array1;
        array1.push_back(pTriangle1); 
        array1.push_back(pTriangle2);
        array1.push_back(pTriangle3);

        array_1d<double,3> point = QEF::QEF_point(*pVoxel,array1);

        KRATOS_CHECK_NEAR(point[0], 0.75, 1e-8);
        KRATOS_CHECK_NEAR(point[1], 0.5, 1e-8);
        KRATOS_CHECK_NEAR(point[2], 0.0, 1e-8);        
    }

     KRATOS_TEST_CASE_IN_SUITE(QEF1dof, KratosCoreFastSuite) {
        //A voxel crossed by a straight plane with only 2 nodes inside the volume (not good case approximation)
        std::vector<double> distances{1, 1, -1, -1, -1, -1, -1, -1};   
        GeometryPtrType pVoxel = GenerateHexahedra3D8(distances);

        //Generate the intersecting triangles
        std::vector<std::vector<double>> triangle1{{-1,0.5,-0.95},{-0.95,0.5,-1.05},{-1.05,0.5,-1.05}}; 
        std::vector<std::vector<double>> triangle2{{-1,-0.95,0},{-0.95,-1.05,0},{-1.05,-1.05,0}}; 
        std::vector<std::vector<double>> triangle3{{1,-0.95,0},{1.05,-1.05,0},{0.95,-1.05,0}}; 
        std::vector<std::vector<double>> triangle4{{1,0.5,-0.95},{1.05,0.5,-1.05},{0.95,0.5,-1.05}}; 
        
        GeometryPtrType pTriangle1 = GenerateTriangle3D3(triangle1);
        GeometryPtrType pTriangle2 = GenerateTriangle3D3(triangle2);
        GeometryPtrType pTriangle3 = GenerateTriangle3D3(triangle3);
        GeometryPtrType pTriangle4 = GenerateTriangle3D3(triangle4);

        GeometryArrayType array1;
        array1.push_back(pTriangle1); 
        array1.push_back(pTriangle2);
        array1.push_back(pTriangle3);
        array1.push_back(pTriangle4);

        array_1d<double,3> point = QEF::QEF_point(*pVoxel,array1);

        KRATOS_CHECK_NEAR(point[0], 0.0, 1e-8);
        KRATOS_CHECK_NEAR(point[1], 0.5, 1e-8);
        KRATOS_CHECK_NEAR(point[2], 0.0, 1e-8);        
    }

    KRATOS_TEST_CASE_IN_SUITE(QEF2dof, KratosCoreFastSuite) {
        //A voxel crossed by a straight plane with 4 nodes inside the volume
        std::vector<double> distances{1, 1, -1, -1, 1, 1, -1, -1};   
        GeometryPtrType pVoxel = GenerateHexahedra3D8(distances);

        //Generate the intersecting triangles
        std::vector<std::vector<double>> triangle1{{1,0.5,1.05},{1.05,0.5,0.95},{0.95,0.5,0.95}};  
        std::vector<std::vector<double>> triangle2{{1,0.5,-0.95},{1.05,0.5,-1.05},{0.95,0.5,-1.05}}; 
        std::vector<std::vector<double>> triangle3{{-1,0.5,-0.95},{-0.95,0.5,-1.05},{-1.05,0.5,-1.05}}; 
        std::vector<std::vector<double>> triangle4{{-1,0.5,1.05},{-0.95,0.5,0.95},{-1.05,0.5,0.95}}; 
        
        GeometryPtrType pTriangle1 = GenerateTriangle3D3(triangle1);
        GeometryPtrType pTriangle2 = GenerateTriangle3D3(triangle2);
        GeometryPtrType pTriangle3 = GenerateTriangle3D3(triangle3);
        GeometryPtrType pTriangle4 = GenerateTriangle3D3(triangle4);

        GeometryArrayType array1;
        array1.push_back(pTriangle1); 
        array1.push_back(pTriangle2);
        array1.push_back(pTriangle3);
        array1.push_back(pTriangle4);

        array_1d<double,3> point = QEF::QEF_point(*pVoxel,array1);

        KRATOS_CHECK_NEAR(point[0], 0.0, 1e-8);
        KRATOS_CHECK_NEAR(point[1], 0.5, 1e-8);
        KRATOS_CHECK_NEAR(point[2], 0.0, 1e-8);        
    } 

    KRATOS_TEST_CASE_IN_SUITE(QEFdoublePlain, KratosCoreFastSuite) {
        //A voxel crossed by two straight planes (enclosing half of the volume) with no nodes inside the volume
        std::vector<double> distances{-1, -1, -1, -1, -1, -1, -1, -1};   
        GeometryPtrType pVoxel = GenerateHexahedra3D8(distances);

        //Generate the intersecting triangles
        std::vector<std::vector<double>> triangle1{{1,0.5,1.05},{1.05,0.5,0.95},{0.95,0.5,0.95}};  
        std::vector<std::vector<double>> triangle2{{1,0.5,-0.95},{1.05,0.5,-1.05},{0.95,0.5,-1.05}}; 
        std::vector<std::vector<double>> triangle3{{-1,0.5,-0.95},{-0.95,0.5,-1.05},{-1.05,0.5,-1.05}}; 
        std::vector<std::vector<double>> triangle4{{-1,0.5,1.05},{-0.95,0.5,0.95},{-1.05,0.5,0.95}}; 
        std::vector<std::vector<double>> triangle5{{1,-0.5,1.05},{1.05,-0.5,0.95},{0.95,-0.5,0.95}};  
        std::vector<std::vector<double>> triangle6{{1,-0.5,-0.95},{1.05,-0.5,-1.05},{0.95,-0.5,-1.05}}; 
        std::vector<std::vector<double>> triangle7{{-1,-0.5,-0.95},{-0.95,-0.5,-1.05},{-1.05,-0.5,-1.05}}; 
        std::vector<std::vector<double>> triangle8{{-1,-0.5,1.05},{-0.95,-0.5,0.95},{-1.05,-0.5,0.95}}; 
        
        GeometryPtrType pTriangle1 = GenerateTriangle3D3(triangle1);
        GeometryPtrType pTriangle2 = GenerateTriangle3D3(triangle2);
        GeometryPtrType pTriangle3 = GenerateTriangle3D3(triangle3);
        GeometryPtrType pTriangle4 = GenerateTriangle3D3(triangle4);
        GeometryPtrType pTriangle5 = GenerateTriangle3D3(triangle5);
        GeometryPtrType pTriangle6 = GenerateTriangle3D3(triangle6);
        GeometryPtrType pTriangle7 = GenerateTriangle3D3(triangle7);
        GeometryPtrType pTriangle8 = GenerateTriangle3D3(triangle8);

        GeometryArrayType array1;
        array1.push_back(pTriangle1); 
        array1.push_back(pTriangle2);
        array1.push_back(pTriangle3);
        array1.push_back(pTriangle4);
        array1.push_back(pTriangle5); 
        array1.push_back(pTriangle6);
        array1.push_back(pTriangle7);
        array1.push_back(pTriangle8);

        //Call the point utility
        array_1d<double,3> point = QEF::QEF_point(*pVoxel,array1);

        KRATOS_CHECK_NEAR(point[0], 0.0, 1e-8);
        KRATOS_CHECK_NEAR(point[1], 0.0, 1e-8);
        KRATOS_CHECK_NEAR(point[2], 0.0, 1e-8);  
    } 
    
    KRATOS_TEST_CASE_IN_SUITE(QEFmovedCenter, KratosCoreFastSuite) {
        //A voxel with one node inside the volume but centered at (1,1,1)
        std::vector<double> distances{1, -1, -1, -1, -1, -1, -1, -1};   
        GeometryPtrType pVoxel = GenerateUncenteredHexahedra3D8(distances);

        std::vector<std::vector<double>> triangle1{{0,1.5,0.05},{0.05,1.5,-0.05},{-0.05,1.5,-0.05}};
        std::vector<std::vector<double>> triangle2{{0,0.05,1},{0.05,-0.05,1},{-0.05,-0.05,1}};
        std::vector<std::vector<double>> triangle3{{1.75,0.05,0},{1.75,-0.05,0.05},{1.75,-0.05,-0.05}};

        GeometryPtrType pTriangle1 = GenerateTriangle3D3(triangle1);
        GeometryPtrType pTriangle2 = GenerateTriangle3D3(triangle2);
        GeometryPtrType pTriangle3 = GenerateTriangle3D3(triangle3);

        GeometryArrayType array1;
        array1.push_back(pTriangle1);
        array1.push_back(pTriangle2);
        array1.push_back(pTriangle3);

        array_1d<double,3> point = QEF::QEF_point(*pVoxel,array1);

        KRATOS_CHECK_NEAR(point[0], 1.75, 1e-8);
        KRATOS_CHECK_NEAR(point[1], 1.5, 1e-8);
        KRATOS_CHECK_NEAR(point[2], 1.0, 1e-8);     
    }

    KRATOS_TEST_CASE_IN_SUITE(QEF2dofMovedCenter, KratosCoreFastSuite) {
        std::vector<double> distances{1, 1, -1, -1, 1, 1, -1, -1};   
        GeometryPtrType pVoxel = GenerateUncenteredHexahedra3D8(distances);

        std::vector<std::vector<double>> triangle1{{2,1.5,2.05},{2.05,1.5,1.95},{1.95,1.5,1.95}};  
        std::vector<std::vector<double>> triangle2{{2,1.5,0.05},{2.05,1.5,-0.05},{1.95,1.5,-0.05}}; 
        std::vector<std::vector<double>> triangle3{{0,1.5,0.05},{0.05,1.5,-0.05},{-0.05,1.5,-0.05}}; 
        std::vector<std::vector<double>> triangle4{{0,1.5,2.05},{0.05,1.5,1.95},{-0.05,1.5,1.95}}; 
        
        GeometryPtrType pTriangle1 = GenerateTriangle3D3(triangle1);
        GeometryPtrType pTriangle2 = GenerateTriangle3D3(triangle2);
        GeometryPtrType pTriangle3 = GenerateTriangle3D3(triangle3);
        GeometryPtrType pTriangle4 = GenerateTriangle3D3(triangle4);

        GeometryArrayType array1;
        array1.push_back(pTriangle1); 
        array1.push_back(pTriangle2);
        array1.push_back(pTriangle3);
        array1.push_back(pTriangle4);

        array_1d<double,3> point = QEF::QEF_point(*pVoxel,array1);

        KRATOS_CHECK_NEAR(point[0], 1.0, 1e-8);
        KRATOS_CHECK_NEAR(point[1], 1.5, 1e-8);
        KRATOS_CHECK_NEAR(point[2], 1.0, 1e-8);  
    }

}  // namespace Testing.
}  // namespace Kratos.