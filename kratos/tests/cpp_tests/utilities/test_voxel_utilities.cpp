//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//               Kratos default license: kratos/license.txt
//
//  Main authors:    Ariadna Cortes
//
//

// Project includes
#include "containers/model.h"
#include "includes/element.h"
#include "geometries/triangle_3d_3.h"
#include "includes/checks.h"
#include "testing/testing.h"
#include "utilities/voxel_utilities.h"

namespace Kratos {
namespace Testing {
namespace {

    typedef Node<3> NodeType;
    typedef Node<3>::Pointer NodePtrType;
    typedef Geometry<NodeType> GeometryType;
    typedef GeometryType::Pointer GeometryPtrType;
    typedef GeometryType::GeometriesArrayType GeometryArrayType;
    typedef GeometryType::PointsArrayType PointsArrayType;

    /*this functions are creating a model for each element used in the test in order to avoid trouble with repeated Ids.
    This is not optimal but since it is just a test (in real cases, funcions will be used with well-defined models)...
    */
    GeometryPtrType GenerateHexahedra3D8() 
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
        return pElement->pGetGeometry();
    }

    void SetDistances(GeometryPtrType& pVoxel, const std::vector<double>& rDistances) {
        //Add the distances
        PointsArrayType nodes = pVoxel->Points();   
        
        for (int i = 0; i < 8; i++) {
            nodes[i].FastGetSolutionStepValue(DISTANCE) = rDistances[i];
        }
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

    GeometryPtrType GenerateQuadrilateral3D4(std::vector<std::vector<double>>& rNodes, const std::vector<double>& rDistances)
    {
        Model my_model;
        ModelPart& Quadrilater = my_model.CreateModelPart("Quadrilater"); 
        Quadrilater.AddNodalSolutionStepVariable(DISTANCE);  

        Quadrilater.CreateNewNode(1, rNodes[0][0], rNodes[0][1], rNodes[0][2]);
        Quadrilater.CreateNewNode(2, rNodes[1][0], rNodes[1][1], rNodes[1][2]);
        Quadrilater.CreateNewNode(3, rNodes[2][0], rNodes[2][1], rNodes[2][2]);
        Quadrilater.CreateNewNode(4, rNodes[3][0], rNodes[3][1], rNodes[3][2]);
        Properties::Pointer p_properties_1(new Properties(0)); 
        Condition::Pointer pElement = Quadrilater.CreateNewCondition("SurfaceCondition3D4N", 1, {1, 2, 3, 4}, p_properties_1);
        GeometryPtrType pQuadrilater = pElement->pGetGeometry();

        PointsArrayType nodes = pQuadrilater->Points();   
        for (int i = 0; i < 4; i++) {
            nodes[i].FastGetSolutionStepValue(DISTANCE) = rDistances[i];
        }
        return pQuadrilater;
    } 
}  //unnamed namespace

    /******************************************************************************************************
     ******************************************************************************************************/

    KRATOS_TEST_CASE_IN_SUITE(VolumeInsideVoxelNodes, KratosCoreFastSuite) 
    {
        //Generate the HEXAHEDRA3D8
        std::vector<double> distances{1, -1, -1, -1, -1, -1, -1, -1,};   
        GeometryPtrType pVoxel = GenerateHexahedra3D8();
        SetDistances(pVoxel, distances);

        std::vector<double> distances2{1, -1, -0.5, -1, 8, -1, -23, 1,}; 
        GeometryPtrType pVoxel2 = GenerateHexahedra3D8();
        SetDistances(pVoxel2, distances2);

        //Call the volume utility
        double volume = VoxelUtilities::NodesApproximation(*pVoxel); 
        double volume2 = VoxelUtilities::NodesApproximation(*pVoxel2); 

        //Expected output of the function
        const double ExpectedVolume1 = 0.125;
        const double ExpectedVolume2 = 0.375;
        
        KRATOS_CHECK_EQUAL(volume, ExpectedVolume1);
        KRATOS_CHECK_EQUAL(volume2, ExpectedVolume2);
    }

    /******************************************************************************************************
     ******************************************************************************************************/

    KRATOS_TEST_CASE_IN_SUITE(VolumeInsideVoxelEdgesPortion, KratosCoreFastSuite) 
    {
        //Both nodes of an edge are outside
        std::vector<double> distances{-1, -1, -1, -1, -1, -1, -1, -1};   
        GeometryPtrType pVoxel = GenerateHexahedra3D8();
        SetDistances(pVoxel, distances);
        GeometryArrayType array;
        double volume = VoxelUtilities::EdgesPortionApproximation(*pVoxel,array); 
        double ExpectedVolume = 0; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.001);

        //Both nodes of the edge are outside and there is an intersection point (meaning it is tangential)
        std::vector<std::vector<double>> triangle1{{0.5,-1,-0.95},{0.5,-0.95,-1.05},{0.5,-1.05,-1.05}};
        GeometryPtrType pTriangle1 = GenerateTriangle3D3(triangle1);
        array.push_back(pTriangle1);
        volume = VoxelUtilities::EdgesPortionApproximation(*pVoxel,array); 
        ExpectedVolume = 0; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.001);

        /* One node is outside and the other inside but there are no intersection points (meaning the node is a
        tangential point) */
        distances = {1, -1, -1, -1, -1, -1, -1, -1};   
        SetDistances(pVoxel, distances);        array = {};
        volume = VoxelUtilities::EdgesPortionApproximation(*pVoxel,array); 
        ExpectedVolume = 0; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.001);

        // One node is outside and the other inside and there is one intersection point
        distances = {1, -1, -1, -1, -1, -1, -1, -1};   
        SetDistances(pVoxel, distances);
        array.push_back(pTriangle1);
        volume = VoxelUtilities::EdgesPortionApproximation(*pVoxel,array); 
        ExpectedVolume = 0.75/12; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.001);

        /* One node is outside and the other inside and there are two intersection point (meaning one is
        a tangential point) */
        distances = {1, -1, -1, -1, -1, -1, -1, -1};   
        SetDistances(pVoxel, distances);
        std::vector<std::vector<double>> triangle2{{0,-1,-0.95},{0,-0.95,-1.05},{0,-1.05,-1.05}};
        GeometryPtrType pTriangle2 = GenerateTriangle3D3(triangle2);
        array.push_back(pTriangle2);
        volume = VoxelUtilities::EdgesPortionApproximation(*pVoxel,array); 
        ExpectedVolume = 0.5/12; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.001);

        //Both nodes are outside and there are 2 intersection points
        distances = {-1, -1, -1, -1, -1, -1, -1, -1};   
        SetDistances(pVoxel, distances);
        volume = VoxelUtilities::EdgesPortionApproximation(*pVoxel,array); 
        ExpectedVolume = 0.25/12; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.001);

        //Both nodes are outside and there are 3 intersection points (meaning one is tangential)
        volume = VoxelUtilities::EdgesPortionApproximation(*pVoxel,array); 
        std::vector<std::vector<double>> triangle3{{-0.5,-1,-0.95},{-0.5,-0.95,-1.05},{-0.5,-1.05,-1.05}};
        GeometryPtrType pTriangle3 = GenerateTriangle3D3(triangle3);
        array.push_back(pTriangle3);
        ExpectedVolume = 0.25/12; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.001);

        // One node is outside and the other inside and there are three intersection point 
        distances = {1, -1, -1, -1, -1, -1, -1, -1};   
        SetDistances(pVoxel, distances);
        volume = VoxelUtilities::EdgesPortionApproximation(*pVoxel,array); 
        ExpectedVolume = 0.5/12; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.001);

        //Both nodes are inside
        distances = {1, 1, -1, -1, -1, -1, -1, -1};   
        SetDistances(pVoxel, distances);
        array = {};
        volume = VoxelUtilities::EdgesPortionApproximation(*pVoxel,array); 
        ExpectedVolume = 1.0/12; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.001);

        //Both nodes are inside and there is an intersection point (tangential)
        array.push_back(pTriangle1);
        volume = VoxelUtilities::EdgesPortionApproximation(*pVoxel,array); 
        ExpectedVolume = 1.0/12; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.001);

        //Both nodes are inside and there are 2 intersection points
        array.push_back(pTriangle2);
        volume = VoxelUtilities::EdgesPortionApproximation(*pVoxel,array); 
        ExpectedVolume = 0.75/12; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.001);

        //Both nodes are inside and there are 3 intersection points (meaning one is tangential)
        array.push_back(pTriangle3);
        volume = VoxelUtilities::EdgesPortionApproximation(*pVoxel,array); 
        ExpectedVolume = 0.5/12; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.001);
    }

    //A couple of test regarding the real efficiency and limitations of this method

    KRATOS_TEST_CASE_IN_SUITE(VolumeInsideVoxelEdgesPortion2, KratosCoreFastSuite) {
        std::vector<double> distances{1, -1, -1, -1, -1, -1, -1, -1};   
        GeometryPtrType pVoxel = GenerateHexahedra3D8();
        SetDistances(pVoxel, distances);

        //Generate the intersecting triangles
        std::vector<std::vector<double>> triangle1{{-1,0.5,-0.95},{-0.95,0.5,-1.05},{-1.05,0.5,-1.05}}; 
        std::vector<std::vector<double>> triangle2{{-1,-0.95,0},{-0.95,-1.05,0},{-1.05,-1.05,0}}; 
        std::vector<std::vector<double>> triangle3{{-0.5,-0.95,-1},{-0.5,-1.05,-0.95},{-0.5,-1.05,-1.05}}; 
        
        GeometryPtrType pTriangle1 = GenerateTriangle3D3(triangle1);
        GeometryPtrType pTriangle2 = GenerateTriangle3D3(triangle2);
        GeometryPtrType pTriangle3 = GenerateTriangle3D3(triangle3);

        GeometryArrayType array1;
        array1.push_back(pTriangle1); 
        array1.push_back(pTriangle2);
        array1.push_back(pTriangle3);

        //Call the volume utility
        double volume = VoxelUtilities::EdgesPortionApproximation(*pVoxel,array1);
        const double ExpectedVolume = 1.0/8; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01);      
    }
    
    KRATOS_TEST_CASE_IN_SUITE(VolumeInsideVoxelEdgesPortion3, KratosCoreFastSuite) 
    {
        //A voxel crossed by a straight plane with only 2 nodes inside th volume (good case approximation)
        std::vector<double> distances{1, 1, -1, -1, -1, -1, -1, -1};   
        GeometryPtrType pVoxel = GenerateHexahedra3D8();
        SetDistances(pVoxel, distances);

        //Generate the intersecting triangles
        std::vector<std::vector<double>> triangle1{{-1,0.5,-0.95},{-0.95,0.5,-1.05},{-1.05,0.5,-1.05}}; 
        std::vector<std::vector<double>> triangle2{{-1,-0.95,0.99},{-0.95,-1.05,0.99},{-1.05,-1.05,0.99}}; 
        std::vector<std::vector<double>> triangle3{{1,-0.95,0.99},{1.05,-1.05,0.99},{0.95,-1.05,0.99}}; 
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

        //Call the volume utility
        double volume = VoxelUtilities::EdgesPortionApproximation(*pVoxel,array1);
        const double ExpectedVolume = 3.0/8; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01); 
        /*in this case the volume expected to be returned by the method approximation is actually very close 
        to the real volume of the test case, but this will normally not happen with this method*/
    }
     KRATOS_TEST_CASE_IN_SUITE(VolumeInsideVoxelEdgesPortion4, KratosCoreFastSuite) 
     {
        //A voxel crossed by a straight plane with only 2 nodes inside the volume (not good case approximation)
        std::vector<double> distances{1, 1, -1, -1, -1, -1, -1, -1};   
        GeometryPtrType pVoxel = GenerateHexahedra3D8();
        SetDistances(pVoxel, distances);

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

        //Call the volume utility
        double volume = VoxelUtilities::EdgesPortionApproximation(*pVoxel,array1);
        const double ExpectedVolume = 0.291; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01); 
        /*in this case the volume returned by the method approximation is not clse to the real volume of the 
        test case, since we would expect a real volume circa 0.1875 */
    }

    KRATOS_TEST_CASE_IN_SUITE(VolumeInsideVoxelEdgesPortion5, KratosCoreFastSuite) 
    {
        //A voxel crossed by a straight plane with 4 nodes inside the volume
        std::vector<double> distances{1, 1, -1, -1, 1, 1, -1, -1};   
        GeometryPtrType pVoxel = GenerateHexahedra3D8();
        SetDistances(pVoxel, distances);

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

        //Call the volume utility
        double volume = VoxelUtilities::EdgesPortionApproximation(*pVoxel,array1);
        const double ExpectedVolume = 7.0/12; //0.5833
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01); 
        /*in this case the volume returned by the method approximation is not close to the real volume of the 
        test case, since we would expect a real volume circa 0.75 */
    } 

    KRATOS_TEST_CASE_IN_SUITE(VolumeInsideVoxelEdgesPortion6, KratosCoreFastSuite) 
    {
        //A voxel crossed by two straight planes (enclosing half of the volume) with no nodes inside the volume
        std::vector<double> distances{-1, -1, -1, -1, -1, -1, -1, -1};   
        GeometryPtrType pVoxel = GenerateHexahedra3D8();
        SetDistances(pVoxel, distances);

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

        //Call the volume utility
        double volume = VoxelUtilities::EdgesPortionApproximation(*pVoxel,array1);
        const double ExpectedVolume = 2.0/12; //0.1666
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01); 
        /*in this case the volume returned by the method approximation is not close to the real volume of the 
        test case, since we would expect a real volume circa 0.5 */
    } 

    /******************************************************************************************************
     ******************************************************************************************************/

    KRATOS_TEST_CASE_IN_SUITE(VolumeInsideVoxelEdgesPortionQuadrilateral, KratosCoreFastSuite)
    {
        std::vector<std::vector<double>> quad{{1,1,-1},{-1,1,-1},{-1,-1,1}, {1,-1,1}};  
        std::vector<double> distances{-1, -1, -1, -1}; 
        GeometryPtrType pFace = GenerateQuadrilateral3D4(quad,distances);

        std::vector<std::vector<double>> triangle1{{-0.5,-1,0.95},{-0.5,-0.95,1.05},{-0.5,-1.05,1.05}};
        GeometryPtrType pTriangle1 = GenerateTriangle3D3(triangle1);

        GeometryArrayType array1;
        array1.push_back(pTriangle1);
        double volume = VoxelUtilities::EdgesPortionApproximation(*pFace,array1);
        const double ExpectedVolume = 0; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01);
    }

    KRATOS_TEST_CASE_IN_SUITE(VolumeInsideVoxelEdgesPortionQuadrilateral2, KratosCoreFastSuite) 
    {
        std::vector<std::vector<double>> quad{{1,1,-1},{-1,1,-1},{-1,-1,1}, {1,-1,1}};  
        std::vector<double> distances{1, -1, -1, -1}; 
        GeometryPtrType pFace = GenerateQuadrilateral3D4(quad,distances);
        
        std::vector<std::vector<double>> triangle1{{0.5,1,-0.95},{0.5,0.95,-1.05},{0.5,1.05,-1.05}};
        GeometryPtrType pTriangle1 = GenerateTriangle3D3(triangle1);

        GeometryArrayType array1;
        array1.push_back(pTriangle1);
        double volume = VoxelUtilities::EdgesPortionApproximation(*pFace,array1);
        const double ExpectedVolume = 0.25/4; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01);
    }

     /******************************************************************************************************
     *****************************************************************************************************/

    KRATOS_TEST_CASE_IN_SUITE(VolumeInsideVoxelFaceArea, KratosCoreFastSuite) 
    { 
        //Quadrilater with one node inside the volume (different cases)
        std::vector<std::vector<double>> quad{{1,1,0},{-1,1,0},{-1,-1,0}, {1,-1,0}};  
        std::vector<double> distances{1, -1, -1, -1}; 
        GeometryPtrType pFace = GenerateQuadrilateral3D4(quad,distances);
        
        std::vector<std::vector<double>> triangle1{{0,1,0.05},{0,0.95,-0.05},{0,1.05,-0.05}};
        std::vector<std::vector<double>> triangle2{{1,0.5,0.05},{0.95,0.5,-0.05},{1.05,0.5,-0.05}};
        std::vector<std::vector<double>> triangle3{{0,-1,0.05},{0,-1.05,-0.05},{0,-0.95,-0.05}};
        std::vector<std::vector<double>> triangle4{{0.95,1,0.05},{0.95,0.95,-0.05},{0.95,1.05,-0.05}};
        std::vector<std::vector<double>> triangle5{{0.95,-1,0.05},{0.95,-1.05,-0.05},{0.95,-0.95,-0.05}};
        GeometryPtrType pTriangle1 = GenerateTriangle3D3(triangle1);
        GeometryPtrType pTriangle2 = GenerateTriangle3D3(triangle2);
        GeometryPtrType pTriangle3 = GenerateTriangle3D3(triangle3);
        GeometryPtrType pTriangle4 = GenerateTriangle3D3(triangle4);
        GeometryPtrType pTriangle5 = GenerateTriangle3D3(triangle5);

        GeometryArrayType array1;
        array1.push_back(pTriangle1);
        array1.push_back(pTriangle2);
        double volume = VoxelUtilities::FaceArea(*pFace,array1);
        double ExpectedVolume = 1.0/16; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01);

        //Two nodes inside the volume
        distances = {-1, 1, 1, -1}; 
        pFace = GenerateQuadrilateral3D4(quad,distances);

        array1.clear();
        array1.push_back(pTriangle1);
        array1.push_back(pTriangle3);
        volume = VoxelUtilities::FaceArea(*pFace,array1);
        ExpectedVolume = 2.0/4; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01);

        array1.clear();
        array1.push_back(pTriangle4);
        array1.push_back(pTriangle5);
        volume = VoxelUtilities::FaceArea(*pFace,array1);
        ExpectedVolume = 0.975; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01);

        //Three nodes inside the volume
        distances = {-1, 1, 1, 1}; 
        pFace = GenerateQuadrilateral3D4(quad,distances);

        array1.clear();
        array1.push_back(pTriangle2);
        array1.push_back(pTriangle1);
        volume = VoxelUtilities::FaceArea(*pFace,array1);
        ExpectedVolume = 15.0/16; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01);

        //Full quadrilateral with no intersections
        distances = {1, 1, 1, 1}; 
        pFace = GenerateQuadrilateral3D4(quad,distances);

        array1.clear();
        volume = VoxelUtilities::FaceArea(*pFace,array1);
        ExpectedVolume = 1.0; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01);

        //Empty quadrilateral
        distances = {-1, -1, -1, -1}; 
        pFace = GenerateQuadrilateral3D4(quad,distances);
        volume = VoxelUtilities::FaceArea(*pFace,array1);
        ExpectedVolume = 0.0; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01);
    }

    /******************************************************************************************************
     ******************************************************************************************************/

    KRATOS_TEST_CASE_IN_SUITE(VolumeInsideVoxelFaceArea0, KratosCoreFastSuite) 
    { 
        GeometryArrayType array1;
        std::vector<std::vector<double>> quad{{1,1,0},{-1,1,0},{-1,-1,0}, {1,-1,0}};  
       
        //Full quadrilateral with no intersections
        std::vector<double> distances{1, 1, 1, 1}; 
        GeometryPtrType  pFace = GenerateQuadrilateral3D4(quad,distances);

        double volume = VoxelUtilities::FaceArea(*pFace,array1);
        double ExpectedVolume = 1.0; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01);

        //Empty quadrilateral
        distances = {-1, -1, -1, -1}; 
        pFace = GenerateQuadrilateral3D4(quad,distances);
        volume = VoxelUtilities::FaceArea(*pFace,array1);
        ExpectedVolume = 0.0; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01);
    }

    KRATOS_TEST_CASE_IN_SUITE(VolumeInsideVoxelFaceArea1, KratosCoreFastSuite) 
    { 
        //Quadrilater with one node inside the volume
        std::vector<std::vector<double>> quad{{1,1,0},{-1,1,0},{-1,-1,0}, {1,-1,0}};  
        std::vector<double> distances{1, -1, -1, -1}; 
        GeometryPtrType pFace = GenerateQuadrilateral3D4(quad,distances);
        
        std::vector<std::vector<double>> triangle1{{0,1,0.05},{0,0.95,-0.05},{0,1.05,-0.05}};
        std::vector<std::vector<double>> triangle2{{1,0.5,0.05},{0.95,0.5,-0.05},{1.05,0.5,-0.05}};

        GeometryPtrType pTriangle1 = GenerateTriangle3D3(triangle1);
        GeometryPtrType pTriangle2 = GenerateTriangle3D3(triangle2);

        GeometryArrayType array1;
        array1.push_back(pTriangle1);
        array1.push_back(pTriangle2);
        double volume = VoxelUtilities::FaceArea(*pFace,array1);
        double ExpectedVolume = 1.0/16; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01);
    }

    KRATOS_TEST_CASE_IN_SUITE(VolumeInsideVoxelFaceArea2, KratosCoreFastSuite) 
    { 
        //Quadrilater with one node inside the volume
        std::vector<std::vector<double>> quad{{1,1,0},{-1,1,0},{-1,-1,0}, {1,-1,0}};  
        std::vector<double> distances{1, -1, -1, -1}; 
        GeometryPtrType pFace = GenerateQuadrilateral3D4(quad,distances);
        
        std::vector<std::vector<double>> triangle1{{-0.99,1,0.05},{-0.99,0.95,-0.05},{-0.99,1.05,-0.05}};
        std::vector<std::vector<double>> triangle2{{1,-0.99,0.05},{0.95,-0.99,-0.05},{1.05,-0.99,-0.05}};

        GeometryPtrType pTriangle1 = GenerateTriangle3D3(triangle1);
        GeometryPtrType pTriangle2 = GenerateTriangle3D3(triangle2);

        GeometryArrayType array1;
        array1.push_back(pTriangle1);
        array1.push_back(pTriangle2);
        double volume = VoxelUtilities::FaceArea(*pFace,array1);
        double ExpectedVolume = 0.495; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.001);
    }

     KRATOS_TEST_CASE_IN_SUITE(VolumeInsideVoxelFaceArea3, KratosCoreFastSuite) 
     { 
        //Quadrilater with one node inside the volume
        std::vector<std::vector<double>> quad{{1,1,0},{-1,1,0},{-1,-1,0}, {1,-1,0}};          
        std::vector<std::vector<double>> triangle1{{0,1,0.05},{0,0.95,-0.05},{0,1.05,-0.05}};
        std::vector<std::vector<double>> triangle3{{0,-1,0.05},{0,-1.05,-0.05},{0,-0.95,-0.05}};

        GeometryPtrType pTriangle1 = GenerateTriangle3D3(triangle1);
        GeometryPtrType pTriangle3 = GenerateTriangle3D3(triangle3);

        //Two nodes inside the volume
        std::vector<double> distances{-1, 1, 1, -1}; 
         GeometryPtrType pFace = GenerateQuadrilateral3D4(quad,distances);

        GeometryArrayType array1;
        array1.push_back(pTriangle1);
        array1.push_back(pTriangle3);
        double volume = VoxelUtilities::FaceArea(*pFace,array1);
        double ExpectedVolume = 2.0/4; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01);
    }

     KRATOS_TEST_CASE_IN_SUITE(VolumeInsideVoxelFaceArea4, KratosCoreFastSuite) 
     { 
        //Two nodes inside the volume
        std::vector<std::vector<double>> quad{{1,1,0},{-1,1,0},{-1,-1,0}, {1,-1,0}};  
        std::vector<std::vector<double>> triangle4{{0.95,1,0.05},{0.95,0.95,-0.05},{0.95,1.05,-0.05}};
        std::vector<std::vector<double>> triangle5{{0.95,-1,0.05},{0.95,-1.05,-0.05},{0.95,-0.95,-0.05}};
        GeometryPtrType pTriangle4 = GenerateTriangle3D3(triangle4);
        GeometryPtrType pTriangle5 = GenerateTriangle3D3(triangle5);

        std::vector<double> distances{-1, 1, 1, -1}; 
         GeometryPtrType pFace = GenerateQuadrilateral3D4(quad,distances);

        GeometryArrayType array1;
        array1.push_back(pTriangle4);
        array1.push_back(pTriangle5);
        double volume = VoxelUtilities::FaceArea(*pFace,array1);
        double ExpectedVolume = 0.975; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01);
    }

    KRATOS_TEST_CASE_IN_SUITE(VolumeInsideVoxelFaceArea5, KratosCoreFastSuite) 
    { 
        //Two nodes inside the volume
        std::vector<std::vector<double>> quad{{1,1,0},{-1,1,0},{-1,-1,0}, {1,-1,0}};  
        std::vector<std::vector<double>> triangle4{{0.99999,1,0.05},{0.9999,0.95,-0.05},{0.9999,1.05,-0.05}};
        std::vector<std::vector<double>> triangle5{{-0.5,-1,0.05},{-0.5,-1.05,-0.05},{-0.5,-0.95,-0.05}};
        GeometryPtrType pTriangle4 = GenerateTriangle3D3(triangle4);
        GeometryPtrType pTriangle5 = GenerateTriangle3D3(triangle5);

        std::vector<double> distances{-1, 1, 1, -1}; 
         GeometryPtrType pFace = GenerateQuadrilateral3D4(quad,distances);

        GeometryArrayType array1;
        array1.push_back(pTriangle4);
        array1.push_back(pTriangle5);
        double volume = VoxelUtilities::FaceArea(*pFace,array1);
        double ExpectedVolume = 0.625; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.001);
    }

    KRATOS_TEST_CASE_IN_SUITE(VolumeInsideVoxelFaceArea6, KratosCoreFastSuite) 
    { 
        //Two nodes inside the volume
        std::vector<std::vector<double>> quad{{1,1,0},{-1,1,0},{-1,-1,0}, {1,-1,0}};  
        std::vector<std::vector<double>> triangle1{{0,1,0.05},{0,0.95,-0.05},{0,1.05,-0.05}};
        std::vector<std::vector<double>> triangle2{{-1,0.5,0.05},{-1.05,0.5,-0.05},{-0.95,0.5,-0.05}};
        std::vector<std::vector<double>> triangle3{{0,-1,0.05},{0,-1.05,-0.05},{0,-0.95,-0.05}};
        std::vector<std::vector<double>> triangleaux{{1,0.5,0.05},{0.95,0.5,-0.05},{1.05,0.5,-0.05}};

        GeometryPtrType pTriangle1 = GenerateTriangle3D3(triangle1);
        GeometryPtrType pTriangle2 = GenerateTriangle3D3(triangle2);
        GeometryPtrType pTriangle3 = GenerateTriangle3D3(triangle3);
        GeometryPtrType pTriangleaux = GenerateTriangle3D3(triangleaux);

        std::vector<double> distances{1, -1, 1, -1}; 
         GeometryPtrType pFace = GenerateQuadrilateral3D4(quad,distances);

        GeometryArrayType array1;
        array1.push_back(pTriangle1);
        array1.push_back(pTriangle2);
        array1.push_back(pTriangle3);
        array1.push_back(pTriangleaux);
        double volume = VoxelUtilities::FaceArea(*pFace,array1);
        double ExpectedVolume = 0.25; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.001);
    }

     KRATOS_TEST_CASE_IN_SUITE(VolumeInsideVoxelFaceArea7, KratosCoreFastSuite) 
     { 
        //three nodes inside the volume
        std::vector<std::vector<double>> quad{{1,1,0},{-1,1,0},{-1,-1,0}, {1,-1,0}};  
        std::vector<std::vector<double>> triangle1{{0,1,0.05},{0,0.95,-0.05},{0,1.05,-0.05}};
        std::vector<std::vector<double>> triangle2{{1,0.5,0.05},{0.95,0.5,-0.05},{1.05,0.5,-0.05}};
        GeometryPtrType pTriangle1 = GenerateTriangle3D3(triangle1);
        GeometryPtrType pTriangle2 = GenerateTriangle3D3(triangle2);

        GeometryArrayType array1;
        array1.push_back(pTriangle1);
        array1.push_back(pTriangle2);
        std::vector<double> distances{-1, 1, 1, 1}; 
        GeometryPtrType pFace = GenerateQuadrilateral3D4(quad,distances);

        double volume = VoxelUtilities::FaceArea(*pFace,array1);
        double ExpectedVolume = 15.0/16; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01);
    }

    KRATOS_TEST_CASE_IN_SUITE(VolumeInsideVoxelFaceArea8, KratosCoreFastSuite) 
    { 
        //QUADRILATERAL with one node inside the volume 
        std::vector<std::vector<double>> quad{{5,3,0},{0,3,0},{0,0,0}, {5,0,0}};  
        std::vector<double> distances{1, -1, -1, -1}; 
        GeometryPtrType pFace = GenerateQuadrilateral3D4(quad,distances);
        
        std::vector<std::vector<double>> triangle1{{1,3,0.05},{1,2.95,-0.05},{1,3.05,-0.05}};
        std::vector<std::vector<double>> triangle2{{5,1.5,0.05},{4.95,1.5,-0.05},{5.05,1.5,-0.05}};

        GeometryPtrType pTriangle1 = GenerateTriangle3D3(triangle1);
        GeometryPtrType pTriangle2 = GenerateTriangle3D3(triangle2);

        GeometryArrayType array1;
        array1.push_back(pTriangle1);
        array1.push_back(pTriangle2);
        double volume = VoxelUtilities::FaceArea(*pFace,array1);
        double ExpectedVolume = 0.2; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01);
    }

    KRATOS_TEST_CASE_IN_SUITE(VolumeInsideVoxelFaceArea9, KratosCoreFastSuite) 
    { 
        //Quadrilater with one node inside the volume 
        std::vector<std::vector<double>> quad{{1,1,0},{-1,1,0},{-3,-3,0}, {1,-1,0}};  
        std::vector<double> distances{1, -1, -1, -1}; 
        GeometryPtrType pFace = GenerateQuadrilateral3D4(quad,distances);
        
        std::vector<std::vector<double>> triangle1{{0,1,0.05},{0,0.95,-0.05},{0,1.05,-0.05}};
        std::vector<std::vector<double>> triangle2{{1,0.5,0.05},{0.95,0.5,-0.05},{1.05,0.5,-0.05}};

        GeometryPtrType pTriangle1 = GenerateTriangle3D3(triangle1);
        GeometryPtrType pTriangle2 = GenerateTriangle3D3(triangle2);

        GeometryArrayType array1;
        array1.push_back(pTriangle1);
        array1.push_back(pTriangle2);
        double volume = VoxelUtilities::FaceArea(*pFace,array1);
        double ExpectedVolume = 1.0/32; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.001);
    }

     KRATOS_TEST_CASE_IN_SUITE(VolumeInsideVoxelFaceArea10, KratosCoreFastSuite)
     { 
        //Two nodes inside the volume
        std::vector<std::vector<double>> quad{{2,1,0},{-1,1,0},{-1,-1,0}, {1.5,-1,0}};  
        std::vector<std::vector<double>> triangle4{{0.95,1,0.05},{0.95,0.95,-0.05},{0.95,1.05,-0.05}};
        std::vector<std::vector<double>> triangle5{{0.95,-1,0.05},{0.95,-1.05,-0.05},{0.95,-0.95,-0.05}};
        GeometryPtrType pTriangle4 = GenerateTriangle3D3(triangle4);
        GeometryPtrType pTriangle5 = GenerateTriangle3D3(triangle5);

        std::vector<double> distances{-1, 1, 1, -1}; 
         GeometryPtrType pFace = GenerateQuadrilateral3D4(quad,distances);

        GeometryArrayType array1;
        array1.push_back(pTriangle4);
        array1.push_back(pTriangle5);
        double volume = VoxelUtilities::FaceArea(*pFace,array1);
        double ExpectedVolume = 0.709; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.001);
    }

    KRATOS_TEST_CASE_IN_SUITE(VolumeInsideVoxelFaceArea11, KratosCoreFastSuite) 
    { 
        //Two nodes inside the volume
        std::vector<std::vector<double>> quad{{2,1,0},{-2,1,0},{-1.5,-1,0}, {1.5,-1,0}};  
        std::vector<std::vector<double>> triangle4{{0.95,1,0.05},{0.95,0.95,-0.05},{0.95,1.05,-0.05}};
        std::vector<std::vector<double>> triangle5{{0.95,-1,0.05},{0.95,-1.05,-0.05},{0.95,-0.95,-0.05}};
        GeometryPtrType pTriangle4 = GenerateTriangle3D3(triangle4);
        GeometryPtrType pTriangle5 = GenerateTriangle3D3(triangle5);

        std::vector<double> distances{-1, 1, 1, -1}; 
         GeometryPtrType pFace = GenerateQuadrilateral3D4(quad,distances);

        GeometryArrayType array1;
        array1.push_back(pTriangle4);
        array1.push_back(pTriangle5);
        double volume = VoxelUtilities::FaceArea(*pFace,array1);
        double ExpectedVolume = 0.772; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.001);
    }

}  // namespace Testing.
}  // namespace Kratos.