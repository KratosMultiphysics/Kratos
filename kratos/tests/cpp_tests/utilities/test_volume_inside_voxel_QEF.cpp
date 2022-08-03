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
#include "containers/model.h"
#include "includes/element.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/triangle_3d_3.h"
#include "includes/checks.h"
#include "testing/testing.h"
#include "utilities/volume_inside_voxel_QEF.h"

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
    GeometryPtrType QEFVolumeGenerateHexahedra3D8(const std::vector<double>& rDistances) 
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

    //rNodes is a 3*3 matrix representin de (x,y,z) coordinates of each of the 3 triangle nodes
    GeometryPtrType QEFVolumeGenerateTriangle3D3(std::vector<std::vector<double>>& rNodes)
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

    //Basic cases tests (all nodes inside/ouside)
    KRATOS_TEST_CASE_IN_SUITE(QEFSimpleNodesFull, KratosCoreFastSuite) {
        //A voxel crossed by a straight plane with only 2 nodes inside th volume (Bad case approximation)
        std::vector<double> distances{1, 1, 1, 1, 1, 1, 1, 1};   
        GeometryPtrType pVoxel = QEFVolumeGenerateHexahedra3D8(distances);

        GeometryArrayType array1;

        //Call the volume utility
        double volume = VolumeInsideVoxelQEF::SimpleNodesQEFApproximation(*pVoxel,array1);
        const double ExpectedVolume = 1.0; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01); 
        /*in this case the volume expected to be returned by the method approximation is actually very close 
        to the real volume of the test case, but this will normally not happen*/
    }

    KRATOS_TEST_CASE_IN_SUITE(QEFSimpleNodesEmpty, KratosCoreFastSuite) {
        //A voxel crossed by a straight plane with only 2 nodes inside th volume (Bad case approximation)
        std::vector<double> distances{1, 1, 1, 1, 1, 1, 1, 1};   
        GeometryPtrType pVoxel = QEFVolumeGenerateHexahedra3D8(distances);

        GeometryArrayType array1;

        //Call the volume utility
        double volume = VolumeInsideVoxelQEF::SimpleNodesQEFApproximation(*pVoxel,array1);
        const double ExpectedVolume = 1.0; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01); 
        /*in this case the volume expected to be returned by the method approximation is actually very close 
        to the real volume of the test case, but this will normally not happen*/
    }

    //A couple of test regarding the real efficiency and limitations of this method

    KRATOS_TEST_CASE_IN_SUITE(QEFSimpleNodes0, KratosCoreFastSuite) {
        std::vector<double> distances{1, -1, -1, -1, -1, -1, -1, -1};   
        GeometryPtrType pVoxel = QEFVolumeGenerateHexahedra3D8(distances);

        //Generate the intersecting triangles
        std::vector<std::vector<double>> triangle1{{-1,0.5,-0.95},{-0.95,0.5,-1.05},{-1.05,0.5,-1.05}}; 
        std::vector<std::vector<double>> triangle2{{-1,-0.95,0},{-0.95,-1.05,0},{-1.05,-1.05,0}}; 
        std::vector<std::vector<double>> triangle3{{-0.5,-0.95,-1},{-0.5,-1.05,-0.95},{-0.5,-1.05,-1.05}}; 
        
        GeometryPtrType pTriangle1 = QEFVolumeGenerateTriangle3D3(triangle1);
        GeometryPtrType pTriangle2 = QEFVolumeGenerateTriangle3D3(triangle2);
        GeometryPtrType pTriangle3 = QEFVolumeGenerateTriangle3D3(triangle3);

        GeometryArrayType array1;
        array1.push_back(pTriangle1); 
        array1.push_back(pTriangle2);
        array1.push_back(pTriangle3);

        //Call the volume utility
        double volume = VolumeInsideVoxelQEF::SimpleNodesQEFApproximation(*pVoxel,array1);
        const double ExpectedVolume = 1.0/8; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01);      
    }

     KRATOS_TEST_CASE_IN_SUITE(QEFSimpleNodes1, KratosCoreFastSuite) {
        //A voxel crossed by a straight plane with only 2 nodes inside the volume (ood case approximation)
        std::vector<double> distances{1, 1, -1, -1, -1, -1, -1, -1};   
        GeometryPtrType pVoxel = QEFVolumeGenerateHexahedra3D8(distances);

        //Generate the intersecting triangles
        std::vector<std::vector<double>> triangle1{{-1,0.0,-0.95},{-0.95,0.0,-1.05},{-1.05,0.0,-1.05}}; 
        std::vector<std::vector<double>> triangle2{{-1,-0.95,0},{-0.95,-1.05,0},{-1.05,-1.05,0}}; 
        std::vector<std::vector<double>> triangle3{{1,-0.95,0},{1.05,-1.05,0},{0.95,-1.05,0}}; 
        std::vector<std::vector<double>> triangle4{{1,0.0,-0.95},{1.05,0.0,-1.05},{0.95,0.0,-1.05}}; 
        
        GeometryPtrType pTriangle1 = QEFVolumeGenerateTriangle3D3(triangle1);
        GeometryPtrType pTriangle2 = QEFVolumeGenerateTriangle3D3(triangle2);
        GeometryPtrType pTriangle3 = QEFVolumeGenerateTriangle3D3(triangle3);
        GeometryPtrType pTriangle4 = QEFVolumeGenerateTriangle3D3(triangle4);

        GeometryArrayType array1;
        array1.push_back(pTriangle1); 
        array1.push_back(pTriangle2);
        array1.push_back(pTriangle3);
        array1.push_back(pTriangle4);

        //Call the volume utility
        double volume = VolumeInsideVoxelQEF::SimpleNodesQEFApproximation(*pVoxel,array1);
        const double ExpectedVolume = 0.25; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01); 
        /*in this case the volume returned by the method approximation exactly the same to the real volume of the 
        test case, since this is a perfect choosen case*/
    }
    
    KRATOS_TEST_CASE_IN_SUITE(QEFSimpleNodes2, KratosCoreFastSuite) {
        //A voxel crossed by a straight plane with only 2 nodes inside th volume (Bad case approximation)
        std::vector<double> distances{1, 1, -1, -1, -1, -1, -1, -1};   
        GeometryPtrType pVoxel = QEFVolumeGenerateHexahedra3D8(distances);

        //Generate the intersecting triangles
        std::vector<std::vector<double>> triangle1{{-1,0.5,-0.95},{-0.95,0.5,-1.05},{-1.05,0.5,-1.05}}; 
        std::vector<std::vector<double>> triangle2{{-1,-0.95,0.99},{-0.95,-1.05,0.99},{-1.05,-1.05,0.99}}; 
        std::vector<std::vector<double>> triangle3{{1,-0.95,0.99},{1.05,-1.05,0.99},{0.95,-1.05,0.99}}; 
        std::vector<std::vector<double>> triangle4{{1,0.5,-0.95},{1.05,0.5,-1.05},{0.95,0.5,-1.05}}; 
        
        GeometryPtrType pTriangle1 = QEFVolumeGenerateTriangle3D3(triangle1);
        GeometryPtrType pTriangle2 = QEFVolumeGenerateTriangle3D3(triangle2);
        GeometryPtrType pTriangle3 = QEFVolumeGenerateTriangle3D3(triangle3);
        GeometryPtrType pTriangle4 = QEFVolumeGenerateTriangle3D3(triangle4);

        GeometryArrayType array1;
        array1.push_back(pTriangle1); 
        array1.push_back(pTriangle2);
        array1.push_back(pTriangle3);
        array1.push_back(pTriangle4);

        //Call the volume utility
        double volume = VolumeInsideVoxelQEF::SimpleNodesQEFApproximation(*pVoxel,array1);
        const double ExpectedVolume = 3.0/8; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01); 
         /*in this case the volume returned by the method approximation is not close to the real volume of the 
        test case, since we would expect a real volume circa 0.75 */
    }

    KRATOS_TEST_CASE_IN_SUITE(QEFSimpleNodes3, KratosCoreFastSuite) {
        //A voxel crossed by a straight plane with 4 nodes inside the volume
        std::vector<double> distances{1, 1, -1, -1, 1, 1, -1, -1};   
        GeometryPtrType pVoxel = QEFVolumeGenerateHexahedra3D8(distances);

        //Generate the intersecting triangles
        std::vector<std::vector<double>> triangle1{{1,0.5,1.05},{1.05,0.5,0.95},{0.95,0.5,0.95}};  
        std::vector<std::vector<double>> triangle2{{1,0.5,-0.95},{1.05,0.5,-1.05},{0.95,0.5,-1.05}}; 
        std::vector<std::vector<double>> triangle3{{-1,0.5,-0.95},{-0.95,0.5,-1.05},{-1.05,0.5,-1.05}}; 
        std::vector<std::vector<double>> triangle4{{-1,0.5,1.05},{-0.95,0.5,0.95},{-1.05,0.5,0.95}}; 
        
        GeometryPtrType pTriangle1 = QEFVolumeGenerateTriangle3D3(triangle1);
        GeometryPtrType pTriangle2 = QEFVolumeGenerateTriangle3D3(triangle2);
        GeometryPtrType pTriangle3 = QEFVolumeGenerateTriangle3D3(triangle3);
        GeometryPtrType pTriangle4 = QEFVolumeGenerateTriangle3D3(triangle4);

        GeometryArrayType array1;
        array1.push_back(pTriangle1); 
        array1.push_back(pTriangle2);
        array1.push_back(pTriangle3);
        array1.push_back(pTriangle4);

        //Call the volume utility
        double volume = VolumeInsideVoxelQEF::SimpleNodesQEFApproximation(*pVoxel,array1);
        const double ExpectedVolume = 7.0/12; //0.5833
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01); 
        /*in this case the volume returned by the method approximation is not close to the real volume of the 
        test case, since we would expect a real volume circa 0.75 */
    } 

    KRATOS_TEST_CASE_IN_SUITE(QEFSimpleNodes4, KratosCoreFastSuite) {
        //A voxel crossed by two straight planes (enclosing half of the volume) with no nodes inside the volume
        std::vector<double> distances{-1, -1, -1, -1, -1, -1, -1, -1};   
        GeometryPtrType pVoxel = QEFVolumeGenerateHexahedra3D8(distances);

        //Generate the intersecting triangles
        std::vector<std::vector<double>> triangle1{{1,0.5,1.05},{1.05,0.5,0.95},{0.95,0.5,0.95}};  
        std::vector<std::vector<double>> triangle2{{1,0.5,-0.95},{1.05,0.5,-1.05},{0.95,0.5,-1.05}}; 
        std::vector<std::vector<double>> triangle3{{-1,0.5,-0.95},{-0.95,0.5,-1.05},{-1.05,0.5,-1.05}}; 
        std::vector<std::vector<double>> triangle4{{-1,0.5,1.05},{-0.95,0.5,0.95},{-1.05,0.5,0.95}}; 
        std::vector<std::vector<double>> triangle5{{1,-0.5,1.05},{1.05,-0.5,0.95},{0.95,-0.5,0.95}};  
        std::vector<std::vector<double>> triangle6{{1,-0.5,-0.95},{1.05,-0.5,-1.05},{0.95,-0.5,-1.05}}; 
        std::vector<std::vector<double>> triangle7{{-1,-0.5,-0.95},{-0.95,-0.5,-1.05},{-1.05,-0.5,-1.05}}; 
        std::vector<std::vector<double>> triangle8{{-1,-0.5,1.05},{-0.95,-0.5,0.95},{-1.05,-0.5,0.95}}; 
        
        GeometryPtrType pTriangle1 = QEFVolumeGenerateTriangle3D3(triangle1);
        GeometryPtrType pTriangle2 = QEFVolumeGenerateTriangle3D3(triangle2);
        GeometryPtrType pTriangle3 = QEFVolumeGenerateTriangle3D3(triangle3);
        GeometryPtrType pTriangle4 = QEFVolumeGenerateTriangle3D3(triangle4);
        GeometryPtrType pTriangle5 = QEFVolumeGenerateTriangle3D3(triangle5);
        GeometryPtrType pTriangle6 = QEFVolumeGenerateTriangle3D3(triangle6);
        GeometryPtrType pTriangle7 = QEFVolumeGenerateTriangle3D3(triangle7);
        GeometryPtrType pTriangle8 = QEFVolumeGenerateTriangle3D3(triangle8);

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
        double volume = VolumeInsideVoxelQEF::SimpleNodesQEFApproximation(*pVoxel,array1);
        const double ExpectedVolume = 0; //no nodes inside
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01); 
        /*in this case the volume returned by the method approximation is not close to the real volume of the 
        test case, since we would expect a real volume circa 0.5 */
    } 

    /******************************************************************************************************
     ******************************************************************************************************/

    //Basic cases tests (all nodes inside/ouside)
    KRATOS_TEST_CASE_IN_SUITE(QEFApproximationsFull, KratosCoreFastSuite) {
        //A voxel crossed by a straight plane with only 2 nodes inside th volume (Bad case approximation)
        std::vector<double> distances{1, 1, 1, 1, 1, 1, 1, 1};   
        GeometryPtrType pVoxel = QEFVolumeGenerateHexahedra3D8(distances);

        GeometryArrayType array1;

        //Call the volume utility
        double volume = VolumeInsideVoxelQEF::FacesPortionQEFApproximation(*pVoxel,array1);
        const double ExpectedVolume = 1.0; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01); 
        volume = VolumeInsideVoxelQEF::GeometricalQEFApproximation(*pVoxel,array1);
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01); 
        /*in this case the volume expected to be returned by the method approximation is actually very close 
        to the real volume of the test case, but this will normally not happen*/
    }

    KRATOS_TEST_CASE_IN_SUITE(QEFApproximationsEmpty, KratosCoreFastSuite) {
        //A voxel crossed by a straight plane with only 2 nodes inside th volume (Bad case approximation)
        std::vector<double> distances{1, 1, 1, 1, 1, 1, 1, 1};   
        GeometryPtrType pVoxel = QEFVolumeGenerateHexahedra3D8(distances);

        GeometryArrayType array1;

        //Call the volume utility
        double volume = VolumeInsideVoxelQEF::FacesPortionQEFApproximation(*pVoxel,array1);
        const double ExpectedVolume = 1.0; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01); 
        volume = VolumeInsideVoxelQEF::GeometricalQEFApproximation(*pVoxel,array1);
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01); 
        /*in this case the volume expected to be returned by the method approximation is actually very close 
        to the real volume of the test case, but this will normally not happen*/
    }

    //A couple of test regarding the real efficiency and limitations of this method
    KRATOS_TEST_CASE_IN_SUITE(QEFApproximations1, KratosCoreFastSuite) {
        //voxel crossed by a plane with only 1 node inside the volume
        std::vector<double> distances{1, -1, -1, -1, -1, -1, -1, -1};   
        GeometryPtrType pVoxel = QEFVolumeGenerateHexahedra3D8(distances);

        //Generate the intersecting triangles
        std::vector<std::vector<double>> triangle1{{-1,0.5,-0.95},{-0.95,0.5,-1.05},{-1.05,0.5,-1.05}}; 
        std::vector<std::vector<double>> triangle2{{-1,-0.95,0},{-0.95,-1.05,0},{-1.05,-1.05,0}}; 
        std::vector<std::vector<double>> triangle3{{-0.5,-0.95,-1},{-0.5,-1.05,-0.95},{-0.5,-1.05,-1.05}}; 
        
        GeometryPtrType pTriangle1 = QEFVolumeGenerateTriangle3D3(triangle1);
        GeometryPtrType pTriangle2 = QEFVolumeGenerateTriangle3D3(triangle2);
        GeometryPtrType pTriangle3 = QEFVolumeGenerateTriangle3D3(triangle3);

        GeometryArrayType array1;
        array1.push_back(pTriangle1); 
        array1.push_back(pTriangle2); 
        array1.push_back(pTriangle3);

        //Call the Faces Portion utility
        double volume = VolumeInsideVoxelQEF::FacesPortionQEFApproximation(*pVoxel,array1);
        const double ExpectedVolume = 0.11458; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01);      
    }

    KRATOS_TEST_CASE_IN_SUITE(QEFApproximations2, KratosCoreFastSuite) {
        //A voxel crossed by a straight plane with only 2 nodes inside the volume (ood case approximation)
        std::vector<double> distances{1, 1, -1, -1, -1, -1, -1, -1};   
        GeometryPtrType pVoxel = QEFVolumeGenerateHexahedra3D8(distances);

        //Generate the intersecting triangles
        std::vector<std::vector<double>> triangle1{{-1,0.0,-0.95},{-0.95,0.0,-1.05},{-1.05,0.0,-1.05}}; 
        std::vector<std::vector<double>> triangle2{{-1,-0.95,0},{-0.95,-1.05,0},{-1.05,-1.05,0}}; 
        std::vector<std::vector<double>> triangle3{{1,-0.95,0},{1.05,-1.05,0},{0.95,-1.05,0}}; 
        std::vector<std::vector<double>> triangle4{{1,0.0,-0.95},{1.05,0.0,-1.05},{0.95,0.0,-1.05}}; 
        
        GeometryPtrType pTriangle1 = QEFVolumeGenerateTriangle3D3(triangle1);
        GeometryPtrType pTriangle2 = QEFVolumeGenerateTriangle3D3(triangle2);
        GeometryPtrType pTriangle3 = QEFVolumeGenerateTriangle3D3(triangle3);
        GeometryPtrType pTriangle4 = QEFVolumeGenerateTriangle3D3(triangle4);

        GeometryArrayType array1;
        array1.push_back(pTriangle1); 
        array1.push_back(pTriangle2);
        array1.push_back(pTriangle3);
        array1.push_back(pTriangle4);

        //Call the volume utility
        double volume = VolumeInsideVoxelQEF::FacesPortionQEFApproximation(*pVoxel,array1);
        double ExpectedVolume = 0.25; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01); 
        volume = VolumeInsideVoxelQEF::GeometricalQEFApproximation(*pVoxel,array1);
        ExpectedVolume = 0.2083; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01); 
        /*in this case the volume returned by the method approximation exactly the same to the real volume of the 
        test case, since this is a perfect choosen case*/
    }

    KRATOS_TEST_CASE_IN_SUITE(QEFApproximations3, KratosCoreFastSuite) {
        //A voxel crossed by a straight plane with only 2 nodes inside th volume 
        std::vector<double> distances{1, 1, -1, -1, -1, -1, -1, -1};   
        GeometryPtrType pVoxel = QEFVolumeGenerateHexahedra3D8(distances);

         //Generate the intersecting triangles
        std::vector<std::vector<double>> triangle1{{-1,0.9,-0.95},{-0.95,0.9,-1.05},{-1.05,0.9,-1.05}}; 
        std::vector<std::vector<double>> triangle2{{-1,-0.95,0},{-0.95,-1.05,0},{-1.05,-1.05,0}}; 
        std::vector<std::vector<double>> triangle3{{1,-0.95,0},{1.05,-1.05,0},{0.95,-1.05,0}}; 
        std::vector<std::vector<double>> triangle4{{1,0.9,-0.95},{1.05,0.9,-1.05},{0.95,0.9,-1.05}}; 
        
        GeometryPtrType pTriangle1 = QEFVolumeGenerateTriangle3D3(triangle1);
        GeometryPtrType pTriangle2 = QEFVolumeGenerateTriangle3D3(triangle2);
        GeometryPtrType pTriangle3 = QEFVolumeGenerateTriangle3D3(triangle3);
        GeometryPtrType pTriangle4 = QEFVolumeGenerateTriangle3D3(triangle4);

        GeometryArrayType array1;
        array1.push_back(pTriangle1); 
        array1.push_back(pTriangle2);
        array1.push_back(pTriangle3);
        array1.push_back(pTriangle4);

        //Call the volume utility
        double volume = VolumeInsideVoxelQEF::FacesPortionQEFApproximation(*pVoxel,array1);
        double ExpectedVolume = 0.4; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01); 
        volume = VolumeInsideVoxelQEF::GeometricalQEFApproximation(*pVoxel,array1);
        ExpectedVolume = 0.4333; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01); 
    }

    KRATOS_TEST_CASE_IN_SUITE(QEFApproximations4, KratosCoreFastSuite) {
        //A voxel crossed by a straight plane with 4 nodes inside the volume
        std::vector<double> distances{1, 1, -1, -1, 1, 1, -1, -1};   
        GeometryPtrType pVoxel = QEFVolumeGenerateHexahedra3D8(distances);

        //Generate the intersecting triangles
        std::vector<std::vector<double>> triangle1{{1,0.5,1.05},{1.05,0.5,0.95},{0.95,0.5,0.95}};  
        std::vector<std::vector<double>> triangle2{{1,0.5,-0.95},{1.05,0.5,-1.05},{0.95,0.5,-1.05}}; 
        std::vector<std::vector<double>> triangle3{{-1,0.5,-0.95},{-0.95,0.5,-1.05},{-1.05,0.5,-1.05}}; 
        std::vector<std::vector<double>> triangle4{{-1,0.5,1.05},{-0.95,0.5,0.95},{-1.05,0.5,0.95}}; 
        
        GeometryPtrType pTriangle1 = QEFVolumeGenerateTriangle3D3(triangle1);
        GeometryPtrType pTriangle2 = QEFVolumeGenerateTriangle3D3(triangle2);
        GeometryPtrType pTriangle3 = QEFVolumeGenerateTriangle3D3(triangle3);
        GeometryPtrType pTriangle4 = QEFVolumeGenerateTriangle3D3(triangle4);

        GeometryArrayType array1;
        array1.push_back(pTriangle1); 
        array1.push_back(pTriangle2);
        array1.push_back(pTriangle3);
        array1.push_back(pTriangle4);

        //Call the volume utility
        double volume = VolumeInsideVoxelQEF::FacesPortionQEFApproximation(*pVoxel,array1);
        double ExpectedVolume = 0.6667; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01); 
        volume = VolumeInsideVoxelQEF::GeometricalQEFApproximation(*pVoxel,array1);
        ExpectedVolume = 0.75; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01); 
        /*in this case the volume returned by the method approximation is quite close to the real volume of the 
        test case, since we would expect a real volume of this case to be circa 0.75 */
    } 

    KRATOS_TEST_CASE_IN_SUITE(QEFApproximations5, KratosCoreFastSuite) {
        //A voxel crossed by two straight planes (enclosing half of the volume) with no nodes inside the volume
        std::vector<double> distances{-1, -1, -1, -1, -1, -1, -1, -1};   
        GeometryPtrType pVoxel = QEFVolumeGenerateHexahedra3D8(distances);

        //Generate the intersecting triangles
        std::vector<std::vector<double>> triangle1{{1,0.5,1.05},{1.05,0.5,0.95},{0.95,0.5,0.95}};  
        std::vector<std::vector<double>> triangle2{{1,0.5,-0.95},{1.05,0.5,-1.05},{0.95,0.5,-1.05}}; 
        std::vector<std::vector<double>> triangle3{{-1,0.5,-0.95},{-0.95,0.5,-1.05},{-1.05,0.5,-1.05}}; 
        std::vector<std::vector<double>> triangle4{{-1,0.5,1.05},{-0.95,0.5,0.95},{-1.05,0.5,0.95}}; 
        std::vector<std::vector<double>> triangle5{{1,-0.5,1.05},{1.05,-0.5,0.95},{0.95,-0.5,0.95}};  
        std::vector<std::vector<double>> triangle6{{1,-0.5,-0.95},{1.05,-0.5,-1.05},{0.95,-0.5,-1.05}}; 
        std::vector<std::vector<double>> triangle7{{-1,-0.5,-0.95},{-0.95,-0.5,-1.05},{-1.05,-0.5,-1.05}}; 
        std::vector<std::vector<double>> triangle8{{-1,-0.5,1.05},{-0.95,-0.5,0.95},{-1.05,-0.5,0.95}}; 
        
        GeometryPtrType pTriangle1 = QEFVolumeGenerateTriangle3D3(triangle1);
        GeometryPtrType pTriangle2 = QEFVolumeGenerateTriangle3D3(triangle2);
        GeometryPtrType pTriangle3 = QEFVolumeGenerateTriangle3D3(triangle3);
        GeometryPtrType pTriangle4 = QEFVolumeGenerateTriangle3D3(triangle4);
        GeometryPtrType pTriangle5 = QEFVolumeGenerateTriangle3D3(triangle5);
        GeometryPtrType pTriangle6 = QEFVolumeGenerateTriangle3D3(triangle6);
        GeometryPtrType pTriangle7 = QEFVolumeGenerateTriangle3D3(triangle7);
        GeometryPtrType pTriangle8 = QEFVolumeGenerateTriangle3D3(triangle8);

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
        double volume = VolumeInsideVoxelQEF::FacesPortionQEFApproximation(*pVoxel,array1);
        double ExpectedVolume = 0.1666; //no nodes inside
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01); 
        volume = VolumeInsideVoxelQEF::GeometricalQEFApproximation(*pVoxel,array1);
        ExpectedVolume = 0.3333; 
        KRATOS_CHECK_NEAR(volume, ExpectedVolume, 0.01); 
        /*in this case the volume returned by the method approximation is not close to the real volume of the 
        test case, since we would expect a real volume circa 0.5. Anyways, it is better that the previous method */
    } 

} //namespace testing
} //namespace kratos