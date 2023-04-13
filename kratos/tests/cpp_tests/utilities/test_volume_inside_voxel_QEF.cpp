//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/license.txt
//
//  Main authors:    Ariadna Cortes
//

// Project includes
#include "containers/model.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/triangle_3d_3.h"
#include "includes/checks.h"
#include "testing/testing.h"
#include "utilities/Volume_inside_voxel_QEF.h"

namespace Kratos {
namespace Testing {

    typedef Geometry<Node<3>> GeometryType;
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

    /******************************************************************************************************
     ******************************************************************************************************/

    //Basic cases tests (all nodes inside/ouside)
    KRATOS_TEST_CASE_IN_SUITE(QEFApproximationsFull, KratosCoreFastSuite) {
        //A voxel with all nodes inside the Volume
        std::vector<double> Distances{1, 1, 1, 1, 1, 1, 1, 1};   
        GeometryPtrType pVoxel = QEFVolumeGenerateHexahedra3D8(Distances);

        GeometryArrayType Array1;

        //Call the Volume utility
        double Volume = VolumeInsideVoxelQEF::SimpleNodesQEFApproximation(*pVoxel,Array1);
        double ExpectedVolume = 1.0; 
        KRATOS_CHECK_NEAR(Volume, ExpectedVolume, 0.01); 

        Volume = VolumeInsideVoxelQEF::VolumeQEFApproximation(*pVoxel,Array1);
        KRATOS_CHECK_NEAR(Volume, ExpectedVolume, 0.01); 
    }

    KRATOS_TEST_CASE_IN_SUITE(QEFApproximationsEmpty, KratosCoreFastSuite) {
        //A voxel with all nodes outside the Volume
        std::vector<double> distances{-1, -1, -1, -1, -1, -1, -1, -1};   
        GeometryPtrType pVoxel = QEFVolumeGenerateHexahedra3D8(distances);

        GeometryArrayType Array1;

        //Call the Volume utility
        double Volume = VolumeInsideVoxelQEF::SimpleNodesQEFApproximation(*pVoxel,Array1);
        double ExpectedVolume = 0.0; 
        KRATOS_CHECK_NEAR(Volume, ExpectedVolume, 0.01); 

        Volume = VolumeInsideVoxelQEF::VolumeQEFApproximation(*pVoxel,Array1);
        KRATOS_CHECK_NEAR(Volume, ExpectedVolume, 0.01); 
    }

    //A couple of test regarding the real efficiency and limitations of this method

    KRATOS_TEST_CASE_IN_SUITE(QEFApproximations1, KratosCoreFastSuite) {
        //voxel crossed by a plane with only 1 node inside the Volume
        std::vector<double> distances{1, -1, -1, -1, -1, -1, -1, -1};   
        GeometryPtrType pVoxel = QEFVolumeGenerateHexahedra3D8(distances);

        //Generate the intersecting triangles
        std::vector<std::vector<double>> triangle1{{-1,0.5,-0.95},{-0.95,0.5,-1.05},{-1.05,0.5,-1.05}}; 
        std::vector<std::vector<double>> triangle2{{-1,-0.95,0},{-0.95,-1.05,0},{-1.05,-1.05,0}}; 
        std::vector<std::vector<double>> triangle3{{-0.5,-0.95,-1},{-0.5,-1.05,-0.95},{-0.5,-1.05,-1.05}}; 
        
        GeometryPtrType pTriangle1 = QEFVolumeGenerateTriangle3D3(triangle1);
        GeometryPtrType pTriangle2 = QEFVolumeGenerateTriangle3D3(triangle2);
        GeometryPtrType pTriangle3 = QEFVolumeGenerateTriangle3D3(triangle3);

        GeometryArrayType Array1;
        Array1.push_back(pTriangle1); 
        Array1.push_back(pTriangle2); 
        Array1.push_back(pTriangle3);


        //Call the Volume utility
        double Volume = VolumeInsideVoxelQEF::SimpleNodesQEFApproximation(*pVoxel,Array1);
        double ExpectedVolume = 0.125; 
        KRATOS_CHECK_NEAR(Volume, ExpectedVolume, 0.01);

        Volume = VolumeInsideVoxelQEF::VolumeQEFApproximation(*pVoxel,Array1);
        ExpectedVolume = 0.046875; 
        KRATOS_CHECK_NEAR(Volume, ExpectedVolume, 0.01); 
        //Note: the real expected Volume assumed in this case was circa 0.09375     
    }

    KRATOS_TEST_CASE_IN_SUITE(QEFApproximations2, KratosCoreFastSuite) {
        //A voxel crossed by a straight plane with only 2 nodes inside the Volume
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

        GeometryArrayType Array1;
        Array1.push_back(pTriangle1); 
        Array1.push_back(pTriangle2);
        Array1.push_back(pTriangle3);
        Array1.push_back(pTriangle4);

        //Call the Volume utility
        double Volume = VolumeInsideVoxelQEF::SimpleNodesQEFApproximation(*pVoxel,Array1);
        double ExpectedVolume = 0.25; 
        KRATOS_CHECK_NEAR(Volume, ExpectedVolume, 0.01); 

        Volume = VolumeInsideVoxelQEF::VolumeQEFApproximation(*pVoxel,Array1);
        ExpectedVolume = 0.2083; 
        KRATOS_CHECK_NEAR(Volume, ExpectedVolume, 0.01); 
        /*in this case the Volume returned by the first two methods is exactly the same as the real Volume of the 
        test case (0.25), since this is a perfect choosen case*/
    }

    KRATOS_TEST_CASE_IN_SUITE(QEFApproximations3, KratosCoreFastSuite) {
        //A voxel crossed by a straight plane with only 2 nodes inside th Volume (extreme case)
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

        GeometryArrayType Array1;
        Array1.push_back(pTriangle1); 
        Array1.push_back(pTriangle2);
        Array1.push_back(pTriangle3);
        Array1.push_back(pTriangle4);

        //Call the Volume utility
        double Volume = VolumeInsideVoxelQEF::SimpleNodesQEFApproximation(*pVoxel,Array1);
        double ExpectedVolume = 3.0/8; 
        KRATOS_CHECK_NEAR(Volume, ExpectedVolume, 0.01);

        Volume = VolumeInsideVoxelQEF::VolumeQEFApproximation(*pVoxel,Array1);
        ExpectedVolume = 0.622; 
        KRATOS_CHECK_NEAR(Volume, ExpectedVolume, 0.01); 
        //Note: the real expected Volume assumed in this case was circa 0.74     
    }

    KRATOS_TEST_CASE_IN_SUITE(QEFApproximations4, KratosCoreFastSuite) {
        //A voxel crossed by a straight plane with 4 nodes inside the Volume
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

        GeometryArrayType Array1;
        Array1.push_back(pTriangle1); 
        Array1.push_back(pTriangle2);
        Array1.push_back(pTriangle3);
        Array1.push_back(pTriangle4);

        //Call the Volume utility
        double Volume = VolumeInsideVoxelQEF::SimpleNodesQEFApproximation(*pVoxel,Array1);
        double ExpectedVolume = 7.0/12; //0.5833
        KRATOS_CHECK_NEAR(Volume, ExpectedVolume, 0.01); 

        Volume = VolumeInsideVoxelQEF::VolumeQEFApproximation(*pVoxel,Array1);
        ExpectedVolume = 0.75; 
        KRATOS_CHECK_NEAR(Volume, ExpectedVolume, 0.01); 
        /*in this case the Volume returned by third method is exactly the expected real Volume of the 
        test case, since we would expect a real Volume of this case to be circa 0.75 */
    } 

    KRATOS_TEST_CASE_IN_SUITE(QEFApproximations5, KratosCoreFastSuite) {
        //A voxel crossed by a straight plane with 4 nodes inside the Volume (extreme case)
        std::vector<double> distances{1, 1, -1, -1, 1, 1, -1, -1};   
        GeometryPtrType pVoxel = QEFVolumeGenerateHexahedra3D8(distances);

        //Generate the intersecting triangles
        std::vector<std::vector<double>> triangle1{{1,0.98,1.05},{1.05,0.98,0.95},{0.95,0.98,0.95}};  
        std::vector<std::vector<double>> triangle2{{1,0.98,-0.95},{1.05,0.98,-1.05},{0.95,0.98,-1.05}}; 
        std::vector<std::vector<double>> triangle3{{-1,0.98,-0.95},{-0.95,0.98,-1.05},{-1.05,0.98,-1.05}}; 
        std::vector<std::vector<double>> triangle4{{-1,0.98,1.05},{-0.95,0.98,0.95},{-1.05,0.98,0.95}}; 
        
        GeometryPtrType pTriangle1 = QEFVolumeGenerateTriangle3D3(triangle1);
        GeometryPtrType pTriangle2 = QEFVolumeGenerateTriangle3D3(triangle2);
        GeometryPtrType pTriangle3 = QEFVolumeGenerateTriangle3D3(triangle3);
        GeometryPtrType pTriangle4 = QEFVolumeGenerateTriangle3D3(triangle4);

        GeometryArrayType Array1;
        Array1.push_back(pTriangle1); 
        Array1.push_back(pTriangle2);
        Array1.push_back(pTriangle3);
        Array1.push_back(pTriangle4);

        //Call the Volume utility
        double Volume = VolumeInsideVoxelQEF::SimpleNodesQEFApproximation(*pVoxel,Array1);
        double ExpectedVolume = 0.6633; 
        KRATOS_CHECK_NEAR(Volume, ExpectedVolume, 0.01); 

        Volume = VolumeInsideVoxelQEF::VolumeQEFApproximation(*pVoxel,Array1);
        ExpectedVolume = 0.99; 
        KRATOS_CHECK_NEAR(Volume, ExpectedVolume, 0.01); 
        /*in this case the Volume returned by third method is exactly the expected real Volume of the 
        test case, since we would expect a real Volume of this case to be circa 0.99 */
    } 

    KRATOS_TEST_CASE_IN_SUITE(QEFApproximations6, KratosCoreFastSuite) {
        //A voxel crossed by two straight planes (enclosing half of the Volume) with no nodes inside the Volume
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

        GeometryArrayType Array1;
        Array1.push_back(pTriangle1); 
        Array1.push_back(pTriangle2);
        Array1.push_back(pTriangle3);
        Array1.push_back(pTriangle4);
        Array1.push_back(pTriangle5); 
        Array1.push_back(pTriangle6);
        Array1.push_back(pTriangle7);
        Array1.push_back(pTriangle8);

        //Call the Volume utility
        double Volume = VolumeInsideVoxelQEF::SimpleNodesQEFApproximation(*pVoxel,Array1);
        double ExpectedVolume = 0; //no nodes inside
        KRATOS_CHECK_NEAR(Volume, ExpectedVolume, 0.01); 

        Volume = VolumeInsideVoxelQEF::VolumeQEFApproximation(*pVoxel,Array1);
        ExpectedVolume = 0.3333; 
        KRATOS_CHECK_NEAR(Volume, ExpectedVolume, 0.01); 
        /*in this case the Volume returned by the method approximation is not close to the real Volume of the 
        test case, since we would expect a real Volume circa 0.5. Anyways, it is better that the previous methods */
    } 

    //A couple of test showing most common use cases (exact approximation)

    KRATOS_TEST_CASE_IN_SUITE(QEFApproximationsUseCase1, KratosCoreFastSuite) {
        //voxel crossed by a plane with only 1 node inside the Volume
        std::vector<double> distances{1, -1, -1, -1, -1, -1, -1, -1};   
        GeometryPtrType pVoxel = QEFVolumeGenerateHexahedra3D8(distances);

        //Generate the intersecting triangles
        std::vector<std::vector<double>> triangle1{{-1,-0.05,-0.95},{-0.95,0,-1.05},{-1.05,0.1,-1.05}}; 
        std::vector<std::vector<double>> triangle2{{-1,-0.95,-0.05},{-0.95,-1.05,0},{-1.05,-1.05,0.1}}; 
        std::vector<std::vector<double>> triangle3{{-0.05,-0.95,-1},{0,-1.05,-0.95},{0.1,-1.05,-1.05}}; 
        
        GeometryPtrType pTriangle1 = QEFVolumeGenerateTriangle3D3(triangle1);
        GeometryPtrType pTriangle2 = QEFVolumeGenerateTriangle3D3(triangle2);
        GeometryPtrType pTriangle3 = QEFVolumeGenerateTriangle3D3(triangle3);

        GeometryArrayType Array1;
        Array1.push_back(pTriangle1); 
        Array1.push_back(pTriangle2); 
        Array1.push_back(pTriangle3);

        static array_1d<double,3> QEF = QuadraticErrorFunction::QuadraticErrorFunctionPoint(*pVoxel, Array1);
        KRATOS_CHECK_NEAR(QEF[0],-0.6667,0.01);
        KRATOS_CHECK_NEAR(QEF[1],-0.6667,0.01);
        KRATOS_CHECK_NEAR(QEF[2],-0.6667,0.01);

        double Volume = VolumeInsideVoxelQEF::VolumeQEFApproximation(*pVoxel,Array1);
        double ExpectedVolume = 0.02083; //no nodes inside
        KRATOS_CHECK_NEAR(Volume, ExpectedVolume, 0.001);
        //Exact expected result  
    }

    KRATOS_TEST_CASE_IN_SUITE(QEFApproximationsUseCase2, KratosCoreFastSuite) {
        //A voxel crossed by a straight plane with only 2 nodes inside the Volume
        std::vector<double> distances{1, 1, -1, -1, -1, -1, -1, -1};   
        GeometryPtrType pVoxel = QEFVolumeGenerateHexahedra3D8(distances);

        //Generate the intersecting triangles
        std::vector<std::vector<double>> triangle1{{-1,-0.05,-0.95},{-0.95,0.05,-1.05},{-1.05,0.05,-1.05}}; 
        std::vector<std::vector<double>> triangle2{{-1,-0.95,-0.05},{-0.95,-1.05,0.05},{-1.05,-1.05,0.05}}; 
        std::vector<std::vector<double>> triangle3{{1,-0.95,-0.05},{1.05,-1.05,0.05},{0.95,-1.05,0.05}}; 
        std::vector<std::vector<double>> triangle4{{1,-0.05,-0.95},{1.05,0.05,-1.05},{0.95,0.05,-1.05}}; 
        
        GeometryPtrType pTriangle1 = QEFVolumeGenerateTriangle3D3(triangle1);
        GeometryPtrType pTriangle2 = QEFVolumeGenerateTriangle3D3(triangle2);
        GeometryPtrType pTriangle3 = QEFVolumeGenerateTriangle3D3(triangle3);
        GeometryPtrType pTriangle4 = QEFVolumeGenerateTriangle3D3(triangle4);

        GeometryArrayType Array1;
        Array1.push_back(pTriangle1); 
        Array1.push_back(pTriangle2);
        Array1.push_back(pTriangle3);
        Array1.push_back(pTriangle4);

        double Volume = VolumeInsideVoxelQEF::VolumeQEFApproximation(*pVoxel,Array1);
        double ExpectedVolume = 0.125; 
        KRATOS_CHECK_NEAR(Volume, ExpectedVolume, 0.001); 
        //Exact expected result
    }

    KRATOS_TEST_CASE_IN_SUITE(QEFApproximationsUseCase3, KratosCoreFastSuite) {
        //A voxel crossed by a shifted plane with 4 nodes inside the Volume
        std::vector<double> distances{1, 1, 1, 1, -1, -1, -1, -1};   
        GeometryPtrType pVoxel = QEFVolumeGenerateHexahedra3D8(distances);

        //Generate the intersecting triangles
        std::vector<std::vector<double>> triangle1{{-1,-0.95+2,-0.025-0.5},{-0.95,-1.05+2,0.025-0.5},{-1.05,-1.05+2,0.025-0.5}}; 
        std::vector<std::vector<double>> triangle2{{-1,-0.95,-0.025+0.5},{-0.95,-1.05,0.025+0.5},{-1.05,-1.05,0.025+0.5}}; 
        std::vector<std::vector<double>> triangle3{{1,-0.95,-0.025+0.5},{1.05,-1.05,0.025+0.5},{0.95,-1.05,0.025+0.5}}; 
        std::vector<std::vector<double>> triangle4{{1,-0.95+2,-0.025-0.5},{1.05,-1.05+2,0.025-0.5},{0.95,-1.05+2,0.025-0.5}}; 
        
        GeometryPtrType pTriangle1 = QEFVolumeGenerateTriangle3D3(triangle1);
        GeometryPtrType pTriangle2 = QEFVolumeGenerateTriangle3D3(triangle2);
        GeometryPtrType pTriangle3 = QEFVolumeGenerateTriangle3D3(triangle3);
        GeometryPtrType pTriangle4 = QEFVolumeGenerateTriangle3D3(triangle4);

        GeometryArrayType Array1;
        Array1.push_back(pTriangle1); 
        Array1.push_back(pTriangle2);
        Array1.push_back(pTriangle3);
        Array1.push_back(pTriangle4);

        double Volume = VolumeInsideVoxelQEF::VolumeQEFApproximation(*pVoxel,Array1);
        double ExpectedVolume = 0.5; 
        KRATOS_CHECK_NEAR(Volume, ExpectedVolume, 0.001);  
        //Exact expected result
    }

    KRATOS_TEST_CASE_IN_SUITE(QEFApproximationsUseCase4, KratosCoreFastSuite) {
        //voxel crossed by a plane with only 1 node inside the Volume
        std::vector<double> distances{-1, 1, 1, 1, 1, 1, 1, 1};   
        GeometryPtrType pVoxel = QEFVolumeGenerateHexahedra3D8(distances);

        //Generate the intersecting triangles
        std::vector<std::vector<double>> triangle1{{-1,-0.05,-0.95},{-0.95,0,-1.05},{-1.05,0.1,-1.05}}; 
        std::vector<std::vector<double>> triangle2{{-1,-0.95,-0.05},{-0.95,-1.05,0},{-1.05,-1.05,0.1}}; 
        std::vector<std::vector<double>> triangle3{{-0.05,-0.95,-1},{0,-1.05,-0.95},{0.1,-1.05,-1.05}}; 
        
        GeometryPtrType pTriangle1 = QEFVolumeGenerateTriangle3D3(triangle1);
        GeometryPtrType pTriangle2 = QEFVolumeGenerateTriangle3D3(triangle2);
        GeometryPtrType pTriangle3 = QEFVolumeGenerateTriangle3D3(triangle3);

        GeometryArrayType Array1;
        Array1.push_back(pTriangle1); 
        Array1.push_back(pTriangle2); 
        Array1.push_back(pTriangle3);

        double Volume = VolumeInsideVoxelQEF::VolumeQEFApproximation(*pVoxel,Array1);
        double ExpectedVolume = 1 - 0.02083; 
        KRATOS_CHECK_NEAR(Volume, ExpectedVolume, 0.001);
        //Exact expected result
    } 
    
    KRATOS_TEST_CASE_IN_SUITE(QEFApproximationsUseCase5, KratosCoreFastSuite) {
        //voxel crossed by a plane with only 1 node inside the Volume (intermediate case)
        std::vector<double> distances{1, -1, -1, -1, -1, -1, -1, -1};   
        GeometryPtrType pVoxel = QEFVolumeGenerateHexahedra3D8(distances);

        //Generate the intersecting triangles
        std::vector<std::vector<double>> triangle1{{-1,-0.05+0.2,-0.95},{-0.95,0+0.2,-1.05},{-1.05,0.1+0.2,-1.05}}; 
        std::vector<std::vector<double>> triangle2{{-1,-0.95,-0.05+0.2},{-0.95,-1.05,0+0.2},{-1.05,-1.05,0.1+0.2}}; 
        std::vector<std::vector<double>> triangle3{{-0.05+0.2,-0.95,-1},{0+0.2,-1.05,-0.95},{0.1+0.2,-1.05,-1.05}}; 
        
        GeometryPtrType pTriangle1 = QEFVolumeGenerateTriangle3D3(triangle1);
        GeometryPtrType pTriangle2 = QEFVolumeGenerateTriangle3D3(triangle2);
        GeometryPtrType pTriangle3 = QEFVolumeGenerateTriangle3D3(triangle3);

        GeometryArrayType Array1;
        Array1.push_back(pTriangle1); 
        Array1.push_back(pTriangle2); 
        Array1.push_back(pTriangle3);

        double Volume = VolumeInsideVoxelQEF::VolumeQEFApproximation(*pVoxel,Array1);
        double ExpectedVolume = 0.036; 
        KRATOS_CHECK_NEAR(Volume, ExpectedVolume, 0.001);
        //The expected result for this case would be circa 0.036
    }

    KRATOS_TEST_CASE_IN_SUITE(QEFApproximationsUseCase6, KratosCoreFastSuite) {
        //voxel crossed by a plane with only 1 node inside the Volume (extreme case)
        std::vector<double> distances{1, -1, -1, -1, -1, -1, -1, -1};   
        GeometryPtrType pVoxel = QEFVolumeGenerateHexahedra3D8(distances);

        //Generate the intersecting triangles
        std::vector<std::vector<double>> triangle1{{-1,-0.05+0.999,-0.95},{-0.95,0+0.999,-1.05},{-1.05,0.1+0.999,-1.05}}; 
        std::vector<std::vector<double>> triangle2{{-1,-0.95,-0.05+0.999},{-0.95,-1.05,0+0.999},{-1.05,-1.05,0.1+0.999}}; 
        std::vector<std::vector<double>> triangle3{{-0.05+0.999,-0.95,-1},{0+0.999,-1.05,-0.95},{0.1+0.999,-1.05,-1.05}}; 
        
        GeometryPtrType pTriangle1 = QEFVolumeGenerateTriangle3D3(triangle1);
        GeometryPtrType pTriangle2 = QEFVolumeGenerateTriangle3D3(triangle2);
        GeometryPtrType pTriangle3 = QEFVolumeGenerateTriangle3D3(triangle3);

        GeometryArrayType Array1;
        Array1.push_back(pTriangle1); 
        Array1.push_back(pTriangle2); 
        Array1.push_back(pTriangle3);

        static array_1d<double,3> QEF = QuadraticErrorFunction::QuadraticErrorFunctionPoint(*pVoxel, Array1);
        KRATOS_CHECK_NEAR(QEF[0],-0.333,0.01);
        KRATOS_CHECK_NEAR(QEF[1],-0.333,0.01);
        KRATOS_CHECK_NEAR(QEF[2],-0.333,0.01);

        double Volume = VolumeInsideVoxelQEF::VolumeQEFApproximation(*pVoxel,Array1);
        double ExpectedVolume = 0.1667; 
        KRATOS_CHECK_NEAR(Volume, ExpectedVolume, 0.001);
        //The expected result for this case would be 1/6 = 0.1667
    }

} //namespace testing
} //namespace kratos