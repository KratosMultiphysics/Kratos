//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// Project includes
#include "containers/model.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/triangle_3d_3.h"
#include "testing/testing.h"
#include "utilities/intersection_utilities.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(VolumeInsideVoxelOctal,KratosCoreFastSuite) {

    typedef Node<3> NodeType;
    typedef Node<3>::Pointer NodePtrType;
    typedef Hexahedra3D8<NodeType> HexaGeometryType;
    typedef HexaGeometryType::Pointer HexaGeometryPtrType;
    typedef Triangle3D3<NodeType> TriGeometryType;
    typedef TriGeometryType::Pointer TriGeometryPtrType;
    typedef TriGeometryType::GeometriesArrayType TriGeometryArray;
    typedef HexaGeometryType::PointsArrayType HexaPointsArray;

    //Define our testing variables

    //Generate the HEXAHEDRA3D8
    Model my_model;
    ModelPart &voxel = my_model.CreateModelPart("Voxel");  
    voxel.AddNodalSolutionStepVariable(DISTANCE);  

    voxel.CreateNewNode(1, -0.5, -0.5, -0.5);
    voxel.CreateNewNode(2,  0.5, -0.5, -0.5);
    voxel.CreateNewNode(3, 0.5,  0.5, -0.5);
    voxel.CreateNewNode(4, -0.5,  0.5, -0.5);
    voxel.CreateNewNode(5, -0.5, -0.5,  0.5);
    voxel.CreateNewNode(6, 0.5, -0.5,  0.5);
    voxel.CreateNewNode(7, 0.5,  0.5,  0.5);
    voxel.CreateNewNode(8, -0.5,  0.5,  0.5); 
    Properties::Pointer p_properties_0(new Properties(0));
    voxel.CreateNewElement("Voxel", 1, {1, 2, 3, 4, 5, 6, 7, 8}, p_properties_0); 
    //HexaGeometryPtrType pVoxel =  voxel.pGetGeometry(1);  

    //Generate some TRIANGLE3D3 and array of them 
    ModelPart& triangles = my_model.CreateModelPart("Triangles");
    
    triangles.CreateNewNode(1, 0.0, 0.0, 0.0);
    triangles.CreateNewNode(2, 0.0, 0.0, 0.0);
    triangles.CreateNewNode(3, 0.0, 0.0, 0.0);
    triangles.CreateNewNode(4, 0.0, 0.0, 0.0);
    triangles.CreateNewNode(5, 0.0, 0.0, 0.0);
    triangles.CreateNewNode(6, 0.0, 0.0, 0.0);
    triangles.CreateNewNode(7, 0.0, 0.0, 0.0);
    triangles.CreateNewNode(8, 0.0, 0.0, 0.0);
    triangles.CreateNewNode(9, 0.0, 0.0, 0.0); 

    Properties::Pointer p_properties_1(new Properties(0));  /*
    TriGeometryPtrType pTriangle1 = triangles.CreateNewElement("Triangle1", 1, {1, 2, 3}, p_properties_1);
    TriGeometryPtrType pTriangle2 = triangles.CreateNewElement("Triangle2", 2, {4, 5, 6}, p_properties_1);
    TriGeometryPtrType pTriangle3 = triangles.CreateNewElement("Triangle3", 3, {7, 8, 9}, p_properties_1);

    //Add the distances for this test
    HexaPointsArray nodes = pVoxel->Points();
    const array_1d<double,1> distances(8);              //Hardcoded to 8 point-geometry
    distances = {1, -1, -1, -1, -1, -1, -1, -1,};       
    
    for (int i = 0; i < 8; i++) {
        nodes[i].FastGetSolutionStepValue(DISTANCE) = distances[i];
    }

    //Select the triangles for this test
    TriGeometryArray triArray;
    triArray.push_back(pTriangle1);
    triArray.push_back(pTriangle2);
    triArray.push_back(pTriangle3);

    //Call the volume utility
    const double volume = VolumeInsideVoxelUtility::OctalApproximation(pVoxel, triArray);

    //Expected output of the function
    const double real_volume = 0.125;
    
    KRATOS_CHECK_EQUAL(volume, real_volume);
    //KRATOS_CHECK_NEAR(volume, real_volume, 0.1);
    */
}
}  // namespace Testing.
}  // namespace Kratos.
