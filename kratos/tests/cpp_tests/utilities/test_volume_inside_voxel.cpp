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
//#include "geometries/node.h"
//#include "geometries/hexahedra_3d_8.h"
#include "geometries/triangle_3d_3.h"
#include "testing/testing.h"
#include "utilities/intersection_utilities.h"

namespace Kratos {
namespace Testing {
/*
    //Define our testing variables
    typedef Node<3> NodeType;
    typedef Node<3>::Pointer NodePtrType;
    typedef Hexahedra3D8<NodeType> HexaGeometryType;
    typedef HexaGeometryType::Pointer HexaGeometryPtrType;

    NodeType n1();
    NodePtrType p_n1= Kratos::make_shared<NodeType>(1,0.0,0.0,0.0);

    p_n1->AddDof(DISTANCE);*/

    KRATOS_TEST_CASE_IN_SUITE(VolumeInsideVoxelOctal,KratosCoreFastSuite) 
    {

    }

}  // namespace Testing.
}  // namespace Kratos.
