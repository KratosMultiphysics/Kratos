//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "utilities/tessellation_utilities/delaunator_utilities.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(CreateTriangleMeshFromNodes, KratosCoreFastSuite)
{
    Model this_model;
    ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);

    // First we create the nodes
    r_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
    r_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
    r_model_part.CreateNewNode(3, 1.0 , 1.0 , 0.0);
    r_model_part.CreateNewNode(4, 0.0 , 1.0 , 0.0);
    r_model_part.CreateNewNode(5, 2.0 , 0.0 , 0.0);
    r_model_part.CreateNewNode(6, 2.0 , 1.0 , 0.0);

    DelaunatorUtilities::CreateTriangleMeshFromNodes(r_model_part);

    KRATOS_EXPECT_TRUE(r_model_part.Elements().size() == 4);
    
    // Check area
    double area = 0.0;
    for (auto& r_elem : r_model_part.Elements()) {
        area += r_elem.GetGeometry().Area();
    }
    
    KRATOS_EXPECT_RELATIVE_NEAR(area, 2.0, 1.0e-12);
}

}   // namespace Testing
}  // namespace Kratos.