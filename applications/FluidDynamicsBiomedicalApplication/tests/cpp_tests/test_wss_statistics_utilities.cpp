//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:   Ruben Zorrilla
//                  Jordi Cotela
//                  Eduardo Soudah
//

#include "testing/testing.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/cfd_variables.h"

#include "custom_utilities/wss_statistics_utilities.h"

namespace Kratos {
namespace Testing {

void TetrahedraModelPartForWSSTests(ModelPart& rModelPart) {

    rModelPart.AddNodalSolutionStepVariable(REACTION);
    rModelPart.SetBufferSize(3);
    Properties::Pointer p_properties = rModelPart.CreateNewProperties(0);

    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
    rModelPart.CreateNewNode(4, 0.0, 0.0, 1.0);
    std::vector<ModelPart::IndexType> element_nodes{1, 2, 3, 4};
    rModelPart.CreateNewElement("Element3D4N", 1, element_nodes, p_properties);

    // Nodal data
    Element& r_element = *(rModelPart.ElementsBegin());
    Geometry< Node<3> >& r_geometry = r_element.GetGeometry();

    constexpr double reaction = 1.0;

    for (unsigned int i = 0; i < 4; i++) {
        Node<3>& r_node = r_geometry[i];
        r_node.FastGetSolutionStepValue(REACTION_X) = reaction;
        r_node.FastGetSolutionStepValue(REACTION_Y) = reaction;
        r_node.FastGetSolutionStepValue(REACTION_Z) = reaction;
    }
}

KRATOS_TEST_CASE_IN_SUITE(WSSUtilities3DWSS, FluidDynamicsApplicationFastSuite) {
    Model model;
    ModelPart& ModelPart = model.CreateModelPart("TestPart");
    TetrahedraModelPartForWSSTests(ModelPart);

    wssutilities::CalculateWSS(rModelPart)

    std::vector< array_1d<double,3> > WSSTangentialStress;
    ModelPart.NodesBegin()->GetValueNodes(WSS_TANGENTIAL_STRESS,WSSTangentialStress,ModelPart.GetProcessInfo());

    KRATOS_CHECK_EQUAL(WSSTangentialStress.size(),4);
    for (unsigned int i = 0; i < WSSTangentialStress.size(); i++) {
        KRATOS_CHECK_NEAR(WSSTangentialStress[i][0], 0.0,1e-6); //Change numbers.
        KRATOS_CHECK_NEAR(WSSTangentialStress[i][1], 2.0,1e-6);
        KRATOS_CHECK_NEAR(WSSTangentialStress[i][2], 0.0,1e-6);
    }
}

}
}