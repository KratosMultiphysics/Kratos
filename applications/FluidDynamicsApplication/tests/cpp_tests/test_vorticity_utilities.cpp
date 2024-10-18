//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//
//

// Project includes
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/cfd_variables.h"

// Application includes
#include "custom_elements/fractional_step.h"
#include "custom_utilities/vorticity_utilities.h"
#include "tests/cpp_tests/fluid_dynamics_fast_suite.h"

namespace Kratos {
namespace Testing {

void TriangleModelPartForVorticityTests(ModelPart& rModelPart) {

    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    rModelPart.SetBufferSize(3);
    Properties::Pointer p_properties = rModelPart.CreateNewProperties(0);

    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> element_nodes{1, 2, 3};
    rModelPart.CreateNewElement("FractionalStep2D3N", 1, element_nodes, p_properties);

    // Nodal data
    Element& r_element = *(rModelPart.ElementsBegin());
    Geometry< Node >& r_geometry = r_element.GetGeometry();

    constexpr double omega = 1.0;
    constexpr double center[2] = {1.0,-1.0};

    for (unsigned int i = 0; i < 3; i++) {
        Node& r_node = r_geometry[i];
        double dx = r_node.X() - center[0];
        double dy = r_node.Y() - center[1];
        double r = std::sqrt(dx*dx + dy*dy);
        double cosine = dx/r;
        double sine = dy/r;

        r_node.FastGetSolutionStepValue(VELOCITY_X) = omega*r*sine;
        r_node.FastGetSolutionStepValue(VELOCITY_Y) = -omega*r*cosine;
    }
}

void TetrahedraModelPartForVorticityTests(ModelPart& rModelPart) {

    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    rModelPart.SetBufferSize(3);
    Properties::Pointer p_properties = rModelPart.CreateNewProperties(0);

    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
    rModelPart.CreateNewNode(4, 0.0, 0.0, 1.0);
    std::vector<ModelPart::IndexType> element_nodes{1, 2, 3, 4};
    rModelPart.CreateNewElement("FractionalStep3D4N", 1, element_nodes, p_properties);

    // Nodal data
    Element& r_element = *(rModelPart.ElementsBegin());
    Geometry< Node >& r_geometry = r_element.GetGeometry();

    constexpr double omega = 1.0;
    constexpr double center[3] = {1.0,0.0,-1.0};

    for (unsigned int i = 0; i < 4; i++) {
        Node& r_node = r_geometry[i];
        double dx = r_node.X() - center[0];
        //double dy = r_node.Y() - center[1];
        double dz = r_node.Z() - center[2];
        double r = std::sqrt(dx*dx + dz*dz);
        double cosine = dx/r;
        double sine = dz/r;

        r_node.FastGetSolutionStepValue(VELOCITY_X) = omega*r*sine;
        r_node.FastGetSolutionStepValue(VELOCITY_Z) = -omega*r*cosine;
    }
}

KRATOS_TEST_CASE_IN_SUITE(VorticityUtilities2DQValue, FluidDynamicsApplicationFastSuite) {
    Model model;
    ModelPart& ModelPart = model.CreateModelPart("TestPart");
    TriangleModelPartForVorticityTests(ModelPart);

    std::vector<double> QValues;
    ModelPart.ElementsBegin()->CalculateOnIntegrationPoints(Q_VALUE,QValues,ModelPart.GetProcessInfo());

    KRATOS_EXPECT_EQ(QValues.size(),3);
    for (unsigned int i = 0; i < QValues.size(); i++) {
        KRATOS_EXPECT_NEAR(QValues[i],1.0,1e-6);
    }
}

KRATOS_TEST_CASE_IN_SUITE(VorticityUtilities2DVorticityMagnitude, FluidDynamicsApplicationFastSuite) {
    Model model;
    ModelPart& ModelPart = model.CreateModelPart("TestPart");
    TriangleModelPartForVorticityTests(ModelPart);

    std::vector<double> VorticityMagnitudes;
    ModelPart.ElementsBegin()->CalculateOnIntegrationPoints(VORTICITY_MAGNITUDE,VorticityMagnitudes,ModelPart.GetProcessInfo());

    KRATOS_EXPECT_EQ(VorticityMagnitudes.size(),3);
    for (unsigned int i = 0; i < VorticityMagnitudes.size(); i++) {
        KRATOS_EXPECT_NEAR(VorticityMagnitudes[i],2.0,1e-6);
    }
}

KRATOS_TEST_CASE_IN_SUITE(VorticityUtilities2DVorticity, FluidDynamicsApplicationFastSuite) {
    Model model;
    ModelPart& ModelPart = model.CreateModelPart("TestPart");
    TriangleModelPartForVorticityTests(ModelPart);

    std::vector< array_1d<double,3> > Vorticities;
    ModelPart.ElementsBegin()->CalculateOnIntegrationPoints(VORTICITY,Vorticities,ModelPart.GetProcessInfo());

    KRATOS_EXPECT_EQ(Vorticities.size(),3);
    for (unsigned int i = 0; i < Vorticities.size(); i++) {
        KRATOS_EXPECT_NEAR(Vorticities[i][0], 0.0,1e-6);
        KRATOS_EXPECT_NEAR(Vorticities[i][1], 0.0,1e-6);
        KRATOS_EXPECT_NEAR(Vorticities[i][2],-2.0,1e-6);
    }
}


KRATOS_TEST_CASE_IN_SUITE(VorticityUtilities3DQValue, FluidDynamicsApplicationFastSuite) {
    Model model;
    ModelPart& ModelPart = model.CreateModelPart("TestPart");
    TetrahedraModelPartForVorticityTests(ModelPart);

    std::vector<double> QValues;
    ModelPart.ElementsBegin()->CalculateOnIntegrationPoints(Q_VALUE,QValues,ModelPart.GetProcessInfo());

    KRATOS_EXPECT_EQ(QValues.size(),4);
    for (unsigned int i = 0; i < QValues.size(); i++) {
        KRATOS_EXPECT_NEAR(QValues[i],1.0,1e-6);
    }
}

KRATOS_TEST_CASE_IN_SUITE(VorticityUtilities3DVorticityMagnitude, FluidDynamicsApplicationFastSuite) {
    Model model;
    ModelPart& ModelPart = model.CreateModelPart("TestPart");
    TetrahedraModelPartForVorticityTests(ModelPart);

    std::vector<double> VorticityMagnitudes;
    ModelPart.ElementsBegin()->CalculateOnIntegrationPoints(VORTICITY_MAGNITUDE,VorticityMagnitudes,ModelPart.GetProcessInfo());

    KRATOS_EXPECT_EQ(VorticityMagnitudes.size(),4);
    for (unsigned int i = 0; i < VorticityMagnitudes.size(); i++) {
        KRATOS_EXPECT_NEAR(VorticityMagnitudes[i],2.0,1e-6);
    }
}

KRATOS_TEST_CASE_IN_SUITE(VorticityUtilities3DVorticity, FluidDynamicsApplicationFastSuite) {
    Model model;
    ModelPart& ModelPart = model.CreateModelPart("TestPart");
    TetrahedraModelPartForVorticityTests(ModelPart);

    std::vector< array_1d<double,3> > Vorticities;
    ModelPart.ElementsBegin()->CalculateOnIntegrationPoints(VORTICITY,Vorticities,ModelPart.GetProcessInfo());

    KRATOS_EXPECT_EQ(Vorticities.size(),4);
    for (unsigned int i = 0; i < Vorticities.size(); i++) {
        KRATOS_EXPECT_NEAR(Vorticities[i][0], 0.0,1e-6);
        KRATOS_EXPECT_NEAR(Vorticities[i][1], 2.0,1e-6);
        KRATOS_EXPECT_NEAR(Vorticities[i][2], 0.0,1e-6);
    }
}

}
}
