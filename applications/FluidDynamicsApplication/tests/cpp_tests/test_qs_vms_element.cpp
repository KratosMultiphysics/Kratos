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
//                   Ruben Zorrilla
//

#include <iomanip> // for std::setprecision

// Project includes
#include "testing/testing.h"
#include "includes/model_part.h"
#include "includes/cfd_variables.h"

// Application includes
#include "custom_constitutive/newtonian_2d_law.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(QSVMS2D4N, FluidDynamicsApplicationFastSuite)
{
    Model model;
    unsigned int buffer_size = 2;
    ModelPart& model_part = model.CreateModelPart("Main",buffer_size);

    // Variables addition
    model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    model_part.AddNodalSolutionStepVariable(ACCELERATION);
    model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    model_part.AddNodalSolutionStepVariable(ADVPROJ);
    model_part.AddNodalSolutionStepVariable(DIVPROJ);

    // Process info creation
    double delta_time = 0.1;
    model_part.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
    model_part.GetProcessInfo().SetValue(DELTA_TIME, delta_time);

    // Set the element properties
    Properties::Pointer p_properties = model_part.pGetProperties(0);
    p_properties->SetValue(DENSITY, 1000.0);
    p_properties->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
    ConstitutiveLaw::Pointer pConsLaw(new Newtonian2DLaw());
    p_properties->SetValue(CONSTITUTIVE_LAW, pConsLaw);

    // Geometry creation
    model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    model_part.CreateNewNode(2, 1.0, 0.1, 0.0);
    model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    model_part.CreateNewNode(4, 1.0, 0.9, 0.0);

    for (ModelPart::NodeIterator it_node=model_part.NodesBegin(); it_node<model_part.NodesEnd(); ++it_node){
        it_node->AddDof(VELOCITY_X,REACTION_X);
        it_node->AddDof(VELOCITY_Y,REACTION_Y);
        it_node->AddDof(VELOCITY_Z,REACTION_Z);
        it_node->AddDof(PRESSURE,REACTION_WATER_PRESSURE);
    }

    std::vector<ModelPart::IndexType> element_nodes {1, 2, 4, 3};
    model_part.CreateNewElement("QSVMS2D4N", 1, element_nodes, p_properties);

    // Loop starts at 1 because you need one less clone than time steps (JC)
    for (unsigned int i = 1; i < buffer_size; i++) {
        model_part.CloneTimeStep(i * delta_time);
    }

    // Define the nodal values
    Matrix reference_velocity(4,2);
    reference_velocity(0,0) = 0.0; reference_velocity(0,1) = 0.1;
    reference_velocity(1,0) = 0.1; reference_velocity(1,1) = 0.2;
    reference_velocity(2,0) = 0.2; reference_velocity(2,1) = 0.3;
    reference_velocity(3,0) = 0.3; reference_velocity(3,1) = 0.4;


    Geometry<Node<3>>& r_geometry = model_part.ElementsBegin()->GetGeometry();


    for(unsigned int i=0; i<4; i++){
        r_geometry[i].FastGetSolutionStepValue(PRESSURE)    = 0.0;
        r_geometry[i].FastGetSolutionStepValue(PRESSURE, 1) = 0.0;
        for(unsigned int k=0; k<2; k++){
            r_geometry[i].FastGetSolutionStepValue(VELOCITY)[k]    = reference_velocity(i,k);
            r_geometry[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9*reference_velocity(i,k);
            r_geometry[i].FastGetSolutionStepValue(MESH_VELOCITY)[k]    = 0.0;
            r_geometry[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
        }
    }

    // RHS and LHS
    Vector RHS = ZeroVector(12);
    Matrix LHS = ZeroMatrix(12,12);

    std::vector<double> output = {-2.665425819,-1.87894198,-0.02477280423,-10.27651236,-5.037560437,-0.05013494554,-19.87147169,-19.57971097,-0.06709466598,-14.8532568,-21.17045328,-0.05799758425}; // QSVMS2D4N

    for (ModelPart::ElementIterator i = model_part.ElementsBegin(); i != model_part.ElementsEnd(); i++) {
        i->Initialize(); // Initialize constitutive law
        i->Check(model_part.GetProcessInfo());
        i->CalculateLocalVelocityContribution(LHS, RHS, model_part.GetProcessInfo());

        //std::cout << i->Info() << std::setprecision(10) << std::endl;
        //KRATOS_WATCH(RHS);

        for (unsigned int j = 0; j < output.size(); j++) {
            KRATOS_CHECK_NEAR(RHS[j], output[j], 1e-6);
        }
    }
}

}  // namespace Testing
}  // namespace Kratos