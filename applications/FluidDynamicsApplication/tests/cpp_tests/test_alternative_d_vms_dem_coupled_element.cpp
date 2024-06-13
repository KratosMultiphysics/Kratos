//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Joaquin Gonzalez-Usua
//

#include <iomanip> // for std::setprecision

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/cfd_variables.h"
#include "fluid_dynamics_application_variables.h"

// Application includes
#include "custom_constitutive/newtonian_2d_law.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(AlternativeDVMSDEMCoupled2D4N, FluidDynamicsApplicationFastSuite)
{
    Model model;
    unsigned int buffer_size = 2;
    unsigned int Dim = 2;
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
    model_part.AddNodalSolutionStepVariable(REACTION);
    model_part.AddNodalSolutionStepVariable(REACTION_WATER_PRESSURE);
    model_part.AddNodalSolutionStepVariable(MASS_SOURCE);
    model_part.AddNodalSolutionStepVariable(FLUID_FRACTION);
    model_part.AddNodalSolutionStepVariable(FLUID_FRACTION_RATE);
    model_part.AddNodalSolutionStepVariable(FLUID_FRACTION_GRADIENT);
    model_part.AddNodalSolutionStepVariable(PERMEABILITY);

    // Process info creation
    double delta_time = 0.1;
    model_part.GetProcessInfo().SetValue(DELTA_TIME, delta_time);

    // Set the element properties
    Properties::Pointer p_properties = model_part.CreateNewProperties(0);
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

        double& r_fluid_fraction = it_node->FastGetSolutionStepValue(FLUID_FRACTION);
        r_fluid_fraction = 1.0;
        Matrix& r_permeability = it_node->FastGetSolutionStepValue(PERMEABILITY);
        r_permeability = ZeroMatrix(Dim, Dim);
        for (unsigned int d = 0; d < Dim; ++d){
            r_permeability(d,d) = 0.0;
        }
    }

    std::vector<ModelPart::IndexType> element_nodes {1, 2, 4, 3};
    model_part.CreateNewElement("AlternativeDVMSDEMCoupled2D4N", 1, element_nodes, p_properties);

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


    Geometry<Node>& r_geometry = model_part.ElementsBegin()->GetGeometry();


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

    std::vector<double> output = {2.885342223,3.817976395,-0.05954391544,-3.764678983,2.452222894,-0.04746809897,-5.675304236,-5.32981538,-0.04009093622,1.144992949,-6.350031955,-0.05289704936}; // AlternativeDVMSDEMCoupled2D4N

    for (ModelPart::ElementIterator i = model_part.ElementsBegin(); i != model_part.ElementsEnd(); i++) {
        const auto& r_process_info = model_part.GetProcessInfo();
        i->Initialize(r_process_info); // Initialize constitutive law
        const auto& rElem = *i;
        rElem.Check(r_process_info);
        i->InitializeNonLinearIteration(r_process_info);
        i->CalculateLocalVelocityContribution(LHS, RHS, r_process_info);

        // std::cout << i->Info() << std::setprecision(10) << std::endl;
        // KRATOS_WATCH(RHS);

        for (unsigned int j = 0; j < output.size(); j++) {
            KRATOS_EXPECT_NEAR(RHS[j], output[j], 1e-5);
        }
    }
    double porosity = 0.5;
    for (ModelPart::NodeIterator it_node=model_part.NodesBegin(); it_node<model_part.NodesEnd(); ++it_node){
        double& r_fluid_fraction = it_node->FastGetSolutionStepValue(FLUID_FRACTION);
        r_fluid_fraction = porosity;
        Matrix& r_permeability = it_node->FastGetSolutionStepValue(PERMEABILITY);
        r_permeability = ZeroMatrix(Dim, Dim);
    }

    for (ModelPart::ElementIterator i = model_part.ElementsBegin(); i != model_part.ElementsEnd(); i++) {
        const auto& r_process_info = model_part.GetProcessInfo();
        i->Initialize(r_process_info); // Initialize constitutive law
        const auto& rElem = *i;
        rElem.Check(r_process_info);
        i->InitializeNonLinearIteration(r_process_info);
        i->CalculateLocalVelocityContribution(LHS, RHS, r_process_info);

        for (unsigned int j = 0; j < output.size(); j++) {
            KRATOS_CHECK_NEAR(RHS[j], porosity*output[j], 1e-5);
        }
    }
}

}  // namespace Testing
}  // namespace Kratos