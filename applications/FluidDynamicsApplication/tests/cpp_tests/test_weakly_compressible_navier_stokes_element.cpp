//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes
#include <iomanip> // for std::setprecision

// External includes


// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/cfd_variables.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "custom_constitutive/newtonian_2d_law.h"

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(WeaklyCompressibleNavierStokes2D3N, FluidDynamicsApplicationFastSuite)
{
    Model model;
    std::size_t buffer_size = 3;
    auto& r_model_part = model.CreateModelPart("TestModelPart", buffer_size);

    // Variables addition
    r_model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    r_model_part.AddNodalSolutionStepVariable(DENSITY);
    r_model_part.AddNodalSolutionStepVariable(PRESSURE);
    r_model_part.AddNodalSolutionStepVariable(VELOCITY);
    r_model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    r_model_part.AddNodalSolutionStepVariable(ACCELERATION);
    r_model_part.AddNodalSolutionStepVariable(REACTION);
    r_model_part.AddNodalSolutionStepVariable(REACTION_WATER_PRESSURE);

    // ProcessInfo container fill
    double delta_time = 0.1;
    r_model_part.GetProcessInfo().SetValue(DYNAMIC_TAU, 1.0);
    r_model_part.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
    Vector bdf_coefs(3);
    bdf_coefs[0] = 3.0 / (2.0 * delta_time);
    bdf_coefs[1] = -2.0 / delta_time;
    bdf_coefs[2] = 0.5 * delta_time;
    r_model_part.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);

    // Set the element properties
    auto p_properties = r_model_part.CreateNewProperties(0);
    p_properties->SetValue(DENSITY, 1.0);
    p_properties->SetValue(DYNAMIC_VISCOSITY, 1.0e-3);
    auto p_cons_law = Kratos::make_shared<Newtonian2DLaw>();
    p_properties->SetValue(CONSTITUTIVE_LAW, p_cons_law);

    // Element creation
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 0.1, 0.0);
    r_model_part.CreateNewNode(3, 1.0, 0.9, 0.0);

    for (auto it_node = r_model_part.NodesBegin(); it_node < r_model_part.NodesEnd(); ++it_node){
        it_node->AddDof(VELOCITY_X,REACTION_X);
        it_node->AddDof(VELOCITY_Y,REACTION_Y);
        it_node->AddDof(VELOCITY_Z,REACTION_Z);
        it_node->AddDof(PRESSURE,REACTION_WATER_PRESSURE);
    }

    std::vector<ModelPart::IndexType> element_nodes {1, 2, 3};
    auto p_elem = r_model_part.CreateNewElement("WeaklyCompressibleNavierStokes2D3N", 1, element_nodes, p_properties);

    // Define and set the nodal values
    Matrix reference_velocity(4,2);
    reference_velocity(0,0) = 0.0; reference_velocity(0,1) = 0.1;
    reference_velocity(1,0) = 0.1; reference_velocity(1,1) = 0.2;
    reference_velocity(2,0) = 0.2; reference_velocity(2,1) = 0.3;

    auto& r_geometry = r_model_part.ElementsBegin()->GetGeometry();
    for (std::size_t i=0; i<3; i++){
        r_geometry[i].SetValue(SOUND_VELOCITY, 1.0e3);
        r_geometry[i].FastGetSolutionStepValue(DENSITY) = 1.0;
        r_geometry[i].FastGetSolutionStepValue(PRESSURE) = 0.0;
        r_geometry[i].FastGetSolutionStepValue(PRESSURE, 1) = 0.0;
        r_geometry[i].FastGetSolutionStepValue(PRESSURE, 2) = 0.0;
        for (std::size_t k=0; k<2; k++){
            r_geometry[i].FastGetSolutionStepValue(VELOCITY)[k] = reference_velocity(i,k);
            r_geometry[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9*reference_velocity(i,k);
            r_geometry[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.7*reference_velocity(i,k);
        }
    }

    // Set the elemental values
    p_elem->SetValue(RESISTANCE, 1.0e3);

    // Calculate RHS and LHS
    Vector RHS = ZeroVector(9);
    Matrix LHS = ZeroMatrix(9,9);
    const auto& r_process_info = r_model_part.GetProcessInfo();
    p_elem->Initialize(r_process_info);
    p_elem->CalculateLocalSystem(LHS, RHS, r_process_info);

    // std::cout << p_elem->Info() << std::setprecision(12) << std::endl;
    // KRATOS_WATCH(RHS)
    // KRATOS_WATCH(row(LHS,0))

    // Check values
    const std::vector<double> rhs_ref = {34.2921172632, -19.1612962435, -0.0212935839184, -58.7701017724, 31.2334272413, -0.0186565395118, -8.37866202262, -77.7742128742, -0.0450498765698}; // WeaklyCompressibleNavierStokes2D3N
    const std::vector<double> lhs_0_ref = {255.745969647, 0, 0.155863699867, -197.183916289, 250.033529091, 0.105861405575, 52.887840474, -250.033529091, 0.13527448821};                     // WeaklyCompressibleNavierStokes2D3N
    KRATOS_EXPECT_VECTOR_NEAR(RHS, rhs_ref, 1.0e-10)
    KRATOS_EXPECT_VECTOR_NEAR(row(LHS,0), lhs_0_ref, 1.0e-8)
}

}  // namespace Kratos::Testing