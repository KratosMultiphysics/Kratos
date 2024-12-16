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
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/cfd_variables.h"

// Application includes
#include "tests/cpp_tests/fluid_dynamics_fast_suite.h"


namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(LowMachNavierStokes2D3N, FluidDynamicsApplicationFastSuite)
{
    Model model;
    unsigned int buffer_size = 3;
    auto& r_model_part = model.CreateModelPart("TestModelPart", buffer_size);

    // Variables addition
    r_model_part.AddNodalSolutionStepVariable(PRESSURE);
    r_model_part.AddNodalSolutionStepVariable(VELOCITY);
    r_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    r_model_part.AddNodalSolutionStepVariable(REACTION);
    r_model_part.AddNodalSolutionStepVariable(REACTION_WATER_PRESSURE);
    r_model_part.AddNodalSolutionStepVariable(REACTION_FLUX);
    r_model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    r_model_part.AddNodalSolutionStepVariable(DENSITY);
    r_model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    r_model_part.AddNodalSolutionStepVariable(HEAT_FLUX);

    // ProcessInfo container fill
    double delta_time = 0.1;
    r_model_part.GetProcessInfo().SetValue(DYNAMIC_TAU, 1.0);
    r_model_part.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
    Vector bdf_coefs(3);
    bdf_coefs[0] = 3.0 / (2.0 * delta_time);
    bdf_coefs[1] = -2.0 / delta_time;
    bdf_coefs[2] = 0.5 * delta_time;
    r_model_part.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);

    // Set thermodynamic pressure to the atmosferic one (open flows case)
    const double p_th = 101325;
    r_model_part.GetProcessInfo().SetValue(PRESSURE, p_th);

    // Set the element properties
    const double c_p = 1.0e-3;
    const double gamma = 1.4e0;
    auto p_properties = r_model_part.CreateNewProperties(0);
    p_properties->SetValue(SPECIFIC_HEAT, c_p);
    p_properties->SetValue(CONDUCTIVITY, 1.0e-3);
    p_properties->SetValue(DYNAMIC_VISCOSITY, 1.0e-3);
    p_properties->SetValue(HEAT_CAPACITY_RATIO, gamma);
    auto p_cons_law = Kratos::make_shared<Newtonian2DLaw>();
    p_properties->SetValue(CONSTITUTIVE_LAW, p_cons_law);

    // Element creation
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 0.1, 0.0);
    r_model_part.CreateNewNode(3, 1.0, 0.9, 0.0);

    for (auto it_node = r_model_part.NodesBegin(); it_node < r_model_part.NodesEnd(); ++it_node){
        it_node->AddDof(PRESSURE, REACTION_WATER_PRESSURE);
        it_node->AddDof(VELOCITY_X, REACTION_X);
        it_node->AddDof(VELOCITY_Y, REACTION_Y);
        it_node->AddDof(VELOCITY_Z, REACTION_Z);
        it_node->AddDof(TEMPERATURE,REACTION_FLUX);
    }

    std::vector<ModelPart::IndexType> element_nodes {1, 2, 3};
    auto p_elem = r_model_part.CreateNewElement("LowMachNavierStokes2D3N", 1, element_nodes, p_properties);

    // Define and set the nodal values
    Matrix reference_velocity(3,2);
    reference_velocity(0,0) = 0.0; reference_velocity(0,1) = 0.1;
    reference_velocity(1,0) = 0.1; reference_velocity(1,1) = 0.2;
    reference_velocity(2,0) = 0.2; reference_velocity(2,1) = 0.3;

    auto& r_geometry = r_model_part.ElementsBegin()->GetGeometry();
    for(unsigned int i = 0; i < 3; ++i){
        const unsigned int id = r_geometry[i].Id();
        const double temp = 273.0 * (id / 3.0);
        r_geometry[i].FastGetSolutionStepValue(PRESSURE) = 0.0;
        r_geometry[i].FastGetSolutionStepValue(PRESSURE, 1) = 0.0;
        r_geometry[i].FastGetSolutionStepValue(PRESSURE, 2) = 0.0;
        r_geometry[i].FastGetSolutionStepValue(TEMPERATURE) = temp;
        r_geometry[i].FastGetSolutionStepValue(TEMPERATURE, 1) = 283.0 * (id / 3.0);
        r_geometry[i].FastGetSolutionStepValue(TEMPERATURE, 2) = 293.0 * (id / 3.0);
        for(unsigned int k = 0; k < 2; ++k){
            r_geometry[i].FastGetSolutionStepValue(VELOCITY)[k] = reference_velocity(i,k);
            r_geometry[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9 * reference_velocity(i,k);
            r_geometry[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.7 * reference_velocity(i,k);
        }
        r_geometry[i].FastGetSolutionStepValue(DENSITY) = gamma * p_th / c_p / (gamma - 1.0) / temp;
    }

    // Calculate RHS and LHS
    Vector RHS = ZeroVector(12);
    Matrix LHS = ZeroMatrix(12,12);
    const auto& r_process_info = r_model_part.GetProcessInfo();
    p_elem->Initialize(r_process_info);

    KRATOS_WATCH("Before LocalSystem")
    p_elem->CalculateLocalSystem(LHS, RHS, r_process_info);
    KRATOS_WATCH("After LocalSystem")

    std::cout << p_elem->Info() << std::setprecision(12) << std::endl;
    KRATOS_WATCH(RHS)
    KRATOS_WATCH(row(LHS,0))

    // Check values
    // const std::vector<double> rhs_ref = {0.101208740426, 0.064297119485, -0.254190386408, 0.0115211432149, 0.137376654428, -0.291181505801, 0.164573097634, 0.174512378122, -0.264759715687}; // AxisymmetricNavierStokes2D3N
    // const std::vector<double> lhs_0_ref = {1.26303343879, -0.048478920937, 0.214700065233, 0.627949407008, 0.0236778793736, 0.224465336895, 1.17170148289, -0.145410069083, 0.398592638829};  // AxisymmetricNavierStokes2D3N
    // KRATOS_EXPECT_VECTOR_NEAR(RHS, rhs_ref, 1.0e-10)
    // KRATOS_EXPECT_VECTOR_NEAR(row(LHS,0), lhs_0_ref, 1.0e-10)
}

}  // namespace Kratos::Testing