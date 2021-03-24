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
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/cfd_variables.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "custom_constitutive/newtonian_2d_law.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(EmbeddedElementDiscontinuous2D3N, FluidDynamicsApplicationFastSuite)
{
    Model model;
    ModelPart& model_part = model.CreateModelPart("Main",3);

    // Variables addition
    model_part.AddNodalSolutionStepVariable(DENSITY);
    model_part.AddNodalSolutionStepVariable(REACTION);
    model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    model_part.AddNodalSolutionStepVariable(DISTANCE);
    model_part.AddNodalSolutionStepVariable(ACCELERATION);
    model_part.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);
    model_part.AddNodalSolutionStepVariable(REACTION_WATER_PRESSURE);

    // For VMS comparison
    model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    model_part.AddNodalSolutionStepVariable(ADVPROJ);
    model_part.AddNodalSolutionStepVariable(DIVPROJ);

    // Process info creation
    const double delta_time = 0.1;
    const double sound_velocity = 1.0e+3;
    model_part.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
    model_part.GetProcessInfo().SetValue(SOUND_VELOCITY, sound_velocity);
    model_part.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
    Vector bdf_coefs(3);
    bdf_coefs[0] = 3.0/(2.0*delta_time);
    bdf_coefs[1] = -2.0/delta_time;
    bdf_coefs[2] = 0.5*delta_time;
    model_part.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);

    // Set the element properties
    const double rho = 1000.0;
    Properties::Pointer p_properties = model_part.CreateNewProperties(0);
    p_properties->SetValue(DENSITY, rho);
    p_properties->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
    ConstitutiveLaw::Pointer pConsLaw(new Newtonian2DLaw());
    p_properties->SetValue(CONSTITUTIVE_LAW, pConsLaw);

    // Geometry creation
    model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    model_part.CreateNewNode(3, 0.0, 1.0, 0.0);

    for (ModelPart::NodeIterator it_node=model_part.NodesBegin(); it_node<model_part.NodesEnd(); ++it_node){
        it_node->AddDof(VELOCITY_X,REACTION_X);
        it_node->AddDof(VELOCITY_Y,REACTION_Y);
        it_node->AddDof(VELOCITY_Z,REACTION_Z);
        it_node->AddDof(PRESSURE,REACTION_WATER_PRESSURE);
    }

    std::vector<ModelPart::IndexType> element_nodes {1, 2, 3};
    model_part.CreateNewElement("EmbeddedWeaklyCompressibleNavierStokesDiscontinuous2D3N", 1, element_nodes, p_properties);
    model_part.CreateNewElement("EmbeddedQSVMSDiscontinuous2D3N", 2, element_nodes, p_properties);

    const auto& r_process_info = model_part.GetProcessInfo();
    for (auto &r_elem : model_part.Elements()) {
        r_elem.Initialize(r_process_info); // Initialize the element to initialize the constitutive law
        const auto& r_const_elem = r_elem;
        r_const_elem.Check(r_process_info); // Otherwise the constitutive law is not seen here
    }

    // Define the nodal values
    Matrix reference_velocity(3,2);
    reference_velocity(0,0) = 0.0; reference_velocity(0,1) = 0.1;
    reference_velocity(1,0) = 0.1; reference_velocity(1,1) = 0.2;
    reference_velocity(2,0) = 0.2; reference_velocity(2,1) = 0.3;

    Element::Pointer p_element = model_part.pGetElement(1);

    for(unsigned int i=0; i<3; i++){
        p_element->GetGeometry()[i].SetValue(SOUND_VELOCITY, sound_velocity); // Required for the weakly compressible instance
        p_element->GetGeometry()[i].FastGetSolutionStepValue(DENSITY) = rho; // Required for the nodal-based density instances
        p_element->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE)    = 0.0;
        p_element->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE, 1) = 0.0;
        p_element->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE, 2) = 0.0;
        for(unsigned int k=0; k<2; k++){
            p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k]    = reference_velocity(i,k);
            p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9*reference_velocity(i,k);
            p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.75*reference_velocity(i,k);
            p_element->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k]    = 0.0;
            p_element->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
            p_element->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;
        }
    }

    // RHS and LHS
    Vector RHS = ZeroVector(9);
    Matrix LHS = ZeroMatrix(9,9);

    std::vector< std::vector<double> > output_uncut(6);
    output_uncut[0] = {-0.6361617846,8.819948812,-0.6557582459,67.57989341,174.5435981,0.1308775154,110.444523,215.3506723,0.3748807306}; // EmbeddedWeaklyCompressibleNavierStokesDiscontinuous
    output_uncut[1] = {-21.81650306,-40.75920676,-0.6557581669,54.90454836,132.1891487,0.1308774929,90.0369547,179.8200581,0.374880674}; // EmbeddedQSVMSDiscontinuous
    int counter = 0;

    // Test Uncut element
    array_1d<double,3> elem_dist;
    elem_dist[0] = 1.0;
    elem_dist[1] = 1.0;
    elem_dist[2] = 1.0;
    for (auto it_elem = model_part.ElementsBegin(); it_elem != model_part.ElementsEnd(); ++it_elem) {
        it_elem->SetValue(ELEMENTAL_DISTANCES, elem_dist);
    }

    for (ModelPart::ElementIterator i = model_part.ElementsBegin(); i != model_part.ElementsEnd(); i++) {
        i->CalculateLocalSystem(LHS, RHS, r_process_info);

        // std::cout << i->Info() << std::setprecision(10) << std::endl;
        // KRATOS_WATCH(RHS);

        for (unsigned int j = 0; j < RHS.size(); j++) {
            KRATOS_CHECK_NEAR(RHS[j], output_uncut[counter][j], 1e-6);
        }

        counter++;
    }

    std::vector< std::vector<double> > output_cut(6);
    output_cut[0] = {18.84223125,59.39337635,-0.4453265312,49.82664899,169.5161544,0.3953265312,32.91666657,-24.94122968,-0.1}; // EmbeddedWeaklyCompressibleNavierStokesDiscontinuous
    output_cut[1] = {3.777844188,12.25294162,-0.4453264623,42.19438001,129.199604,0.3953264623,32.91666657,-24.94122968,-0.1}; // EmbeddedQSVMSDiscontinuous
    counter = 0;

    // Test cut element
    elem_dist[0] = -1.0;
    elem_dist[1] = -1.0;
    elem_dist[2] =  0.5;
    for (auto it_elem = model_part.ElementsBegin(); it_elem != model_part.ElementsEnd(); ++it_elem) {
        it_elem->SetValue(ELEMENTAL_DISTANCES, elem_dist);
    }

    model_part.GetProcessInfo().SetValue(SLIP_LENGTH, 0.0);
    model_part.GetProcessInfo().SetValue(PENALTY_COEFFICIENT, 0.1);
    for (ModelPart::ElementIterator i = model_part.ElementsBegin(); i != model_part.ElementsEnd(); i++) {
        i->Set(SLIP, false);
        i->CalculateLocalSystem(LHS, RHS, r_process_info);

        // std::cout << i->Info() << std::setprecision(10) << std::endl;
        // KRATOS_WATCH(RHS);

        for (unsigned int j = 0; j < RHS.size(); j++) {
            KRATOS_CHECK_NEAR(RHS[j], output_cut[counter][j], 1e-6);
        }

        counter++;
    }

    std::vector< std::vector<double> > output_slip_cut(6);
    output_slip_cut[0] = {18.84227218,59.39341828,-0.4453265312,49.82660608,169.5161125,0.3953265312,32.91666667,-24.94122968,-0.1}; // EmbeddedWeaklyCompressibleNavierStokesDiscontinuous
    output_slip_cut[1] = {3.777885122,12.25298355,-0.4453264623,42.1943371,129.199562,0.3953264623,32.91666667,-24.94122968,-0.1}; // EmbeddedQSVMSDiscontinuous
    counter = 0;

    // Test slip cut element
    model_part.GetProcessInfo().SetValue(SLIP_LENGTH, 1.0e+08);
    model_part.GetProcessInfo().SetValue(PENALTY_COEFFICIENT, 0.1);
    for (ModelPart::ElementIterator i = model_part.ElementsBegin(); i != model_part.ElementsEnd(); ++i) {
        i->Set(SLIP, true);
        i->CalculateLocalSystem(LHS, RHS, r_process_info);

        // std::cout << i->Info() << std::setprecision(10) << std::endl;
        // KRATOS_WATCH(RHS);

        for (unsigned int j = 0; j < RHS.size(); j++) {
            KRATOS_CHECK_NEAR(RHS[j], output_slip_cut[counter][j], 1e-6);
        }

        counter++;
    }
}

}  // namespace Testing
}  // namespace Kratos
