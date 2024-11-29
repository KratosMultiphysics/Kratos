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
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/cfd_variables.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "custom_constitutive/newtonian_2d_law.h"
#include "tests/cpp_tests/fluid_dynamics_fast_suite.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(EmbeddedElement2D3N, FluidDynamicsApplicationFastSuite)
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
    model_part.CreateNewElement("EmbeddedNavierStokes2D3N", 1, element_nodes, p_properties);
    model_part.CreateNewElement("EmbeddedWeaklyCompressibleNavierStokes2D3N", 2, element_nodes, p_properties);
    model_part.CreateNewElement("NavierStokes2D3N", 3, element_nodes, p_properties);
    model_part.CreateNewElement("WeaklyCompressibleNavierStokes2D3N", 4, element_nodes, p_properties);
    model_part.CreateNewElement("TimeIntegratedQSVMS2D3N", 5, element_nodes, p_properties);
    model_part.CreateNewElement("EmbeddedQSVMS2D3N", 6, element_nodes, p_properties);

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
    output_uncut[0] = {7.044040063, 26.56103979, -0.5033880823, 60.93071102, 153.3702507, 0.08535256462, 95.7560183, 186.7293649, 0.2680355177}; // EmbeddedNavierStokes
    output_uncut[1] = {-0.6361617846, 8.819948812, -0.6557582459, 67.57989341, 174.5435981, 0.1308775154, 110.444523, 215.3506723, 0.3748807306}; // EmbeddedFluidElement
    output_uncut[2] = {7.044040063, 26.56103979, -0.5033880823, 60.93071102, 153.3702507, 0.08535256462, 95.7560183, 186.7293649, 0.2680355177};  // NavierStokes
    output_uncut[3] = {-0.6361617846, 8.819948812, -0.6557582459, 67.57989341, 174.5435981, 0.1308775154, 110.444523, 215.3506723, 0.3748807306}; // WeaklyCompressibleNavierStokes2D3N
    output_uncut[4] = {-21.81650306, -40.75920676, -0.6557581669, 54.90454836, 132.1891487, 0.1308774929, 90.0369547, 179.8200581, 0.374880674};  // TimeIntegratedQSVMS
    output_uncut[5] = {-21.81650306, -40.75920676, -0.6557581669, 54.90454836, 132.1891487, 0.1308774929, 90.0369547, 179.8200581, 0.374880674};  // EmbeddedQSVMS
    int counter = 0;

    // Test Uncut element
    p_element->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = 1.0;
    p_element->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = 1.0;
    p_element->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = 1.0;

    for (ModelPart::ElementIterator i = model_part.ElementsBegin(); i != model_part.ElementsEnd(); i++) {
        i->CalculateLocalSystem(LHS, RHS, r_process_info);

        // std::cout << i->Info() << std::setprecision(10) << std::endl;
        // KRATOS_WATCH(RHS);

        for (unsigned int j = 0; j < RHS.size(); j++) {
            KRATOS_EXPECT_NEAR(RHS[j], output_uncut[counter][j], 1e-6);
        }

        counter++;
    }

    std::vector< std::vector<double> > output_cut(6);
    output_cut[0] = {-0.008024691358, -0.01358024691, -0.05463909243, -0.008641975309, -0.01419753086, 0.01767964152, 28.23027063, 46.48478011, 0.02029278424}; // EmbeddedNavierStokes
    output_cut[1] = {-0.008024691358, -0.01358024691, -0.07247223359, -0.008641975309, -0.01419753086, 0.02427807437, 31.53550661, 51.62714426, 0.03152749255}; // EmbeddedFluidElement
    output_cut[2] = {7.044040063,26.56103979,-0.5033880823,60.93071102,153.3702507,0.08535256462,95.7560183,186.7293649,0.2680355177}; // NavierStokes
    output_cut[3] = {-0.6361617846,8.819948812,-0.6557582459,67.57989341,174.5435981,0.1308775154,110.444523,215.3506723,0.3748807306}; // WeaklyCompressibleNavierStokes2D3N
    output_cut[4] = {-21.81650306,-40.75920676,-0.6557581669,54.90454836,132.1891487,0.1308774929,90.0369547,179.8200581,0.374880674}; // TimeIntegratedQSVMS
    output_cut[5] = {-0.008024691358, -0.01358024691, -0.07247222729, -0.008641975309, -0.01419753086, 0.02427807205, 25.42462911, 42.1803488, 0.03152748858}; // EmbeddedQSVMS
    counter = 0;

    // Test cut element
    p_element->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
    p_element->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = -1.0;
    p_element->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) =  0.5;

    for (ModelPart::ElementIterator i = model_part.ElementsBegin(); i != model_part.ElementsEnd(); i++) {
        i->Set(SLIP, false);
        i->CalculateLocalSystem(LHS, RHS, r_process_info);

        //std::cout << i->Info() << std::setprecision(10) << std::endl;
        //KRATOS_WATCH(RHS);

        for (unsigned int j = 0; j < RHS.size(); j++) {
            KRATOS_EXPECT_NEAR(RHS[j], output_cut[counter][j], 1e-6);
        }

        counter++;
    }

    std::vector< std::vector<double> > output_slip_cut(6);
    output_slip_cut[0] = {-3.97751477, -753.7076317, -0.06821933934, 5.148358206, -768.9308712, 0.003482110658, 28.23027124, -3005.776109, -0.03526277131}; // EmbeddedNavierStokes
    output_slip_cut[1] = {-6.153635179,-1010.492889,-0.0860524805,5.998772503,-1033.104416,0.01008054351,31.53550721,-4039.007447,-0.02402806301};  // EmbeddedFluidElement
    output_slip_cut[2] = {7.044040063,26.56103979,-0.5033880823,60.93071102,153.3702507,0.08535256462,95.7560183,186.7293649,0.2680355177}; // NavierStokes
    output_slip_cut[3] = {-0.6361617846,8.819948812,-0.6557582459,67.57989341,174.5435981,0.1308775154,110.444523,215.3506723,0.3748807306}; // WeaklyCompressibleNavierStokes2D3N
    output_slip_cut[4] = {-21.81650306,-40.75920676,-0.6557581669,54.90454836,132.1891487,0.1308774929,90.0369547,179.8200581,0.374880674}; // TimeIntegratedQSVMS
    output_slip_cut[5] = {-7.936602916,-1012.910524,-0.0860524742,6.053639874,-1034.587123,0.01008054118,25.42462971,-4048.454242,-0.02402806698}; // EmbeddedQSVMS
    counter = 0;

    // Test slip cut element
    model_part.GetProcessInfo().SetValue(SLIP_LENGTH, 1.0e+08);
    model_part.GetProcessInfo().SetValue(PENALTY_COEFFICIENT, 10.0);
    for (ModelPart::ElementIterator i = model_part.ElementsBegin(); i != model_part.ElementsEnd(); ++i) {
        i->Set(SLIP, true);
        i->CalculateLocalSystem(LHS, RHS, r_process_info);

        //std::cout << i->Info() << std::setprecision(10) << std::endl;
        //KRATOS_WATCH(RHS);

        for (unsigned int j = 0; j < RHS.size(); j++) {
            KRATOS_EXPECT_NEAR(RHS[j], output_slip_cut[counter][j], 1e-6);
        }

        counter++;
    }

    std::vector< std::vector<double> > output_embedded_velocity;
    output_embedded_velocity.resize(6);
    output_embedded_velocity[0] = {0.0475308641975, 0.0975308641975, -0.0546390924304, 0.0469135802469, 0.0969135802469, 0.0176796415227, 16436.8507492, 33828.9387065, 0.020292784241}; // EmbeddedNavierStokes
    output_embedded_velocity[1] = {0.0475308641975,0.0975308641975,-0.0724722335903,0.0469135802469,0.0969135802469,0.0242780743742,29532.6149603,60789.7777059,0.0315274925494}; // EmbeddedFluidElement
    output_embedded_velocity[2] = {7.044040063,26.56103979,-0.5033880823,60.93071102,153.3702507,0.08535256462,95.7560183,186.7293649,0.2680355177}; // NavierStokes
    output_embedded_velocity[3] = {-0.636161784635,8.8199488124,-0.655758245947,67.5798934101,174.543598083,0.130877515367,110.444522985,215.350672279,0.37488073058}; // WeaklyCompressibleNavierStokes2D3N
    output_embedded_velocity[4] = {-21.81650306,-40.75920676,-0.6557581669,54.90454836,132.1891487,0.1308774929,90.0369547,179.8200581,0.374880674}; // TimeIntegratedQSVMS
    output_embedded_velocity[5] = {0.0475308641975,0.0975308641975,-0.0724722272894,0.0469135802469,0.0969135802469,0.0242780720459,29526.5040828,60780.3309105,0.0315274885768}; // EmbeddedQSVMS
    counter = 0;

    // Test cut element with embedded velocity
    array_1d<double, 3> embedded_vel;
    embedded_vel(0) = 1.0;
    embedded_vel(1) = 2.0;
    embedded_vel(2) = 0.0;

    // Set the nodal EMBEDDED_VELOCITY
    for (auto &r_node : model_part.Nodes()) {
        r_node.SetValue(EMBEDDED_VELOCITY, embedded_vel);
    }

    for (ModelPart::ElementIterator i = model_part.ElementsBegin(); i != model_part.ElementsEnd(); i++) {
        i->Set(SLIP, false);
        i->CalculateLocalSystem(LHS, RHS, r_process_info);

        //std::cout << i->Info() << std::setprecision(12) << std::endl;
        //KRATOS_WATCH(RHS);

        for (unsigned int j = 0; j < RHS.size(); j++) {
            KRATOS_EXPECT_NEAR(RHS[j], output_embedded_velocity[counter][j], 1e-6);
        }

        counter++;
    }

    std::vector< std::vector<double> > output_slip_embedded_velocity(6);
    output_slip_embedded_velocity[0] = {-3.97751921427, 5350.81415787, 0.0428917717671, 5.14836265024, 5335.59090954, 0.11459322177, 28.2302712431, 21412.3110048, 0.40918167313}; // EmbeddedNavierStokes
    output_slip_embedded_velocity[1] = {-6.15363962331,7161.74321338,0.0250586306072,5.99877694774,7156.06872216,0.121191654621,31.5355072101,28683.8110085,0.420416381438}; // EmbeddedFluidElement
    output_slip_embedded_velocity[2] = {7.044040063,26.56103979,-0.5033880823,60.93071102,153.3702507,0.08535256462,95.7560183,186.7293649,0.2680355177}; // NavierStokes
    output_slip_embedded_velocity[3] = {-0.636161784635,8.8199488124,-0.655758245947,67.5798934101,174.543598083,0.130877515367,110.444522985,215.350672279,0.37488073058}; // WeaklyCompressibleNavierStokes2D3N
    output_slip_embedded_velocity[4] = {-21.8165030581,-40.7592067601,-0.655758166948,54.9045483625,132.189148699,0.130877492936,90.0369546956,179.820058061,0.374880674012}; // TimeIntegratedQSVMS
    output_slip_embedded_velocity[5] = {-7.93660736037,7159.32557879,0.0250586369081,6.05364431844,7154.58601556,0.121191652293,25.4246297086,28674.3642131,0.420416377466}; // EmbeddedQSVMS
    counter = 0;

    // Set the nodal EMBEDDED_VELOCITY
    for (auto &r_node : model_part.Nodes()) {
        r_node.SetValue(EMBEDDED_VELOCITY, embedded_vel);
    }

    // Test slip cut element with embedded velocity
    model_part.GetProcessInfo().SetValue(SLIP_LENGTH, 1.0e+08);
    model_part.GetProcessInfo().SetValue(PENALTY_COEFFICIENT, 10.0);
    for (ModelPart::ElementIterator i = model_part.ElementsBegin(); i != model_part.ElementsEnd(); i++) {
        i->Set(SLIP, true);
        i->CalculateLocalSystem(LHS, RHS, r_process_info);

        //std::cout << i->Info() << std::setprecision(12) << std::endl;
        //KRATOS_WATCH(RHS);

        for (unsigned int j = 0; j < RHS.size(); j++) {
            KRATOS_EXPECT_NEAR(RHS[j], output_slip_embedded_velocity[counter][j], 1e-6);
        }

        counter++;
    }
}

}  // namespace Testing
}  // namespace Kratos
