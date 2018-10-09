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

KRATOS_TEST_CASE_IN_SUITE(EmbeddedElement2D3N, FluidDynamicsApplicationFastSuite)
{
    ModelPart model_part("Main");
    model_part.SetBufferSize(3);

    // Variables addition
    model_part.AddNodalSolutionStepVariable(DENSITY); // TODO: To be removed once the element migration is finally finished
    model_part.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY); // TODO: To be removed once the element migration is finally finished
    model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    model_part.AddNodalSolutionStepVariable(DYNAMIC_TAU);
    model_part.AddNodalSolutionStepVariable(SOUND_VELOCITY);
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    model_part.AddNodalSolutionStepVariable(DISTANCE);
    model_part.AddNodalSolutionStepVariable(ACCELERATION);

    // For VMS comparison
    model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    model_part.AddNodalSolutionStepVariable(ADVPROJ);
    model_part.AddNodalSolutionStepVariable(DIVPROJ);


    // Process info creation
    double delta_time = 0.1;
    model_part.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
    model_part.GetProcessInfo().SetValue(SOUND_VELOCITY, 1.0e+3);
    model_part.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
    Vector bdf_coefs(3);
    bdf_coefs[0] = 3.0/(2.0*delta_time);
    bdf_coefs[1] = -2.0/delta_time;
    bdf_coefs[2] = 0.5*delta_time;
    model_part.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);

    // Set the element properties
    Properties::Pointer p_properties = model_part.pGetProperties(0);
    p_properties->SetValue(DENSITY, 1000.0);
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
    model_part.CreateNewElement("EmbeddedSymbolicNavierStokes2D3N", 2, element_nodes, p_properties);
    model_part.CreateNewElement("NavierStokes2D3N", 3, element_nodes, p_properties);
    model_part.CreateNewElement("SymbolicNavierStokes2D3N", 4, element_nodes, p_properties);
    model_part.CreateNewElement("TimeIntegratedQSVMS2D3N", 5, element_nodes, p_properties);
    model_part.CreateNewElement("EmbeddedQSVMS2D3N", 6, element_nodes, p_properties);

    // Define the nodal values
    Matrix reference_velocity(3,2);
    reference_velocity(0,0) = 0.0; reference_velocity(0,1) = 0.1;
    reference_velocity(1,0) = 0.1; reference_velocity(1,1) = 0.2;
    reference_velocity(2,0) = 0.2; reference_velocity(2,1) = 0.3;

    Element::Pointer p_element = model_part.pGetElement(1);

    for(unsigned int i=0; i<3; i++){
        p_element->GetGeometry()[i].FastGetSolutionStepValue(DENSITY) = p_properties->GetValue(DENSITY); // TODO: To be removed once the element migration is finally finished
        p_element->GetGeometry()[i].FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = p_properties->GetValue(DYNAMIC_VISCOSITY); // TODO: To be removed once the element migration is finally finished
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
    output_uncut[0] = {7.044040063,26.56103979,-0.5033880823,60.93071102,153.3702507,0.08535256462,95.7560183,186.7293649,0.2680355177}; // EmbeddedNavierStokes
    output_uncut[1] = {7.044040063,26.56103979,-0.5033880823,60.93071102,153.3702507,0.08535256462,95.7560183,186.7293649,0.2680355177}; // EmbeddedFluidElement
    output_uncut[2] = {7.044040063,26.56103979,-0.5033880823,60.93071102,153.3702507,0.08535256462,95.7560183,186.7293649,0.2680355177}; // NavierStokes
    output_uncut[3] = {7.044040063,26.56103979,-0.5033880823,60.93071102,153.3702507,0.08535256462,95.7560183,186.7293649,0.2680355177}; // SymbolicNavierStokes2D3N
    output_uncut[4] = {-21.81650306,-40.75920676,-0.6557581669,54.90454836,132.1891487,0.1308774929,90.0369547,179.8200581,0.374880674}; // TimeIntegratedQSVMS
    output_uncut[5] = {-21.81650306,-40.75920676,-0.6557581669,54.90454836,132.1891487,0.1308774929,90.0369547,179.8200581,0.374880674}; // EmbeddedQSVMS
    int counter = 0;

    // Test Uncut element
    p_element->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = 1.0;
    p_element->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = 1.0;
    p_element->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = 1.0;

    for (ModelPart::ElementIterator i = model_part.ElementsBegin(); i != model_part.ElementsEnd(); i++) {
        i->Initialize(); // Initialize the element to initialize the constitutive law
        i->Check(model_part.GetProcessInfo()); // Otherwise the constitutive law is not seen here
        i->CalculateLocalSystem(LHS, RHS, model_part.GetProcessInfo());

        //std::cout << i->Info() << std::setprecision(10) << std::endl;
        //KRATOS_WATCH(RHS);

        for (unsigned int j = 0; j < RHS.size(); j++) {
            KRATOS_CHECK_NEAR(RHS[j], output_uncut[counter][j], 1e-6);
        }

        counter++;
    }

    std::vector< std::vector<double> > output_cut(6);
    output_cut[0] = {-0.008024691358,-0.01358024691,-0.05463909243,-0.008641975309,-0.01419753086,0.01767964152,-2867.408637,-4779.580067,0.02029278424}; // EmbeddedNavierStokes
    output_cut[1] = {-0.008024691358,-0.01358024691,-0.05463909243,-0.008641975309,-0.01419753086,0.01767964152,-2867.408636,-4779.580066,0.02029278424}; // EmbeddedFluidElement
    output_cut[2] = {7.044040063,26.56103979,-0.5033880823,60.93071102,153.3702507,0.08535256462,95.7560183,186.7293649,0.2680355177}; // NavierStokes
    output_cut[3] = {7.044040063,26.56103979,-0.5033880823,60.93071102,153.3702507,0.08535256462,95.7560183,186.7293649,0.2680355177}; // SymbolicNavierStokes2D3N
    output_cut[4] = {-21.81650306,-40.75920676,-0.6557581669,54.90454836,132.1891487,0.1308774929,90.0369547,179.8200581,0.374880674}; // TimeIntegratedQSVMS
    output_cut[5] = {-0.008024691358,-0.01358024691,-0.07247222729,-0.008641975309,-0.01419753086,0.02427807205,-5132.699253,-8554.692789,0.03152748858}; // EmbeddedQSVMS
    counter = 0;

    // Test cut element
    p_element->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
    p_element->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = -1.0;
    p_element->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) =  0.5;
    
    for (ModelPart::ElementIterator i = model_part.ElementsBegin(); i != model_part.ElementsEnd(); i++) {
        i->Set(SLIP, false);
        i->Initialize(); // Initialize the element to initialize the constitutive law
        i->CalculateLocalSystem(LHS, RHS, model_part.GetProcessInfo());

        //std::cout << i->Info() << std::setprecision(10) << std::endl;
        //KRATOS_WATCH(RHS);
        
        for (unsigned int j = 0; j < RHS.size(); j++) {
            KRATOS_CHECK_NEAR(RHS[j], output_cut[counter][j], 1e-6);
        }

        counter++;
    }

    std::vector< std::vector<double> > output_slip_cut(6);
    output_slip_cut[0] = {-3.97751477,-753.7076317,-0.06821933934,5.148358206,-768.9308712,0.003482110658,28.23027124,-3005.776109,-0.03526277131}; // EmbeddedNavierStokes
    output_slip_cut[1] = {-3.977514436,-753.7076315,-0.06821933934,5.148358539,-768.930871,0.003482110658,28.23027258,-3005.776108,-0.03526277131}; // EmbeddedFluidElement
    output_slip_cut[2] = {7.044040063,26.56103979,-0.5033880823,60.93071102,153.3702507,0.08535256462,95.7560183,186.7293649,0.2680355177}; // NavierStokes
    output_slip_cut[3] = {7.044040063,26.56103979,-0.5033880823,60.93071102,153.3702507,0.08535256462,95.7560183,186.7293649,0.2680355177}; // SymbolicNavierStokes2D3N
    output_slip_cut[4] = {-21.81650306,-40.75920676,-0.6557581669,54.90454836,132.1891487,0.1308774929,90.0369547,179.8200581,0.374880674}; // TimeIntegratedQSVMS
    output_slip_cut[5] = {-7.936602583,-1004.669448,-0.0860524742,6.053640207,-1023.809413,0.01008054118,25.42463104,-4010.41667,-0.02402806698}; // EmbeddedQSVMS
    counter = 0;

    // Test slip cut element
    for (ModelPart::ElementIterator i = model_part.ElementsBegin(); i != model_part.ElementsEnd(); ++i) {
        i->Set(SLIP, true);
        i->Initialize(); // Initialize the element to initialize the constitutive law
        i->CalculateLocalSystem(LHS, RHS, model_part.GetProcessInfo());

        //std::cout << i->Info() << std::setprecision(10) << std::endl;
        //KRATOS_WATCH(RHS);
        
        for (unsigned int j = 0; j < RHS.size(); j++) {
            KRATOS_CHECK_NEAR(RHS[j], output_slip_cut[counter][j], 1e-6);
        }

        counter++;
    }

    std::vector< std::vector<double> > output_embedded_velocity;
    output_embedded_velocity.resize(6);
    output_embedded_velocity[0] = {0.0475308641975,0.0975308641975,-0.0546390924304,0.0469135802469,0.0969135802469,0.0176796415227,16436.8507492,33828.9387065,0.020292784241}; // EmbeddedNavierStokes
    output_embedded_velocity[1] = {0.0475308641975,0.0975308641975,-0.0546390924304,0.0469135802469,0.0969135802469,0.0176796415227,16436.8507505,33828.9387074,0.020292784241}; // EmbeddedFluidElement
    output_embedded_velocity[2] = {7.044040063,26.56103979,-0.5033880823,60.93071102,153.3702507,0.08535256462,95.7560183,186.7293649,0.2680355177}; // NavierStokes
    output_embedded_velocity[3] = {7.044040063,26.56103979,-0.5033880823,60.93071102,153.3702507,0.08535256462,95.7560183,186.7293649,0.2680355177}; // SymbolicNavierStokes2D3N
    output_embedded_velocity[4] = {-21.81650306,-40.75920676,-0.6557581669,54.90454836,132.1891487,0.1308774929,90.0369547,179.8200581,0.374880674}; // TimeIntegratedQSVMS
    output_embedded_velocity[5] = {0.0475308641975,0.0975308641975,-0.0724722272894,0.0469135802469,0.0969135802469,0.0242780720459,29254.7933008,60220.2923181,0.0315274885768}; // EmbeddedQSVMS
    counter = 0;

    // Test cut element with embedded velocity
    array_1d<double, 3> embedded_vel;
    embedded_vel(0) = 1.0;
    embedded_vel(1) = 2.0;
    embedded_vel(2) = 0.0;
    
    for (ModelPart::ElementIterator i = model_part.ElementsBegin(); i != model_part.ElementsEnd(); i++) {
        i->Set(SLIP, false);
        i->SetValue(EMBEDDED_VELOCITY, embedded_vel);
        i->Initialize(); // Initialize the element to initialize the constitutive law
        i->CalculateLocalSystem(LHS, RHS, model_part.GetProcessInfo());

        //std::cout << i->Info() << std::setprecision(12) << std::endl;
        //KRATOS_WATCH(RHS);
        
        for (unsigned int j = 0; j < RHS.size(); j++) {
            KRATOS_CHECK_NEAR(RHS[j], output_embedded_velocity[counter][j], 1e-6);
        }

        counter++;
    }

    std::vector< std::vector<double> > output_slip_embedded_velocity(6);
    output_slip_embedded_velocity[0] = {-3.97751921427,5350.81415787,0.0428917717671,5.14836265024,5335.59090954,0.11459322177,28.2302712431,21412.3110048,0.40918167313}; // EmbeddedNavierStokes
    output_slip_embedded_velocity[1] = {-3.97751888093,5350.8141581,0.0428917717671,5.14836298358,5335.59090976,0.11459322177,28.2302725764,21412.3110057,0.40918167313}; // EmbeddedFluidElement
    output_slip_embedded_velocity[2] = {7.044040063,26.56103979,-0.5033880823,60.93071102,153.3702507,0.08535256462,95.7560183,186.7293649,0.2680355177}; // NavierStokes
    output_slip_embedded_velocity[3] = {7.044040063,26.56103979,-0.5033880823,60.93071102,153.3702507,0.08535256462,95.7560183,186.7293649,0.2680355177}; // SymbolicNavierStokes2D3N
    output_slip_embedded_velocity[4] = {-21.8165030581,-40.7592067601,-0.655758166948,54.9045483625,132.189148699,0.130877492936,90.0369546956,179.820058061,0.374880674012}; // TimeIntegratedQSVMS
    output_slip_embedded_velocity[5] = {-7.93660702703,7100.52460316,0.0250586369081,6.05364465178,7081.38462929,0.121191652293,25.4246310419,28410.3594893,0.420416377466}; // EmbeddedQSVMS
    counter = 0;

    // Test slip cut element with embedded velocity   
    for (ModelPart::ElementIterator i = model_part.ElementsBegin(); i != model_part.ElementsEnd(); i++) {
        i->Set(SLIP, true);
        i->SetValue(EMBEDDED_VELOCITY, embedded_vel);
        i->Initialize(); // Initialize the element to initialize the constitutive law
        i->CalculateLocalSystem(LHS, RHS, model_part.GetProcessInfo());

        //std::cout << i->Info() << std::setprecision(12) << std::endl;
        //KRATOS_WATCH(RHS);
        
        for (unsigned int j = 0; j < RHS.size(); j++) {
            KRATOS_CHECK_NEAR(RHS[j], output_slip_embedded_velocity[counter][j], 1e-6);
        }

        counter++;
    }
}

}  // namespace Testing
}  // namespace Kratos