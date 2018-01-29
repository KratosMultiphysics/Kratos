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

#include <iomanip> // for std::setprecision

// Project includes
#include "testing/testing.h"
#include "includes/model_part.h"
#include "includes/cfd_variables.h"

#include "custom_utilities/fluid_element_data.h"
#include "custom_constitutive/newtonian_2d_law.h"

namespace Kratos {
namespace Testing {

// Auxiliary classes to test FluidData member functions
class TestNodalVectorData : public FluidElementData<2, 3, true> {
public:
    NodalVectorData Velocity;
    NodalVectorData Velocity_OldStep1;

    void Initialize(
        const Element& rElement, const ProcessInfo& rProcessInfo) override {
        this->FillFromNodalData(Velocity, VELOCITY, rElement.GetGeometry());
        this->FillFromHistoricalNodalData(
            Velocity_OldStep1, VELOCITY, rElement.GetGeometry(), 1);
    }

    static int Check(const Element& rElement, const ProcessInfo& rProcessInfo) {
        const Geometry<Node<3> >& r_geometry = rElement.GetGeometry();
        KRATOS_CHECK_VARIABLE_KEY(VELOCITY);

        for (unsigned int i = 0; i < 3; i++) {
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_geometry[i]);
        }

        return 0;
    }
};

class TestNodalScalarData : public FluidElementData<2, 3, true> {
public:
    NodalScalarData Pressure;
    NodalScalarData Pressure_OldStep1;

    void Initialize(
        const Element& rElement, const ProcessInfo& rProcessInfo) override {
        this->FillFromNodalData(Pressure, PRESSURE, rElement.GetGeometry());
        this->FillFromHistoricalNodalData(
            Pressure_OldStep1, PRESSURE, rElement.GetGeometry(), 1);
    }

    static int Check(const Element& rElement, const ProcessInfo& rProcessInfo) {
        const Geometry<Node<3> >& r_geometry = rElement.GetGeometry();
        KRATOS_CHECK_VARIABLE_KEY(PRESSURE);

        for (unsigned int i = 0; i < 3; i++) {
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE, r_geometry[i]);
        }

        return 0;
    }
};

class TestElementData : public FluidElementData<2, 3, true> {
public:
    double CSmagorinsky;

    void Initialize(
        const Element& rElement, const ProcessInfo& rProcessInfo) override {
        this->FillFromElementData(CSmagorinsky, C_SMAGORINSKY, rElement);
    }

    static int Check(const Element& rElement, const ProcessInfo& rProcessInfo) {
        KRATOS_CHECK_VARIABLE_KEY(C_SMAGORINSKY);
        return 0;
    }
};

class TestPropertiesData : public FluidElementData<2, 3, true> {
public:
    double KinematicViscosity;

    void Initialize(
        const Element& rElement, const ProcessInfo& rProcessInfo) override {
        this->FillFromProperties(
            KinematicViscosity, KINEMATIC_VISCOSITY, rElement);
    }

    static int Check(const Element& rElement, const ProcessInfo& rProcessInfo) {
        KRATOS_CHECK_VARIABLE_KEY(KINEMATIC_VISCOSITY);
        return 0;
    }
};

class TestProcessInfoData : public FluidElementData<2, 3, true> {
public:
    double DeltaTime;
    int UseOSS;

    void Initialize(
        const Element& rElement, const ProcessInfo& rProcessInfo) override {
        this->FillFromProcessInfo(DeltaTime, DELTA_TIME, rProcessInfo);
        this->FillFromProcessInfo(UseOSS, OSS_SWITCH, rProcessInfo);
    }

    static int Check(const Element& rElement, const ProcessInfo& rProcessInfo) {
        KRATOS_CHECK_VARIABLE_KEY(DELTA_TIME);
        KRATOS_CHECK_VARIABLE_KEY(OSS_SWITCH);
        return 0;
    }
};

void FluidElementDataTestEmptyModelPart(
    ModelPart& rModelPart, double DeltaTime, unsigned int BufferSize) {
    rModelPart.SetBufferSize(BufferSize);
    Properties::Pointer p_properties = rModelPart.pGetProperties(0);

    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> element_nodes{1, 2, 3};
    rModelPart.CreateNewElement("QSVMS2D3N", 1, element_nodes, p_properties);

    // Loop starts at 1 because you need one less clone than time steps (JC)
    for (unsigned int i = 1; i < BufferSize; i++) {
        rModelPart.CloneTimeStep(i * DeltaTime);
    }
}

void FluidElementDataTestCompleteModelPart(
    ModelPart& rModelPart, double DeltaTime, unsigned int BufferSize) {
    
    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(PRESSURE);
    rModelPart.SetBufferSize(BufferSize);
    Properties::Pointer p_properties = rModelPart.pGetProperties(0);

    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> element_nodes{1, 2, 3};
    rModelPart.CreateNewElement("QSVMS2D3N", 1, element_nodes, p_properties);

    // Nodal data
    Element& r_element = *(rModelPart.ElementsBegin());
    Geometry< Node<3> >& r_geometry = r_element.GetGeometry();

    for (unsigned int i = 0; i < 3; i++) {
        Node<3>& r_node = r_geometry[i];
        r_node.FastGetSolutionStepValue(PRESSURE) = 10.0 * r_node.Id();
        r_node.FastGetSolutionStepValue(VELOCITY_X) = r_node.Id() + 5.0;
    }

    // Element data
    r_element.SetValue(C_SMAGORINSKY,0.16);
    
    // ProcessInfo
    rModelPart.GetProcessInfo().SetValue(OSS_SWITCH,1);

    // Loop starts at 1 because you need one less clone than time steps
    for (unsigned int i = 1; i < BufferSize; i++) {
        rModelPart.CloneTimeStep(i * DeltaTime);

        for (unsigned int j = 0; j < 3; j++) {
            Node<3>& r_node = r_geometry[j];
            r_node.FastGetSolutionStepValue(PRESSURE) += float(i);
            r_node.FastGetSolutionStepValue(VELOCITY_Y) = r_node.Id() + i;
        }
    }
}
    
KRATOS_TEST_CASE_IN_SUITE(FluidElementDataRead, FluidDynamicsApplicationFastSuite) {

    TestNodalScalarData nodal_scalar_data;
    TestNodalVectorData nodal_vector_data;
    TestElementData element_data;
    TestPropertiesData properties_data;
    TestProcessInfoData process_info_data;

    ModelPart full_model_part("Test Full");

    constexpr double DeltaTime = 0.1;
    FluidElementDataTestCompleteModelPart(full_model_part,DeltaTime,2);
    Element& r_element = *(full_model_part.ElementsBegin());
    ProcessInfo& r_process_info = full_model_part.GetProcessInfo();

    nodal_scalar_data.Initialize(r_element,r_process_info);
    KRATOS_CHECK_EQUAL(nodal_scalar_data.Pressure[0], 11.0);
    KRATOS_CHECK_EQUAL(nodal_scalar_data.Pressure[1], 21.0);
    KRATOS_CHECK_EQUAL(nodal_scalar_data.Pressure[2], 31.0);
    KRATOS_CHECK_EQUAL(nodal_scalar_data.Pressure_OldStep1[0], 10.0);
    KRATOS_CHECK_EQUAL(nodal_scalar_data.Pressure_OldStep1[1], 20.0);
    KRATOS_CHECK_EQUAL(nodal_scalar_data.Pressure_OldStep1[2], 30.0);

    nodal_vector_data.Initialize(r_element,r_process_info);
    KRATOS_CHECK_EQUAL(nodal_vector_data.Velocity(0,1), 2.0);
    KRATOS_CHECK_EQUAL(nodal_vector_data.Velocity(1,1), 3.0);
    KRATOS_CHECK_EQUAL(nodal_vector_data.Velocity(2,1), 4.0);
    KRATOS_CHECK_EQUAL(nodal_vector_data.Velocity_OldStep1(0,0), 6.0);
    KRATOS_CHECK_EQUAL(nodal_vector_data.Velocity_OldStep1(1,0), 7.0);
    KRATOS_CHECK_EQUAL(nodal_vector_data.Velocity_OldStep1(2,0), 8.0);

    element_data.Initialize(r_element,r_process_info);
    KRATOS_CHECK_EQUAL(element_data.CSmagorinsky, r_element.GetValue(C_SMAGORINSKY));

    properties_data.Initialize(r_element,r_process_info);
    KRATOS_CHECK_EQUAL(properties_data.KinematicViscosity, r_element.GetProperties().GetValue(KINEMATIC_VISCOSITY));

    process_info_data.Initialize(r_element,r_process_info);
    KRATOS_CHECK_EQUAL(process_info_data.UseOSS, r_process_info.GetValue(OSS_SWITCH));
    KRATOS_CHECK_EQUAL(process_info_data.DeltaTime, r_process_info.GetValue(DELTA_TIME));
}

KRATOS_TEST_CASE_IN_SUITE(FluidElementDataCheck, FluidDynamicsApplicationFastSuite) {

    ModelPart empty_model_part("Test Empty");

    constexpr double DeltaTime = 0.1;
    FluidElementDataTestEmptyModelPart(empty_model_part,DeltaTime,1);

    Element& r_element = *(empty_model_part.ElementsBegin());
    ProcessInfo& r_process_info = empty_model_part.GetProcessInfo();

    // historical data container should not work with variables not added to model part
    int out;
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        out = TestNodalScalarData::Check(r_element, r_process_info),
        "Missing PRESSURE variable in solution step data for node 1.");
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
    out = TestNodalVectorData::Check(r_element, r_process_info),
        "Missing VELOCITY variable in solution step data for node 1.");

    // Other containers can work with non-initialized variables, but should return 0 
    out = TestElementData::Check(r_element,r_process_info);
    KRATOS_CHECK_EQUAL(out,0);
    out = TestPropertiesData::Check(r_element,r_process_info);
    KRATOS_CHECK_EQUAL(out,0);
    out = TestProcessInfoData::Check(r_element,r_process_info);
    KRATOS_CHECK_EQUAL(out,0);

    TestElementData element_data;
    TestPropertiesData properties_data;
    TestProcessInfoData process_info_data;
    
    element_data.Initialize(r_element,r_process_info);
    KRATOS_CHECK_EQUAL(element_data.CSmagorinsky, 0.0);

    properties_data.Initialize(r_element,r_process_info);
    KRATOS_CHECK_EQUAL(properties_data.KinematicViscosity, 0.0);

    process_info_data.Initialize(r_element,r_process_info);
    KRATOS_CHECK_EQUAL(process_info_data.UseOSS, 0.0);
    KRATOS_CHECK_EQUAL(process_info_data.DeltaTime, 0.0);
}

KRATOS_TEST_CASE_IN_SUITE(EmbeddedElement2D3N, FluidDynamicsApplicationFastSuite)
{
    ModelPart model_part("Main");
    model_part.SetBufferSize(3);

    // Variables addition
    model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    model_part.AddNodalSolutionStepVariable(DENSITY);
    model_part.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
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

    // Define the nodal values
    Matrix reference_velocity(3,2);
    reference_velocity(0,0) = 0.0; reference_velocity(0,1) = 0.1;
    reference_velocity(1,0) = 0.1; reference_velocity(1,1) = 0.2;
    reference_velocity(2,0) = 0.2; reference_velocity(2,1) = 0.3;

    Element::Pointer p_element = model_part.pGetElement(1);

    // Set the nodal DENSITY and DYNAMIC_VISCOSITY values
    for (ModelPart::NodeIterator it_node=model_part.NodesBegin(); it_node<model_part.NodesEnd(); ++it_node){
        it_node->FastGetSolutionStepValue(DENSITY) = p_properties->GetValue(DENSITY);
        it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = p_properties->GetValue(DYNAMIC_VISCOSITY);
    }

    for(unsigned int i=0; i<3; i++){
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

    std::vector< std::vector<double> > output_uncut;
    output_uncut.resize(5);
    output_uncut[0] = {-1.818177294,17.69882244,-0.5033880823,69.79292838,153.3702507,0.08535256462,95.7560183,195.5915822,0.2680355177}; // EmbeddedNavierStokes
    output_uncut[1] = {-12.52607405,-3.069963453,-0.6557582459,79.46980568,174.5435981,0.1308775154,110.444523,227.2405845,0.3748807306}; // EmbeddedFluidElement
    output_uncut[2] = {-1.818177294,17.69882244,-0.5033880823,69.79292838,153.3702507,0.08535256462,95.7560183,195.5915822,0.2680355177}; // NavierStokes
    output_uncut[3] = {-12.52607405,-3.069963453,-0.6557582459,79.46980568,174.5435981,0.1308775154,110.444523,227.2405845,0.3748807306}; // SymbolicNavierStokes2D3N
    output_uncut[4] = {-21.81650306,-40.75920676,-0.6557581669,54.90454836,132.1891487,0.1308774929,90.0369547,179.8200581,0.374880674}; // TimeIntegratedQSVMS
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

    std::vector< std::vector<double> > output_cut;
    output_cut.resize(5);
    output_cut[0] = {-0.008024691358,-0.01358024691,-0.05463909243,-0.008641975309,-0.01419753086,0.01767964152,-2867.408637,-4778.200165,0.02029278424}; // EmbeddedNavierStokes
    output_cut[1] = {-0.008024691358,-0.01358024691,-0.07247223359,-0.008641975309,-0.01419753086,0.02427807437,-5126.588376,-8543.394662,0.03152749255}; // EmbeddedFluidElement
    output_cut[2] = {-1.818177294,17.69882244,-0.5033880823,69.79292838,153.3702507,0.08535256462,95.7560183,195.5915822,0.2680355177}; // NavierStokes
    output_cut[3] = {-12.52607405,-3.069963453,-0.6557582459,79.46980568,174.5435981,0.1308775154,110.444523,227.2405845,0.3748807306}; // SymbolicNavierStokes2D3N
    output_cut[4] = {-21.81650306,-40.75920676,-0.6557581669,54.90454836,132.1891487,0.1308774929,90.0369547,179.8200581,0.374880674}; // TimeIntegratedQSVMS
    counter = 0;

    // Test cut element
    p_element->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
    p_element->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = -1.0;
    p_element->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) =  0.5;
    
    for (ModelPart::ElementIterator i = model_part.ElementsBegin(); i != model_part.ElementsEnd(); i++) {
        i->Initialize(); // Initialize the element to initialize the constitutive law
        i->CalculateLocalSystem(LHS, RHS, model_part.GetProcessInfo());

        //std::cout << i->Info() << std::setprecision(10) << std::endl;
        //KRATOS_WATCH(RHS);
        
        for (unsigned int j = 0; j < RHS.size(); j++) {
            KRATOS_CHECK_NEAR(RHS[j], output_cut[counter][j], 1e-6);
        }

        counter++;
    }

    std::vector< std::vector<double> > output_embedded_velocity;
    output_embedded_velocity.resize(5);
    output_embedded_velocity[0] = {0.0475308641975,0.0975308641975,-0.0546390924304,0.0469135802469,0.0969135802469,0.0176796415227,16436.8507492,33830.318608,0.020292784241}; // EmbeddedNavierStokes
    output_embedded_velocity[1] = {0.0475308641975,0.0975308641975,-0.0724722335903,0.0469135802469,0.0969135802469,0.0242780743742,29260.904177,60231.5904448,0.0315274925494}; // EmbeddedFluidElement
    output_embedded_velocity[2] = {-1.818177294,17.69882244,-0.5033880823,69.79292838,153.3702507,0.08535256462,95.7560183,195.5915822,0.2680355177}; // NavierStokes
    output_embedded_velocity[3] = {-12.52607405,-3.069963453,-0.6557582459,79.46980568,174.5435981,0.1308775154,110.444523,227.2405845,0.3748807306}; // SymbolicNavierStokes2D3N
    output_embedded_velocity[4] = {-21.81650306,-40.75920676,-0.6557581669,54.90454836,132.1891487,0.1308774929,90.0369547,179.8200581,0.374880674}; // TimeIntegratedQSVMS
    counter = 0;

    // Test cut element with embedded velocity
    array_1d<double, 3> embedded_vel;
    embedded_vel(0) = 1.0;
    embedded_vel(1) = 2.0;
    embedded_vel(2) = 0.0;
    
    for (ModelPart::ElementIterator i = model_part.ElementsBegin(); i != model_part.ElementsEnd(); i++) {
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
}

}  // namespace Testing
}  // namespace Kratos