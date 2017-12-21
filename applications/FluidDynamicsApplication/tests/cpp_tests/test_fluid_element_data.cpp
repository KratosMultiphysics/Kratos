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

    p_element->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = 1.0;
    p_element->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = 1.0;
    p_element->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = 1.0;
    /*p_element->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
    p_element->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = -1.0;
    p_element->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) =  0.5;
    array_1d<double, 3> embedded_vel;
    embedded_vel(0) = 1.0;
    embedded_vel(1) = 2.0;
    embedded_vel(2) = 0.0;
    p_element->SetValue(EMBEDDED_VELOCITY, embedded_vel);*/

    // Compute RHS and LHS
    Vector RHS = ZeroVector(9);
    Matrix LHS = ZeroMatrix(9,9);

    for (ModelPart::ElementIterator i = model_part.ElementsBegin(); i != model_part.ElementsEnd(); i++) {
        i->Initialize(); // Initialize the element to initialize the constitutive law
        i->CalculateLocalSystem(LHS, RHS, model_part.GetProcessInfo());

        KRATOS_WATCH(RHS);
    }
/*
    Element::Pointer p_element2 = model_part.pGetElement(2);
    p_element2->Initialize(); // Initialize the element to initialize the constitutive law
    p_element2->CalculateLocalSystem(LHS, RHS, model_part.GetProcessInfo());
    KRATOS_WATCH(RHS);

    Element::Pointer p_element3 = model_part.pGetElement(3);
    p_element3->Initialize(); // Initialize the element to initialize the constitutive law
    p_element3->CalculateLocalSystem(LHS, RHS, model_part.GetProcessInfo());
    KRATOS_WATCH(RHS);

    Element::Pointer p_element4 = model_part.pGetElement(4);
    p_element4->Initialize(); // Initialize the element to initialize the constitutive law
    p_element4->CalculateLocalSystem(LHS, RHS, model_part.GetProcessInfo());
    KRATOS_WATCH(RHS);*/
}

}  // namespace Testing
}  // namespace Kratos