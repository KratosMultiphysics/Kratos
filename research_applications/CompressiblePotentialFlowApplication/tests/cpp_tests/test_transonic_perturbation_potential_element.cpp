//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Inigo Lopez
//
//

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "compressible_potential_flow_application_variables.h"
#include "fluid_dynamics_application_variables.h"
#include "custom_elements/transonic_perturbation_potential_flow_element.h"
#include "custom_utilities/potential_flow_utilities.h"
#include "processes/find_nodal_neighbours_process.h"

namespace Kratos {
namespace Testing {

typedef ModelPart::IndexType IndexType;
typedef ModelPart::NodeIterator NodeIteratorType;

void GenerateTransonicPerturbationElement(ModelPart& rModelPart) {
    // Variables addition
    rModelPart.AddNodalSolutionStepVariable(VELOCITY_POTENTIAL);
    rModelPart.AddNodalSolutionStepVariable(AUXILIARY_VELOCITY_POTENTIAL);

    // Set the element properties
    Properties::Pointer pElemProp = rModelPart.CreateNewProperties(0);
    rModelPart.GetProcessInfo()[FREE_STREAM_DENSITY] = 1.225;
    rModelPart.GetProcessInfo()[FREE_STREAM_MACH] = 0.6;
    rModelPart.GetProcessInfo()[HEAT_CAPACITY_RATIO] = 1.4;
    rModelPart.GetProcessInfo()[SOUND_VELOCITY] = 340.3;
    rModelPart.GetProcessInfo()[MACH_LIMIT] = 1.73205080756887729;
    rModelPart.GetProcessInfo()[CRITICAL_MACH] = 0.99;
    rModelPart.GetProcessInfo()[UPWIND_FACTOR_CONSTANT] = 1.0;

    BoundedVector<double, 3> free_stream_velocity = ZeroVector(3);
    free_stream_velocity(0) = rModelPart.GetProcessInfo().GetValue(FREE_STREAM_MACH) * rModelPart.GetProcessInfo().GetValue(SOUND_VELOCITY);
    rModelPart.GetProcessInfo()[FREE_STREAM_VELOCITY] = free_stream_velocity;

    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> elemNodes{1, 2, 3};
    rModelPart.CreateNewElement("TransonicPerturbationPotentialFlowElement2D3N", 1, elemNodes, pElemProp);
}

void GenerateTransonicPerturbationUpwindElement(ModelPart& rModelPart) {
    // Variables addition
    // Set the element properties
    Properties::Pointer pElemProp = rModelPart.CreateNewProperties(1);

    // Geometry creation
    rModelPart.CreateNewNode(4, 0.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> elemNodes{1, 3, 4};
    rModelPart.CreateNewElement("TransonicPerturbationPotentialFlowElement2D3N", 2, elemNodes, pElemProp);
}

void AssignPotentialsToNormalTransonicPerturbationElement(Element::Pointer pElement)
{
    std::array<double, 3> potential{1.0, 100.0, 150.0};

    for (unsigned int i = 0; i < 3; i++)
        pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = potential[i];
}

void AssignPerturbationPotentialsToTransonicElement(Element& rElement, const std::array<double, 3> rPotential) {
    for (unsigned int i = 0; i < 3; i++){
        rElement.GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = rPotential[i];
    }
}

/** Checks the TransonicPerturbationPotentialFlowElement.
 * Checks the RHS computation.
 */
KRATOS_TEST_CASE_IN_SUITE(TransonicPerturbationPotentialFlowElementRHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    p_element->SetFlags(INLET);

    AssignPotentialsToNormalTransonicPerturbationElement(p_element);

    // Compute RHS
    Vector RHS = ZeroVector(4);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateRightHandSide(RHS, r_current_process_info);

    std::vector<double> reference{146.2643261263345,-122.1426284341492,-24.12169769218525};

    KRATOS_EXPECT_VECTOR_NEAR(RHS, reference, 1e-13);
}

/** Checks the TransonicPerturbationPotentialFlowElement.
 * Checks the inlet RHS computation.
 */
KRATOS_TEST_CASE_IN_SUITE(TransonicPerturbationPotentialFlowInletElementRHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->Initialize(r_current_process_info);

    AssignPotentialsToNormalTransonicPerturbationElement(p_element);

    // Compute RHS
    Vector RHS = ZeroVector(3);

    p_element->CalculateRightHandSide(RHS, r_current_process_info);

    std::vector<double> reference{146.2643261263345,-122.1426284341492,-24.12169769218525};

    KRATOS_EXPECT_VECTOR_NEAR(RHS, reference, 1e-13);
}

/** Checks the TransonicPerturbationPotentialFlowElement.
 * Checks the LHS computation.
 */
KRATOS_TEST_CASE_IN_SUITE(TransonicPerturbationPotentialFlowElementLHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    p_element->SetFlags(INLET);

    AssignPotentialsToNormalTransonicPerturbationElement(p_element);

    // Compute LHS
    Matrix LHS = ZeroMatrix(4, 4);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateLeftHandSide(LHS, r_current_process_info);

    std::array<double, 16> reference{ 0.061142784644415527,-0.1306215050744058, 0.06947872042999037, 0.0,
                                     -0.1306215050744058, 0.6710758508914103,-0.5404543458170046, 0.0,
                                      0.06947872042999037,-0.5404543458170046,0.4709756253870142, 0.0,
                                      0.0, 0.0, 0.0, 0.0};

    for (unsigned int i = 0; i < LHS.size1(); i++) {
        for (unsigned int j = 0; j < LHS.size2(); j++) {
            KRATOS_EXPECT_RELATIVE_NEAR(LHS(i, j), reference[i * 4 + j], 1e-15);
        }
    }
}

/** Checks the TransonicPerturbationPotentialFlowElement.
 * Checks the Supersonic LHS computation.
 */
KRATOS_TEST_CASE_IN_SUITE(TransonicPerturbationPotentialFlowSupersonicElementLHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement(model_part);
    GenerateTransonicPerturbationUpwindElement(model_part);

    Element::Pointer p_element = model_part.pGetElement(1);
    Element::Pointer pUpwindElement = model_part.pGetElement(2);

    FindNodalNeighboursProcess find_nodal_neighbours_process(model_part);
    find_nodal_neighbours_process.Execute();

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->Initialize(r_current_process_info);
    pUpwindElement->SetFlags(INLET);

    std::array<double, 3> high_potential{1.0, 200.0, 100.0};  // node id order 23 74 55
    std::array<double, 3> low_potential{1.0, 100.0, 150.0};   // node id order 23 55 67
    // mach number 1.92516
    AssignPerturbationPotentialsToTransonicElement(*p_element, high_potential);
    // mach number 0.39943
    AssignPerturbationPotentialsToTransonicElement(*pUpwindElement, low_potential);

    for (auto& r_node : model_part.Nodes()){
        r_node.AddDof(VELOCITY_POTENTIAL);
    }

    Element::DofsVectorType CurrentElementalDofList;
    p_element->GetDofList(CurrentElementalDofList, r_current_process_info);

    Element::DofsVectorType UpwindElementalDofList;
    pUpwindElement->GetDofList(UpwindElementalDofList, r_current_process_info);

    std::vector<int> current_ids{23, 74, 55};
    std::vector<int> upwind_ids{23, 55, 67};
    for (int i = 0; i < 3; i++) {
        CurrentElementalDofList[i]->SetEquationId(current_ids[i]);
        UpwindElementalDofList[2]->SetEquationId(upwind_ids[i]);
    }

    // // Compute LHS
    Matrix LHS = ZeroMatrix(4, 4);

    p_element->CalculateLeftHandSide(LHS, r_current_process_info);

    std::array<double, 16> reference{ 0.37708651121240516,-0.54915594944726343,0.17743192602659938,-0.0053624877917411414,
                                     -0.35664864389216139,1.0308840557120569,-0.68092794767744869,0.0066925358575532206,
                                     -0.020437867320243799,-0.48172810626479345,0.50349602165084928,-0.0013300480658120792,
                                      0.0,0.0,0.0,0.0};

    for (unsigned int i = 0; i < LHS.size1(); i++) {
        for (unsigned int j = 0; j < LHS.size2(); j++) {
            KRATOS_EXPECT_RELATIVE_NEAR(LHS(i, j), reference[i * 4 + j], 1e-15);
        }
    }
}

/** Checks the TransonicPerturbationPotentialFlowElement.
 * Checks the Supersonic LHS computation.
 */
KRATOS_TEST_CASE_IN_SUITE(TransonicPerturbationPotentialFlowSupersonicDeceleratingElementLHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement(model_part);
    GenerateTransonicPerturbationUpwindElement(model_part);

    Element::Pointer p_element = model_part.pGetElement(1);
    Element::Pointer pUpwindElement = model_part.pGetElement(2);

    FindNodalNeighboursProcess find_nodal_neighbours_process(model_part);
    find_nodal_neighbours_process.Execute();

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->Initialize(r_current_process_info);
    pUpwindElement->SetFlags(INLET);

    std::array<double, 3> high_potential{10.0, 200.0, 350.0};  // node id order 23 74 55
    std::array<double, 3> higher_potential{10.0, 350.0, 150.0};   // node id order 23 55 67
    // mach number 1.92516
    AssignPerturbationPotentialsToTransonicElement(*p_element, high_potential);
    // mach number 2.0898
    AssignPerturbationPotentialsToTransonicElement(*pUpwindElement, higher_potential);


    for (auto& r_node : model_part.Nodes()){
        r_node.AddDof(VELOCITY_POTENTIAL);
    }

    Element::DofsVectorType CurrentElementalDofList;
    p_element->GetDofList(CurrentElementalDofList, r_current_process_info);

    Element::DofsVectorType UpwindElementalDofList;
    pUpwindElement->GetDofList(UpwindElementalDofList, r_current_process_info);

    std::vector<int> current_ids{23, 74, 55};
    std::vector<int> upwind_ids{23, 55, 67};
    for (int i = 0; i < 3; i++) {
        CurrentElementalDofList[i]->SetEquationId(current_ids[i]);
        UpwindElementalDofList[2]->SetEquationId(upwind_ids[i]);
    }

    // // Compute LHS
    Matrix LHS = ZeroMatrix(4, 4);

    p_element->CalculateLeftHandSide(LHS, r_current_process_info);

    std::array<double, 16> reference{ -0.054669771246690091,-0.14888540383260493,0.40255671949894734,-0.19900154441965234,
                                      -0.083557300280464597,0.51822478070562072,-0.55794161182804491,0.12327413140288881,
                                       0.13822707152715466,-0.36933937687301577,0.15538489232909758,0.075727413016763542,
                                       0.0,0.0,0.0,0.0};

    for (unsigned int i = 0; i < LHS.size1(); i++) {
        for (unsigned int j = 0; j < LHS.size2(); j++) {
            KRATOS_EXPECT_RELATIVE_NEAR(LHS(i, j), reference[i * 4 + j], 1e-15);
        }
    }
}

/** Checks the TransonicPerturbationPotentialFlowElement.
 * Checks the Supersonic RHS computation.
 */
KRATOS_TEST_CASE_IN_SUITE(TransonicPerturbationPotentialFlowSupersonicElementRHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement(model_part);
    GenerateTransonicPerturbationUpwindElement(model_part);

    Element::Pointer p_element = model_part.pGetElement(1);
    Element::Pointer pUpwindElement = model_part.pGetElement(2);

    FindNodalNeighboursProcess find_nodal_neighbours_process(model_part);
    find_nodal_neighbours_process.Execute();

    std::array<double, 3> high_potential{1.0, 200.0, 100.0};  // node id order 23 74 55
    std::array<double, 3> low_potential{1.0, 100.0, 150.0};   // node id order 23 55 67
    // mach number 1.92516
    AssignPerturbationPotentialsToTransonicElement(*p_element, high_potential);
    // mach number 0.3999
    AssignPerturbationPotentialsToTransonicElement(*pUpwindElement, low_potential);


    for (auto& r_node : model_part.Nodes()){
        r_node.AddDof(VELOCITY_POTENTIAL);
    }

    Element::DofsVectorType CurrentElementalDofList;
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->GetDofList(CurrentElementalDofList, r_current_process_info);

    Element::DofsVectorType UpwindElementalDofList;
    pUpwindElement->GetDofList(UpwindElementalDofList, r_current_process_info);

    std::vector<int> current_ids{23, 74, 55};
    std::vector<int> upwind_ids{23, 55, 67};
    for (int i = 0; i < 3; i++) {
        CurrentElementalDofList[i]->SetEquationId(current_ids[i]);
        UpwindElementalDofList[i]->SetEquationId(upwind_ids[i]);
    }

    p_element->Initialize(r_current_process_info);
    pUpwindElement->SetFlags(INLET);

    // Compute RHS
    Vector RHS = ZeroVector(4);

    p_element->CalculateRightHandSide(RHS, r_current_process_info);

    std::vector<double> reference{185.25633340652948,-231.20512387394589,45.948790467416408,0.0};

    KRATOS_EXPECT_VECTOR_NEAR(RHS, reference, 1e-15);
}

/** Checks the TransonicPerturbationPotentialFlowElement.
 * Checks the LHS inlet computation.
 */
KRATOS_TEST_CASE_IN_SUITE(TransonicPerturbationPotentialFlowInletElementLHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->Initialize(r_current_process_info);

    p_element->SetFlags(INLET);

    AssignPotentialsToNormalTransonicPerturbationElement(p_element);

    // Compute LHS
    Matrix LHS = ZeroMatrix(4, 4);

    p_element->CalculateLeftHandSide(LHS, r_current_process_info);

    std::array<double, 9> reference{ 0.061142784644415527,-0.1306215050744058, 0.06947872042999037,
                                     -0.1306215050744058, 0.67107585089141042,-0.5404543458170046,
                                      0.06947872042999037,-0.5404543458170046,0.4709756253870142};

    for (unsigned int i = 0; i < LHS.size1(); i++) {
        for (unsigned int j = 0; j < LHS.size2(); j++) {
            KRATOS_EXPECT_NEAR(LHS(i, j), reference[i * 3 + j], 1e-16);
        }
    }
}

/** Checks the TransonicPerturbationPotentialFlowElement.
 * Tests the LHS computation.
 */
KRATOS_TEST_CASE_IN_SUITE(PingTransonicPerturbationPotentialFlowElementLHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    p_element->SetFlags(INLET);

    AssignPotentialsToNormalTransonicPerturbationElement(p_element);

    Vector RHS_original = ZeroVector(number_of_nodes);
    Matrix LHS_original = ZeroMatrix(number_of_nodes, number_of_nodes);
    Matrix LHS_finite_diference = ZeroMatrix(number_of_nodes, number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(number_of_nodes, number_of_nodes);

    // Compute original RHS and LHS
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateLocalSystem(LHS_original, RHS_original, r_current_process_info);

    double delta = 1e-3;
    for(unsigned int i = 0; i < number_of_nodes; i++){
        // Pinging
        p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += delta;

        Vector RHS_pinged = ZeroVector(number_of_nodes);
        Matrix LHS_pinged = ZeroMatrix(number_of_nodes, number_of_nodes);
        // Compute pinged LHS and RHS
        p_element->CalculateLocalSystem(LHS_pinged, RHS_pinged, r_current_process_info);

        for(unsigned int k = 0; k < number_of_nodes; k++){
            // Compute the finite difference estimate of the sensitivity
            LHS_finite_diference( k, i) = -(RHS_pinged(k)-RHS_original(k)) / delta;
            // Compute the average of the original and pinged analytic sensitivities
            LHS_analytical( k, i) = 0.5 * (LHS_original(k,i) + LHS_pinged(k,i));
        }

        // Unpinging
        p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= delta;
    }

    KRATOS_EXPECT_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

/** Checks the TransonicPerturbationPotentialFlowElement.
 * Tests the LHS computation.
 */
KRATOS_TEST_CASE_IN_SUITE(PingTransonicPerturbationPotentialFlowSupersonicElementLHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement(model_part);
    GenerateTransonicPerturbationUpwindElement(model_part);

    Element::Pointer p_element = model_part.pGetElement(1);
    Element::Pointer pUpwindElement = model_part.pGetElement(2);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    FindNodalNeighboursProcess find_nodal_neighbours_process(model_part);
    find_nodal_neighbours_process.Execute();

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->Initialize(r_current_process_info);
    pUpwindElement->SetFlags(INLET);

    std::array<double, 3> high_potential{1.0, 200.0, 100.0};  // node id order 23 74 55
    std::array<double, 3> low_potential{1.0, 100.0, 150.0};   // node id order 23 55 67
    // mach number 1.92516
    AssignPerturbationPotentialsToTransonicElement(*p_element, high_potential);
    // mach number 0.39943
    AssignPerturbationPotentialsToTransonicElement(*pUpwindElement, low_potential);

    for (auto& r_node : model_part.Nodes()){
        r_node.AddDof(VELOCITY_POTENTIAL);
    }

    Element::DofsVectorType CurrentElementalDofList;
    p_element->GetDofList(CurrentElementalDofList, r_current_process_info);

    Element::DofsVectorType UpwindElementalDofList;
    pUpwindElement->GetDofList(UpwindElementalDofList, r_current_process_info);

    std::vector<int> current_ids{23, 74, 55}; // 1 2 3
    std::vector<int> upwind_ids{23, 55, 67};  // 1 3 4
    for (int i = 0; i < 3; i++) {
        CurrentElementalDofList[i]->SetEquationId(current_ids[i]);
        UpwindElementalDofList[i]->SetEquationId(upwind_ids[i]);
    }

    Vector RHS_original = ZeroVector(number_of_nodes + 1);
    Matrix LHS_original = ZeroMatrix(number_of_nodes + 1, number_of_nodes + 1);
    Matrix LHS_finite_diference = ZeroMatrix(number_of_nodes + 1, number_of_nodes + 1);
    Matrix LHS_analytical = ZeroMatrix(number_of_nodes + 1, number_of_nodes + 1);

    // Compute original RHS and LHS
    p_element->CalculateLocalSystem(LHS_original, RHS_original, r_current_process_info);

    double delta = 1e-3;
    for(unsigned int i = 0; i < 4; i++){
        // Pinging
        if (i < 3) {
            p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += delta;
        }
        else {
            pUpwindElement->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += delta;
        }

        Vector RHS_pinged = ZeroVector(number_of_nodes + 1);
        Matrix LHS_pinged = ZeroMatrix(number_of_nodes + 1, number_of_nodes + 1);
        // Compute pinged LHS and RHS
        p_element->CalculateLocalSystem(LHS_pinged, RHS_pinged, r_current_process_info);

        for(unsigned int k = 0; k < number_of_nodes + 1; k++){
            // Compute the finite difference estimate of the sensitivity
            LHS_finite_diference( k, i) = -(RHS_pinged(k)-RHS_original(k)) / delta;
            // Compute the average of the original and pinged analytic sensitivities
            LHS_analytical( k, i) = 0.5 * (LHS_original(k,i) + LHS_pinged(k,i));
        }

        // Unpinging
        if (i < 3) {
            p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= delta;
        }
        else {
            pUpwindElement->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= delta;
        }
    }

    KRATOS_EXPECT_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(TransonicPerturbationPotentialFlowElementEquationId, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement(model_part);
    GenerateTransonicPerturbationUpwindElement(model_part);

    FindNodalNeighboursProcess find_nodal_neighbours_process(model_part);
    find_nodal_neighbours_process.Execute();

    Element::Pointer pCurrentElement = model_part.pGetElement(1);
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    pCurrentElement->Initialize(r_current_process_info);

    for (auto& r_node : model_part.Nodes()){
        r_node.AddDof(VELOCITY_POTENTIAL);
    }

    Element::DofsVectorType CurrentElementalDofList;
    pCurrentElement->GetDofList(CurrentElementalDofList, r_current_process_info);

    std::vector<int> ids{23, 74, 55};
    for (int i = 0; i < 3; i++) {
        CurrentElementalDofList[i]->SetEquationId(ids[i]);
    }

    // upwind element equation id
    Element::Pointer pUpwindElement = model_part.pGetElement(2);

    pUpwindElement->GetGeometry()[2].AddDof(VELOCITY_POTENTIAL);

    Element::DofsVectorType UpwindElementalDofList;
    pUpwindElement->GetDofList(UpwindElementalDofList, r_current_process_info);

    UpwindElementalDofList[2]->SetEquationId(67);

    // make and check equation ids
    Element::EquationIdVectorType EquationIdVector;
    pCurrentElement->EquationIdVector(EquationIdVector, r_current_process_info);

    std::vector<double> reference{23.0, 74.0, 55.0, 67.0};
    KRATOS_EXPECT_VECTOR_NEAR(EquationIdVector, reference, 1e-15);
}

BoundedVector<double,3> AssignDistancesToPerturbationTransonicElement()
{
    BoundedVector<double,3> distances;
    distances(0) = 1.0;
    distances(1) = -1.0;
    distances(2) = -1.0;
    return distances;
}

void AssignPotentialsToWakeTransonicPerturbationElement(Element::Pointer pElement, const array_1d<double, 3>& rDistances, const std::array<double, 6>& rPotential)
{
    for (unsigned int i = 0; i < 3; i++){
        if (rDistances(i) > 0.0)
            pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = rPotential[i];
        else
            pElement->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) = rPotential[i];
    }
    for (unsigned int i = 0; i < 3; i++){
        if (rDistances(i) < 0.0)
            pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = rPotential[i+3];
        else
            pElement->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) = rPotential[i+3];
    }
}

KRATOS_TEST_CASE_IN_SUITE(WakeTransonicPerturbationPotentialFlowElementRHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    BoundedVector<double,3> distances = AssignDistancesToPerturbationTransonicElement();

    p_element->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    p_element->GetValue(WAKE) = true;

    const std::array<double, 6> potential{1.0, 101.0, 150.0, 6.0, 105.0, 155.0};
    AssignPotentialsToWakeTransonicPerturbationElement(p_element, distances, potential);

    // Compute RHS and LHS
    Vector RHS = ZeroVector(6);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateRightHandSide(RHS, r_current_process_info);

    std::vector<double> reference{146.392649744264,-1.225,0.6125,-0.6125,-122.1426284341492,-24.12169769218525};

    KRATOS_EXPECT_VECTOR_NEAR(RHS, reference, 1e-13);
}

KRATOS_TEST_CASE_IN_SUITE(WakeTransonicPerturbationPotentialFlowElementLHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    BoundedVector<double,3> distances = AssignDistancesToPerturbationTransonicElement();

    p_element->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    p_element->GetValue(WAKE) = true;

    const std::array<double, 6> potential{1.0, 100.0, 150.0, 6.0, 105.0, 155.0};
    AssignPotentialsToWakeTransonicPerturbationElement(p_element, distances, potential);

    // Compute LHS
    Matrix LHS = ZeroMatrix(6, 6);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateLeftHandSide(LHS, r_current_process_info);

    // Check the LHS values
    std::array<double,36> reference{0.06114278464441553,-0.1306215050744058,
    0.06947872042999036,0,0,0,-0.6125,1.225,-0.6125,0.6125,-1.225,0.6125,0,
    -0.6125,0.6125,-0,0.6125,-0.6125,-0.6125,0.6125,-0,0.6125,-0.6125,0,0,0,
    0,-0.1306215050744058,0.6710758508914104,-0.5404543458170046,0,0,0,
    0.06947872042999036,-0.5404543458170046,0.4709756253870142};

    for (unsigned int i = 0; i < LHS.size1(); i++) {
        for (unsigned int j = 0; j < LHS.size2(); j++) {
            KRATOS_EXPECT_NEAR(LHS(i, j), reference[6 * i + j], 1e-16);
        }
    }
}

} // namespace Testing
} // namespace Kratos.
