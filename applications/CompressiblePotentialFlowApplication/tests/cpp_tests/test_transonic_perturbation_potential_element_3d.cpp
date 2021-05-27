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
#include "tests/cpp_tests/test_utilities.h"

namespace Kratos {
namespace Testing {

typedef ModelPart::IndexType IndexType;
typedef ModelPart::NodeIterator NodeIteratorType;

void GenerateTransonicPerturbationElement3D(ModelPart& rModelPart) {
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
    free_stream_velocity(0) = rModelPart.GetProcessInfo().GetValue(FREE_STREAM_MACH) *
                              rModelPart.GetProcessInfo().GetValue(SOUND_VELOCITY);
    rModelPart.GetProcessInfo()[FREE_STREAM_VELOCITY] = free_stream_velocity;

    BoundedVector<double, 3> free_stream_velocity_direction = ZeroVector(3);
    free_stream_velocity_direction(0) = 1.0;
    rModelPart.GetProcessInfo()[FREE_STREAM_VELOCITY_DIRECTION] = free_stream_velocity_direction;

    BoundedVector<double, 3> wake_normal = ZeroVector(3);
    wake_normal(2) = 1.0;
    rModelPart.GetProcessInfo()[WAKE_NORMAL] = wake_normal;

    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, -0.2, -0.2);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 0.1, 1.0, 0.0);
    rModelPart.CreateNewNode(4, -0.1, 0.0, 1.0);
    std::vector<ModelPart::IndexType> elemNodes{1, 2, 3, 4};
    rModelPart.CreateNewElement(
        "TransonicPerturbationPotentialFlowElement3D4N", 1, elemNodes, pElemProp);
}

void GenerateTransonicPerturbationUpwindElement3D(ModelPart& rModelPart) {
    // Variables addition
    // Set the element properties
    Properties::Pointer pElemProp = rModelPart.CreateNewProperties(1);

    // Geometry creation
    rModelPart.CreateNewNode(5, -1.0, 0.0, 0.0);
    std::vector<ModelPart::IndexType> elemNodes{5, 1, 3, 4};
    rModelPart.CreateNewElement("TransonicPerturbationPotentialFlowElement3D4N", 2, elemNodes, pElemProp);
}

void AssignPerturbationPotentialsToTransonicElement3D(Element& rElement, const std::array<double, 4> rPotential) {
    for (unsigned int i = 0; i < 4; i++){
        rElement.GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = rPotential[i];
    }
}

void ComputeElementalSensitivitiesTransonic3D(Matrix& rLHS_finite_diference, Matrix& rLHS_analytical, const std::array<double, 4> rPotential, const std::array<double, 4> rPotentialUpwind){
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement3D(model_part);
    GenerateTransonicPerturbationUpwindElement3D(model_part);

    Element::Pointer pElement = model_part.pGetElement(1);
    Element::Pointer pUpwindElement = model_part.pGetElement(2);
    const unsigned int number_of_nodes = pElement->GetGeometry().size();

    FindNodalNeighboursProcess find_nodal_neighbours_process(model_part);
    find_nodal_neighbours_process.Execute();

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    pElement->Initialize(r_current_process_info);
    pUpwindElement->SetFlags(INLET);

    AssignPerturbationPotentialsToTransonicElement3D(*pElement, rPotential);
    AssignPerturbationPotentialsToTransonicElement3D(*pUpwindElement, rPotentialUpwind);

    for (auto& r_node : model_part.Nodes()){
        r_node.AddDof(VELOCITY_POTENTIAL);
    }

    Element::DofsVectorType CurrentElementalDofList;
    pElement->GetDofList(CurrentElementalDofList, r_current_process_info);

    Element::DofsVectorType UpwindElementalDofList;
    pUpwindElement->GetDofList(UpwindElementalDofList, r_current_process_info);

    std::vector<int> current_ids{23, 74, 55, 35}; // 1 2 3 4
    std::vector<int> upwind_ids{87, 23, 55, 35};  // 5 1 3 4
    for (unsigned int i = 0; i < number_of_nodes; i++) {
        CurrentElementalDofList[i]->SetEquationId(current_ids[i]);
        UpwindElementalDofList[i]->SetEquationId(upwind_ids[i]);
    }

    // Compute original RHS and LHS
    Vector RHS_original = ZeroVector(number_of_nodes);
    Matrix LHS_original = ZeroMatrix(number_of_nodes, number_of_nodes);
    pElement->CalculateLocalSystem(LHS_original, RHS_original, r_current_process_info);

    double delta = 1e-3;
    for(unsigned int i = 0; i < number_of_nodes+1; i++){
        // Pinging
        if (i < number_of_nodes) {
            pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += delta;
        }
        else {
            pUpwindElement->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += delta;
        }

        PotentialFlowTestUtilities::ComputeElementalSensitivitiesMatrixRow(model_part, delta, i, LHS_original, RHS_original, rLHS_finite_diference, rLHS_analytical);

        // Unpinging
        if (i < number_of_nodes) {
            pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= delta;
        }
        else {
            pUpwindElement->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= delta;
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(PingTransonicPerturbationPotentialFlowElementLHS3D,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement3D(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    p_element->SetFlags(INLET);

    std::array<double, 4> potential{1.39572, 110.69275, 221.1549827, 304.284736};

    Matrix LHS_finite_diference = ZeroMatrix(number_of_nodes, number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(number_of_nodes, number_of_nodes);

    PotentialFlowTestUtilities::ComputeElementalSensitivities<4>(
        model_part, LHS_finite_diference, LHS_analytical, potential);

    KRATOS_CHECK_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(PingTransonicPerturbationPotentialFlowSupersonicElementLHS3D,
                          CompressiblePotentialApplicationFastSuite)
{
    const double potential_1 = 1.386736;
    const double potential_2 = 210.6927598;
    const double potential_3 = 221.1549827;
    const double potential_4 = 304.2847368;
    const double potential_5 = -3.349587;
    // mach number squared = 1.35
    std::array<double, 4> high_potential{potential_1, potential_2, potential_3, potential_4};
    // mach number squared = 0.95
    std::array<double, 4> low_potential{potential_5, potential_1, potential_3, potential_4};

    const unsigned int number_of_nodes = 4;
    Matrix LHS_finite_diference = ZeroMatrix(number_of_nodes + 1, number_of_nodes + 1);
    Matrix LHS_analytical = ZeroMatrix(number_of_nodes + 1, number_of_nodes + 1);

    ComputeElementalSensitivitiesTransonic3D(LHS_finite_diference, LHS_analytical, high_potential, low_potential);

    KRATOS_CHECK_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(PingTransonicPerturbationPotentialFlowSupersonicElementLHSClamped3D, CompressiblePotentialApplicationFastSuite) {
    const double potential_1 = 1.386736;
    const double potential_2 = 610.6927598;
    const double potential_3 = 221.1549827;
    const double potential_4 = 304.2847368;
    const double potential_5 = -300.349587;
    // mach number squared = 3.0
    std::array<double, 4> high_potential{potential_1, potential_2, potential_3, potential_4};
    // mach number squared = 3.0
    std::array<double, 4> low_potential{potential_5, potential_1, potential_3, potential_4};

    const unsigned int number_of_nodes = 4;
    Matrix LHS_finite_diference = ZeroMatrix(number_of_nodes + 1, number_of_nodes + 1);
    Matrix LHS_analytical = ZeroMatrix(number_of_nodes + 1, number_of_nodes + 1);

    ComputeElementalSensitivitiesTransonic3D(LHS_finite_diference, LHS_analytical, high_potential, low_potential);

    KRATOS_CHECK_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(PingTransonicPerturbationPotentialFlowSupersonicDeceleratingElementLHS3D, CompressiblePotentialApplicationFastSuite) {
    const double potential_1 = 1.386736;
    const double potential_2 = 210.6927598;
    const double potential_3 = 221.1549827;
    const double potential_4 = 304.2847368;
    const double potential_5 = -103.349587;
    // mach number squared = 1.35
    std::array<double, 4> high_potential{potential_1, potential_2, potential_3, potential_4};
    // mach number squared = 1.79
    std::array<double, 4> low_potential{potential_5, potential_1, potential_3, potential_4};

    const unsigned int number_of_nodes = 4;
    Matrix LHS_finite_diference = ZeroMatrix(number_of_nodes + 1, number_of_nodes + 1);
    Matrix LHS_analytical = ZeroMatrix(number_of_nodes + 1, number_of_nodes + 1);

    ComputeElementalSensitivitiesTransonic3D(LHS_finite_diference, LHS_analytical, high_potential, low_potential);

    KRATOS_CHECK_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(PingTransonicPerturbationPotentialFlowSupersonicDeceleratingElementLHSClamped3D, CompressiblePotentialApplicationFastSuite) {
    const double potential_1 = 1.386736;
    const double potential_2 = 210.6927598;
    const double potential_3 = 221.1549827;
    const double potential_4 = 304.2847368;
    const double potential_5 = -300.349587;
    // mach number squared = 1.35
    std::array<double, 4> high_potential{potential_1, potential_2, potential_3, potential_4};
    // mach number squared = 3.0
    std::array<double, 4> low_potential{potential_5, potential_1, potential_3, potential_4};

    const unsigned int number_of_nodes = 4;
    Matrix LHS_finite_diference = ZeroMatrix(number_of_nodes + 1, number_of_nodes + 1);
    Matrix LHS_analytical = ZeroMatrix(number_of_nodes + 1, number_of_nodes + 1);

    ComputeElementalSensitivitiesTransonic3D(LHS_finite_diference, LHS_analytical, high_potential, low_potential);

    KRATOS_CHECK_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(TransonicPerturbationPotentialFlowElementEquationId3D,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement3D(model_part);
    GenerateTransonicPerturbationUpwindElement3D(model_part);

    FindNodalNeighboursProcess find_nodal_neighbours_process(model_part);
    find_nodal_neighbours_process.Execute();

    Element::Pointer pCurrentElement = model_part.pGetElement(1);
    Element::Pointer pUpwindElement = model_part.pGetElement(2);
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    pCurrentElement->Initialize(r_current_process_info);

    for (auto& r_node : model_part.Nodes()){
        r_node.AddDof(VELOCITY_POTENTIAL);
    }

    Element::DofsVectorType CurrentElementalDofList;
    pCurrentElement->GetDofList(CurrentElementalDofList, r_current_process_info);

    Element::DofsVectorType UpwindElementalDofList;
    pUpwindElement->GetDofList(UpwindElementalDofList, r_current_process_info);

    std::vector<int> current_ids{23, 74, 55, 35}; // 1 2 3 4
    std::vector<int> upwind_ids{87, 23, 55, 35};  // 5 1 3 4
    for (unsigned int i = 0; i < 4; i++) {
        CurrentElementalDofList[i]->SetEquationId(current_ids[i]);
        UpwindElementalDofList[i]->SetEquationId(upwind_ids[i]);
    }

    // make and check equation ids
    Element::EquationIdVectorType EquationIdVector;
    pCurrentElement->EquationIdVector(EquationIdVector, r_current_process_info);

    std::vector<double> reference{23.0, 74.0, 55.0, 35.0, 87.0};
    KRATOS_CHECK_VECTOR_NEAR(EquationIdVector, reference, 1e-15);
}

KRATOS_TEST_CASE_IN_SUITE(PingWakeTransonicPerturbationPotentialFlowElementLHS3D,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement3D(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    std::array<double, 8> potential{1.39572, 110.69275, 121.1549827, 104.284736,
                                    2.39572, 46.69275,  100.1549827, 102.284736};

    Matrix LHS_finite_diference = ZeroMatrix(2 * number_of_nodes, 2 * number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(2 * number_of_nodes, 2 * number_of_nodes);

    PotentialFlowTestUtilities::ComputeWakeElementalSensitivities<4>(
        model_part, LHS_finite_diference, LHS_analytical, potential);

    KRATOS_CHECK_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(PingWakeTransonicPerturbationPotentialFlowElementLHS3DClamping,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement3D(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    std::array<double, 8> potential{1.39572, 357.69275, 321.1549827, 304.284736,
                                    2.39572, 346.69275, 200.1549827, 302.284736};

    Matrix LHS_finite_diference = ZeroMatrix(2 * number_of_nodes, 2 * number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(2 * number_of_nodes, 2 * number_of_nodes);

    PotentialFlowTestUtilities::ComputeWakeElementalSensitivities<4>(
        model_part, LHS_finite_diference, LHS_analytical, potential);

    KRATOS_CHECK_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(PingWakeStructureTransonicPerturbationPotentialFlowElementLHS3D,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement3D(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    p_element->Set(STRUCTURE);
    p_element->GetGeometry()[number_of_nodes-1].SetValue(TRAILING_EDGE, true);

    std::array<double, 8> potential{1.39572, 110.69275, 121.1549827, 104.284736,
                                    2.39572,  46.69275, 100.1549827, 102.284736};

    Matrix LHS_finite_diference = ZeroMatrix(2 * number_of_nodes, 2 * number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(2 * number_of_nodes, 2 * number_of_nodes);

    PotentialFlowTestUtilities::ComputeWakeElementalSensitivities<4>(
        model_part, LHS_finite_diference, LHS_analytical, potential);

    KRATOS_CHECK_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(PingWakeStructureTransonicPerturbationPotentialFlowElementLHS3DClamping,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement3D(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    p_element->Set(STRUCTURE);
    p_element->GetGeometry()[number_of_nodes-1].SetValue(TRAILING_EDGE, true);

    std::array<double, 8> potential{1.39572, 357.69275, 321.1549827, 304.284736,
                                    2.39572, 346.69275, 200.1549827, 302.284736};

    Matrix LHS_finite_diference = ZeroMatrix(2 * number_of_nodes, 2 * number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(2 * number_of_nodes, 2 * number_of_nodes);

    PotentialFlowTestUtilities::ComputeWakeElementalSensitivities<4>(
        model_part, LHS_finite_diference, LHS_analytical, potential);

    KRATOS_CHECK_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(TransonicPerturbationPotentialFlowElementRHS3D,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement3D(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    p_element->SetFlags(INLET);

    std::array<double, 4> potential{1.39572, 110.69275, 121.1549827, 104.284736};
    PotentialFlowTestUtilities::AssignPotentialsToNormalElement<4>(*p_element, potential);

    // Compute RHS
    Vector RHS = ZeroVector(3);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateRightHandSide(RHS, r_current_process_info);

    std::vector<double> reference{71.66991905097665, -64.11826564927853,
                                  -3.932086180475159, -3.619567221222969};

    KRATOS_CHECK_VECTOR_NEAR(RHS, reference, 1e-13);
}

KRATOS_TEST_CASE_IN_SUITE(TransonicPerturbationPotentialFlowInletElementRHS3D,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement3D(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->Initialize(r_current_process_info);

    std::array<double, 4> potential{1.39572, 110.69275, 121.1549827, 104.284736};
    PotentialFlowTestUtilities::AssignPotentialsToNormalElement<4>(*p_element, potential);

    // Compute RHS
    Vector RHS = ZeroVector(3);

    p_element->CalculateRightHandSide(RHS, r_current_process_info);

    std::vector<double> reference{71.66991905097665, -64.11826564927853,
                                  -3.932086180475159, -3.619567221222969};

    KRATOS_CHECK_VECTOR_NEAR(RHS, reference, 1e-13);
}

KRATOS_TEST_CASE_IN_SUITE(TransonicPerturbationPotentialFlowInletElementLHS3D,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement3D(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    p_element->SetFlags(INLET);

    std::array<double, 4> potential{1.39572, 110.69275, 121.1549827, 104.284736};
    PotentialFlowTestUtilities::AssignPotentialsToNormalElement<4>(*p_element,potential);

    // Compute LHS
    Matrix LHS = ZeroMatrix(3, 3);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateLeftHandSide(LHS, r_current_process_info);

    std::vector<double> reference{
        0.1376306838326365,   0.0248823393106424,  -0.06452385209048064,
        -0.09798917105279824, 0.0248823393106424,  0.06159497346557513,
        -0.0664293742338969,  -0.0200479385423207, -0.06452385209048064,
        -0.0664293742338969,  0.1825781102929537,  -0.0516248839685762,
        -0.09798917105279824, -0.0200479385423207, -0.0516248839685762,
        0.1696619935636952};

    for (unsigned int i = 0; i < LHS.size1(); i++) {
        for (unsigned int j = 0; j < LHS.size2(); j++) {
            KRATOS_CHECK_NEAR(LHS(i, j), reference[i * 4 + j], 1e-16);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(TransonicPerturbationPotentialFlowSupersonicElementLHS3D,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement3D(model_part);
    GenerateTransonicPerturbationUpwindElement3D(model_part);

    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();
    Element::Pointer pUpwindElement = model_part.pGetElement(2);

    FindNodalNeighboursProcess find_nodal_neighbours_process(model_part);
    find_nodal_neighbours_process.Execute();

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->Initialize(r_current_process_info);
    pUpwindElement->SetFlags(INLET);

    const double potential_1 = 1.386736;
    const double potential_2 = 210.6927598;
    const double potential_3 = 221.1549827;
    const double potential_4 = 304.2847368;
    const double potential_5 = -3.349587;
    // mach number squared = 1.35
    std::array<double, 4> high_potential{potential_1, potential_2, potential_3, potential_4};
    PotentialFlowTestUtilities::AssignPotentialsToNormalElement<4>(*p_element,high_potential);
    // mach number squared = 0.95
    std::array<double, 4> low_potential{potential_5, potential_1, potential_3, potential_4};
    PotentialFlowTestUtilities::AssignPotentialsToNormalElement<4>(*pUpwindElement,low_potential);

    for (auto& r_node : model_part.Nodes()){
        r_node.AddDof(VELOCITY_POTENTIAL);
    }

    Element::DofsVectorType CurrentElementalDofList;
    p_element->GetDofList(CurrentElementalDofList, r_current_process_info);

    Element::DofsVectorType UpwindElementalDofList;
    pUpwindElement->GetDofList(UpwindElementalDofList, r_current_process_info);

    std::vector<int> current_ids{23, 74, 55, 35}; // 1 2 3 4
    std::vector<int> upwind_ids{87, 23, 55, 35};  // 5 1 3 4
    for (unsigned int i = 0; i < number_of_nodes; i++) {
        CurrentElementalDofList[i]->SetEquationId(current_ids[i]);
        UpwindElementalDofList[i]->SetEquationId(upwind_ids[i]);
    }

    // Compute LHS
    Matrix LHS = ZeroMatrix(3, 3);

    p_element->CalculateLeftHandSide(LHS, r_current_process_info);

    std::vector<double> reference{
        0.07266558291387809,-0.02095045152455358,0.0206984047580022,0.07804597981741059,-0.1504595159647373,
        -0.001147820720733481,0.09455275283978169,-0.08611110193404595,-0.1091917784084459,0.1018979482234437,
        -0.04216803465169104,-0.04244906075669839,0.1188411475887801,-0.04247590493843136,0.008251852758040613,
        -0.02934972754145357,-0.03115324055852974,-0.05342845041273643,0.07362170352946668,0.04030971498325305,
        0,-0,-0,-0,0};

    for (unsigned int i = 0; i < LHS.size1(); i++) {
        for (unsigned int j = 0; j < LHS.size2(); j++) {
            KRATOS_CHECK_NEAR(LHS(i, j), reference[i * 5 + j], 1e-16);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(TransonicPerturbationPotentialFlowSupersonicElementLHS3DClamping,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement3D(model_part);
    GenerateTransonicPerturbationUpwindElement3D(model_part);

    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();
    Element::Pointer pUpwindElement = model_part.pGetElement(2);

    FindNodalNeighboursProcess find_nodal_neighbours_process(model_part);
    find_nodal_neighbours_process.Execute();

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->Initialize(r_current_process_info);
    pUpwindElement->SetFlags(INLET);

    const double potential_1 = 1.386736;
    const double potential_2 = 610.6927598;
    const double potential_3 = 221.1549827;
    const double potential_4 = 304.2847368;
    const double potential_5 = -300.349587;
    // mach number squared = 3.0
    std::array<double, 4> high_potential{potential_1, potential_2, potential_3, potential_4};
    PotentialFlowTestUtilities::AssignPotentialsToNormalElement<4>(*p_element,high_potential);
    // mach number squared = 3.0
    std::array<double, 4> low_potential{potential_5, potential_1, potential_3, potential_4};
    PotentialFlowTestUtilities::AssignPotentialsToNormalElement<4>(*pUpwindElement,low_potential);

    for (auto& r_node : model_part.Nodes()){
        r_node.AddDof(VELOCITY_POTENTIAL);
    }

    Element::DofsVectorType CurrentElementalDofList;
    p_element->GetDofList(CurrentElementalDofList, r_current_process_info);

    Element::DofsVectorType UpwindElementalDofList;
    pUpwindElement->GetDofList(UpwindElementalDofList, r_current_process_info);

    std::vector<int> current_ids{23, 74, 55, 35}; // 1 2 3 4
    std::vector<int> upwind_ids{87, 23, 55, 35};  // 5 1 3 4
    for (unsigned int i = 0; i < number_of_nodes; i++) {
        CurrentElementalDofList[i]->SetEquationId(current_ids[i]);
        UpwindElementalDofList[i]->SetEquationId(upwind_ids[i]);
    }

    // Compute LHS
    Matrix LHS = ZeroMatrix(3, 3);

    p_element->CalculateLeftHandSide(LHS, r_current_process_info);

    std::vector<double> reference{
        0.1618269140113864,-0.07651948119478803,-0.03515180648724157,-0.05015562632935684,-0,
        -0.07651948119478803,0.1071272736727032,-0.02580657012843831,-0.004801222349476892,0,
        -0.03515180648724157,-0.02580657012843831,0.084492939739455,-0.02353456312377514,-0,
        -0.05015562632935684,-0.004801222349476892,-0.02353456312377514,0.07849141180260888,0,
        0,0,0,0,0};

    for (unsigned int i = 0; i < LHS.size1(); i++) {
        for (unsigned int j = 0; j < LHS.size2(); j++) {
            KRATOS_CHECK_NEAR(LHS(i, j), reference[i * 5 + j], 1e-16);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(TransonicPerturbationPotentialFlowSupersonicDeceleratingElementLHS3D,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement3D(model_part);
    GenerateTransonicPerturbationUpwindElement3D(model_part);

    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();
    Element::Pointer pUpwindElement = model_part.pGetElement(2);

    FindNodalNeighboursProcess find_nodal_neighbours_process(model_part);
    find_nodal_neighbours_process.Execute();

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->Initialize(r_current_process_info);
    pUpwindElement->SetFlags(INLET);

    const double potential_1 = 1.386736;
    const double potential_2 = 210.6927598;
    const double potential_3 = 221.1549827;
    const double potential_4 = 304.2847368;
    const double potential_5 = -103.349587;
    // mach number squared = 1.35
    std::array<double, 4> high_potential{potential_1, potential_2, potential_3, potential_4};
    PotentialFlowTestUtilities::AssignPotentialsToNormalElement<4>(*p_element,high_potential);
    // mach number squared = 1.79
    std::array<double, 4> low_potential{potential_5, potential_1, potential_3, potential_4};
    PotentialFlowTestUtilities::AssignPotentialsToNormalElement<4>(*pUpwindElement,low_potential);

    for (auto& r_node : model_part.Nodes()){
        r_node.AddDof(VELOCITY_POTENTIAL);
    }

    Element::DofsVectorType CurrentElementalDofList;
    p_element->GetDofList(CurrentElementalDofList, r_current_process_info);

    Element::DofsVectorType UpwindElementalDofList;
    pUpwindElement->GetDofList(UpwindElementalDofList, r_current_process_info);

    std::vector<int> current_ids{23, 74, 55, 35}; // 1 2 3 4
    std::vector<int> upwind_ids{87, 23, 55, 35};  // 5 1 3 4
    for (unsigned int i = 0; i < number_of_nodes; i++) {
        CurrentElementalDofList[i]->SetEquationId(current_ids[i]);
        UpwindElementalDofList[i]->SetEquationId(upwind_ids[i]);
    }

    // Compute LHS
    Matrix LHS = ZeroMatrix(3, 3);

    p_element->CalculateLeftHandSide(LHS, r_current_process_info);

    std::vector<double> reference{
        0.04648677446666163,0.01195966508506741,0.03759047222406028,0.115287595888341,-0.2113245076641303,
        0.006233937351858579,0.05496313396997898,-0.08203039475694361,-0.1222851333951752,0.1431184568302813,
        -0.03251200493727034,-0.03486786947275742,0.09208529523935334,-0.03629537384910351,0.0115899530197779,
        -0.02020870688124987,-0.03205492958228898,-0.04764537270647004,0.04329291135593773,0.05661609781407116,
        0,-0,-0,-0,0};

    for (unsigned int i = 0; i < LHS.size1(); i++) {
        for (unsigned int j = 0; j < LHS.size2(); j++) {
            KRATOS_CHECK_NEAR(LHS(i, j), reference[i * 5 + j], 1e-16);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(TransonicPerturbationPotentialFlowSupersonicDeceleratingElementLHS3DClamping,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement3D(model_part);
    GenerateTransonicPerturbationUpwindElement3D(model_part);

    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();
    Element::Pointer pUpwindElement = model_part.pGetElement(2);

    FindNodalNeighboursProcess find_nodal_neighbours_process(model_part);
    find_nodal_neighbours_process.Execute();

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->Initialize(r_current_process_info);
    pUpwindElement->SetFlags(INLET);

    const double potential_1 = 1.386736;
    const double potential_2 = 210.6927598;
    const double potential_3 = 221.1549827;
    const double potential_4 = 304.2847368;
    const double potential_5 = -300.349587;
    // mach number squared = 1.35
    std::array<double, 4> high_potential{potential_1, potential_2, potential_3, potential_4};
    PotentialFlowTestUtilities::AssignPotentialsToNormalElement<4>(*p_element,high_potential);
    // mach number squared = 3.0
    std::array<double, 4> low_potential{potential_5, potential_1, potential_3, potential_4};
    PotentialFlowTestUtilities::AssignPotentialsToNormalElement<4>(*pUpwindElement,low_potential);

    for (auto& r_node : model_part.Nodes()){
        r_node.AddDof(VELOCITY_POTENTIAL);
    }

    Element::DofsVectorType CurrentElementalDofList;
    p_element->GetDofList(CurrentElementalDofList, r_current_process_info);

    Element::DofsVectorType UpwindElementalDofList;
    pUpwindElement->GetDofList(UpwindElementalDofList, r_current_process_info);

    std::vector<int> current_ids{23, 74, 55, 35}; // 1 2 3 4
    std::vector<int> upwind_ids{87, 23, 55, 35};  // 5 1 3 4
    for (unsigned int i = 0; i < number_of_nodes; i++) {
        CurrentElementalDofList[i]->SetEquationId(current_ids[i]);
        UpwindElementalDofList[i]->SetEquationId(upwind_ids[i]);
    }

    // Compute LHS
    Matrix LHS = ZeroMatrix(3, 3);

    p_element->CalculateLeftHandSide(LHS, r_current_process_info);

    std::vector<double> reference{
        0.0396349216873926,0.0094190248662884,-0.03098012588944438,-0.01807382066423663,-0,
        0.0094190248662884,0.05425051039991628,-0.03340847095809917,-0.03026106430810551,0,
        -0.03098012588944438,-0.03340847095809917,0.09221337106454636,-0.02782477421700282,0,
        -0.01807382066423663,-0.03026106430810551,-0.02782477421700282,0.07615965918934496,0,
        0,-0,-0,-0,0};

    for (unsigned int i = 0; i < LHS.size1(); i++) {
        for (unsigned int j = 0; j < LHS.size2(); j++) {
            KRATOS_CHECK_NEAR(LHS(i, j), reference[i * 5 + j], 1e-16);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(TransonicPerturbationPotentialFlowSupersonicElementRHS3D,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement3D(model_part);
    GenerateTransonicPerturbationUpwindElement3D(model_part);

    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();
    Element::Pointer pUpwindElement = model_part.pGetElement(2);

    FindNodalNeighboursProcess find_nodal_neighbours_process(model_part);
    find_nodal_neighbours_process.Execute();

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->Initialize(r_current_process_info);
    pUpwindElement->SetFlags(INLET);

    const double potential_1 = 1.386736;
    const double potential_2 = 210.6927598;
    const double potential_3 = 221.1549827;
    const double potential_4 = 304.2847368;
    const double potential_5 = -3.349587;
    // mach number squared = 1.35
    std::array<double, 4> high_potential{potential_1, potential_2, potential_3, potential_4};
    PotentialFlowTestUtilities::AssignPotentialsToNormalElement<4>(*p_element,high_potential);
    // mach number squared = 0.95
    std::array<double, 4> low_potential{potential_5, potential_1, potential_3, potential_4};
    PotentialFlowTestUtilities::AssignPotentialsToNormalElement<4>(*pUpwindElement,low_potential);

    for (auto& r_node : model_part.Nodes()){
        r_node.AddDof(VELOCITY_POTENTIAL);
    }

    Element::DofsVectorType CurrentElementalDofList;
    p_element->GetDofList(CurrentElementalDofList, r_current_process_info);

    Element::DofsVectorType UpwindElementalDofList;
    pUpwindElement->GetDofList(UpwindElementalDofList, r_current_process_info);

    std::vector<int> current_ids{23, 74, 55, 35}; // 1 2 3 4
    std::vector<int> upwind_ids{87, 23, 55, 35};  // 5 1 3 4
    for (unsigned int i = 0; i < number_of_nodes; i++) {
        CurrentElementalDofList[i]->SetEquationId(current_ids[i]);
        UpwindElementalDofList[i]->SetEquationId(upwind_ids[i]);
    }

    // Compute RHS
    Vector RHS = ZeroVector(4);

    p_element->CalculateRightHandSide(RHS, r_current_process_info);

    std::vector<double> reference{78.83234736892321,-53.38880960120657,-4.323508014019686,-21.12002975369695,0};

    KRATOS_CHECK_VECTOR_NEAR(RHS, reference, 1e-14);
}

KRATOS_TEST_CASE_IN_SUITE(WakeTransonicPerturbationPotentialFlowElementRHS3D,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement3D(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    BoundedVector<double,4> distances = PotentialFlowTestUtilities::AssignDistancesToElement<4>();

    p_element->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    p_element->GetValue(WAKE) = true;

    std::array<double, 8> potential{1.39572, 110.69275, 121.1549827, 104.284736,
                                    2.39572,  46.69275, 100.1549827, 102.284736};
    PotentialFlowTestUtilities::AssignPotentialsToWakeElement<4>(*p_element,distances,potential);

    // Compute RHS and LHS
    Vector RHS = ZeroVector(6);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateRightHandSide(RHS, r_current_process_info);

    std::vector<double> reference{11.25952380952381, -14.46333333333333,
                                  2.251904761904762, -3.619567221222969,
                                  68.65551596318301, -58.62766030853704,
                                  -4.30462713896052, -0.9519047619047626};

    KRATOS_CHECK_VECTOR_NEAR(RHS, reference, 1e-13);
}

KRATOS_TEST_CASE_IN_SUITE(WakeStructureTransonicPerturbationPotentialFlowElementRHS3D,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement3D(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    BoundedVector<double,4> distances = PotentialFlowTestUtilities::AssignDistancesToElement<4>();

    p_element->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    p_element->GetValue(WAKE) = true;
    p_element->Set(STRUCTURE);
    p_element->GetGeometry()[number_of_nodes-1].SetValue(TRAILING_EDGE, true);

    std::array<double, 8> potential{1.39572, 110.69275, 121.1549827, 104.284736,
                                    2.39572,  46.69275, 100.1549827, 102.284736};
    PotentialFlowTestUtilities::AssignPotentialsToWakeElement<4>(*p_element,distances,potential);

    // Compute RHS and LHS
    Vector RHS = ZeroVector(6);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateRightHandSide(RHS, r_current_process_info);

    std::vector<double> reference{11.25952380952381, -14.46333333333333,
                                  2.251904761904762, -0.4524459026528712,
                                  68.65551596318301, -58.62766030853704,
                                  -4.30462713896052, -5.007824951224748};

    KRATOS_CHECK_VECTOR_NEAR(RHS, reference, 1e-13);
}

KRATOS_TEST_CASE_IN_SUITE(WakeTransonicPerturbationPotentialFlowElementLHS3D,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement3D(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    BoundedVector<double,4> distances = PotentialFlowTestUtilities::AssignDistancesToElement<4>();

    p_element->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    p_element->GetValue(WAKE) = true;

    std::array<double, 8> potential{1.39572, 110.69275, 121.1549827, 104.284736,
                                    2.39572,  46.69275, 100.1549827, 102.284736};
    PotentialFlowTestUtilities::AssignPotentialsToWakeElement<4>(*p_element,distances,potential);

    // Compute LHS
    Matrix LHS = ZeroMatrix(6, 6);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateLeftHandSide(LHS, r_current_process_info);

    // Check the LHS values
    std::vector<double> reference{0.263095238095238,-0.185,0.05261904761904762,-0.1307142857142857,
    -0.263095238095238,0.185,-0.05261904761904762,0.1307142857142857,-0.185,0.2356666666666666,
    -0.03699999999999999,-0.01366666666666666,0.185,-0.2356666666666666,0.03699999999999999,
    0.01366666666666666,0.05261904761904762,-0.03699999999999999,0.01052380952380952,
    -0.02614285714285715,-0.05261904761904762,0.03699999999999999,-0.01052380952380952,
    0.02614285714285715,-0.09798917105279824,-0.0200479385423207,-0.0516248839685762,
    0.1696619935636952,0,0,0,0,0,0,0,0,0.2513125629816388,-0.05867500741543832,-0.07898156406683947,
    -0.1136559914993611,0,0,0,0,-0.05867500741543832,0.1557535754651047,-0.07370192884392225,
    -0.0233766392057441,0,0,0,0,-0.07898156406683947,-0.07370192884392225,0.2130140967141911,
    -0.06033060380342937,0.1307142857142857,0.01366666666666666,0.02614285714285715,
    -0.1705238095238095,-0.1307142857142857,-0.01366666666666666,-0.02614285714285715,
    0.1705238095238095};

    for (unsigned int i = 0; i < LHS.size1(); i++)
    {
        for (unsigned int j = 0; j < LHS.size2(); j++)
        {
            KRATOS_CHECK_NEAR(LHS(i, j), reference[8 * i + j], 1e-16);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(WakeTransonicPerturbationPotentialFlowElementLHS3DClamping,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement3D(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    BoundedVector<double,4> distances = PotentialFlowTestUtilities::AssignDistancesToElement<4>();

    p_element->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    p_element->GetValue(WAKE) = true;

    std::array<double, 8> potential{1.39572, 357.69275, 321.1549827, 304.284736,
                                    2.39572, 346.69275, 200.1549827, 302.284736};
    PotentialFlowTestUtilities::AssignPotentialsToWakeElement<4>(*p_element,distances,potential);

    // Compute LHS
    Matrix LHS = ZeroMatrix(6, 6);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateLeftHandSide(LHS, r_current_process_info);

    // Check the LHS values
    std::vector<double> reference{0.263095238095238,-0.185,0.05261904761904762,
    -0.1307142857142857,-0.263095238095238,0.185,-0.05261904761904762,
    0.1307142857142857,-0.185,0.2356666666666666,-0.03699999999999999,
    -0.01366666666666666,0.185,-0.2356666666666666,0.03699999999999999,
    0.01366666666666666,0.05261904761904762,-0.03699999999999999,
    0.01052380952380952,-0.02614285714285715,-0.05261904761904762,
    0.03699999999999999,-0.01052380952380952,0.02614285714285715,
    -0.05015562632935684,-0.004801222349476892,-0.02353456312377514,
    0.07849141180260888,0,0,0,0,0,0,0,0,0.1618269140113864,-0.07651948119478803,
    -0.03515180648724157,-0.05015562632935684,0,0,0,0,-0.07651948119478803,
    0.1071272736727032,-0.02580657012843831,-0.004801222349476892,0,0,0,0,
    -0.03515180648724157,-0.02580657012843831,0.084492939739455,
    -0.02353456312377514,0.1307142857142857,0.01366666666666666,
    0.02614285714285715,-0.1705238095238095,-0.1307142857142857,
    -0.01366666666666666,-0.02614285714285715,0.1705238095238095};

    for (unsigned int i = 0; i < LHS.size1(); i++)
    {
        for (unsigned int j = 0; j < LHS.size2(); j++)
        {
            KRATOS_CHECK_NEAR(LHS(i, j), reference[8 * i + j], 1e-16);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(WakeStructureTransonicPerturbationPotentialFlowElementLHS3D,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement3D(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    BoundedVector<double,4> distances = PotentialFlowTestUtilities::AssignDistancesToElement<4>();

    p_element->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    p_element->GetValue(WAKE) = true;
    p_element->Set(STRUCTURE);
    p_element->GetGeometry()[number_of_nodes-1].SetValue(TRAILING_EDGE, true);

    std::array<double, 8> potential{1.39572, 110.69275, 121.1549827, 104.284736,
                                    2.39572,  46.69275, 100.1549827, 102.284736};
    PotentialFlowTestUtilities::AssignPotentialsToWakeElement<4>(*p_element,distances,potential);

    // Compute LHS
    Matrix LHS = ZeroMatrix(6, 6);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateLeftHandSide(LHS, r_current_process_info);

    // Check the LHS values
    std::vector<double> reference{0.263095238095238,-0.185,0.05261904761904762,-0.1307142857142857,
    -0.263095238095238,0.185,-0.05261904761904762,0.1307142857142857,-0.185,0.2356666666666666,
    -0.03699999999999999,-0.01366666666666666,0.185,-0.2356666666666666,0.03699999999999999,
    0.01366666666666666,0.05261904761904762,-0.03699999999999999,0.01052380952380952,
    -0.02614285714285715,-0.05261904761904762,0.03699999999999999,-0.01052380952380952,
    0.02614285714285715,-0.01224864638159978,-0.002505992317790087,-0.006453110496072026,
    0.02120774919546189,0,0,0,0,0,0,0,0,0.2513125629816388,-0.05867500741543832,
    -0.07898156406683947,-0.1136559914993611,0,0,0,0,-0.05867500741543832,0.1557535754651047,
    -0.07370192884392225,-0.0233766392057441,0,0,0,0,-0.07898156406683947,-0.07370192884392225,
    0.2130140967141911,-0.06033060380342937,0,0,0,0,-0.09944899256194101,-0.0204545593050261,
    -0.05278927832800071,0.1726928301949679};

    for (unsigned int i = 0; i < LHS.size1(); i++)
    {
        for (unsigned int j = 0; j < LHS.size2(); j++)
        {
            KRATOS_CHECK_NEAR(LHS(i, j), reference[8 * i + j], 1e-16);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(WakeStructureTransonicPerturbationPotentialFlowElementLHS3DClamping,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement3D(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    BoundedVector<double,4> distances = PotentialFlowTestUtilities::AssignDistancesToElement<4>();

    p_element->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    p_element->GetValue(WAKE) = true;
    p_element->Set(STRUCTURE);
    p_element->GetGeometry()[number_of_nodes-1].SetValue(TRAILING_EDGE, true);

    std::array<double, 8> potential{1.39572, 357.69275, 321.1549827, 304.284736,
                                    2.39572, 346.69275, 200.1549827, 302.284736};
    PotentialFlowTestUtilities::AssignPotentialsToWakeElement<4>(*p_element,distances,potential);

    // Compute LHS
    Matrix LHS = ZeroMatrix(6, 6);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateLeftHandSide(LHS, r_current_process_info);

    // Check the LHS values
    std::vector<double> reference{0.263095238095238,-0.185,0.05261904761904762,
    -0.1307142857142857,-0.263095238095238,0.185,-0.05261904761904762,
    0.1307142857142857,-0.185,0.2356666666666666,-0.03699999999999999,
    -0.01366666666666666,0.185,-0.2356666666666666,0.03699999999999999,
    0.01366666666666666,0.05261904761904762,-0.03699999999999999,
    0.01052380952380952,-0.02614285714285715,-0.05261904761904762,
    0.03699999999999999,-0.01052380952380952,0.02614285714285715,
    -0.006269453291169606,-0.0006001527936846115,-0.002941820390471893,
    0.009811426475326112,0,0,0,0,0,0,0,0,0.1618269140113864,
    -0.07651948119478803,-0.03515180648724157,-0.05015562632935684,0,0,0,0,
    -0.07651948119478803,0.1071272736727032,-0.02580657012843831,
    -0.004801222349476892,0,0,0,0,-0.03515180648724157,-0.02580657012843831,
    0.084492939739455,-0.02353456312377514,0,0,0,0,-0.04388617303818725,
    -0.004201069555792282,-0.02059274273330326,0.0686799853272828};

    for (unsigned int i = 0; i < LHS.size1(); i++)
    {
        for (unsigned int j = 0; j < LHS.size2(); j++)
        {
            KRATOS_CHECK_NEAR(LHS(i, j), reference[8 * i + j], 1e-16);
        }
    }
}

} // namespace Testing
} // namespace Kratos.
