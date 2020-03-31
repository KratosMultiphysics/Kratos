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
#include "processes/find_nodal_neighbours_process.h"
#include "custom_utilities/potential_flow_utilities.h"

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
    rModelPart.GetProcessInfo()[MACH_LIMIT] = 0.99;
    rModelPart.GetProcessInfo()[UPWINDING_FACTOR_CONSTANT] = 1.0;

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

void GenerateTestingTransonicUpstreamElement(ModelPart& rModelPart) {
    // Create extra node
    rModelPart.CreateNewNode(4, 0.0, 1.0, 0.0);
    // Nodes Ids
    std::vector<ModelPart::IndexType> upstream_elemNodes{1, 3, 4};
    Properties::Pointer pElemProp = rModelPart.pGetProperties(0);
    rModelPart.CreateNewElement("TransonicPerturbationPotentialFlowElement2D3N", 2, upstream_elemNodes, pElemProp);
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
    Element::Pointer pElement = model_part.pGetElement(1);
    AssignPotentialsToNormalTransonicPerturbationElement(pElement);

    GenerateTestingTransonicUpstreamElement(model_part);
    Element::Pointer pUpstreamElement = model_part.pGetElement(2);
    std::array<double, 3> upstream_potential{1.0, 150.0, 90.0};
    AssignPerturbationPotentialsToTransonicElement(*pUpstreamElement, upstream_potential);

    FindNodalNeighboursProcess find_nodal_neighbours_process(model_part);
    find_nodal_neighbours_process.Execute();
    pElement->Initialize(model_part.GetProcessInfo());

    // Compute RHS
    Vector RHS = ZeroVector(3);

    pElement->CalculateRightHandSide(RHS, model_part.GetProcessInfo());

    std::vector<double> reference{146.2643261263345,-122.1426284341492,-24.12169769218525, 0.0};

    KRATOS_CHECK_VECTOR_NEAR(RHS, reference, 1e-13);
}

/** Checks the TransonicPerturbationPotentialFlowElement.
 * Checks the RHS computation.
 */
KRATOS_TEST_CASE_IN_SUITE(TransonicPerturbationPotentialFlowElementRHSSupersonicAccelerating, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);
    std::array<double, 3> potential{1.0, 120.0, 180.0};
    AssignPerturbationPotentialsToTransonicElement(*pElement, potential);

    GenerateTestingTransonicUpstreamElement(model_part);
    Element::Pointer pUpstreamElement = model_part.pGetElement(2);
    std::array<double, 3> upstream_potential{1.0, 180.0, 90.0};
    AssignPerturbationPotentialsToTransonicElement(*pUpstreamElement, upstream_potential);

    FindNodalNeighboursProcess find_nodal_neighbours_process(model_part);
    find_nodal_neighbours_process.Execute();
    pElement->Initialize(model_part.GetProcessInfo());

    // Compute RHS
    Vector RHS = ZeroVector(3);

    pElement->CalculateRightHandSide(RHS, model_part.GetProcessInfo());

    std::vector<double> reference{146.7051200733924,-119.4685732437509,-27.23654682964152, 0.0};

    KRATOS_CHECK_VECTOR_NEAR(RHS, reference, 1e-13);
}

/** Checks the TransonicPerturbationPotentialFlowElement.
 * Checks the RHS computation.
 */
KRATOS_TEST_CASE_IN_SUITE(TransonicPerturbationPotentialFlowElementRHSSupersonicDecelerating, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);
    std::array<double, 3> potential{1.0, 120.0, 180.0};
    AssignPerturbationPotentialsToTransonicElement(*pElement, potential);

    GenerateTestingTransonicUpstreamElement(model_part);
    Element::Pointer pUpstreamElement = model_part.pGetElement(2);
    std::array<double, 3> upstream_potential{1.0, 180.0, 51.0};
    AssignPerturbationPotentialsToTransonicElement(*pUpstreamElement, upstream_potential);

    FindNodalNeighboursProcess find_nodal_neighbours_process(model_part);
    find_nodal_neighbours_process.Execute();
    pElement->Initialize(model_part.GetProcessInfo());

    // Compute RHS
    Vector RHS = ZeroVector(3);

    pElement->CalculateRightHandSide(RHS, model_part.GetProcessInfo());

    std::vector<double> reference{145.5366000056261,-118.5169948309941,-27.01960517463199, 0.0};

    KRATOS_CHECK_VECTOR_NEAR(RHS, reference, 1e-13);
}

/** Checks the TransonicPerturbationPotentialFlowElement.
 * Checks the RHS computation.
 */
KRATOS_TEST_CASE_IN_SUITE(TransonicPerturbationPotentialFlowElementInlet, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);
    AssignPotentialsToNormalTransonicPerturbationElement(pElement);
    pElement->Initialize(model_part.GetProcessInfo());

    // Compute RHS
    Vector RHS = ZeroVector(3);

    pElement->CalculateRightHandSide(RHS, model_part.GetProcessInfo());

    std::vector<double> reference{146.2643261263345,-122.1426284341492,-24.12169769218525};

    KRATOS_CHECK_VECTOR_NEAR(RHS, reference, 1e-13);
}

/** Checks the TransonicPerturbationPotentialFlowElement.
 * Checks the LHS computation.
 */
KRATOS_TEST_CASE_IN_SUITE(TransonicPerturbationPotentialFlowElementLHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);
    AssignPotentialsToNormalTransonicPerturbationElement(pElement);

    GenerateTestingTransonicUpstreamElement(model_part);
    Element::Pointer pUpstreamElement = model_part.pGetElement(2);
    std::array<double, 3> upstream_potential{1.0, 150.0, 90.0};
    AssignPerturbationPotentialsToTransonicElement(*pUpstreamElement, upstream_potential);

    FindNodalNeighboursProcess find_nodal_neighbours_process(model_part);
    find_nodal_neighbours_process.Execute();
    pElement->Initialize(model_part.GetProcessInfo());

    // Compute LHS
    Matrix LHS = ZeroMatrix(3, 3);

    pElement->CalculateLeftHandSide(LHS, model_part.GetProcessInfo());

    std::array<double, 16> reference{ 0.06114278464441542,-0.1306215050744058, 0.06947872042999037, 0.0,
                                     -0.1306215050744058, 0.6710758508914103,-0.5404543458170046, 0.0,
                                      0.06947872042999037,-0.5404543458170046,0.4709756253870142, 0.0,
                                      0.0, 0.0, 0.0, 0.0};

    for (unsigned int i = 0; i < LHS.size1(); i++) {
        for (unsigned int j = 0; j < LHS.size2(); j++) {
            KRATOS_CHECK_NEAR(LHS(i, j), reference[i * 4 + j], 1e-16);
        }
    }
}

/** Checks the TransonicPerturbationPotentialFlowElement.
 * Checks the LHS computation.
 */
KRATOS_TEST_CASE_IN_SUITE(TransonicPerturbationPotentialFlowElementLHSSuperSonicAccelerating, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);
    std::array<double, 3> potential{1.0, 120.0, 180.0};
    AssignPerturbationPotentialsToTransonicElement(*pElement, potential);

    GenerateTestingTransonicUpstreamElement(model_part);
    Element::Pointer pUpstreamElement = model_part.pGetElement(2);
    std::array<double, 3> upstream_potential{1.0, 180.0, 90.0};
    AssignPerturbationPotentialsToTransonicElement(*pUpstreamElement, upstream_potential);

    for (auto& r_node : model_part.Nodes()){
        r_node.AddDof(VELOCITY_POTENTIAL);
    }

    Element::DofsVectorType ElementalDofList, UpstreamElementalDofList;
    pElement->GetDofList(ElementalDofList, model_part.GetProcessInfo());
    pUpstreamElement->GetDofList(UpstreamElementalDofList, model_part.GetProcessInfo());

    std::vector<int> ids{23, 74, 55};
    std::vector<int> upstream_ids{23, 55, 67};
    for (int i = 0; i < 3; i++){
        ElementalDofList[i]->SetEquationId(ids[i]);
        UpstreamElementalDofList[i]->SetEquationId(upstream_ids[i]);
    }

    FindNodalNeighboursProcess find_nodal_neighbours_process(model_part);
    find_nodal_neighbours_process.Execute();
    pElement->Initialize(model_part.GetProcessInfo());

    // Compute LHS
    Matrix LHS = ZeroMatrix(3, 3);

    pElement->CalculateLeftHandSide(LHS, model_part.GetProcessInfo());

    std::array<double, 16> reference{ 0.08311820046089136,-0.1594672105670851, 0.09759158478435792,-0.02124257467816418,
                                     -0.1519635956028805,  0.6680804733587455,-0.5334156611081742,  0.01729878335230908,
                                      0.06884539514198909,-0.5086132627916604,0.4358240763238163,0.003943791325855098,
                                      0.0, 0.0, 0.0, 0.0};

    for (unsigned int i = 0; i < LHS.size1(); i++) {
        for (unsigned int j = 0; j < LHS.size2(); j++) {
            KRATOS_CHECK_NEAR(LHS(i, j), reference[i * 4 + j], 1e-15);
        }
    }
}

/** Checks the TransonicPerturbationPotentialFlowElement.
 * Tests the LHS computation.
 */
KRATOS_TEST_CASE_IN_SUITE(PingTransonicPerturbationPotentialFlowElementLHSSuperSonicAccelerating, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

GenerateTransonicPerturbationElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);
    const unsigned int number_of_nodes = pElement->GetGeometry().size();
    std::array<double, 3> potential{1.0, 120.0, 180.0};
    AssignPerturbationPotentialsToTransonicElement(*pElement, potential);

    GenerateTestingTransonicUpstreamElement(model_part);
    Element::Pointer pUpstreamElement = model_part.pGetElement(2);
    std::array<double, 3> upstream_potential{1.0, 180.0, 90.0};
    AssignPerturbationPotentialsToTransonicElement(*pUpstreamElement, upstream_potential);

    for (auto& r_node : model_part.Nodes()){
        r_node.AddDof(VELOCITY_POTENTIAL);
    }

    Element::DofsVectorType ElementalDofList, UpstreamElementalDofList;
    pElement->GetDofList(ElementalDofList, model_part.GetProcessInfo());
    pUpstreamElement->GetDofList(UpstreamElementalDofList, model_part.GetProcessInfo());

    std::vector<int> ids{23, 74, 55};
    std::vector<int> upstream_ids{23, 55, 67};
    for (int i = 0; i < 3; i++){
        ElementalDofList[i]->SetEquationId(ids[i]);
        UpstreamElementalDofList[i]->SetEquationId(upstream_ids[i]);
    }

    FindNodalNeighboursProcess find_nodal_neighbours_process(model_part);
    find_nodal_neighbours_process.Execute();
    pElement->Initialize(model_part.GetProcessInfo());

    Vector RHS_original = ZeroVector(number_of_nodes);
    Matrix LHS_original = ZeroMatrix(number_of_nodes, number_of_nodes);
    Matrix LHS_finite_diference = ZeroMatrix(number_of_nodes + 1, number_of_nodes + 1);
    Matrix LHS_analytical = ZeroMatrix(number_of_nodes + 1, number_of_nodes + 1);

    // Compute original RHS and LHS
    pElement->CalculateLocalSystem(LHS_original, RHS_original, model_part.GetProcessInfo());
    // KRATOS_WATCH(LHS_original)
    // KRATOS_WATCH(RHS_original)

    double delta = 1e-3;
    for(unsigned int i = 0; i < number_of_nodes; i++){
        // Pinging
        pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += delta;

        Vector RHS_pinged = ZeroVector(number_of_nodes + 1);
        Matrix LHS_pinged = ZeroMatrix(number_of_nodes + 1, number_of_nodes + 1);
        // Compute pinged LHS and RHS
        pElement->CalculateLocalSystem(LHS_pinged, RHS_pinged, model_part.GetProcessInfo());

        for(unsigned int k = 0; k < number_of_nodes + 1; k++){
            // Compute the finite difference estimate of the sensitivity
            LHS_finite_diference( k, i) = -(RHS_pinged(k)-RHS_original(k)) / delta;
            // Compute the average of the original and pinged analytic sensitivities
            LHS_analytical( k, i) = 0.5 * (LHS_original(k,i) + LHS_pinged(k,i));
        }

        // Unpinging
        pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= delta;
    }

    std::cout.precision(16);
    // KRATOS_WATCH(LHS_analytical)
    // KRATOS_WATCH(LHS_finite_diference)
    for (unsigned int i = 0; i < LHS_finite_diference.size1(); i++) {
        for (unsigned int j = 0; j < LHS_finite_diference.size2(); j++) {
            KRATOS_CHECK_NEAR(LHS_finite_diference(i,j), LHS_analytical(i,j), 1e-10);
        }
    }
}

/** Checks the TransonicPerturbationPotentialFlowElement.
 * Tests the LHS computation.
 */
KRATOS_TEST_CASE_IN_SUITE(PingTransonicPerturbationPotentialFlowElementLHSSuperSonicAccelerating2, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    double maximum_potential = 280.0;
    double minimum_potential = 0.0;

    GenerateTransonicPerturbationElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);
    const unsigned int number_of_nodes = pElement->GetGeometry().size();
    std::array<double, 3> potential{minimum_potential, maximum_potential, maximum_potential};
    AssignPerturbationPotentialsToTransonicElement(*pElement, potential);

    GenerateTestingTransonicUpstreamElement(model_part);
    Element::Pointer pUpstreamElement = model_part.pGetElement(2);
    std::array<double, 3> upstream_potential{minimum_potential, maximum_potential, minimum_potential};
    AssignPerturbationPotentialsToTransonicElement(*pUpstreamElement, upstream_potential);

    for (auto& r_node : model_part.Nodes()){
        r_node.AddDof(VELOCITY_POTENTIAL);
    }

    Element::DofsVectorType ElementalDofList, UpstreamElementalDofList;
    pElement->GetDofList(ElementalDofList, model_part.GetProcessInfo());
    pUpstreamElement->GetDofList(UpstreamElementalDofList, model_part.GetProcessInfo());

    std::vector<int> ids{23, 74, 55};
    std::vector<int> upstream_ids{23, 55, 67};
    for (int i = 0; i < 3; i++){
        ElementalDofList[i]->SetEquationId(ids[i]);
        UpstreamElementalDofList[i]->SetEquationId(upstream_ids[i]);
    }

    FindNodalNeighboursProcess find_nodal_neighbours_process(model_part);
    find_nodal_neighbours_process.Execute();
    pElement->Initialize(model_part.GetProcessInfo());

    // double local_mach_number =
    //     PotentialFlowUtilities::ComputePerturbationLocalMachNumber<2, 3>(
    //         *pElement, model_part.GetProcessInfo());


    // double local_mach_number_upstream =
    //     PotentialFlowUtilities::ComputePerturbationLocalMachNumber<2, 3>(
    //         *pUpstreamElement, model_part.GetProcessInfo());

    // KRATOS_WATCH(local_mach_number)
    // KRATOS_WATCH(local_mach_number_upstream)

    Vector RHS_original = ZeroVector(number_of_nodes);
    Matrix LHS_original = ZeroMatrix(number_of_nodes, number_of_nodes);
    Matrix LHS_finite_diference = ZeroMatrix(number_of_nodes + 1, number_of_nodes + 1);
    Matrix LHS_analytical = ZeroMatrix(number_of_nodes + 1, number_of_nodes + 1);

    // Compute original RHS and LHS
    pElement->CalculateLocalSystem(LHS_original, RHS_original, model_part.GetProcessInfo());
    // KRATOS_WATCH(LHS_original)
    // KRATOS_WATCH(RHS_original)

    maximum_potential = 500.0;
    minimum_potential = 0.0;


    potential = {minimum_potential, maximum_potential, maximum_potential};
    AssignPerturbationPotentialsToTransonicElement(*pElement, potential);

    upstream_potential= {minimum_potential, maximum_potential, minimum_potential};
    AssignPerturbationPotentialsToTransonicElement(*pUpstreamElement, upstream_potential);

    // local_mach_number =
    //     PotentialFlowUtilities::ComputePerturbationLocalMachNumber<2, 3>(
    //         *pElement, model_part.GetProcessInfo());


    // local_mach_number_upstream =
    //     PotentialFlowUtilities::ComputePerturbationLocalMachNumber<2, 3>(
    //         *pUpstreamElement, model_part.GetProcessInfo());

    // KRATOS_WATCH(local_mach_number)
    // KRATOS_WATCH(local_mach_number_upstream)

    // Compute original RHS and LHS
    pElement->CalculateLocalSystem(LHS_original, RHS_original, model_part.GetProcessInfo());
    // KRATOS_WATCH(LHS_original)
    // KRATOS_WATCH(RHS_original)

    // double delta = 1e-3;
    // for(unsigned int i = 0; i < number_of_nodes; i++){
    //     // Pinging
    //     pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += delta;

    //     Vector RHS_pinged = ZeroVector(number_of_nodes + 1);
    //     Matrix LHS_pinged = ZeroMatrix(number_of_nodes + 1, number_of_nodes + 1);
    //     // Compute pinged LHS and RHS
    //     pElement->CalculateLocalSystem(LHS_pinged, RHS_pinged, model_part.GetProcessInfo());

    //     for(unsigned int k = 0; k < number_of_nodes + 1; k++){
    //         // Compute the finite difference estimate of the sensitivity
    //         LHS_finite_diference( k, i) = -(RHS_pinged(k)-RHS_original(k)) / delta;
    //         // Compute the average of the original and pinged analytic sensitivities
    //         LHS_analytical( k, i) = 0.5 * (LHS_original(k,i) + LHS_pinged(k,i));
    //     }

    //     // Unpinging
    //     pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= delta;
    // }

    // std::cout.precision(16);
    // KRATOS_WATCH(LHS_analytical)
    // KRATOS_WATCH(LHS_finite_diference)
    // for (unsigned int i = 0; i < LHS_finite_diference.size1(); i++) {
    //     for (unsigned int j = 0; j < LHS_finite_diference.size2(); j++) {
    //         KRATOS_CHECK_NEAR(LHS_finite_diference(i,j), LHS_analytical(i,j), 1e-10);
    //     }
    // }
}

/** Checks the TransonicPerturbationPotentialFlowElement.
 * Checks the LHS computation.
 */
KRATOS_TEST_CASE_IN_SUITE(TransonicPerturbationPotentialFlowElementLHSSuperSonicDecelerating, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);
    std::array<double, 3> potential{1.0, 120.0, 180.0};
    AssignPerturbationPotentialsToTransonicElement(*pElement, potential);

    GenerateTestingTransonicUpstreamElement(model_part);
    Element::Pointer pUpstreamElement = model_part.pGetElement(2);
    std::array<double, 3> upstream_potential{1.0, 180.0, 51.0};
    AssignPerturbationPotentialsToTransonicElement(*pUpstreamElement, upstream_potential);

    for (auto& r_node : model_part.Nodes()){
        r_node.AddDof(VELOCITY_POTENTIAL);
    }

    Element::DofsVectorType ElementalDofList, UpstreamElementalDofList;
    pElement->GetDofList(ElementalDofList, model_part.GetProcessInfo());
    pUpstreamElement->GetDofList(UpstreamElementalDofList, model_part.GetProcessInfo());

    std::vector<int> ids{23, 74, 55};
    std::vector<int> upstream_ids{23, 55, 67};
    for (int i = 0; i < 3; i++){
        ElementalDofList[i]->SetEquationId(ids[i]);
        UpstreamElementalDofList[i]->SetEquationId(upstream_ids[i]);
    }

    FindNodalNeighboursProcess find_nodal_neighbours_process(model_part);
    find_nodal_neighbours_process.Execute();
    pElement->Initialize(model_part.GetProcessInfo());

    // Compute LHS
    Matrix LHS = ZeroMatrix(3, 3);

    pElement->CalculateLeftHandSide(LHS, model_part.GetProcessInfo());

    std::array<double, 16> reference{ 0.03486940937962629,-0.1220506719403047, 0.1570735000967481,-0.06989223753606974,
                                     -0.112001164475407,0.6333235380330079,-0.5782387640357952,0.0569163904781943,
                                      0.07713175509578073,-0.5112728660927033,0.4211652639390471,0.01297584705787544,
                                      0.0, 0.0, 0.0, 0.0};

    for (unsigned int i = 0; i < LHS.size1(); i++) {
        for (unsigned int j = 0; j < LHS.size2(); j++) {
            KRATOS_CHECK_NEAR(LHS(i, j), reference[i * 4 + j], 1e-15);
        }
    }
}

/** Checks the TransonicPerturbationPotentialFlowElement.
 * Tests the LHS computation.
 */
KRATOS_TEST_CASE_IN_SUITE(PingTransonicPerturbationPotentialFlowElementLHSSuperSonicDecelerating, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

GenerateTransonicPerturbationElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);
    const unsigned int number_of_nodes = pElement->GetGeometry().size();
    std::array<double, 3> potential{1.0, 120.0, 180.0};
    AssignPerturbationPotentialsToTransonicElement(*pElement, potential);

    GenerateTestingTransonicUpstreamElement(model_part);
    Element::Pointer pUpstreamElement = model_part.pGetElement(2);
    std::array<double, 3> upstream_potential{1.0, 180.0, 51.0};
    AssignPerturbationPotentialsToTransonicElement(*pUpstreamElement, upstream_potential);

    for (auto& r_node : model_part.Nodes()){
        r_node.AddDof(VELOCITY_POTENTIAL);
    }

    Element::DofsVectorType ElementalDofList, UpstreamElementalDofList;
    pElement->GetDofList(ElementalDofList, model_part.GetProcessInfo());
    pUpstreamElement->GetDofList(UpstreamElementalDofList, model_part.GetProcessInfo());

    std::vector<int> ids{23, 74, 55};
    std::vector<int> upstream_ids{23, 55, 67};
    for (int i = 0; i < 3; i++){
        ElementalDofList[i]->SetEquationId(ids[i]);
        UpstreamElementalDofList[i]->SetEquationId(upstream_ids[i]);
    }

    FindNodalNeighboursProcess find_nodal_neighbours_process(model_part);
    find_nodal_neighbours_process.Execute();
    pElement->Initialize(model_part.GetProcessInfo());

    Vector RHS_original = ZeroVector(number_of_nodes);
    Matrix LHS_original = ZeroMatrix(number_of_nodes, number_of_nodes);
    Matrix LHS_finite_diference = ZeroMatrix(number_of_nodes + 1, number_of_nodes + 1);
    Matrix LHS_analytical = ZeroMatrix(number_of_nodes + 1, number_of_nodes + 1);

    // Compute original RHS and LHS
    pElement->CalculateLocalSystem(LHS_original, RHS_original, model_part.GetProcessInfo());
    KRATOS_WATCH(LHS_original)
    KRATOS_WATCH(RHS_original)

    double delta = 1e-3;
    for(unsigned int i = 0; i < number_of_nodes; i++){
        // Pinging
        pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += delta;

        Vector RHS_pinged = ZeroVector(number_of_nodes + 1);
        Matrix LHS_pinged = ZeroMatrix(number_of_nodes + 1, number_of_nodes + 1);
        // Compute pinged LHS and RHS
        pElement->CalculateLocalSystem(LHS_pinged, RHS_pinged, model_part.GetProcessInfo());

        for(unsigned int k = 0; k < number_of_nodes + 1; k++){
            // Compute the finite difference estimate of the sensitivity
            LHS_finite_diference( k, i) = -(RHS_pinged(k)-RHS_original(k)) / delta;
            // Compute the average of the original and pinged analytic sensitivities
            LHS_analytical( k, i) = 0.5 * (LHS_original(k,i) + LHS_pinged(k,i));
        }

        // Unpinging
        pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= delta;
    }

    std::cout.precision(16);
    KRATOS_WATCH(LHS_analytical)
    KRATOS_WATCH(LHS_finite_diference)
    for (unsigned int i = 0; i < LHS_finite_diference.size1(); i++) {
        for (unsigned int j = 0; j < LHS_finite_diference.size2(); j++) {
            KRATOS_CHECK_NEAR(LHS_finite_diference(i,j), LHS_analytical(i,j), 1e-10);
        }
    }
}

/** Checks the TransonicPerturbationPotentialFlowElement.
 * Checks the LHS computation.
 */
KRATOS_TEST_CASE_IN_SUITE(TransonicPerturbationPotentialFlowElementLHSInlet, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);
    AssignPotentialsToNormalTransonicPerturbationElement(pElement);
    pElement->Initialize(model_part.GetProcessInfo());

    // Compute LHS
    Matrix LHS = ZeroMatrix(3, 3);

    pElement->CalculateLeftHandSide(LHS, model_part.GetProcessInfo());

    std::array<double, 16> reference{ 0.06114278464441542,-0.1306215050744058, 0.06947872042999037,
                                     -0.1306215050744058, 0.6710758508914103,-0.5404543458170046,
                                      0.06947872042999037,-0.5404543458170046,0.4709756253870142};

    for (unsigned int i = 0; i < LHS.size1(); i++) {
        for (unsigned int j = 0; j < LHS.size2(); j++) {
            KRATOS_CHECK_NEAR(LHS(i, j), reference[i * 3 + j], 1e-16);
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
    Element::Pointer pElement = model_part.pGetElement(1);
    pElement->Initialize(model_part.GetProcessInfo());
    const unsigned int number_of_nodes = pElement->GetGeometry().size();

    AssignPotentialsToNormalTransonicPerturbationElement(pElement);

    Vector RHS_original = ZeroVector(number_of_nodes);
    Matrix LHS_original = ZeroMatrix(number_of_nodes, number_of_nodes);
    Matrix LHS_finite_diference = ZeroMatrix(number_of_nodes, number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(number_of_nodes, number_of_nodes);

    // Compute original RHS and LHS
    pElement->CalculateLocalSystem(LHS_original, RHS_original, model_part.GetProcessInfo());

    double delta = 1e-3;
    for(unsigned int i = 0; i < number_of_nodes; i++){
        // Pinging
        pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += delta;

        Vector RHS_pinged = ZeroVector(number_of_nodes);
        Matrix LHS_pinged = ZeroMatrix(number_of_nodes, number_of_nodes);
        // Compute pinged LHS and RHS
        pElement->CalculateLocalSystem(LHS_pinged, RHS_pinged, model_part.GetProcessInfo());

        for(unsigned int k = 0; k < number_of_nodes; k++){
            // Compute the finite difference estimate of the sensitivity
            LHS_finite_diference( k, i) = -(RHS_pinged(k)-RHS_original(k)) / delta;
            // Compute the average of the original and pinged analytic sensitivities
            LHS_analytical( k, i) = 0.5 * (LHS_original(k,i) + LHS_pinged(k,i));
        }

        // Unpinging
        pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= delta;
    }

    for (unsigned int i = 0; i < LHS_finite_diference.size1(); i++) {
        for (unsigned int j = 0; j < LHS_finite_diference.size2(); j++) {
            KRATOS_CHECK_NEAR(LHS_finite_diference(i,j), LHS_analytical(i,j), 1e-10);
        }
    }
}

BoundedVector<double,3> AssignDistancesToPerturbationTransonicElement()
{
    BoundedVector<double,3> distances;
    distances(0) = 1.0;
    distances(1) = -1.0;
    distances(2) = -1.0;
    return distances;
}

void AssignPotentialsToWakeTransonicPerturbationElement(Element::Pointer pElement, const array_1d<double, 3>& rDistances)
{
    const std::array<double, 3> potential{1.0, 100.0, 150.0};

    for (unsigned int i = 0; i < 3; i++){
        if (rDistances(i) > 0.0)
            pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = potential[i];
        else
            pElement->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) = potential[i];
    }
    for (unsigned int i = 0; i < 3; i++){
        if (rDistances(i) < 0.0)
            pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = potential[i]+5;
        else
            pElement->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) = potential[i]+5;
    }
}

KRATOS_TEST_CASE_IN_SUITE(WakeTransonicPerturbationPotentialFlowElementRHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);

    BoundedVector<double,3> distances = AssignDistancesToPerturbationTransonicElement();

    pElement->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    pElement->GetValue(WAKE) = true;

    AssignPotentialsToWakeTransonicPerturbationElement(pElement, distances);

    // Compute RHS and LHS
    Vector RHS = ZeroVector(6);

    pElement->CalculateRightHandSide(RHS, model_part.GetProcessInfo());

    std::vector<double> reference{146.2643261263345,-0,-0,-0,-122.1426284341492,-24.12169769218525};

    KRATOS_CHECK_VECTOR_NEAR(RHS, reference, 1e-13);
}

KRATOS_TEST_CASE_IN_SUITE(WakeTransonicPerturbationPotentialFlowElementLHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);

    BoundedVector<double,3> distances = AssignDistancesToPerturbationTransonicElement();

    pElement->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    pElement->GetValue(WAKE) = true;

    AssignPotentialsToWakeTransonicPerturbationElement(pElement, distances);

    // Compute LHS
    Matrix LHS = ZeroMatrix(6, 6);

    pElement->CalculateLeftHandSide(LHS, model_part.GetProcessInfo());

    // Check the LHS values
    std::array<double,36> reference{0.06114278464441542,-0.1306215050744058,0.06947872042999037,0,0,0,
    -0.1306215050744058,0.6710758508914103,-0.5404543458170046,0.1306215050744058,-0.6710758508914103,0.5404543458170046,
    0.06947872042999037,-0.5404543458170046,0.4709756253870142,-0.06947872042999037,0.5404543458170046,-0.4709756253870142,
    -0.06114278464441542,0.1306215050744058,-0.06947872042999037,0.06114278464441542,-0.1306215050744058,0.06947872042999037,
    0,0,0,-0.1306215050744058,0.6710758508914103,-0.5404543458170046,
    0,0,0,0.06947872042999037,-0.5404543458170046,0.4709756253870142};

    for (unsigned int i = 0; i < LHS.size1(); i++) {
        for (unsigned int j = 0; j < LHS.size2(); j++) {
            KRATOS_CHECK_NEAR(LHS(i, j), reference[6 * i + j], 1e-16);
        }
    }
}

/** Checks the TransonicPotentialFlowElement element.
* Checks the EquationIdVector.
*/
KRATOS_TEST_CASE_IN_SUITE(TransonicPotentialFlowElementEquationIdVector, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);
    std::array<double, 3> potential{1.0, 120.0, 180.0};
    AssignPerturbationPotentialsToTransonicElement(*pElement, potential);

    GenerateTestingTransonicUpstreamElement(model_part);
    Element::Pointer pUpstreamElement = model_part.pGetElement(2);
    std::array<double, 3> upstream_potential{1.0, 180.0, 90.0};
    AssignPerturbationPotentialsToTransonicElement(*pUpstreamElement, upstream_potential);

    FindNodalNeighboursProcess find_nodal_neighbours_process(model_part);
    find_nodal_neighbours_process.Execute();
    pElement->Initialize(model_part.GetProcessInfo());

    for (auto& r_node : model_part.Nodes()){
        r_node.AddDof(VELOCITY_POTENTIAL);
    }

    Element::DofsVectorType ElementalDofList, UpstreamElementalDofList;
    pElement->GetDofList(ElementalDofList, model_part.GetProcessInfo());
    pUpstreamElement->GetDofList(UpstreamElementalDofList, model_part.GetProcessInfo());

    std::vector<int> ids{23, 74, 55};
    std::vector<int> upstream_ids{23, 55, 67};
    for (int i = 0; i < 3; i++){
        ElementalDofList[i]->SetEquationId(ids[i]);
        UpstreamElementalDofList[i]->SetEquationId(upstream_ids[i]);
    }

    Element::EquationIdVectorType EquationIdVector;
    pElement->EquationIdVector(EquationIdVector, model_part.GetProcessInfo());

    std::vector<double> reference{23, 74, 55, 67};
    KRATOS_CHECK_VECTOR_NEAR(EquationIdVector, reference, 1e-16);
}

/** Checks the TransonicPotentialFlowElement element.
* Checks the EquationIdVector for the Wake.
*/
KRATOS_TEST_CASE_IN_SUITE(TransonicPotentialFlowElementEquationIdVectorWake, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);
    pElement->SetValue(WAKE, true);
    BoundedVector<double, 3> distances = AssignDistancesToPerturbationTransonicElement();
    pElement->SetValue(WAKE_ELEMENTAL_DISTANCES, distances);

    GenerateTestingTransonicUpstreamElement(model_part);

    FindNodalNeighboursProcess find_nodal_neighbours_process(model_part);
    find_nodal_neighbours_process.Execute();
    pElement->Initialize(model_part.GetProcessInfo());

    for (unsigned int i = 0; i < 3; i++) {
        pElement->GetGeometry()[i].AddDof(VELOCITY_POTENTIAL);
        pElement->GetGeometry()[i].AddDof(AUXILIARY_VELOCITY_POTENTIAL);
    }

    Element::DofsVectorType ElementalDofList;
    pElement->GetDofList(ElementalDofList, model_part.GetProcessInfo());

    std::vector<int> ids{15, 366, 24, 98, 103, 254};
    for (int i = 0; i < 6; i++){
        ElementalDofList[i]->SetEquationId(ids[i]);
    }

    Element::EquationIdVectorType EquationIdVector;
    pElement->EquationIdVector(EquationIdVector, model_part.GetProcessInfo());

    std::vector<double> reference{15, 366, 24, 98, 103, 254};
    KRATOS_CHECK_VECTOR_NEAR(EquationIdVector, reference, 1e-16);
}

} // namespace Testing
} // namespace Kratos.
