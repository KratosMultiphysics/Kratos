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
#include "tests/cpp_tests/compressible_potential_flow_fast_suite.h"
#include "compressible_potential_flow_application_variables.h"
#include "fluid_dynamics_application_variables.h"
#include "custom_elements/compressible_perturbation_potential_flow_element.h"
#include "custom_utilities/potential_flow_utilities.h"

namespace Kratos {
namespace Testing {

typedef ModelPart::IndexType IndexType;
typedef ModelPart::NodeIterator NodeIteratorType;

void GenerateCompressiblePerturbationElement(ModelPart& rModelPart) {
    // Variables addition
    rModelPart.AddNodalSolutionStepVariable(VELOCITY_POTENTIAL);
    rModelPart.AddNodalSolutionStepVariable(AUXILIARY_VELOCITY_POTENTIAL);

    // Set the element properties
    Properties::Pointer pElemProp = rModelPart.CreateNewProperties(0);
    rModelPart.GetProcessInfo()[FREE_STREAM_DENSITY] = 1.225;
    rModelPart.GetProcessInfo()[FREE_STREAM_MACH] = 0.6;
    rModelPart.GetProcessInfo()[HEAT_CAPACITY_RATIO] = 1.4;
    rModelPart.GetProcessInfo()[SOUND_VELOCITY] = 340.3;
    rModelPart.GetProcessInfo()[MACH_LIMIT] = 0.94;

    BoundedVector<double, 3> free_stream_velocity = ZeroVector(3);
    free_stream_velocity(0) = rModelPart.GetProcessInfo().GetValue(FREE_STREAM_MACH) * rModelPart.GetProcessInfo().GetValue(SOUND_VELOCITY);
    rModelPart.GetProcessInfo()[FREE_STREAM_VELOCITY] = free_stream_velocity;

    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> elemNodes{1, 2, 3};
    rModelPart.CreateNewElement("CompressiblePerturbationPotentialFlowElement2D3N", 1, elemNodes, pElemProp);
}

void AssignPotentialsToNormalCompressiblePerturbationElement(Element::Pointer pElement, const std::array<double, 3> rPotential)
{
    for (unsigned int i = 0; i < 3; i++)
        pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = rPotential[i];
}

BoundedVector<double,3> AssignDistancesToPerturbationCompressibleElement()
{
    BoundedVector<double,3> distances;
    distances(0) = 1.0;
    distances(1) = -1.0;
    distances(2) = -1.0;
    return distances;
}

void AssignPotentialsToWakeCompressiblePerturbationElement(Element::Pointer pElement, const array_1d<double, 3>& rDistances, const std::array<double, 6>& rPotential)
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

void PrintTestWakeMatrixPretty(Matrix& rMatrix){
    std::cout.precision(5);
    std::cout << std::scientific;
    std::cout << std::showpos;
    std::cout << std::endl;
    for(unsigned int row = 0; row < rMatrix.size1(); ++row){
        for(unsigned int column = 0; column < rMatrix.size2(); column++){
            if(column == 2 || column == 5){
                std::cout << " " << rMatrix(row, column) << " |";
            }
            else{
                std::cout << " " << rMatrix(row, column) << " ";
            }
        }

        std::cout << std::endl;

        if(row ==2|| row == 5){
            for(unsigned int j = 0; j < 14*rMatrix.size1(); j++){
            std::cout << "_" ;
            }
            std::cout << " " << std::endl;
        }
        else{
            for(unsigned int i = 0; i < 3; i++){
                for(unsigned int j = 0; j < 14*3; j++){
                    std::cout << " " ;
                }
                if(i != 2){
                    std::cout << "|" ;
                }
            }
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void PrintTestElementInfo(ModelPart& rModelPart){
    const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
    Element::Pointer p_element = rModelPart.pGetElement(1);
    array_1d<double, 2> perturbed_velocity = PotentialFlowUtilities::ComputePerturbedVelocity<2,3>(*p_element, r_current_process_info);
    const double local_velocity_squared = inner_prod(perturbed_velocity, perturbed_velocity);
    const double local_mach_squared = PotentialFlowUtilities::ComputeLocalMachNumberSquared<2, 3>(perturbed_velocity, r_current_process_info);
    const double max_velocity_squared = PotentialFlowUtilities::ComputeMaximumVelocitySquared<2, 3>(r_current_process_info);

    std::cout.precision(16);
    KRATOS_WATCH(perturbed_velocity)
    KRATOS_WATCH(std::sqrt(max_velocity_squared))
    KRATOS_WATCH(std::sqrt(local_velocity_squared))
    KRATOS_WATCH(local_mach_squared)
}


/** Checks the CompressiblePerturbationPotentialFlowElement.
 * Checks the RHS computation.
 */
KRATOS_TEST_CASE_IN_SUITE(CompressiblePerturbationPotentialFlowElementRHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    std::array<double, 3> potential{1.0, 20.0, 50.0};
    AssignPotentialsToNormalCompressiblePerturbationElement(p_element, potential);

    // Compute RHS
    Vector RHS = ZeroVector(3);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateRightHandSide(RHS, r_current_process_info);

    std::vector<double> reference{131.4361747323354,-113.768439084114,-17.66773564822145};

    KRATOS_EXPECT_VECTOR_NEAR(RHS, reference, 1e-13);
}

KRATOS_TEST_CASE_IN_SUITE(CompressiblePerturbationPotentialFlowElementRHSClamping, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    std::array<double, 3> potential{1.0, 220.0, 250.0};
    AssignPotentialsToNormalCompressiblePerturbationElement(p_element, potential);

    // Compute RHS
    Vector RHS = ZeroVector(3);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateRightHandSide(RHS, r_current_process_info);

    std::vector<double> reference{205.3219372530133,-190.7662916232804,-14.55564562973297};

    KRATOS_EXPECT_VECTOR_NEAR(RHS, reference, 1e-13);
}

// /** Checks the CompressiblePerturbationPotentialFlowElement.
//  * Checks the LHS computation.
//  */
KRATOS_TEST_CASE_IN_SUITE(CompressiblePerturbationPotentialFlowElementLHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    std::array<double, 3> potential{1.0, 20.0, 50.0};
    AssignPotentialsToNormalCompressiblePerturbationElement(p_element, potential);

    // Compute LHS
    Matrix LHS = ZeroMatrix(3, 3);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateLeftHandSide(LHS, r_current_process_info);

    std::array<double, 9> reference{ 0.3316096612183497,-0.3661980912374865,0.03458843001913683,
                                    -0.3661980912374865,0.9850616437217248,-0.6188635524842382,
                                    0.03458843001913683,-0.6188635524842382,0.5842751224651014};

    for (unsigned int i = 0; i < LHS.size1(); i++) {
        for (unsigned int j = 0; j < LHS.size2(); j++) {
            KRATOS_EXPECT_NEAR(LHS(i, j), reference[i * 3 + j], 1e-16);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(CompressiblePerturbationPotentialFlowElementLHSClamping, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    std::array<double, 3> potential{1.0, 220.0, 250.0};
    AssignPotentialsToNormalCompressiblePerturbationElement(p_element, potential);

    // Compute LHS
    Matrix LHS = ZeroMatrix(3, 3);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateLeftHandSide(LHS, r_current_process_info);

    std::array<double, 9> reference{ 0.4851881876577658,-0.4851881876577658,0,
                                    -0.4851881876577658,0.9703763753155316,-0.4851881876577658,
                                    0,-0.4851881876577658,0.4851881876577658};

    for (unsigned int i = 0; i < LHS.size1(); i++) {
        for (unsigned int j = 0; j < LHS.size2(); j++) {
            KRATOS_EXPECT_NEAR(LHS(i, j), reference[i * 3 + j], 1e-16);
        }
    }
}

void ComputeElementalSensitivitiesMatrixRow(ModelPart& rModelPart, double delta, unsigned int row, Matrix& rLHS_original, Vector& rRHS_original, Matrix& rLHS_finite_diference, Matrix& rLHS_analytical){
    Element::Pointer p_element = rModelPart.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    // Compute pinged LHS and RHS
    Vector RHS_pinged = ZeroVector(number_of_nodes);
    Matrix LHS_pinged = ZeroMatrix(number_of_nodes, number_of_nodes);
    const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
    p_element->CalculateLocalSystem(LHS_pinged, RHS_pinged, r_current_process_info);

    for(unsigned int k = 0; k < rLHS_original.size2(); k++){
        // Compute the finite difference estimate of the sensitivity
        rLHS_finite_diference( k, row) = -(RHS_pinged(k)-rRHS_original(k)) / delta;
        // Compute the average of the original and pinged analytic sensitivities
        rLHS_analytical( k, row) = 0.5 * (rLHS_original(k,row) + LHS_pinged(k,row));
    }

}

void ComputeElementalSensitivities(ModelPart& rModelPart, Matrix& rLHS_finite_diference, Matrix& rLHS_analytical, const std::array<double, 3> rPotential){
    Element::Pointer p_element = rModelPart.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    AssignPotentialsToNormalCompressiblePerturbationElement(p_element, rPotential);

    // Compute original RHS and LHS
    Vector RHS_original = ZeroVector(number_of_nodes);
    Matrix LHS_original = ZeroMatrix(number_of_nodes, number_of_nodes);
    const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
    p_element->CalculateLocalSystem(LHS_original, RHS_original, r_current_process_info);

    double delta = 1e-3;
    for(unsigned int i = 0; i < number_of_nodes; i++){
        // Pinging
        p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += delta;

        ComputeElementalSensitivitiesMatrixRow(rModelPart, delta, i, LHS_original, RHS_original, rLHS_finite_diference, rLHS_analytical);

        // Unpinging
        p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= delta;
    }
}

void ComputeWakeElementalSensitivities(ModelPart& rModelPart, Matrix& rLHS_finite_diference, Matrix& rLHS_analytical, const std::array<double, 6> rPotential){
    Element::Pointer p_element = rModelPart.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    BoundedVector<double,3> distances = AssignDistancesToPerturbationCompressibleElement();
    p_element->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    p_element->GetValue(WAKE) = true;

    AssignPotentialsToWakeCompressiblePerturbationElement(p_element, distances, rPotential);

    // Compute original RHS and LHS
    Vector RHS_original = ZeroVector(2*number_of_nodes);
    Matrix LHS_original = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);
    const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
    p_element->CalculateLocalSystem(LHS_original, RHS_original, r_current_process_info);

    double delta = 1e-3;
    for(unsigned int i = 0; i < 2*number_of_nodes; i++){
        if(i < number_of_nodes){
            // Pinging
            if (distances(i) > 0.0)
                p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += delta;
            else
                p_element->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) += delta;

            ComputeElementalSensitivitiesMatrixRow(rModelPart, delta, i, LHS_original, RHS_original, rLHS_finite_diference, rLHS_analytical);

            // Unpinging
            if (distances(i) > 0.0)
                p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= delta;
            else
                p_element->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) -= delta;
        }
        else{
            // Pinging
            if (distances(i-number_of_nodes) > 0.0)
                p_element->GetGeometry()[i-number_of_nodes].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) += delta;
            else
                p_element->GetGeometry()[i-number_of_nodes].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += delta;

            ComputeElementalSensitivitiesMatrixRow(rModelPart, delta, i, LHS_original, RHS_original, rLHS_finite_diference, rLHS_analytical);

            // Unpinging
            if (distances(i-number_of_nodes) > 0.0)
                p_element->GetGeometry()[i-number_of_nodes].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) -= delta;
            else
                p_element->GetGeometry()[i-number_of_nodes].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= delta;
        }
    }
}

/** Checks the CompressiblePerturbationPotentialFlowElement.
 * Tests the LHS computation.
 */
KRATOS_TEST_CASE_IN_SUITE(PingCompressiblePerturbationPotentialFlowElementLHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    std::array<double, 3> potential{1.0, 20.0, 50.0};

    Matrix LHS_finite_diference = ZeroMatrix(number_of_nodes, number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(number_of_nodes, number_of_nodes);

    ComputeElementalSensitivities(model_part, LHS_finite_diference, LHS_analytical, potential);

    KRATOS_EXPECT_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(PingCompressiblePerturbationPotentialFlowElementLHSClamping, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    std::array<double, 3> potential{1.2495, 94.1948, 182.149583};

    Matrix LHS_finite_diference = ZeroMatrix(number_of_nodes, number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(number_of_nodes, number_of_nodes);

    ComputeElementalSensitivities(model_part, LHS_finite_diference, LHS_analytical, potential);

    KRATOS_EXPECT_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(WakeCompressiblePerturbationPotentialFlowElementRHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    BoundedVector<double,3> distances = AssignDistancesToPerturbationCompressibleElement();

    p_element->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    p_element->GetValue(WAKE) = true;

    const std::array<double, 6> potential{1.0, 31.0, 150.0, 6.0, 75.0, 55.0};
    AssignPotentialsToWakeCompressiblePerturbationElement(p_element, distances, potential);

    // Compute RHS and LHS
    Vector RHS = ZeroVector(6);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateRightHandSide(RHS, r_current_process_info);

    std::vector<double> reference{127.1146544469925,109.025,-85.1375,23.8875,-154.8303022595422,10.56213263248122};

    KRATOS_EXPECT_VECTOR_NEAR(RHS, reference, 1e-13);
}

KRATOS_TEST_CASE_IN_SUITE(WakeCompressiblePerturbationPotentialFlowElementRHSClamping, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    BoundedVector<double,3> distances = AssignDistancesToPerturbationCompressibleElement();

    p_element->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    p_element->GetValue(WAKE) = true;

    const std::array<double, 6> potential{1.0, 151.0, 190.0, 6.0, 165.0, 195.0};
    AssignPotentialsToWakeCompressiblePerturbationElement(p_element, distances, potential);

    // Compute RHS and LHS
    Vector RHS = ZeroVector(6);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateRightHandSide(RHS, r_current_process_info);

    std::vector<double> reference{171.8439523046275,11.025,-5.5125,5.5125,-161.6550003638144,-14.55564562973297};

    KRATOS_EXPECT_VECTOR_NEAR(RHS, reference, 1e-13);
}

KRATOS_TEST_CASE_IN_SUITE(WakeStructureCompressiblePerturbationPotentialFlowElementRHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    BoundedVector<double,3> distances = AssignDistancesToPerturbationCompressibleElement();

    p_element->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    p_element->GetValue(WAKE) = true;
    p_element->Set(STRUCTURE);
    p_element->GetGeometry()[number_of_nodes-1].SetValue(TRAILING_EDGE, true);

    const std::array<double, 6> potential{1.0, 31.0, 150.0, 6.0, 75.0, 55.0};
    AssignPotentialsToWakeCompressiblePerturbationElement(p_element, distances, potential);

    // Compute RHS and LHS
    Vector RHS = ZeroVector(6);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateRightHandSide(RHS, r_current_process_info);

    std::vector<double> reference{127.1146544469925,109.025,-16.14852237508765,23.8875,-154.8303022595422,7.921599474360912};

    KRATOS_EXPECT_VECTOR_NEAR(RHS, reference, 1e-13);
}

KRATOS_TEST_CASE_IN_SUITE(WakeStructureCompressiblePerturbationPotentialFlowElementRHSClamping, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    BoundedVector<double,3> distances = AssignDistancesToPerturbationCompressibleElement();

    p_element->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    p_element->GetValue(WAKE) = true;
    p_element->Set(STRUCTURE);
    p_element->GetGeometry()[number_of_nodes-1].SetValue(TRAILING_EDGE, true);

    const std::array<double, 6> potential{1.0, 151.0, 190.0, 6.0, 165.0, 195.0};
    AssignPotentialsToWakeCompressiblePerturbationElement(p_element, distances, potential);

    // Compute RHS and LHS
    Vector RHS = ZeroVector(6);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateRightHandSide(RHS, r_current_process_info);

    std::vector<double> reference{171.8439523046275,11.025,-4.730584829663217,5.5125,-161.6550003638144,-10.91673422229973};

    KRATOS_EXPECT_VECTOR_NEAR(RHS, reference, 1e-13);
}

KRATOS_TEST_CASE_IN_SUITE(WakeCompressiblePerturbationPotentialFlowElementLHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    BoundedVector<double,3> distances = AssignDistancesToPerturbationCompressibleElement();

    p_element->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    p_element->GetValue(WAKE) = true;

    const std::array<double, 6> potential{1.0, 31.0, 150.0, 6.0, 75.0, 55.0};
    AssignPotentialsToWakeCompressiblePerturbationElement(p_element, distances, potential);

    // Compute LHS
    Matrix LHS = ZeroMatrix(6, 6);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateLeftHandSide(LHS, r_current_process_info);

    // Check the LHS values
    std::array<double,36> reference{0.2730300317674839,-0.4101190902695764,0.1370890585020925,0,0,0,-0.6125,1.225,-0.6125,0.6125,-1.225,0.6125,0,-0.6125,0.6125,-0,0.6125,-0.6125,-0.6125,0.6125,-0,0.6125,-0.6125,0,0,0,0,-0.1405503745970895,0.6402833143676492,-0.4997329397705597,0,0,0,-0.02643811017306578,-0.4997329397705597,0.5261710499436255};

    for (unsigned int i = 0; i < LHS.size1(); i++) {
        for (unsigned int j = 0; j < LHS.size2(); j++) {
            KRATOS_EXPECT_NEAR(LHS(i, j), reference[6 * i + j], 1e-16);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(WakeCompressiblePerturbationPotentialFlowElementLHSClamping, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    BoundedVector<double,3> distances = AssignDistancesToPerturbationCompressibleElement();

    p_element->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    p_element->GetValue(WAKE) = true;

    const std::array<double, 6> potential{1.0, 151.0, 190.0, 6.0, 165.0, 195.0};
    AssignPotentialsToWakeCompressiblePerturbationElement(p_element, distances, potential);

    // Compute LHS
    Matrix LHS = ZeroMatrix(6, 6);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateLeftHandSide(LHS, r_current_process_info);

    // Check the LHS values
    std::array<double,36> reference{0.4851881876577658,-0.4851881876577658,0,0,0,0,-0.6125,1.225,-0.6125,0.6125,
    -1.225,0.6125,0,-0.6125,0.6125,-0,0.6125,-0.6125,-0.6125,0.6125,-0,0.6125,-0.6125,0,0,0,0,-0.4851881876577658,0.9703763753155316,-0.4851881876577658,0,0,0,0,-0.4851881876577658,0.4851881876577658};

    for (unsigned int i = 0; i < LHS.size1(); i++) {
        for (unsigned int j = 0; j < LHS.size2(); j++) {
            KRATOS_EXPECT_NEAR(LHS(i, j), reference[6 * i + j], 1e-16);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(WakeStructureCompressiblePerturbationPotentialFlowElementLHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    BoundedVector<double,3> distances = AssignDistancesToPerturbationCompressibleElement();

    p_element->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    p_element->GetValue(WAKE) = true;
    p_element->Set(STRUCTURE);
    p_element->GetGeometry()[number_of_nodes-1].SetValue(TRAILING_EDGE, true);

    const std::array<double, 6> potential{1.0, 31.0, 150.0, 6.0, 75.0, 55.0};
    AssignPotentialsToWakeCompressiblePerturbationElement(p_element, distances, potential);

    // Compute LHS
    Matrix LHS = ZeroMatrix(6, 6);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateLeftHandSide(LHS, r_current_process_info);

    // Check the LHS values
    std::array<double,36> reference{0.2730300317674839,-0.4101190902695764,0.1370890585020925,0,0,0,-0.6125,1.225,
    -0.6125,0.6125,-1.225,0.6125,0.03427226462552314,-0.1525584723345968,0.1182862077090737,0,0,0,-0.6125,0.6125,
    -0,0.6125,-0.6125,0,0,0,0,-0.1405503745970895,0.6402833143676492,-0.4997329397705597,0,0,0,-0.01982858262979934,-0.3747997048279197,0.3946282874577191};

    for (unsigned int i = 0; i < LHS.size1(); i++) {
        for (unsigned int j = 0; j < LHS.size2(); j++) {
            KRATOS_EXPECT_NEAR(LHS(i, j), reference[6 * i + j], 1e-16);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(WakeStructureCompressiblePerturbationPotentialFlowElementLHSClamping, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    BoundedVector<double,3> distances = AssignDistancesToPerturbationCompressibleElement();

    p_element->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    p_element->GetValue(WAKE) = true;
    p_element->Set(STRUCTURE);
    p_element->GetGeometry()[number_of_nodes-1].SetValue(TRAILING_EDGE, true);

    const std::array<double, 6> potential{1.0, 151.0, 190.0, 6.0, 165.0, 195.0};
    AssignPotentialsToWakeCompressiblePerturbationElement(p_element, distances, potential);

    // Compute LHS
    Matrix LHS = ZeroMatrix(6, 6);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateLeftHandSide(LHS, r_current_process_info);

    // Check the LHS values
    std::array<double,36> reference{0.4851881876577658,-0.4851881876577658,0,0,0,0,-0.6125,1.225,-0.6125,0.6125,
    -1.225,0.6125,0,-0.1212970469144414,0.1212970469144414,0,0,0,-0.6125,0.6125,-0,0.6125,-0.6125,0,0,0,0,
    -0.4851881876577658,0.9703763753155316,-0.4851881876577658,0,0,0,0,-0.3638911407433243,0.3638911407433243};

    for (unsigned int i = 0; i < LHS.size1(); i++) {
        for (unsigned int j = 0; j < LHS.size2(); j++) {
            KRATOS_EXPECT_NEAR(LHS(i, j), reference[6 * i + j], 1e-16);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(PingWakeCompressiblePerturbationPotentialFlowElementLHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    const std::array<double, 6> potential{1.0, 40.0, 35.0, 6.0, 26.0, 14.0};

    Matrix LHS_finite_diference = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);

    ComputeWakeElementalSensitivities(model_part, LHS_finite_diference, LHS_analytical, potential);

    KRATOS_EXPECT_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(PingWakeCompressiblePerturbationPotentialFlowElementLHSClamping, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    const std::array<double, 6> potential{1.285837, 170.29384, 135.1583, 6.0, 196.345, 114.0};

    Matrix LHS_finite_diference = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);

    ComputeWakeElementalSensitivities(model_part, LHS_finite_diference, LHS_analytical, potential);

    KRATOS_EXPECT_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(PingWakeStructureCompressiblePerturbationPotentialFlowElementLHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    p_element->Set(STRUCTURE);
    p_element->GetGeometry()[number_of_nodes-1].SetValue(TRAILING_EDGE, true);

    const std::array<double, 6> potential{1.285837, 30.29384, 35.1583, 6.0, 46.345, 64.0};

    Matrix LHS_finite_diference = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);

    ComputeWakeElementalSensitivities(model_part, LHS_finite_diference, LHS_analytical, potential);

    PrintTestElementInfo(model_part);

    KRATOS_EXPECT_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(PingWakeStructureCompressiblePerturbationPotentialFlowElementLHSClamping, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    p_element->Set(STRUCTURE);
    p_element->GetGeometry()[number_of_nodes-1].SetValue(TRAILING_EDGE, true);

    const std::array<double, 6> potential{1.285837, 170.29384, 135.1583, 6.0, 196.345, 114.0};

    Matrix LHS_finite_diference = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);

    ComputeWakeElementalSensitivities(model_part, LHS_finite_diference, LHS_analytical, potential);

    PrintTestElementInfo(model_part);

    KRATOS_EXPECT_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

} // namespace Testing
} // namespace Kratos.
