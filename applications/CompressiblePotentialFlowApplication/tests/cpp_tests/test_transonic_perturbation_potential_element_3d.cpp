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
    rModelPart.GetProcessInfo()[MACH_LIMIT] = 0.94;
    rModelPart.GetProcessInfo()[MACH_SQUARED_LIMIT] = 3.0;
    rModelPart.GetProcessInfo()[CRITICAL_MACH] = 0.99;
    rModelPart.GetProcessInfo()[UPWIND_FACTOR_CONSTANT] = 1.0;

    BoundedVector<double, 3> free_stream_velocity = ZeroVector(3);
    free_stream_velocity(0) = rModelPart.GetProcessInfo().GetValue(FREE_STREAM_MACH) * rModelPart.GetProcessInfo().GetValue(SOUND_VELOCITY);
    rModelPart.GetProcessInfo()[FREE_STREAM_VELOCITY] = free_stream_velocity;

    const double velocity_squared_limit = PotentialFlowUtilities::ComputeMaximumVelocitySquared<2,3>(rModelPart.GetProcessInfo());
    rModelPart.GetProcessInfo()[VELOCITY_SQUARED_LIMIT] = velocity_squared_limit;

    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, -0.2, -0.2);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 0.1, 1.0, 0.0);
    rModelPart.CreateNewNode(4, -0.1, 0.0, 1.0);
    std::vector<ModelPart::IndexType> elemNodes{1, 2, 3, 4};
    rModelPart.CreateNewElement("TransonicPerturbationPotentialFlowElement3D4N", 1, elemNodes, pElemProp);
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

void PrintTestMatrixPretty3D(Matrix& rMatrix){
    std::cout.precision(5);
    std::cout << std::scientific;
    std::cout << std::showpos;
    std::cout << std::endl;
    for(unsigned int row = 0; row < rMatrix.size1(); ++row){
        for(unsigned int column = 0; column < rMatrix.size2(); column++){
            if(column == rMatrix.size1()/2-1 || column == rMatrix.size1()-1){
                std::cout << " " << rMatrix(row, column) << " |";
            }
            else{
                std::cout << " " << rMatrix(row, column) << " ";
            }
        }

        std::cout << std::endl;

        if(row ==rMatrix.size1()/2-1|| row == rMatrix.size1()-1){
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

void PrintTransonicTestElementInfo3D(ModelPart& rModelPart, array_1d<double, 3>& perturbed_velocity){
    const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

    const double max_velocity_squared = PotentialFlowUtilities::ComputeMaximumVelocitySquared<2, 3>(r_current_process_info);
    const double local_velocity_squared = inner_prod(perturbed_velocity, perturbed_velocity);
    const double local_mach_squared = PotentialFlowUtilities::ComputeLocalMachNumberSquared<2, 3>(perturbed_velocity, r_current_process_info);

    std::cout.precision(16);
    std::string element_info = "\n\n\nElement Info:";
    KRATOS_WATCH(element_info)
    KRATOS_WATCH(perturbed_velocity)
    KRATOS_WATCH(std::sqrt(max_velocity_squared))
    KRATOS_WATCH(std::sqrt(local_velocity_squared))
    KRATOS_WATCH(local_mach_squared)
}

void ComputeElementalSensitivitiesMatrixRowTransonic3D(ModelPart& rModelPart, double delta, unsigned int row, Matrix& rLHS_original, Vector& rRHS_original, Matrix& rLHS_finite_diference, Matrix& rLHS_analytical){
    Element::Pointer pElement = rModelPart.pGetElement(1);
    const unsigned int number_of_nodes = pElement->GetGeometry().size();

    // Compute pinged LHS and RHS
    Vector RHS_pinged = ZeroVector(number_of_nodes);
    Matrix LHS_pinged = ZeroMatrix(number_of_nodes, number_of_nodes);
    const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
    pElement->CalculateLocalSystem(LHS_pinged, RHS_pinged, r_current_process_info);

    for(unsigned int k = 0; k < rLHS_original.size2(); k++){
        // Compute the finite difference estimate of the sensitivity
        rLHS_finite_diference( k, row) = -(RHS_pinged(k)-rRHS_original(k)) / delta;
        // Compute the average of the original and pinged analytic sensitivities
        rLHS_analytical( k, row) = 0.5 * (rLHS_original(k,row) + LHS_pinged(k,row));
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

        ComputeElementalSensitivitiesMatrixRowTransonic3D(model_part, delta, i, LHS_original, RHS_original, rLHS_finite_diference, rLHS_analytical);

        // Unpinging
        if (i < number_of_nodes) {
            pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= delta;
        }
        else {
            pUpwindElement->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= delta;
        }
    }


    // array_1d<double, 3> perturbed_velocity = PotentialFlowUtilities::ComputePerturbedVelocity<3,4>(*pElement, r_current_process_info);
    // PrintTransonicTestElementInfo3D(model_part, perturbed_velocity);
    // array_1d<double, 3> upwind_perturbed_velocity = PotentialFlowUtilities::ComputePerturbedVelocity<3,4>(*pUpwindElement, r_current_process_info);
    // PrintTransonicTestElementInfo3D(model_part, upwind_perturbed_velocity);
}

/** Checks the TransonicPerturbationPotentialFlowElement.
 * Tests the LHS computation.
 */
KRATOS_TEST_CASE_IN_SUITE(PingTransonicPerturbationPotentialFlowElementLHS3D, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement3D(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);
    const unsigned int number_of_nodes = pElement->GetGeometry().size();

    pElement->SetFlags(INLET);

    std::array<double, 4> potential{1.39572, 110.69275, 221.1549827, 304.284736};
    AssignPerturbationPotentialsToTransonicElement3D(*pElement, potential);

    Vector RHS_original = ZeroVector(number_of_nodes);
    Matrix LHS_original = ZeroMatrix(number_of_nodes, number_of_nodes);
    Matrix LHS_finite_diference = ZeroMatrix(number_of_nodes, number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(number_of_nodes, number_of_nodes);

    // Compute original RHS and LHS
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    pElement->CalculateLocalSystem(LHS_original, RHS_original, r_current_process_info);

    double delta = 1e-3;
    for(unsigned int i = 0; i < number_of_nodes; i++){
        // Pinging
        pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += delta;

        Vector RHS_pinged = ZeroVector(number_of_nodes);
        Matrix LHS_pinged = ZeroMatrix(number_of_nodes, number_of_nodes);
        // Compute pinged LHS and RHS
        pElement->CalculateLocalSystem(LHS_pinged, RHS_pinged, r_current_process_info);

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

/** Checks the TransonicPerturbationPotentialFlowElement.
 * Tests the LHS computation.
 */
KRATOS_TEST_CASE_IN_SUITE(PingTransonicPerturbationPotentialFlowSupersonicElementLHS3D, CompressiblePotentialApplicationFastSuite) {
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

    for (unsigned int i = 0; i < LHS_finite_diference.size1(); i++) {
        for (unsigned int j = 0; j < LHS_finite_diference.size2(); j++) {
            KRATOS_CHECK_NEAR(LHS_finite_diference(i,j), LHS_analytical(i,j), 1e-10);
        }
    }
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

    for (unsigned int i = 0; i < LHS_finite_diference.size1(); i++) {
        for (unsigned int j = 0; j < LHS_finite_diference.size2(); j++) {
            KRATOS_CHECK_NEAR(LHS_finite_diference(i,j), LHS_analytical(i,j), 1e-10);
        }
    }
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

    for (unsigned int i = 0; i < LHS_finite_diference.size1(); i++) {
        for (unsigned int j = 0; j < LHS_finite_diference.size2(); j++) {
            KRATOS_CHECK_NEAR(LHS_finite_diference(i,j), LHS_analytical(i,j), 1e-10);
        }
    }
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

    for (unsigned int i = 0; i < LHS_finite_diference.size1(); i++) {
        for (unsigned int j = 0; j < LHS_finite_diference.size2(); j++) {
            KRATOS_CHECK_NEAR(LHS_finite_diference(i,j), LHS_analytical(i,j), 1e-10);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(TransonicPerturbationPotentialFlowElementEquationId3D, CompressiblePotentialApplicationFastSuite) {
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

BoundedVector<double,4> AssignDistancesToPerturbationTransonicElement3D()
{
    BoundedVector<double,4> distances;
    distances(0) = -1.0;
    distances(1) = -1.0;
    distances(2) = -1.0;
    distances(3) = 1.0;
    return distances;
}

void AssignPotentialsToWakeTransonicPerturbationElement3D(Element::Pointer pElement, const array_1d<double, 4>& rDistances, const std::array<double, 8> rPotential)
{
    for (unsigned int i = 0; i < 4; i++){
        if (rDistances(i) > 0.0)
            pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = rPotential[i];
        else
            pElement->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) = rPotential[i];
    }
    for (unsigned int i = 0; i < 4; i++){
        if (rDistances(i) < 0.0)
            pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = rPotential[i+4];
        else
            pElement->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) = rPotential[i+4];
    }
}

void ComputeWakeElementalSensitivities3D(ModelPart& rModelPart, Matrix& rLHS_finite_diference, Matrix& rLHS_analytical, const std::array<double, 8> rPotential){
    Element::Pointer p_element = rModelPart.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    BoundedVector<double,4> distances = AssignDistancesToPerturbationTransonicElement3D();
    p_element->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    p_element->GetValue(WAKE) = true;

    AssignPotentialsToWakeTransonicPerturbationElement3D(p_element, distances, rPotential);

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

            ComputeElementalSensitivitiesMatrixRowTransonic3D(rModelPart, delta, i, LHS_original, RHS_original, rLHS_finite_diference, rLHS_analytical);

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

            ComputeElementalSensitivitiesMatrixRowTransonic3D(rModelPart, delta, i, LHS_original, RHS_original, rLHS_finite_diference, rLHS_analytical);

            // Unpinging
            if (distances(i-number_of_nodes) > 0.0)
                p_element->GetGeometry()[i-number_of_nodes].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) -= delta;
            else
                p_element->GetGeometry()[i-number_of_nodes].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= delta;
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(PingWakeTransonicPerturbationPotentialFlowElementLHS3D, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement3D(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    const std::array<double, 8>& potential{1.39572, 110.69275, 121.1549827, 104.284736, 2.39572, 46.69275, 100.1549827, 102.284736};

    Matrix LHS_finite_diference = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);

    ComputeWakeElementalSensitivities3D(model_part, LHS_finite_diference, LHS_analytical, potential);

    // PrintTestMatrixPretty(LHS_analytical);
    // PrintTestMatrixPretty(LHS_finite_diference);

    for (unsigned int i = 0; i < LHS_finite_diference.size1(); i++) {
        for (unsigned int j = 0; j < LHS_finite_diference.size2(); j++) {
            if(std::abs(LHS_finite_diference(i,j)-LHS_analytical(i,j)) > 1e-10){
                KRATOS_WATCH(i)
                KRATOS_WATCH(j)
                std::cout.precision(16);
                KRATOS_WATCH(std::abs(LHS_finite_diference(i,j)-LHS_analytical(i,j)))
            }
            KRATOS_CHECK_NEAR(LHS_finite_diference(i,j), LHS_analytical(i,j), 1e-10);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(PingWakeTransonicPerturbationPotentialFlowElementLHSClamping3D, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement3D(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    std::array<double, 8> potential{1.39572, 357.69275, 321.1549827, 304.284736, 2.39572, 346.69275, 200.1549827, 302.284736};

    Matrix LHS_finite_diference = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);

    ComputeWakeElementalSensitivities3D(model_part, LHS_finite_diference, LHS_analytical, potential);

    // PrintTestMatrixPretty(LHS_analytical);
    // PrintTestMatrixPretty(LHS_finite_diference);

    for (unsigned int i = 0; i < LHS_finite_diference.size1(); i++) {
        for (unsigned int j = 0; j < LHS_finite_diference.size2(); j++) {
            if(std::abs(LHS_finite_diference(i,j)-LHS_analytical(i,j)) > 1e-10){
                KRATOS_WATCH(i)
                KRATOS_WATCH(j)
                std::cout.precision(16);
                KRATOS_WATCH(std::abs(LHS_finite_diference(i,j)-LHS_analytical(i,j)))
            }
            KRATOS_CHECK_NEAR(LHS_finite_diference(i,j), LHS_analytical(i,j), 1e-10);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(PingWakeStructureTransonicPerturbationPotentialFlowElementLHS3D, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement3D(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    p_element->Set(STRUCTURE);
    p_element->GetGeometry()[number_of_nodes-1].SetValue(TRAILING_EDGE, true);

    const std::array<double, 8>& potential{1.39572, 110.69275, 121.1549827, 104.284736, 2.39572, 46.69275, 100.1549827, 102.284736};

    Matrix LHS_finite_diference = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);

    ComputeWakeElementalSensitivities3D(model_part, LHS_finite_diference, LHS_analytical, potential);

    // PrintTestMatrixPretty(LHS_analytical);
    // PrintTestMatrixPretty(LHS_finite_diference);

    for (unsigned int i = 0; i < LHS_finite_diference.size1(); i++) {
        for (unsigned int j = 0; j < LHS_finite_diference.size2(); j++) {
            if(std::abs(LHS_finite_diference(i,j)-LHS_analytical(i,j)) > 1e-10){
                KRATOS_WATCH(i)
                KRATOS_WATCH(j)
                std::cout.precision(16);
                KRATOS_WATCH(std::abs(LHS_finite_diference(i,j)-LHS_analytical(i,j)))
            }
            KRATOS_CHECK_NEAR(LHS_finite_diference(i,j), LHS_analytical(i,j), 1e-10);
        }
    }
    // const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    // array_1d<double, 3> perturbed_velocity = PotentialFlowUtilities::ComputePerturbedVelocity<3,4>(*p_element, r_current_process_info);
    // PrintTransonicTestElementInfo3D(model_part, perturbed_velocity);
    // array_1d<double, 3> upwind_perturbed_velocity = PotentialFlowUtilities::ComputePerturbedVelocityLowerElement<3,4>(*p_element, r_current_process_info);
    // PrintTransonicTestElementInfo3D(model_part, upwind_perturbed_velocity);
}

KRATOS_TEST_CASE_IN_SUITE(PingWakeStructureTransonicPerturbationPotentialFlowElementLHSClamping3D, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransonicPerturbationElement3D(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    p_element->Set(STRUCTURE);
    p_element->GetGeometry()[number_of_nodes-1].SetValue(TRAILING_EDGE, true);

    std::array<double, 8> potential{1.39572, 357.69275, 321.1549827, 304.284736, 2.39572, 346.69275, 200.1549827, 302.284736};

    Matrix LHS_finite_diference = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);

    ComputeWakeElementalSensitivities3D(model_part, LHS_finite_diference, LHS_analytical, potential);

    // PrintTestMatrixPretty(LHS_analytical);
    // PrintTestMatrixPretty(LHS_finite_diference);

    for (unsigned int i = 0; i < LHS_finite_diference.size1(); i++) {
        for (unsigned int j = 0; j < LHS_finite_diference.size2(); j++) {
            if(std::abs(LHS_finite_diference(i,j)-LHS_analytical(i,j)) > 1e-10){
                KRATOS_WATCH(i)
                KRATOS_WATCH(j)
                std::cout.precision(16);
                KRATOS_WATCH(std::abs(LHS_finite_diference(i,j)-LHS_analytical(i,j)))
            }
            KRATOS_CHECK_NEAR(LHS_finite_diference(i,j), LHS_analytical(i,j), 1e-10);
        }
    }
    // const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    // array_1d<double, 3> perturbed_velocity = PotentialFlowUtilities::ComputePerturbedVelocity<3,4>(*p_element, r_current_process_info);
    // PrintTransonicTestElementInfo3D(model_part, perturbed_velocity);
    // array_1d<double, 3> upwind_perturbed_velocity = PotentialFlowUtilities::ComputePerturbedVelocityLowerElement<3,4>(*p_element, r_current_process_info);
    // PrintTransonicTestElementInfo3D(model_part, upwind_perturbed_velocity);
}

} // namespace Testing
} // namespace Kratos.
