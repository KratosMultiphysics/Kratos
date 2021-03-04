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
#include "custom_elements/compressible_perturbation_potential_flow_element.h"
#include "custom_utilities/potential_flow_utilities.h"

namespace Kratos {
namespace Testing {

typedef ModelPart::IndexType IndexType;
typedef ModelPart::NodeIterator NodeIteratorType;

void GenerateCompressiblePerturbationElement3D(ModelPart& rModelPart) {
    // Variables addition
    rModelPart.AddNodalSolutionStepVariable(VELOCITY_POTENTIAL);
    rModelPart.AddNodalSolutionStepVariable(AUXILIARY_VELOCITY_POTENTIAL);

    // Set the element properties
    Properties::Pointer pElemProp = rModelPart.CreateNewProperties(0);
    rModelPart.GetProcessInfo()[FREE_STREAM_DENSITY] = 1.225;
    rModelPart.GetProcessInfo()[FREE_STREAM_MACH] = 0.6;
    rModelPart.GetProcessInfo()[HEAT_CAPACITY_RATIO] = 1.4;
    rModelPart.GetProcessInfo()[SOUND_VELOCITY] = 340.3;
    // rModelPart.GetProcessInfo()[MACH_LIMIT] = 0.99;
    rModelPart.GetProcessInfo()[MACH_SQUARED_LIMIT] = 0.8836;

    BoundedVector<double, 3> free_stream_velocity = ZeroVector(3);
    free_stream_velocity(0) = rModelPart.GetProcessInfo().GetValue(FREE_STREAM_MACH) * rModelPart.GetProcessInfo().GetValue(SOUND_VELOCITY);
    rModelPart.GetProcessInfo()[FREE_STREAM_VELOCITY] = free_stream_velocity;

    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, -0.2, -0.2);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 0.1, 1.0, 0.0);
    rModelPart.CreateNewNode(4, -0.1, 0.0, 1.0);
    std::vector<ModelPart::IndexType> elemNodes{1, 2, 3, 4};
    rModelPart.CreateNewElement("CompressiblePerturbationPotentialFlowElement3D4N", 1, elemNodes, pElemProp);
}

void AssignPotentialsToNormalCompressiblePerturbationElement(Element::Pointer pElement)
{
    std::array<double, 3> potential{1.0, 100.0, 150.0};

    for (unsigned int i = 0; i < 3; i++)
        pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = potential[i];
}

void AssignPerturbationPotentialsToNormalCompressibleElement3D(Element& rElement, const std::array<double, 4> rPotential) {
    for (unsigned int i = 0; i < 4; i++){
        rElement.GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = rPotential[i];
    }
}

// /** Checks the CompressiblePerturbationPotentialFlowElement.
//  * Checks the RHS computation.
//  */
// KRATOS_TEST_CASE_IN_SUITE(CompressiblePerturbationPotentialFlowElementRHS3D, CompressiblePotentialApplicationFastSuite) {
//     Model this_model;
//     ModelPart& model_part = this_model.CreateModelPart("Main", 3);

//     GenerateCompressiblePerturbationElement3D(model_part);
//     Element::Pointer pElement = model_part.pGetElement(1);

//     std::array<double, 4> potential{1.39572, 143.39275, 151.1549827, 134.284736};
//     AssignPerturbationPotentialsToNormalCompressibleElement3D(*pElement, potential);

//     // Compute RHS
//     Vector RHS = ZeroVector(4);

//     const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
//     pElement->CalculateRightHandSide(RHS, r_current_process_info);

//     std::vector<double> reference{72.53100745113115,-39.92932281025485,-17.27378025102953,-15.32790438984678};

//     KRATOS_CHECK_VECTOR_NEAR(RHS, reference, 1e-13);
// }

// /** Checks the CompressiblePerturbationPotentialFlowElement.
//  * Checks the LHS computation.
//  */
// KRATOS_TEST_CASE_IN_SUITE(CompressiblePerturbationPotentialFlowElementLHS3D, CompressiblePotentialApplicationFastSuite) {
//     Model this_model;
//     ModelPart& model_part = this_model.CreateModelPart("Main", 3);

//     GenerateCompressiblePerturbationElement3D(model_part);
//     Element::Pointer pElement = model_part.pGetElement(1);

//     std::array<double, 4> potential{1.39572, 143.39275, 151.1549827, 134.284736};
//     AssignPerturbationPotentialsToNormalCompressibleElement3D(*pElement, potential);

//     // Compute LHS
//     Matrix LHS = ZeroMatrix(4, 4);

//     const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
//     pElement->CalculateLeftHandSide(LHS, r_current_process_info);

//     std::array<double, 16> reference{ -0.1488791034771964,0.1571111809330265,0.002522732543031606,-0.01075480999886159,
//                                        0.1571111809330265,-0.03464650206168052,-0.06488707484410745,-0.05757760402723854,
//                                        0.002522732543031606,-0.06488707484410745,0.08727292600217176,-0.02490858370109592,
//                                       -0.01075480999886159,-0.05757760402723854,-0.02490858370109592,0.09324099772719605};

//     std::cout.precision(16);

//     KRATOS_WATCH(LHS)

//     for (unsigned int i = 0; i < LHS.size1(); i++) {
//         for (unsigned int j = 0; j < LHS.size2(); j++) {
//             KRATOS_CHECK_NEAR(LHS(i, j), reference[i * 4 + j], 1e-16);
//         }
//     }
// }

/** Checks the CompressiblePerturbationPotentialFlowElement.
 * Tests the LHS computation.
 */
KRATOS_TEST_CASE_IN_SUITE(PingCompressiblePerturbationPotentialFlowElementLHS3D, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement3D(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);
    const unsigned int number_of_nodes = pElement->GetGeometry().size();

    std::array<double, 4> potential{1.39572, 110.69275, 221.1549827, 304.284736};
    AssignPerturbationPotentialsToNormalCompressibleElement3D(*pElement, potential);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    // array_1d<double, 3> perturbed_velocity = PotentialFlowUtilities::ComputePerturbedVelocity<3,4>(*pElement, r_current_process_info);
    // const double local_velocity_squared = inner_prod(perturbed_velocity, perturbed_velocity);
    // const double local_mach_squared = PotentialFlowUtilities::ComputeLocalMachNumberSquared<3, 4>(perturbed_velocity, r_current_process_info);
    // const double max_velocity_squared = PotentialFlowUtilities::ComputeMaximumVelocitySquared<3, 4>(r_current_process_info);

    Vector RHS_original = ZeroVector(number_of_nodes);
    Matrix LHS_original = ZeroMatrix(number_of_nodes, number_of_nodes);
    Matrix LHS_finite_diference = ZeroMatrix(number_of_nodes, number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(number_of_nodes, number_of_nodes);

    // Compute original RHS and LHS
    //const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
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

    // std::cout.precision(16);

    // KRATOS_WATCH(perturbed_velocity)
    // KRATOS_WATCH(std::sqrt(max_velocity_squared))
    // KRATOS_WATCH(std::sqrt(local_velocity_squared))
    // KRATOS_WATCH(local_mach_squared)
    // KRATOS_WATCH(LHS_analytical)
    // KRATOS_WATCH(LHS_finite_diference)

    for (unsigned int i = 0; i < LHS_finite_diference.size1(); i++) {
        for (unsigned int j = 0; j < LHS_finite_diference.size2(); j++) {
            KRATOS_CHECK_NEAR(LHS_finite_diference(i,j), LHS_analytical(i,j), 1e-10);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(PingCompressiblePerturbationPotentialFlowElementLHS3DClamping, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement3D(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);
    const unsigned int number_of_nodes = pElement->GetGeometry().size();

    std::array<double, 4> potential{1.39572, 117.69275, 221.1549827, 304.284736};
    AssignPerturbationPotentialsToNormalCompressibleElement3D(*pElement, potential);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();

    Vector RHS_original = ZeroVector(number_of_nodes);
    Matrix LHS_original = ZeroMatrix(number_of_nodes, number_of_nodes);
    Matrix LHS_finite_diference = ZeroMatrix(number_of_nodes, number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(number_of_nodes, number_of_nodes);

    // Compute original RHS and LHS
    //const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
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

BoundedVector<double,4> AssignDistancesToPerturbationCompressibleElement3D()
{
    BoundedVector<double,4> distances;
    distances(0) = -1.0;
    distances(1) = -1.0;
    distances(2) = -1.0;
    distances(3) = 1.0;
    return distances;
}

void AssignPotentialsToWakeCompressiblePerturbationElement3D(Element::Pointer pElement, const array_1d<double, 4>& rDistances,  const std::array<double, 8> rPotential)
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

// KRATOS_TEST_CASE_IN_SUITE(WakeCompressiblePerturbationPotentialFlowElementRHS, CompressiblePotentialApplicationFastSuite) {
//     Model this_model;
//     ModelPart& model_part = this_model.CreateModelPart("Main", 3);

//     GenerateCompressiblePerturbationElement3D(model_part);
//     Element::Pointer pElement = model_part.pGetElement(1);

//     BoundedVector<double,4> distances = AssignDistancesToPerturbationCompressibleElement3D();

//     pElement->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
//     pElement->GetValue(WAKE) = true;

//     std::array<double, 8> potential{1.39572, 66.69275, 101.1549827, 104.284736, 2.39572, 46.69275, 100.1549827, 102.284736};
//     AssignPotentialsToWakeCompressiblePerturbationElement3D(pElement, distances, potential);

//     // Compute RHS and LHS
//     Vector RHS = ZeroVector(4);

//     const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
//     pElement->CalculateRightHandSide(RHS, r_current_process_info);

//     // std::cout.precision(16);
//     // KRATOS_WATCH(RHS)

//     std::vector<double> reference{4.141666666666665,-5.798333333333331,0.8283333333333329,-5.68648586161461,65.67375589700852,-56.08141746367978,-4.117673984175311,-0.8283333333333329};

//     KRATOS_CHECK_VECTOR_NEAR(RHS, reference, 1e-13);
// }

KRATOS_TEST_CASE_IN_SUITE(PingWakeCompressiblePerturbationPotentialFlowElementLHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement3D(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);
    const unsigned int number_of_nodes = pElement->GetGeometry().size();

    BoundedVector<double,4> distances = AssignDistancesToPerturbationCompressibleElement3D();

    pElement->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    pElement->GetValue(WAKE) = true;

    std::array<double, 8> potential{1.39572, 110.69275, 121.1549827, 104.284736, 2.39572, 46.69275, 100.1549827, 102.284736};
    AssignPotentialsToWakeCompressiblePerturbationElement3D(pElement, distances, potential);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    // array_1d<double, 3> perturbed_velocity = PotentialFlowUtilities::ComputePerturbedVelocity<3,4>(*pElement, r_current_process_info);
    // const double local_velocity_squared = inner_prod(perturbed_velocity, perturbed_velocity);
    // const double local_mach_squared = PotentialFlowUtilities::ComputeLocalMachNumberSquared<3, 4>(perturbed_velocity, r_current_process_info);
    // const double max_velocity_squared = PotentialFlowUtilities::ComputeMaximumVelocitySquared<3, 4>(r_current_process_info);

    Vector RHS_original = ZeroVector(2*number_of_nodes);
    Matrix LHS_original = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);
    Matrix LHS_finite_diference = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);

    // Compute original RHS and LHS
    //const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    pElement->CalculateLocalSystem(LHS_original, RHS_original, r_current_process_info);
    // KRATOS_WATCH(RHS_original)

    double delta = 1e-3;
    for(unsigned int i = 0; i < 2*number_of_nodes; i++){
        if(i < number_of_nodes){
            // Pinging
            if (distances(i) > 0.0)
                pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += delta;
            else
                pElement->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) += delta;

            Vector RHS_pinged = ZeroVector(number_of_nodes);
            Matrix LHS_pinged = ZeroMatrix(number_of_nodes, number_of_nodes);
            // Compute pinged LHS and RHS
            pElement->CalculateLocalSystem(LHS_pinged, RHS_pinged, r_current_process_info);

            for(unsigned int k = 0; k < 2*number_of_nodes; k++){
                // Compute the finite difference estimate of the sensitivity
                LHS_finite_diference( k, i) = -(RHS_pinged(k)-RHS_original(k)) / delta;
                // Compute the average of the original and pinged analytic sensitivities
                LHS_analytical( k, i) = 0.5 * (LHS_original(k,i) + LHS_pinged(k,i));
            }

            // Unpinging
            if (distances(i) > 0.0)
                pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= delta;
            else
                pElement->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) -= delta;
        }
        else{
            // Pinging
            if (distances(i-4) > 0.0)
                pElement->GetGeometry()[i-4].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) += delta;
            else
                pElement->GetGeometry()[i-4].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += delta;

            Vector RHS_pinged = ZeroVector(number_of_nodes);
            Matrix LHS_pinged = ZeroMatrix(number_of_nodes, number_of_nodes);
            // Compute pinged LHS and RHS
            pElement->CalculateLocalSystem(LHS_pinged, RHS_pinged, r_current_process_info);

            for(unsigned int k = 0; k < 2*number_of_nodes; k++){
                // Compute the finite difference estimate of the sensitivity
                LHS_finite_diference( k, i) = -(RHS_pinged(k)-RHS_original(k)) / delta;
                // Compute the average of the original and pinged analytic sensitivities
                LHS_analytical( k, i) = 0.5 * (LHS_original(k,i) + LHS_pinged(k,i));
            }

            // Unpinging
            if (distances(i-4) > 0.0)
                pElement->GetGeometry()[i-4].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) -= delta;
            else
                pElement->GetGeometry()[i-4].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= delta;
        }
    }

    // std::cout.precision(16);

    // KRATOS_WATCH(perturbed_velocity)
    // KRATOS_WATCH(std::sqrt(max_velocity_squared))
    // KRATOS_WATCH(std::sqrt(local_velocity_squared))
    // KRATOS_WATCH(local_mach_squared)

    // std::cout.precision(5);
    // std::cout << std::scientific;
    // std::cout << std::showpos;
    // std::cout << std::endl;
    // for(unsigned int row = 0; row < LHS_analytical.size1(); ++row){
    //     for(unsigned int column = 0; column < LHS_analytical.size2(); column++){
    //         if(column == 3 || column == 7){
    //             std::cout << " " << LHS_analytical(row, column) << " |";
    //         }
    //         else{
    //             std::cout << " " << LHS_analytical(row, column) << " ";
    //         }
    //     }

    //     std::cout << std::endl;

    //     if(row ==3|| row == 7){
    //         for(unsigned int j = 0; j < 14*LHS_analytical.size1(); j++){
    //         std::cout << "_" ;
    //         }
    //         std::cout << " " << std::endl;
    //     }
    //     else{
    //         for(unsigned int i = 0; i < 3; i++){
    //             for(unsigned int j = 0; j < 14*4; j++){
    //                 std::cout << " " ;
    //             }
    //             if(i != 2){
    //                 std::cout << "|" ;
    //             }
    //         }
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;

    // for(unsigned int row = 0; row < LHS_finite_diference.size1(); ++row){
    //     for(unsigned int column = 0; column < LHS_finite_diference.size2(); column++){
    //         if(column == 3 || column == 7){
    //             std::cout << " " << LHS_finite_diference(row, column) << " |";
    //         }
    //         else{
    //             std::cout << " " << LHS_finite_diference(row, column) << " ";
    //         }
    //     }

    //     std::cout << std::endl;

    //     if(row ==3|| row == 7){
    //         for(unsigned int j = 0; j < 14*LHS_finite_diference.size1(); j++){
    //         std::cout << "_" ;
    //         }
    //         std::cout << " " << std::endl;
    //     }
    //     else{
    //         for(unsigned int i = 0; i < 3; i++){
    //             for(unsigned int j = 0; j < 14*4; j++){
    //                 std::cout << " " ;
    //             }
    //             if(i != 2){
    //                 std::cout << "|" ;
    //             }
    //         }
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;


    for (unsigned int i = 0; i < LHS_finite_diference.size1(); i++) {
        for (unsigned int j = 0; j < LHS_finite_diference.size2(); j++) {
            if(std::abs(LHS_finite_diference(i,j)-LHS_analytical(i,j)) > 1e-10){
                KRATOS_WATCH(i)
                KRATOS_WATCH(j)
                KRATOS_WATCH(std::abs(LHS_finite_diference(i,j)-LHS_analytical(i,j)))
            }
            KRATOS_CHECK_NEAR(LHS_finite_diference(i,j), LHS_analytical(i,j), 1e-10);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(PingWakeCompressiblePerturbationPotentialFlowElementLHSClamping, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement3D(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);
    const unsigned int number_of_nodes = pElement->GetGeometry().size();

    BoundedVector<double,4> distances = AssignDistancesToPerturbationCompressibleElement3D();

    pElement->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    pElement->GetValue(WAKE) = true;

    std::array<double, 8> potential{1.39572, 117.69275, 121.1549827, 104.284736, 2.39572, 46.69275, 100.1549827, 102.284736};
    AssignPotentialsToWakeCompressiblePerturbationElement3D(pElement, distances, potential);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();

    Vector RHS_original = ZeroVector(2*number_of_nodes);
    Matrix LHS_original = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);
    Matrix LHS_finite_diference = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);

    // Compute original RHS and LHS
    //const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    pElement->CalculateLocalSystem(LHS_original, RHS_original, r_current_process_info);
    // KRATOS_WATCH(RHS_original)

    double delta = 1e-3;
    for(unsigned int i = 0; i < 2*number_of_nodes; i++){
        if(i < number_of_nodes){
            // Pinging
            if (distances(i) > 0.0)
                pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += delta;
            else
                pElement->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) += delta;

            Vector RHS_pinged = ZeroVector(number_of_nodes);
            Matrix LHS_pinged = ZeroMatrix(number_of_nodes, number_of_nodes);
            // Compute pinged LHS and RHS
            pElement->CalculateLocalSystem(LHS_pinged, RHS_pinged, r_current_process_info);

            for(unsigned int k = 0; k < 2*number_of_nodes; k++){
                // Compute the finite difference estimate of the sensitivity
                LHS_finite_diference( k, i) = -(RHS_pinged(k)-RHS_original(k)) / delta;
                // Compute the average of the original and pinged analytic sensitivities
                LHS_analytical( k, i) = 0.5 * (LHS_original(k,i) + LHS_pinged(k,i));
            }

            // Unpinging
            if (distances(i) > 0.0)
                pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= delta;
            else
                pElement->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) -= delta;
        }
        else{
            // Pinging
            if (distances(i-4) > 0.0)
                pElement->GetGeometry()[i-4].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) += delta;
            else
                pElement->GetGeometry()[i-4].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += delta;

            Vector RHS_pinged = ZeroVector(number_of_nodes);
            Matrix LHS_pinged = ZeroMatrix(number_of_nodes, number_of_nodes);
            // Compute pinged LHS and RHS
            pElement->CalculateLocalSystem(LHS_pinged, RHS_pinged, r_current_process_info);

            for(unsigned int k = 0; k < 2*number_of_nodes; k++){
                // Compute the finite difference estimate of the sensitivity
                LHS_finite_diference( k, i) = -(RHS_pinged(k)-RHS_original(k)) / delta;
                // Compute the average of the original and pinged analytic sensitivities
                LHS_analytical( k, i) = 0.5 * (LHS_original(k,i) + LHS_pinged(k,i));
            }

            // Unpinging
            if (distances(i-4) > 0.0)
                pElement->GetGeometry()[i-4].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) -= delta;
            else
                pElement->GetGeometry()[i-4].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= delta;
        }
    }

    for (unsigned int i = 0; i < LHS_finite_diference.size1(); i++) {
        for (unsigned int j = 0; j < LHS_finite_diference.size2(); j++) {
            KRATOS_CHECK_NEAR(LHS_finite_diference(i,j), LHS_analytical(i,j), 1e-10);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(PingWakeStructureCompressiblePerturbationPotentialFlowElementLHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement3D(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);
    const unsigned int number_of_nodes = pElement->GetGeometry().size();

    BoundedVector<double,4> distances = AssignDistancesToPerturbationCompressibleElement3D();

    pElement->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    pElement->GetValue(WAKE) = true;
    pElement->Set(STRUCTURE);
    pElement->GetGeometry()[number_of_nodes-1].SetValue(TRAILING_EDGE, true);

    std::array<double, 8> potential{1.39572, 110.69275, 121.1549827, 104.284736, 2.39572, 46.69275, 100.1549827, 102.284736};
    AssignPotentialsToWakeCompressiblePerturbationElement3D(pElement, distances, potential);

    Vector RHS_original = ZeroVector(2*number_of_nodes);
    Matrix LHS_original = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);
    Matrix LHS_finite_diference = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);

    // Compute original RHS and LHS
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    pElement->CalculateLocalSystem(LHS_original, RHS_original, r_current_process_info);
    // KRATOS_WATCH(RHS_original)

    double delta = 1e-3;
    for(unsigned int i = 0; i < 2*number_of_nodes; i++){
        if(i < number_of_nodes){
            // Pinging
            if (distances(i) > 0.0)
                pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += delta;
            else
                pElement->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) += delta;

            Vector RHS_pinged = ZeroVector(number_of_nodes);
            Matrix LHS_pinged = ZeroMatrix(number_of_nodes, number_of_nodes);
            // Compute pinged LHS and RHS
            pElement->CalculateLocalSystem(LHS_pinged, RHS_pinged, r_current_process_info);

            for(unsigned int k = 0; k < 2*number_of_nodes; k++){
                // Compute the finite difference estimate of the sensitivity
                LHS_finite_diference( k, i) = -(RHS_pinged(k)-RHS_original(k)) / delta;
                // Compute the average of the original and pinged analytic sensitivities
                LHS_analytical( k, i) = 0.5 * (LHS_original(k,i) + LHS_pinged(k,i));
            }

            // Unpinging
            if (distances(i) > 0.0)
                pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= delta;
            else
                pElement->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) -= delta;
        }
        else{
            // Pinging
            if (distances(i-4) > 0.0)
                pElement->GetGeometry()[i-4].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) += delta;
            else
                pElement->GetGeometry()[i-4].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += delta;

            Vector RHS_pinged = ZeroVector(number_of_nodes);
            Matrix LHS_pinged = ZeroMatrix(number_of_nodes, number_of_nodes);
            // Compute pinged LHS and RHS
            pElement->CalculateLocalSystem(LHS_pinged, RHS_pinged, r_current_process_info);

            for(unsigned int k = 0; k < 2*number_of_nodes; k++){
                // Compute the finite difference estimate of the sensitivity
                LHS_finite_diference( k, i) = -(RHS_pinged(k)-RHS_original(k)) / delta;
                // Compute the average of the original and pinged analytic sensitivities
                LHS_analytical( k, i) = 0.5 * (LHS_original(k,i) + LHS_pinged(k,i));
            }

            // Unpinging
            if (distances(i-4) > 0.0)
                pElement->GetGeometry()[i-4].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) -= delta;
            else
                pElement->GetGeometry()[i-4].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= delta;
        }
    }

    for (unsigned int i = 0; i < LHS_finite_diference.size1(); i++) {
        for (unsigned int j = 0; j < LHS_finite_diference.size2(); j++) {
            KRATOS_CHECK_NEAR(LHS_finite_diference(i,j), LHS_analytical(i,j), 1e-10);
        }
    }
}

void ComputeWakeElementalSensitivities(ModelPart& rModelPart, Matrix& rLHS_finite_diference, Matrix& rLHS_analytical, const std::array<double, 8> rPotential){
    Element::Pointer pElement = rModelPart.pGetElement(1);
    const unsigned int number_of_nodes = pElement->GetGeometry().size();

    BoundedVector<double,4> distances = AssignDistancesToPerturbationCompressibleElement3D();
    pElement->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    pElement->GetValue(WAKE) = true;

    AssignPotentialsToWakeCompressiblePerturbationElement3D(pElement, distances, rPotential);

    // Compute original RHS and LHS
    Vector RHS_original = ZeroVector(2*number_of_nodes);
    Matrix LHS_original = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);
    const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
    pElement->CalculateLocalSystem(LHS_original, RHS_original, r_current_process_info);

    double delta = 1e-3;
    for(unsigned int i = 0; i < 2*number_of_nodes; i++){
        if(i < number_of_nodes){
            // Pinging
            if (distances(i) > 0.0)
                pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += delta;
            else
                pElement->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) += delta;

            Vector RHS_pinged = ZeroVector(number_of_nodes);
            Matrix LHS_pinged = ZeroMatrix(number_of_nodes, number_of_nodes);
            // Compute pinged LHS and RHS
            pElement->CalculateLocalSystem(LHS_pinged, RHS_pinged, r_current_process_info);

            for(unsigned int k = 0; k < 2*number_of_nodes; k++){
                // Compute the finite difference estimate of the sensitivity
                rLHS_finite_diference( k, i) = -(RHS_pinged(k)-RHS_original(k)) / delta;
                // Compute the average of the original and pinged analytic sensitivities
                rLHS_analytical( k, i) = 0.5 * (LHS_original(k,i) + LHS_pinged(k,i));
            }

            // Unpinging
            if (distances(i) > 0.0)
                pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= delta;
            else
                pElement->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) -= delta;
        }
        else{
            // Pinging
            if (distances(i-4) > 0.0)
                pElement->GetGeometry()[i-4].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) += delta;
            else
                pElement->GetGeometry()[i-4].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += delta;

            Vector RHS_pinged = ZeroVector(number_of_nodes);
            Matrix LHS_pinged = ZeroMatrix(number_of_nodes, number_of_nodes);
            // Compute pinged LHS and RHS
            pElement->CalculateLocalSystem(LHS_pinged, RHS_pinged, r_current_process_info);

            for(unsigned int k = 0; k < 2*number_of_nodes; k++){
                // Compute the finite difference estimate of the sensitivity
                rLHS_finite_diference( k, i) = -(RHS_pinged(k)-RHS_original(k)) / delta;
                // Compute the average of the original and pinged analytic sensitivities
                rLHS_analytical( k, i) = 0.5 * (LHS_original(k,i) + LHS_pinged(k,i));
            }

            // Unpinging
            if (distances(i-4) > 0.0)
                pElement->GetGeometry()[i-4].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) -= delta;
            else
                pElement->GetGeometry()[i-4].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= delta;
        }
    }

}

KRATOS_TEST_CASE_IN_SUITE(PingWakeStructureCompressiblePerturbationPotentialFlowElementLHSClamping, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement3D(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);
    const unsigned int number_of_nodes = pElement->GetGeometry().size();

    pElement->Set(STRUCTURE);
    pElement->GetGeometry()[number_of_nodes-1].SetValue(TRAILING_EDGE, true);

    std::array<double, 8> potential{1.39572, 117.69275, 121.1549827, 104.284736, 2.39572, 146.69275, 100.1549827, 102.284736};

    Matrix LHS_finite_diference = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);

    ComputeWakeElementalSensitivities(model_part, LHS_finite_diference, LHS_analytical, potential);

    for (unsigned int i = 0; i < LHS_finite_diference.size1(); i++) {
        for (unsigned int j = 0; j < LHS_finite_diference.size2(); j++) {
            KRATOS_CHECK_NEAR(LHS_finite_diference(i,j), LHS_analytical(i,j), 1e-10);
        }
    }
}

} // namespace Testing
} // namespace Kratos.
