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
    rModelPart.GetProcessInfo()[MACH_LIMIT] = 0.99;

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

void AssignPotentialsToNormalCompressiblePerturbationElement(Element::Pointer pElement)
{
    std::array<double, 3> potential({1.0, 100.0, 150.0});

    for (unsigned int i = 0; i < 3; i++)
        pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = potential[i];
}

/** Checks the CompressiblePerturbationPotentialFlowElement.
 * Checks the RHS computation.
 */
KRATOS_TEST_CASE_IN_SUITE(CompressiblePerturbationPotentialFlowElementRHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);

    AssignPotentialsToNormalCompressiblePerturbationElement(pElement);

    // Compute RHS
    Vector RHS = ZeroVector(3);

    pElement->CalculateRightHandSide(RHS, model_part.GetProcessInfo());

    std::vector<double> reference({146.264326, -122.142628, -24.1216977});

    KRATOS_CHECK_VECTOR_NEAR(RHS, reference, 1e-6);
}

/** Checks the CompressiblePerturbationPotentialFlowElement.
 * Checks the LHS computation.
 */
KRATOS_TEST_CASE_IN_SUITE(CompressiblePerturbationPotentialFlowElementLHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);

    AssignPotentialsToNormalCompressiblePerturbationElement(pElement);

    // Compute LHS
    Matrix LHS = ZeroMatrix(3, 3);

    pElement->CalculateLeftHandSide(LHS, model_part.GetProcessInfo());

    std::array<double, 9> reference({0.0611427846, -0.130621505, 0.0694787204,
                                     -0.130621505,0.671075851,-0.540454346,
                                     0.0694787204,-0.540454346,0.470975625});

    for (unsigned int i = 0; i < LHS.size1(); i++) {
        for (unsigned int j = 0; j < LHS.size2(); j++) {
            KRATOS_CHECK_NEAR(LHS(i, j), reference[i * 3 + j], 1e-6);
        }
    }
}

BoundedVector<double,3> AssignDistancesToPerturbationCompressibleElement()
{
    BoundedVector<double,3> distances;
    distances(0) = 1.0;
    distances(1) = -1.0;
    distances(2) = -1.0;
    return distances;
}

void AssignPotentialsToWakeCompressiblePerturbationElement(Element::Pointer pElement, const array_1d<double, 3>& rDistances)
{
    const std::array<double, 3> potential({1.0, 100.0, 150.0});

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

KRATOS_TEST_CASE_IN_SUITE(WakeCompressiblePerturbationPotentialFlowElementRHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);

    BoundedVector<double,3> distances = AssignDistancesToPerturbationCompressibleElement();

    pElement->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    pElement->GetValue(WAKE) = true;

    AssignPotentialsToWakeCompressiblePerturbationElement(pElement, distances);

    // Compute RHS and LHS
    Vector RHS = ZeroVector(6);

    pElement->CalculateRightHandSide(RHS, model_part.GetProcessInfo());

    std::vector<double> reference({146.264326, 0, 0, 0, -122.142628, -24.1216977});

    KRATOS_CHECK_VECTOR_NEAR(RHS, reference, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(WakeCompressiblePerturbationPotentialFlowElementLHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);

    BoundedVector<double,3> distances = AssignDistancesToPerturbationCompressibleElement();

    pElement->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    pElement->GetValue(WAKE) = true;

    AssignPotentialsToWakeCompressiblePerturbationElement(pElement, distances);

    // Compute LHS
    Matrix LHS = ZeroMatrix(6, 6);

    pElement->CalculateLeftHandSide(LHS, model_part.GetProcessInfo());
    std::cout.precision(9);
    KRATOS_WATCH(LHS)

    // Check the LHS values
    std::array<double,36> reference({0.0611427846,-0.130621505,0.0694787204,0,0,0,
                                  -0.130621505,0.671075851,-0.540454346,0.130621505,-0.671075851,0.540454346,
                                  0.0694787204,-0.540454346,0.470975625,-0.0694787204,0.540454346,-0.470975625,
                                  -0.0611427846,0.130621505,-0.0694787204,0.0611427846,-0.130621505,0.0694787204,
                                  0,0,0,-0.130621505,0.671075851,-0.540454346,
                                  0,0,0,0.0694787204,-0.540454346,0.470975625});

    for (unsigned int i = 0; i < LHS.size1(); i++) {
        for (unsigned int j = 0; j < LHS.size2(); j++) {
            KRATOS_CHECK_NEAR(LHS(i, j), reference[6 * i + j], 1e-6);
        }
    }
}

} // namespace Testing
} // namespace Kratos.
