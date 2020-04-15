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
#include "custom_utilities/potential_flow_utilities.h"
#include "compressible_potential_flow_application_variables.h"
#include "fluid_dynamics_application_variables.h"
#include "custom_elements/compressible_potential_flow_element.h"

namespace Kratos {
namespace Testing {

typedef ModelPart::IndexType IndexType;
typedef ModelPart::NodeIterator NodeIteratorType;

void GenerateTestingElement(ModelPart& rModelPart) {
    // Variables addition
    rModelPart.AddNodalSolutionStepVariable(VELOCITY_POTENTIAL);
    rModelPart.AddNodalSolutionStepVariable(AUXILIARY_VELOCITY_POTENTIAL);

    // Set the element properties
    Properties::Pointer pElemProp = rModelPart.CreateNewProperties(0);

    rModelPart.GetProcessInfo()[FREE_STREAM_DENSITY] = 1.225;
    rModelPart.GetProcessInfo()[FREE_STREAM_MACH] = 0.6;
    rModelPart.GetProcessInfo()[HEAT_CAPACITY_RATIO] = 1.4;
    rModelPart.GetProcessInfo()[SOUND_VELOCITY] = 340.0;

    BoundedVector<double, 3> free_stream_velocity = ZeroVector(3);
    free_stream_velocity(0) = rModelPart.GetProcessInfo().GetValue(FREE_STREAM_MACH) *
                              rModelPart.GetProcessInfo().GetValue(SOUND_VELOCITY);
    rModelPart.GetProcessInfo()[FREE_STREAM_VELOCITY] = free_stream_velocity;

    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> elemNodes{1, 2, 3};
    rModelPart.CreateNewElement("CompressiblePotentialFlowElement2D3N", 1, elemNodes, pElemProp);
}

void AssignPotentialsToElement(Element& rElement) {
    // Define the nodal values
    std::array<double, 3> potential{0.0, 150.0, 350.0};

    for (unsigned int i = 0; i < 3; i++)
        rElement.GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) =
            potential[i];
}

void AssignPerturbationPotentialsToElement(Element& rElement) {
    // Define the nodal values
    std::array<double, 3> potential{1.0, 100.0, 150.0};

    for (unsigned int i = 0; i < 3; i++)
        rElement.GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) =
            potential[i];
}

// Checks the function ComputeLocalSpeedOfSound from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeLocalSpeedOfSound, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTestingElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);

    AssignPotentialsToElement(*pElement);
    const double local_speed_of_sound =
        PotentialFlowUtilities::ComputeLocalSpeedOfSound<2, 3>(
            *pElement, model_part.GetProcessInfo());
    KRATOS_CHECK_NEAR(local_speed_of_sound, 333.801138, 1e-6);
}

// Checks the function ComputeLocalMachNumber from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeLocalMachNumber, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTestingElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);

    AssignPotentialsToElement(*pElement);
    const double local_mach_number =
        PotentialFlowUtilities::ComputeLocalMachNumber<2, 3>(
            *pElement, model_part.GetProcessInfo());
    KRATOS_CHECK_NEAR(local_mach_number, 0.748948914, 1e-6);
}

// Checks the function ComputeIncompressiblePerturbationPressureCoefficient from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputePerturbationIncompressiblePressureCoefficient, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTestingElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);

    AssignPerturbationPotentialsToElement(*pElement);
    const double pressure_coefficient =
        PotentialFlowUtilities::ComputePerturbationIncompressiblePressureCoefficient<2, 3>(
            *pElement, model_part.GetProcessInfo());

    KRATOS_CHECK_NEAR(pressure_coefficient, -1.266171664744329, 1e-15);
}

// Checks the function ComputePerturbationCompressiblePressureCoefficient from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputePerturbationCompressiblePressureCoefficient, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTestingElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);

    AssignPerturbationPotentialsToElement(*pElement);
    const double pressure_coefficient =
        PotentialFlowUtilities::ComputePerturbationCompressiblePressureCoefficient<2, 3>(
            *pElement, model_part.GetProcessInfo());

    KRATOS_CHECK_NEAR(pressure_coefficient, -1.128385779511008, 1e-15);
}

// Checks the function ComputePerturbationLocalSpeedOfSound from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputePerturbationLocalSpeedOfSound, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTestingElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);

    AssignPerturbationPotentialsToElement(*pElement);
    const double local_speed_of_sound =
        PotentialFlowUtilities::ComputePerturbationLocalSpeedOfSound<2, 3>(
            *pElement, model_part.GetProcessInfo());

    KRATOS_CHECK_NEAR(local_speed_of_sound, 324.1317633309022, 1e-13);
}

// Checks the function ComputePerturbationLocalMachNumber from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputePerturbationLocalMachNumber, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTestingElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);

    AssignPerturbationPotentialsToElement(*pElement);
    const double local_mach_number =
        PotentialFlowUtilities::ComputePerturbationLocalMachNumber<2, 3>(
            *pElement, model_part.GetProcessInfo());

    KRATOS_CHECK_NEAR(local_mach_number, 0.9474471158469713, 1e-16);
}

// Checks the function ComputeLocalMachSquaredDerivative from the utilities, transonic local Mach number
KRATOS_TEST_CASE_IN_SUITE(ComputeLocalMachSquaredDerivativeTransonicMach, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    model_part.GetProcessInfo()[FREE_STREAM_DENSITY] = 1.225;
    model_part.GetProcessInfo()[FREE_STREAM_MACH] = 0.6;
    model_part.GetProcessInfo()[HEAT_CAPACITY_RATIO] = 1.4;
    model_part.GetProcessInfo()[SOUND_VELOCITY] = 340.0;

    BoundedVector<double, 3> free_stream_velocity = ZeroVector(3);
    free_stream_velocity(0) = model_part.GetProcessInfo().GetValue(FREE_STREAM_MACH) *
                              model_part.GetProcessInfo().GetValue(SOUND_VELOCITY);
    model_part.GetProcessInfo()[FREE_STREAM_VELOCITY] = free_stream_velocity;

    array_1d<double, 2> velocity(2, 0.0);
    velocity[0] = 68.0 * sqrt(67.0/3.0);

    double sq_local_mach_number = PotentialFlowUtilities::ComputeLocalMachNumberSquared<2, 3>(velocity, model_part.GetProcessInfo());
    double local_mach_number = sqrt(sq_local_mach_number);

    auto mach_derivative = PotentialFlowUtilities::ComputeLocalMachSquaredDerivative<2, 3>(velocity,
                local_mach_number, model_part.GetProcessInfo());

    double reference_derivative = 1.16201001910861e-05;

    KRATOS_CHECK_NEAR(mach_derivative, reference_derivative, 1e-16);
}

// Checks the function ComputeLocalMachSquaredDerivative from the utilities, supersonic local Mach number
KRATOS_TEST_CASE_IN_SUITE(ComputeLocalMachSquaredDerivativeSupersonicMach, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    model_part.GetProcessInfo()[FREE_STREAM_DENSITY] = 1.225;
    model_part.GetProcessInfo()[FREE_STREAM_MACH] = 0.6;
    model_part.GetProcessInfo()[HEAT_CAPACITY_RATIO] = 1.4;
    model_part.GetProcessInfo()[SOUND_VELOCITY] = 340.0;

    BoundedVector<double, 3> free_stream_velocity = ZeroVector(3);
    free_stream_velocity(0) = model_part.GetProcessInfo().GetValue(FREE_STREAM_MACH) *
                              model_part.GetProcessInfo().GetValue(SOUND_VELOCITY);
    model_part.GetProcessInfo()[FREE_STREAM_VELOCITY] = free_stream_velocity;

    array_1d<double, 2> velocity(2, 0.0);
    velocity[0] = 272.0 * sqrt(134.0/21.0);

    double sq_local_mach_number = PotentialFlowUtilities::ComputeLocalMachNumberSquared<2, 3>(velocity, model_part.GetProcessInfo());
    double local_mach_number = sqrt(sq_local_mach_number);

    double mach_derivative = PotentialFlowUtilities::ComputeLocalMachSquaredDerivative<2, 3>(velocity,
                local_mach_number, model_part.GetProcessInfo());
                
    double reference_derivative = 2.06579558952642e-05;

    KRATOS_CHECK_NEAR(mach_derivative, reference_derivative, 1e-16);
}

// Checks the function ComputeMaximumVelocitySquared from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeMaximumVelocitySquared, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    model_part.GetProcessInfo()[FREE_STREAM_DENSITY] = 1.225;
    model_part.GetProcessInfo()[FREE_STREAM_MACH] = 0.6;
    model_part.GetProcessInfo()[HEAT_CAPACITY_RATIO] = 1.4;
    model_part.GetProcessInfo()[SOUND_VELOCITY] = 340.0;

    BoundedVector<double, 3> free_stream_velocity = ZeroVector(3);
    free_stream_velocity(0) = model_part.GetProcessInfo().GetValue(FREE_STREAM_MACH) *
                              model_part.GetProcessInfo().GetValue(SOUND_VELOCITY);
    model_part.GetProcessInfo()[FREE_STREAM_VELOCITY] = free_stream_velocity;

    // Max local Mach number = sqrt(3.0), hard coded in ComputeMaximumVelocitySquared

    double reference_sq_max_velocity = 232356.0000000000000000;

    double sq_max_velocity = PotentialFlowUtilities::ComputeMaximumVelocitySquared<2, 3>(model_part.GetProcessInfo());

    KRATOS_CHECK_NEAR(sq_max_velocity, reference_sq_max_velocity, 1e-10);
}

// Checks ComputeLocalSpeedofSoundSquared in utilities, local velocity that should be clamped
KRATOS_TEST_CASE_IN_SUITE(ComputeLocalSpeedofSoundSquared, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    model_part.GetProcessInfo()[FREE_STREAM_DENSITY] = 1.225;
    model_part.GetProcessInfo()[FREE_STREAM_MACH] = 0.6;
    model_part.GetProcessInfo()[HEAT_CAPACITY_RATIO] = 1.4;
    model_part.GetProcessInfo()[SOUND_VELOCITY] = 340.0;

    BoundedVector<double, 3> free_stream_velocity = ZeroVector(3);
    free_stream_velocity(0) = model_part.GetProcessInfo().GetValue(FREE_STREAM_MACH) *
                              model_part.GetProcessInfo().GetValue(SOUND_VELOCITY);
    model_part.GetProcessInfo()[FREE_STREAM_VELOCITY] = free_stream_velocity;

    array_1d<double, 2> velocity(2, 0.0);
    velocity[0] = 350000.0;
    
    double sq_local_speed_sound = PotentialFlowUtilities::ComputeLocalSpeedofSoundSquared<2,3>(velocity, model_part.GetProcessInfo());

    double reference_sq_local_speed_sound = 77452.0000000000;

    KRATOS_CHECK_NEAR(sq_local_speed_sound, reference_sq_local_speed_sound, 1e-16);
}

// Checks ComputeLocalMachNumberSquared in utilities, local velocity that should be clamped
KRATOS_TEST_CASE_IN_SUITE(ComputeLocalMachNumberSquared, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    model_part.GetProcessInfo()[FREE_STREAM_DENSITY] = 1.225;
    model_part.GetProcessInfo()[FREE_STREAM_MACH] = 0.6;
    model_part.GetProcessInfo()[HEAT_CAPACITY_RATIO] = 1.4;
    model_part.GetProcessInfo()[SOUND_VELOCITY] = 340.0;

    BoundedVector<double, 3> free_stream_velocity = ZeroVector(3);
    free_stream_velocity(0) = model_part.GetProcessInfo().GetValue(FREE_STREAM_MACH) *
                              model_part.GetProcessInfo().GetValue(SOUND_VELOCITY);
    model_part.GetProcessInfo()[FREE_STREAM_VELOCITY] = free_stream_velocity;

    array_1d<double, 2> velocity(2, 0.0);
    velocity[0] = 350000.0;

    double sq_local_mach_number = PotentialFlowUtilities::ComputeLocalMachNumberSquared<2, 3>(velocity, model_part.GetProcessInfo());

    double reference_sq_local_mach_number = 3.0;

    KRATOS_CHECK_NEAR(sq_local_mach_number, reference_sq_local_mach_number, 1e-10);
}

} // namespace Testing
} // namespace Kratos.
