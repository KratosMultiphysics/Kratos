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
    rModelPart.GetProcessInfo()[MACH_LIMIT] = 1.73205080756887729;

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

void AssignFreeStreamValues(ModelPart& rModelPart) {
    rModelPart.GetProcessInfo()[FREE_STREAM_DENSITY] = 1.225;
    rModelPart.GetProcessInfo()[FREE_STREAM_MACH] = 0.6;
    rModelPart.GetProcessInfo()[HEAT_CAPACITY_RATIO] = 1.4;
    rModelPart.GetProcessInfo()[SOUND_VELOCITY] = 340.0;
    rModelPart.GetProcessInfo()[MACH_LIMIT] = 1.73205080756887729;
    rModelPart.GetProcessInfo()[CRITICAL_MACH] = 0.99;
    rModelPart.GetProcessInfo()[UPWIND_FACTOR_CONSTANT] = 1.0;

    BoundedVector<double, 3> free_stream_velocity = ZeroVector(3);
    free_stream_velocity(0) = rModelPart.GetProcessInfo().GetValue(FREE_STREAM_MACH) *
                              rModelPart.GetProcessInfo().GetValue(SOUND_VELOCITY);
    rModelPart.GetProcessInfo()[FREE_STREAM_VELOCITY] = free_stream_velocity;
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

void AssignCustomPerturbationPotentialsToElement(Element& rElement, const std::array<double, 3> rPotential)
{
    for (unsigned int i = 0; i < 3; i++)
        rElement.GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = rPotential[i];
}

KRATOS_TEST_CASE_IN_SUITE(ComputePerturbedVelocity, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTestingElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    AssignPerturbationPotentialsToElement(*p_element);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    array_1d<double, 2> perturbed_velocity = PotentialFlowUtilities::ComputePerturbedVelocity<2,3>(*p_element, r_current_process_info);

    KRATOS_EXPECT_RELATIVE_NEAR(perturbed_velocity[0], 303.0, 1e-15);
    KRATOS_EXPECT_RELATIVE_NEAR(perturbed_velocity[1], 50.0, 1e-15);
}

// checks the function ComputeVelocityMagnitude from utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeVelocityMagnitude, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    AssignFreeStreamValues(model_part);

    // velocity corresponding to squared mach number of 3.0
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    const double local_velocity_squared = PotentialFlowUtilities::ComputeVelocityMagnitude<2, 3>(3.0, r_current_process_info);

    const double reference_velocity_squared = 232356.0;

    KRATOS_EXPECT_RELATIVE_NEAR(local_velocity_squared, reference_velocity_squared, 1e-15);
}

// Checks the function ComputeMaximumVelocitySquared from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeMaximumVelocitySquared, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    AssignFreeStreamValues(model_part);

    const double reference_max_velocity_squared = 232356.0;

    // Max local Mach number = sqrt(3.0), from MACH_LIMIT
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    const double max_velocity_squared = PotentialFlowUtilities::ComputeMaximumVelocitySquared<2, 3>(r_current_process_info);

    KRATOS_EXPECT_RELATIVE_NEAR(max_velocity_squared, reference_max_velocity_squared, 1e-15);
}

// Checks the function ComputeVacuumVelocitySquared from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeVacuumVelocitySquared, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    AssignFreeStreamValues(model_part);

    const double reference_max_velocity_squared = 619616.0;

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    const double vacuum_velocity_squared = PotentialFlowUtilities::ComputeVacuumVelocitySquared(r_current_process_info);

    KRATOS_EXPECT_RELATIVE_NEAR(vacuum_velocity_squared, reference_max_velocity_squared, 1e-15);
}

// Checks the function ComputeLocalSpeedOfSound from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeLocalSpeedOfSound, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTestingElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    AssignPotentialsToElement(*p_element);
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    const double local_speed_of_sound =
        PotentialFlowUtilities::ComputeLocalSpeedOfSound<2, 3>(
            *p_element, r_current_process_info);
    KRATOS_EXPECT_NEAR(local_speed_of_sound, 333.801138, 1e-6);
}

// Checks ComputeLocalSpeedofSoundSquared in utilities, local velocity that should be clamped
KRATOS_TEST_CASE_IN_SUITE(ComputeLocalSpeedofSoundSquared, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    AssignFreeStreamValues(model_part);

    const double local_mach_number_squared = 3.0;

    // velocity corresponding to mach number 3.0
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    const double local_velocity_squared = PotentialFlowUtilities::ComputeVelocityMagnitude<2, 3>(local_mach_number_squared, r_current_process_info);

    array_1d<double, 2> velocity(2, 0.0);
    velocity[0] = std::sqrt(local_velocity_squared);

    const double local_speed_sound_squared = PotentialFlowUtilities::ComputeLocalSpeedofSoundSquared<2,3>(velocity, r_current_process_info);

    const double reference_local_speed_sound_squared = 77452.0;

    KRATOS_EXPECT_RELATIVE_NEAR(local_speed_sound_squared, reference_local_speed_sound_squared, 1e-15);
}

// Checks the function ComputeLocalMachNumber from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeLocalMachNumber, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTestingElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    AssignPotentialsToElement(*p_element);
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    const double local_mach_number =
        PotentialFlowUtilities::ComputeLocalMachNumber<2, 3>(
            *p_element, r_current_process_info);
    KRATOS_EXPECT_NEAR(local_mach_number, 0.748948914, 1e-6);
}

// Checks ComputeLocalMachNumberSquared in utilities, local velocity that should be clamped
KRATOS_TEST_CASE_IN_SUITE(ComputeLocalMachNumberSquared, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    AssignFreeStreamValues(model_part);

    const double reference_local_mach_squared = 3.0;

    // velocity corresponding to mach number 3.0
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    const double local_velocity_squared = PotentialFlowUtilities::ComputeVelocityMagnitude<2, 3>(reference_local_mach_squared, r_current_process_info);

    array_1d<double, 2> velocity(2, 0.0);
    velocity[0] = std::sqrt(local_velocity_squared);

    // computes mach number with clamping
    const double local_mach_squared = PotentialFlowUtilities::ComputeLocalMachNumberSquared<2, 3>(velocity, r_current_process_info);

    KRATOS_EXPECT_RELATIVE_NEAR(local_mach_squared, reference_local_mach_squared, 1e-15);
}

// Checks the function ComputeDerivativeLocalMachSquaredWRTVelocitySquared from the utilities, transonic local Mach number
KRATOS_TEST_CASE_IN_SUITE(ComputeDerivativeLocalMachSquaredWRTVelocitySquaredTransonicMach, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    AssignFreeStreamValues(model_part);

    const double local_mach_number_squared = 1.0;

    // velocity corresponding to mach number 1.0
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    const double local_velocity_squared = PotentialFlowUtilities::ComputeVelocityMagnitude<2, 3>(local_mach_number_squared, r_current_process_info);

    array_1d<double, 2> velocity(2, 0.0);
    velocity[0] = std::sqrt(local_velocity_squared);

    // performs clamping on mach number
    const double local_mach_squared = PotentialFlowUtilities::ComputeLocalMachNumberSquared<2, 3>(velocity, r_current_process_info);

    const double mach_derivative = PotentialFlowUtilities::ComputeDerivativeLocalMachSquaredWRTVelocitySquared<2, 3>(velocity,
                local_mach_squared, r_current_process_info);

    const double reference_derivative = 1.1620100191086091883e-05;

    KRATOS_EXPECT_RELATIVE_NEAR(mach_derivative, reference_derivative, 1e-16);
}

// Checks the function ComputeDerivativeLocalMachSquaredWRTVelocitySquared from the utilities, supersonic local Mach number
KRATOS_TEST_CASE_IN_SUITE(ComputeDerivativeLocalMachSquaredWRTVelocitySquaredSupersonicMach, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    AssignFreeStreamValues(model_part);

    const double local_mach_number_squared = 3.0;

    // velocity corresponding to mach number 3.0
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    const double local_velocity_squared = PotentialFlowUtilities::ComputeVelocityMagnitude<2, 3>(local_mach_number_squared, r_current_process_info);

    array_1d<double, 2> velocity(2, 0.0);
    velocity[0] = std::sqrt(local_velocity_squared);

    // performs clamping on mach number
    const double local_mach_squared = PotentialFlowUtilities::ComputeLocalMachNumberSquared<2, 3>(velocity, r_current_process_info);

    const double mach_derivative = PotentialFlowUtilities::ComputeDerivativeLocalMachSquaredWRTVelocitySquared<2, 3>(velocity,
                local_mach_squared, r_current_process_info);

    const double reference_derivative = 2.0657955895264170124e-05;

    KRATOS_EXPECT_RELATIVE_NEAR(mach_derivative, reference_derivative, 1e-16);
}

// Checks the function ComputeIncompressiblePerturbationPressureCoefficient from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputePerturbationIncompressiblePressureCoefficient, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTestingElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    AssignPerturbationPotentialsToElement(*p_element);
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    const double pressure_coefficient =
        PotentialFlowUtilities::ComputePerturbationIncompressiblePressureCoefficient<2, 3>(
            *p_element, r_current_process_info);

    KRATOS_EXPECT_NEAR(pressure_coefficient, -1.266171664744329, 1e-15);
}

// Checks the function ComputePerturbationCompressiblePressureCoefficient from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputePerturbationCompressiblePressureCoefficient, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTestingElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    AssignPerturbationPotentialsToElement(*p_element);
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    const double pressure_coefficient =
        PotentialFlowUtilities::ComputePerturbationCompressiblePressureCoefficient<2, 3>(
            *p_element, r_current_process_info);

    KRATOS_EXPECT_NEAR(pressure_coefficient, -1.128385779511008, 1e-15);
}

// Checks the function ComputePerturbationCompressiblePressureCoefficient from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputePerturbationCompressiblePressureCoefficientClamped, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTestingElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    const std::array<double, 3>& potential{1.0, 733.13764, 929.1948};
    AssignCustomPerturbationPotentialsToElement(*p_element, potential);
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    const double pressure_coefficient =
        PotentialFlowUtilities::ComputePerturbationCompressiblePressureCoefficient<2, 3>(
            *p_element, r_current_process_info);

    const double reference_pressure_coefficient = -2.99132886677513;
    const double tolerance = 1e-15;

    KRATOS_ERROR_IF(!(std::abs(pressure_coefficient - reference_pressure_coefficient) < tolerance))
        << "Check failed because pressure_coefficient  = " << pressure_coefficient <<
        " is not near to reference_pressure_coefficient = " << reference_pressure_coefficient <<
        " within the tolerance " << tolerance << std::endl;
}

// Checks the function ComputePerturbationLocalSpeedOfSound from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputePerturbationLocalSpeedOfSound, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTestingElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    AssignPerturbationPotentialsToElement(*p_element);
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    const double local_speed_of_sound =
        PotentialFlowUtilities::ComputePerturbationLocalSpeedOfSound<2, 3>(
            *p_element, r_current_process_info);

    KRATOS_EXPECT_NEAR(local_speed_of_sound, 324.1317633309022, 1e-13);
}

// Checks the function ComputePerturbationLocalMachNumber from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputePerturbationLocalMachNumber, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTestingElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    AssignPerturbationPotentialsToElement(*p_element);
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    const double local_mach_number =
        PotentialFlowUtilities::ComputePerturbationLocalMachNumber<2, 3>(
            *p_element, r_current_process_info);

    KRATOS_EXPECT_NEAR(local_mach_number, 0.9474471158469713, 1e-16);
}

// tests the function ComputeUpwindFactor from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeUpwindFactor, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    AssignFreeStreamValues(model_part);

    const double local_mach_squared = 3.0;

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    const double upwind_factor = PotentialFlowUtilities::ComputeUpwindFactor<2, 3>(local_mach_squared, r_current_process_info);

    KRATOS_EXPECT_RELATIVE_NEAR(upwind_factor, 0.6733, 1e-15);
}

// tests the function SelectMaxUpwindFactor from the utilities
KRATOS_TEST_CASE_IN_SUITE(SelectMaxUpwindFactor, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    AssignFreeStreamValues(model_part);

    const double local_mach_number_squared = 3.0;

    // velocity corresponding to mach number sqrt(3.0)
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    const double local_velocity_squared = PotentialFlowUtilities::ComputeVelocityMagnitude<2, 3>(local_mach_number_squared, r_current_process_info);

    array_1d<double, 2> current_velocity(2, 0.0);
    current_velocity[0] = std::sqrt(local_velocity_squared);

    const double upwind_mach_number_squared = 0.7 * 0.7;

    // velocity corresponding to mach number 0.7
    const double upwind_velocity_squared = PotentialFlowUtilities::ComputeVelocityMagnitude<2, 3>(upwind_mach_number_squared, r_current_process_info);

    array_1d<double, 2> upwind_velocity(2, 0.0);
    upwind_velocity[0] = std::sqrt(upwind_velocity_squared);

    // find max of current element upwind factor, upwind element upwind factor, and 0
    const double upwind_factor = PotentialFlowUtilities::SelectMaxUpwindFactor<2, 3>(current_velocity, upwind_velocity, r_current_process_info);

    KRATOS_EXPECT_RELATIVE_NEAR(upwind_factor, 0.6733, 1e-15);
}

// tests the function ComputeUpwindFactorCase from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeUpwindFactorCase0, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    AssignFreeStreamValues(model_part);

    array_1d<double, 3> upwind_factor_options(3, 0.0);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    upwind_factor_options[1] = PotentialFlowUtilities::ComputeUpwindFactor<2, 3>(0.35, r_current_process_info);
    upwind_factor_options[2] = PotentialFlowUtilities::ComputeUpwindFactor<2, 3>(0.49, r_current_process_info);

    const auto upwind_factor_case = PotentialFlowUtilities::ComputeUpwindFactorCase<2,3>(upwind_factor_options);

    KRATOS_EXPECT_RELATIVE_NEAR(upwind_factor_case, 0.0, 1e-15);
}

// tests the function ComputeUpwindFactorCase from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeUpwindFactorCase1, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    AssignFreeStreamValues(model_part);

    array_1d<double, 3> upwind_factor_options(3, 0.0);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    upwind_factor_options[1] = PotentialFlowUtilities::ComputeUpwindFactor<2, 3>(3.0, r_current_process_info);
    upwind_factor_options[2] = PotentialFlowUtilities::ComputeUpwindFactor<2, 3>(0.49, r_current_process_info);

    const auto upwind_factor_case = PotentialFlowUtilities::ComputeUpwindFactorCase<2,3>(upwind_factor_options);

    KRATOS_EXPECT_RELATIVE_NEAR(upwind_factor_case, 1.0, 1e-15);
}

// tests the function ComputeUpwindFactorCase from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeUpwindFactorCase2, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    AssignFreeStreamValues(model_part);

    array_1d<double, 3> upwind_factor_options(3, 0.0);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    upwind_factor_options[1] = PotentialFlowUtilities::ComputeUpwindFactor<2, 3>(1.3, r_current_process_info);
    upwind_factor_options[2] = PotentialFlowUtilities::ComputeUpwindFactor<2, 3>(3.0, r_current_process_info);

    const auto upwind_factor_case = PotentialFlowUtilities::ComputeUpwindFactorCase<2,3>(upwind_factor_options);

    KRATOS_EXPECT_RELATIVE_NEAR(upwind_factor_case, 2.0, 1e-15);
}

// tests the function ComputeUpwindFactorCase from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeUpwindFactorCaseSubsonicElement, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    AssignFreeStreamValues(model_part);

    array_1d<double, 3> upwind_factor_options(3, 0.0);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    upwind_factor_options[1] = PotentialFlowUtilities::ComputeUpwindFactor<2, 3>(0.49, r_current_process_info);
    upwind_factor_options[2] = PotentialFlowUtilities::ComputeUpwindFactor<2, 3>(3.0, r_current_process_info);

    const auto upwind_factor_case = PotentialFlowUtilities::ComputeUpwindFactorCase<2,3>(upwind_factor_options);

    KRATOS_EXPECT_RELATIVE_NEAR(upwind_factor_case, 0.0, 1e-15);
}

// tests the function ComputeUpwindFactorDerivativeWRTMachSquared from utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeUpwindFactorDerivativeWRTMachSquared, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    AssignFreeStreamValues(model_part);

    const double local_mach_number_squared = 3.0;

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    const double upwind_factor_derivative = PotentialFlowUtilities::ComputeUpwindFactorDerivativeWRTMachSquared<2,3>(local_mach_number_squared,r_current_process_info);

    KRATOS_EXPECT_RELATIVE_NEAR(upwind_factor_derivative, 0.1089, 1e-15);
}

// tests the function ComputeUpwindFactorDerivativeWRTVelocitySquared from utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeUpwindFactorDerivativeWRTVelocitySquared, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    AssignFreeStreamValues(model_part);

    const double local_mach_number_squared = 3.0;

    // velocity corresponding to mach number sqrt(3.0)
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    const double local_velocity_squared = PotentialFlowUtilities::ComputeVelocityMagnitude<2, 3>(local_mach_number_squared, r_current_process_info);

    array_1d<double, 2> current_velocity(2, 0.0);
    current_velocity[0] = std::sqrt(local_velocity_squared);

    const double upwind_factor_derivative = PotentialFlowUtilities::ComputeUpwindFactorDerivativeWRTVelocitySquared<2,3>(current_velocity, r_current_process_info);

    const double mach_derivative_ref = 2.0657955895264163348e-05;
    const double upwind_factor_ref = 0.1089;

    KRATOS_EXPECT_RELATIVE_NEAR(upwind_factor_derivative, mach_derivative_ref * upwind_factor_ref, 1e-15);
}

// tests the function ComputeDensity from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeDensity, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    AssignFreeStreamValues(model_part);

    const double local_mach_number_squared = 3.0;

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    const double density = PotentialFlowUtilities::ComputeDensity<2, 3>(local_mach_number_squared, r_current_process_info);

    KRATOS_EXPECT_RELATIVE_NEAR(density, 0.450114595263459, 1e-15);
}

// tests the function ComputeDensity from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeUpwindedDensity, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    AssignFreeStreamValues(model_part);

    const double local_mach_number_squared = 3.0;

    // velocity corresponding to mach number sqrt(3.0)
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    const double local_velocity_squared = PotentialFlowUtilities::ComputeVelocityMagnitude<2, 3>(local_mach_number_squared, r_current_process_info);

    array_1d<double, 2> current_velocity(2, 0.0);
    current_velocity[0] = std::sqrt(local_velocity_squared);

    const double upwind_mach_number_squared = 0.7 * 0.7;

    // velocity corresponding to mach number 0.7
    const double upwind_velocity_squared = PotentialFlowUtilities::ComputeVelocityMagnitude<2, 3>(upwind_mach_number_squared, r_current_process_info);

    array_1d<double, 2> upwind_velocity(2, 0.0);
    upwind_velocity[0] = std::sqrt(upwind_velocity_squared);

    const double upwinded_density = PotentialFlowUtilities::ComputeUpwindedDensity<2,3>(current_velocity, upwind_velocity, r_current_process_info);

    KRATOS_EXPECT_RELATIVE_NEAR(upwinded_density, 0.92388212928098, 1e-15);
}

// tests the function ComputeDensityDerivativeWRTVelocitySquared from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeDensityDerivativeWRTVelocitySquared, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    AssignFreeStreamValues(model_part);

    const double local_mach_number_squared = 3.0;

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    const double density_derivative = PotentialFlowUtilities::ComputeDensityDerivativeWRTVelocitySquared<2,3>(local_mach_number_squared, r_current_process_info);

    KRATOS_EXPECT_RELATIVE_NEAR(density_derivative, -2.905764830239754E-06, 1e-15);
}

// tests the function ComputeUpwindedDensityDerivativeWRTVelocitySquaredSupersonicAccelerating from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeUpwindedDensityDerivativeWRTVelocitySquaredSupersonicAccelerating, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    AssignFreeStreamValues(model_part);

    const double local_mach_number_squared = 3.0;

    // velocity corresponding to mach number sqrt(3.0)
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    const double local_velocity_squared = PotentialFlowUtilities::ComputeVelocityMagnitude<2, 3>(local_mach_number_squared, r_current_process_info);

    array_1d<double, 2> current_velocity(2, 0.0);
    current_velocity[0] = std::sqrt(local_velocity_squared);

    const double upwind_mach_number_squared = 0.7 * 0.7;

    const double density_derivative_accel = PotentialFlowUtilities::ComputeUpwindedDensityDerivativeWRTVelocitySquaredSupersonicAccelerating<2,3>(
        current_velocity, local_mach_number_squared, upwind_mach_number_squared, model_part.GetProcessInfo());

    const double reference = 6.33653798760679499e-07;

    KRATOS_EXPECT_RELATIVE_NEAR(density_derivative_accel, reference, 1e-13);
}

// tests the function ComputeUpwindedDensityDerivativeWRTVelocitySquaredSupersonicDeaccelerating from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeUpwindedDensityDerivativeWRTVelocitySquaredSupersonicDeaccelerating, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    AssignFreeStreamValues(model_part);

    const double local_mach_number_squared = 1.3;

    const double upwind_mach_number_squared = 3.0;

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    const double density_derivative_deaccel = PotentialFlowUtilities::ComputeUpwindedDensityDerivativeWRTVelocitySquaredSupersonicDeaccelerating<2,3>(
        local_mach_number_squared, upwind_mach_number_squared, r_current_process_info);

    KRATOS_EXPECT_RELATIVE_NEAR(density_derivative_deaccel, -1.3584190201217436E-06, 1e-15);
}

// tests the function ComputeUpwindedDensityDerivativeWRTUpwindVelocitySquaredSupersonicAccelerating from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeUpwindedDensityDerivativeWRTUpwindVelocitySquaredSupersonicAccelerating, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    AssignFreeStreamValues(model_part);

    const double local_mach_number_squared = 3.0;

    const double upwind_mach_number_squared = 0.7 * 0.7;

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    const double density_derivative_accel = PotentialFlowUtilities::ComputeUpwindedDensityDerivativeWRTUpwindVelocitySquaredSupersonicAccelerating<2,3>(
        local_mach_number_squared, upwind_mach_number_squared, r_current_process_info);

    KRATOS_EXPECT_RELATIVE_NEAR(density_derivative_accel, -3.441482308103857E-06, 1e-15);
}

// tests the function ComputeUpwindedDensityDerivativeWRTUpwindVelocitySquaredSupersonicDeaccelerating from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeUpwindedDensityDerivativeWRTUpwindVelocitySquaredSupersonicDeaccelerating, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    AssignFreeStreamValues(model_part);

    const double local_mach_number_squared = 1.3;

    const double upwind_mach_number_squared = 3.0;

    // velocity corresponding to mach number sqrt(3.0)
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    const double upwind_velocity_squared = PotentialFlowUtilities::ComputeVelocityMagnitude<2, 3>(upwind_mach_number_squared, r_current_process_info);

    array_1d<double, 2> upwind_velocity(2, 0.0);
    upwind_velocity[0] = std::sqrt(upwind_velocity_squared);

    const double density_derivative_deaccel = PotentialFlowUtilities::ComputeUpwindedDensityDerivativeWRTUpwindVelocitySquaredSupersonicDeaccelerating<2,3>(
        upwind_velocity, local_mach_number_squared, upwind_mach_number_squared, r_current_process_info);

    KRATOS_EXPECT_RELATIVE_NEAR(density_derivative_deaccel, -2.7838255012122669E-06, 1e-15);
}

} // namespace Testing
} // namespace Kratos.
