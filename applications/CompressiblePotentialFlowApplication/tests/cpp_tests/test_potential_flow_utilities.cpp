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
    rModelPart.GetProcessInfo()[MACH_LIMIT] = 0.99;
    rModelPart.GetProcessInfo()[UPWINDING_FACTOR_CONSTANT] = 1.0;

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

void GenerateTestingUpstreamElement(ModelPart& rModelPart) {
    // Create extra node
    rModelPart.CreateNewNode(4, 0.0, 1.0, 0.0);
    // Nodes Ids
    std::vector<ModelPart::IndexType> upstream_elemNodes{1, 3, 4};
    Properties::Pointer pElemProp = rModelPart.pGetProperties(0);
    rModelPart.CreateNewElement("CompressiblePotentialFlowElement2D3N", 2, upstream_elemNodes, pElemProp);
}

void AssignPotentialsToElement(Element& rElement) {
    // Define the nodal values
    std::array<double, 3> potential{0.0, 150.0, 350.0};

    for (unsigned int i = 0; i < 3; i++)
        rElement.GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) =
            potential[i];
}

void AssignPerturbationPotentialsToElement(Element& rElement, const std::array<double, 3> rPotential) {
    for (unsigned int i = 0; i < 3; i++){
        rElement.GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = rPotential[i];
    }
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

    std::array<double, 3> potential{1.0, 100.0, 150.0};
    AssignPerturbationPotentialsToElement(*pElement, potential);

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

    std::array<double, 3> potential{1.0, 100.0, 150.0};
    AssignPerturbationPotentialsToElement(*pElement, potential);

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

    std::array<double, 3> potential{1.0, 100.0, 150.0};
    AssignPerturbationPotentialsToElement(*pElement, potential);

    const double local_speed_of_sound =
        PotentialFlowUtilities::ComputePerturbationLocalSpeedOfSound<2, 3>(
            *pElement, model_part.GetProcessInfo());

    KRATOS_CHECK_NEAR(local_speed_of_sound, 324.1317633309022, 1e-13);
}

// Checks the function ComputePerturbationLocalSpeedOfSound from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeMaximumVelocitySquared, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTestingElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);

    std::array<double, 3> potential{1.0, 100.0, 150.0};
    AssignPerturbationPotentialsToElement(*pElement, potential);

    std::cout.precision(16);
    const double maximum_velocity_squared =
        PotentialFlowUtilities::ComputeMaximumVelocitySquared<2, 3>(
            *pElement, model_part.GetProcessInfo());

    KRATOS_CHECK_NEAR(maximum_velocity_squared, 232356.00000000003, 1e-13);
}

// Checks the function ComputePerturbationLocalMachNumber from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputePerturbationLocalMachNumber, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTestingElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);

    std::array<double, 3> potential{1.0, 100.0, 150.0};
    AssignPerturbationPotentialsToElement(*pElement, potential);

    const double local_mach_number =
        PotentialFlowUtilities::ComputePerturbationLocalMachNumber<2, 3>(
            *pElement, model_part.GetProcessInfo());

    KRATOS_CHECK_NEAR(local_mach_number, 0.9474471158469713, 1e-16);
}

// Checks the function ComputeUpwindFactor from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeUpwindFactor, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTestingElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);

    std::array<double, 3> potential{1.0, 100.0, 150.0};
    AssignPerturbationPotentialsToElement(*pElement, potential);

    const double upwind_factor = PotentialFlowUtilities::ComputeUpwindFactor<2, 3>(
        *pElement, model_part.GetProcessInfo());

    KRATOS_CHECK_NEAR(upwind_factor, -0.09184360071679243, 1e-16);
}

// Checks the function ComputeSwitchingOperatorSubsonic from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeSwitchingOperatorSubsonic, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTestingElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);
    std::array<double, 3> potential{1.0, 100.0, 150.0};
    AssignPerturbationPotentialsToElement(*pElement, potential);

    GenerateTestingUpstreamElement(model_part);
    Element::Pointer pUpstreamElement = model_part.pGetElement(2);
    std::array<double, 3> upstream_potential{1.0, 150.0, 51.0};
    AssignPerturbationPotentialsToElement(*pUpstreamElement, upstream_potential);

    const double upwind_factor = PotentialFlowUtilities::ComputeSwitchingOperator<2, 3>(
        *pElement, *pUpstreamElement, model_part.GetProcessInfo());

    KRATOS_CHECK_NEAR(upwind_factor, 0.0, 1e-16);
}

// Checks the function ComputeSwitchingOperator from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeSwitchingOperatorSupersonicAccelerating, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTestingElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);
    std::array<double, 3> potential{1.0, 120.0, 180.0};
    AssignPerturbationPotentialsToElement(*pElement, potential);

    GenerateTestingUpstreamElement(model_part);
    Element::Pointer pUpstreamElement = model_part.pGetElement(2);
    std::array<double, 3> upstream_potential{1.0, 180.0, 90.0};
    AssignPerturbationPotentialsToElement(*pUpstreamElement, upstream_potential);

    const double upwind_factor = PotentialFlowUtilities::ComputeSwitchingOperator<2, 3>(
        *pElement, *pUpstreamElement, model_part.GetProcessInfo());

    KRATOS_CHECK_NEAR(upwind_factor, 0.07067715127537522, 1e-16);
}

// Checks the function ComputeSwitchingOperator from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeSwitchingOperatorSupersonicDecelerating, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTestingElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);
    std::array<double, 3> potential{1.0, 120.0, 180.0};
    AssignPerturbationPotentialsToElement(*pElement, potential);

    GenerateTestingUpstreamElement(model_part);
    Element::Pointer pUpstreamElement = model_part.pGetElement(2);
    std::array<double, 3> upstream_potential{1.0, 180.0, 51.0};
    AssignPerturbationPotentialsToElement(*pUpstreamElement, upstream_potential);

    const double upwind_factor = PotentialFlowUtilities::ComputeSwitchingOperator<2, 3>(
        *pElement, *pUpstreamElement, model_part.GetProcessInfo());

    KRATOS_CHECK_NEAR(upwind_factor, 0.1248655818465636, 1e-16);
}

// Checks the function ComputePerturbationDensity from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputePerturbationDensity, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTestingElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);
    std::array<double, 3> potential{1.0, 100.0, 150.0};
    AssignPerturbationPotentialsToElement(*pElement, potential);

    const double density = PotentialFlowUtilities::ComputePerturbationDensity<2, 3>(
        *pElement, model_part.GetProcessInfo());

    KRATOS_CHECK_NEAR(density, 0.9646048979014981, 1e-16);
}

// Checks the function ComputeUpwindDensity from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeUpwindDensitySubsonic, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTestingElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);
    std::array<double, 3> potential{1.0, 100.0, 150.0};
    AssignPerturbationPotentialsToElement(*pElement, potential);

    GenerateTestingUpstreamElement(model_part);
    Element::Pointer pUpstreamElement = model_part.pGetElement(2);
    std::array<double, 3> upstream_potential{1.0, 150.0, 51.0};
    AssignPerturbationPotentialsToElement(*pUpstreamElement, upstream_potential);

    const double upwind_density = PotentialFlowUtilities::ComputeUpwindDensity<2, 3>(
        *pElement, *pUpstreamElement, model_part.GetProcessInfo());

    KRATOS_CHECK_NEAR(upwind_density, 0.9646048979014981, 1e-16);
}

// Checks the function ComputeUpwindDensity from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeUpwindDensitySupersonicAccelerating, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTestingElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);
    std::array<double, 3> potential{1.0, 120.0, 180.0};
    AssignPerturbationPotentialsToElement(*pElement, potential);

    GenerateTestingUpstreamElement(model_part);
    Element::Pointer pUpstreamElement = model_part.pGetElement(2);
    std::array<double, 3> upstream_potential{1.0, 180.0, 90.0};
    AssignPerturbationPotentialsToElement(*pUpstreamElement, upstream_potential);

    const double upwind_density = PotentialFlowUtilities::ComputeUpwindDensity<2, 3>(
        *pElement, *pUpstreamElement, model_part.GetProcessInfo());

    KRATOS_CHECK_NEAR(upwind_density, 0.9076084702345625, 1e-16);
}

// Checks the function ComputeUpwindDensity from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeUpwindDensitySupersonicDecelerating, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTestingElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);
    std::array<double, 3> potential{1.0, 120.0, 180.0};
    AssignPerturbationPotentialsToElement(*pElement, potential);

    GenerateTestingUpstreamElement(model_part);
    Element::Pointer pUpstreamElement = model_part.pGetElement(2);
    std::array<double, 3> upstream_potential{1.0, 180.0, 51.0};
    AssignPerturbationPotentialsToElement(*pUpstreamElement, upstream_potential);

    const double upwind_density = PotentialFlowUtilities::ComputeUpwindDensity<2, 3>(
        *pElement, *pUpstreamElement, model_part.GetProcessInfo());

    KRATOS_CHECK_NEAR(upwind_density, 0.9003057182942017, 1e-16);
}

// Checks the function ComputeDerivativeUpwindFactorWRTSquareMachNumber from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeDerivativeUpwindFactorWRTSquareMachNumber, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTestingElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);
    std::array<double, 3> potential{1.0, 120.0, 180.0};
    AssignPerturbationPotentialsToElement(*pElement, potential);

    const double DUpwindFactor_DMachNumberSquared =
        PotentialFlowUtilities::ComputeDerivativeUpwindFactorWRTMachNumberSquared<2, 3>(
            *pElement, model_part.GetProcessInfo());

    KRATOS_CHECK_NEAR(DUpwindFactor_DMachNumberSquared, 0.8811763668622098, 1e-16);
}

// Checks the function ComputeDerivativeMachNumberSquaredWRTVelocitySquared from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeDerivativeMachNumberSquaredWRTVelocitySquared, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTestingElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);
    std::array<double, 3> potential{1.0, 120.0, 180.0};
    AssignPerturbationPotentialsToElement(*pElement, potential);

    const array_1d<double, 3> free_stream_velocity = model_part.GetProcessInfo()[FREE_STREAM_VELOCITY];
    array_1d<double, 2> velocity = PotentialFlowUtilities::ComputeVelocity<2,3>(*pElement);
    for (unsigned int i = 0; i < 2; i++){
        velocity[i] += free_stream_velocity[i];
    }

    const double DMachNumberSquared_DVelocitySquared =
        PotentialFlowUtilities::ComputeDerivativeMachNumberSquaredWRTVelocitySquared<2, 3>(
            *pElement, model_part.GetProcessInfo());

    KRATOS_CHECK_NEAR(DMachNumberSquared_DVelocitySquared, 1.1832700207409171 * 1e-5, 1e-16);
}

// Checks the function ComputeUpwindDensity from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeDrhoDphiSupersonicAccelerating, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTestingElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);
    std::array<double, 3> potential{1.0, 120.0, 180.0};
    AssignPerturbationPotentialsToElement(*pElement, potential);

    GenerateTestingUpstreamElement(model_part);
    Element::Pointer pUpstreamElement = model_part.pGetElement(2);
    std::array<double, 3> upstream_potential{1.0, 180.0, 90.0};
    AssignPerturbationPotentialsToElement(*pUpstreamElement, upstream_potential);

    const BoundedVector<double, 3> DrhoDphi = PotentialFlowUtilities::ComputeDrhoDphiSupersonicAccelerating<2, 3>(
        *pElement, *pUpstreamElement, model_part.GetProcessInfo());

    std::vector<double> reference{0.002237981723350806,-0.001822257564214433,-0.0004157241591363726};

    KRATOS_CHECK_VECTOR_NEAR(DrhoDphi, reference, 1e-16);
}

// Checks the function ComputeUpwindDensity from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeDrhoDphiUpSupersonicAccelerating, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTestingElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);
    std::array<double, 3> potential{1.0, 120.0, 180.0};
    AssignPerturbationPotentialsToElement(*pElement, potential);

    GenerateTestingUpstreamElement(model_part);
    Element::Pointer pUpstreamElement = model_part.pGetElement(2);
    std::array<double, 3> upstream_potential{1.0, 180.0, 90.0};
    AssignPerturbationPotentialsToElement(*pUpstreamElement, upstream_potential);

    const BoundedVector<double, 3> DrhoDphiUp = PotentialFlowUtilities::ComputeDrhoDphiUpSupersonicAccelerating<2, 3>(
        *pElement, *pUpstreamElement, model_part.GetProcessInfo());

    std::vector<double> reference{5.774518723453986e-05, -0.0001907537645725249, 0.0001330085773379851};

    KRATOS_CHECK_VECTOR_NEAR(DrhoDphiUp, reference, 1e-16);
}

// Checks the function ComputeUpwindDensity from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeDrhoDphiSupersonicDecelerating, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTestingElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);
    std::array<double, 3> potential{1.0, 120.0, 180.0};
    AssignPerturbationPotentialsToElement(*pElement, potential);

    GenerateTestingUpstreamElement(model_part);
    Element::Pointer pUpstreamElement = model_part.pGetElement(2);
    std::array<double, 3> upstream_potential{1.0, 180.0, 51.0};
    AssignPerturbationPotentialsToElement(*pUpstreamElement, upstream_potential);

    const BoundedVector<double, 3> DrhoDphi = PotentialFlowUtilities::ComputeDrhoDphiSupersonicDecelerating<2, 3>(
        *pElement, *pUpstreamElement, model_part.GetProcessInfo());

    std::vector<double> reference{0.002494998902346807,-0.002031531613985171,-0.0004634672883616359};

    KRATOS_CHECK_VECTOR_NEAR(DrhoDphi, reference, 1e-16);
}

// Checks the function ComputeUpwindDensity from the utilities
KRATOS_TEST_CASE_IN_SUITE(ComputeDrhoDphiUpSupersonicDecelerating, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTestingElement(model_part);
    Element::Pointer pElement = model_part.pGetElement(1);
    std::array<double, 3> potential{1.0, 120.0, 180.0};
    AssignPerturbationPotentialsToElement(*pElement, potential);

    GenerateTestingUpstreamElement(model_part);
    Element::Pointer pUpstreamElement = model_part.pGetElement(2);
    std::array<double, 3> upstream_potential{1.0, 180.0, 51.0};
    AssignPerturbationPotentialsToElement(*pUpstreamElement, upstream_potential);

    const BoundedVector<double, 3> DrhoDphiUp = PotentialFlowUtilities::ComputeDrhoDphiUpSupersonicDecelerating<2, 3>(
        *pElement, *pUpstreamElement, model_part.GetProcessInfo());

    std::vector<double> reference{7.680874010502327e-05,-0.0005115462090994551,0.0004347374689944318};

    KRATOS_CHECK_VECTOR_NEAR(DrhoDphiUp, reference, 1e-16);
}

} // namespace Testing
} // namespace Kratos.
