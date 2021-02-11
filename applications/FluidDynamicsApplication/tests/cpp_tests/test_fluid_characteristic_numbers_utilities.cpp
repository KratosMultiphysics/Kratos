//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/cfd_variables.h"
#include "utilities/element_size_calculator.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "custom_utilities/fluid_characteristic_numbers_utilities.h"

namespace Kratos {
namespace Testing  {

namespace {

void TestFluidCharacteristicNumberInitializeModelPart(
    ModelPart& rModelPart,
    const double DeltaTime)
{

    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    Properties::Pointer p_properties = rModelPart.CreateNewProperties(0);

    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 2.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
    rModelPart.CreateNewNode(4, 1.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> element_nodes_1{1, 2, 3};
    std::vector<ModelPart::IndexType> element_nodes_2{3, 2, 4};
    rModelPart.CreateNewElement("Element2D3N", 1, element_nodes_1, p_properties);
    rModelPart.CreateNewElement("Element2D3N", 2, element_nodes_2, p_properties);

    // Set a fake current delta time
    rModelPart.GetProcessInfo().SetValue(DELTA_TIME, DeltaTime);

    // Set nodal data
    for (auto& rNode : rModelPart.Nodes()) {
        rNode.FastGetSolutionStepValue(VELOCITY_X) = rNode.Id() * rNode.X();
        rNode.FastGetSolutionStepValue(VELOCITY_Y) = rNode.Id() * rNode.Y();
    }
}

void TestFluidCharacteristicNumberElementSet(ModelPart& rModelPart)
{

    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    Properties::Pointer p_properties = rModelPart.CreateNewProperties(0);
    p_properties->SetValue(DENSITY, 0.1);
    p_properties->SetValue(CONDUCTIVITY, 2.0);
    p_properties->SetValue(SPECIFIC_HEAT, 1.0);
    p_properties->SetValue(DYNAMIC_VISCOSITY, 1.0);

    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 2.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> element_nodes_1{1, 2, 3};
    rModelPart.CreateNewElement("Element2D3N", 1, element_nodes_1, p_properties);

    // Set nodal data
    for (auto& rNode : rModelPart.Nodes()) {
        rNode.SetValue(ARTIFICIAL_DYNAMIC_VISCOSITY, rNode.Id());
        rNode.SetValue(ARTIFICIAL_CONDUCTIVITY, 2.0 * rNode.Id());
        rNode.FastGetSolutionStepValue(VELOCITY_X) = rNode.Id() * rNode.X();
        rNode.FastGetSolutionStepValue(VELOCITY_Y) = rNode.Id() * rNode.Y();
    }
}

} // namespace internals

KRATOS_TEST_CASE_IN_SUITE(FluidCharacteristicNumbersCalculateLocalCFL, FluidDynamicsApplicationFastSuite)
{
    // Set the current delta time to calculate the CFL number
    const double current_dt = 1.0e-1;

    // Create the test model part
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("TestModelPart");
    TestFluidCharacteristicNumberInitializeModelPart(r_model_part, current_dt);

    // Calculate the CFL number for each element
    FluidCharacteristicNumbersUtilities::CalculateLocalCFL(r_model_part);

    // Check results
    const double tolerance = 2.0e-6;
    KRATOS_CHECK_NEAR(r_model_part.GetElement(1).GetValue(CFL_NUMBER), 0.186339, tolerance);
    KRATOS_CHECK_NEAR(r_model_part.GetElement(2).GetValue(CFL_NUMBER), 0.792324, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(FluidCharacteristicNumbersCalculateElementPrandtlNumber, FluidDynamicsApplicationFastSuite)
{
    // Create the test element
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("TestModelPart");
    TestFluidCharacteristicNumberElementSet(r_model_part);

    // Calculate the element Prantdl number
    const double prandtl_number = FluidCharacteristicNumbersUtilities::CalculateElementPrandtlNumber<true>(r_model_part.GetElement(1));

    // Check results
    const double tolerance = 1.0e-8;
    KRATOS_CHECK_NEAR(prandtl_number, 0.5, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(FluidCharacteristicNumbersCalculateElementPecletNumbers, FluidDynamicsApplicationFastSuite)
{
    // Create the test element
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("TestModelPart");
    TestFluidCharacteristicNumberElementSet(r_model_part);

    // Calculate the element Peclet numbers
    std::function<double(const Geometry<Node<3>>&)> avg_elem_function = ElementSizeCalculator<2,3>::AverageElementSize;
    const auto peclet_numbers = FluidCharacteristicNumbersUtilities::CalculateElementPecletNumbers<true, false>(
        r_model_part.GetElement(1),
        avg_elem_function);

    // Check results
    const double tolerance = 1.0e-8;
    KRATOS_CHECK_NEAR(std::get<0>(peclet_numbers), 0.055555555555, tolerance);
    KRATOS_CHECK_NEAR(std::get<1>(peclet_numbers), 0.027777777777, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(FluidCharacteristicNumbersCalculateElementThermalPecletNumber, FluidDynamicsApplicationFastSuite)
{
    // Create the test element
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("TestModelPart");
    TestFluidCharacteristicNumberElementSet(r_model_part);

    // Calculate the element Peclet numbers
    std::function<double(const Geometry<Node<3>>&)> avg_elem_function = ElementSizeCalculator<2,3>::AverageElementSize;
    const double k_peclet_number = FluidCharacteristicNumbersUtilities::CalculateElementThermalPecletNumber<true, false>(
        r_model_part.GetElement(1),
        avg_elem_function);

    // Check results
    const double tolerance = 1.0e-8;
    KRATOS_CHECK_NEAR(k_peclet_number, 0.027777777777, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(FluidCharacteristicNumbersCalculateElementViscousPecletNumber, FluidDynamicsApplicationFastSuite)
{
    // Create the test element
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("TestModelPart");
    TestFluidCharacteristicNumberElementSet(r_model_part);

    // Calculate the element Peclet numbers
    std::function<double(const Geometry<Node<3>>&)> avg_elem_function = ElementSizeCalculator<2,3>::AverageElementSize;
    const double mu_peclet_number = FluidCharacteristicNumbersUtilities::CalculateElementViscousPecletNumber<true, false>(
        r_model_part.GetElement(1),
        avg_elem_function);

    // Check results
    const double tolerance = 1.0e-8;
    KRATOS_CHECK_NEAR(mu_peclet_number, 0.055555555555, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(FluidCharacteristicNumbersCalculateElementFourierNumbers, FluidDynamicsApplicationFastSuite)
{
    // Set the current delta time to calculate the Fourier number
    const double current_dt = 1.0e-1;

    // Create the test element
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("TestModelPart");
    TestFluidCharacteristicNumberElementSet(r_model_part);

    // Calculate the element Peclet numbers
    std::function<double(const Geometry<Node<3>>&)> min_elem_function = ElementSizeCalculator<2,3>::MinimumElementSize;
    const auto fourier_numbers = FluidCharacteristicNumbersUtilities::CalculateElementFourierNumbers<true, false>(
        r_model_part.GetElement(1),
        min_elem_function,
        current_dt);

    // Check results
    const double tolerance = 1.0e-8;
    KRATOS_CHECK_NEAR(std::get<0>(fourier_numbers), 3.75, tolerance);
    KRATOS_CHECK_NEAR(std::get<1>(fourier_numbers), 7.5, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(FluidCharacteristicNumbersCalculateElementThermalFourierNumber, FluidDynamicsApplicationFastSuite)
{
    // Set the current delta time to calculate the Fourier number
    const double current_dt = 1.0e-1;
    
    // Create the test element
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("TestModelPart");
    TestFluidCharacteristicNumberElementSet(r_model_part);

    // Calculate the element Peclet numbers
    std::function<double(const Geometry<Node<3>>&)> min_elem_function = ElementSizeCalculator<2,3>::MinimumElementSize;
    const double thermal_fourier_number = FluidCharacteristicNumbersUtilities::CalculateElementThermalFourierNumber<true, false>(
        r_model_part.GetElement(1),
        min_elem_function,
        current_dt);

    // Check results
    const double tolerance = 1.0e-8;
    KRATOS_CHECK_NEAR(thermal_fourier_number, 7.5, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(FluidCharacteristicNumbersCalculateElementViscousFourierNumber, FluidDynamicsApplicationFastSuite)
{
    // Set the current delta time to calculate the Fourier number
    const double current_dt = 1.0e-1;

    // Create the test element
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("TestModelPart");
    TestFluidCharacteristicNumberElementSet(r_model_part);

    // Calculate the element Peclet numbers
    std::function<double(const Geometry<Node<3>>&)> min_elem_function = ElementSizeCalculator<2,3>::MinimumElementSize;
    const double viscous_fourier_number = FluidCharacteristicNumbersUtilities::CalculateElementViscousFourierNumber<true, false>(
        r_model_part.GetElement(1),
        min_elem_function,
        current_dt);

    // Check results
    const double tolerance = 1.0e-8;
    KRATOS_CHECK_NEAR(viscous_fourier_number, 3.75, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(FluidCharacteristicNumbersCalculateElementMachNumber, FluidDynamicsApplicationFastSuite)
{
    // Create the test element
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("TestModelPart");
    TestFluidCharacteristicNumberElementSet(r_model_part);

    // Set SOUND_VELOCITY
    for (auto& rNode: r_model_part.Nodes()) {
        rNode.SetValue(SOUND_VELOCITY, 2.0 * norm_2(rNode.FastGetSolutionStepValue(VELOCITY)));
    }

    // Calculate the element Mach number
    const double mach_number = FluidCharacteristicNumbersUtilities::CalculateElementMachNumber(r_model_part.GetElement(1));

    // Check results
    const double tolerance = 1.0e-6;
    KRATOS_CHECK_NEAR(mach_number, 0.357143, tolerance);
}

}
}