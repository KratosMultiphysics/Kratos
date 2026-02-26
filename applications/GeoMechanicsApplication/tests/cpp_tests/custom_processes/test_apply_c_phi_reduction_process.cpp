// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra
//                   Gennady Markelov
//
#include "containers/model.h"
#include "custom_processes/apply_c_phi_reduction_process.h"
#include "geo_mechanics_application_variables.h"
#include "geometries/triangle_2d_3.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/stub_linear_elastic_law.h"
#include <boost/numeric/ublas/assignment.hpp>

using namespace Kratos;

namespace
{
ModelPart& CreateModelPartWithElements(Model& rModel)
{
    auto& result       = rModel.CreateModelPart("dummy");
    auto  p_properties = result.CreateNewProperties(0);

    const auto p_node1 = make_intrusive<Node>(1, 0.0, 0.0, 0.0);
    const auto p_node2 = make_intrusive<Node>(2, 1.0, 0.0, 0.0);
    const auto p_node3 = make_intrusive<Node>(3, 1.0, 1.0, 0.0);
    const auto p_node4 = make_intrusive<Node>(4, 0.0, 1.0, 0.0);

    result.AddElement(make_intrusive<Element>(
        1, std::make_shared<Triangle2D3<Node>>(p_node1, p_node2, p_node3), p_properties));
    result.AddElement(make_intrusive<Element>(
        3, std::make_shared<Triangle2D3<Node>>(p_node1, p_node3, p_node4), p_properties));

    return result;
}

ModelPart& PrepareCPhiTestModelPart(Model& rModel)
{
    auto& result = CreateModelPartWithElements(rModel);

    auto& r_model_part_properties = result.GetProperties(0);
    auto  p_dummy_law             = std::make_shared<Testing::StubLinearElasticLaw>();

    r_model_part_properties.SetValue(CONSTITUTIVE_LAW, p_dummy_law);
    Vector umat_parameters(6);
    umat_parameters <<= 10000000, 0.2, 10.0, 25.0, 25.0, 1000;
    r_model_part_properties.SetValue(UMAT_PARAMETERS, umat_parameters);
    r_model_part_properties.SetValue(INDEX_OF_UMAT_C_PARAMETER, 3);
    r_model_part_properties.SetValue(INDEX_OF_UMAT_PHI_PARAMETER, 4);

    return result;
}

void CheckReducedCPhi(const ModelPart& rModelPart, double COrig, double PhiOrig, double ReductionFactor)
{
    for (Element& rElement : rModelPart.Elements()) {
        const auto& element_properties     = rElement.GetProperties();
        const auto& umat_properties_vector = element_properties.GetValue(UMAT_PARAMETERS);
        const auto  c_index   = element_properties.GetValue(INDEX_OF_UMAT_C_PARAMETER) - 1;
        const auto  phi_index = element_properties.GetValue(INDEX_OF_UMAT_PHI_PARAMETER) - 1;

        KRATOS_EXPECT_DOUBLE_EQ(umat_properties_vector(c_index), ReductionFactor * COrig);

        const double phi_rad = MathUtils<>::DegreesToRadians(PhiOrig);
        const double tan_phi = std::tan(phi_rad);
        KRATOS_EXPECT_DOUBLE_EQ(std::tan(MathUtils<>::DegreesToRadians(umat_properties_vector(phi_index))),
                                ReductionFactor * tan_phi);
    }
}

} // namespace

namespace Kratos::Testing
{
KRATOS_TEST_CASE_IN_SUITE(CheckCAndPhiReducedAfterCallingApplyCPhiReductionProcess, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model       model;
    const auto& r_model_part = PrepareCPhiTestModelPart(model);

    const auto                parameters = Parameters{R"({"model_part_name" : "dummy"})"};
    ApplyCPhiReductionProcess process{model, parameters};
    process.ExecuteInitializeSolutionStep();

    CheckReducedCPhi(r_model_part, 10.0, 25.0, 0.9);
}

KRATOS_TEST_CASE_IN_SUITE(CheckCAndPhiTwiceReducedAfterCallingApplyCPhiReductionProcessTwice,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model       model;
    const auto& r_model_part = PrepareCPhiTestModelPart(model);

    const auto                parameters = Parameters{R"({"model_part_name" : "dummy"})"};
    ApplyCPhiReductionProcess process{model, parameters};
    process.ExecuteInitializeSolutionStep();
    process.ExecuteFinalizeSolutionStep();
    process.ExecuteInitializeSolutionStep();

    CheckReducedCPhi(r_model_part, 10.0, 25.0, 0.8);
}

KRATOS_TEST_CASE_IN_SUITE(CheckFailureUmatInputsApplyCPhiReductionProcess, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    auto& r_model_part = CreateModelPartWithElements(model);

    auto& r_model_part_properties = r_model_part.GetProperties(0);
    auto  p_dummy_law             = std::make_shared<Testing::StubLinearElasticLaw>();
    r_model_part_properties.SetValue(CONSTITUTIVE_LAW, p_dummy_law);
    const auto parameters = Parameters{R"({"model_part_name" : "dummy"})"};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        (ApplyCPhiReductionProcess{model, parameters}.Check()),
        "UMAT_PARAMETERS does not exist in the model part property with Id 0.")

    Vector umat_parameters(6);
    umat_parameters <<= 10000000, 0.2, 10.0, 25.0, 25.0, 1000;
    r_model_part_properties.SetValue(UMAT_PARAMETERS, umat_parameters);

    // checking settings for Phi
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        (ApplyCPhiReductionProcess{model, parameters}.Check()),
        "INDEX_OF_UMAT_PHI_PARAMETER does not exist in the model part property with Id 0.")

    r_model_part_properties.SetValue(INDEX_OF_UMAT_PHI_PARAMETER, 0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN((ApplyCPhiReductionProcess{model, parameters}.Check()),
                                      "INDEX_OF_UMAT_PHI_PARAMETER in the model part property with "
                                      "Id 0 has an invalid value: 0 is out of the range [1, 6].")

    r_model_part_properties.SetValue(INDEX_OF_UMAT_PHI_PARAMETER, 7);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN((ApplyCPhiReductionProcess{model, parameters}.Check()),
                                      "INDEX_OF_UMAT_PHI_PARAMETER in the model part property with "
                                      "Id 0 has an invalid value: 7 is out of the range [1, 6].")

    r_model_part_properties.SetValue(INDEX_OF_UMAT_PHI_PARAMETER, 4);

    // checking settings for c
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        (ApplyCPhiReductionProcess{model, parameters}.Check()),
        "INDEX_OF_UMAT_C_PARAMETER does not exist in the model part property with Id 0.")

    r_model_part_properties.SetValue(INDEX_OF_UMAT_C_PARAMETER, 0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN((ApplyCPhiReductionProcess{model, parameters}.Check()),
                                      "INDEX_OF_UMAT_C_PARAMETER in the model part property with "
                                      "Id 0 has an invalid value: 0 is out of the range [1, 6].")

    r_model_part_properties.SetValue(INDEX_OF_UMAT_C_PARAMETER, 7);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN((ApplyCPhiReductionProcess{model, parameters}.Check()),
                                      "INDEX_OF_UMAT_C_PARAMETER in the model part property with "
                                      "Id 0 has an invalid value: 7 is out of the range [1, 6].")

    // checking Phi value
    umat_parameters(3) = -0.0001;
    r_model_part_properties.SetValue(UMAT_PARAMETERS, umat_parameters);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN((ApplyCPhiReductionProcess{model, parameters}.ExecuteInitializeSolutionStep()), "Friction angle Phi in the model part property with Id 0 has an invalid value: -0.0001 is out of range [0,90] (degrees).")

    umat_parameters(3) = 90.0001;
    r_model_part_properties.SetValue(UMAT_PARAMETERS, umat_parameters);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN((ApplyCPhiReductionProcess{model, parameters}.ExecuteInitializeSolutionStep()), "Friction angle Phi in the model part property with Id 0 has an invalid value: 90.0001 is out of range [0,90] (degrees).")

    umat_parameters(3) = 25.0;
    r_model_part_properties.SetValue(UMAT_PARAMETERS, umat_parameters);

    // checking c value
    r_model_part_properties.SetValue(INDEX_OF_UMAT_C_PARAMETER, 3);
    umat_parameters(2) = -0.00001;
    r_model_part_properties.SetValue(UMAT_PARAMETERS, umat_parameters);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        (ApplyCPhiReductionProcess{model, parameters}.ExecuteInitializeSolutionStep()),
        "Cohesion C out of range: -1e-05")
}

KRATOS_TEST_CASE_IN_SUITE(CheckFailureEmptyModelPartApplyCPhiReductionProcess, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    model.CreateModelPart("dummy");
    const auto parameters = Parameters{R"({"model_part_name" : "dummy"})"};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        (ApplyCPhiReductionProcess{model, parameters}.Check()),
        "None of the provided model parts contains at least one element. A c-phi reduction "
        "analysis requires at least one element.")
}

KRATOS_TEST_CASE_IN_SUITE(CheckReturnsZeroForValidModelPartApplyCPhiReductionProcess, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    PrepareCPhiTestModelPart(model);
    const auto parameters = Parameters{R"({"model_part_name_list" : ["dummy"]})"};

    ApplyCPhiReductionProcess process{model, parameters};
    KRATOS_CHECK_EQUAL(process.Check(), 0);
}

KRATOS_TEST_CASE_IN_SUITE(CheckFailureNegativeReductionFactorApplyCPhiReductionProcess,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    PrepareCPhiTestModelPart(model);

    const auto parameters = Parameters{R"({"model_part_name_list" : ["dummy"]})"};

    ApplyCPhiReductionProcess process{model, parameters};
    for (size_t i = 0; i < 9; ++i) {
        process.ExecuteInitializeSolutionStep();
        process.ExecuteFinalizeSolutionStep();
    }
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        process.ExecuteInitializeSolutionStep(),
        "Reduction factor should not drop below 0.01, calculation stopped.");
}

KRATOS_TEST_CASE_IN_SUITE(CheckFailureTooSmallReductionIncrementApplyCPhiReductionProcess,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    auto& r_model_part = PrepareCPhiTestModelPart(model);
    // Number of cycles should be larger than one to trigger the step halving of the
    // reduction increment
    r_model_part.GetProcessInfo().GetValue(NUMBER_OF_CYCLES) = 10;
    const auto parameters = Parameters{R"({"model_part_name_list" : ["dummy"]})"};

    ApplyCPhiReductionProcess process{model, parameters};
    for (size_t i = 0; i < 6; ++i) {
        process.ExecuteInitializeSolutionStep();
        process.ExecuteFinalizeSolutionStep();
    }

    // After halving the initial reduction increment (0.1) for the seventh time, it becomes
    // 0.00078125, so it should throw an exception
    // The final safety factor can be calculated as (1.0 / (1 - 0.1 * 0.5 - 0.1 * 0.5^2 ... + 0.1*0.5^6)) = 1.10919
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(process.ExecuteInitializeSolutionStep(),
                                      "Reduction increment should not drop below 0.001, "
                                      "calculation stopped. Final safety factor = 1.10919");
}

KRATOS_TEST_CASE_IN_SUITE(CheckInfoApplyCPhiReductionProcess, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    model.CreateModelPart("foo");
    const auto                      parameters = Parameters{R"({"model_part_name" : "foo"})"};
    const ApplyCPhiReductionProcess process{model, parameters};

    // Act & Assert
    KRATOS_EXPECT_EQ(process.Info(), "ApplyCPhiReductionProcess");
}

} // namespace Kratos::Testing
