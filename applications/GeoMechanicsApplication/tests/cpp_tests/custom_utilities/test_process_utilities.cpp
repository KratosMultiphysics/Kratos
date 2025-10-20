// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//

#include "custom_utilities/process_utilities.h"
#include "testing/testing.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(GetModelPartsFromSettings_SingleModelPart, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    model.CreateModelPart("Main");

    Parameters settings(R"(
        {
            "model_part_name": "Main"
        })");

    auto model_parts = ProcessUtilities::GetModelPartsFromSettings(model, settings, "TestProcess");

    KRATOS_CHECK_EQUAL(model_parts.size(), 1);
    KRATOS_CHECK_EQUAL(model_parts[0].get().Name(), "Main");
}

KRATOS_TEST_CASE_IN_SUITE(GetModelPartsFromSettings_ListOfModelParts, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    model.CreateModelPart("Part1");
    model.CreateModelPart("Part2");

    Parameters settings(R"(
        {
            "model_part_name_list": ["Part1", "Part2"]
        })");

    auto model_parts = ProcessUtilities::GetModelPartsFromSettings(model, settings, "TestProcess");

    KRATOS_CHECK_EQUAL(model_parts.size(), 2);
    KRATOS_CHECK_EQUAL(model_parts[0].get().Name(), "Part1");
    KRATOS_CHECK_EQUAL(model_parts[1].get().Name(), "Part2");
}

KRATOS_TEST_CASE_IN_SUITE(GetModelPartsFromSettings_BothParametersPresent_Throws, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;

    Parameters settings(R"(
        {
            "model_part_name": "Part1",
            "model_part_name_list": ["Part2"]
        })");

    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        ProcessUtilities::GetModelPartsFromSettings(model, settings, "TestProcess"),
        "The parameters 'model_part_name' and 'model_part_name_list' are mutually exclusive");
}

KRATOS_TEST_CASE_IN_SUITE(GetModelPartsFromSettings_MissingParameters_Throws, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model      model;
    Parameters settings("{}");

    KRATOS_CHECK_EXCEPTION_IS_THROWN(ProcessUtilities::GetModelPartsFromSettings(model, settings, "TestProcess"),
                                     "Please specify 'model_part_name' or 'model_part_name_list'");
}

KRATOS_TEST_CASE_IN_SUITE(GetModelPartsFromSettings_EmptyList_Throws, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model      model;
    Parameters settings(R"(
        {
            "model_part_name_list": []
        })");

    KRATOS_CHECK_EXCEPTION_IS_THROWN(ProcessUtilities::GetModelPartsFromSettings(model, settings, "TestProcess"), "The parameters 'model_part_name_list' needs to contain at least one model part name for TestProcess");
}

} // namespace Kratos::Testing