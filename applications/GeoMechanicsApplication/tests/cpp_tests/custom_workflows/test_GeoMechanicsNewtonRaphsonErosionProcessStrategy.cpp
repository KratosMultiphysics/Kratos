// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Jonathan Nuttall
//

// System includes
#include <map>

/* Project includes */
#include "containers/model.h"
#include "custom_workflows/dgeoflow.h"
#include "geo_mechanics_application_variables.h"
#include "includes/kratos_parameters.h"
#include "tests/cpp_tests/flow_stubs.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

namespace Kratos::Testing
{

using namespace Kratos;

class ParametrizedDGeoFlowTests : public ::testing::TestWithParam<std::string>
{
};

TEST_P(ParametrizedDGeoFlowTests, ErosionProcessStrategy)
{
    // Arrange
    auto working_directory =
        "./applications/GeoMechanicsApplication/tests/test_compare_sellmeijer/" + GetParam();
    auto project_file = "ProjectParameters.json";

    auto execute = KratosExecute();

    std::vector<double> progress_report_values;
    auto                reportProgress = [&progress_report_values](double Progress) {
        progress_report_values.push_back(Progress);
    };

    std::vector<std::string> reported_textual_progress;
    auto report_textual_progress = [&reported_textual_progress](const char* pMessage) {
        reported_textual_progress.emplace_back(pMessage);
    };

    const KratosExecute::CriticalHeadInfo  critical_head_info(3, 4, 0.1);
    const KratosExecute::CallBackFunctions call_back_functions(
        &flow_stubs::emptyLog, reportProgress, report_textual_progress, &flow_stubs::emptyCancel);

    // Act
    const auto status = execute.ExecuteFlowAnalysis(working_directory, project_file, critical_head_info,
                                                    "PorousDomain.Left_head", call_back_functions);

    // Assert
    KRATOS_EXPECT_EQ(status, 0);
    KRATOS_EXPECT_EQ(reported_textual_progress.front(), "Calculating head level 3m (1/12)");
    KRATOS_EXPECT_EQ(reported_textual_progress.back(), "Calculating head level 3.8m (9/12)");
    KRATOS_EXPECT_EQ(reported_textual_progress.size(), 9);
    KRATOS_EXPECT_DOUBLE_EQ(progress_report_values.front(), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(progress_report_values.back(), 0.75);
    KRATOS_EXPECT_EQ(progress_report_values.size(), 10);
}

KRATOS_TEST_CASE_IN_SUITE(DGeoFlowParsesSeepageAndNormalFluxLoads, KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("PorousDomain");
    r_model_part.AddNodalSolutionStepVariable(NORMAL_FLUID_FLUX);
    r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION_X);
    r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION_Y);
    r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION_Z);
    r_model_part.CreateSubModelPart("flux_boundary");
    r_model_part.CreateSubModelPart("seepage_boundary");
    r_model_part.CreateSubModelPart("gravity");

    Parameters project_parameters(R"(
    {
        "processes": {
            "constraints_process_list": [],
            "loads_process_list": [
                {
                    "process_name": "ApplyScalarConstraintTableProcess",
                    "Parameters": {
                        "model_part_name": "PorousDomain.flux_boundary",
                        "variable_name": "NORMAL_FLUID_FLUX",
                        "value": -5.0,
                        "table": 1
                    }
                },
                {
                    "process_name": "ApplySeepageBoundaryProcess",
                    "Parameters": {
                        "model_part_name": "PorousDomain.seepage_boundary"
                    }
                },
                {
                    "process_name": "ApplyVectorConstraintTableProcess",
                    "Parameters": {
                        "model_part_name": "PorousDomain.gravity",
                        "variable_name": "VOLUME_ACCELERATION"
                    }
                }
            ]
        }
    })");

    KratosExecute execute;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(execute.ParseProcesses(r_model_part, project_parameters),
                                      "DGeoFlow only supports constant NORMAL_FLUID_FLUX values.");

    project_parameters["processes"]["loads_process_list"][0]["Parameters"]["table"].SetInt(0);
    execute.ParseProcesses(r_model_part, project_parameters);
}

KRATOS_TEST_CASE_IN_SUITE(DGeoFlowParsesSeepageConstraints, KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("PorousDomain");
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION_X);
    r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION_Y);
    r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION_Z);
    r_model_part.CreateSubModelPart("fixed_head");
    r_model_part.CreateSubModelPart("gravity");

    Parameters project_parameters(R"(
    {
        "processes": {
            "constraints_process_list": [
                {
                    "Parameters": {
                        "model_part_name": "PorousDomain.fixed_head",
                        "fluid_pressure_type": "Uniform",
                        "value": 0.0
                    }
                },
                {
                    "process_name": "ApplySeepageBoundaryProcess",
                    "Parameters": {
                        "model_part_name": "PorousDomain.seepage_boundary"
                    }
                }
            ],
            "loads_process_list": [
                {
                    "process_name": "ApplyVectorConstraintTableProcess",
                    "Parameters": {
                        "model_part_name": "PorousDomain.gravity",
                        "variable_name": "VOLUME_ACCELERATION"
                    }
                }
            ]
        }
    })");

    KratosExecute execute;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        execute.ParseProcesses(r_model_part, project_parameters),
        "ApplySeepageBoundaryProcess: \"PorousDomain.seepage_boundary\" does not correspond to any model part");

    r_model_part.CreateSubModelPart("seepage_boundary");
    KratosExecute valid_execute;
    valid_execute.ParseProcesses(r_model_part, project_parameters);
}

INSTANTIATE_TEST_SUITE_P(KratosGeoMechanicsIntegrationSuite,
                         ParametrizedDGeoFlowTests,
                         ::testing::Values("HeightAquiferD10L30line"));

} // namespace Kratos::Testing
