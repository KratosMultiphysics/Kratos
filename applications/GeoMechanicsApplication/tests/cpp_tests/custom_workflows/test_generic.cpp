// KRATOS___
//     //   ) )
//    //         ___      ___
//  //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//

#include "custom_workflows/generic.h"

#include <string>

#include "custom_utilities/json_process_info_parser.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/stub_input_utility.h"
#include "tests/cpp_tests/stub_time_loop_executor.h"

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(GenericWorkflowDispatchesNativeProcessesFromProjectParameters,
                          KratosGeoMechanicsFastSuite)
{
    const std::string project_parameters = R"(
    {
        "problem_data": {
            "start_time": 0.0,
            "end_time": 1.0
        },
        "solver_settings": {
            "model_part_name": "test",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "mesh"
            },
            "problem_domain_sub_model_part_list": [],
            "time_stepping": {
                "time_step": 1.0
            },
            "number_cycles": 1,
            "reduction_factor": 0.5,
            "increase_factor": 2.0,
            "min_iterations": 1,
            "max_iterations": 2,
            "reset_totals": false
        },
        "processes": {
            "constraints_process_list": [
                {
                    "process_name": "ApplyConstantHydrostaticPressureProcess",
                    "Parameters": {
                        "model_part_name": "test",
                        "fluid_pressure_type": "Uniform",
                        "value": 12.3
                    }
                }
            ],
            "loads_process_list": [
                {
                    "process_name": "GeoApplyConstantScalarValueProcess",
                    "Parameters": {
                        "model_part_name": "test",
                        "variable_name": "NORMAL_FLUID_FLUX",
                        "value": 2.0
                    }
                }
            ],
            "auxiliary_process_list": [
                {
                    "process_name": "ApplySeepageBoundaryProcess",
                    "Parameters": {
                        "model_part_name": "test"
                    }
                }
            ]
        }
    })";

    KratosGeneric generic{std::make_unique<StubInputUtility>(project_parameters),
                          std::make_unique<JsonProcessInfoParser>(),
                          std::make_unique<StubTimeLoopExecutor>()};
    std::string log;
    const auto status = generic.Run({}, {}, [&log](const char* pMessage) { log = pMessage; },
                                    [](double) {}, [](const char*) {}, [] { return false; });

    KRATOS_EXPECT_EQ(status, 1);
    KRATOS_EXPECT_NE(log.find("strategy_type"), std::string::npos);
}

} // namespace Kratos::Testing
