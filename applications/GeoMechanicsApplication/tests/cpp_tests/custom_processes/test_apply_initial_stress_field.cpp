// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#include "custom_processes/apply_initial_stress_field.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(ApplyInitialStressFieldProcessCanBeCreated, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Parameters parameters(R"({
        "model_part_name": "test_model_part",
        "variable_name": "CAUCHY_STRESS_VECTOR",
        "value": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    })");

    Model model;
    ModelPart& r_model_part = model.CreateModelPart("test_model_part");

    ApplyInitialStressField process(r_model_part, parameters);


}

} // namespace Kratos::Testing
