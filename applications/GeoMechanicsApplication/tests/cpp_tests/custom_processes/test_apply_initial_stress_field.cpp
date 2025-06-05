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
#include "tests/cpp_tests/test_utilities/model_setup_utilities.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(ApplyInitialUniformStressFieldProcessAppliesStressesToElementsInModelParts,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Parameters parameters(R"({
        "model_part_name": "test_model_part",
        "variable_name": "CAUCHY_STRESS_VECTOR",
        "value": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    })");

    Model model;
    auto& rModelPart = ModelSetupUtilities::CreateModelPartWithASingle2D3NElement(model);

    ApplyInitialUniformStressField process(rModelPart, parameters);
    process.ExecuteInitialize();

    std::vector<Vector> actual_stresses;

    auto p_element = rModelPart.Elements()[0];

    rModelPart.GetElement(1).CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, actual_stresses,
                                                          rModelPart.GetProcessInfo());
    KRATOS_EXPECT_EQ(actual_stresses.size(), 3);

    for (const auto& stress : actual_stresses) {
        Vector expected_stress(6);
        expected_stress[0] = 1.0;
        expected_stress[1] = 2.0;
        expected_stress[2] = 3.0;
        expected_stress[3] = 4.0;
        expected_stress[4] = 5.0;
        expected_stress[5] = 6.0;
        KRATOS_EXPECT_VECTOR_NEAR(stress, expected_stress, 1e-6);
    }
}

} // namespace Kratos::Testing
