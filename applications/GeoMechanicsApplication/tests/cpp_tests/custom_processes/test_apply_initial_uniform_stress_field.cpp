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

#include "custom_constitutive/incremental_linear_elastic_law.h"
#include "custom_constitutive/plane_strain.h"
#include "custom_constitutive/three_dimensional.h"
#include "custom_processes/apply_initial_uniform_stress_field.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities/model_setup_utilities.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(ApplyInitialUniformStressFieldProcessAppliesStressesToPlaneStrainElementsInModelParts,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Note that when creating the process using python, these parameters
    // should also contain the "model_part_name" field.
    Parameters parameters(R"({
        "value": [1.0, 2.0, 3.0, 4.0]
    })");

    Model model;
    auto& rModelPart = ModelSetupUtilities::CreateModelPartWithASingle2D3NElement(model);
    rModelPart.GetElement(1).GetProperties()[CONSTITUTIVE_LAW] =
        std::make_shared<GeoIncrementalLinearElasticLaw>(std::make_unique<PlaneStrain>());

    ApplyInitialUniformStressField process(rModelPart, parameters);
    process.ExecuteInitialize();

    std::vector<Vector> actual_stresses;
    rModelPart.GetElement(1).CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, actual_stresses,
                                                          rModelPart.GetProcessInfo());
    KRATOS_EXPECT_EQ(actual_stresses.size(), 3);

    const std::vector<double> expected_stress = {1.0, 2.0, 3.0, 4.0};
    for (const auto& stress : actual_stresses) {
        KRATOS_EXPECT_VECTOR_NEAR(stress, expected_stress, 1e-6);
    }
}

KRATOS_TEST_CASE_IN_SUITE(ApplyInitialUniformStressFieldProcessAppliesStressesToThreeDimensionalElementsInModelParts,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Parameters parameters(R"({
        "value": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    })");

    Model model;
    auto& rModelPart = ModelSetupUtilities::CreateModelPartWithASingle3D4NElement(model);
    rModelPart.GetElement(1).GetProperties()[CONSTITUTIVE_LAW] =
        std::make_shared<GeoIncrementalLinearElasticLaw>(std::make_unique<ThreeDimensional>());
    ApplyInitialUniformStressField process(rModelPart, parameters);
    process.ExecuteInitialize();

    std::vector<Vector> actual_stresses;

    rModelPart.GetElement(1).CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, actual_stresses,
                                                          rModelPart.GetProcessInfo());
    KRATOS_EXPECT_EQ(actual_stresses.size(), 4);

    const std::vector<double> expected_stress = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    for (const auto& stress : actual_stresses) {
        KRATOS_EXPECT_VECTOR_NEAR(stress, expected_stress, 1e-6);
    }
}

KRATOS_TEST_CASE_IN_SUITE(ApplyInitialUniformStressFieldProcessAppliesStressesToThreeDimensionalDiffOrderElementsInModelParts,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Parameters parameters(R"({
        "value": [1.0, 2.0, 3.0, 4.0]
    })");

    Model model;
    auto& rModelPart = ModelSetupUtilities::CreateModelPartWithASingle2D6NDiffOrderElement(model);
    rModelPart.GetElement(1).GetProperties()[CONSTITUTIVE_LAW] =
        std::make_shared<GeoIncrementalLinearElasticLaw>(std::make_unique<PlaneStrain>());
    ApplyInitialUniformStressField process(rModelPart, parameters);
    process.ExecuteInitialize();

    std::vector<Vector> actual_stresses;

    rModelPart.GetElement(1).CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, actual_stresses,
                                                          rModelPart.GetProcessInfo());
    KRATOS_EXPECT_EQ(actual_stresses.size(), 3);

    const std::vector<double> expected_stress = {1.0, 2.0, 3.0, 4.0};
    for (const auto& stress : actual_stresses) {
        KRATOS_EXPECT_VECTOR_NEAR(stress, expected_stress, 1e-6);
    }
}

KRATOS_TEST_CASE_IN_SUITE(ApplyInitialUniformStressFieldProcessThrowsUponConstructionWhenValuesAreMissing,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Parameters parameters;
    Model      model;
    auto&      rModelPart = ModelSetupUtilities::CreateModelPartWithASingle2D3NElement(model);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(ApplyInitialUniformStressField process(rModelPart, parameters),
                                      "Getting a value that does not exist. entry string : value");
}

KRATOS_TEST_CASE_IN_SUITE(ApplyInitialUniformStressFieldThrowsUponConstructionWhenStressVectorSizeIsIncorrectForPlaneStrain,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Parameters parameters(R"({
        "value": [1.0, 2.0, 3.0]
    })");

    Model model;
    auto& rModelPart = ModelSetupUtilities::CreateModelPartWithASingle2D3NElement(model);
    rModelPart.GetElement(1).GetProperties()[CONSTITUTIVE_LAW] =
        std::make_shared<GeoIncrementalLinearElasticLaw>(std::make_unique<PlaneStrain>());
    rModelPart.GetElement(1).Initialize(rModelPart.GetProcessInfo());
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        ApplyInitialUniformStressField(rModelPart, parameters),
        "The size of the input stress vector for applying a uniform initial stress field must "
        "match the strain size of the constitutive law, which is 4, but is 3 for element 1 in model part "
        "'Main'. Please check the process parameters.");
}

} // namespace Kratos::Testing
