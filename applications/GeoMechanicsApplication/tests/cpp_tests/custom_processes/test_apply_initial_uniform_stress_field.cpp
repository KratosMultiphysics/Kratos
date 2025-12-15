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

#include "containers/model.h"
#include "custom_constitutive/incremental_linear_elastic_law.h"
#include "custom_constitutive/plane_strain.h"
#include "custom_constitutive/three_dimensional.h"
#include "custom_processes/apply_initial_uniform_stress_field.h"
#include "includes/expect.h"
#include "test_setup_utilities/model_setup_utilities.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h"
#include "tests/cpp_tests/test_utilities.h"

namespace Kratos::Testing
{
TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel,
       ApplyInitialUniformStressFieldProcessAppliesStressesToPlaneStrainElementsInModelPart)
{
    // Note that when creating the process using python, these parameters
    // should also contain the "model_part_name" field.
    Parameters parameters(R"({
        "value": [1.0, 2.0, 3.0, 4.0]
    })");

    Model model;
    auto& r_model_part = ModelSetupUtilities::CreateModelPartWithASingle2D3NElement(model);
    r_model_part.GetElement(1).GetProperties()[CONSTITUTIVE_LAW] =
        std::make_shared<GeoIncrementalLinearElasticLaw>(std::make_unique<PlaneStrain>());

    ApplyInitialUniformStressField process(r_model_part, parameters);
    process.ExecuteInitialize();

    std::vector<Vector> actual_stresses;
    r_model_part.GetElement(1).CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, actual_stresses,
                                                            r_model_part.GetProcessInfo());
    KRATOS_EXPECT_EQ(actual_stresses.size(), 3);

    const std::vector<double> expected_stress = {1.0, 2.0, 3.0, 4.0};
    for (const auto& stress : actual_stresses) {
        KRATOS_EXPECT_VECTOR_NEAR(stress, expected_stress, Defaults::absolute_tolerance);
    }
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel,
       ApplyInitialUniformStressFieldProcessAppliesStressesToThreeDimensionalElementsInModelPart)
{
    Parameters parameters(R"({
        "value": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    })");

    Model model;
    auto& r_model_part = ModelSetupUtilities::CreateModelPartWithASingle3D4NElement(model);
    r_model_part.GetElement(1).GetProperties()[CONSTITUTIVE_LAW] =
        std::make_shared<GeoIncrementalLinearElasticLaw>(std::make_unique<ThreeDimensional>());
    ApplyInitialUniformStressField process(r_model_part, parameters);
    process.ExecuteInitialize();

    std::vector<Vector> actual_stresses;

    r_model_part.GetElement(1).CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, actual_stresses,
                                                            r_model_part.GetProcessInfo());
    KRATOS_EXPECT_EQ(actual_stresses.size(), 4);

    const std::vector<double> expected_stress = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    for (const auto& stress : actual_stresses) {
        KRATOS_EXPECT_VECTOR_NEAR(stress, expected_stress, Defaults::absolute_tolerance);
    }
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel,
       ApplyInitialUniformStressFieldProcessAppliesStressesToPlaneStrainDiffOrderElementsInModelPart)
{
    Parameters parameters(R"({
        "value": [1.0, 2.0, 3.0, 4.0]
    })");

    Model model;
    auto& r_model_part = ModelSetupUtilities::CreateModelPartWithASingle2D6NDiffOrderElement(model);
    r_model_part.GetElement(1).GetProperties()[CONSTITUTIVE_LAW] =
        std::make_shared<GeoIncrementalLinearElasticLaw>(std::make_unique<PlaneStrain>());
    ApplyInitialUniformStressField process(r_model_part, parameters);
    process.ExecuteInitialize();

    std::vector<Vector> actual_stresses;

    r_model_part.GetElement(1).CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, actual_stresses,
                                                            r_model_part.GetProcessInfo());
    KRATOS_EXPECT_EQ(actual_stresses.size(), 3);

    const std::vector<double> expected_stress = {1.0, 2.0, 3.0, 4.0};
    for (const auto& stress : actual_stresses) {
        KRATOS_EXPECT_VECTOR_NEAR(stress, expected_stress, Defaults::absolute_tolerance);
    }
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel,
       ApplyInitialUniformStressFieldProcessThrowsUponConstructionWhenValuesAreMissing)
{
    Parameters parameters;
    Model      model;
    auto&      r_model_part = ModelSetupUtilities::CreateModelPartWithASingle2D3NElement(model);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(ApplyInitialUniformStressField(r_model_part, parameters),
                                      "Getting a value that does not exist. entry string : value");
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel,
       ApplyInitialUniformStressFieldThrowsUponConstructionWhenStressVectorSizeIsIncorrectForPlaneStrain)
{
    Parameters parameters(R"({
        "value": [1.0, 2.0, 3.0]
    })");

    Model model;
    auto& r_model_part = ModelSetupUtilities::CreateModelPartWithASingle2D3NElement(model);
    r_model_part.GetElement(1).GetProperties()[CONSTITUTIVE_LAW] =
        std::make_shared<GeoIncrementalLinearElasticLaw>(std::make_unique<PlaneStrain>());
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        ApplyInitialUniformStressField(r_model_part, parameters),
        "The size of the input stress vector for applying a uniform initial stress field must "
        "match the strain size of the constitutive law, which is 4, but is 3 for element 1 in "
        "model part "
        "'Main'. Please check the process parameters.");
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CheckInfoApplyInitialUniformStressField)
{
    // Arrange
    Model model;
    auto& r_empty_model_part = model.CreateModelPart("foo");
    const ApplyInitialUniformStressField process(r_empty_model_part, {R"({"value": [1.0, 2.0, 3.0, 4.0]})"});

    // Act & assert
    KRATOS_EXPECT_EQ(process.Info(), "ApplyInitialUniformStressField");
}

} // namespace Kratos::Testing
