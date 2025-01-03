// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//
#include "containers/model.h"
#include "custom_processes/apply_k0_procedure_process.h"
#include "geo_mechanics_application_variables.h"
#include "includes/element.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/stub_linear_elastic_law.h"
#include "tests/cpp_tests/test_utilities.h"
#include <custom_constitutive/incremental_linear_elastic_law.h>

#include <boost/numeric/ublas/assignment.hpp>
#include <gmock/gmock.h>

namespace
{

using namespace Kratos;

class MockConstitutiveLaw : public GeoIncrementalLinearElasticLaw
{
public:
    MOCK_METHOD(std::size_t, WorkingSpaceDimension, (), (override));
};

class StubElement : public Element
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(StubElement);

    void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                      std::vector<Vector>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override
    {
        rOutput = {mStressVector};
    }

    using Element::CalculateOnIntegrationPoints;

    void SetValuesOnIntegrationPoints(const Variable<Vector>&    rVariable,
                                      const std::vector<Vector>& rValues,
                                      const ProcessInfo&         rCurrentProcessInfo) override
    {
        mStressVector = rValues[0];
    }

    using Element::SetValuesOnIntegrationPoints;

private:
    Vector mStressVector;
};

Vector ApplyK0ProcedureOnStubElement(const Properties::Pointer& rProperties, const Vector& rInitialStressVector)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("main");
    auto  p_element    = make_intrusive<StubElement>();
    p_element->SetProperties(rProperties);
    r_model_part.AddElement(p_element);

    p_element->SetValuesOnIntegrationPoints(CAUCHY_STRESS_VECTOR, {rInitialStressVector},
                                            r_model_part.GetProcessInfo());

    const auto              k0_settings = Parameters{};
    ApplyK0ProcedureProcess process{r_model_part, k0_settings};

    // Act
    process.ExecuteFinalizeSolutionStep();

    // Assert
    std::vector<Vector> actual_stress_vector;
    p_element->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, actual_stress_vector,
                                            r_model_part.GetProcessInfo());

    return actual_stress_vector[0];
}

} // namespace

using namespace Kratos;

namespace
{

ModelPart& PrepareTestModelPart(Model& rModel)
{
    auto& result                  = rModel.CreateModelPart("dummy");
    auto  p_dummy_law             = std::make_shared<Testing::StubLinearElasticLaw>();
    auto  p_model_part_properties = result.pGetProperties(0);
    p_model_part_properties->SetValue(CONSTITUTIVE_LAW, p_dummy_law);

    auto p_element = make_intrusive<StubElement>();
    p_element->SetProperties(p_model_part_properties);
    result.AddElement(p_element);

    return result;
}

bool ElementConsidersDiagonalEntriesOnlyAndNoShear(const Element& rElement)
{
    auto p_constitutive_law =
        dynamic_cast<const GeoLinearElasticLaw*>(rElement.GetProperties().GetValue(CONSTITUTIVE_LAW).get());
    return p_constitutive_law && p_constitutive_law->GetConsiderDiagonalEntriesOnlyAndNoShear();
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(AllElementsConsiderDiagonalEntriesOnlyAndNoShearWhenUseStandardProcedureFlagIsNotDefined,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    auto& r_model_part = PrepareTestModelPart(model);

    Parameters k0_settings; // 'use_standard_procedure' is not defined, assume it to be true

    ApplyK0ProcedureProcess process{r_model_part, k0_settings};
    process.ExecuteInitialize();

    KRATOS_EXPECT_TRUE(ElementConsidersDiagonalEntriesOnlyAndNoShear(r_model_part.Elements()[0]))
}

KRATOS_TEST_CASE_IN_SUITE(AllElementsConsiderDiagonalEntriesOnlyAndNoShearWhenUsingStandardProcedure,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    auto& r_model_part = PrepareTestModelPart(model);

    Parameters k0_settings{R"({"use_standard_procedure": true})"};

    ApplyK0ProcedureProcess process{r_model_part, k0_settings};
    process.ExecuteInitialize();

    KRATOS_EXPECT_TRUE(ElementConsidersDiagonalEntriesOnlyAndNoShear(r_model_part.Elements()[0]))
}

KRATOS_TEST_CASE_IN_SUITE(NoneOfElementsConsiderDiagonalEntriesOnlyAndNoShearWhenNotUsingStandardProcedure,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    auto& r_model_part = PrepareTestModelPart(model);

    Parameters k0_settings{R"({"use_standard_procedure": false})"};

    ApplyK0ProcedureProcess process{r_model_part, k0_settings};
    process.ExecuteInitialize();

    KRATOS_EXPECT_FALSE(ElementConsidersDiagonalEntriesOnlyAndNoShear(r_model_part.Elements()[0]))
}

KRATOS_TEST_CASE_IN_SUITE(UseStandardProcedureFlagIsInEffectDuringProcessExecutionOnly,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    auto& r_model_part = PrepareTestModelPart(model);

    Parameters k0_settings{R"({"use_standard_procedure": true})"};

    ApplyK0ProcedureProcess process{r_model_part, k0_settings};
    process.ExecuteInitialize(); // start considering diagonal entries only and no shear
    process.ExecuteFinalize();   // stop considering diagonal entries only and no shear

    KRATOS_EXPECT_FALSE(ElementConsidersDiagonalEntriesOnlyAndNoShear(r_model_part.Elements()[0]))
}

KRATOS_TEST_CASE_IN_SUITE(K0ProcedureIsAppliedCorrectlyWithK0_NC, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto p_properties = std::make_shared<Properties>();
    p_properties->SetValue(K0_NC, 0.5);
    p_properties->SetValue(K0_MAIN_DIRECTION, 1);
    Vector initial_stress_vector{4};
    initial_stress_vector <<= 0.0, -10.0, 0.0, 27.0;

    // Act
    const auto actual_stress_vector = ApplyK0ProcedureOnStubElement(p_properties, initial_stress_vector);

    // Assert
    Vector expected_stress_vector{4};
    expected_stress_vector <<= -5.0, -10.0, -5.0, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(actual_stress_vector, expected_stress_vector, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(K0ProcedureIsAppliedCorrectlyWithK0_NC_3D, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto p_properties = std::make_shared<Properties>();
    p_properties->SetValue(K0_NC, 0.5);
    p_properties->SetValue(K0_MAIN_DIRECTION, 2);
    Vector initial_stress_vector{6};
    initial_stress_vector <<= 0.0, -10.0, -10.0, 27.0, 10.0, 5.0;

    // Act
    const auto actual_stress_vector = ApplyK0ProcedureOnStubElement(p_properties, initial_stress_vector);

    // Assert
    Vector expected_stress_vector{6};
    expected_stress_vector <<= -5.0, -5.0, -10.0, 0.0, 0.0, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(actual_stress_vector, expected_stress_vector, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(K0ProcedureIsAppliedCorrectlyWithPhi, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto p_properties = std::make_shared<Properties>();
    p_properties->SetValue(INDEX_OF_UMAT_PHI_PARAMETER, 1);
    Vector umat_parameters{1};
    umat_parameters[0] = 30.0;
    p_properties->SetValue(UMAT_PARAMETERS, umat_parameters);
    p_properties->SetValue(K0_MAIN_DIRECTION, 1);

    Vector initial_stress_vector{4};
    initial_stress_vector <<= 0.0, -10.0, 0.0, 27.0;

    // Act
    const auto actual_stress_vector = ApplyK0ProcedureOnStubElement(p_properties, initial_stress_vector);

    // Assert
    Vector expected_stress_vector{4};
    expected_stress_vector <<= -5.0, -10.0, -5.0, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(actual_stress_vector, expected_stress_vector, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(K0ProcedureIsAppliedCorrectlyWithPhi_3D, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto p_properties = std::make_shared<Properties>();
    p_properties->SetValue(INDEX_OF_UMAT_PHI_PARAMETER, 1);
    Vector umat_parameters{1};
    umat_parameters[0] = 30.0;
    p_properties->SetValue(UMAT_PARAMETERS, umat_parameters);
    p_properties->SetValue(K0_MAIN_DIRECTION, 2);

    Vector initial_stress_vector{6};
    initial_stress_vector <<= 0.0, -10.0, -10.0, 27.0, 10.0, 5.0;

    // Act
    const auto actual_stress_vector = ApplyK0ProcedureOnStubElement(p_properties, initial_stress_vector);

    // Assert
    Vector expected_stress_vector{6};
    expected_stress_vector <<= -5.0, -5.0, -10.0, 0.0, 0.0, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(actual_stress_vector, expected_stress_vector, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(K0ProcedureIsAppliedCorrectlyWithK0_NCandOCR, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto p_properties = std::make_shared<Properties>();
    p_properties->SetValue(K0_NC, 0.5);
    p_properties->SetValue(K0_MAIN_DIRECTION, 1);
    p_properties->SetValue(OCR, 1.5);
    Vector initial_stress_vector{4};
    initial_stress_vector <<= 0.0, -10.0, 0.0, 27.0;

    // Act
    const auto actual_stress_vector = ApplyK0ProcedureOnStubElement(p_properties, initial_stress_vector);

    // Assert
    Vector expected_stress_vector{4};
    expected_stress_vector <<= -7.5, -10.0, -7.5, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(actual_stress_vector, expected_stress_vector, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(K0ProcedureIsAppliedCorrectlyWithK0_NCandOCR_3D, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto p_properties = std::make_shared<Properties>();
    p_properties->SetValue(K0_NC, 0.5);
    p_properties->SetValue(K0_MAIN_DIRECTION, 2);
    p_properties->SetValue(OCR, 1.5);
    Vector initial_stress_vector{6};
    initial_stress_vector <<= 0.0, -10.0, -10.0, 27.0, 10.0, 5.0;

    // Act
    const auto actual_stress_vector = ApplyK0ProcedureOnStubElement(p_properties, initial_stress_vector);

    // Assert
    Vector expected_stress_vector{6};
    expected_stress_vector <<= -7.5, -7.5, -10.0, 0.0, 0.0, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(actual_stress_vector, expected_stress_vector, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(K0ProcedureIsAppliedCorrectlyWithK0_NCandOCRandNu_UR, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto p_properties = std::make_shared<Properties>();
    p_properties->SetValue(K0_NC, 0.5);
    p_properties->SetValue(K0_MAIN_DIRECTION, 1);
    p_properties->SetValue(OCR, 1.5);
    p_properties->SetValue(POISSON_UNLOADING_RELOADING, 0.25);
    Vector initial_stress_vector{4};
    initial_stress_vector <<= 0.0, -10.0, 0.0, 27.0;

    // Act
    const auto actual_stress_vector = ApplyK0ProcedureOnStubElement(p_properties, initial_stress_vector);

    // Assert
    Vector expected_stress_vector{4};
    expected_stress_vector <<= -7 * 10.0 / 12.0, -10.0, -7 * 10.0 / 12.0, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(actual_stress_vector, expected_stress_vector, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(K0ProcedureIsAppliedCorrectlyWithK0_NCandOCRandNu_UR_3D, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto p_properties = std::make_shared<Properties>();
    p_properties->SetValue(K0_NC, 0.5);
    p_properties->SetValue(K0_MAIN_DIRECTION, 2);
    p_properties->SetValue(OCR, 1.5);
    p_properties->SetValue(POISSON_UNLOADING_RELOADING, 0.25);
    Vector initial_stress_vector{6};
    initial_stress_vector <<= 0.0, -10.0, -10.0, 27.0, 10.0, 5.0;

    // Act
    const auto actual_stress_vector = ApplyK0ProcedureOnStubElement(p_properties, initial_stress_vector);

    // Assert
    Vector expected_stress_vector{6};
    expected_stress_vector <<= -7 * 10.0 / 12.0, -7 * 10.0 / 12.0, -10.0, 0.0, 0.0, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(actual_stress_vector, expected_stress_vector, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(K0ProcedureIsAppliedCorrectlyWithK0_NCandPOP, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto p_properties = std::make_shared<Properties>();
    p_properties->SetValue(K0_NC, 0.5);
    p_properties->SetValue(K0_MAIN_DIRECTION, 1);
    p_properties->SetValue(POP, 50.0);
    Vector initial_stress_vector{4};
    initial_stress_vector <<= 0.0, -10.0, 0.0, 27.0;

    // Act
    const auto actual_stress_vector = ApplyK0ProcedureOnStubElement(p_properties, initial_stress_vector);

    // Assert
    Vector expected_stress_vector{4};
    expected_stress_vector <<= -30.0, -10.0, -30.0, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(actual_stress_vector, expected_stress_vector, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(K0ProcedureIsAppliedCorrectlyWithK0_NCandPOP_3D, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto p_properties = std::make_shared<Properties>();
    p_properties->SetValue(K0_NC, 0.5);
    p_properties->SetValue(K0_MAIN_DIRECTION, 2);
    p_properties->SetValue(POP, 50.0);
    Vector initial_stress_vector{6};
    initial_stress_vector <<= 0.0, -10.0, -10.0, 27.0, 10.0, 5.0;

    // Act
    const auto actual_stress_vector = ApplyK0ProcedureOnStubElement(p_properties, initial_stress_vector);

    // Assert
    Vector expected_stress_vector{6};
    expected_stress_vector <<= -30.0, -30.0, -10.0, 0.0, 0.0, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(actual_stress_vector, expected_stress_vector, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(K0ProcedureIsAppliedCorrectlyWithK0_NCandPOPandNu_UR, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto p_properties = std::make_shared<Properties>();
    p_properties->SetValue(K0_NC, 0.5);
    p_properties->SetValue(K0_MAIN_DIRECTION, 1);
    p_properties->SetValue(POP, 50.0);
    p_properties->SetValue(POISSON_UNLOADING_RELOADING, 0.25);
    Vector initial_stress_vector{4};
    initial_stress_vector <<= 0.0, -10.0, 0.0, 27.0;

    // Act
    const auto actual_stress_vector = ApplyK0ProcedureOnStubElement(p_properties, initial_stress_vector);

    // Assert
    Vector expected_stress_vector{4};
    expected_stress_vector <<= -30.0 + 50.0 / 3.0, -10.0, -30.0 + 50.0 / 3.0, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(actual_stress_vector, expected_stress_vector, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(K0ProcedureIsAppliedCorrectlyWithK0_NCandPOPandNu_UR_3D, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto p_properties = std::make_shared<Properties>();
    p_properties->SetValue(K0_NC, 0.5);
    p_properties->SetValue(K0_MAIN_DIRECTION, 2);
    p_properties->SetValue(POP, 50.0);
    p_properties->SetValue(POISSON_UNLOADING_RELOADING, 0.25);
    Vector initial_stress_vector{6};
    initial_stress_vector <<= 0.0, -10.0, -10.0, 27.0, 10.0, 5.0;

    // Act
    const auto actual_stress_vector = ApplyK0ProcedureOnStubElement(p_properties, initial_stress_vector);

    // Assert
    Vector expected_stress_vector{6};
    expected_stress_vector <<= -30.0 + 50.0 / 3.0, -30.0 + 50.0 / 3.0, -10.0, 0.0, 0.0, 0.0;
    KRATOS_EXPECT_VECTOR_NEAR(actual_stress_vector, expected_stress_vector, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(K0ProcedureChecksIfProcessHasCorrectMaterialData, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto& r_model_part = model.CreateModelPart("main");
    auto  p_element    = make_intrusive<StubElement>();
    p_element->SetId(1);
    p_element->SetProperties(std::make_shared<Properties>());
    r_model_part.AddElement(p_element);

    const auto              k0_settings = Parameters{};
    ApplyK0ProcedureProcess process{r_model_part, k0_settings};

    auto mock_constitutive_law = std::make_shared<MockConstitutiveLaw>();
    p_element->GetProperties().SetValue(CONSTITUTIVE_LAW, mock_constitutive_law);
    EXPECT_CALL(*mock_constitutive_law, WorkingSpaceDimension()).WillRepeatedly(testing::Return(2));

    // Act & Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(process.Check(),
                                      "K0_MAIN_DIRECTION is not defined for element 1.");

    p_element->GetProperties().SetValue(K0_MAIN_DIRECTION, 4);

    EXPECT_CALL(*mock_constitutive_law, WorkingSpaceDimension()).WillRepeatedly(testing::Return(1));
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(process.Check(), "dimension should be 2 or 3 for element 1.")

    EXPECT_CALL(*mock_constitutive_law, WorkingSpaceDimension()).WillRepeatedly(testing::Return(2));
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(process.Check(),
                                      "K0_MAIN_DIRECTION should be 0 or 1 for element 1.")

    EXPECT_CALL(*mock_constitutive_law, WorkingSpaceDimension()).WillRepeatedly(testing::Return(3));
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(process.Check(),
                                      "K0_MAIN_DIRECTION should be 0, 1 or 2 for element 1.")

    p_element->GetProperties().SetValue(K0_MAIN_DIRECTION, 1);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        process.Check(),
        "Insufficient material data for K0 procedure process for element 1. No K0_NC, "
        "(INDEX_OF_UMAT_PHI_PARAMETER and UMAT_PARAMETERS) or (K0_VALUE_XX, _YY and _ZZ found).")
    p_element->GetProperties().SetValue(K0_VALUE_XX, -0.5);
    p_element->GetProperties().SetValue(K0_VALUE_YY, -0.5);
    p_element->GetProperties().SetValue(K0_VALUE_ZZ, -0.5);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        process.Check(), "K0_VALUE_XX (-0.5) should be in the range [0.0,-> for element 1.")
    p_element->GetProperties().SetValue(K0_VALUE_XX, 0.5);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        process.Check(), "K0_VALUE_YY (-0.5) should be in the range [0.0,-> for element 1.")
    p_element->GetProperties().SetValue(K0_VALUE_YY, 0.5);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        process.Check(), "K0_VALUE_ZZ (-0.5) should be in the range [0.0,-> for element 1.")
    p_element->GetProperties().SetValue(K0_VALUE_ZZ, 0.5);

    p_element->GetProperties().SetValue(POISSON_UNLOADING_RELOADING, 0.75);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        process.Check(),
        "POISSON_UNLOADING_RELOADING (0.75) is not in range [-1.0, 0.5> for element 1.")
    p_element->GetProperties().SetValue(POISSON_UNLOADING_RELOADING, 0.25);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        process.Check(), "Insufficient material data for K0 procedure process for "
                         "element 1. Poisson unloading-reloading, OCR and POP functionality cannot "
                         "be combined with K0_VALUE_XX, _YY and _ZZ.")

    p_element->GetProperties().Erase(K0_VALUE_XX);
    p_element->GetProperties().Erase(K0_VALUE_YY);
    p_element->GetProperties().Erase(K0_VALUE_ZZ);

    p_element->GetProperties().SetValue(INDEX_OF_UMAT_PHI_PARAMETER, 2);
    Vector umat_parameters{1};
    umat_parameters[0] = -30.0;
    p_element->GetProperties().SetValue(UMAT_PARAMETERS, umat_parameters);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        process.Check(),
        "INDEX_OF_UMAT_PHI_PARAMETER (2) is not in range 1, size of UMAT_PARAMETERS for element 1")
    p_element->GetProperties().SetValue(INDEX_OF_UMAT_PHI_PARAMETER, 1);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(process.Check(),
                                      "Phi (-30) should be between 0 and 90 degrees for element 1.")

    umat_parameters[0] = 30.0;
    p_element->GetProperties().SetValue(UMAT_PARAMETERS, umat_parameters);
    p_element->GetProperties().SetValue(OCR, 0.5);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(process.Check(),
                                      "OCR (0.5) should be in the range [1.0,-> for element 1.")

    p_element->GetProperties().Erase(OCR);
    p_element->GetProperties().SetValue(POP, -100.0);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(process.Check(),
                                      "POP (-100) should be in the range [0.0,-> for element 1.")

    p_element->GetProperties().Erase(POP);
    p_element->GetProperties().Erase(INDEX_OF_UMAT_PHI_PARAMETER);
    p_element->GetProperties().Erase(UMAT_PARAMETERS);
    p_element->GetProperties().SetValue(K0_NC, -0.5);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(process.Check(),
                                      "K0_NC (-0.5) should be in the range [0.0,-> for element 1.")

    p_element->GetProperties().SetValue(K0_NC, 0.5);
    KRATOS_EXPECT_EQ(process.Check(), 0);
}

KRATOS_TEST_CASE_IN_SUITE(K0ProcedureChecksIfModelPartHasElements, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    auto& r_modelpart = model.CreateModelPart("dummy");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN((ApplyK0ProcedureProcess{r_modelpart, {}}.Check()),
                                      "ApplyK0ProcedureProces has no elements in modelpart dummy")
}

} // namespace Kratos::Testing
