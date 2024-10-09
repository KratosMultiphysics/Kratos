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
#include "geo_mechanics_fast_suite.h"
#include "geometries/quadrilateral_2d_4.h"
#include "processes/structured_mesh_generator_process.h"
#include "stub_linear_elastic_law.h"
#include "test_utilities.h"

#include <boost/algorithm/cxx11/all_of.hpp>
#include <boost/algorithm/cxx11/none_of.hpp>
#include <boost/numeric/ublas/assignment.hpp>

namespace
{

using namespace Kratos;

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
    auto& result = rModel.CreateModelPart("dummy");

    // Set up the test model part mesh
    auto                   p_point_1 = Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0);
    auto                   p_point_2 = Kratos::make_intrusive<Node>(2, 0.0, 1.0, 0.0);
    auto                   p_point_3 = Kratos::make_intrusive<Node>(3, 1.0, 1.0, 0.0);
    auto                   p_point_4 = Kratos::make_intrusive<Node>(4, 1.0, 0.0, 0.0);
    Quadrilateral2D4<Node> domain_geometry(p_point_1, p_point_2, p_point_3, p_point_4);

    Parameters mesher_parameters(R"({
        "number_of_divisions": 2,
        "element_name": "Element2D3N",
        "condition_name": "LineCondition",
        "create_skin_sub_model_part": true
    })");
    StructuredMeshGeneratorProcess(domain_geometry, result, mesher_parameters).Execute();

    auto  p_dummy_law             = std::make_shared<Testing::StubLinearElasticLaw>();
    auto& r_model_part_properties = result.GetProperties(0);
    r_model_part_properties.SetValue(CONSTITUTIVE_LAW, p_dummy_law);

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
                          KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = PrepareTestModelPart(model);

    Parameters k0_settings; // 'use_standard_procedure' is not defined, assume it to be true

    ApplyK0ProcedureProcess process{r_model_part, k0_settings};
    process.ExecuteInitialize();

    KRATOS_EXPECT_TRUE(boost::algorithm::all_of(r_model_part.Elements(), ElementConsidersDiagonalEntriesOnlyAndNoShear))
}

KRATOS_TEST_CASE_IN_SUITE(AllElementsConsiderDiagonalEntriesOnlyAndNoShearWhenUsingStandardProcedure,
                          KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = PrepareTestModelPart(model);

    Parameters k0_settings{R"({"use_standard_procedure": true})"};

    ApplyK0ProcedureProcess process{r_model_part, k0_settings};
    process.ExecuteInitialize();

    KRATOS_EXPECT_TRUE(boost::algorithm::all_of(r_model_part.Elements(), ElementConsidersDiagonalEntriesOnlyAndNoShear))
}

KRATOS_TEST_CASE_IN_SUITE(NoneOfElementsConsiderDiagonalEntriesOnlyAndNoShearWhenNotUsingStandardProcedure,
                          KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = PrepareTestModelPart(model);

    Parameters k0_settings{R"({"use_standard_procedure": false})"};

    ApplyK0ProcedureProcess process{r_model_part, k0_settings};
    process.ExecuteInitialize();

    KRATOS_EXPECT_TRUE(boost::algorithm::none_of(r_model_part.Elements(), ElementConsidersDiagonalEntriesOnlyAndNoShear))
}

KRATOS_TEST_CASE_IN_SUITE(UseStandardProcedureFlagIsInEffectDuringProcessExecutionOnly, KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = PrepareTestModelPart(model);

    Parameters k0_settings{R"({"use_standard_procedure": true})"};

    ApplyK0ProcedureProcess process{r_model_part, k0_settings};
    process.ExecuteInitialize(); // start considering diagonal entries only and no shear
    process.ExecuteFinalize();   // stop considering diagonal entries only and no shear

    KRATOS_EXPECT_TRUE(boost::algorithm::none_of(r_model_part.Elements(), ElementConsidersDiagonalEntriesOnlyAndNoShear))
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

KRATOS_TEST_CASE_IN_SUITE(K0ProcedureIsAppliedCorrectlyWithPhi, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto p_properties = std::make_shared<Properties>();
    p_properties->SetValue(INDEX_OF_UMAT_PHI_PARAMETER, 1);
    p_properties->SetValue(NUMBER_OF_UMAT_PARAMETERS, 1);
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

    // Act & Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(process.Check(),
                                      "K0_MAIN_DIRECTION is not defined for element 1.");

    p_element->GetProperties().SetValue(K0_MAIN_DIRECTION, 4);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(process.Check(),
                                      "K0_MAIN_DIRECTION should be 0 or 1 for element 1.")

    p_element->GetProperties().SetValue(K0_MAIN_DIRECTION, 1);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        process.Check(),
        "Insufficient material data for K0 procedure process for element 1. No K0_NC, "
        "(INDEX_OF_UMAT_PHI_PARAMETER, NUMBER_OF_UMAT_PARAMETERS and "
        "UMAT_PARAMETERS) or (K0_VALUE_XX, _YY and _ZZ found).");
    p_element->GetProperties().SetValue(K0_VALUE_XX, 0.5);
    p_element->GetProperties().SetValue(K0_VALUE_YY, 0.5);
    p_element->GetProperties().SetValue(K0_VALUE_ZZ, 0.5);

    p_element->GetProperties().SetValue(POISSON_UNLOADING_RELOADING, 0.75);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        process.Check(),
        "POISSON_UNLOADING_RELOADING (0.75) is not in range [-1.0, 0.5> for element 1.");
    p_element->GetProperties().SetValue(POISSON_UNLOADING_RELOADING, 0.25);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        process.Check(), "Insufficient material data for K0 procedure process for "
                         "element 1. Poisson unloading-reloading, OCR and POP functionality cannot "
                         "be combined with K0_VALUE_XX, _YY and _ZZ.");

    p_element->GetProperties().Erase(K0_VALUE_XX);
    p_element->GetProperties().Erase(K0_VALUE_YY);
    p_element->GetProperties().Erase(K0_VALUE_ZZ);

    p_element->GetProperties().SetValue(INDEX_OF_UMAT_PHI_PARAMETER, 1);
    p_element->GetProperties().SetValue(NUMBER_OF_UMAT_PARAMETERS, 0);
    Vector umat_parameters{1};
    umat_parameters[0] = -30.0;
    p_element->GetProperties().SetValue(UMAT_PARAMETERS, umat_parameters);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(process.Check(),
                                      "INDEX_OF_UMAT_PHI_PARAMETER (1) is not in range 1, "
                                      "NUMBER_OF_UMAT_PARAMETERS (0) for element 1.");

    p_element->GetProperties().SetValue(NUMBER_OF_UMAT_PARAMETERS, 1);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(process.Check(),
                                      "Phi should be between 0 and 90 degrees for element 1.");

    umat_parameters[0] = 30.0;
    p_element->GetProperties().SetValue(UMAT_PARAMETERS, umat_parameters);
    KRATOS_EXPECT_EQ(process.Check(), 0);
}

} // namespace Kratos::Testing
