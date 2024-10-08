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
#include "containers/variable.h"
#include "custom_processes/reset_displacement_process.h"
#include "geometries/triangle_2d_3.h"
#include "includes/constitutive_law.h"
#include "includes/element.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

#include <boost/numeric/ublas/assignment.hpp>

namespace
{

using namespace Kratos;

class StubConstitutiveLaw : public ConstitutiveLaw
{
public:
    SizeType GetStrainSize() const override { return 4; }
};

class StubElement : public Element
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(StubElement);

    StubElement(IndexType NewId, const GeometryType::Pointer& pGeometry) : Element(NewId, pGeometry)
    {
        mConstitutiveLaws = std::vector<ConstitutiveLaw::Pointer>(3, make_shared<StubConstitutiveLaw>());
    }

    void CalculateOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>&,
                                      std::vector<ConstitutiveLaw::Pointer>& rOutput,
                                      const ProcessInfo&) override
    {
        rOutput = mConstitutiveLaws;
    }

    void SetValuesOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>&,
                                      const std::vector<ConstitutiveLaw::Pointer>& rValues,
                                      const ProcessInfo&) override
    {
        mConstitutiveLaws = rValues;
    }

    void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                      std::vector<Vector>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override
    {
        if (rVariable == PK2_STRESS_VECTOR) rOutput = mIntegrationPointVectors;
    }

    using Element::CalculateOnIntegrationPoints;

    void SetValuesOnIntegrationPoints(const Variable<Vector>&    rVariable,
                                      const std::vector<Vector>& rValues,
                                      const ProcessInfo&) override
    {
        if (rVariable == PK2_STRESS_VECTOR) mIntegrationPointVectors = rValues;
    }

    using Element::SetValuesOnIntegrationPoints;

private:
    std::vector<Vector>                   mIntegrationPointVectors = {};
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLaws        = {};
};

ModelPart& CreateModelPartWithAStubElement(Model& rModel)
{
    auto& model_part = rModel.CreateModelPart("MainModelPart");
    auto  node_1     = model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto  node_2     = model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    auto  node_3     = model_part.CreateNewNode(3, 1.0, 1.0, 0.0);

    auto geometry = std::make_shared<Triangle2D3<Node>>(node_1, node_2, node_3);
    model_part.AddElement(make_intrusive<StubElement>(1, geometry));

    return model_part;
}
} // namespace

namespace Kratos::Testing
{
KRATOS_TEST_CASE_IN_SUITE(ResetDisplacementProcess_SetsInitialStressOfConstitutiveLaws,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model  model;
    auto&  model_part = CreateModelPartWithAStubElement(model);
    Vector initial_stress_vector(4);
    initial_stress_vector <<= 1.0, 2.0, 3.0, 4.0;

    constexpr auto number_of_integration_points = 3;
    const auto     dummy_process_info           = ProcessInfo{};
    model_part.Elements()[1].SetValuesOnIntegrationPoints(
        PK2_STRESS_VECTOR, std::vector<Vector>(number_of_integration_points, initial_stress_vector),
        dummy_process_info);

    const auto               dummy_parameters = Parameters{};
    ResetDisplacementProcess reset_displacement_process(model_part, dummy_parameters);
    reset_displacement_process.ExecuteInitialize();

    std::vector<ConstitutiveLaw::Pointer> constitutive_laws;
    model_part.Elements()[1].CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, constitutive_laws,
                                                          model_part.GetProcessInfo());

    KRATOS_EXPECT_EQ(constitutive_laws.size(), number_of_integration_points);
    for (const auto& constitutive_law : constitutive_laws) {
        KRATOS_EXPECT_VECTOR_NEAR(constitutive_law->GetInitialState().GetInitialStressVector(),
                                  initial_stress_vector, 1e-12)
        KRATOS_EXPECT_VECTOR_NEAR(constitutive_law->GetInitialState().GetInitialStrainVector(),
                                  Vector{ZeroVector{4}}, 1e-12)
    }
}

KRATOS_TEST_CASE_IN_SUITE(ResetDisplacementProcess_ThrowsInCheck_WhenModelIsNotRestarted,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    auto& model_part                          = model.CreateModelPart("Main");
    model_part.GetProcessInfo()[IS_RESTARTED] = false;

    const auto               dummy_parameters = Parameters{};
    ResetDisplacementProcess process(model_part, dummy_parameters);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        process.Check(), "The IS_RESTARTED flag must be set to true in the ProcessInfo of the "
                         "model part. Please use the \"rest\" option for the model input type")

    model_part.GetProcessInfo()[IS_RESTARTED] = true;
    KRATOS_EXPECT_EQ(process.Check(), 0);
}

KRATOS_TEST_CASE_IN_SUITE(ResetDisplacementProcess_ThrowsInExecuteInitialize_WhenConstitutiveLawsCannotBeRetrieved,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    auto& model_part = CreateModelPartWithAStubElement(model);

    const auto               dummy_parameters = Parameters{};
    ResetDisplacementProcess process(model_part, dummy_parameters);

    const auto dummy_process_info          = ProcessInfo{};
    const auto empty_constitutive_law_list = std::vector<ConstitutiveLaw::Pointer>{};
    model_part.Elements()[1].SetValuesOnIntegrationPoints(
        CONSTITUTIVE_LAW, empty_constitutive_law_list, dummy_process_info);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        process.ExecuteInitialize(),
        "The constitutive laws on the integration points could not be retrieved for element 1")
}

KRATOS_TEST_CASE_IN_SUITE(ResetDisplacementProcess_ThrowsInExecuteInitialize_WhenStressVectorsCannotBeRetrieved,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    auto& model_part = CreateModelPartWithAStubElement(model);

    const auto               dummy_parameters = Parameters{};
    ResetDisplacementProcess process(model_part, dummy_parameters);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        process.ExecuteInitialize(),
        "The stress vectors on the integration points could not be retrieved for element 1")
}

KRATOS_TEST_CASE_IN_SUITE(ResetDisplacementProcess_ThrowsInExecuteInitialize_WhenStressVectorsAndConstitutiveLawsAreNotTheSameSize,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model      model;
    auto&      model_part         = CreateModelPartWithAStubElement(model);
    const auto dummy_process_info = ProcessInfo{};
    model_part.Elements()[1].SetValuesOnIntegrationPoints(
        PK2_STRESS_VECTOR, std::vector<Vector>(2, ScalarVector(4, 1.0)), dummy_process_info);

    const auto               dummy_parameters = Parameters{};
    ResetDisplacementProcess process(model_part, dummy_parameters);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(process.ExecuteInitialize(),
                                      "Number of retrieved stress vectors (2) does not match the "
                                      "number of constitutive laws (3) for element 1")
}

} // namespace Kratos::Testing