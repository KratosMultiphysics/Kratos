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

#include "custom_processes/reset_displacement_process.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

#include <boost/numeric/ublas/assignment.hpp>

namespace
{

using namespace Kratos;

class StubConstitutiveLaw : public ConstitutiveLaw
{
};

class StubElement : public Element
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(StubElement);

    StubElement(IndexType NewId, const GeometryType::Pointer& pGeometry) : Element(NewId, pGeometry)
    {
        mConstitutiveLaws = std::vector<ConstitutiveLaw::Pointer>(3, make_shared<StubConstitutiveLaw>());
    }

    void CalculateOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>& rVariable,
                                      std::vector<ConstitutiveLaw::Pointer>&    rOutput,
                                      const ProcessInfo& rCurrentProcessInfo) override
    {
        rOutput = mConstitutiveLaws;
    }

    void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                      std::vector<Vector>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override
    {
        rOutput = mIntegrationPointVectors;
    }

    using Element::CalculateOnIntegrationPoints;

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
    auto element  = Kratos::make_intrusive<StubElement>(1, geometry);
    model_part.AddElement(element);

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
    dynamic_cast<StubElement&>(model_part.Elements()[1]).mIntegrationPointVectors =
        std::vector<Vector>(number_of_integration_points, initial_stress_vector);

    ResetDisplacementProcess reset_displacement_process(model_part, {});
    reset_displacement_process.ExecuteInitialize();

    std::vector<ConstitutiveLaw::Pointer> constitutive_laws;
    model_part.Elements()[1].CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, constitutive_laws,
                                                          model_part.GetProcessInfo());

    KRATOS_EXPECT_EQ(constitutive_laws.size(), number_of_integration_points);
    for (const auto& constitutive_law : constitutive_laws) {
        KRATOS_EXPECT_VECTOR_NEAR(constitutive_law->GetInitialState().GetInitialStressVector(),
                                  initial_stress_vector, 1e-12)
    }
}

KRATOS_TEST_CASE_IN_SUITE(ResetDisplacementProcess_ThrowsInCheck_WhenModelIsNotRestarted,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    auto& model_part                          = model.CreateModelPart("Main");
    model_part.GetProcessInfo()[IS_RESTARTED] = false;

    ResetDisplacementProcess process(model_part, {});

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

    ResetDisplacementProcess process(model_part, {});
    dynamic_cast<StubElement&>(model_part.Elements()[1]).mConstitutiveLaws.clear();
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        process.ExecuteInitialize(),
        "The constitutive laws on the integration points could not be retrieved for element 1")
}

KRATOS_TEST_CASE_IN_SUITE(ResetDisplacementProcess_ThrowsInExecuteInitialize_WhenStressVectorsCannotBeRetrieved,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    auto& model_part = CreateModelPartWithAStubElement(model);

    ResetDisplacementProcess process(model_part, {});

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        process.ExecuteInitialize(),
        "The stress vectors on the integration points could not be retrieved for element 1")
}

KRATOS_TEST_CASE_IN_SUITE(ResetDisplacementProcess_ThrowsInExecuteInitialize_WhenStressVectorsAndConstitutiveLawsAreNotTheSameSize,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    auto& model_part = CreateModelPartWithAStubElement(model);
    dynamic_cast<StubElement&>(model_part.Elements()[1]).mIntegrationPointVectors =
        std::vector<Vector>(2, ScalarVector(4, 1.0));

    ResetDisplacementProcess process(model_part, {});

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(process.ExecuteInitialize(),
                                      "Number of retrieved stress vectors does not match the "
                                      "number of constitutive laws for element 1")
}

} // namespace Kratos::Testing