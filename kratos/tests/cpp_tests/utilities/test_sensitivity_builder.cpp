//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:
//

// External includes
#include <array>
#include <limits>
#include <string>

// Project includes
#include "containers/model.h"
#include "geometries/triangle_2d_3.h"
#include "includes/element.h"
#include "response_functions/adjoint_response_function.h"
#include "testing/testing.h"
#include "utilities/sensitivity_builder.h"
#include "utilities/variable_utils.h"

namespace
{
namespace test_sensitivity_builder
{
using namespace Kratos;
template <class TElement>
class TestAdjoint : public TElement
{
public:
    typedef Kratos::intrusive_ptr<TestAdjoint> Pointer;
    typedef Kratos::unique_ptr<TestAdjoint> UniquePointer;


    static typename TestAdjoint::Pointer Create(std::size_t NewId, const PointerVector<Node<3>>& rNodes)
    {
        Geometry<Node<3>>::Pointer p_geom =
            Kratos::make_shared<Triangle2D3<Node<3>>>(rNodes);
        return Kratos::make_intrusive<TestAdjoint>(NewId, p_geom);
    }

    TestAdjoint(std::size_t NewId, Geometry<Node<3>>::Pointer pGeom)
        : TElement(NewId, pGeom)
    {
    }

    void GetValuesVector(Vector& rValues, int) const override
    {
        rValues.resize(3, false);
        rValues[0] = 12.;
        rValues[1] = -6.;
        rValues[2] = 2.;
    }

    void CalculateSensitivityMatrix(const Variable<double>& rDesignVariable,
                                    Matrix& rOutput,
                                    const ProcessInfo& rProcessInfo) override
    {
        if (rDesignVariable == NORMAL_SENSITIVITY)
        {
            rOutput = ZeroMatrix(3, 3);
            rOutput(0, 0) = rOutput(1, 1) = rOutput(2, 2) = 2.;
        }
        else if (rDesignVariable == SCALAR_SENSITIVITY)
        {
            rOutput = ZeroMatrix(1, 3);
            rOutput(0, 0) = 2.;
        }
        else
        {
            KRATOS_ERROR << "Invalid variable: " << rDesignVariable << std::endl;
        }
    }

    void CalculateSensitivityMatrix(const Variable<array_1d<double, 3>>& rDesignVariable,
                                    Matrix& rOutput,
                                    const ProcessInfo& rProcessInfo) override
    {
        if (rDesignVariable == SHAPE_SENSITIVITY)
        {
            rOutput = ZeroMatrix(6, 3);
            rOutput(0, 0) = rOutput(1, 1) = rOutput(2, 2) = 2.;
            rOutput(3, 0) = rOutput(4, 1) = rOutput(5, 2) = 2.;
        }
        else
        {
            KRATOS_ERROR << "Invalid variable: " << rDesignVariable << std::endl;
        }
    }
};

using TestAdjointElement = TestAdjoint<Element>;
using TestAdjointCondition = TestAdjoint<Condition>;

class TestResponseFunction : public AdjointResponseFunction
{
public:
    void CalculatePartialSensitivity(Element& rAdjointElement,
                                     const Variable<double>& rVariable,
                                     const Matrix& rSensitivityMatrix,
                                     Vector& rSensitivityGradient,
                                     const ProcessInfo& rProcessInfo) override
    {
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        rSensitivityGradient[0] = 1.;
    }

    void CalculatePartialSensitivity(Condition& rAdjointCondition,
                                     const Variable<double>& rVariable,
                                     const Matrix& rSensitivityMatrix,
                                     Vector& rSensitivityGradient,
                                     const ProcessInfo& rProcessInfo) override
    {
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
    }

    void CalculatePartialSensitivity(Element& rAdjointElement,
                                     const Variable<array_1d<double, 3>>& rVariable,
                                     const Matrix& rSensitivityMatrix,
                                     Vector& rSensitivityGradient,
                                     const ProcessInfo& rProcessInfo) override
    {
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        rSensitivityGradient[0] = 1.;
    }

    void CalculatePartialSensitivity(Condition& rAdjointCondition,
                                     const Variable<array_1d<double, 3>>& rVariable,
                                     const Matrix& rSensitivityMatrix,
                                     Vector& rSensitivityGradient,
                                     const ProcessInfo& rProcessInfo) override
    {
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
    }

    double CalculateValue(ModelPart& rModelPart) override
    {
        return 0.;
    };
};

PointerVector<Node<3>> GetNodes(ModelPart& rModelPart, const std::vector<unsigned>& rIds)
{
    PointerVector<Node<3>> nodes;
    for (auto id : rIds)
    {
        nodes.push_back(rModelPart.pGetNode(id));
    }
    return nodes;
};

ModelPart& CreateModelPartWithTestElements(Model& rModel)
{
    auto& r_model_part = rModel.CreateModelPart("test");
    r_model_part.AddNodalSolutionStepVariable(NORMAL_SENSITIVITY);
    r_model_part.AddNodalSolutionStepVariable(SHAPE_SENSITIVITY);
    r_model_part.AddNodalSolutionStepVariable(SCALAR_SENSITIVITY);
    r_model_part.CreateNewNode(1, 1., 0., 0.);
    r_model_part.CreateNewNode(2, 3., 0., 0.);
    r_model_part.CreateNewNode(3, 5., 0., 0.);
    r_model_part.CreateNewNode(4, 2., 1., 0.);
    r_model_part.CreateNewNode(5, 4., 1., 0.);
    r_model_part.CreateNewNode(6, 3., 2., 0.);
    r_model_part.AddElement(
        TestAdjointElement::Create(1, GetNodes(r_model_part, {1, 2, 4})));
    r_model_part.AddElement(
        TestAdjointElement::Create(2, GetNodes(r_model_part, {2, 5, 4})));
    r_model_part.AddElement(
        TestAdjointElement::Create(3, GetNodes(r_model_part, {2, 3, 5})));
    r_model_part.AddElement(
        TestAdjointElement::Create(4, GetNodes(r_model_part, {4, 5, 6})));
    r_model_part.SetBufferSize(1);
    VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, true,
                                             r_model_part.Nodes());
    VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, true,
                                             r_model_part.Elements());
    return r_model_part;
}

ModelPart& CreateModelPartWithTestConditions(Model& rModel)
{
    auto& r_model_part = rModel.CreateModelPart("test");
    r_model_part.AddNodalSolutionStepVariable(NORMAL_SENSITIVITY);
    r_model_part.AddNodalSolutionStepVariable(SHAPE_SENSITIVITY);
    r_model_part.AddNodalSolutionStepVariable(SCALAR_SENSITIVITY);
    r_model_part.CreateNewNode(1, 1., 0., 0.);
    r_model_part.CreateNewNode(2, 3., 0., 0.);
    r_model_part.CreateNewNode(3, 5., 0., 0.);
    r_model_part.CreateNewNode(4, 2., 1., 0.);
    r_model_part.CreateNewNode(5, 4., 1., 0.);
    r_model_part.CreateNewNode(6, 3., 2., 0.);
    r_model_part.AddCondition(
        TestAdjointCondition::Create(1, GetNodes(r_model_part, {1, 2, 4})));
    r_model_part.AddCondition(
        TestAdjointCondition::Create(2, GetNodes(r_model_part, {2, 5, 4})));
    r_model_part.AddCondition(
        TestAdjointCondition::Create(3, GetNodes(r_model_part, {2, 3, 5})));
    r_model_part.AddCondition(
        TestAdjointCondition::Create(4, GetNodes(r_model_part, {4, 5, 6})));
    r_model_part.SetBufferSize(1);
    VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, true,
                                             r_model_part.Nodes());
    VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, true,
                                             r_model_part.Conditions());
    return r_model_part;
}

// Test number magnitudes are approx. 1-150 -> choosing tolerance approx. 10 * <magnitude> * epsilon
constexpr double tolerance = std::numeric_limits<double>::epsilon() * 1500.;

} // namespace test_sensitivity_builder
} // namespace

namespace Kratos
{
namespace Testing
{
KRATOS_TEST_CASE_IN_SUITE(SensitivityBuilder_CalculateNodalSolutionStepSensitivities_Double,
                          KratosSensitivityTestSuite)
{
    using namespace test_sensitivity_builder;
    Model model;
    auto& model_part = CreateModelPartWithTestElements(model);
    TestResponseFunction response_function;
    SensitivityBuilder::CalculateNodalSolutionStepSensitivities(
        {"NORMAL_SENSITIVITY"}, model_part, response_function, 4.);
    auto normal_sensitivity = [&model_part](std::size_t i) -> double {
        return model_part.GetNode(i).FastGetSolutionStepValue(NORMAL_SENSITIVITY);
    };
    KRATOS_CHECK_NEAR(normal_sensitivity(1), 100., tolerance);
    KRATOS_CHECK_NEAR(normal_sensitivity(2), 152., tolerance);
    KRATOS_CHECK_NEAR(normal_sensitivity(3), -48., tolerance);
    KRATOS_CHECK_NEAR(normal_sensitivity(4), 132., tolerance);
    KRATOS_CHECK_NEAR(normal_sensitivity(5), -80., tolerance);
    KRATOS_CHECK_NEAR(normal_sensitivity(6), 16., tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(SensitivityBuilder_CalculateNodalSolutionStepSensitivities_Array1d,
                          KratosSensitivityTestSuite)
{
    using namespace test_sensitivity_builder;
    Model model;
    auto& model_part = CreateModelPartWithTestElements(model);
    TestResponseFunction response_function;
    SensitivityBuilder::CalculateNodalSolutionStepSensitivities(
        {"SHAPE_SENSITIVITY"}, model_part, response_function, 4.);
    auto CheckNode = [](Node<3>& rNode, std::array<double, 2> RefValue) {
        const double val_x = rNode.FastGetSolutionStepValue(SHAPE_SENSITIVITY_X);
        const double val_y = rNode.FastGetSolutionStepValue(SHAPE_SENSITIVITY_Y);
        KRATOS_CHECK_NEAR(val_x, RefValue[0], tolerance);
        KRATOS_CHECK_NEAR(val_y, RefValue[1], tolerance);
    };
    CheckNode(model_part.GetNode(1), {{100., -48.}});
    CheckNode(model_part.GetNode(2), {{216., 0.}});
    CheckNode(model_part.GetNode(3), {{16., 96.}});
    CheckNode(model_part.GetNode(4), {{4., -16.}});
    CheckNode(model_part.GetNode(5), {{-16., 208.}});
    CheckNode(model_part.GetNode(6), {{-48., 16.}});
}

KRATOS_TEST_CASE_IN_SUITE(SensitivityBuilder_CalculateNodalSolutionStepSensitivities_InactiveNodes,
                          KratosSensitivityTestSuite)
{
    using namespace test_sensitivity_builder;
    Model model;
    auto& model_part = CreateModelPartWithTestElements(model);
    for (auto& r_node : model_part.GetElement(4).GetGeometry())
    {
        r_node.SetValue(UPDATE_SENSITIVITIES, false);
    }
    TestResponseFunction response_function;
    SensitivityBuilder::CalculateNodalSolutionStepSensitivities(
        {"NORMAL_SENSITIVITY"}, model_part, response_function, 4.);
    auto normal_sensitivity = [&model_part](std::size_t i) -> double {
        return model_part.GetNode(i).FastGetSolutionStepValue(NORMAL_SENSITIVITY);
    };
    KRATOS_CHECK_NEAR(normal_sensitivity(1), 100., tolerance);
    KRATOS_CHECK_NEAR(normal_sensitivity(2), 152., tolerance);
    KRATOS_CHECK_NEAR(normal_sensitivity(3), -48., tolerance);
    KRATOS_CHECK_NEAR(normal_sensitivity(4), 0., tolerance);
    KRATOS_CHECK_NEAR(normal_sensitivity(5), 0., tolerance);
    KRATOS_CHECK_NEAR(normal_sensitivity(6), 0., tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(SensitivityBuilder_CalculateNodalSolutionStepSensitivities_Unsupported,
                          KratosSensitivityTestSuite)
{
    using namespace test_sensitivity_builder;
    Model model;
    auto& model_part = CreateModelPartWithTestElements(model);
    TestResponseFunction response_function;
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        SensitivityBuilder::CalculateNodalSolutionStepSensitivities(
            {"UPDATE_SENSITIVITIES"}, model_part, response_function, 4.);
        , "Unsupported variable: UPDATE_SENSITIVITIES");
}

KRATOS_TEST_CASE_IN_SUITE(SensitivityBuilder_CalculateNodalSolutionStepSensitivities_Condition,
                          KratosSensitivityTestSuite)
{
    using namespace test_sensitivity_builder;
    Model model;
    auto& model_part = CreateModelPartWithTestConditions(model);
    TestResponseFunction response_function;
    SensitivityBuilder::CalculateNodalSolutionStepSensitivities(
        {"NORMAL_SENSITIVITY"}, model_part, response_function, 4.);
    auto normal_sensitivity = [&model_part](std::size_t i) -> double {
        return model_part.GetNode(i).FastGetSolutionStepValue(NORMAL_SENSITIVITY);
    };
    KRATOS_CHECK_NEAR(normal_sensitivity(1), 96., tolerance);
    KRATOS_CHECK_NEAR(normal_sensitivity(2), 144., tolerance);
    KRATOS_CHECK_NEAR(normal_sensitivity(3), -48., tolerance);
    KRATOS_CHECK_NEAR(normal_sensitivity(4), 128., tolerance);
    KRATOS_CHECK_NEAR(normal_sensitivity(5), -80., tolerance);
    KRATOS_CHECK_NEAR(normal_sensitivity(6), 16., tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(SensitivityBuilder_CalculateNonHistoricalSensitivities_Double,
                          KratosSensitivityTestSuite)
{
    using namespace test_sensitivity_builder;
    Model model;
    auto& model_part = CreateModelPartWithTestElements(model);
    VariableUtils().SetNonHistoricalVariable(
        SCALAR_SENSITIVITY, SCALAR_SENSITIVITY.Zero(), model_part.Elements());
    model_part.GetElement(4).SetValue(UPDATE_SENSITIVITIES, false);
    TestResponseFunction response_function;
    SensitivityBuilder::CalculateNonHistoricalSensitivities(
        {"SCALAR_SENSITIVITY"}, model_part.Elements(), response_function,
        model_part.GetProcessInfo(), 4.);
    KRATOS_CHECK_NEAR(model_part.GetElement(1).GetValue(SCALAR_SENSITIVITY), 100., tolerance);
    KRATOS_CHECK_NEAR(model_part.GetElement(2).GetValue(SCALAR_SENSITIVITY), 100., tolerance);
    KRATOS_CHECK_NEAR(model_part.GetElement(3).GetValue(SCALAR_SENSITIVITY), 100., tolerance);
    KRATOS_CHECK_NEAR(model_part.GetElement(4).GetValue(SCALAR_SENSITIVITY), 0., tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(SensitivityBuilder_CalculateNonHistoricalSensitivities_Unsupported,
                          KratosSensitivityTestSuite)
{
    using namespace test_sensitivity_builder;
    Model model;
    auto& model_part = CreateModelPartWithTestElements(model);
    VariableUtils().SetNonHistoricalVariable(
        SCALAR_SENSITIVITY, SCALAR_SENSITIVITY.Zero(), model_part.Elements());
    TestResponseFunction response_function;
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        SensitivityBuilder::CalculateNonHistoricalSensitivities(
            {"UPDATE_SENSITIVITIES"}, model_part.Elements(), response_function,
            model_part.GetProcessInfo(), 4.);
        , "Unsupported variable: UPDATE_SENSITIVITIES");
}

KRATOS_TEST_CASE_IN_SUITE(SensitivityBuilder_CalculateNonHistoricalSensitivities_Condition,
                          KratosSensitivityTestSuite)
{
    using namespace test_sensitivity_builder;
    Model model;
    auto& model_part = CreateModelPartWithTestConditions(model);
    VariableUtils().SetNonHistoricalVariable(
        SCALAR_SENSITIVITY, SCALAR_SENSITIVITY.Zero(), model_part.Conditions());
    model_part.GetCondition(4).SetValue(UPDATE_SENSITIVITIES, false);
    TestResponseFunction response_function;
    SensitivityBuilder::CalculateNonHistoricalSensitivities(
        {"SCALAR_SENSITIVITY"}, model_part.Conditions(), response_function,
        model_part.GetProcessInfo(), 4.);
    KRATOS_CHECK_NEAR(model_part.GetCondition(1).GetValue(SCALAR_SENSITIVITY), 96., tolerance);
    KRATOS_CHECK_NEAR(model_part.GetCondition(2).GetValue(SCALAR_SENSITIVITY), 96., tolerance);
    KRATOS_CHECK_NEAR(model_part.GetCondition(3).GetValue(SCALAR_SENSITIVITY), 96., tolerance);
    KRATOS_CHECK_NEAR(model_part.GetCondition(4).GetValue(SCALAR_SENSITIVITY), 0., tolerance);
}

} // namespace Testing
} // namespace Kratos.