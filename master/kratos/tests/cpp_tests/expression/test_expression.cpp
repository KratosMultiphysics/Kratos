//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <algorithm>
#include <limits>
#include <vector>
#include <numeric>

// External includes

// Project includes
#include "containers/model.h"
#include "expression/expression_io_utils.h"
#include "expression/integration_point_expression_io.h"
#include "expression/c_array_expression_io.h"
#include "includes/condition.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/smart_pointers.h"
#include "testing/testing.h"

namespace Kratos::Testing
{

namespace ExpressionIntegrationPointIOTestUtilities
{

template<class TEntityType>
class DummyEntity : public TEntityType {
public:
    ///@name Type Definitions
    ///@{

    using BaseType = TEntityType;

    using IndexType = typename BaseType::IndexType;

    using NodesArrayType = typename BaseType::NodesArrayType;

    using GeometryType = typename BaseType::GeometryType;

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(DummyEntity);

    ///@}
    ///@name Life cycle
    ///@{

    DummyEntity(IndexType NewId)
        : BaseType(NewId)
    {
    }

    ~DummyEntity() override = default;

    ///@}
    ///@name Public operations
    ///@{

    void CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rProcessInfo) override
    {
        if (rOutput.size() != mIPDouble.size()) {
            rOutput.resize(mIPDouble.size());
        }

        std::copy(mIPDouble.begin(), mIPDouble.end(), rOutput.begin());
    }

    void CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 3>>& rVariable,
        std::vector<array_1d<double, 3>>& rOutput,
        const ProcessInfo& rProcessInfo) override
    {
        if (rOutput.size() != mIPArray3.size()) {
            rOutput.resize(mIPArray3.size());
        }

        std::copy(mIPArray3.begin(), mIPArray3.end(), rOutput.begin());
    }

    void CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rOutput,
        const ProcessInfo& rProcessInfo) override
    {
        if (mIPVector.size() == 0) return;

        if (rOutput.size() != mIPVector.size()) {
            rOutput.resize(mIPVector.size(), Vector(mIPVector[0].size()));
        }

        std::copy(mIPVector.begin(), mIPVector.end(), rOutput.begin());
    }

    void CalculateOnIntegrationPoints(
        const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rOutput,
        const ProcessInfo& rProcessInfo) override
    {
        if (mIPMatrix.size() == 0) return;

        if (rOutput.size() != mIPMatrix.size()) {
            rOutput.resize(mIPMatrix.size(), Matrix(mIPMatrix[0].size1(), mIPMatrix[0].size2()));
        }

        std::copy(mIPMatrix.begin(), mIPMatrix.end(), rOutput.begin());
    }

    void SetValuesOnIntegrationPoints(
        const Variable<double>& rVariable,
        const std::vector<double>& rValues,
        const ProcessInfo& rProcessInfo) override
    {
        if (mIPDouble.size() != rValues.size()) {
            mIPDouble.resize(rValues.size());
        }

        std::copy(rValues.begin(), rValues.end(), mIPDouble.begin());
    }

    void SetValuesOnIntegrationPoints(
        const Variable<array_1d<double, 3>>& rVariable,
        const std::vector<array_1d<double, 3>>& rValues,
        const ProcessInfo& rProcessInfo) override
    {
        if (mIPArray3.size() != rValues.size()) {
            mIPArray3.resize(rValues.size());
        }

        std::copy(rValues.begin(), rValues.end(), mIPArray3.begin());
    }

    void SetValuesOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        const std::vector<Vector>& rValues,
        const ProcessInfo& rProcessInfo) override
    {
        if (rValues.size() == 0) return;

        if (mIPVector.size() != rValues.size()) {
            mIPVector.resize(rValues.size(), Vector(rValues[0].size()));
        }

        std::copy(rValues.begin(), rValues.end(), mIPVector.begin());
    }

    void SetValuesOnIntegrationPoints(
        const Variable<Matrix>& rVariable,
        const std::vector<Matrix>& rValues,
        const ProcessInfo& rProcessInfo) override
    {
        if (rValues.size() == 0) return;

        if (mIPMatrix.size() != rValues.size()) {
            mIPMatrix.resize(rValues.size(), Matrix(rValues[0].size1(), rValues[1].size2()));
        }

        std::copy(rValues.begin(), rValues.end(), mIPMatrix.begin());
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    std::vector<double> mIPDouble;

    std::vector<array_1d<double, 3>> mIPArray3;

    std::vector<Vector> mIPVector;

    std::vector<Matrix> mIPMatrix;

    ///@}
};

template<class TContainerType>
void FillData(
    TContainerType& rContainer,
    const ProcessInfo& rProcessInfo)
{
    IndexType number_of_integration_points = 10;
    for (auto& r_element : rContainer) {
        std::vector<double> double_values(10);
        std::vector<array_1d<double, 3>> array_3_values(10);
        std::vector<Vector> vector_values(10, Vector(12));
        std::vector<Matrix> matrix_3_values(10, Matrix(4, 3));

        for (IndexType i = 0; i < number_of_integration_points; ++i) {
            double_values[i] = r_element.Id() + i + 1;
            array_3_values[i][0] = r_element.Id() + i + 1;
            array_3_values[i][1] = r_element.Id() + i + 2;
            array_3_values[i][2] = r_element.Id() + i + 3;
            for (IndexType j = 0; j < 12; ++j) {
                vector_values[i][j] = r_element.Id() + i + j;
                matrix_3_values[i].data()[j] = r_element.Id() + i + j;
            }
        }

        r_element.SetValuesOnIntegrationPoints(PRESSURE, double_values, rProcessInfo);
        r_element.SetValuesOnIntegrationPoints(VELOCITY, array_3_values, rProcessInfo);
        r_element.SetValuesOnIntegrationPoints(INITIAL_STRAIN, vector_values, rProcessInfo);
        r_element.SetValuesOnIntegrationPoints(NORMAL_SHAPE_DERIVATIVE, matrix_3_values, rProcessInfo);
    }
}

ModelPart& SetUpModelPart(Model& rModel)
{
    auto& r_model_part = rModel.CreateModelPart("test");
    auto p_element_1 = Kratos::make_intrusive<DummyEntity<Element>>(1);
    auto p_element_2 = Kratos::make_intrusive<DummyEntity<Element>>(2);
    auto p_condition_1 = Kratos::make_intrusive<DummyEntity<Condition>>(3);
    auto p_condition_2 = Kratos::make_intrusive<DummyEntity<Condition>>(4);

    r_model_part.AddElement(p_element_1);
    r_model_part.AddElement(p_element_2);
    r_model_part.AddCondition(p_condition_1);
    r_model_part.AddCondition(p_condition_2);

    FillData(r_model_part.Conditions(), r_model_part.GetProcessInfo());
    FillData(r_model_part.Elements(), r_model_part.GetProcessInfo());

    return r_model_part;
}

ModelPart& SetUpEmptyModelPart(Model& rModel)
{
    auto& r_model_part = rModel.CreateModelPart("test");
    return r_model_part;
}

template<class TContainerType>
void ExecuteReadTest(
    TContainerType& rContainer,
    ModelPart& rModelPart,
    const IndexType NumberOfIntegrationPoints)
{
    ContainerExpression<TContainerType> double_exp(rModelPart);
    IntegrationPointExpressionIO::Read(double_exp, &PRESSURE);
    if (NumberOfIntegrationPoints > 0) {
        std::vector<std::size_t> double_size = {NumberOfIntegrationPoints};
        KRATOS_EXPECT_EQ(double_exp.GetItemShape(), double_size);
    } else {
        KRATOS_EXPECT_EQ(double_exp.GetItemShape().size(), 0);
    }

    ContainerExpression<TContainerType> array3_exp(rModelPart);
    IntegrationPointExpressionIO::Read(array3_exp, &VELOCITY);
    if (NumberOfIntegrationPoints > 0) {
        std::vector<std::size_t> array3_size = {NumberOfIntegrationPoints, 3};
        KRATOS_EXPECT_EQ(array3_exp.GetItemShape(), array3_size);
    } else {
        KRATOS_EXPECT_EQ(array3_exp.GetItemShape().size(), 0);
    }

    ContainerExpression<TContainerType> vector_exp(rModelPart);
    IntegrationPointExpressionIO::Read(vector_exp, &INITIAL_STRAIN);
    if (NumberOfIntegrationPoints > 0) {
        std::vector<std::size_t> vector_size = {NumberOfIntegrationPoints, 12};
        KRATOS_EXPECT_EQ(vector_exp.GetItemShape(), vector_size);
    } else {
        KRATOS_EXPECT_EQ(vector_exp.GetItemShape().size(), 0);
    }

    ContainerExpression<TContainerType> matrix_exp(rModelPart);
    IntegrationPointExpressionIO::Read(matrix_exp, &NORMAL_SHAPE_DERIVATIVE);
    if (NumberOfIntegrationPoints > 0) {
        std::vector<std::size_t> matrix_size = {NumberOfIntegrationPoints, 4, 3};
        KRATOS_EXPECT_EQ(matrix_exp.GetItemShape(), matrix_size);
    } else {
        KRATOS_EXPECT_EQ(matrix_exp.GetItemShape().size(), 0);
    }

    for (IndexType entity_index = 0; entity_index < rContainer.size(); ++entity_index) {
        const auto& r_entity = *(rContainer.begin() + entity_index);
        const auto entity_id = r_entity.Id();
        for (IndexType i_gauss = 0; i_gauss < NumberOfIntegrationPoints; ++i_gauss) {
            KRATOS_EXPECT_EQ(double_exp.GetExpression().Evaluate(entity_index, entity_index * NumberOfIntegrationPoints, i_gauss), entity_id + i_gauss + 1);

            KRATOS_EXPECT_EQ(array3_exp.GetExpression().Evaluate(entity_index, entity_index * NumberOfIntegrationPoints * 3, i_gauss * 3 + 0), entity_id + i_gauss + 1);
            KRATOS_EXPECT_EQ(array3_exp.GetExpression().Evaluate(entity_index, entity_index * NumberOfIntegrationPoints * 3, i_gauss * 3 + 1), entity_id + i_gauss + 2);
            KRATOS_EXPECT_EQ(array3_exp.GetExpression().Evaluate(entity_index, entity_index * NumberOfIntegrationPoints * 3, i_gauss * 3 + 2), entity_id + i_gauss + 3);

            for (IndexType i_comp = 0; i_comp < 12; ++i_comp) {
                KRATOS_EXPECT_EQ(vector_exp.GetExpression().Evaluate(entity_index, entity_index * NumberOfIntegrationPoints * 12, i_gauss * 12 + i_comp), entity_id + i_gauss + i_comp);
                KRATOS_EXPECT_EQ(matrix_exp.GetExpression().Evaluate(entity_index, entity_index * NumberOfIntegrationPoints * 12, i_gauss * 12 + i_comp), entity_id + i_gauss + i_comp);
            }
        }
    }
}

template<class TContainerType>
void ExecuteWriteTest(
    TContainerType& rContainer,
    ModelPart& rModelPart)
{
    IndexType number_of_entities = rContainer.size();
    IndexType number_of_integration_points = 10;

    std::vector<double> values(number_of_integration_points * number_of_entities * 12);
    std::iota(values.begin(), values.end(), 1);

    ContainerExpression<TContainerType> double_exp(rModelPart);
    std::vector<int> double_shape{10};
    CArrayExpressionIO::Read(double_exp, &values[0], number_of_entities, &double_shape[0], 1);
    IntegrationPointExpressionIO::Write(double_exp, &PRESSURE);

    ContainerExpression<TContainerType> array3_exp(rModelPart);
    std::vector<int> array3_shape{10, 3};
    CArrayExpressionIO::Read(array3_exp, &values[0], number_of_entities, &array3_shape[0], 2);
    IntegrationPointExpressionIO::Write(array3_exp, &VELOCITY);

    ContainerExpression<TContainerType> vector_exp(rModelPart);
    std::vector<int> vector_shape{10, 12};
    CArrayExpressionIO::Read(vector_exp, &values[0], number_of_entities, &vector_shape[0], 2);
    IntegrationPointExpressionIO::Write(vector_exp, &INITIAL_STRAIN);

    ContainerExpression<TContainerType> matrix_exp(rModelPart);
    std::vector<int> matrix_shape{10, 4, 3};
    CArrayExpressionIO::Read(matrix_exp, &values[0], number_of_entities, &matrix_shape[0], 3);
    IntegrationPointExpressionIO::Write(matrix_exp, &NORMAL_SHAPE_DERIVATIVE);

    for (IndexType entity_index = 0; entity_index < rContainer.size(); ++entity_index) {
        auto& r_entity = *(rContainer.begin() + entity_index);

        std::vector<double> double_values, ref_double_values(number_of_integration_points);
        std::iota(ref_double_values.begin(), ref_double_values.end(), entity_index * number_of_integration_points + 1);
        r_entity.CalculateOnIntegrationPoints(PRESSURE, double_values, rModelPart.GetProcessInfo());
        KRATOS_EXPECT_EQ(double_values, ref_double_values);

        std::vector<array_1d<double, 3>> array3_values, ref_array3_values(number_of_integration_points);
        for (IndexType i = 0; i < number_of_integration_points; ++i) {
            ref_array3_values[i][0] = entity_index * number_of_integration_points * 3 + i * 3 + 1;
            ref_array3_values[i][1] = entity_index * number_of_integration_points * 3 + i * 3 + 2;
            ref_array3_values[i][2] = entity_index * number_of_integration_points * 3 + i * 3 + 3;
        }
        r_entity.CalculateOnIntegrationPoints(VELOCITY, array3_values, rModelPart.GetProcessInfo());
        KRATOS_EXPECT_EQ(array3_values, ref_array3_values);

        std::vector<Vector> vector_values, ref_vector_values(number_of_integration_points, Vector(12));
        for (IndexType i = 0; i < number_of_integration_points; ++i) {
            for (IndexType j = 0; j < 12; ++j) {
                ref_vector_values[i][j] = entity_index * number_of_integration_points * 12 + i * 12 + j + 1;
            }
        }
        r_entity.CalculateOnIntegrationPoints(INITIAL_STRAIN, vector_values, rModelPart.GetProcessInfo());
        KRATOS_EXPECT_EQ(vector_values.size(), ref_vector_values.size());
        for (IndexType i = 0; i < number_of_integration_points; ++i) {
            KRATOS_EXPECT_VECTOR_EQ(vector_values[i], ref_vector_values[i]);
        }

        std::vector<Matrix> matrix_values, ref_matrix_values(number_of_integration_points, Matrix(4, 3));
        for (IndexType i = 0; i < number_of_integration_points; ++i) {
            for (IndexType j = 0; j < 12; ++j) {
                *(&(ref_matrix_values[i](0, 0)) + j) = entity_index * number_of_integration_points * 12 + i * 12 + j + 1;
            }
        }
        r_entity.CalculateOnIntegrationPoints(NORMAL_SHAPE_DERIVATIVE, matrix_values, rModelPart.GetProcessInfo());
        KRATOS_EXPECT_EQ(matrix_values.size(), ref_matrix_values.size());
        for (IndexType i = 0; i < number_of_integration_points; ++i) {
            KRATOS_EXPECT_MATRIX_EQ(matrix_values[i], ref_matrix_values[i]);
        }
    }
}

} // namespace ExpressionIntegrationPointIOTestUtilities

KRATOS_TEST_CASE_IN_SUITE(ExpressionIntegrationPointIOReadConditions, KratosCoreFastSuite)
{
    Model model;
    auto& r_model_part = ExpressionIntegrationPointIOTestUtilities::SetUpModelPart(model);
    ExpressionIntegrationPointIOTestUtilities::ExecuteReadTest(r_model_part.Conditions(), r_model_part, 10);
}

KRATOS_TEST_CASE_IN_SUITE(ExpressionIntegrationPointIOReadElements, KratosCoreFastSuite)
{
    Model model;
    auto& r_model_part = ExpressionIntegrationPointIOTestUtilities::SetUpModelPart(model);
    ExpressionIntegrationPointIOTestUtilities::ExecuteReadTest(r_model_part.Elements(), r_model_part, 10);
}

KRATOS_TEST_CASE_IN_SUITE(ExpressionIntegrationPointIOReadEmptyConditions, KratosCoreFastSuite)
{
    Model model;
    auto& r_model_part = ExpressionIntegrationPointIOTestUtilities::SetUpEmptyModelPart(model);
    ExpressionIntegrationPointIOTestUtilities::ExecuteReadTest(r_model_part.Conditions(), r_model_part, 0);
}

KRATOS_TEST_CASE_IN_SUITE(ExpressionIntegrationPointIOReadEmptyElements, KratosCoreFastSuite)
{
    Model model;
    auto& r_model_part = ExpressionIntegrationPointIOTestUtilities::SetUpEmptyModelPart(model);
    ExpressionIntegrationPointIOTestUtilities::ExecuteReadTest(r_model_part.Elements(), r_model_part, 0);
}

KRATOS_TEST_CASE_IN_SUITE(ExpressionIntegrationPointIOWriteConditions, KratosCoreFastSuite)
{
    Model model;
    auto& r_model_part = ExpressionIntegrationPointIOTestUtilities::SetUpModelPart(model);
    ExpressionIntegrationPointIOTestUtilities::ExecuteWriteTest(r_model_part.Conditions(), r_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(ExpressionIntegrationPointIOWriteElements, KratosCoreFastSuite)
{
    Model model;
    auto& r_model_part = ExpressionIntegrationPointIOTestUtilities::SetUpModelPart(model);
    ExpressionIntegrationPointIOTestUtilities::ExecuteWriteTest(r_model_part.Elements(), r_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(ExpressionIntegrationPointIOWriteEmptyConditions, KratosCoreFastSuite)
{
    Model model;
    auto& r_model_part = ExpressionIntegrationPointIOTestUtilities::SetUpEmptyModelPart(model);
    ExpressionIntegrationPointIOTestUtilities::ExecuteWriteTest(r_model_part.Conditions(), r_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(ExpressionIntegrationPointIOWriteEmptyElements, KratosCoreFastSuite)
{
    Model model;
    auto& r_model_part = ExpressionIntegrationPointIOTestUtilities::SetUpEmptyModelPart(model);
    ExpressionIntegrationPointIOTestUtilities::ExecuteWriteTest(r_model_part.Elements(), r_model_part);
}

}  // namespace Kratos::Testing.

