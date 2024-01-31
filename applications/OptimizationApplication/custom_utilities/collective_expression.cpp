//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

// System includes
#include <sstream>
#include <string>
#include <vector>
#include <variant>

// Project includes
#include "expression/arithmetic_operators.h"

// Application includes

// Include base h
#include "collective_expression.h"

namespace Kratos {

///@name Kratos Classes
///@{

CollectiveExpression::CollectiveExpression(const std::vector<CollectiveExpressionType>& rContainerVariableDataHolderPointersList)
{
    for (const auto& p_container_variable_data_holder : rContainerVariableDataHolderPointersList){
        this->Add(p_container_variable_data_holder);
    }
}

CollectiveExpression::CollectiveExpression(const CollectiveExpression& rOther)
{
    for (const auto& p_container_variable_data_holder : rOther.mExpressionPointersList) {
        std::visit([&](const auto& v) {
            mExpressionPointersList.push_back(v->Clone());
        }, p_container_variable_data_holder);
    }
}

CollectiveExpression& CollectiveExpression::operator=(const CollectiveExpression& rOther)
{
    mExpressionPointersList.clear();
    for (const auto& r_container_variable_data_holder : rOther.mExpressionPointersList) {
        std::visit([this](const auto& v) {
            mExpressionPointersList.push_back(v);
        }, r_container_variable_data_holder);
    }
    return *this;
}

CollectiveExpression CollectiveExpression::Clone()  const
{
    CollectiveExpression result;
    for (const auto& p_container_variable_data_holder : mExpressionPointersList) {
        std::visit([&result](const auto& v) {
            result.Add(v->Clone());
        }, p_container_variable_data_holder);
    }
    return result;
}

void CollectiveExpression::Add(const CollectiveExpressionType& pVariableDataHolder)
{
    std::visit([&](const auto& v) {
        mExpressionPointersList.push_back(v);
    }, pVariableDataHolder);
}

void CollectiveExpression::Add(const CollectiveExpression& rCollectiveVariableDataHolder)
{
    for (const auto& p_container_variable_data_holder : rCollectiveVariableDataHolder.mExpressionPointersList) {
        std::visit([&](const auto& v) {
            mExpressionPointersList.push_back(v);
        }, p_container_variable_data_holder);
    }
}

void CollectiveExpression::Clear()
{
    mExpressionPointersList.clear();
}

IndexType CollectiveExpression::GetCollectiveFlattenedDataSize() const
{
    IndexType size = 0;
    for (const auto& p_container_variable_data_holder : mExpressionPointersList) {
        std::visit([&size](const auto& v) {
            size += v->GetContainer().size() * v->GetItemComponentCount();
        }, p_container_variable_data_holder);
    }

    return size;
}

std::vector<CollectiveExpression::CollectiveExpressionType> CollectiveExpression::GetContainerExpressions()
{
    return mExpressionPointersList;
}

std::vector<CollectiveExpression::CollectiveExpressionType> CollectiveExpression::GetContainerExpressions() const
{
    return mExpressionPointersList;
}

bool CollectiveExpression::IsCompatibleWith(const CollectiveExpression& rOther) const
{
    if (mExpressionPointersList.size() != rOther.mExpressionPointersList.size()) {
        return false;
    }

    bool is_compatible = true;

    for (IndexType i = 0; i < mExpressionPointersList.size(); ++i) {
        const auto& r_other_expression = rOther.mExpressionPointersList[i];
        std::visit([&r_other_expression, &is_compatible](const auto& v) {
            using v_type = std::decay_t<decltype(v)>;
            auto* ptr = std::get_if<v_type>(&r_other_expression);
            is_compatible = is_compatible && (ptr != nullptr) && ((*ptr)->GetContainer().size() == v->GetContainer().size());
        }, mExpressionPointersList[i]);
    }

    return is_compatible;
}

std::string CollectiveExpression::Info() const
{
    std::stringstream msg;

    msg << "CollectiveExpression contains following data holders:\n";

    for (const auto& p_container_variable_data_holder : mExpressionPointersList) {
        std::visit([&msg](const auto& v) { msg << "\t" << *v; }, p_container_variable_data_holder);
    }

    return msg.str();
}

#define KRATOS_DEFINE_BINARY_COLLECTIVE_EXPRESSION_OPERATOR(OPERATOR_NAME)                    \
    CollectiveExpression OPERATOR_NAME(const CollectiveExpression& rLeft, const double Right) \
    {                                                                                         \
        KRATOS_TRY                                                                            \
                                                                                              \
        auto result = rLeft;                                                                  \
        auto r_list_of_container_expressions = result.GetContainerExpressions();              \
        for (IndexType i = 0; i < r_list_of_container_expressions.size(); ++i) {              \
            std::visit(                                                                       \
                [Right](auto& pResult) {                                                      \
                    pResult->SetExpression(OPERATOR_NAME(pResult->pGetExpression(), Right));  \
                },                                                                            \
                r_list_of_container_expressions[i]);                                          \
        }                                                                                     \
        return result;                                                                        \
                                                                                              \
        KRATOS_CATCH("")                                                                      \
    }                                                                                         \
                                                                                              \
    CollectiveExpression OPERATOR_NAME(const double Left, const CollectiveExpression& rRight) \
    {                                                                                         \
        KRATOS_TRY                                                                            \
                                                                                              \
        auto result = rRight;                                                                 \
        auto r_list_of_container_expressions = result.GetContainerExpressions();              \
        for (IndexType i = 0; i < r_list_of_container_expressions.size(); ++i) {              \
            std::visit(                                                                       \
                [Left](auto& pResult) {                                                       \
                    pResult->SetExpression(OPERATOR_NAME(pResult->pGetExpression(), Left));   \
                },                                                                            \
                r_list_of_container_expressions[i]);                                          \
        }                                                                                     \
        return result;                                                                        \
                                                                                              \
        KRATOS_CATCH("")                                                                      \
    }                                                                                         \
                                                                                              \
    CollectiveExpression OPERATOR_NAME(const CollectiveExpression& rLeft,                     \
                                       const CollectiveExpression& rRight)                    \
    {                                                                                         \
        KRATOS_TRY                                                                            \
                                                                                              \
        KRATOS_ERROR_IF_NOT(rLeft.IsCompatibleWith(rRight))                                   \
            << "Unsupported collective variable data holders provided for "                   \
               "\""                                                                           \
            << #OPERATOR_NAME << "\"."                                                        \
            << "\nLeft operand : " << rLeft << "\nRight operand: " << rRight                  \
            << std::endl;                                                                     \
                                                                                              \
        auto result = rLeft;                                                                  \
        auto r_list_of_container_expressions = result.GetContainerExpressions();              \
        const auto& r_right_container_expressions = rRight.GetContainerExpressions();         \
        for (IndexType i = 0; i < r_list_of_container_expressions.size(); ++i) {              \
            std::visit(                                                                       \
                [&r_right_container_expressions, i](auto& pResult) {                          \
                    auto p_right = std::get<std::decay_t<decltype(pResult)>>(                 \
                        r_right_container_expressions[i]);                                    \
                    pResult->SetExpression(OPERATOR_NAME(                                     \
                        pResult->pGetExpression(), p_right->pGetExpression()));               \
                },                                                                            \
                r_list_of_container_expressions[i]);                                          \
        }                                                                                     \
        return result;                                                                        \
                                                                                              \
        KRATOS_CATCH("");                                                                     \
    }

#define KRATOS_DEFINE_UNARY_COLLECTIVE_EXPRESSION_OPERATOR(OPERATOR_NAME, EXPRESSION_OPERATOR_NAME) \
    CollectiveExpression& CollectiveExpression::OPERATOR_NAME(const double Value)                   \
    {                                                                                               \
        KRATOS_TRY                                                                                  \
                                                                                                    \
        for (auto& p_expression : mExpressionPointersList) {                                        \
            std::visit(                                                                             \
                [Value](auto& pContainerExpression) {                                               \
                    pContainerExpression->SetExpression(EXPRESSION_OPERATOR_NAME(                   \
                        pContainerExpression->pGetExpression(), Value));                            \
                },                                                                                  \
                p_expression);                                                                      \
        }                                                                                           \
                                                                                                    \
        return *this;                                                                               \
                                                                                                    \
        KRATOS_CATCH("");                                                                           \
    }                                                                                               \
                                                                                                    \
    CollectiveExpression& CollectiveExpression::OPERATOR_NAME(const CollectiveExpression& rOther)   \
    {                                                                                               \
        KRATOS_TRY                                                                                  \
                                                                                                    \
        KRATOS_ERROR_IF_NOT(IsCompatibleWith(rOther))                                               \
            << "Unsupported collective variable data holders provided for "                         \
               "\""                                                                                 \
            << #OPERATOR_NAME << "\"."                                                              \
            << "\nLeft operand : " << *this << "\nRight operand: " << rOther                        \
            << std::endl;                                                                           \
                                                                                                    \
        const auto& r_other_epxressions_list = rOther.GetContainerExpressions();                    \
        for (IndexType i = 0; i < mExpressionPointersList.size(); ++i) {                            \
            std::visit(                                                                             \
                [&r_other_epxressions_list, i](auto& pContainerExpression) {                        \
                    auto p_other =                                                                  \
                        std::get<std::decay_t<decltype(pContainerExpression)>>(                     \
                            r_other_epxressions_list[i]);                                           \
                                                                                                    \
                    pContainerExpression->SetExpression(EXPRESSION_OPERATOR_NAME(                   \
                        pContainerExpression->pGetExpression(), p_other->pGetExpression()));        \
                },                                                                                  \
                mExpressionPointersList[i]);                                                        \
        }                                                                                           \
                                                                                                    \
        return *this;                                                                               \
                                                                                                    \
        KRATOS_CATCH("");                                                                           \
    }

KRATOS_DEFINE_BINARY_COLLECTIVE_EXPRESSION_OPERATOR(operator+)
KRATOS_DEFINE_BINARY_COLLECTIVE_EXPRESSION_OPERATOR(operator-)
KRATOS_DEFINE_BINARY_COLLECTIVE_EXPRESSION_OPERATOR(operator*)
KRATOS_DEFINE_BINARY_COLLECTIVE_EXPRESSION_OPERATOR(operator/)

KRATOS_DEFINE_UNARY_COLLECTIVE_EXPRESSION_OPERATOR(operator+=, operator+)
KRATOS_DEFINE_UNARY_COLLECTIVE_EXPRESSION_OPERATOR(operator-=, operator-)
KRATOS_DEFINE_UNARY_COLLECTIVE_EXPRESSION_OPERATOR(operator*=, operator*)
KRATOS_DEFINE_UNARY_COLLECTIVE_EXPRESSION_OPERATOR(operator/=, operator/)

#undef KRATOS_DEFINE_BINARY_COLLECTIVE_EXPRESSION_OPERATOR
#undef KRATOS_DEFINE_UNARY_COLLECTIVE_EXPRESSION_OPERATOR

} // namespace Kratos