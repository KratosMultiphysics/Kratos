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
#include "containers/container_expression/expressions/binary/binary_expression.h"
#include "containers/container_expression/expressions/literal/literal_expression.h"

// Application includes

// Include base h
#include "collective_expressions.h"

namespace Kratos {

///@name Kratos Classes
///@{

CollectiveExpressions::CollectiveExpressions(const std::vector<CollectiveExpressionType>& rContainerVariableDataHolderPointersList)
{
    for (const auto& p_container_variable_data_holder : rContainerVariableDataHolderPointersList){
        this->Add(p_container_variable_data_holder);
    }
}

CollectiveExpressions::CollectiveExpressions(const CollectiveExpressions& rOther)
{
    for (const auto& p_container_variable_data_holder : rOther.mExpressionPointersList) {
        std::visit([&](const auto& v) {
            mExpressionPointersList.push_back(v->Clone());
        }, p_container_variable_data_holder);
    }
}

CollectiveExpressions CollectiveExpressions::Clone()  const
{
    CollectiveExpressions result;
    for (const auto& p_container_variable_data_holder : mExpressionPointersList) {
        std::visit([&](const auto& v) {
            result.Add(v->Clone());
        }, p_container_variable_data_holder);
    }
    return result;
}

void CollectiveExpressions::SetToZero()
{
    for (const auto& p_container_variable_data_holder : mExpressionPointersList) {
        std::visit([&](const auto& v) {
            v->SetDataToZero();
        }, p_container_variable_data_holder);
    }
}

void CollectiveExpressions::Add(const CollectiveExpressionType& pVariableDataHolder)
{
    std::visit([&](const auto& v) {
        mExpressionPointersList.push_back(v);
    }, pVariableDataHolder);
}

void CollectiveExpressions::Add(const CollectiveExpressions& rCollectiveVariableDataHolder)
{
    for (const auto& p_container_variable_data_holder : rCollectiveVariableDataHolder.mExpressionPointersList) {
        std::visit([&](const auto& v) {
            mExpressionPointersList.push_back(v);
        }, p_container_variable_data_holder);
    }
}

void CollectiveExpressions::Clear()
{
    mExpressionPointersList.clear();
}

std::vector<CollectiveExpressions::CollectiveExpressionType> CollectiveExpressions::GetContainerExpressions()
{
    return mExpressionPointersList;
}

std::vector<CollectiveExpressions::CollectiveExpressionType> CollectiveExpressions::GetContainerExpressions() const
{
    return mExpressionPointersList;
}

bool CollectiveExpressions::IsCompatibleWith(const CollectiveExpressions& rOther) const
{
    if (mExpressionPointersList.size() != rOther.mExpressionPointersList.size()) {
        return false;
    }

    bool is_compatible = true;

    for (IndexType i = 0; i < mExpressionPointersList.size(); ++i) {
        std::visit([&](auto&& v) {
            using v_type = std::decay_t<decltype(v)>;
            auto* ptr = std::get_if<v_type>(&(rOther.mExpressionPointersList[i]));
            is_compatible = is_compatible && (ptr != nullptr) && (&((*ptr)->GetModelPart()) == &(v->GetModelPart()));
        }, mExpressionPointersList[i]);
    }

    return is_compatible;
}

std::string CollectiveExpressions::Info() const
{
    std::stringstream msg;

    msg << "CollectiveExpressions contains following data holders:\n";

    for (const auto& p_container_variable_data_holder : mExpressionPointersList) {
        std::visit([&msg](auto&& v) { msg << "\t" << *v; }, p_container_variable_data_holder);
    }

    return msg.str();
}

CollectiveExpressions CollectiveExpressions::operator+(const CollectiveExpressions& rOther) const
{
    KRATOS_ERROR_IF_NOT(IsCompatibleWith(rOther))
        << "Unsupported collective variable data holders provided for \"+\" operation."
        << "\nLeft operand : " << *this << "\nRight operand: " << rOther << std::endl;

    CollectiveExpressions result(rOther);
    for (IndexType i = 0; i < mExpressionPointersList.size(); ++i) {
        std::visit([&](auto&& v) {
            using v_type = std::decay_t<decltype(v)>;
            v->operator+=(*(std::get<v_type>(mExpressionPointersList[i])));
        }, result.mExpressionPointersList[i]);
    }

    return result;
}

CollectiveExpressions& CollectiveExpressions::operator+=(const CollectiveExpressions& rOther)
{
    KRATOS_ERROR_IF_NOT(IsCompatibleWith(rOther))
        << "Unsupported collective variable data holders provided for \"+=\" operation."
        << "\nLeft operand : " << *this << "\nRight operand: " << rOther << std::endl;

    for (IndexType i = 0; i < mExpressionPointersList.size(); ++i) {
        std::visit([&](auto&& v) {
            using v_type = std::decay_t<decltype(v)>;
            v->operator+=(*(std::get<v_type>(rOther.mExpressionPointersList[i])));
        }, mExpressionPointersList[i]);
    }

    return *this;
}

CollectiveExpressions CollectiveExpressions::operator+(const double rOther) const
{
    CollectiveExpressions result(*this);
    for (IndexType i = 0; i < mExpressionPointersList.size(); ++i) {
        std::visit([&](auto&& v) {
            v->operator+=(rOther);
        }, result.mExpressionPointersList[i]);
    }

    return result;
}

CollectiveExpressions& CollectiveExpressions::operator+=(const double rOther)
{
    for (IndexType i = 0; i < mExpressionPointersList.size(); ++i) {
        std::visit([&](auto&& v) {
            v->operator+=(rOther);
        }, mExpressionPointersList[i]);
    }

    return *this;
}

CollectiveExpressions CollectiveExpressions::operator-(const CollectiveExpressions& rOther) const
{
    KRATOS_ERROR_IF_NOT(IsCompatibleWith(rOther))
        << "Unsupported collective variable data holders provided for \"-\" operation."
        << "\nLeft operand : " << *this << "\nRight operand: " << rOther << std::endl;

    CollectiveExpressions result(*this);
    for (IndexType i = 0; i < mExpressionPointersList.size(); ++i) {
        std::visit([&](auto&& v) {
            using v_type = std::decay_t<decltype(v)>;
            v->operator-=(*(std::get<v_type>(rOther.mExpressionPointersList[i])));
        }, result.mExpressionPointersList[i]);
    }

    return result;
}

CollectiveExpressions& CollectiveExpressions::operator-=(const CollectiveExpressions& rOther)
{
    KRATOS_ERROR_IF_NOT(IsCompatibleWith(rOther))
        << "Unsupported collective variable data holders provided for \"-=\" operation."
        << "\nLeft operand : " << *this << "\nRight operand: " << rOther << std::endl;

    for (IndexType i = 0; i < mExpressionPointersList.size(); ++i) {
        std::visit([&](auto&& v) {
            using v_type = std::decay_t<decltype(v)>;
            v->operator-=(*(std::get<v_type>(rOther.mExpressionPointersList[i])));
        }, mExpressionPointersList[i]);
    }

    return *this;
}

CollectiveExpressions CollectiveExpressions::operator-(const double rOther) const
{
    CollectiveExpressions result(*this);
    for (IndexType i = 0; i < mExpressionPointersList.size(); ++i) {
        std::visit([&](auto&& v) {
            v->operator-=(rOther);
        }, result.mExpressionPointersList[i]);
    }

    return result;
}

CollectiveExpressions& CollectiveExpressions::operator-=(const double rOther)
{
    for (IndexType i = 0; i < mExpressionPointersList.size(); ++i) {
        std::visit([&](auto&& v) {
            v->operator-=(rOther);
        }, mExpressionPointersList[i]);
    }

    return *this;
}

CollectiveExpressions CollectiveExpressions::operator*(const CollectiveExpressions& rOther) const
{
    KRATOS_ERROR_IF_NOT(IsCompatibleWith(rOther))
        << "Unsupported collective variable data holders provided for \"-\" operation."
        << "\nLeft operand : " << *this << "\nRight operand: " << rOther << std::endl;

    CollectiveExpressions result(*this);
    for (IndexType i = 0; i < mExpressionPointersList.size(); ++i) {
        std::visit([&](auto&& v) {
            using v_type = std::decay_t<decltype(v)>;
            v->operator*=(*(std::get<v_type>(rOther.mExpressionPointersList[i])));
        }, result.mExpressionPointersList[i]);
    }

    return result;
}

CollectiveExpressions& CollectiveExpressions::operator*=(const CollectiveExpressions& rOther)
{
    KRATOS_ERROR_IF_NOT(IsCompatibleWith(rOther))
        << "Unsupported collective variable data holders provided for \"-=\" operation."
        << "\nLeft operand : " << *this << "\nRight operand: " << rOther << std::endl;

    for (IndexType i = 0; i < mExpressionPointersList.size(); ++i) {
        std::visit([&](auto&& v) {
            using v_type = std::decay_t<decltype(v)>;
            v->operator*=(*(std::get<v_type>(rOther.mExpressionPointersList[i])));
        }, mExpressionPointersList[i]);
    }

    return *this;
}

CollectiveExpressions CollectiveExpressions::operator*(const double rOther) const
{
    CollectiveExpressions result(*this);
    for (IndexType i = 0; i < mExpressionPointersList.size(); ++i) {
        std::visit([&](auto&& v) {
            v->operator*=(rOther);
        }, result.mExpressionPointersList[i]);
    }

    return result;
}

CollectiveExpressions& CollectiveExpressions::operator*=(const double rOther)
{
    for (IndexType i = 0; i < mExpressionPointersList.size(); ++i) {
        std::visit([&](auto&& v) {
            v->operator*=(rOther);
        }, mExpressionPointersList[i]);
    }

    return *this;
}

CollectiveExpressions CollectiveExpressions::operator/(const CollectiveExpressions& rOther) const
{
    KRATOS_ERROR_IF_NOT(IsCompatibleWith(rOther))
        << "Unsupported collective variable data holders provided for \"-\" operation."
        << "\nLeft operand : " << *this << "\nRight operand: " << rOther << std::endl;

    CollectiveExpressions result(*this);
    for (IndexType i = 0; i < mExpressionPointersList.size(); ++i) {
        std::visit([&](auto&& v) {
            using v_type = std::decay_t<decltype(v)>;
            v->operator/=(*(std::get<v_type>(rOther.mExpressionPointersList[i])));
        }, result.mExpressionPointersList[i]);
    }

    return result;
}

CollectiveExpressions& CollectiveExpressions::operator/=(const CollectiveExpressions& rOther)
{
    KRATOS_ERROR_IF_NOT(IsCompatibleWith(rOther))
        << "Unsupported collective variable data holders provided for \"-=\" operation."
        << "\nLeft operand : " << *this << "\nRight operand: " << rOther << std::endl;

    for (IndexType i = 0; i < mExpressionPointersList.size(); ++i) {
        std::visit([&](auto&& v) {
            using v_type = std::decay_t<decltype(v)>;
            v->operator/=(*(std::get<v_type>(rOther.mExpressionPointersList[i])));
        }, mExpressionPointersList[i]);
    }

    return *this;
}

CollectiveExpressions CollectiveExpressions::operator/(const double rOther) const
{
    CollectiveExpressions result(*this);
    for (IndexType i = 0; i < mExpressionPointersList.size(); ++i) {
        std::visit([&](auto&& v) {
            v->operator/=(rOther);
        }, result.mExpressionPointersList[i]);
    }

    return result;
}

CollectiveExpressions& CollectiveExpressions::operator/=(const double rOther)
{
    for (IndexType i = 0; i < mExpressionPointersList.size(); ++i) {
        std::visit([&](auto&& v) {
            v->operator/=(rOther);
        }, mExpressionPointersList[i]);
    }

    return *this;
}

CollectiveExpressions CollectiveExpressions::Pow(const CollectiveExpressions& rOther) const
{
    KRATOS_ERROR_IF_NOT(IsCompatibleWith(rOther))
        << "Unsupported collective variable data holders provided for \"-\" operation."
        << "\nLeft operand : " << *this << "\nRight operand: " << rOther << std::endl;

    CollectiveExpressions result(*this);
    for (IndexType i = 0; i < mExpressionPointersList.size(); ++i) {
        std::visit([&](auto& v) {
            using v_type = std::decay_t<decltype(v)>;
            v->SetExpression(
                BinaryExpression<BinaryOperations::Power>::Create(
                    v->pGetExpression(),
                    std::get<v_type>(rOther.mExpressionPointersList[i])->pGetExpression()));
        }, result.mExpressionPointersList[i]);
    }

    return result;
}

CollectiveExpressions CollectiveExpressions::Pow(const double Value) const
{
    CollectiveExpressions result(*this);
    for (IndexType i = 0; i < mExpressionPointersList.size(); ++i) {
        std::visit([&](auto&& v) {
            v->SetExpression(
                BinaryExpression<BinaryOperations::Power>::Create(
                    v->pGetExpression(),
                    LiteralExpression<double>::Create(Value)));
        }, result.mExpressionPointersList[i]);
    }
    return result;
}

} // namespace Kratos