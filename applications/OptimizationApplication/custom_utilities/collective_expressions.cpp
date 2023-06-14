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
#include "expression/binary_expression.h"
#include "expression/literal_expression.h"

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
        std::visit([&result](const auto& v) {
            result.Add(v->Clone());
        }, p_container_variable_data_holder);
    }
    return result;
}

void CollectiveExpressions::SetToZero()
{
    for (const auto& p_container_variable_data_holder : mExpressionPointersList) {
        std::visit([](const auto& v) {
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

void CollectiveExpressions::Read(
    double const* pBegin,
    int const* NumberOfEntities,
    int const** pListShapeBegin,
    int const* ShapeSizes,
    const int NumberOfContainers)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(NumberOfContainers > 0 && static_cast<IndexType>(NumberOfContainers) == this->mExpressionPointersList.size())
        << "Number of containers mismatch. [ Input number of containers = " << NumberOfContainers
        << ", CollectiveExpressions number of containers = "
        << this->mExpressionPointersList.size() << " ].\n";

    for (const auto& p_container_variable_data_holder : mExpressionPointersList) {
        std::visit([&pBegin, &pListShapeBegin, &ShapeSizes, &NumberOfEntities](const auto& v) {
            v->Read(pBegin, *NumberOfEntities, *pListShapeBegin, *ShapeSizes);

            // now offset everything
            pBegin += v->GetContainer().size() * v->GetItemComponentCount();
            ++pListShapeBegin;
            ++ShapeSizes;
            ++NumberOfEntities;
        }, p_container_variable_data_holder);
    }

    KRATOS_CATCH("");
}

void CollectiveExpressions::Read(const VariableTypes& rVariable)
{
    KRATOS_TRY

    for (const auto& p_container_variable_data_holder : mExpressionPointersList) {
        std::visit([](const auto& v, const auto pVariable) {
            v->Read(*pVariable);
        }, p_container_variable_data_holder, rVariable);
    }

    KRATOS_CATCH("");
}

void CollectiveExpressions::Read(const std::vector<VariableTypes>& rVariables)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(rVariables.size() == mExpressionPointersList.size())
        << "Variables and container expressions size mismatch. [ Provided number of variables: "
        << rVariables.size() << ", number of ContainerExpressions = "
        << mExpressionPointersList.size() << " ].\n";

    for (IndexType i = 0; i < rVariables.size(); ++i) {
        std::visit([](const auto& v, const auto pVariable) {
            v->Read(*pVariable);
        }, mExpressionPointersList[i], rVariables[i]);
    }

    KRATOS_CATCH("");
}

void CollectiveExpressions::MoveFrom(
    double* pBegin,
    int const* NumberOfEntities,
    int const** pListShapeBegin,
    int const* ShapeSizes,
    const int NumberOfContainers)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(NumberOfContainers > 0 && static_cast<IndexType>(NumberOfContainers) == this->mExpressionPointersList.size())
        << "Number of containers mismatch. [ Input number of containers = " << NumberOfContainers
        << ", CollectiveExpressions number of containers = "
        << this->mExpressionPointersList.size() << " ].\n";

    for (const auto& p_container_variable_data_holder : mExpressionPointersList) {
        std::visit([&pBegin, &pListShapeBegin, &ShapeSizes, &NumberOfEntities](const auto& v) {
            v->MoveFrom(pBegin, *NumberOfEntities, *pListShapeBegin, *ShapeSizes);

            // now offset everything
            pBegin += v->GetContainer().size() * v->GetItemComponentCount();
            ++pListShapeBegin;
            ++ShapeSizes;
            ++NumberOfEntities;
        }, p_container_variable_data_holder);
    }

    KRATOS_CATCH("");
}

void CollectiveExpressions::Evaluate(
    double* pBegin,
    const int Size) const
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(Size > 0 && static_cast<IndexType>(Size) == this->GetCollectiveFlattenedDataSize())
        << "The size of the double vector does not match with the required "
           "collective expression size. [ "
           "Size = "
        << Size << ", collective expression data size = "
        << this->GetCollectiveFlattenedDataSize() << " ].\n";

    for (const auto& p_container_variable_data_holder : mExpressionPointersList) {
        std::visit([&pBegin](const auto& v) {
            // get the shape of the container expression.
            const auto& r_shape = v->GetItemShape();

            // transform unsigned Index type shape to signed int shape.
            std::vector<int> shape(r_shape.size());
            std::transform(r_shape.begin(), r_shape.end(), shape.begin(), [](const IndexType Value) -> int { return Value; });

            // get the number of entities in the container.
            const auto number_of_entities = v->GetContainer().size();

            // evaluate the expression and put the result in a continuous array starting with pBegin.
            v->Evaluate(pBegin, number_of_entities, shape.data(), shape.size());

            // increase the offset to place the evaluated values of the next container expression correctly.
            pBegin += number_of_entities * v->GetItemComponentCount();
        }, p_container_variable_data_holder);
    }

    KRATOS_CATCH("");
}

void CollectiveExpressions::Evaluate(const VariableTypes& rVariable)
{
    KRATOS_TRY

    for (const auto& p_container_variable_data_holder : mExpressionPointersList) {
        std::visit([](const auto& v, const auto pVariable) {
            v->Evaluate(*pVariable);
        }, p_container_variable_data_holder, rVariable);
    }

    KRATOS_CATCH("");
}

void CollectiveExpressions::Evaluate(const std::vector<VariableTypes>& rVariables)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(rVariables.size() == mExpressionPointersList.size())
        << "Variables and container expressions size mismatch. [ Provided number of variables: "
        << rVariables.size() << ", number of ContainerExpressions = "
        << mExpressionPointersList.size() << " ].\n";

    for (IndexType i = 0; i < rVariables.size(); ++i) {
        std::visit([](const auto& v, const auto pVariable) {
            v->Evaluate(*pVariable);
        }, mExpressionPointersList[i], rVariables[i]);
    }

    KRATOS_CATCH("");
}

IndexType CollectiveExpressions::GetCollectiveFlattenedDataSize() const
{
    IndexType size = 0;
    for (const auto& p_container_variable_data_holder : mExpressionPointersList) {
        std::visit([&size](const auto& v) {
            size += v->GetContainer().size() * v->GetItemComponentCount();
        }, p_container_variable_data_holder);
    }

    return size;
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
        const auto& r_other_expression = rOther.mExpressionPointersList[i];
        std::visit([&r_other_expression, &is_compatible](const auto& v) {
            using v_type = std::decay_t<decltype(v)>;
            auto* ptr = std::get_if<v_type>(&r_other_expression);
            is_compatible = is_compatible && (ptr != nullptr) && ((*ptr)->GetContainer().size() == v->GetContainer().size());
        }, mExpressionPointersList[i]);
    }

    return is_compatible;
}

std::string CollectiveExpressions::Info() const
{
    std::stringstream msg;

    msg << "CollectiveExpressions contains following data holders:\n";

    for (const auto& p_container_variable_data_holder : mExpressionPointersList) {
        std::visit([&msg](const auto& v) { msg << "\t" << *v; }, p_container_variable_data_holder);
    }

    return msg.str();
}

CollectiveExpressions CollectiveExpressions::operator+(const CollectiveExpressions& rOther) const
{
    KRATOS_ERROR_IF_NOT(IsCompatibleWith(rOther))
        << "Unsupported collective variable data holders provided for \"+\" operation."
        << "\nLeft operand : " << *this << "\nRight operand: " << rOther << std::endl;

    CollectiveExpressions result(*this);
    for (IndexType i = 0; i < mExpressionPointersList.size(); ++i) {
        const auto& r_other_expression = rOther.mExpressionPointersList[i];
        std::visit([&r_other_expression](const auto& v) {
            using v_type = std::decay_t<decltype(v)>;
            v->operator+=(*(std::get<v_type>(r_other_expression)));
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
        const auto& r_other_expression = rOther.mExpressionPointersList[i];
        std::visit([&r_other_expression](const auto& v) {
            using v_type = std::decay_t<decltype(v)>;
            v->operator+=(*(std::get<v_type>(r_other_expression)));
        }, mExpressionPointersList[i]);
    }

    return *this;
}

CollectiveExpressions CollectiveExpressions::operator+(const double Value) const
{
    CollectiveExpressions result(*this);
    for (IndexType i = 0; i < mExpressionPointersList.size(); ++i) {
        std::visit([Value](const auto& v) {
            v->operator+=(Value);
        }, result.mExpressionPointersList[i]);
    }

    return result;
}

CollectiveExpressions& CollectiveExpressions::operator+=(const double Value)
{
    for (IndexType i = 0; i < mExpressionPointersList.size(); ++i) {
        std::visit([Value](const auto& v) {
            v->operator+=(Value);
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
        const auto& r_other_expression = rOther.mExpressionPointersList[i];
        std::visit([&r_other_expression](const auto& v) {
            using v_type = std::decay_t<decltype(v)>;
            v->operator-=(*(std::get<v_type>(r_other_expression)));
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
        const auto& r_other_expression = rOther.mExpressionPointersList[i];
        std::visit([&r_other_expression](const auto& v) {
            using v_type = std::decay_t<decltype(v)>;
            v->operator-=(*(std::get<v_type>(r_other_expression)));
        }, mExpressionPointersList[i]);
    }

    return *this;
}

CollectiveExpressions CollectiveExpressions::operator-(const double Value) const
{
    CollectiveExpressions result(*this);
    for (IndexType i = 0; i < mExpressionPointersList.size(); ++i) {
        std::visit([Value](const auto& v) {
            v->operator-=(Value);
        }, result.mExpressionPointersList[i]);
    }

    return result;
}

CollectiveExpressions& CollectiveExpressions::operator-=(const double Value)
{
    for (IndexType i = 0; i < mExpressionPointersList.size(); ++i) {
        std::visit([Value](const auto& v) {
            v->operator-=(Value);
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
        const auto& r_other_expression = rOther.mExpressionPointersList[i];
        std::visit([&r_other_expression](const auto& v) {
            using v_type = std::decay_t<decltype(v)>;
            v->operator*=(*(std::get<v_type>(r_other_expression)));
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
        const auto& r_other_expression = rOther.mExpressionPointersList[i];
        std::visit([&r_other_expression](const auto& v) {
            using v_type = std::decay_t<decltype(v)>;
            v->operator*=(*(std::get<v_type>(r_other_expression)));
        }, mExpressionPointersList[i]);
    }

    return *this;
}

CollectiveExpressions CollectiveExpressions::operator*(const double Value) const
{
    CollectiveExpressions result(*this);
    for (IndexType i = 0; i < mExpressionPointersList.size(); ++i) {
        std::visit([Value](const auto& v) {
            v->operator*=(Value);
        }, result.mExpressionPointersList[i]);
    }

    return result;
}

CollectiveExpressions& CollectiveExpressions::operator*=(const double Value)
{
    for (IndexType i = 0; i < mExpressionPointersList.size(); ++i) {
        std::visit([Value](const auto& v) {
            v->operator*=(Value);
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
        const auto& r_other_expression = rOther.mExpressionPointersList[i];
        std::visit([&r_other_expression](const auto& v) {
            using v_type = std::decay_t<decltype(v)>;
            v->operator/=(*(std::get<v_type>(r_other_expression)));
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
        const auto& r_other_expression = rOther.mExpressionPointersList[i];
        std::visit([&r_other_expression](const auto& v) {
            using v_type = std::decay_t<decltype(v)>;
            v->operator/=(*(std::get<v_type>(r_other_expression)));
        }, mExpressionPointersList[i]);
    }

    return *this;
}

CollectiveExpressions CollectiveExpressions::operator/(const double Value) const
{
    CollectiveExpressions result(*this);
    for (IndexType i = 0; i < mExpressionPointersList.size(); ++i) {
        std::visit([Value](const auto& v) {
            v->operator/=(Value);
        }, result.mExpressionPointersList[i]);
    }

    return result;
}

CollectiveExpressions& CollectiveExpressions::operator/=(const double Value)
{
    for (IndexType i = 0; i < mExpressionPointersList.size(); ++i) {
        std::visit([Value](const auto& v) {
            v->operator/=(Value);
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
        const auto& r_other_expression = rOther.mExpressionPointersList[i];
        std::visit([&r_other_expression](const auto& v) {
            using v_type = std::decay_t<decltype(v)>;
            v->SetExpression(
                BinaryExpression<BinaryOperations::Power>::Create(
                    v->pGetExpression(),
                    std::get<v_type>(r_other_expression)->pGetExpression()));
        }, result.mExpressionPointersList[i]);
    }

    return result;
}

CollectiveExpressions CollectiveExpressions::Pow(const double Value) const
{
    CollectiveExpressions result(*this);
    for (IndexType i = 0; i < mExpressionPointersList.size(); ++i) {
        std::visit([Value](const auto& v) {
            v->SetExpression(
                BinaryExpression<BinaryOperations::Power>::Create(
                    v->pGetExpression(),
                    LiteralExpression<double>::Create(Value, v->GetExpression().NumberOfEntities())));
        }, result.mExpressionPointersList[i]);
    }
    return result;
}

} // namespace Kratos