//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: HDF5Application/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

// System includes
#include <sstream>
#include <string>
#include <vector>
#include <variant>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

// Application includes

// Include base h
#include "collective_variable_data_holder.h"

namespace Kratos {

///@name Kratos Classes
///@{

CollectiveVariableDataHolder::CollectiveVariableDataHolder(const std::vector<ContainerVariableDataHolderPointerVariantType>& rContainerVariableDataHolderPointersList)
{
    for (const auto& p_container_variable_data_holder : rContainerVariableDataHolderPointersList){
        this->AddVariableDataHolder(p_container_variable_data_holder);
    }
}

CollectiveVariableDataHolder::CollectiveVariableDataHolder(const CollectiveVariableDataHolder& rOther)
{
    for (const auto& p_container_variable_data_holder : rOther.mContainerVariableDataHolderPointersList) {
        std::visit([&](auto&& v) {
            mContainerVariableDataHolderPointersList.push_back(v->Clone());
        }, p_container_variable_data_holder);
    }
}

CollectiveVariableDataHolder CollectiveVariableDataHolder::CloneWithDataInitializedToZero()  const
{
    CollectiveVariableDataHolder result;
    for (const auto& p_container_variable_data_holder : mContainerVariableDataHolderPointersList) {
        std::visit([&](auto&& v) {
            result.AddVariableDataHolder(v->CloneWithDataInitializedToZero());
        }, p_container_variable_data_holder);
    }
    return result;
}

void CollectiveVariableDataHolder::AddVariableDataHolder(const ContainerVariableDataHolderPointerVariantType& pVariableDataHolder)
{
    std::visit([&](auto&& v) {
        mContainerVariableDataHolderPointersList.push_back(v);
    }, pVariableDataHolder);
}

std::vector<CollectiveVariableDataHolder::ContainerVariableDataHolderPointerVariantType> CollectiveVariableDataHolder::GetVariableDataHolders()
{
    return mContainerVariableDataHolderPointersList;
}

std::vector<CollectiveVariableDataHolder::ContainerVariableDataHolderPointerVariantType> CollectiveVariableDataHolder::GetVariableDataHolders() const
{
    return mContainerVariableDataHolderPointersList;
}

bool CollectiveVariableDataHolder::IsCompatibleWith(const CollectiveVariableDataHolder& rOther) const
{
    if (mContainerVariableDataHolderPointersList.size() != rOther.mContainerVariableDataHolderPointersList.size()) {
        return false;
    }

    bool is_compatible = true;

    for (IndexType i = 0; i < mContainerVariableDataHolderPointersList.size(); ++i) {
        std::visit([&](auto&& v) {
            using v_type = std::decay_t<decltype(v)>;
            auto* ptr = std::get_if<v_type>(&(rOther.mContainerVariableDataHolderPointersList[i]));
            is_compatible = is_compatible && (ptr != nullptr);
        }, mContainerVariableDataHolderPointersList[i]);
    }

    return is_compatible;
}

std::string CollectiveVariableDataHolder::Info() const
{
    std::stringstream msg;

    msg << "CollectiveVariableDataHolder contains following data holders:\n";

    for (const auto& p_container_variable_data_holder : mContainerVariableDataHolderPointersList) {
        std::visit([&msg](auto&& v) { msg << "\t" << *v; }, p_container_variable_data_holder);
    }

    return msg.str();
}

CollectiveVariableDataHolder CollectiveVariableDataHolder::operator+(const CollectiveVariableDataHolder& rOther) const
{
    KRATOS_ERROR_IF_NOT(IsCompatibleWith(rOther))
        << "Unsupported collective variable data holders provided for \"+\" operation."
        << "\nLeft operand : " << *this << "\nRight operand: " << rOther << std::endl;

    CollectiveVariableDataHolder result(rOther);
    for (IndexType i = 0; i < mContainerVariableDataHolderPointersList.size(); ++i) {
        std::visit([&](auto&& v) {
            using v_type = std::decay_t<decltype(v)>;
            v->operator+=(*(std::get<v_type>(mContainerVariableDataHolderPointersList[i])));
        }, result.mContainerVariableDataHolderPointersList[i]);
    }

    return result;
}

CollectiveVariableDataHolder& CollectiveVariableDataHolder::operator+=(const CollectiveVariableDataHolder& rOther)
{
    KRATOS_ERROR_IF_NOT(IsCompatibleWith(rOther))
        << "Unsupported collective variable data holders provided for \"+=\" operation."
        << "\nLeft operand : " << *this << "\nRight operand: " << rOther << std::endl;

    for (IndexType i = 0; i < mContainerVariableDataHolderPointersList.size(); ++i) {
        std::visit([&](auto&& v) {
            using v_type = std::decay_t<decltype(v)>;
            v->operator+=(*(std::get<v_type>(rOther.mContainerVariableDataHolderPointersList[i])));
        }, mContainerVariableDataHolderPointersList[i]);
    }

    return *this;
}

CollectiveVariableDataHolder CollectiveVariableDataHolder::operator+(const double rOther) const
{
    CollectiveVariableDataHolder result(*this);
    for (IndexType i = 0; i < mContainerVariableDataHolderPointersList.size(); ++i) {
        std::visit([&](auto&& v) {
            v->operator+=(rOther);
        }, result.mContainerVariableDataHolderPointersList[i]);
    }

    return result;
}

CollectiveVariableDataHolder& CollectiveVariableDataHolder::operator+=(const double rOther)
{
    for (IndexType i = 0; i < mContainerVariableDataHolderPointersList.size(); ++i) {
        std::visit([&](auto&& v) {
            v->operator+=(rOther);
        }, mContainerVariableDataHolderPointersList[i]);
    }

    return *this;
}

CollectiveVariableDataHolder CollectiveVariableDataHolder::operator-(const CollectiveVariableDataHolder& rOther) const
{
    KRATOS_ERROR_IF_NOT(IsCompatibleWith(rOther))
        << "Unsupported collective variable data holders provided for \"-\" operation."
        << "\nLeft operand : " << *this << "\nRight operand: " << rOther << std::endl;

    CollectiveVariableDataHolder result(*this);
    for (IndexType i = 0; i < mContainerVariableDataHolderPointersList.size(); ++i) {
        std::visit([&](auto&& v) {
            using v_type = std::decay_t<decltype(v)>;
            v->operator-=(*(std::get<v_type>(rOther.mContainerVariableDataHolderPointersList[i])));
        }, result.mContainerVariableDataHolderPointersList[i]);
    }

    return result;
}

CollectiveVariableDataHolder& CollectiveVariableDataHolder::operator-=(const CollectiveVariableDataHolder& rOther)
{
    KRATOS_ERROR_IF_NOT(IsCompatibleWith(rOther))
        << "Unsupported collective variable data holders provided for \"-=\" operation."
        << "\nLeft operand : " << *this << "\nRight operand: " << rOther << std::endl;

    for (IndexType i = 0; i < mContainerVariableDataHolderPointersList.size(); ++i) {
        std::visit([&](auto&& v) {
            using v_type = std::decay_t<decltype(v)>;
            v->operator-=(*(std::get<v_type>(rOther.mContainerVariableDataHolderPointersList[i])));
        }, mContainerVariableDataHolderPointersList[i]);
    }

    return *this;
}

CollectiveVariableDataHolder CollectiveVariableDataHolder::operator-(const double rOther) const
{
    CollectiveVariableDataHolder result(*this);
    for (IndexType i = 0; i < mContainerVariableDataHolderPointersList.size(); ++i) {
        std::visit([&](auto&& v) {
            v->operator-=(rOther);
        }, result.mContainerVariableDataHolderPointersList[i]);
    }

    return result;
}

CollectiveVariableDataHolder& CollectiveVariableDataHolder::operator-=(const double rOther)
{
    for (IndexType i = 0; i < mContainerVariableDataHolderPointersList.size(); ++i) {
        std::visit([&](auto&& v) {
            v->operator-=(rOther);
        }, mContainerVariableDataHolderPointersList[i]);
    }

    return *this;
}

CollectiveVariableDataHolder CollectiveVariableDataHolder::operator*(const double rOther) const
{
    CollectiveVariableDataHolder result(*this);
    for (IndexType i = 0; i < mContainerVariableDataHolderPointersList.size(); ++i) {
        std::visit([&](auto&& v) {
            v->operator*=(rOther);
        }, result.mContainerVariableDataHolderPointersList[i]);
    }

    return result;
}

CollectiveVariableDataHolder& CollectiveVariableDataHolder::operator*=(const double rOther)
{
    for (IndexType i = 0; i < mContainerVariableDataHolderPointersList.size(); ++i) {
        std::visit([&](auto&& v) {
            v->operator*=(rOther);
        }, mContainerVariableDataHolderPointersList[i]);
    }

    return *this;
}

CollectiveVariableDataHolder CollectiveVariableDataHolder::operator/(const double rOther) const
{
    CollectiveVariableDataHolder result(*this);
    for (IndexType i = 0; i < mContainerVariableDataHolderPointersList.size(); ++i) {
        std::visit([&](auto&& v) {
            v->operator/=(rOther);
        }, result.mContainerVariableDataHolderPointersList[i]);
    }

    return result;
}

CollectiveVariableDataHolder& CollectiveVariableDataHolder::operator/=(const double rOther)
{
    for (IndexType i = 0; i < mContainerVariableDataHolderPointersList.size(); ++i) {
        std::visit([&](auto&& v) {
            v->operator/=(rOther);
        }, mContainerVariableDataHolderPointersList[i]);
    }

    return *this;
}

CollectiveVariableDataHolder CollectiveVariableDataHolder::operator^(const double rOther) const
{
    CollectiveVariableDataHolder result(*this);
    for (IndexType i = 0; i < mContainerVariableDataHolderPointersList.size(); ++i) {
        std::visit([&](auto&& v) {
            v->operator^=(rOther);
        }, result.mContainerVariableDataHolderPointersList[i]);
    }

    return result;
}

CollectiveVariableDataHolder& CollectiveVariableDataHolder::operator^=(const double rOther)
{
    for (IndexType i = 0; i < mContainerVariableDataHolderPointersList.size(); ++i) {
        std::visit([&](auto&& v) {
            v->operator^=(rOther);
        }, mContainerVariableDataHolderPointersList[i]);
    }

    return *this;
}

} // namespace Kratos