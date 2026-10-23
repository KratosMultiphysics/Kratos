//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Daniel Diez
//

// System includes
#include <algorithm>
#include <vector>

// External includes

// Project includes
#include "containers/data_value_container.h"

namespace Kratos {

KRATOS_CREATE_LOCAL_FLAG(DataValueContainer, OVERWRITE_OLD_VALUES, 0);


DataValueContainer& DataValueContainer::operator=(const DataValueContainer& rOther)
{
    if (this != &rOther) {
        Clear();

        for(const_iterator i = rOther.mData.begin() ; i != rOther.mData.end() ; ++i)
            mData.push_back(ValueType(i->first, i->first->Clone(i->second)));
    } // if this != &rOther

    return *this;
}


void DataValueContainer::Merge(
    const DataValueContainer& rOther,
    const Flags Options
    )
{
    const bool overwrite_values = Options.Is(OVERWRITE_OLD_VALUES);

    if (overwrite_values) {
        for (const_iterator i = rOther.mData.begin(); i != rOther.mData.end(); ++i) {
            bool variable_already_exist = false;
                for (iterator j = mData.begin(); j != mData.end(); ++j) {
                    if (i->first == j->first) {
                        variable_already_exist = true;
                        j->first->Delete(j->second);
                        j->second = i->first->Clone(i->second);
                    }
                }
            if (!variable_already_exist) {
                mData.push_back(ValueType(i->first, i->first->Clone(i->second)));
            }
        }
    } else {
        for (const_iterator i = rOther.mData.begin(); i != rOther.mData.end(); ++i) {
            bool variable_already_exist = false;
                for (iterator j = mData.begin(); j != mData.end(); ++j) {
                    if (i->first == j->first) {
                        variable_already_exist = true;
                    }
                }
            if (!variable_already_exist) {
                mData.push_back(ValueType(i->first, i->first->Clone(i->second)));
            }
        }
    }
}

void DataValueContainer::save(Serializer& rSerializer) const
{
    KRATOS_TRY

    std::size_t size = mData.size();
    rSerializer.save("Size", size);
    for (std::size_t i = 0; i < size; i++) {
        rSerializer.save("Variable Name", mData[i].first->Name());
        mData[i].first->Save(rSerializer, mData[i].second);
    }

    KRATOS_CATCH("")
}

void DataValueContainer::load(Serializer& rSerializer)
{
    KRATOS_TRY

    std::string name;
    if (!rSerializer.IsDataOnly()) {
        // an in-place load allocates every value anew, so delete the current ones first
        Clear();
        std::size_t size;
        rSerializer.load("Size", size);
        mData.resize(size);
        for (std::size_t i = 0; i < size; i++) {
            rSerializer.load("Variable Name", name);
            mData[i].first = KratosComponents<VariableData>::pGet(name);
            mData[i].first->Allocate(&(mData[i].second));
            mData[i].first->Load(rSerializer, mData[i].second);
        }
    } else {
        // DataOnly loads into the existing values, so they are matched by variable (not by position, which
        // may differ or hold another type); values not in the archive are deleted to restore the saved state
        std::size_t size;
        rSerializer.load("Size", size);
        std::vector<bool> is_loaded(mData.size(), false);
        for (std::size_t i = 0; i < size; i++) {
            rSerializer.load("Variable Name", name);
            const auto p_variable = KratosComponents<VariableData>::pGet(name);
            const auto it = std::find_if(mData.begin(), mData.end(), IndexCheck(p_variable->SourceKey()));
            KRATOS_ERROR_IF(it == mData.end())
                << "DataOnly loading restores values into an existing container, but the variable " << name
                << " is not in it.\n";
            p_variable->Load(rSerializer, it->second);
            is_loaded[it - mData.begin()] = true;
        }

        std::size_t number_of_kept_values = 0;
        for (std::size_t i = 0; i < mData.size(); ++i) {
            if (is_loaded[i]) {
                mData[number_of_kept_values++] = mData[i];
            } else {
                mData[i].first->Delete(mData[i].second);
            }
        }
        mData.resize(number_of_kept_values);
    }

    KRATOS_CATCH("")
}

}  // namespace Kratos.