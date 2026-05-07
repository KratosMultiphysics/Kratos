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

    std::size_t size;
    rSerializer.load("Size", size);
    mData.resize(size);
    std::string name;
    for (std::size_t i = 0; i < size; i++) {
        rSerializer.load("Variable Name", name);
        mData[i].first = KratosComponents<VariableData>::pGet(name);
        if (!rSerializer.IsDataOnly()) {
            mData[i].first->Allocate(&(mData[i].second));
        }
        mData[i].first->Load(rSerializer, mData[i].second);
    }

    KRATOS_CATCH("")
}

}  // namespace Kratos.