//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand, Daniel Diez
//                   
//


// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "containers/data_value_container.h"


namespace Kratos
{
    KRATOS_CREATE_LOCAL_FLAG(DataValueContainer, OVERWRITE_OLD_VALUES, 0);


    void DataValueContainer::Merge(
        const DataValueContainer& rOther,
        Flags Options)
    {
        const bool overwrite_values = Options.Is(OVERWRITE_OLD_VALUES);

        if (overwrite_values) {
            for (const_iterator i = rOther.mData.begin(); i != rOther.mData.end(); ++i) {
            bool variable_already_exist = false;
                for (iterator j = mData.begin(); j != mData.end(); ++j) {
                    if (i->first == j->first) {
                        variable_already_exist = true;
                        j->second = i->first->Clone(i->second);
                    }
                }
            if (!variable_already_exist) mData.push_back(ValueType(i->first, i->first->Clone(i->second)));
            }
        }

        else {
            for (const_iterator i = rOther.mData.begin(); i != rOther.mData.end(); ++i) {
            bool variable_already_exist = false;
                for (iterator j = mData.begin(); j != mData.end(); ++j) {
                    if (i->first == j->first) {
                        variable_already_exist = true;
                    }
                }
            if (!variable_already_exist) mData.push_back(ValueType(i->first, i->first->Clone(i->second)));
            }
        }
    }

}  // namespace Kratos.


