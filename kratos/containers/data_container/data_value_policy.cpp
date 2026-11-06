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
//  Collaborator:    Vicente Mataix Ferrandiz
//
//

// System includes

// External includes

// Project includes
#include "containers/data_container/data_value_policy.h"

namespace Kratos
{

DataValuePolicyBase::DataValuePolicyBase(std::size_t SizeOfData) : mSizeOfData(SizeOfData)
{
    KRATOS_DEBUG_ERROR_IF_NOT(SizeOfData > 0) << "Size must be greater than zero." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t DataValuePolicyBase::SizeOfData() const
{
    return mSizeOfData;
}

/***********************************************************************************/
/***********************************************************************************/

bool DataValuePolicyBase::IsSparse() const
{
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

DataAccessor<int> DataValuePolicyBase::GetSparseIndexAccessor() const
{
    KRATOS_ERROR << "Calling base class DataValuePolicyBase::GetSparseIndexAccessor" << std::endl;
}

} // namespace Kratos
