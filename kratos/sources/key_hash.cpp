//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/key_hash.h"
#include "containers/variable_data.h"

namespace Kratos
{

HashType VariableHasher::operator()(const VariableData& rVariable) const
{
    return rVariable.Key();
}

/***********************************************************************************/
/***********************************************************************************/

bool VariableComparator::operator()(
    const VariableData& rFirst,
    const VariableData& rSecond
    ) const
{
    return rFirst.Key() == rSecond.Key();
}

/***********************************************************************************/
/***********************************************************************************/

HashType pVariableHasher::operator()(const VariableData* pVariable) const
{
    return pVariable->Key();
}

/***********************************************************************************/
/***********************************************************************************/

bool pVariableComparator::operator()(
    const VariableData* pFirst,
    const VariableData* pSecond
    ) const
{
    return pFirst->Key() == pSecond->Key();
}

}  // namespace Kratos