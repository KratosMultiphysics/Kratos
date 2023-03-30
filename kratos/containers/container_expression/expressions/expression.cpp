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

// Project includes

// Include base h
#include "expression.h"

namespace Kratos {

std::size_t Expression::GetFlattenedSize() const
{
    IndexType result = 1;
    for (const auto v : this->GetShape()) {
        result *= v;
    }
    return result;
}

} // namespace Kratos