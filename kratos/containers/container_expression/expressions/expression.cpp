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
#include <numeric>

// Project includes

// Include base h
#include "expression.h"

namespace Kratos {

std::size_t Expression::GetFlattenedShapeSize() const
{
    const auto& r_shape = this->GetShape();
    return std::accumulate(
        r_shape.begin(),
        r_shape.end(), 1UL,
        [](const auto V1, const auto V2) { return V1 * V2; });
}

} // namespace Kratos