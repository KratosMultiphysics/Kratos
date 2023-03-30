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

void intrusive_ptr_add_ref(const Expression* x)
{
    x->mReferenceCounter.fetch_add(1, std::memory_order_relaxed);
}

void intrusive_ptr_release(const Expression* x)
{
    if (x->mReferenceCounter.fetch_sub(1, std::memory_order_release) == 1) {
        std::atomic_thread_fence(std::memory_order_acquire);
        delete x;
    }
}

} // namespace Kratos