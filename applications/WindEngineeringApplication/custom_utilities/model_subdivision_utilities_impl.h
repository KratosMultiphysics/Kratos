//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

#ifndef KRATOS_WIND_MODEL_SUBDIVISION_UTILITIES_IMPL
#define KRATOS_WIND_MODEL_SUBDIVISION_UTILITIES_IMPL

// post-included in "model_subdivision_utilities.h"

// Core includes
#include "utilities/math_utils.h"

// STL includes
#include <mutex> // for std::lock_guard


namespace Kratos
{
namespace Wind
{


inline void ModelSubdivisionUtilities::ThreadSafeIndexSet::Push(ModelPart* pSubModelPart,
                                                                ModelSubdivisionUtilities::ThreadSafeIndexSet::IndexType value)
{
    KRATOS_TRY

    // Push index to the corresponding container of the model part pointer
    // Note: the case where pSubModelPart==NULL is valid, and means the value will not be pushed anywhere
    if (pSubModelPart) {
        auto& rPair = mIndexSets[pSubModelPart];
        std::lock_guard<LockType> lock(*rPair.second); // replace with std::scoped_lock (c++17)
        rPair.first.push_back(value);
    }

    KRATOS_CATCH("");
}


template <class TFunction>
inline void ModelSubdivisionUtilities::ThreadSafeIndexSet::Apply(TFunction function)
{
    KRATOS_TRY

    for (auto& rPair : mIndexSets) {
        (rPair.first->*function)(rPair.second.first, 0);
    }

    KRATOS_CATCH("")
}


inline bool ModelSubdivisionUtilities::Slab::IsBelow(const array_1d<double,3>& rPoint) const
{
    return mBottomPlane.IsOnPositiveSide(rPoint, mIsOpen);
}

inline bool ModelSubdivisionUtilities::Slab::IsAbove(const array_1d<double,3>& rPoint) const
{
    return mTopPlane.IsOnPositiveSide(rPoint, mIsOpen);
}

inline bool ModelSubdivisionUtilities::Slab::IsInside(const array_1d<double,3>& rPoint) const
{
    return (!IsBelow(rPoint)) && (!IsAbove(rPoint));
}


inline bool ModelSubdivisionUtilities::Slab::Plane::IsOnPositiveSide(const array_1d<double,3>& rPoint, bool open) const
{
    const double product = MathUtils<double>::Dot(rPoint - mReferencePoint, mNormal);
    if (product < 0) {
        return false;
    }
    else if (0 < product) {
        return true;
    }
    else if (open) {
        return true;
    }
    else {
        return false;
    }
}


inline const array_1d<double,3>& ModelSubdivisionUtilities::Slab::Normal() const
{
    return mTopPlane.mNormal;
}


} // namespace Wind
} // namespace Kratos

#endif