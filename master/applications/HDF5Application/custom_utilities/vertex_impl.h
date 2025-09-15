//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: HDF5Application/license.txt
//
//  Main author:     Máté Kelemen
//

#ifndef KRATOS_HDF5_APPLICATION_VERTEX_IMPL_H
#define KRATOS_HDF5_APPLICATION_VERTEX_IMPL_H

// Application includes
#include "hdf5_application_define.h"


namespace Kratos
{
namespace HDF5
{
namespace Detail
{


template <class TValue>
struct ZeroInitialized
{
    static TValue Get()
    {
        return TValue(0);
    }
};


template <class TValue, std::size_t ArraySize>
struct ZeroInitialized<array_1d<TValue, ArraySize>>
{
    static array_1d<TValue, ArraySize> Get()
    {
        return ZeroVector(ArraySize);
    }
};


template <class TValue>
struct ZeroInitialized<Matrix<TValue>>
{
    static Matrix<TValue> Get()
    {
        return ZeroMatrix();
    }
};


inline std::size_t Vertex::GetID() const
{
    return mID;
}


template <class TValue>
inline TValue Vertex::GetValue(const Variable<TValue>& rVariable) const
{
    KRATOS_TRY

    KRATOS_ERROR_IF(!mpContainingElement.get()) << "attempt to interpolate on a non-located vertex";

    const auto& r_geometry = mpContainingElement->GetGeometry();
    auto value = ZeroInitialized<TValue>::Get();

    // Interpolate variable
    for (std::size_t i_node=0; i_node<r_geometry.size(); ++i_node) {
        value += mShapeFunctionValues[i_node] * mpVariableGetter->GetValue(r_geometry.GetPoint(i_node), rVariable);
    }

    return value;

    KRATOS_CATCH("");
}


} // namespace Detail
} // namespace HDF5
} // namespace Kratos


#endif