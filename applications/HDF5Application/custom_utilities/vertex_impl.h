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

#pragma once

// Application includes
#include "hdf5_application_define.h"


namespace Kratos
{
namespace HDF5
{
namespace Detail
{

inline std::size_t Vertex::GetID() const
{
    return mID;
}


template <class TValue, class TVariableGetterType>
inline TValue Vertex::GetValue(
    const Variable<TValue>& rVariable,
    const TVariableGetterType& rVariableGetter) const
{
    KRATOS_TRY

    KRATOS_ERROR_IF(!mpContainingElement.get()) << "attempt to interpolate on a non-located vertex";

    const auto& r_geometry = mpContainingElement->GetGeometry();

    char dummy{};

    TValue value = rVariableGetter.GetValue(r_geometry[0], rVariable, dummy) * mShapeFunctionValues[0];

    // Interpolate variable
    for (std::size_t i_node = 1; i_node < r_geometry.size(); ++i_node) {
        value += rVariableGetter.GetValue(r_geometry[i_node], rVariable, dummy) * mShapeFunctionValues[i_node];
    }

    return value;

    KRATOS_CATCH("");
}


} // namespace Detail
} // namespace HDF5
} // namespace Kratos