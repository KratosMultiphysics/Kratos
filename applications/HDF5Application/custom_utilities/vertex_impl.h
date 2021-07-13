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


template <class TLocator>
bool Vertex::Locate(TLocator&& rLocator)
{
    mShapeFunctionValues.clear();
    mEntityID = rLocator.FindElement(*this, mShapeFunctionValues);

    // Return false if the location failed
    return -1 < mEntityID; 
}


template <class TValue>
TValue Vertex::GetValue(const Variable<TValue>& rVariable) const
{
    KRATOS_TRY

    const auto& rGeometry = this->GetLocatedGeometry();
    auto value = ZeroInitialized<TValue>::Get();

    if (mIsHistorical) {
        for (std::size_t i_node=0; i_node<rGeometry.size(); ++i_node) {
            value += (TValue)(mShapeFunctionValues[i_node] * rGeometry.GetPoint(i_node).GetSolutionStepValue(rVariable));
        }
    }
    else {
        for (std::size_t i_node=0; i_node<rGeometry.size(); ++i_node) {
            value += (TValue)(mShapeFunctionValues[i_node] * rGeometry.GetPoint(i_node).GetValue(rVariable));
        }
    }

    return value;

    KRATOS_CATCH("");
}


} // namespace Detail
} // namespace HDF5
} // namespace Kratos


#endif