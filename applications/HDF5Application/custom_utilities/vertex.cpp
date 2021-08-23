//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: HDF5Application/license.txt
//
//  Main author:     Máté Kelemen
//

// Internal includes
#include "vertex.h"


namespace Kratos
{
namespace HDF5
{
namespace Detail
{


Vertex::Vertex(const array_1d<double,3>& rPosition,
               const PointLocatorAdaptor& rLocator,
               bool isHistorical)
    : Point(rPosition),
      mpContainingElement(rLocator.FindElement(*this))
{
    KRATOS_TRY

    if (isHistorical) {
        mpVariableGetter = NodalVariableGetter::UniquePointer(new HistoricalVariableGetter);
    }
    else {
        mpVariableGetter = NodalVariableGetter::UniquePointer(new NonHistoricalVariableGetter);
    }

    if (mpContainingElement.get()) {
        mpContainingElement->GetGeometry().IsInside(*this, mLocalCoordinates);
    }

    KRATOS_CATCH("");
}


Vertex::Vertex()
    : Point(array_1d<double,3>{{0.0, 0.0, 0.0}}),
      mpContainingElement()
{
}


Vertex::Pointer Vertex::MakeShared(const array_1d<double,3>& rPosition,
                                   const PointLocatorAdaptor& rLocator,
                                   bool isHistorical)
{
    return Vertex::Pointer(new Vertex(rPosition, rLocator, isHistorical));
}


bool Vertex::IsLocated() const
{
    return mpContainingElement.get();
}


} // namespace Detail
} // namespace HDF5
} // namespace Kratos