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


// STL includes
#include <limits>


namespace Kratos
{
namespace HDF5
{
namespace Detail
{


Vertex::Vertex(const array_1d<double,3>& rPosition,
               const PointLocatorAdaptor& rLocator,
               std::size_t id,
               bool isHistorical)
    : Point(rPosition),
      mID(id),
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
      mID(std::numeric_limits<std::size_t>::max()),
      mpContainingElement()
{
}


Vertex::Pointer Vertex::MakeShared(const array_1d<double,3>& rPosition,
                                   const PointLocatorAdaptor& rLocator,
                                   std::size_t id,
                                   bool isHistorical)
{
    return Vertex::Pointer(new Vertex(rPosition, rLocator, id, isHistorical));
}


bool Vertex::IsLocated() const
{
    return mpContainingElement.get();
}


} // namespace Detail
} // namespace HDF5
} // namespace Kratos