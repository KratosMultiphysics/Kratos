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
               std::size_t id)
    : Point(rPosition),
      mID(id),
      mpContainingElement(rLocator.FindElement(*this))
{
    KRATOS_TRY

    // Get shape function values if the containing element was found
    if (mpContainingElement.get()) {
        Point local_coordinates;
        const auto& r_geometry = mpContainingElement->GetGeometry();
        r_geometry.PointLocalCoordinates(local_coordinates, *this);
        r_geometry.ShapeFunctionsValues(mShapeFunctionValues, local_coordinates);
    }

    KRATOS_CATCH("");
}


Vertex::Vertex()
    : Point(array_1d<double,3>{{0.0, 0.0, 0.0}}),
      mID(std::numeric_limits<std::size_t>::max()),
      mpContainingElement()
{
    KRATOS_ERROR << "Call to default constructor is forbidden";
}


Vertex::Pointer Vertex::MakeShared(const array_1d<double,3>& rPosition,
                                   const PointLocatorAdaptor& rLocator,
                                   std::size_t id)
{
    return Vertex::Pointer(new Vertex(rPosition, rLocator, id));
}


bool Vertex::IsLocated() const
{
    return mpContainingElement.get();
}


} // namespace Detail
} // namespace HDF5
} // namespace Kratos