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


Vertex::Vertex()
    : Point(array_1d<double,3>{{0.0, 0.0, 0.0}}),
      mpContainingElement()
{
}


bool Vertex::IsLocated() const
{
    return mpContainingElement.get();
}


} // namespace Detail
} // namespace HDF5
} // namespace Kratos