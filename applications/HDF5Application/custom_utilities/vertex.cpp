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

// Internal includes
#include "vertex.h"


namespace Kratos
{
namespace HDF5
{
namespace Detail
{


Vertex::Vertex(const array_1d<double,3>& rPosition,
               const ModelPart& rModelPart,
               bool isHistorical)
    : Point(rPosition),
      mpModelPart(&rModelPart),
      mEntityID(-1),
      mIsHistorical(isHistorical)
{
}


Vertex::Vertex()
    : Point(array_1d<double,3>{{0.0, 0.0, 0.0}}),
      mpModelPart(nullptr),
      mEntityID(-1),
      mIsHistorical(true)
{
}


const Element::GeometryType& Vertex::GetLocatedGeometry() const
{
    KRATOS_TRY

    // Check whether the location was performed and successful
    KRATOS_ERROR_IF(mEntityID < 0)
        << "Vertex was not located!" 
        << " Possible reasons are: Vertex::Locate was not called, or failed";

    // Get the geometry of the located entity for finding its associated nodes
    const auto& rGeometry = mpModelPart->GetElement(mEntityID).GetGeometry();
    const std::size_t number_of_nodes = rGeometry.size(); 

    // Check number of nodes and shape functions
    KRATOS_ERROR_IF(number_of_nodes != mShapeFunctionValues.size())
        << "Number of nodes in the located entity (" << number_of_nodes << ") "
        << "does not match the number of shape function values (" << mShapeFunctionValues.size() << ")";

    return rGeometry;

    KRATOS_CATCH("");
}


} // namespace Detail
} // namespace HDF5
} // namespace Kratos