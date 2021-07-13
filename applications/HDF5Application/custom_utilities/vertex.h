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

#ifndef KRATOS_HDF5_APPLICATION_VERTEX_H
#define KRATOS_HDF5_APPLICATION_VERTEX_H

// Core includes
#include "includes/model_part.h"


namespace Kratos
{
namespace HDF5
{
namespace Detail
{


class Vertex : public Point
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(Vertex);

    Vertex(const array_1d<double,3>& rPosition,
           const ModelPart& rModelPart,
           bool isHistorical = true);

    // Required by PointerVectorSet
    Vertex();

    /** Locate the point and update member variables
     *  @details the locator must expose TLocator::FindElement(const Point&, Vector&, ...)
     */
    template <class TLocator>
    bool Locate(TLocator&& rLocator);

    /// Interpolate the requested Variable
    template <class TValue>
    TValue GetValue(const Variable<TValue>& rVariable) const;

private:
    /// Check whether @ref{Vertex::Locate} was called and successful and return the located Element's geometry
    const Element::GeometryType& GetLocatedGeometry() const;

    const ModelPart* mpModelPart;

    int mEntityID;

    const bool mIsHistorical;

    Kratos::Vector mShapeFunctionValues;
};


using VertexContainerType = PointerVectorSet<Vertex,
                                             IndexedObject,
                                             std::less<typename IndexedObject::result_type>,
                                             std::equal_to<typename IndexedObject::result_type>,
                                             Vertex::Pointer,
                                             std::vector<Vertex::Pointer>>;


} // namespace Detail
} // namespace HDF5
} // namespace Kratos

// Template implementations
#include "vertex_impl.h"

#endif