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

#ifndef KRATOS_HDF5_APPLICATION_VERTEX_H
#define KRATOS_HDF5_APPLICATION_VERTEX_H

//  Project includes
#include "nodal_variable_getter.h"

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

    /** Constructor
     *
     *  @details a search using the provided locator is used to find which element contains this vertex.
     *  If the search fails, @ref{GetValue} should not be called and will throw an exception.
     *  The search result can be checked with @ref{IsLocated}.
     *
     *  @param rPosition position attribute
     *  @param rLocator point locator exposing a FindElement member (see the notes for details)
     *
     *  @note the TLocator type must expose a member function for locating the element containing this vertex with the following signature:
     *  const Element* TLocator::FindElement(const array_1d<double,3>&)
     */
    template <class TLocator>
    Vertex(const array_1d<double,3>& rPosition,
           TLocator&& rLocator,
           bool isHistorical = true);

    /** Default constructor required by PointerVectorSet
     *  @note this constructor does not locate the vertex and should
     *  not be called.
     */
    Vertex();

    /** Interpolate the requested variable
     *  @param rVariable variable to interpolate
     *  @note throws an exception if the vertex was not located successfully upon construction
     */
    template <class TValue>
    TValue GetValue(const Variable<TValue>& rVariable) const;

    /** Check whether the vertex was located successfully on construction
     *  @details calling @ref{GetValue} or @ref{GetSolutionStepValue} will result
     *  in exceptions if IsLocated returns false.
     */
    bool IsLocated() const;

private:
    const Element* mpContainingElement;

    NodalVariableGetter::UniquePointer mpVariableGetter;
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