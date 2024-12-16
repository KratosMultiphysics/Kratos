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
#include "vertex_utilities.h"

// Core includes
#include "includes/model_part.h"


namespace Kratos
{
namespace HDF5
{
namespace Detail
{


class KRATOS_API(HDF5_APPLICATION) Vertex : public Point
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(Vertex);

    Vertex(Vertex&& rOther) = default;

    Vertex(const Vertex& rOther) = default;

    /** Constructor
     *
     *  @details a search using the provided locator is used to find which element contains this vertex.
     *  If the search fails, @ref{GetValue} should not be called and will throw an exception.
     *  The search result can be checked with @ref{IsLocated}.
     *
     *  @param rPosition position attribute
     *  @param rLocator point locator exposing a FindElement member
     *  @param id vertex identifier (should be unique, but no checks are performed)
     *  @param isHistorical decides whether historical or non-historical nodal variables
     *  are used for interpolation.
     */
    Vertex(const array_1d<double,3>& rPosition,
           const PointLocatorAdaptor& rLocator,
           std::size_t id,
           bool isHistorical = true);

    /** Default constructor required by PointerVector
     *  @note this constructor does not locate the vertex and must
     *  not be called.
     */
    Vertex();

    /** Pointer constructor for python
     *
     *  @details a search using the provided locator is used to find which element contains this vertex.
     *  If the search fails, @ref{GetValue} should not be called and will throw an exception.
     *  The search result can be checked with @ref{IsLocated}.
     *
     *  @param rPosition position attribute
     *  @param rLocator point locator exposing a FindElement member
     *  @param id vertex identifier (should be unique, but no checks are performed)
     *  @param isHistorical decides whether historical or non-historical nodal variables
     */
    static Vertex::Pointer MakeShared(const array_1d<double,3>& rPosition,
                                      const PointLocatorAdaptor& rLocator,
                                      std::size_t id,
                                      bool isHistorical);

    std::size_t GetID() const;

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
    std::size_t mID;

    const Element::WeakPointer mpContainingElement;

    const NodalVariableGetter* mpVariableGetter;

    Kratos::Vector mShapeFunctionValues;

    static const NodalVariableGetter::UniquePointer mpHistoricalVariableGetter;

    static const NodalVariableGetter::UniquePointer mpNonHistoricalVariableGetter;
};


using VertexContainerType = PointerVector<Vertex,
                                          Vertex::Pointer,
                                          std::vector<Vertex::Pointer>>;


} // namespace Detail
} // namespace HDF5
} // namespace Kratos

// Template implementations
#include "vertex_impl.h"

#endif