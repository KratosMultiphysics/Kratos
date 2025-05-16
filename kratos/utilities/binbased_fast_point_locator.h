//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Pablo Becker
//                   Carlos Roig
//                   Vicente Mataix Ferrandiz
//

#pragma once

// Project includes
#include "includes/node.h"
#include "includes/condition.h"
#include "spatial_containers/spatial_containers.h"
#include "spatial_containers/cell.h"
#include "spatial_containers/bins_dynamic_objects.h"

#include "utilities/spatial_containers_configure.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @ingroup KratosCore
 * @brief This class is designed to allow the fast location of MANY points on the top of a 3D mesh.
 * @details The utility relies on the creation of a Bin of objects that allows finding quikly a reduced number of element candidates for the location of a point.
 * The basic idea is to allow finding the element in which a given spatial position sits
 * The user should call the function "UpdateSearchDatabase" to mount the bin and subsequently locate the points as needed
 * @author  Riccardo Rossi <rrossi@cimne.upc.edu>
 * @note The location function is threadsafe, and can be used in OpenMP loops
 * @tparam TDim If we work in a 2D or 3D space
 * @tparam TConfigureType The spatial container
 */
template< SizeType TDim, class TConfigureType = SpatialContainersConfigure<TDim> >
class BinBasedFastPointLocator
{
public:
    /// The configure type
    typedef TConfigureType ConfigureType;

    /// The definition of the different containers
    typedef typename ConfigureType::PointType PointType;
    typedef typename ConfigureType::EntityType EntityType;
    typedef typename ConfigureType::ContainerType ContainerType;
    typedef typename ConfigureType::IteratorType IteratorType;
    typedef typename ConfigureType::ResultContainerType ResultContainerType;
    typedef typename ConfigureType::ResultIteratorType ResultIteratorType;

    // The definition of the bins
    typedef BinsObjectDynamic<ConfigureType> BinsType;
    typedef typename BinsObjectDynamic<ConfigureType>::CoordinateType BinsCoordinateType;
    typedef typename BinsObjectDynamic<ConfigureType>::PointType BinsPointType;

    /// The definition of the geometry
    typedef Geometry<Node> GeometryType;

    /// The size definition
    typedef std::size_t SizeType;

    /// The index definition
    typedef std::size_t IndexType;

    /// Pointer definition of BinBasedFastPointLocator
    KRATOS_CLASS_POINTER_DEFINITION(BinBasedFastPointLocator);

    /**
     * @brief This is the default constructor
     * @param rModelPart The model part of the mesh used in the search
     */
    explicit BinBasedFastPointLocator(ModelPart& rModelPart)
        : mrModelPart(rModelPart)
    {
    }

    /// Destructor.
    virtual ~BinBasedFastPointLocator() = default;

    /// Copy constructor.
    BinBasedFastPointLocator(BinBasedFastPointLocator const& rOther)
        : mrModelPart(rOther.mrModelPart)
    {
        auto paux = typename BinsType::Pointer(new BinsType(*rOther.mpBinsObjectDynamic));
        paux.swap(mpBinsObjectDynamic);
    }

    /**
     * @brief Function to construct or update the search database
     */
    void UpdateSearchDatabase()
    {
        KRATOS_TRY

        // Copy the entities to a new container, as the list will be shuffled during the construction of the tree
        ContainerType entities_array;
        GetContainer(mrModelPart, entities_array);
        IteratorType it_begin = entities_array.begin();
        IteratorType it_end = entities_array.end();

        auto paux = typename BinsType::Pointer(new BinsType(it_begin, it_end));
        paux.swap(mpBinsObjectDynamic);

        KRATOS_CATCH("")
    }

    /**
     * @brief Function to construct or update the search database
     * @param CellSize The current size of the cell used for search
     */
    void UpdateSearchDatabaseAssignedSize(const BinsCoordinateType CellSize)
    {
        KRATOS_TRY

        // Copy the entities to a new container, as the list will be shuffled duringthe construction of the tree
        ContainerType entities_array;
        GetContainer(mrModelPart, entities_array);
        IteratorType it_begin = entities_array.begin();
        IteratorType it_end = entities_array.end();

        auto paux = typename BinsType::Pointer(new BinsType(it_begin, it_end, CellSize));
        paux.swap(mpBinsObjectDynamic);

        KRATOS_CATCH("")
    }

    /**
     * @brief This function should find the element into which a given node is located
     * and return a pointer to the element and the vector containing the
     * shape functions that define the postion within the element
     * @param rCoordinates The vector containign the coordinates of the point to be searched
     * @param rNShapeFunction The vector containing the shape function of the located point
     * @param pEntity The pointer to the element containing the located point
     * @param ItResultBegin The iterator of the search
     * @param MaxNumberOfResults The max number of results to be considered
     * @param Tolerance The tolerance considered on the search
     * @return If "false" is devolved the element is not found
     * @note this function is threadsafe and can be used within OpenMP loops
     * @warning This is legacy version (using array instead of vector for shape function)
     */
    KRATOS_DEPRECATED_MESSAGE("This is legacy version (using array instead of vector for shape function)") bool FindPointOnMesh(
        const array_1d<double, 3 >& rCoordinates,
        array_1d<double, TDim + 1 >& rNShapeFunction,
        typename EntityType::Pointer& pEntity,
        ResultIteratorType ItResultBegin,
        const SizeType MaxNumberOfResults = 1000,
        const double Tolerance = 1.0e-5
        )
    {
        // Ask to the container for the list of candidate entities
        SizeType results_found = mpBinsObjectDynamic->SearchObjectsInCell(BinsPointType{rCoordinates}, ItResultBegin, MaxNumberOfResults);

        if (results_found > 0) {
            // Loop over the candidate entities and check if the particle falls within
            for (IndexType i = 0; i < results_found; i++) {
                GeometryType& r_geom = (*(ItResultBegin + i))->GetGeometry();

                // Find local position
                array_1d<double, 3> point_local_coordinates;
                Vector shape_function;
                const bool is_found = LocalIsInside(r_geom, rCoordinates, point_local_coordinates, Tolerance);
                r_geom.ShapeFunctionsValues(shape_function, point_local_coordinates);
                noalias(rNShapeFunction) = shape_function;

                if (is_found) {
                    pEntity = (*(ItResultBegin + i));
                    return true;
                }
            }
        }

        // Not found case
        pEntity = nullptr;
        return false;
    }

    /**
     * @brief This function should find the element into which a given node is located
     * and return a pointer to the element and the vector containing the
     * shape functions that define the postion within the element
     * @param rCoordinates The vector containign the coordinates of the point to be searched
     * @param rNShapeFunction The vector containing the shape function of the located point
     * @param pEntity The pointer to the element containing the located point
     * @param ItResultBegin The iterator of the search
     * @param MaxNumberOfResults The max number of results to be considered
     * @param Tolerance The tolerance considered on the search
     * @return If "false" is devolved the element is not found
     * @note this function is threadsafe and can be used within OpenMP loops
     */
    bool FindPointOnMesh(
        const array_1d<double, 3 >& rCoordinates,
        Vector& rNShapeFunction,
        typename EntityType::Pointer& pEntity,
        ResultIteratorType ItResultBegin,
        const SizeType MaxNumberOfResults = 1000,
        const double Tolerance = 1.0e-5
        )
    {
        // Ask to the container for the list of candidate entities
        const SizeType results_found = mpBinsObjectDynamic->SearchObjectsInCell(BinsPointType{rCoordinates}, ItResultBegin, MaxNumberOfResults);

        if (results_found > 0) {
            // Loop over the candidate entities and check if the particle falls within
            for (IndexType i = 0; i < static_cast<IndexType>(results_found); i++) {

                GeometryType& r_geom = (*(ItResultBegin + i))->GetGeometry();

                // Find local position
                array_1d<double, 3> point_local_coordinates;
                const bool is_found = LocalIsInside(r_geom, rCoordinates, point_local_coordinates, Tolerance);
                r_geom.ShapeFunctionsValues(rNShapeFunction, point_local_coordinates);

                if (is_found) {
                    pEntity = (*(ItResultBegin + i));
                    return true;
                }
            }
        }

        // Not found case
        pEntity = nullptr;
        return false;
    }

    /**
     * @brief Simplified (less efficient) function to find the element into which a given node is located and return a pointer to the element and the vector containing the shape functions that define the postion within the element
     * @param rCoordinates The vector containign the coordinates of the point to be searched
     * @param rNShapeFunction The vector containing the shape function of the located point
     * @param pEntity The pointer to the element containing the located point
     * @param MaxNumberOfResults The max number of results to be considered
     * @param Tolerance The tolerance considered on the search
     * @return If "false" is devolved the element is not found
     * @note this function is threadsafe and can be used within OpenMP loops
     */
    bool FindPointOnMeshSimplified(
        const array_1d<double, 3 >& rCoordinates,
        Vector& rNShapeFunction,
        typename EntityType::Pointer& pEntity,
        const SizeType MaxNumberOfResults = 1000,
        const double Tolerance = 1.0e-5
        )
    {
        ResultContainerType results(MaxNumberOfResults);

        const bool is_found = FindPointOnMesh(rCoordinates, rNShapeFunction, pEntity, results.begin(), MaxNumberOfResults, Tolerance);

        return is_found;
    }

protected:
    /**
    * @brief Checks if given point in global space coordinates
    *        is inside the geometry boundaries. This function
    *        computes the local coordinates and checks then if
    *        this point lays within the boundaries.
    * @param rPointGlobalCoordinates the global coordinates of the
    *        external point.
    * @param rResult the local coordinates of the point.
    * @param Tolerance the tolerance to the boundary.
    * @return true if the point is inside, false otherwise
    */
    virtual bool LocalIsInside(
        const GeometryType& rGeometry,
        const GeometryType::CoordinatesArrayType& rPointGlobalCoordinates,
        GeometryType::CoordinatesArrayType& rResult,
        const double Tolerance = std::numeric_limits<double>::epsilon()
        ) const
    {
        return rGeometry.IsInside(rPointGlobalCoordinates, rResult, Tolerance);
    }

private:
    ModelPart& mrModelPart; /// The model part containing the mesh for the search

    typename BinsType::Pointer mpBinsObjectDynamic; /// The pointer of the bins used for the search

    /**
     * @brief This operation is defined to the the corresponding container type
     * @param rModelPart The model part to get the element container
     * @param The corresponding element array
     */
    static inline void GetContainer(
        ModelPart& rModelPart,
        PointerVectorSet<Element, IndexedObject>::ContainerType& rContainerArray
        )
    {
        rContainerArray = rModelPart.ElementsArray();
    }

    /**
     * @brief This operation is defined to the the corresponding container type
     * @param rModelPart The model part to get the condition container
     * @param The corresponding condition array
     */
    static inline void GetContainer(
        ModelPart& rModelPart,
        PointerVectorSet<Condition, IndexedObject>::ContainerType& rContainerArray
        )
    {
        rContainerArray = rModelPart.ConditionsArray();
    }
};

} // namespace Kratos.
