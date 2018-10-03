//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                       license: license.txt
//
//  License:          BSD License
//  Main authors:     Riccardo Rossi
//                    Pablo Becker
//                    Carlos Roig
//                    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_BINBASED_FAST_POINT_LOCATOR_INCLUDED )
#define  KRATOS_BINBASED_FAST_POINT_LOCATOR_INCLUDED

// System includes

// External includes


// Project includes
#include "includes/define.h"
#include "includes/node.h"

#include "spatial_containers/spatial_containers.h"
#include "spatial_containers/bounding_box.h"
#include "spatial_containers/cell.h"
#include "spatial_containers/bins_dynamic_objects.h"

#include "utilities/spatial_containers_configure.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    /// The size definition
    typedef std::size_t SizeType;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class BinBasedFastPointLocator
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
    ///@name Type Definitions
    ///@{

    /// The configure type
    typedef TConfigureType ConfigureType;

    /// The definition of the different containers
    typedef typename ConfigureType::PointType PointType;
    typedef typename ConfigureType::EntityType EntityType;
    typedef typename ConfigureType::ContainerType ContainerType;
    typedef typename ConfigureType::IteratorType IteratorType;
    typedef typename ConfigureType::ResultContainerType ResultContainerType;
    typedef typename ConfigureType::ResultIteratorType ResultIteratorType;

    /// The definition of the node
    typedef Node<3> NodeType;

    /// The definition of the geometry
    typedef Geometry<NodeType> GeometryType;

    /// The index definition
    typedef std::size_t IndexType;

    /// Pointer definition of BinBasedFastPointLocator
    KRATOS_CLASS_POINTER_DEFINITION(BinBasedFastPointLocator);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief This is the default constructor
     * @param rModelPart The model part of the mesh used in the search
     */
    BinBasedFastPointLocator(ModelPart& rModelPart)
        : mrModelPart(rModelPart)
    {
    }

    /// Destructor.
    virtual ~BinBasedFastPointLocator() = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Function to construct or update the search database
     */
    void UpdateSearchDatabase()
    {
        KRATOS_TRY

        // Copy the entities to a new container, as the list will be shuffled duringthe construction of the tree
        ContainerType entities_array;
        GetContainer(mrModelPart, entities_array);
        IteratorType it_begin = entities_array.begin();
        IteratorType it_end = entities_array.end();

        auto paux = typename BinsObjectDynamic<ConfigureType>::Pointer(new BinsObjectDynamic<ConfigureType > (it_begin, it_end));
        paux.swap(mpBinsObjectDynamic);

        KRATOS_CATCH("")
    }

    /**
     * @brief Function to construct or update the search database
     * @param CellSize The current size of the cell used for search
     */
    void UpdateSearchDatabaseAssignedSize(double CellSize)
    {
        KRATOS_TRY

        // Copy the entities to a new container, as the list will be shuffled duringthe construction of the tree
        ContainerType entities_array;
        GetContainer(mrModelPart, entities_array);
        IteratorType it_begin = entities_array.begin();
        IteratorType it_end = entities_array.end();

        auto paux = typename BinsObjectDynamic<ConfigureType>::Pointer(new BinsObjectDynamic<ConfigureType > (it_begin, it_end, CellSize));
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
        SizeType results_found = mpBinsObjectDynamic->SearchObjectsInCell(rCoordinates, ItResultBegin, MaxNumberOfResults);

        if (results_found > 0) {
            // Loop over the candidate entities and check if the particle falls within
            for (IndexType i = 0; i < results_found; i++) {
                GeometryType& geom = (*(ItResultBegin + i))->GetGeometry();

                // Find local position
                array_1d<double, 3> point_local_coordinates;
                Vector shape_function;
                const bool is_found = geom.IsInside(rCoordinates, point_local_coordinates, Tolerance);
                geom.ShapeFunctionsValues(shape_function, point_local_coordinates);
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
        const int results_found = mpBinsObjectDynamic->SearchObjectsInCell(rCoordinates, ItResultBegin, MaxNumberOfResults);

        if (results_found > 0) {
            // Loop over the candidate entities and check if the particle falls within
            for (IndexType i = 0; i < static_cast<IndexType>(results_found); i++) {
              
                GeometryType& geom = (*(ItResultBegin + i))->GetGeometry();

                // Find local position
                array_1d<double, 3> point_local_coordinates;
                const bool is_found = geom.IsInside(rCoordinates, point_local_coordinates, Tolerance);
                geom.ShapeFunctionsValues(rNShapeFunction, point_local_coordinates);

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

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{
    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{
    ///@}

private:

    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart; /// The model part containing the mesh for the search

    typename BinsObjectDynamic<ConfigureType>::Pointer mpBinsObjectDynamic; /// The pointer of the bins used for the search

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

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
    
    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}
};
    
} // namespace Kratos.

#endif // KRATOS_BINBASED_FAST_POINT_LOCATOR_INCLUDED  defined


