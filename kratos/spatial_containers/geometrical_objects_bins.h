//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

#pragma once

// System includes
#include <unordered_set>

// External includes

// Project includes
#include "geometries/bounding_box.h"
#include "geometries/point.h"
#include "spatial_containers/spatial_search_result.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

class GeometricalObject; // forward declaration, to be included in the cpp. This is needed to reduce the compilation time. Can be done as we consider the GeometricalObject as a pointer

/**
 * @class GeometricalObjectsBins
 * @ingroup KratosCore
 * @brief A bins container for 3 dimensional GeometricalObject entities.
 * @details It provides efficient search in radius and search nearest methods.
 * All of the geometries should be given at construction time. After
 * constructing the bins the geometries cannot be modified. In case of
 * any modification, the bins should be reconstructed.
 * @author Pooyan Dadvand
*/
class KRATOS_API(KRATOS_CORE) GeometricalObjectsBins
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GeometricalObjectsBins
    KRATOS_CLASS_POINTER_DEFINITION(GeometricalObjectsBins);

    /// The point type definition
    using PointType = Point;

    /// The type of geometrical object to be stored in the bins
    using ObjectType = GeometricalObject;

    /// The type of geometrical object to be stored in the bins
    using CellType = std::vector<GeometricalObject*>;
    using ResultType = SpatialSearchResult<GeometricalObject>;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief The constructor with all geometries to be stored. Please note that all of them should be available at construction time and cannot be modified after.
     * @param GeometricalObjectsBegin The begin iterator of the geometries to be stored
     * @param GeometricalObjectsEnd The end iterator of the geometries to be stored
     * @tparam TIteratorType The type of the iterator
     */
    template<typename TIteratorType>
    GeometricalObjectsBins(
        TIteratorType GeometricalObjectsBegin,
        TIteratorType GeometricalObjectsEnd,
        const double Tolerance = 1e-12
        )
    {
        mTolerance = Tolerance;
        const std::size_t number_of_objects = std::distance(GeometricalObjectsBegin, GeometricalObjectsEnd);
        if (number_of_objects > 0){
            mBoundingBox.Set(GeometricalObjectsBegin->GetGeometry().begin(), GeometricalObjectsBegin->GetGeometry().end());
            for (TIteratorType i_object = GeometricalObjectsBegin ; i_object != GeometricalObjectsEnd ; i_object++){
                mBoundingBox.Extend(i_object->GetGeometry().begin() , i_object->GetGeometry().end());
            }
            mBoundingBox.Extend(Tolerance);
        }
        CalculateCellSize(number_of_objects);
        mCells.resize(GetTotalNumberOfCells());
        AddObjectsToCells(GeometricalObjectsBegin, GeometricalObjectsEnd);
    }

    /**
     * @brief The constructor with all geometries to be stored. Please note that all of them should be available at construction time and cannot be modified after.
     * @param rGeometricalObjectsVector The geometries to be stored
     * @tparam TContainer The container type
     */
    template<typename TContainer>
    GeometricalObjectsBins(
        TContainer& rGeometricalObjectsVector,
        const double Tolerance = 1e-12)
        : GeometricalObjectsBins(rGeometricalObjectsVector.begin(), rGeometricalObjectsVector.end(), Tolerance)
    {
    }

    /// Destructor.
    virtual ~GeometricalObjectsBins(){}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Accessing cell_ijk giving the 3 indices
     * @param I The index in x direction
     * @param J The index in y direction
     * @param K The index in z direction
     * @return The cell_ijk
     */
    CellType& GetCell(
        const std::size_t I,
        const std::size_t J,
        const std::size_t K
        );

    /**
     * @brief Calculating the cell_ijk bounding box
     * @param I The index in x direction
     * @param J The index in y direction
     * @param K The index in z direction
     * @return The bounding box of the cell
     */
    BoundingBox<PointType> GetCellBoundingBox(
        const std::size_t I,
        const std::size_t J,
        const std::size_t K
        );

    /**
     * @brief This method takes a point and finds all of the objects in the given radius to it.
     * @details The result contains the object and also its distance to the point.
     * @param rPoint The point to be checked
     * @param Radius The radius to be checked
     * @param rResults The results of the search
     */
    void SearchInRadius(
        const PointType& rPoint,
        const double Radius,
        std::vector<ResultType>& rResults
        );

    /**
     * @brief This method takes a point and finds all of the objects in the given radius to it (iterative version).
     * @details The result contains the object and also its distance to the point.
     * @param itPointBegin The first point iterator
     * @param itPointEnd The last point iterator
     * @param Radius The radius to be checked
     * @param rResults The results of the search
     * @tparam TPointIteratorType The type of the point iterator
     */
    template<typename TPointIteratorType>
    void SearchInRadius(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        const double Radius,
        std::vector<std::vector<ResultType>>& rResults
        )
    {
        const std::size_t number_of_points = std::distance(itPointBegin, itPointEnd);
        rResults.resize(number_of_points);
        for (auto it_point = itPointBegin ; it_point != itPointEnd ; it_point++){
            SearchInRadius(*it_point, Radius, rResults[it_point - itPointBegin]);
        }
    }

    /**
     * @brief This method takes a point and finds the nearest object to it in a given radius.
     * @details If there are more than one object in the same minimum distance only one is returned
     * If there are no objects in that radius the result will be set to not found.
     * Result contains a flag is the object has been found or not.
     * @param rPoint The point to be checked
     * @param Radius The radius to be checked
     * @return ResultType The result of the search
     */
    ResultType SearchNearestInRadius(
        const PointType& rPoint,
        const double Radius
        );

    /**
     * @brief This method takes a point and finds the nearest object to it in a given radius (iterative version).
     * @details If there are more than one object in the same minimum distance only one is returned
     * If there are no objects in that radius the result will be set to not found.
     * Result contains a flag is the object has been found or not.
     * @param itPointBegin The first point iterator
     * @param itPointEnd The last point iterator
     * @param Radius The radius to be checked
     * @return std::vector<ResultType> The result of the search
     * @tparam TPointIteratorType The type of the point iterator
     */
    template<typename TPointIteratorType>
    std::vector<ResultType> SearchNearestInRadius(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd,
        const double Radius
        )
    {
        // Doing a vector of results
        std::vector<ResultType> results;
        const std::size_t number_of_points = std::distance(itPointBegin, itPointEnd);
        results.resize(number_of_points);
        for (auto it_point = itPointBegin ; it_point != itPointEnd ; it_point++){
            results[it_point - itPointBegin] = SearchNearestInRadius(*it_point, Radius);
        }
        return results;
    }

    /**
     * @brief This method takes a point and finds the nearest object to it.
     * @details If there are more than one object in the same minimum distance only one is returned
     * Result contains a flag is the object has been found or not.
     * @param rPoint The point to be checked
     * @return ResultType The result of the search
    */
    ResultType SearchNearest(const PointType& rPoint);

    /**
     * @brief This method takes a point and finds the nearest object to it (iterative version).
     * @details If there are more than one object in the same minimum distance only one is returned
     * Result contains a flag is the object has been found or not.
     * @param itPointBegin The first point iterator
     * @param itPointEnd The last point iterator
     * @return std::vector<ResultType> The result of the search
     * @tparam TPointIteratorType The type of the point iterator
     */
    template<typename TPointIteratorType>
    std::vector<ResultType> SearchNearest(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd
        )
    {
        // Doing a vector of results
        std::vector<ResultType> results;
        const std::size_t number_of_points = std::distance(itPointBegin, itPointEnd);
        results.resize(number_of_points);
        for (auto it_point = itPointBegin ; it_point != itPointEnd ; it_point++){
            results[it_point - itPointBegin] = SearchNearest(*it_point);
        }
        return results;
    }

    /**
     * @brief This method takes a point and search if it's inside an geometrical object of the domain.
     * @details If it is inside an object, it returns it, and search distance is set to zero.
     * If there is no object, the result will be set to not found.
     * Result contains a flag is the object has been found or not.
     * This method is a simplified and faster method of SearchNearest.
     * @param rPoint The point to be checked
     * @return ResultType The result of the search
     */
    ResultType SearchIsInside(const PointType& rPoint);

    /**
     * @brief This method takes a point and search if it's inside an geometrical object of the domain (iterative version).
     * @details If it is inside an object, it returns it, and search distance is set to zero.
     * If there is no object, the result will be set to not found.
     * Result contains a flag is the object has been found or not.
     * This method is a simplified and faster method of SearchNearest.
     * @param itPointBegin The first point iterator
     * @param itPointEnd The last point iterator
     * @return std::vector<ResultType> The result of the search
     * @tparam TPointIteratorType The type of the point iterator
     */
    template<typename TPointIteratorType>
    std::vector<ResultType> SearchIsInside(
        TPointIteratorType itPointBegin,
        TPointIteratorType itPointEnd
        )
    {
        // Doing a vector of results
        std::vector<ResultType> results;
        const std::size_t number_of_points = std::distance(itPointBegin, itPointEnd);
        results.resize(number_of_points);
        for (auto it_point = itPointBegin ; it_point != itPointEnd ; it_point++){
            results[it_point - itPointBegin] = SearchIsInside(*it_point);
        }
        return results;
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Getting the bins bounding box
     * @return The bounding box of the bins
     */
    const BoundingBox<PointType>& GetBoundingBox() const {
        return mBoundingBox;
    }

    /**
     * @brief return an array with the x,y and z size of the cube
     * @return The size of the cube
     */
    const array_1d<double, 3>& GetCellSizes(){
        return mCellSizes;
    }

    /**
     * @brief returns a 3D array having the number of cells in direction x, y and z
     * @return The number of cells in each direction
     */
    const array_1d<std::size_t, 3>& GetNumberOfCells(){
        return mNumberOfCells;
    }

    /**
     * @brief The total number of cells in the container
     * @return The total number of cells
     */
    std::size_t GetTotalNumberOfCells(){
        return mNumberOfCells[0] * mNumberOfCells[1] * mNumberOfCells[2];
    }

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "GeometricalObjectsBins" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "GeometricalObjectsBins";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

    ///@}
    ///@name Friends
    ///@{

    ///@}
protected:
    ///@name Protected Life Cycle
    ///@{

    /// Default constructor protected.
    GeometricalObjectsBins() = default;

    ///@}
    ///@name Protected Static Member Variables
    ///@{

    static constexpr unsigned int Dimension = 3;    /// The dimension of the problem

    ///@}
    ///@name Protected Member Variables
    ///@{

    BoundingBox<PointType> mBoundingBox;             /// The bounding box of the domain
    array_1d<std::size_t, Dimension> mNumberOfCells; /// The number of cells in each direction
    array_1d<double, 3> mCellSizes;                  /// The size of each cell in each direction
    array_1d<double, 3> mInverseOfCellSize;          /// The inverse of the size of each cell in each direction
    std::vector<CellType> mCells;                    /// The cells of the domain
    double mTolerance;                               /// The tolerance considered

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief This method checks if a point is inside any bounding box of the global bounding boxes
     * @param rCoords The coordinates of the point
     * @return True if the point is inside the bounding box
     */
    bool PointIsInsideBoundingBox(const array_1d<double, 3>& rCoords);

    /**
     * @brief This method checks if a point is inside any bounding box of the global bounding boxes considering a certain tolerance
     * @param rCoords The coordinates of the point
     * @param Tolerance The tolerance
     * @return True if the point is inside the bounding box
     */
    bool PointIsInsideBoundingBoxWithTolerance(
        const array_1d<double, 3>& rCoords,
        const double Tolerance
        );

    /**
     * @brief Calculate the cell sizes to be as equilateral as possible and tries to approximate (roughly) the given number of cells
     * @details This method calculates the cell sizes to be as equilateral as possible and tries to approximate (roughly) the given number of cells
     * @param NumberOfCells The number of cells to be calculated
     */
    void CalculateCellSize(const std::size_t NumberOfCells);

    /**
     * @brief Adding objects to the cells that intersecting with it.
     * @details This method takes a geometrical object and adds it to the cells that intersecting with it.
     * @tparam TIteratorType The type of the iterator of the geometrical objects
     * @param GeometricalObjectsBegin The begining of the geometrical objects
     * @param GeometricalObjectsEnd The end of the geometrical objects
     */
    template<typename TIteratorType>
    void AddObjectsToCells(
        TIteratorType GeometricalObjectsBegin,
        TIteratorType GeometricalObjectsEnd
        )
    {
        for(auto i_geometrical_object = GeometricalObjectsBegin ; i_geometrical_object != GeometricalObjectsEnd ; i_geometrical_object++){
            array_1d<std::size_t, 3> min_position(3,0);
            array_1d<std::size_t, 3> max_position(3,0);
            CalculateMinMaxPositions(i_geometrical_object->GetGeometry(), min_position, max_position);
            for(std::size_t k = min_position[2] ; k < max_position[2] ; k++){
                for(std::size_t j = min_position[1] ; j < max_position[1] ; j++){
                    for(std::size_t i = min_position[0] ; i < max_position[0] ; i++){
                        auto cell_bounding_box = GetCellBoundingBox(i,j,k);
                        if(IsIntersected(i_geometrical_object->GetGeometry(), cell_bounding_box, mTolerance)){
                            GetCell(i,j,k).push_back(&(*i_geometrical_object));
                        }
                    }
                }
            }
        }
    }

    ///@}
private:
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Giving the min and max position of cells intersecting with the bounding box of the geometry.
     * @details This method takes a geometrical object and calculates the min and max position of cells intersecting with the bounding box of the geometry.
     * @tparam TGeometryType The type of the geometrical object
     * @param rGeometry The geometrical object to be checked
     * @param rMinPosition The min position of cells intersecting with the bounding box of the geometry
     * @param rMaxPosition The max position of cells intersecting with the bounding box of the geometry
     * @tparam TGeometryType The type of the geometrical object
     */
    template<typename TGeometryType>
    void CalculateMinMaxPositions(
        const TGeometryType& rGeometry,
        array_1d<std::size_t, 3>& rMinPosition,
        array_1d<std::size_t, 3>& rMaxPosition
        )
    {
        if(rGeometry.empty())
            return;

        BoundingBox<PointType> bounding_box(rGeometry.begin(), rGeometry.end());

        for(unsigned int i = 0; i < 3; i++ ) {
            rMinPosition[i] = CalculatePosition( bounding_box.GetMinPoint()[i], i );
            rMaxPosition[i] = CalculatePosition( bounding_box.GetMaxPoint()[i], i ) + 1;
        }
    }

    /**
     * @brief Calculating the cell position of a given coordinate in respective axis given by dimension
     * @param Coordinate The coordinate to be calculated
     * @param ThisDimension The current dimension index
     * @return std::size_t The cell position of a given coordinate in respective axis given by dimension
     */
    std::size_t CalculatePosition(
        const double Coordinate,
        const int ThisDimension
        ) const;

    /**
     * @brief Expands by the tolerance the geometry bounding box and checks the intersection
     * @details This method expands by the tolerance the geometry bounding box and checks the intersection
     * @tparam TGeometryType The type of the geometrical object
     * @param rGeometry The geometrical object to be checked
     * @param rBox The bounding box of the cell
     * @param ThisTolerance The tolerance to be considered
     * @return true if the geometry bounding box intersects with the cell bounding box
     */
    template<typename TGeometryType>
    static inline bool IsIntersected(
        TGeometryType& rGeometry,
        const BoundingBox<PointType>& rBox,
        const double ThisTolerance
        )
    {
        PointType low_point_tolerance;
        PointType high_point_tolerance;

        for(unsigned int i = 0; i<3; i++) {
            low_point_tolerance[i]  =  rBox.GetMinPoint()[i] - ThisTolerance;
            high_point_tolerance[i] =  rBox.GetMaxPoint()[i] + ThisTolerance;
        }

        return rGeometry.HasIntersection(low_point_tolerance,high_point_tolerance);
    }

    /**
     * @brief Searchs in objects in the given cell for the ones with distance less than given radius + tolerance.
     * @details This method takes a cell and a point and searchs in objects in the given cell for the ones with distance less than given radius + tolerance.
     * @param rCell The cell to be checked
     * @param rPoint The point to be checked
     * @param Radius The radius to be considered
     * @param rResults The results of the search
     */
    void SearchInRadiusInCell(
        const CellType& rCell,
        const PointType& rPoint,
        const double Radius,
        std::unordered_map<GeometricalObject*, double>& rResults
        );

    /**
     * @brief Searchs in objects in the given cell for the nearest one.
     * @details This method takes a cell and a point and searchs in objects in the given cell for the nearest one.
     * @param rCell The cell to be checked
     * @param rPoint The point to be checked
     * @param rResult The result of the search
     * @param MaxRadius The max radius to be considered
     */
    void SearchNearestInCell(
        const CellType& rCell,
        const PointType& rPoint,
        ResultType& rResult,
        const double MaxRadius
        );

    /**
     * @brief Searchs in objects in the given cell for the one inside only.
     * @details This method takes a cell and a point and searchs in objects in the given cell for the one inside only.
     * @param rCell The cell to be checked
     * @param rPoint The point to be checked
     * @param rResult The result of the search
     */
    void SearchIsInsideInCell(
        const CellType& rCell,
        const PointType& rPoint,
        ResultType& rResult
        );

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator deleted.
    GeometricalObjectsBins& operator=(GeometricalObjectsBins const& rOther) = delete;

    /// Copy constructor deleted.
    GeometricalObjectsBins(GeometricalObjectsBins const& rOther) = delete;

    ///@}

}; // Class GeometricalObjectsBins

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

// /// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                 GeometricalObjectsBins& rThis){
//                     return rIStream;
//                 }

// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                 const GeometricalObjectsBins& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);

//     return rOStream;
// }

///@}

///@} addtogroup block

}  // namespace Kratos.