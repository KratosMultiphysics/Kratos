//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//


#if !defined(KRATOS_GEOMETRICAL_OBJECTS_BINS_H_INCLUDED )
#define  KRATOS_GEOMETRICAL_OBJECTS_BINS_H_INCLUDED


// System includes
#include <string>
#include <iostream>
#include <unordered_set>


// External includes


// Project includes
#include "includes/define.h"
#include "geometries/bounding_box.h"
#include "geometries/point.h"
#include "includes/geometrical_object.h"
#include "includes/global_pointer.h"
#include "utilities/geometry_utilities.h"
#include "spatial_containers/spatial_search_result.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/// A bins container for 3 dimensional GeometricalObject entities.
/** It provides efficent search in radius and search nearest methods.
 * All of the geometries should be given at construction time. After 
 * constructing the bins the geometries cannot be modified. In case of
 * any modification, the bins should be reconstructed.  
 * Please note that the current implementation is only for triangles.
 * Addapt it to other geometries is postponed to after adding distance
 * method to the geometry.
*/
class GeometricalObjectsBins
{
    static constexpr int Dimension = 3;
    static constexpr double Tolerance = 1e-12;
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GeometricalObjectsBins
    KRATOS_CLASS_POINTER_DEFINITION(GeometricalObjectsBins);

    using CellType = std::vector<GeometricalObject*>;
    using ResultType = SpatialSearchResult<GeometricalObject>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor deleted.
    GeometricalObjectsBins() = delete;

    /// The constructor with all geometries to be stored. Please note that all of them should be available at construction time and cannot be modified after.
    template<typename TIteratorType>
    GeometricalObjectsBins(TIteratorType GeometricalObjectsBegin, TIteratorType GeometricalObjectsEnd) {
        std::size_t number_of_objects = std::distance(GeometricalObjectsBegin, GeometricalObjectsEnd);
        if(number_of_objects > 0){
            mBoundingBox.Set(GeometricalObjectsBegin->GetGeometry().begin(), GeometricalObjectsBegin->GetGeometry().end());
            for(TIteratorType i_object = GeometricalObjectsBegin ; i_object != GeometricalObjectsEnd ; i_object++){
                mBoundingBox.Extend(i_object->GetGeometry().begin() , i_object->GetGeometry().end());
            }
        }
        mBoundingBox.Extend(Tolerance);
        CalculateCellSize(number_of_objects);
        mCells.resize(GetTotalNumberOfCells());
        AddObjectsToCells(GeometricalObjectsBegin, GeometricalObjectsEnd);
    }


    /// Destructor.
    virtual ~GeometricalObjectsBins(){}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Access
    ///@{

    /// Getting the bins bounding box
    const BoundingBox<Point>& GetBoundingBox() const {
        return mBoundingBox;
    }

    /// return an array with the x,y and z size of the cube
    const array_1d<double, 3>& GetCellSizes(){
        return mCellSizes;
    }

    /// returns a 3D array having the number of cells in direction x,y and z
    const array_1d<std::size_t, 3>& GetNumberOfCells(){
        return mNumberOfCells;
    }

    /// The total number of cells in the container
    std::size_t GetTotalNumberOfCells(){
        return mNumberOfCells[0] * mNumberOfCells[1] * mNumberOfCells[2];
    }

    /// Accessing cell_ijk giving the 3 indices 
    CellType& GetCell(std::size_t I, std::size_t J, std::size_t K){
        KRATOS_DEBUG_ERROR_IF(I > mNumberOfCells[0]) << "Index " << I << " is larger than number of cells in x direction : " << mNumberOfCells[0] << std::endl;
        KRATOS_DEBUG_ERROR_IF(J > mNumberOfCells[1]) << "Index " << J << " is larger than number of cells in y direction : " << mNumberOfCells[1] << std::endl;
        KRATOS_DEBUG_ERROR_IF(K > mNumberOfCells[2]) << "Index " << K << " is larger than number of cells in z direction : " << mNumberOfCells[2] << std::endl;

        const std::size_t index = I + J * mNumberOfCells[0] + K * mNumberOfCells[1] * mNumberOfCells[0];
        return mCells[index];
    }

    /// Calculating the cell_ijk bounding box
    BoundingBox<Point> GetCellBoundingBox(std::size_t I, std::size_t J, std::size_t K){
        KRATOS_DEBUG_ERROR_IF(I > mNumberOfCells[0]) << "Index " << I << " is larger than number of cells in x direction : " << mNumberOfCells[0] << std::endl;
        KRATOS_DEBUG_ERROR_IF(J > mNumberOfCells[1]) << "Index " << J << " is larger than number of cells in y direction : " << mNumberOfCells[1] << std::endl;
        KRATOS_DEBUG_ERROR_IF(K > mNumberOfCells[2]) << "Index " << K << " is larger than number of cells in z direction : " << mNumberOfCells[2] << std::endl;

        BoundingBox<Point> result;
        
        result.GetMinPoint()[0] = mBoundingBox.GetMinPoint()[0] + I * mCellSizes[0];
        result.GetMinPoint()[1] = mBoundingBox.GetMinPoint()[1] + J * mCellSizes[1];
        result.GetMinPoint()[2] = mBoundingBox.GetMinPoint()[2] + K * mCellSizes[2];
         
        result.GetMaxPoint()[0] = mBoundingBox.GetMinPoint()[0] + (I + 1) * mCellSizes[0];
        result.GetMaxPoint()[1] = mBoundingBox.GetMinPoint()[1] + (J + 1) * mCellSizes[1];
        result.GetMaxPoint()[2] = mBoundingBox.GetMinPoint()[2] + (K + 1) * mCellSizes[2];

        return result;
    }

    /** This method takes a point and finds all of the objects in the given radius to it.
     * The result contains the object and also its distance to the point.
    */
    template<typename TPointType>
    void SearchInRadius(TPointType const& ThePoint, double Radius, std::vector<ResultType>& rResults) {
        std::unordered_set<GeometricalObject*> results;

        array_1d< std::size_t, 3 > min_position;
        array_1d< std::size_t, 3 > max_position;

        for(int i = 0; i < 3; i++ ) {
            min_position[i] = CalculatePosition(ThePoint[i] - Radius, i);
            max_position[i] = CalculatePosition(ThePoint[i] + Radius, i) + 1;
        }
        for(std::size_t k = min_position[2] ; k < max_position[2] ; k++){
            for(std::size_t j = min_position[1] ; j < max_position[1] ; j++){
                for(std::size_t i = min_position[0] ; i < max_position[0] ; i++){
                    auto& cell = GetCell(i,j,k);
                    SearchInRadiusInCell(cell, ThePoint, Radius, results);
                }
            }
        }

        rResults.clear();
        for(auto p_object : results){
            rResults.push_back(ResultType(p_object));
        }
    }

    /** This method takes a point and finds the nearest object to it in a given radius.
     * If there are more than one object in the same minimum distance only one is returned
     * If there are no objects in that radius the result will be set to not found.
     * Result contains a flag is the object has been found or not. 
    */
     template<typename TPointType>
    ResultType SearchNearestInRadius(TPointType const& ThePoint, double Radius) {
        ResultType current_result;
        current_result.SetDistance(std::numeric_limits<double>::max());

        double radius_increment = *std::max_element(mCellSizes.begin(), mCellSizes.end());

        array_1d< std::size_t, 3 > min_position;
        array_1d< std::size_t, 3 > max_position;

        for(double current_radius = 0 ; current_radius < Radius + radius_increment ; current_radius += radius_increment){
            current_radius = (current_radius > Radius) ? Radius : current_radius;
            for(int i = 0; i < 3; i++ ) {
                min_position[i] = CalculatePosition(ThePoint[i] - current_radius, i);
                max_position[i] = CalculatePosition(ThePoint[i] + current_radius, i) + 1;
            }

            for(std::size_t k = min_position[2] ; k < max_position[2] ; k++){
                for(std::size_t j = min_position[1] ; j < max_position[1] ; j++){
                    for(std::size_t i = min_position[0] ; i < max_position[0] ; i++){
                        auto& cell = GetCell(i,j,k);
                        SearchNearestInCell(cell, ThePoint, current_result, current_radius);
                    }
                }
            }
            bool all_cells_are_covered = (min_position[0] == 0) && (min_position[1] == 0) && (min_position[2] == 0);
            all_cells_are_covered &= (max_position[0] == mNumberOfCells[0] - 1) && (max_position[1] == mNumberOfCells[1] - 1) && (max_position[2] == mNumberOfCells[2] - 1);

            if(all_cells_are_covered || current_result.IsObjectFound()){
                break;
            }
        }
        return current_result;
    }

    /** This method takes a point and finds the nearest object to it.
     * If there are more than one object in the same minimum distance only one is returned
     * Result contains a flag is the object has been found or not. 
    */
    template<typename TPointType>
    ResultType SearchNearest(TPointType const& ThePoint) {
        ResultType current_result;

        array_1d<double, 3> box_size = mBoundingBox.GetMaxPoint() - mBoundingBox.GetMinPoint();
        double max_radius= *std::max_element(box_size.begin(), box_size.end());

        return SearchNearestInRadius(ThePoint, max_radius);
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

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{
    
    BoundingBox<Point> mBoundingBox;
    array_1d<std::size_t, Dimension> mNumberOfCells;
    array_1d<double, 3>  mCellSizes;
    array_1d<double, 3>  mInverseOfCellSize;
    std::vector<CellType> mCells;


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /// Caclculate the cell sizes to be as equilateral as possible and tries to approximate (roughly) the given number of cells
    void CalculateCellSize(std::size_t NumberOfCells) {
        std::size_t avarage_number_of_cells = static_cast<std::size_t>(std::pow(static_cast<double>(NumberOfCells), 1.00 / Dimension));
        std::array<double, 3> lengths;
        double avarage_length = 0.00;
        for (int i = 0; i < Dimension; i++) {
            lengths[i] = mBoundingBox.GetMaxPoint()[i] - mBoundingBox.GetMinPoint()[i];
            avarage_length += lengths[i];
        }
        avarage_length *= 1.00 / 3.00;

        if (avarage_length < std::numeric_limits<double>::epsilon()) {
            mNumberOfCells = ScalarVector(3, 1);
            return;
        }

        for (int i = 0; i < Dimension; i++) {
            mNumberOfCells[i] = static_cast<std::size_t>(lengths[i] / avarage_length * avarage_number_of_cells) + 1;
            if (mNumberOfCells[i] > 1)
                mCellSizes[i] = lengths[i] / mNumberOfCells[i];
            else
                mCellSizes[i] = avarage_length;

            mInverseOfCellSize[i] = 1.00 / mCellSizes[i];
        }

    }

    /// Adding objects to the cells that intersecting with it.
    template<typename TIteratorType>
    void AddObjectsToCells(TIteratorType GeometricalObjectsBegin, TIteratorType GeometricalObjectsEnd) {
        for(auto i_geometrical_object = GeometricalObjectsBegin ; i_geometrical_object != GeometricalObjectsEnd ; i_geometrical_object++){
            array_1d<std::size_t, 3> min_position(3,0);
            array_1d<std::size_t, 3> max_position(3,0);
            CalculateMinMaxPositions(i_geometrical_object->GetGeometry(), min_position, max_position);
            for(std::size_t k = min_position[2] ; k < max_position[2] ; k++){
                for(std::size_t j = min_position[1] ; j < max_position[1] ; j++){
                    for(std::size_t i = min_position[0] ; i < max_position[0] ; i++){
                        auto cell_bounding_box = GetCellBoundingBox(i,j,k);
                        if(IsIntersected(i_geometrical_object->GetGeometry(), cell_bounding_box, Tolerance)){
                            GetCell(i,j,k).push_back(&(*i_geometrical_object));
                        }
                    }
                }
            }
        }        
    }

    /// Giving the min and max position of cells intersecting with the bounding box of the geometry.
    template<typename TGeometryType>
    void CalculateMinMaxPositions(TGeometryType const& TheGeometry, array_1d< std::size_t, 3 >& MinPosition, array_1d< std::size_t, 3 >& MaxPosition){
        if(TheGeometry.empty())
            return;

        BoundingBox<Point> bounding_box(TheGeometry.begin(), TheGeometry.end());

        for(int i = 0; i < 3; i++ ) {
            MinPosition[ i ] = CalculatePosition( bounding_box.GetMinPoint()[i], i );
            MaxPosition[ i ] = CalculatePosition( bounding_box.GetMaxPoint()[i], i ) + 1;
        }
    }

    /// calculating the cell position of a given coordinate in respective axis given by dimension
    std::size_t CalculatePosition( double Coordinate, int ThisDimension ) const {
        auto distance = Coordinate - mBoundingBox.GetMinPoint()[ ThisDimension ];
        distance = ( distance < 0.00 ) ? 0.00 : distance;
        std::size_t position =
            static_cast< std::size_t >( distance * mInverseOfCellSize[ ThisDimension ] );
        std::size_t result= ( position > mNumberOfCells[ ThisDimension ] - 1 )
                                ? mNumberOfCells[ ThisDimension ] - 1
                                : position;
        return result;
    }
    
    /// Expands by the tolerance the geometry bounding box and checks the intersection
    template<typename TGeometryType>
    static inline bool IsIntersected(TGeometryType& TheGeometry, BoundingBox<Point> Box, const double tolerance)
    {
        Point rLowPointTolerance;
        Point rHighPointTolerance;
        
        for(std::size_t i = 0; i<3; i++)
        {
            rLowPointTolerance[i]  =  Box.GetMinPoint()[i] - tolerance;
            rHighPointTolerance[i] =  Box.GetMaxPoint()[i] + tolerance;
        }
        
        return  TheGeometry.HasIntersection(rLowPointTolerance,rHighPointTolerance);
    }

    /// Searchs in objects in the given cell for the ones with distance less than given radius + tolerance.
    template<typename TPointType>
    void SearchInRadiusInCell(CellType const& TheCell, TPointType const& ThePoint, double Radius, std::unordered_set<GeometricalObject*>& rResults) {
        for(auto p_geometrical_object : TheCell){  
            auto& geometry = p_geometrical_object->GetGeometry();
            // TODO: Change this to new Distance method of the geometry to be more general
            double distance = GeometryUtils::PointDistanceToTriangle3D(
            geometry[0],
            geometry[1],
            geometry[2],
            ThePoint);
            if((Radius + Tolerance) > distance){
                rResults.insert(p_geometrical_object);
            }
        }
    }

   /// Searchs in objects in the given cell for the nearest one.
    template<typename TPointType>
    void SearchNearestInCell(CellType const& TheCell, TPointType const& ThePoint, ResultType& rResult, double MaxRadius) {
        for(auto p_geometrical_object : TheCell){  
            auto& geometry = p_geometrical_object->GetGeometry();
            // TODO: Change this to new Distance method of the geometry to be more general
            double distance = GeometryUtils::PointDistanceToTriangle3D(
            geometry[0],
            geometry[1],
            geometry[2],
            ThePoint);
            if ((distance < rResult.GetDistance()) && (distance < MaxRadius)) {
                rResult.Set(p_geometrical_object);
                rResult.SetDistance(distance);
            }
        }
    }

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


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                GeometricalObjectsBins& rThis){
                    return rIStream;
                }

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const GeometricalObjectsBins& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_GEOMETRICAL_OBJECTS_BINS_H_INCLUDED  defined
