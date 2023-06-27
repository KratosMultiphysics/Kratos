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


// System includes

// External includes

// Project includes
#include "includes/geometrical_object.h"
#include "spatial_containers/geometrical_objects_bins.h"

namespace Kratos
{

GeometricalObjectsBins::CellType& GeometricalObjectsBins::GetCell(
    const std::size_t I,
    const std::size_t J,
    const std::size_t K
    )
{
    KRATOS_DEBUG_ERROR_IF(I > mNumberOfCells[0]) << "Index " << I << " is larger than number of cells in x direction : " << mNumberOfCells[0] << std::endl;
    KRATOS_DEBUG_ERROR_IF(J > mNumberOfCells[1]) << "Index " << J << " is larger than number of cells in y direction : " << mNumberOfCells[1] << std::endl;
    KRATOS_DEBUG_ERROR_IF(K > mNumberOfCells[2]) << "Index " << K << " is larger than number of cells in z direction : " << mNumberOfCells[2] << std::endl;

    const std::size_t index = I + J * mNumberOfCells[0] + K * mNumberOfCells[1] * mNumberOfCells[0];
    return mCells[index];
}

/***********************************************************************************/
/***********************************************************************************/

BoundingBox<Point> GeometricalObjectsBins::GetCellBoundingBox(
    const std::size_t I,
    const std::size_t J,
    const std::size_t K
    )
{
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

/***********************************************************************************/
/***********************************************************************************/

void GeometricalObjectsBins::SearchInRadius(
    const Point& rPoint,
    const double Radius,
    std::vector<ResultType>& rResults
    )
{
    std::unordered_set<GeometricalObject*> results;

    array_1d<std::size_t, Dimension> min_position;
    array_1d<std::size_t, Dimension> max_position;

    for(unsigned int i = 0; i < Dimension; i++ ) {
        min_position[i] = CalculatePosition(rPoint[i] - Radius, i);
        max_position[i] = CalculatePosition(rPoint[i] + Radius, i) + 1;
    }
    for(std::size_t k = min_position[2] ; k < max_position[2] ; k++){
        for(std::size_t j = min_position[1] ; j < max_position[1] ; j++){
            for(std::size_t i = min_position[0] ; i < max_position[0] ; i++){
                auto& r_cell = GetCell(i,j,k);
                SearchInRadiusInCell(r_cell, rPoint, Radius, results);
            }
        }
    }

    rResults.clear();
    for(auto p_object : results){
        rResults.push_back(ResultType(p_object));
    }
}

/***********************************************************************************/
/***********************************************************************************/

GeometricalObjectsBins::ResultType GeometricalObjectsBins::SearchNearestInRadius(
    const Point& rPoint,
    const double Radius
    )
{
    ResultType current_result;
    current_result.SetDistance(std::numeric_limits<double>::max());

    const double radius_increment = *std::max_element(mCellSizes.begin(), mCellSizes.end());

    array_1d<std::size_t, Dimension> min_position;
    array_1d<std::size_t, Dimension> max_position;

    for(double current_radius = radius_increment ; current_radius < Radius + radius_increment ; current_radius += radius_increment){
        current_radius = (current_radius > Radius) ? Radius : current_radius;
        for(unsigned int i = 0; i < Dimension; i++ ) {
            min_position[i] = CalculatePosition(rPoint[i] - current_radius, i);
            max_position[i] = CalculatePosition(rPoint[i] + current_radius, i) + 1;
        }

        for(std::size_t k = min_position[2] ; k < max_position[2] ; k++){
            for(std::size_t j = min_position[1] ; j < max_position[1] ; j++){
                for(std::size_t i = min_position[0] ; i < max_position[0] ; i++){
                    auto& r_cell = GetCell(i,j,k);
                    SearchNearestInCell(r_cell, rPoint, current_result, Radius);
                }
            }
        }
        bool all_cells_are_covered = (min_position[0] == 0) && (min_position[1] == 0) && (min_position[2] == 0);
        all_cells_are_covered &= (max_position[0] == mNumberOfCells[0] - 1) && (max_position[1] == mNumberOfCells[1] - 1) && (max_position[2] == mNumberOfCells[2] - 1);

        const bool object_found_within_current_radius = current_result.IsObjectFound() && (current_result.GetDistance() < current_radius);
        if (all_cells_are_covered || object_found_within_current_radius) {
            break;
        }
    }
    return current_result;
}

/***********************************************************************************/
/***********************************************************************************/

GeometricalObjectsBins::ResultType GeometricalObjectsBins::SearchNearest(const Point& rPoint)
{
    ResultType current_result;

    const array_1d<double, 3> box_size = mBoundingBox.GetMaxPoint() - mBoundingBox.GetMinPoint();
    const double max_radius= *std::max_element(box_size.begin(), box_size.end());

    return SearchNearestInRadius(rPoint, max_radius);
}

/***********************************************************************************/
/***********************************************************************************/

GeometricalObjectsBins::ResultType GeometricalObjectsBins::SearchIsInside(const Point& rPoint)
{
    ResultType current_result;
    current_result.SetDistance(std::numeric_limits<double>::max());

    array_1d<std::size_t, Dimension> position;
    for(unsigned int i = 0; i < Dimension; i++ ) {
        position[i] = CalculatePosition(rPoint[i], i);
    }

    auto& r_cell = GetCell(position[0],position[1],position[2]);
    SearchIsInsideInCell(r_cell, rPoint, current_result);

    return current_result;
}

/***********************************************************************************/
/***********************************************************************************/

void GeometricalObjectsBins::CalculateCellSize(const std::size_t NumberOfCells)
{
    const std::size_t avarage_number_of_cells = static_cast<std::size_t>(std::pow(static_cast<double>(NumberOfCells), 1.00 / Dimension));
    std::array<double, Dimension> lengths;
    double avarage_length = 0.0;
    for (unsigned int i = 0; i < Dimension; i++) {
        lengths[i] = mBoundingBox.GetMaxPoint()[i] - mBoundingBox.GetMinPoint()[i];
        avarage_length += lengths[i];
    }
    avarage_length *= 0.33333333333333333333333333333333;

    if (avarage_length < std::numeric_limits<double>::epsilon()) {
        mNumberOfCells = ScalarVector(3, 1);
        return;
    }

    for (unsigned int i = 0; i < Dimension; i++) {
        mNumberOfCells[i] = static_cast<std::size_t>(lengths[i] / avarage_length * avarage_number_of_cells) + 1;
        if (mNumberOfCells[i] > 1)
            mCellSizes[i] = lengths[i] / mNumberOfCells[i];
        else
            mCellSizes[i] = avarage_length;

        mInverseOfCellSize[i] = 1.00 / mCellSizes[i];
    }

}

/***********************************************************************************/
/***********************************************************************************/

std::size_t GeometricalObjectsBins::CalculatePosition(
    const double Coordinate,
    const int ThisDimension
    ) const
{
    auto distance = Coordinate - mBoundingBox.GetMinPoint()[ ThisDimension ];
    distance = ( distance < 0.0 ) ? 0.0 : distance;
    const std::size_t position = static_cast< std::size_t >( distance * mInverseOfCellSize[ ThisDimension ] );
    const std::size_t result= ( position > mNumberOfCells[ ThisDimension ] - 1 )
                            ? mNumberOfCells[ ThisDimension ] - 1
                            : position;
    return result;
}

/***********************************************************************************/
/***********************************************************************************/

void GeometricalObjectsBins::SearchInRadiusInCell(
    const CellType& rCell,
    const Point& rPoint,
    const double Radius,
    std::unordered_set<GeometricalObject*>& rResults
    )
{
    double distance = 0.0;
    for(auto p_geometrical_object : rCell){
        auto& r_geometry = p_geometrical_object->GetGeometry();
        distance = r_geometry.CalculateDistance(rPoint, Tolerance);
        if((Radius + Tolerance) > distance){
            rResults.insert(p_geometrical_object);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void GeometricalObjectsBins::SearchNearestInCell(
    const CellType& rCell,
    const Point& rPoint,
    ResultType& rResult,
    const double MaxRadius
    )
{
    double distance = 0.0;
    for(auto p_geometrical_object : rCell){
        auto& r_geometry = p_geometrical_object->GetGeometry();
        distance = r_geometry.CalculateDistance(rPoint, Tolerance);
        if ((distance < rResult.GetDistance()) && (distance < MaxRadius)) {
            rResult.Set(p_geometrical_object);
            rResult.SetDistance(distance);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void GeometricalObjectsBins::SearchIsInsideInCell(
    const CellType& rCell,
    const Point& rPoint,
    ResultType& rResult
    )
{
    array_1d<double, 3> point_local_coordinates;
    for(auto p_geometrical_object : rCell){
        auto& r_geometry = p_geometrical_object->GetGeometry();
        const bool is_inside = r_geometry.IsInside(rPoint,point_local_coordinates,Tolerance);
        if (is_inside) {
            rResult.Set(p_geometrical_object);
            rResult.SetDistance(0.0);
            return;
        }
    }
}

}  // namespace Kratos.