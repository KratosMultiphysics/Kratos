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

BoundingBox<GeometricalObjectsBins::PointType> GeometricalObjectsBins::GetCellBoundingBox(
    const std::size_t I,
    const std::size_t J,
    const std::size_t K
    )
{
    KRATOS_DEBUG_ERROR_IF(I > mNumberOfCells[0]) << "Index " << I << " is larger than number of cells in x direction : " << mNumberOfCells[0] << std::endl;
    KRATOS_DEBUG_ERROR_IF(J > mNumberOfCells[1]) << "Index " << J << " is larger than number of cells in y direction : " << mNumberOfCells[1] << std::endl;
    KRATOS_DEBUG_ERROR_IF(K > mNumberOfCells[2]) << "Index " << K << " is larger than number of cells in z direction : " << mNumberOfCells[2] << std::endl;

    BoundingBox<PointType> result;

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
    const PointType& rPoint,
    const double Radius,
    std::vector<ResultType>& rResults
    )
{
    // Clear the results
    rResults.clear();

    // Initialize the current size
    std::size_t current_size = 0;

    // Initialize the candidates
    std::unordered_set<GeometricalObject*> candidates;
    std::unordered_set<GeometricalObject*> to_remove;

    // Initialize the position bounds
    array_1d<std::size_t, Dimension> min_position;
    array_1d<std::size_t, Dimension> max_position;

    // Calculate the position bounds
    for(unsigned int i = 0; i < Dimension; i++ ) {
        min_position[i] = CalculatePosition(rPoint[i] - Radius, i);
        max_position[i] = CalculatePosition(rPoint[i] + Radius, i) + 1;
    }

    // Loop over the cells and gather candidates
    for(std::size_t k = min_position[2]; k < max_position[2]; k++) {
        for(std::size_t j = min_position[1]; j < max_position[1]; j++) {
            for(std::size_t i = min_position[0]; i < max_position[0]; i++) {
                auto& r_cell = GetCell(i, j, k);
                if (IsCellBoundingBoxInsideRadius(i, j, k, rPoint, Radius)) {
                    current_size = rResults.size();
                    rResults.reserve(current_size + r_cell.size());
                    for(auto p_geometrical_object : r_cell) {
                        to_remove.insert(p_geometrical_object);
                        auto& r_geometry = p_geometrical_object->GetGeometry();
                        const double distance = r_geometry.CalculateDistance(rPoint, mTolerance);
                        rResults.push_back(ResultType(p_geometrical_object));
                        rResults.back().SetDistance(distance);
                    }
                } else {
                    for(auto p_geometrical_object : r_cell) {
                        candidates.insert(p_geometrical_object);
                    }
                }
            }
        }
    }

    // Clear the candidates
    for (const auto& r_element : to_remove) {
        candidates.erase(r_element);
    }

    // Loop over the candidates and filter by distance and fill the results
    current_size = rResults.size();
    rResults.reserve(current_size + candidates.size());
    for(auto& p_geometrical_object : candidates) {
        auto& r_geometry = p_geometrical_object->GetGeometry();
        const double distance = r_geometry.CalculateDistance(rPoint, mTolerance);
        if((Radius + mTolerance) > distance) {
            rResults.push_back(ResultType(p_geometrical_object));
            rResults.back().SetDistance(distance);
        }
    }

    // Shrink the results
    rResults.shrink_to_fit();
}

/***********************************************************************************/
/***********************************************************************************/

GeometricalObjectsBins::ResultType GeometricalObjectsBins::SearchNearestInRadius(
    const PointType& rPoint,
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

void GeometricalObjectsBins::SearchInBoundingBox(
    const PointType& rPoint,
    const array_1d<double, 3>& rMinPoint,
    const array_1d<double, 3>& rMaxPoint,
    std::vector<ResultType>& rResults
    )
{
    // Clear the results
    rResults.clear();

    // Define bounding box
    BoundingBox<PointType> bounding_box;
    noalias(bounding_box.GetMinPoint().Coordinates()) = rMinPoint;
    noalias(bounding_box.GetMaxPoint().Coordinates()) = rMaxPoint;

    // Initialize the current size
    std::size_t current_size = 0;

    // Initialize the candidates
    std::unordered_set<GeometricalObject*> candidates;

    // Initialize the position bounds
    array_1d<std::size_t, Dimension> min_position;
    array_1d<std::size_t, Dimension> max_position;

    // Calculate the position bounds
    for(unsigned int i = 0; i < Dimension; i++ ) {
        min_position[i] = CalculatePosition(rMinPoint[i], i);
        max_position[i] = CalculatePosition(rMaxPoint[i], i) + 1;
    }

    // Loop over the cells and gather candidates
    for(std::size_t k = min_position[2]; k < max_position[2]; k++) {
        for(std::size_t j = min_position[1]; j < max_position[1]; j++) {
            for(std::size_t i = min_position[0]; i < max_position[0]; i++) {
                auto& r_cell = GetCell(i, j, k);
                if (IsCellBoundingBoxInsideBoundingBox(i, j, k, rPoint, bounding_box)) {
                    current_size = rResults.size();
                    rResults.reserve(current_size + r_cell.size());
                    for(auto p_geometrical_object : r_cell) {
                        auto& r_geometry = p_geometrical_object->GetGeometry();
                        const double distance = r_geometry.CalculateDistance(rPoint, mTolerance);
                        rResults.push_back(ResultType(p_geometrical_object));
                        rResults.back().SetDistance(distance);
                    }
                } else {
                    for(auto p_geometrical_object : r_cell) {
                        candidates.insert(p_geometrical_object);
                    }
                }
            }
        }
    }

    // Loop over the candidates and filter by distance and fill the results
    const auto& r_min_point = bounding_box.GetMinPoint();
    const auto& r_max_point = bounding_box.GetMaxPoint();
    current_size = rResults.size();
    rResults.reserve(current_size + candidates.size());
    for(auto& p_geometrical_object : candidates) {
        const auto& r_geometry = p_geometrical_object->GetGeometry();
        if(r_geometry.HasIntersection(r_min_point, r_max_point)) { // NOTE: This operation is expensive, we should think a way to avoid it
            const double distance = r_geometry.CalculateDistance(rPoint, mTolerance);
            rResults.push_back(ResultType(p_geometrical_object));
            rResults.back().SetDistance(distance);
        }
    }

    // Shrink the results
    rResults.shrink_to_fit();
}

/***********************************************************************************/
/***********************************************************************************/

GeometricalObjectsBins::ResultType GeometricalObjectsBins::SearchNearest(const PointType& rPoint)
{
    ResultType current_result;

    const array_1d<double, 3> box_size = mBoundingBox.GetMaxPoint() - mBoundingBox.GetMinPoint();
    const double max_radius= *std::max_element(box_size.begin(), box_size.end());

    return SearchNearestInRadius(rPoint, max_radius);
}

/***********************************************************************************/
/***********************************************************************************/

GeometricalObjectsBins::ResultType GeometricalObjectsBins::SearchIsInside(const PointType& rPoint)
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

bool GeometricalObjectsBins::IsCellBoundingBoxInsideRadius(
    const std::size_t I,
    const std::size_t J,
    const std::size_t K,
    const Point& rPoint,
    const double Radius
    )
{
    // Check the indices
    KRATOS_DEBUG_ERROR_IF(I > mNumberOfCells[0]) << "Index " << I << " is larger than number of cells in x direction : " << mNumberOfCells[0] << std::endl;
    KRATOS_DEBUG_ERROR_IF(J > mNumberOfCells[1]) << "Index " << J << " is larger than number of cells in y direction : " << mNumberOfCells[1] << std::endl;
    KRATOS_DEBUG_ERROR_IF(K > mNumberOfCells[2]) << "Index " << K << " is larger than number of cells in z direction : " << mNumberOfCells[2] << std::endl;

    // Get the bounding box points min point
    const array_1d<double, 3>& r_min_point = mBoundingBox.GetMinPoint();

    // Calculate the minimum point in the bounding box
    const double min_point_x = r_min_point[0] + I * mCellSizes[0];
    const double min_point_y = r_min_point[1] + J * mCellSizes[1];
    const double min_point_z = r_min_point[2] + K * mCellSizes[2];

    // Calculate the maximum point in the bounding box
    const double max_point_x = r_min_point[0] + (I + 1) * mCellSizes[0];
    const double max_point_y = r_min_point[1] + (J + 1) * mCellSizes[1];
    const double max_point_z = r_min_point[2] + (K + 1) * mCellSizes[2];

    // Determine the farthest point in the bounding box from the center point
    const double x_farthest = (rPoint[0] < (min_point_x + max_point_x) / 2) ? max_point_x : min_point_x;
    const double y_farthest = (rPoint[1] < (min_point_y + max_point_y) / 2) ? max_point_y : min_point_y;
    const double z_farthest = (rPoint[2] < (min_point_z + max_point_z) / 2) ? max_point_z : min_point_z;

    // Calculate the squared distance to avoid using sqrt
    const double dx = x_farthest - rPoint[0];
    const double dy = y_farthest - rPoint[1];
    const double dz = z_farthest - rPoint[2];
    const double distance_squared = dx * dx + dy * dy + dz * dz;

    // Compare squared distance with squared radius
    return distance_squared <= Radius * Radius;
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometricalObjectsBins::IsCellBoundingBoxInsideBoundingBox(
    const std::size_t I,
    const std::size_t J,
    const std::size_t K,
    const Point& rPoint,
    const BoundingBox<PointType>& rBoundingBox
    )
{
    // Check the indices
    KRATOS_DEBUG_ERROR_IF(I > mNumberOfCells[0]) << "Index " << I << " is larger than number of cells in x direction : " << mNumberOfCells[0] << std::endl;
    KRATOS_DEBUG_ERROR_IF(J > mNumberOfCells[1]) << "Index " << J << " is larger than number of cells in y direction : " << mNumberOfCells[1] << std::endl;
    KRATOS_DEBUG_ERROR_IF(K > mNumberOfCells[2]) << "Index " << K << " is larger than number of cells in z direction : " << mNumberOfCells[2] << std::endl;

    // Get the bounding box points min point
    const array_1d<double, 3>& r_min_point = mBoundingBox.GetMinPoint();

    // Calculate the minimum point in the bounding box
    array_1d<double, 3> min_point;
    min_point[0] = r_min_point[0] + I * mCellSizes[0];
    min_point[1] = r_min_point[1] + J * mCellSizes[1];
    min_point[2] = r_min_point[2] + K * mCellSizes[2];

    // Calculate the maximum point in the bounding box
    array_1d<double, 3> max_point;
    max_point[0] = r_min_point[0] + (I + 1) * mCellSizes[0];
    max_point[1] = r_min_point[1] + (J + 1) * mCellSizes[1];
    max_point[2] = r_min_point[2] + (K + 1) * mCellSizes[2];

    // We do the checks of both points
    return rBoundingBox.IsInside(min_point) && rBoundingBox.IsInside(max_point);
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometricalObjectsBins::PointIsInsideBoundingBox(const array_1d<double, 3>& rCoords)
{
    // Get the bounding box points
    const auto& r_max_point = mBoundingBox.GetMaxPoint();
    const auto& r_min_point = mBoundingBox.GetMinPoint();

    // The Bounding Box check
    if (rCoords[0] < r_max_point[0] && rCoords[0] > r_min_point[0])           // check x-direction
        if (rCoords[1] < r_max_point[1] && rCoords[1] > r_min_point[1])       // check y-direction
            if (rCoords[2] < r_max_point[2] && rCoords[2] > r_min_point[2])   // check z-direction
                return true;
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometricalObjectsBins::PointIsInsideBoundingBoxWithTolerance(
    const array_1d<double, 3>& rCoords,
    const double Tolerance
    )
{
    // Get the bounding box points
    auto max_point = mBoundingBox.GetMaxPoint();
    auto min_point = mBoundingBox.GetMinPoint();

    // Apply Tolerances (only in non zero BB cases)
    const double epsilon = std::numeric_limits<double>::epsilon();
    if (norm_2(max_point) > epsilon && norm_2(min_point) > epsilon) {
        for (unsigned int i=0; i<3; ++i) {
            max_point[i] += Tolerance;
            min_point[i] -= Tolerance;
        }
    }

    // The Bounding Box check
    if (rCoords[0] < max_point[0] && rCoords[0] > min_point[0])           // check x-direction
        if (rCoords[1] < max_point[1] && rCoords[1] > min_point[1])       // check y-direction
            if (rCoords[2] < max_point[2] && rCoords[2] > min_point[2])   // check z-direction
                return true;
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

void GeometricalObjectsBins::CalculateCellSize(const std::size_t NumberOfCells)
{
    const std::size_t average_number_of_cells = static_cast<std::size_t>(std::pow(static_cast<double>(NumberOfCells), 1.00 / Dimension));
    std::array<double, Dimension> lengths;
    double average_length = 0.0;
    for (unsigned int i = 0; i < Dimension; i++) {
        lengths[i] = mBoundingBox.GetMaxPoint()[i] - mBoundingBox.GetMinPoint()[i];
        average_length += lengths[i];
    }
    average_length *= 0.33333333333333333333333333333333;

    if (average_length < std::numeric_limits<double>::epsilon()) {
        mNumberOfCells = ScalarVector(3, 1);
        for (unsigned int i = 0; i < Dimension; i++) {
            mNumberOfCells[i] = 0;
            mCellSizes[i] = 0.0;
            mInverseOfCellSize[i] = std::numeric_limits<double>::max();
        }
        return;
    }

    for (unsigned int i = 0; i < Dimension; i++) {
        mNumberOfCells[i] = static_cast<std::size_t>(lengths[i] / average_length * average_number_of_cells) + 1;
        if (mNumberOfCells[i] > 1) {
            mCellSizes[i] = lengths[i] / mNumberOfCells[i];
        } else {
            mCellSizes[i] = average_length;
        }

        mInverseOfCellSize[i] = 1.0 / mCellSizes[i];
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

void GeometricalObjectsBins::SearchNearestInCell(
    const CellType& rCell,
    const PointType& rPoint,
    ResultType& rResult,
    const double MaxRadius
    )
{
    double distance = 0.0;
    for(auto p_geometrical_object : rCell){
        auto& r_geometry = p_geometrical_object->GetGeometry();
        distance = r_geometry.CalculateDistance(rPoint, mTolerance);
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
    const PointType& rPoint,
    ResultType& rResult
    )
{
    array_1d<double, 3> point_local_coordinates;
    for(auto p_geometrical_object : rCell){
        auto& r_geometry = p_geometrical_object->GetGeometry();
        const bool is_inside = r_geometry.IsInside(rPoint,point_local_coordinates,mTolerance);
        if (is_inside) {
            rResult.Set(p_geometrical_object);
            rResult.SetDistance(0.0);
            return;
        }
    }
}

}  // namespace Kratos.