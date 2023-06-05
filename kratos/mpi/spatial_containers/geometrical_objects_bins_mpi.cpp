//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes
#include "mpi.h"

// Project includes
#include "input_output/vtk_output.h"
#include "mpi/includes/mpi_data_communicator.h"
#include "utilities/search_utilities.h"
#include "mpi/spatial_containers/geometrical_objects_bins_mpi.h"

namespace Kratos
{

BoundingBox<Point> GeometricalObjectsBinsMPI::GetBoundingBox() const
{
    // Generate BB
    BoundingBox<Point> bb;
    auto& r_max = bb.GetMaxPoint();
    auto& r_min = bb.GetMinPoint();

    const auto& r_local_bb = BaseType::GetBoundingBox();
    const auto& r_local_max = r_local_bb.GetMaxPoint();
    const auto& r_local_min = r_local_bb.GetMinPoint();

    // Getting max values
    r_max[0] = mrDataCommunicator.MaxAll(r_local_max[0]);
    r_max[1] = mrDataCommunicator.MaxAll(r_local_max[1]);
    r_max[2] = mrDataCommunicator.MaxAll(r_local_max[2]);

    // Getting min values
    r_min[0] = mrDataCommunicator.MinAll(r_local_min[0]);
    r_min[1] = mrDataCommunicator.MinAll(r_local_min[1]);
    r_min[2] = mrDataCommunicator.MinAll(r_local_min[2]);

    return bb;
}

/***********************************************************************************/
/***********************************************************************************/

void GeometricalObjectsBinsMPI::SearchInRadius(
    const Point& rPoint,
    const double Radius,
    ResultTypeContainer& rResults,
    const bool SyncronizeResults
    )
{
    // Check if the point is inside the set
    if (BaseType::PointIsInsideBoundingBoxWithTolerance(rPoint, Radius)) {
        // Call local search
        BaseType::SearchInRadius(rPoint, Radius, rResults, false);
    }

    // Synchronize if needed
    if (SyncronizeResults) {
        rResults.SynchronizeAll(mrDataCommunicator);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void GeometricalObjectsBinsMPI::SearchNearestInRadius(
    const Point& rPoint,
    const double Radius,
    ResultTypeContainer& rResults,
    const bool SyncronizeResults
    )
{
    // Result to return
    ResultType local_result;

    // Get the rank
    const int current_rank = GetRank();

    // Check if the point is inside the set
    if (BaseType::PointIsInsideBoundingBoxWithTolerance(rPoint, Radius)) {
        // Call local search
        local_result = BaseType::SearchNearestInRadius(rPoint, Radius);
    }

    /* Now sync results between partitions */

    // Get the distance
    const double local_distance = (local_result.IsObjectFound() && local_result.IsDistanceCalculated()) ? local_result.GetDistance() : std::numeric_limits<double>::max();

    // Find the minimum value and the rank that holds it
    struct {
        double value;
        int rank;
    } local_min, global_min;

    local_min.value = local_distance;
    local_min.rank = current_rank;
    MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPIDataCommunicator::GetMPICommunicator(mrDataCommunicator));

    // Get the solution from the computed_rank
    if (global_min.rank == current_rank) {
        // Add the local search
        rResults.AddResult(local_result);
    }

    // Synchronize if needed
    if (SyncronizeResults) {
        rResults.SynchronizeAll(mrDataCommunicator);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void GeometricalObjectsBinsMPI::SearchNearest(
    const Point& rPoint,
    ResultTypeContainer& rResults,
    const bool SyncronizeResults
    )
{
    ResultType current_result;
    const auto bb = GetBoundingBox();
    const array_1d<double, 3> box_size = bb.GetMaxPoint() - bb.GetMinPoint();
    const double max_radius= *std::max_element(box_size.begin(), box_size.end());
    this->SearchNearestInRadius(rPoint, max_radius, rResults, SyncronizeResults);
}

/***********************************************************************************/
/***********************************************************************************/

void GeometricalObjectsBinsMPI::SearchIsInside(
    const Point& rPoint,
    ResultTypeContainer& rResults,
    const bool SyncronizeResults
    )
{
    // Result to return
    ResultType local_result;

    // Get the rank
    const int current_rank = GetRank();

    // Check if the point is inside the set
    int computed_rank = std::numeric_limits<int>::max();
    if (BaseType::PointIsInsideBoundingBox(rPoint)) {
        // Call local search
        local_result = BaseType::SearchIsInside(rPoint);

        // Set current rank
        computed_rank = current_rank;
    }

    // Now sync results between partitions
    //KRATOS_WARNING("GeometricalObjectsBinsMPI.SearchIsInside") << "This assumes that first one of the viable results is the closest one. The algorithm  sets distance to 0, so it is ambiguous which is the closest one." << std::endl;

    // Min rank of all ranks
    computed_rank = mrDataCommunicator.MinAll(computed_rank);

    // Get the solution from the computed_rank
    if (computed_rank == current_rank) {
        // Add the local search
        rResults.AddResult(local_result);
    }

    // Synchronize if needed
    if (SyncronizeResults) {
        rResults.SynchronizeAll(mrDataCommunicator);
    }
}

/***********************************************************************************/
/***********************************************************************************/

int GeometricalObjectsBinsMPI::GetRank() const
{
    return mrDataCommunicator.Rank();
}

/***********************************************************************************/
/***********************************************************************************/

int GeometricalObjectsBinsMPI::GetWorldSize() const
{
    return mrDataCommunicator.Size();
}

/***********************************************************************************/
/***********************************************************************************/

void GeometricalObjectsBinsMPI::InitializeGlobalBoundingBoxes()
{
    // We get the world size
    const int world_size = GetWorldSize();

    // Set up the global bounding boxes
    if (static_cast<int>(mGlobalBoundingBoxes.size()) != 6*world_size) {
        mGlobalBoundingBoxes.resize(6*world_size);
    }

    // Set up the local bounding boxes
    std::vector<double> local_bounding_box(6);
    const auto& r_bb = BaseType::GetBoundingBox();
    const auto& r_max = r_bb.GetMaxPoint();
    const auto& r_min = r_bb.GetMinPoint();
    for (int i = 0; i < 3; ++i) {
        local_bounding_box[2 * i] = r_max[i];
        local_bounding_box[2 * i + 1] = r_min[i];
    }

    MPI_Allgather(local_bounding_box.data(),   6, MPI_DOUBLE,
                  mGlobalBoundingBoxes.data(), 6, MPI_DOUBLE,
                  MPI_COMM_WORLD);
}

/***********************************************************************************/
/***********************************************************************************/

std::vector<int> GeometricalObjectsBinsMPI::RansksPointIsInsideBoundingBox(const array_1d<double, 3>& rCoords)
{
    std::vector<int> ranks;
    const int world_size = GetWorldSize();
    std::array<double, 6> local_bb;
    const auto it_begin = mGlobalBoundingBoxes.begin();
    for (int i = 0; i < world_size; ++i) {
        auto vec_it = it_begin + 6 * i;
        for (unsigned int j = 0; j < 6; ++j, ++vec_it) {
            local_bb[j] = *vec_it;
        }
        if (SearchUtilities::PointIsInsideBoundingBox(local_bb, rCoords)) {
            ranks.push_back(i);
        }
    }

    return ranks;
}

/***********************************************************************************/
/***********************************************************************************/

std::vector<int> GeometricalObjectsBinsMPI::RansksPointIsInsideBoundingBoxWithTolerance(
    const array_1d<double, 3>& rCoords,
    const double Tolerance
    )
{
    std::vector<int> ranks;
    const int world_size = GetWorldSize();
    std::array<double, 6> local_bb;
    std::vector<double> bb_tolerance(mGlobalBoundingBoxes.size());
    SearchUtilities::ComputeBoundingBoxesWithTolerance(mGlobalBoundingBoxes, Tolerance, bb_tolerance);
    const auto it_begin = bb_tolerance.begin();
    for (int i = 0; i < world_size; ++i) {
        auto vec_it = it_begin + 6 * i;
        for (unsigned int j = 0; j < 6; ++j, ++vec_it) {
            local_bb[j] = *vec_it;
        }
        if (SearchUtilities::PointIsInsideBoundingBox(local_bb, rCoords)) {
            ranks.push_back(i);
        }
    }

    return ranks;
}

}  // namespace Kratos.