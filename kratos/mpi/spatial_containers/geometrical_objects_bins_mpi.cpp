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
#include "includes/geometrical_object.h"
#include "input_output/vtk_output.h"
#include "utilities/pointer_communicator.h"
#include "mpi/utilities/mpi_search_utilities.h"
#include "mpi/spatial_containers/geometrical_objects_bins_mpi.h"

namespace Kratos
{

BoundingBox<Point> GeometricalObjectsBinsMPI::GetBoundingBox() const
{
    // Generate BB
    BoundingBox<Point> bb;
    auto& r_max = bb.GetMaxPoint();
    auto& r_min = bb.GetMinPoint();

    const auto& r_local_bb = mLocalGeometricalObjectsBins.GetBoundingBox();
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
    std::vector<ResultType>& rResults
    )
{
    // TODO
    const int current_rank = GetRank();
    std::vector<int> ranks = RansksPointIsInsideBoundingBoxWithTolerance(rPoint.Coordinates(), Radius);
}

/***********************************************************************************/
/***********************************************************************************/

GeometricalObjectsBinsMPI::ResultType GeometricalObjectsBinsMPI::SearchNearestInRadius(
    const Point& rPoint,
    const double Radius
    )
{
    ResultType current_result;
    current_result.SetDistance(std::numeric_limits<double>::max());
    // TODO
    const int current_rank = GetRank();
    std::vector<int> ranks = RansksPointIsInsideBoundingBoxWithTolerance(rPoint.Coordinates(), Radius);
    return current_result;
}

/***********************************************************************************/
/***********************************************************************************/

GeometricalObjectsBinsMPI::ResultType GeometricalObjectsBinsMPI::SearchNearest(const Point& rPoint)
{
    ResultType current_result;
    const auto bb = GetBoundingBox();
    const array_1d<double, 3> box_size = bb.GetMaxPoint() - bb.GetMinPoint();
    const double max_radius= *std::max_element(box_size.begin(), box_size.end());
    return SearchNearestInRadius(rPoint, max_radius);
}

/***********************************************************************************/
/***********************************************************************************/

GeometricalObjectsBinsMPI::ResultType GeometricalObjectsBinsMPI::SearchIsInside(const Point& rPoint)
{
    ResultType current_result;
    // TODO
    const int current_rank = GetRank();
    std::vector<int> ranks = RansksPointIsInsideBoundingBox(rPoint.Coordinates());
    return current_result;
}

/***********************************************************************************/
/***********************************************************************************/

// void GeometricalObjectsBinsMPI::InitializeSearch()
// {
//     KRATOS_TRY;

//     // Reset to zero
//     std::fill(mSendSizes.begin(), mSendSizes.end(), 0);
//     std::fill(mRecvSizes.begin(), mRecvSizes.end(), 0);

//     // Apply tolerance to bounding boxes
//     std::vector<double> bounding_boxes_with_tol;
//     MPISearchUtilities::ComputeBoundingBoxesWithTolerance(mGlobalBoundingBoxes,
//                                                           mRadius,
//                                                           bounding_boxes_with_tol);

//     // // Compute Candidate Partitions and fill the send buffer
//     // MapperUtilities::FillBufferBeforeLocalSearch(mrMapperLocalSystems,
//     //                                              bounding_boxes_with_tol,
//     //                                              GetBufferSizeEstimate(),
//     //                                              mSendBufferDouble,
//     //                                              mSendSizes);

//     // // Copy the local information directly
//     // mRecvBufferDouble[mCommRank] = mSendBufferDouble[mCommRank];

//     // const int err = MPISearchUtilities::ExchangeDataAsync(mSendBufferDouble, mRecvBufferDoubl, , mCommRank, mCommSize, mSendSizes, mRecvSizes);

//     // KRATOS_ERROR_IF_NOT(err == MPI_SUCCESS) << "Error in exchanging the information for "
//     //     << "the construction of the MapperInterfaceInfos in MPI" << std::endl;

//     // // Construct MapperInterfaceInfos
//     // MapperUtilities::CreateMapperInterfaceInfosFromBuffer(mRecvBufferDouble,
//     //                                                       rpRefInterfaceInfo,
//     //                                                       mCommRank,
//     //                                                       mMapperInterfaceInfosContainer);

//     MPI_Barrier(MPI_COMM_WORLD);

//     KRATOS_CATCH("");
// }

// /***********************************************************************************/
// /***********************************************************************************/

// void GeometricalObjectsBinsMPI::FinalizeSearch()
// {
//     KRATOS_TRY;

//     // Reset to zero
//     std::fill(mSendSizes.begin(), mSendSizes.end(), 0);
//     std::fill(mRecvSizes.begin(), mRecvSizes.end(), 0);

//     // FilterInterfaceInfosSuccessfulSearch();

//     // MapperUtilities::FillBufferAfterLocalSearch(mMapperInterfaceInfosContainer,
//     //                                             rpRefInterfaceInfo,
//     //                                             mCommRank,
//     //                                             mSendBufferChar,
//     //                                             mSendSizes);

//     // const int err = ExchangeDataAsync(mSendBufferChar, mRecvBufferChar);

//     // KRATOS_ERROR_IF_NOT(err == MPI_SUCCESS) << "Error in exchanging the "
//     //     << "serialized MapperInterfaceInfos in MPI" << std::endl;

//     // MapperUtilities::DeserializeMapperInterfaceInfosFromBuffer(mRecvBufferChar,
//     //                                                            rpRefInterfaceInfo,
//     //                                                            mCommRank,
//     //                                                            mMapperInterfaceInfosContainer);

//     // AssignInterfaceInfos();

//     MPI_Barrier(MPI_COMM_WORLD);

//     KRATOS_CATCH("");
// }

// void GeometricalObjectsBinsMPI::SynchronizeSearchInRadius(std::vector<ResultType>& rLocalResults)
// {
//     KRATOS_TRY;

//     // A priori some results can be duplicated, therefore we need to  first creeate the list of results to actually syncronize

//     MPI_Barrier(MPI_COMM_WORLD);

//     KRATOS_CATCH("");
// }

// /***********************************************************************************/
// /***********************************************************************************/

// GeometricalObjectsBinsMPI::ResultType GeometricalObjectsBinsMPI::SynchronizeSearchNearestInRadius(ResultType& rLocalResult)
// {
//     KRATOS_TRY;

//     // The same algorithm as SynchronizeSearchNearest, because we consider the closest one of all the local results
//     return SynchronizeSearchNearest(rLocalResult);

//     KRATOS_CATCH("");
// }

// /***********************************************************************************/
// /***********************************************************************************/

// GeometricalObjectsBinsMPI::ResultType GeometricalObjectsBinsMPI::SynchronizeSearchNearest(ResultType& rLocalResult)
// {
//     KRATOS_TRY;

//     // Get the distance
//     const double local_distance = (rLocalResult.IsObjectFound() && rLocalResult.IsDistanceCalculated()) ? rLocalResult.GetDistance() : std::numeric_limits<double>::max();

//     // Find the minimum value and the rank that holds it
//     struct {
//         double value;
//         int rank;
//     } local_min, global_min;

//     local_min.value = local_distance;
//     local_min.rank = mCommRank;
//     MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);

//     std::vector<GlobalPointer<GeometricalObject>> list(1, rLocalResult.Get());
//     auto global_pointer_communicator = GlobalPointerCommunicator<GeometricalObject>(mrDataCommunicator, list.begin(), list.end());

//     // TODO: Pensar en usar el serializer

//     // Get the global minimum
//     auto global_result = global_pointer_communicator.Apply([&](GlobalPointer<GeometricalObject>& r_object) {
//         if (local_min.rank == mCommRank) {
//             auto result = ResultType(r_object.get());
//             result.SetDistance(global_min.value);
//             return result;
//         }
//     });

//     return rLocalResult;

//     KRATOS_CATCH("");
// }

// /***********************************************************************************/
// /***********************************************************************************/

// GeometricalObjectsBinsMPI::ResultType GeometricalObjectsBinsMPI::SynchronizeSearchIsInside(ResultType& rLocalResult)
// {
//     KRATOS_TRY;

//     KRATOS_WARNING("GeometricalObjectsBinsMPI.SynchronizeSearchIsInside") << "This assumes that first one of the viable results is the closest one. The algorithm  sets distance to 0, so it is ambiguous which is the closest one." << std::endl;
//     return SynchronizeSearchNearest(rLocalResult);

//     KRATOS_CATCH("");
// }

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
    const auto& r_bb = mLocalGeometricalObjectsBins.GetBoundingBox();
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
        for (size_t i = 0; i < 6; ++i, ++vec_it) {
            local_bb[i] = *vec_it;
        }
        if (MPISearchUtilities::PointIsInsideBoundingBox(local_bb, rCoords)) {
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
    MPISearchUtilities::ComputeBoundingBoxesWithTolerance(mGlobalBoundingBoxes, Tolerance, bb_tolerance);
    const auto it_begin = bb_tolerance.begin();
    for (int i = 0; i < world_size; ++i) {
        auto vec_it = it_begin + 6 * i;
        for (size_t i = 0; i < 6; ++i, ++vec_it) {
            local_bb[i] = *vec_it;
        }
        if (MPISearchUtilities::PointIsInsideBoundingBox(local_bb, rCoords)) {
            ranks.push_back(i);
        }
    }

    return ranks;
}

}  // namespace Kratos.