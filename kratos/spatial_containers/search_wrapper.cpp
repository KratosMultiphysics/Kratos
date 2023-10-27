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

// Project includes
#include "utilities/search_utilities.h"
#include "spatial_containers/search_wrapper.h"

namespace Kratos
{

template<class TSearchObject>
BoundingBox<Point> SearchWrapper<TSearchObject>::GetBoundingBox() const
{
    // Generate BB
    BoundingBox<Point> bb;
    auto& r_max = bb.GetMaxPoint();
    auto& r_min = bb.GetMinPoint();

    const auto& r_local_bb = mpSearchObject->GetBoundingBox();
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

template<class TSearchObject>
int SearchWrapper<TSearchObject>::GetRank() const
{
    return mrDataCommunicator.Rank();
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject>
int SearchWrapper<TSearchObject>::GetWorldSize() const
{
    return mrDataCommunicator.Size();
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject>
void SearchWrapper<TSearchObject>::InitializeGlobalBoundingBoxes()
{
    // Just executed in MPI
    if (mrDataCommunicator.IsDistributed()) {
        // We get the world size
        const int world_size = GetWorldSize();

        // Set up the global bounding boxes
        if (static_cast<int>(mGlobalBoundingBoxes.size()) != 6*world_size) {
            mGlobalBoundingBoxes.resize(6*world_size);
        }

        // Set up the local bounding boxes
        std::vector<double> local_bounding_box(6);
        const auto& r_bb = mpSearchObject->GetBoundingBox();
        const auto& r_max = r_bb.GetMaxPoint();
        const auto& r_min = r_bb.GetMinPoint();
        for (int i = 0; i < 3; ++i) {
            local_bounding_box[2 * i] = r_max[i];
            local_bounding_box[2 * i + 1] = r_min[i];
        }

        // Gather all bounding boxes
        mrDataCommunicator.AllGather(local_bounding_box, mGlobalBoundingBoxes);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject>
std::vector<int> SearchWrapper<TSearchObject>::RansksPointIsInsideBoundingBox(const array_1d<double, 3>& rCoords)
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

template<class TSearchObject>
std::vector<int> SearchWrapper<TSearchObject>::RansksPointIsInsideBoundingBoxWithTolerance(
    const array_1d<double, 3>& rCoords,
    const double Tolerance
    )
{
    std::vector<int> ranks;
    const int world_size = GetWorldSize();
    std::array<double, 6> local_bb;
    std::vector<double> bb_tolerance(mGlobalBoundingBoxes.size());
    SearchUtilities::ComputeBoundingBoxesWithToleranceCheckingNullBB(mGlobalBoundingBoxes, Tolerance, bb_tolerance);
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

/***********************************************************************************/
/***********************************************************************************/

template<class TSearchObject>
const Parameters SearchWrapper<TSearchObject>::GetDefaultParameters() const 
{
    const Parameters default_parameters = Parameters(R"(
    {
        "allocation_size" : 1000,
        "bucket_size"     : 4
    })" );
    return default_parameters;
}

/***********************************************************************************/
/***********************************************************************************/

template class SearchWrapper<GeometricalObjectsBins>;

}  // namespace Kratos.