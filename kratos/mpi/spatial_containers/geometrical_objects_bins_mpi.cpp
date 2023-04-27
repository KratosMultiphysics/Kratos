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
#include "mpi/spatial_containers/geometrical_objects_bins_mpi.h"

namespace Kratos
{
void GeometricalObjectsBinsMPI::InitializeSearch()
{
    KRATOS_TRY;

    // Reset to zero
    std::fill(mSendSizes.begin(), mSendSizes.end(), 0);
    std::fill(mRecvSizes.begin(), mRecvSizes.end(), 0);

    // Apply tolerance to bounding boxes
    std::vector<double> bounding_boxes_with_tol;
    MPISearchUtilities::ComputeBoundingBoxesWithTolerance(mGlobalBoundingBoxes,
                                                          mRadius,
                                                          bounding_boxes_with_tol);

    // // Compute Candidate Partitions and fill the send buffer
    // MapperUtilities::FillBufferBeforeLocalSearch(mrMapperLocalSystems,
    //                                              bounding_boxes_with_tol,
    //                                              GetBufferSizeEstimate(),
    //                                              mSendBufferDouble,
    //                                              mSendSizes);

    // // Copy the local information directly
    // mRecvBufferDouble[mCommRank] = mSendBufferDouble[mCommRank];

    // const int err = MPISearchUtilities::ExchangeDataAsync(mSendBufferDouble, mRecvBufferDoubl, , mCommRank, mCommSize, mSendSizes, mRecvSizes);

    // KRATOS_ERROR_IF_NOT(err == MPI_SUCCESS) << "Error in exchanging the information for "
    //     << "the construction of the MapperInterfaceInfos in MPI" << std::endl;

    // // Construct MapperInterfaceInfos
    // MapperUtilities::CreateMapperInterfaceInfosFromBuffer(mRecvBufferDouble,
    //                                                       rpRefInterfaceInfo,
    //                                                       mCommRank,
    //                                                       mMapperInterfaceInfosContainer);

    MPI_Barrier(MPI_COMM_WORLD);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void GeometricalObjectsBinsMPI::FinalizeSearch()
{
    KRATOS_TRY;

    // Reset to zero
    std::fill(mSendSizes.begin(), mSendSizes.end(), 0);
    std::fill(mRecvSizes.begin(), mRecvSizes.end(), 0);

    // FilterInterfaceInfosSuccessfulSearch();

    // MapperUtilities::FillBufferAfterLocalSearch(mMapperInterfaceInfosContainer,
    //                                             rpRefInterfaceInfo,
    //                                             mCommRank,
    //                                             mSendBufferChar,
    //                                             mSendSizes);

    // const int err = ExchangeDataAsync(mSendBufferChar, mRecvBufferChar);

    // KRATOS_ERROR_IF_NOT(err == MPI_SUCCESS) << "Error in exchanging the "
    //     << "serialized MapperInterfaceInfos in MPI" << std::endl;

    // MapperUtilities::DeserializeMapperInterfaceInfosFromBuffer(mRecvBufferChar,
    //                                                            rpRefInterfaceInfo,
    //                                                            mCommRank,
    //                                                            mMapperInterfaceInfosContainer);

    // AssignInterfaceInfos();

    MPI_Barrier(MPI_COMM_WORLD);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void GeometricalObjectsBinsMPI::ComputeGlobalBoundingBoxes()
{
    // xmax, xmin, ymax, ymin, zmax, zmin
    std::array<double, 6> local_bounding_box;
    auto& r_min_point = mBoundingBox.GetMinPoint();
    auto& r_max_point = mBoundingBox.GetMaxPoint();
    for (unsigned int i=0; i<3; ++i) {
        local_bounding_box[i*2]   = r_max_point[i];
        local_bounding_box[i*2+1] = r_min_point[i];
    }

    if (static_cast<int>(mGlobalBoundingBoxes.size()) != 6*mCommSize) {
        mGlobalBoundingBoxes.resize(6*mCommSize);
    }

    MPI_Allgather(local_bounding_box.data(),   6, MPI_DOUBLE,
                  mGlobalBoundingBoxes.data(), 6, MPI_DOUBLE,
                  MPI_COMM_WORLD);
}

}  // namespace Kratos.