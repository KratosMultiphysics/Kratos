//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

// System includes
#include <numeric>

// External includes
#include "mpi.h"

// Project includes
#include "containers/model.h"
#include "input_output/vtk_output.h"
#include "interface_communicator_mpi.h"
#include "custom_utilities/mapper_utilities.h"

namespace Kratos
{

using SizeType = std::size_t;

/***********************************************************************************/
/* PUBLIC Methods */
/***********************************************************************************/
InterfaceCommunicatorMPI::InterfaceCommunicatorMPI(ModelPart& rModelPartOrigin,
                                MapperLocalSystemPointerVector& rMapperLocalSystems,
                                Parameters SearchSettings) :
    InterfaceCommunicator(rModelPartOrigin,
                          rMapperLocalSystems,
                          SearchSettings)
{
    // Initialize MPI
    mSearchData.Initialize();
}

/***********************************************************************************/
/* PROTECTED Methods */
/***********************************************************************************/
void InterfaceCommunicatorMPI::InitializeSearch(const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo)
{
    InterfaceCommunicator::InitializeSearch(rpRefInterfaceInfo);

    // Exchange Bounding Boxes => has to be done every time
    ComputeGlobalBoundingBoxes();
}

void InterfaceCommunicatorMPI::InitializeSearchIteration(const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo)
{
    // Reset to zero
    std::fill(mSearchData.SendSizes.begin(), mSearchData.SendSizes.end(), 0);
    std::fill(mSearchData.RecvSizes.begin(), mSearchData.RecvSizes.end(), 0);

    // Apply tolerance to bounding boxes
    std::vector<double> bounding_boxes_with_tol;
    MPISearchUtilities::ComputeBoundingBoxesWithTolerance(mGlobalBoundingBoxes,
                                                          mSearchRadius,
                                                          bounding_boxes_with_tol);

    // Compute Candidate Partitions and fill the send buffer
    MapperUtilities::FillBufferBeforeLocalSearch(mrMapperLocalSystems,
                                                 bounding_boxes_with_tol,
                                                 GetBufferSizeEstimate(),
                                                 mSearchData.SendBufferDouble,
                                                 mSearchData.SendSizes);

    // copy the local information directly
    mSearchData.RecvBufferDouble[mSearchData.CommRank] = mSearchData.SendBufferDouble[mSearchData.CommRank];

    const int err = MPISearchUtilities::ExchangeDataAsync(mSearchData.SendBufferDouble, mSearchData.RecvBufferDouble, mSearchData);

    KRATOS_ERROR_IF_NOT(err == MPI_SUCCESS) << "Error in exchanging the information for "
        << "the construction of the MapperInterfaceInfos in MPI" << std::endl;

    // Construct MapperInterfaceInfos
    MapperUtilities::CreateMapperInterfaceInfosFromBuffer(mSearchData.RecvBufferDouble,
                                                          rpRefInterfaceInfo,
                                                          mSearchData.CommRank,
                                                          mMapperInterfaceInfosContainer);

    MPI_Barrier(MPI_COMM_WORLD);
}

void InterfaceCommunicatorMPI::FinalizeSearchIteration(const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo)
{
    // Reset to zero
    std::fill(mSearchData.SendSizes.begin(), mSearchData.SendSizes.end(), 0);
    std::fill(mSearchData.RecvSizes.begin(), mSearchData.RecvSizes.end(), 0);

    FilterInterfaceInfosSuccessfulSearch();

    MapperUtilities::FillBufferAfterLocalSearch(mMapperInterfaceInfosContainer,
                                                rpRefInterfaceInfo,
                                                mSearchData.CommRank,
                                                mSearchData.SendBufferChar,
                                                mSearchData.SendSizes);

    const int err = MPISearchUtilities::ExchangeDataAsync(mSearchData.SendBufferChar, mSearchData.RecvBufferChar, mSearchData);

    KRATOS_ERROR_IF_NOT(err == MPI_SUCCESS) << "Error in exchanging the "
        << "serialized MapperInterfaceInfos in MPI" << std::endl;

    MapperUtilities::DeserializeMapperInterfaceInfosFromBuffer(mSearchData.RecvBufferChar,
                                                               rpRefInterfaceInfo,
                                                               mSearchData.CommRank,
                                                               mMapperInterfaceInfosContainer);

    AssignInterfaceInfos();

    MPI_Barrier(MPI_COMM_WORLD);
}

/***********************************************************************************/
/* PRIVATE Methods */
/***********************************************************************************/
void InterfaceCommunicatorMPI::ComputeGlobalBoundingBoxes()
{
    const auto local_bounding_box = MapperUtilities::ComputeLocalBoundingBox(mrModelPartOrigin);

    if (static_cast<int>(mGlobalBoundingBoxes.size()) != 6*mSearchData.CommSize) {
        mGlobalBoundingBoxes.resize(6*mSearchData.CommSize);
    }

    MPI_Allgather(local_bounding_box.data(),   6, MPI_DOUBLE,
                  mGlobalBoundingBoxes.data(), 6, MPI_DOUBLE,
                  MPI_COMM_WORLD);

    const bool print_bounding_boxes_to_file = mSearchSettings.Has("print_bounding_boxes_to_file") ? mSearchSettings["print_bounding_boxes_to_file"].GetBool() : false;

    if (print_bounding_boxes_to_file &&
        mrModelPartOrigin.GetCommunicator().GetDataCommunicator().IsDefinedOnThisRank() &&
        mrModelPartOrigin.GetCommunicator().MyPID() == 0) {
        Model tmp_model;
        ModelPart& r_bbox_model_part = tmp_model.CreateModelPart("bounding_box");
        auto props = r_bbox_model_part.CreateNewProperties(0);

        // create hexahedral elements representing the bounding boxes
        for (int i=0; i<mSearchData.CommSize; ++i) {
            const double x_max = mGlobalBoundingBoxes[(i*6)];
            const double x_min = mGlobalBoundingBoxes[(i*6)+1];
            const double y_max = mGlobalBoundingBoxes[(i*6)+2];
            const double y_min = mGlobalBoundingBoxes[(i*6)+3];
            const double z_max = mGlobalBoundingBoxes[(i*6)+4];
            const double z_min = mGlobalBoundingBoxes[(i*6)+5];

            if (x_max < x_min) {
                // the bounding boxes are initialized inverted
                // hence if this condition is true then it means that
                // this partition does not have part of the interface
                continue;
            }

            // create vertices
            r_bbox_model_part.CreateNewNode((i*8),   x_min, y_min, z_min);
            r_bbox_model_part.CreateNewNode((i*8)+1, x_max, y_min, z_min);
            r_bbox_model_part.CreateNewNode((i*8)+2, x_max, y_max, z_min);
            r_bbox_model_part.CreateNewNode((i*8)+3, x_min, y_max, z_min);
            r_bbox_model_part.CreateNewNode((i*8)+4, x_min, y_min, z_max);
            r_bbox_model_part.CreateNewNode((i*8)+5, x_max, y_min, z_max);
            r_bbox_model_part.CreateNewNode((i*8)+6, x_max, y_max, z_max);
            r_bbox_model_part.CreateNewNode((i*8)+7, x_min, y_max, z_max);

            // create hexa
            std::vector<ModelPart::IndexType> bbox_vertex_ids(8);
            std::iota(bbox_vertex_ids.begin(), bbox_vertex_ids.end(), i*8);
            auto p_elem = r_bbox_model_part.CreateNewElement("Element3D8N", i, bbox_vertex_ids, props);
            p_elem->SetValue(PARTITION_INDEX, i);
        }

        const std::string file_name = "MapperMPISearch_BoundingBoxes_" + mrModelPartOrigin.FullName();

        KRATOS_INFO("MPIMapper") << "Printing file with search bounding boxes: " << file_name << ".vtk" << std::endl;

        Parameters vtk_params( R"({
            "file_format"                        : "binary",
            "save_output_files_in_folder"        : true,
            "output_path"                        : "",
            "element_data_value_variables"       : ["PARTITION_INDEX"]
        })");

        if (mSearchSettings.Has("bounding_boxes_file_path")) {
            vtk_params["output_path"].SetString(mSearchSettings["bounding_boxes_file_path"].GetString());
        }

        VtkOutput(r_bbox_model_part, vtk_params).PrintOutput(file_name);
    }
}

}  // namespace Kratos.
