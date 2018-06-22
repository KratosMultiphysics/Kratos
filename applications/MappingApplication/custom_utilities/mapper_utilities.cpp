//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

// System includes

// External includes

// Project includes
#include "mapper_utilities.h"
#include "mapping_application_variables.h"

namespace Kratos
{
namespace MapperUtilities
{

using IndexType = std::size_t;
using SizeType = std::size_t;

void AssignInterfaceEquationIds(Communicator& rModelPartCommunicator)
{
    const int num_nodes_local = rModelPartCommunicator.LocalMesh().NumberOfNodes();

    int num_nodes_accumulated;

    rModelPartCommunicator.ScanSum(num_nodes_local, num_nodes_accumulated);

    const int start_equation_id = num_nodes_accumulated - num_nodes_local;

    const auto nodes_begin = rModelPartCommunicator.LocalMesh().NodesBegin();

    #pragma omp parallel for
    for (int i=0; i<num_nodes_local; ++i)
    {
        // TODO this should be working in omp, not usre though
        ( nodes_begin + i )->SetValue(INTERFACE_EQUATION_ID, start_equation_id + i);
    }

    rModelPartCommunicator.SynchronizeVariable(INTERFACE_EQUATION_ID);
}

double ComputeSearchRadius(ModelPart& rModelPart, int EchoLevel)
{
    double search_safety_factor = 1.2;
    double max_element_size = 0.0;

    int num_conditions_global = ComputeNumberOfConditions(rModelPart);
    int num_elements_global = ComputeNumberOfElements(rModelPart);

    if (num_conditions_global > 0)
        max_element_size = ComputeMaxEdgeLengthLocal(rModelPart.GetCommunicator().LocalMesh().Conditions());
    else if (num_elements_global > 0)
        max_element_size = ComputeMaxEdgeLengthLocal(rModelPart.GetCommunicator().LocalMesh().Elements());
    else
    {
        if (EchoLevel >= 2 && rModelPart.GetCommunicator().MyPID() == 0)
            std::cout << "MAPPER WARNING, no conditions/elements for search radius "
                        << "computations in ModelPart \"" << rModelPart.Name() << "\" found, "
                        << "using nodes (less efficient, bcs search radius will be larger)" << std::endl;
        max_element_size = ComputeMaxEdgeLengthLocal(rModelPart.GetCommunicator().LocalMesh().Nodes());
    }

    rModelPart.GetCommunicator().MaxAll(max_element_size); // Compute the maximum among the partitions
    return max_element_size * search_safety_factor;
}

void CheckInterfaceModelParts(const int CommRank)
{
    // const int num_nodes_origin = MapperUtilities::ComputeNumberOfNodes(mrModelPartOrigin);
    // const int num_conditions_origin = MapperUtilities::ComputeNumberOfConditions(mrModelPartOrigin);
    // const int num_elements_origin = MapperUtilities::ComputeNumberOfElements(mrModelPartOrigin);

    // const int num_nodes_destination = MapperUtilities::ComputeNumberOfNodes(mrModelPartDestination);
    // const int num_conditions_destination = MapperUtilities::ComputeNumberOfConditions(mrModelPartDestination);
    // const int num_elements_destination = MapperUtilities::ComputeNumberOfElements(mrModelPartDestination);

    // // Check if the ModelPart contains entities
    // KRATOS_ERROR_IF(num_nodes_origin + num_conditions_origin + num_elements_origin < 1)
    //     << "Neither Nodes nor Conditions nor Elements found "
    //     << "in the Origin ModelPart" << std::endl;

    // KRATOS_ERROR_IF(num_nodes_destination + num_conditions_destination + num_elements_destination < 1)
    //     << "Neither Nodes nor Conditions nor Elements found "
    //     << "in the Destination ModelPart" << std::endl;

    // // Check if the inpt ModelParts contain both Elements and Conditions
    // // This is NOT possible, bcs the InterfaceObjects are constructed
    // // with whatever exists in the Modelpart (see the InterfaceObjectManagerBase,
    // // function "InitializeInterfaceGeometryObjectManager")
    // KRATOS_ERROR_IF(num_conditions_origin > 0 && num_elements_origin > 0)
    //     << "Origin ModelPart contains both Conditions and Elements "
    //     << "which is not permitted" << std::endl;

    // KRATOS_ERROR_IF(num_conditions_destination > 0 && num_elements_destination > 0)
    //     << "Destination ModelPart contains both Conditions and Elements "
    //     << "which is not permitted" << std::endl;

    // if (mEchoLevel >= 2) {
    //     std::vector<double> model_part_origin_bbox = MapperUtilities::ComputeModelPartBoundingBox(mrModelPartOrigin);
    //     std::vector<double> model_part_destination_bbox = MapperUtilities::ComputeModelPartBoundingBox(mrModelPartDestination);

    //     bool bbox_overlapping = MapperUtilities::ComputeBoundingBoxIntersection(
    //                                                 model_part_origin_bbox,
    //                                                 model_part_destination_bbox);
    //     if(CommRank == 0)
    //     {
    //         if (!bbox_overlapping) {
    //             std::cout << "MAPPER WARNING, the bounding boxes of the "
    //                         << "Modelparts do not overlap! "
    //                         << MapperUtilities::PrintModelPartBoundingBoxes(model_part_origin_bbox,
    //                                                                         model_part_destination_bbox)
    //                         << std::endl;
    //         } else if (mEchoLevel >= 3)
    //         {
    //             std::cout << MapperUtilities::PrintModelPartBoundingBoxes(model_part_origin_bbox,
    //                                                                         model_part_destination_bbox)
    //                         << std::endl;
    //         }
    //     }
    // }
}

std::vector<double> ComputeLocalBoundingBox(ModelPart& rModelPart)
{
    std::vector<double> local_bounding_box {-1e10, 1e10, -1e10, 1e10, -1e10, 1e10}; // initialize "inverted"
    // xmax, xmin,  ymax, ymin,  zmax, zmin

    // loop over all nodes (local and ghost(necessary if conditions have only ghost nodes) )
    for (auto &r_node : rModelPart.Nodes())
    {
        local_bounding_box[0] = std::max(r_node.X(), local_bounding_box[0]);
        local_bounding_box[1] = std::min(r_node.X(), local_bounding_box[1]);
        local_bounding_box[2] = std::max(r_node.Y(), local_bounding_box[2]);
        local_bounding_box[3] = std::min(r_node.Y(), local_bounding_box[3]);
        local_bounding_box[4] = std::max(r_node.Z(), local_bounding_box[4]);
        local_bounding_box[5] = std::min(r_node.Z(), local_bounding_box[5]);
    }
    return local_bounding_box;
}

void ComputeBoundingBoxesWithTolerance(const std::vector<double>& rBoundingBoxes,
                                       const double Tolerance,
                                       std::vector<double>& rBoundingBoxesWithTolerance)
{
    const SizeType size_vec = rBoundingBoxes.size();

    KRATOS_ERROR_IF_NOT(std::fmod(size_vec, 6) == 0)
        << "Bounding Boxes size has to be a multiple of 6!" << std::endl;

    if (rBoundingBoxesWithTolerance.size() != size_vec)
        rBoundingBoxesWithTolerance.resize(size_vec);

    // Apply Tolerances
    for (IndexType i=0; i<size_vec; i+=2)
        rBoundingBoxesWithTolerance[i] = rBoundingBoxes[i] + Tolerance;

    for (IndexType i=1; i<size_vec; i+=2)
        rBoundingBoxesWithTolerance[i] = rBoundingBoxes[i] - Tolerance;
}

std::string BoundingBoxStringStream(const std::vector<double>& rBoundingBox)
{
    // xmax, xmin,  ymax, ymin,  zmax, zmin
    std::stringstream buffer;
    buffer << "[" << rBoundingBox[1] << " "    // xmin
                  << rBoundingBox[3] << " "    // ymin
                  << rBoundingBox[5] << "]|["  // zmin
                  << rBoundingBox[0] << " "    // xmax
                  << rBoundingBox[2] << " "    // ymax
                  << rBoundingBox[4] << "]";   // zmax
    return buffer.str();
}

bool PointIsInsideBoundingBox(const std::vector<double>& rBoundingBox,
                              const Point& rPoint)
{   // The Bounding Box should have some tolerance already!
    if (rPoint.X() < rBoundingBox[0] && rPoint.X() > rBoundingBox[1])   // check x-direction
    {
        if (rPoint.Y() < rBoundingBox[2] && rPoint.Y() > rBoundingBox[3])   // check y-direction
        {
            if (rPoint.Z() < rBoundingBox[4] && rPoint.Z() > rBoundingBox[5])   // check z-direction
                return true;
        }
    }
    return false;
}

void FillBufferBeforeLocalSearch(BufferType& rSendBuffer,
                                 std::vector<int>& rSendSizes,
                                 const std::vector<double>& rBoundingBoxes)
{
    const SizeType comm_size = rSendBuffer.size();
    const SizeType num_local_sys = mpMapperLocalSystems->size();

    std::vector<double> bounding_box(6);

    // #pragma omp parallel for => "bounding_box" has to be threadprivate!
    for (IndexType i_rank=0; i_rank<comm_size; ++i_rank)
    {
        for (IndexType i=0; i<num_local_sys; ++i)
        {
            for (IndexType j=0; j<6; ++j)
                bounding_box[j] = rBoundingBoxes[(i_rank*6) + j]; // retrieve bounding box of partition

            const auto& rp_local_sys = (*mpMapperLocalSystems)[i];

            if (!rp_local_sys->HasInterfaceInfo())
            {
                const auto& r_coords = rp_local_sys->Coordinates();
                if (MapperUtilities::PointIsInsideBoundingBox(bounding_box, r_coords))
                {
                    // These push_backs are threadsafe bcs only one vector is accessed per thread!
                    rSendBuffer[i_rank].push_back(static_cast<double>(i)); // this it the "mSourceLocalSystemIndex" of the MapperInterfaceInfo
                    rSendBuffer[i_rank].push_back(r_coords[0]);
                    rSendBuffer[i_rank].push_back(r_coords[1]);
                    rSendBuffer[i_rank].push_back(r_coords[2]);

                    rSendSizes[i_rank] += 4;
                }
            }
        }
    }
}

void FillBufferAfterLocalSearch()
{
    /*
    Here we copy or move or however the InterfaceInfos that had a successful search to a new vector
    of the same type, which we later pass to the serialiazerhelper
    => on the receiving side we copy it to the the member-var mpInterfaceInfos


    */

}

void CreateMapperInterfaceInfosFromBuffer()
{
    mpMapperInterfaceInfosContainer->clear();
    mpMapperInterfaceInfosContainer->reserve(comm_size);

    array_1d<double, 3> coords;
    // Loop the ranks and construct the MapperInterfaceInfos
    for (int i_rank=0; i_rank<comm_size; ++i_rank)
    {
        KRATOS_DEBUG_ERROR_IF_NOT(std::fmod(recv_sizes[i_rank], 4) == 0) << "Rank " << comm_rank
                << " received wrong size from rank " << i_rank << std::endl;
        int num_objs = recv_sizes[i_rank] / 4;

        if ((*mpMapperInterfaceInfosContainer)[i_rank].size() != num_objs)
            (*mpMapperInterfaceInfosContainer)[i_rank].resize(num_objs);

        for (int j=0; j<num_objs; ++j)
        {
            // retrive data from buffer
            const int local_sys_idx = static_cast<int>(recv_buffer[i_rank][j*4]);
            coords[0] = recv_buffer[i_rank][j*4 + 1];
            coords[1] = recv_buffer[i_rank][j*4 + 2];
            coords[2] = recv_buffer[i_rank][j*4 + 3];
            (*mpMapperInterfaceInfosContainer)[i_rank][j] = rpRefInterfaceInfo->Create(coords, local_sys_idx, i_rank);
        }
    }
}

void MapperInterfaceInfoSerializer::save(Kratos::Serializer& rSerializer) const
{
    const SizeType num_ranks = mrInterfaceInfos.size();
    rSerializer.save("size1", num_ranks);

    for(IndexType i=0; i<num_ranks; ++i)
    {
        const SizeType size_rank = mrInterfaceInfos[i].size();
        rSerializer.save("size2", size_rank);

        for (IndexType j=0; j<size_rank; ++j)
            rSerializer.save("E", *(mrInterfaceInfos[i][j])); // NOT serializing the shared_ptr!
    }
}

void MapperInterfaceInfoSerializer::load(Kratos::Serializer& rSerializer)
{
    mrInterfaceInfos.clear(); // make sure it has no leftovers

    SizeType num_ranks;
    rSerializer.load("size1", num_ranks);
    if (mrInterfaceInfos.size() != num_ranks)
        mrInterfaceInfos.resize(num_ranks);

    for(IndexType i=0; i<num_ranks; ++i)
    {
        SizeType size_rank;
        rSerializer.load("size2", size_rank);

        if (mrInterfaceInfos[i].size() != size_rank)
            mrInterfaceInfos[i].resize(size_rank);

        for (IndexType j=0; j<size_rank; ++j)
        {
            // first we create a new object, then we load its data
            // this is needed bcs of the polymorphic behavior of the InterfaceInfos
            // i.e. in order to create the correct type
            // => the vector contains baseclass-pointers!
            // Jordi I am quite sure that I could get around it by registering it, what do you think... TODO
            // I think doing it manually is more efficient, which I want so I would probably leave it ...
            // The serializer does some nasty casting when pointers are serialized...
            // I could do a benchmark at some point but I highly doubt that the serializer is faster ...
            mrInterfaceInfos[i][j] = mrpRefInterfaceInfo->Create();
            rSerializer.load("E", *(mrInterfaceInfos[i][j])); // NOT serializing the shared_ptr!
        }
    }
}

} // namespace MapperUtilities
} // namespace Kratos.
