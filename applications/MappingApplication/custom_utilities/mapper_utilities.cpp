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
#include "includes/stream_serializer.h"
#include "utilities/parallel_utilities.h"
#include "mapper_utilities.h"
#include "mapping_application_variables.h"

namespace Kratos
{
namespace MapperUtilities
{

typedef std::size_t SizeType;
typedef std::size_t IndexType;

void AssignInterfaceEquationIds(Communicator& rModelPartCommunicator)
{
    const int num_nodes_local = rModelPartCommunicator.LocalMesh().NumberOfNodes();
    int num_nodes_accumulated = rModelPartCommunicator.GetDataCommunicator().ScanSum(num_nodes_local);
    const int start_equation_id = num_nodes_accumulated - num_nodes_local;
    const auto nodes_begin = rModelPartCommunicator.LocalMesh().NodesBegin();

    IndexPartition<unsigned int>(num_nodes_local).for_each(
        [nodes_begin, start_equation_id](unsigned int i){
            (nodes_begin + i)->SetValue(INTERFACE_EQUATION_ID, start_equation_id + i);
        }
    );

    rModelPartCommunicator.SynchronizeNonHistoricalVariable(INTERFACE_EQUATION_ID);
}

template <typename T>
double ComputeMaxEdgeLengthLocal(const T& rEntityContainer)
{
    double max_element_size = 0.0;
    // Loop through each edge of a geometrical entity ONCE
    for (const auto& r_entity : rEntityContainer) {
        for (std::size_t i = 0; i < (r_entity.GetGeometry().size() - 1); ++i) {
            for (std::size_t j = i + 1; j < r_entity.GetGeometry().size(); ++j) {
                double edge_length = ComputeDistance(r_entity.GetGeometry()[i].Coordinates(),
                                                        r_entity.GetGeometry()[j].Coordinates());
                max_element_size = std::max(max_element_size, edge_length);
            }
        }
    }
    return max_element_size;
}

double ComputeSearchRadius(const ModelPart& rModelPart, const int EchoLevel)
{
    static constexpr double search_safety_factor = 1.2;
    double max_element_size = 0.0;

    const auto r_comm = rModelPart.GetCommunicator();

    if (r_comm.GlobalNumberOfConditions() > 0) {
        max_element_size = ComputeMaxEdgeLengthLocal(rModelPart.GetCommunicator().LocalMesh().Conditions());
    }
    else if (r_comm.GlobalNumberOfElements() > 0) {
        max_element_size = ComputeMaxEdgeLengthLocal(rModelPart.GetCommunicator().LocalMesh().Elements());
    }
    else {
        KRATOS_WARNING_IF("Mapper", EchoLevel > 0)
            << "No conditions/elements for computation of search radius found in\n"
            << "ModelPart \"" << rModelPart.Name() << "\", using nodes for computing it\n"
            << "(less efficient, because search radius will be larger)\n"
            << "It is recommended to specify the search-radius manually\n"
            << "through \"search_radius\" in the mapper-settings (~2*element-size)" << std::endl;

        const auto bounding_box = ComputeGlobalBoundingBox(rModelPart);
        const double dx = bounding_box[0] - bounding_box[1];
        const double dy = bounding_box[2] - bounding_box[3];
        const double dz = bounding_box[4] - bounding_box[5];

        const double nominator = std::sqrt((dx*dx) + (dy*dy) + (dz*dz));
        const double denominator = std::sqrt(static_cast<double>(r_comm.GlobalNumberOfNodes()));

        max_element_size = nominator / denominator;
    }

    max_element_size = r_comm.GetDataCommunicator().MaxAll(max_element_size); // Compute the maximum among the partitions
    return max_element_size * search_safety_factor;
}

double ComputeSearchRadius(const ModelPart& rModelPart1, const ModelPart& rModelPart2, const int EchoLevel)
{
    double search_radius = std::max(ComputeSearchRadius(rModelPart1, EchoLevel),
                                    ComputeSearchRadius(rModelPart2, EchoLevel));

    KRATOS_INFO_IF("Mapper", EchoLevel > 0) << "Computed search-radius: "
        << search_radius << std::endl;

    return search_radius;
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

BoundingBoxType ComputeLocalBoundingBox(const ModelPart& rModelPart)
{
    BoundingBoxType local_bounding_box {-1e10, 1e10, -1e10, 1e10, -1e10, 1e10}; // initialize "inverted"
    // xmax, xmin,  ymax, ymin,  zmax, zmin

    // loop over all nodes (local and ghost(necessary if conditions have only ghost nodes) )
    for (auto &r_node : rModelPart.Nodes()) {
        local_bounding_box[0] = std::max(r_node.X(), local_bounding_box[0]);
        local_bounding_box[1] = std::min(r_node.X(), local_bounding_box[1]);
        local_bounding_box[2] = std::max(r_node.Y(), local_bounding_box[2]);
        local_bounding_box[3] = std::min(r_node.Y(), local_bounding_box[3]);
        local_bounding_box[4] = std::max(r_node.Z(), local_bounding_box[4]);
        local_bounding_box[5] = std::min(r_node.Z(), local_bounding_box[5]);
    }
    return local_bounding_box;
}

BoundingBoxType ComputeGlobalBoundingBox(const ModelPart& rModelPart)
{
    BoundingBoxType local_bounding_box = ComputeLocalBoundingBox(rModelPart);

    array_1d<double,3> max_vals;
    array_1d<double,3> min_vals;

    // fill buffers for MPI
    for (std::size_t i=0; i<3; ++i) {
        max_vals[i] = local_bounding_box[i*2];
        min_vals[i] = local_bounding_box[i*2+1];
    }

    // compute global values
    const auto& r_data_comm = rModelPart.GetCommunicator().GetDataCommunicator();
    max_vals = r_data_comm.MaxAll(max_vals);
    min_vals = r_data_comm.MinAll(min_vals);

    BoundingBoxType global_bounding_box;
    // extract information from buffers
    for (std::size_t i=0; i<3; ++i) {
        global_bounding_box[i*2] = max_vals[i];
        global_bounding_box[i*2+1] = min_vals[i];
    }

    return global_bounding_box;
}

void ComputeBoundingBoxesWithTolerance(const std::vector<double>& rBoundingBoxes,
                                       const double Tolerance,
                                       std::vector<double>& rBoundingBoxesWithTolerance)
{
    const SizeType size_vec = rBoundingBoxes.size();

    KRATOS_DEBUG_ERROR_IF_NOT(std::fmod(size_vec, 6) == 0)
        << "Bounding Boxes size has to be a multiple of 6!" << std::endl;

    if (rBoundingBoxesWithTolerance.size() != size_vec) {
        rBoundingBoxesWithTolerance.resize(size_vec);
    }

    // Apply Tolerances
    for (IndexType i=0; i<size_vec; i+=2) {
        rBoundingBoxesWithTolerance[i] = rBoundingBoxes[i] + Tolerance;
    }

    for (IndexType i=1; i<size_vec; i+=2) {
        rBoundingBoxesWithTolerance[i] = rBoundingBoxes[i] - Tolerance;
    }
}

std::string BoundingBoxStringStream(const BoundingBoxType& rBoundingBox)
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

bool PointIsInsideBoundingBox(const BoundingBoxType& rBoundingBox,
                              const array_1d<double, 3>& rCoords)
{   // The Bounding Box should have some tolerance already!
    if (rCoords[0] < rBoundingBox[0] && rCoords[0] > rBoundingBox[1])   // check x-direction
        if (rCoords[1] < rBoundingBox[2] && rCoords[1] > rBoundingBox[3])   // check y-direction
            if (rCoords[2] < rBoundingBox[4] && rCoords[2] > rBoundingBox[5])   // check z-direction
                return true;
    return false;
}

void CreateMapperLocalSystemsFromGeometries(const MapperLocalSystem& rMapperLocalSystemPrototype,
                                            const Communicator& rModelPartCommunicator,
                                            std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rLocalSystems)
{
    const std::size_t num_conditions = rModelPartCommunicator.LocalMesh().NumberOfConditions();
    const auto cond_begin = rModelPartCommunicator.LocalMesh().ConditionsBegin();

    if (rLocalSystems.size() != num_conditions) rLocalSystems.resize(num_conditions);

    IndexPartition<std::size_t>(num_conditions).for_each([&](std::size_t i){
        InterfaceObject::GeometryPointerType p_geom = &((cond_begin+i)->GetGeometry());
        rLocalSystems[i] = rMapperLocalSystemPrototype.Create(p_geom);
    });

    const int num_local_systems = rModelPartCommunicator.GetDataCommunicator().SumAll((int)(rLocalSystems.size())); // int bcs of MPI

    KRATOS_ERROR_IF_NOT(num_local_systems > 0) << "No mapper local systems were created" << std::endl;
}

void SaveCurrentConfiguration(ModelPart& rModelPart)
{
    KRATOS_TRY;

    block_for_each(rModelPart.Nodes(), [&](Node<3>& rNode){
        rNode.SetValue(CURRENT_COORDINATES, rNode.Coordinates());
    });

    KRATOS_CATCH("");
}

void RestoreCurrentConfiguration(ModelPart& rModelPart)
{
    KRATOS_TRY;

    if (rModelPart.NumberOfNodes() > 0) {
        KRATOS_ERROR_IF_NOT(rModelPart.NodesBegin()->Has(CURRENT_COORDINATES)) << "Nodes do not have CURRENT_COORDINATES for restoring the current configuration!" << std::endl;

        block_for_each(rModelPart.Nodes(), [&](Node<3>& rNode){
            noalias(rNode.Coordinates()) = rNode.GetValue(CURRENT_COORDINATES);
            rNode.Data().Erase(CURRENT_COORDINATES);
        });
    }

    KRATOS_CATCH("");
}


void FillBufferBeforeLocalSearch(const MapperLocalSystemPointerVector& rMapperLocalSystems,
                                 const std::vector<double>& rBoundingBoxes,
                                 const SizeType BufferSizeEstimate,
                                 std::vector<std::vector<double>>& rSendBuffer,
                                 std::vector<int>& rSendSizes)
{
    const SizeType comm_size = rSendBuffer.size();

    BoundingBoxType bounding_box;
    KRATOS_DEBUG_ERROR_IF_NOT(std::fmod(rBoundingBoxes.size(), 6) == 0)
        << "Bounding Boxes size has to be a multiple of 6!" << std::endl;

    // #pragma omp parallel for => "bounding_box" has to be threadprivate!
    for (IndexType i_rank=0; i_rank<comm_size; ++i_rank) {
        auto& r_rank_buffer = rSendBuffer[i_rank];
        r_rank_buffer.clear();
        r_rank_buffer.reserve(BufferSizeEstimate);
        rSendSizes[i_rank] = 0;

        for (IndexType i_local_sys=0; i_local_sys<rMapperLocalSystems.size(); ++i_local_sys) {
            for (IndexType j=0; j<6; ++j) {
                bounding_box[j] = rBoundingBoxes[(i_rank*6) + j]; // retrieve bounding box of partition
            }

            const auto& rp_local_sys = rMapperLocalSystems[i_local_sys];

            if (!rp_local_sys->HasInterfaceInfo()) {
                const auto& r_coords = rp_local_sys->Coordinates();
                if (MapperUtilities::PointIsInsideBoundingBox(bounding_box, r_coords)) {
                    // These push_backs are threadsafe bcs only one vector is accessed per thread!
                    r_rank_buffer.push_back(static_cast<double>(i_local_sys)); // this it the "mSourceLocalSystemIndex" of the MapperInterfaceInfo
                    r_rank_buffer.push_back(r_coords[0]);
                    r_rank_buffer.push_back(r_coords[1]);
                    r_rank_buffer.push_back(r_coords[2]);

                    rSendSizes[i_rank] += 4;
                }
            }
        }
    }
}

void CreateMapperInterfaceInfosFromBuffer(const std::vector<std::vector<double>>& rRecvBuffer,
                                          const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
                                          const int CommRank,
                                          MapperInterfaceInfoPointerVectorType& rMapperInterfaceInfosContainer)
{
    const SizeType comm_size = rMapperInterfaceInfosContainer.size();

    KRATOS_DEBUG_ERROR_IF_NOT(rRecvBuffer.size() == comm_size)
        << "Buffer-size mismatch!" << std::endl;

    array_1d<double, 3> coords;
    // Loop the ranks and construct the MapperInterfaceInfos
    for (IndexType i_rank=0; i_rank<comm_size; ++i_rank) {
        const SizeType recv_buffer_size_rank = rRecvBuffer[i_rank].size();
        KRATOS_DEBUG_ERROR_IF_NOT(std::fmod(recv_buffer_size_rank, 4) == 0) << "Rank " << CommRank
            << " received a wrong buffer-size from rank " << i_rank+1 << "!" << std::endl;
        const SizeType num_objs = recv_buffer_size_rank / 4; // 1 index and 3 coordinates

        const auto& r_rank_buffer = rRecvBuffer[i_rank];
        auto& r_interface_infos_rank = rMapperInterfaceInfosContainer[i_rank];
        r_interface_infos_rank.clear();

        if (r_interface_infos_rank.size() != num_objs) {
            r_interface_infos_rank.resize(num_objs);
        }

        for (IndexType j=0; j<num_objs; ++j) {
#ifdef KRATOS_DEBUG
            // with this check we make sure that this field
            // only contains doubles converted from ints
            // e.g. 4.5 is not allowed!
            double int_part;
            double fract_part = std::modf((r_rank_buffer[j*4]+0.1), &int_part);

            KRATOS_ERROR_IF(std::abs(fract_part-0.1) > 1e-12)
                << "Buffer contains a double (" << r_rank_buffer[j*4]
                << ") that was not casted from an int, i.e. it contains a "
                << "fractional part of " << std::abs(fract_part-0.1) << "!" << std::endl;
#endif
            // retrive data from buffer
            const int local_sys_idx = static_cast<IndexType>(r_rank_buffer[j*4]+0.1);
            // 0.1 is added to prevent truncation errors like (int)1.9999 = 1
            coords[0] = r_rank_buffer[j*4 + 1];
            coords[1] = r_rank_buffer[j*4 + 2];
            coords[2] = r_rank_buffer[j*4 + 3];
            r_interface_infos_rank[j] = rpRefInterfaceInfo->Create(coords, local_sys_idx, i_rank);
        }
    }
}

void FillBufferAfterLocalSearch(MapperInterfaceInfoPointerVectorType& rMapperInterfaceInfosContainer,
                                const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
                                const int CommRank,
                                std::vector<std::vector<char>>& rSendBuffer,
                                std::vector<int>& rSendSizes)
{
    const SizeType comm_size = rMapperInterfaceInfosContainer.size();

    KRATOS_DEBUG_ERROR_IF_NOT(rSendSizes.size() == comm_size)
        << "Buffer-size mismatch!" << std::endl;

    for (IndexType i_rank=0; i_rank<comm_size; ++i_rank) {
        if (i_rank != static_cast<IndexType>(CommRank)) {
            MapperUtilities::MapperInterfaceInfoSerializer interface_infos_serializer(
                rMapperInterfaceInfosContainer[i_rank], rpRefInterfaceInfo );

            Kratos::StreamSerializer serializer;
            serializer.save("interface_infos", interface_infos_serializer);

            const auto p_serializer_buffer = dynamic_cast<std::stringstream*>(serializer.pGetBuffer());
            const std::string& stream_str = p_serializer_buffer->str();

            const SizeType send_size = sizeof(char) * (stream_str.size()+1); // +1 fof Null-terminated string

            rSendSizes[i_rank] = send_size;

            auto& r_rank_buffer = rSendBuffer[i_rank];
            r_rank_buffer.clear();

            if (r_rank_buffer.size() != send_size) {
                r_rank_buffer.resize(send_size);
            }

            std::memcpy(r_rank_buffer.data(), stream_str.c_str(), send_size);
        }
    }
}

void DeserializeMapperInterfaceInfosFromBuffer(
    const std::vector<std::vector<char>>& rRecvBuffer,
    const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
    const int CommRank,
    MapperInterfaceInfoPointerVectorType& rMapperInterfaceInfosContainer)
{
    const SizeType comm_size = rMapperInterfaceInfosContainer.size();

    KRATOS_DEBUG_ERROR_IF_NOT(rRecvBuffer.size() == comm_size)
        << "Buffer-size mismatch!" << std::endl;

    for (IndexType i_rank=0; i_rank<comm_size; ++i_rank) {
        if (i_rank != static_cast<IndexType>(CommRank)) {
            Kratos::StreamSerializer serializer;

            const auto p_serializer_buffer = dynamic_cast<std::stringstream*>(serializer.pGetBuffer());
            p_serializer_buffer->write(rRecvBuffer[i_rank].data(), rRecvBuffer[i_rank].size());

            MapperUtilities::MapperInterfaceInfoSerializer interface_infos_serializer(
                rMapperInterfaceInfosContainer[i_rank], rpRefInterfaceInfo );

            serializer.load("interface_infos", interface_infos_serializer);
        }
    }
}

void MapperInterfaceInfoSerializer::save(Kratos::Serializer& rSerializer) const
{
    const SizeType num_infos = mrInterfaceInfos.size();
    rSerializer.save("size", num_infos);

    for (IndexType i=0; i<num_infos; ++i) {
        rSerializer.save("E", *(mrInterfaceInfos[i])); // NOT serializing the shared_ptr!
    }
}

void MapperInterfaceInfoSerializer::load(Kratos::Serializer& rSerializer)
{
    mrInterfaceInfos.clear(); // make sure it has no leftovers

    SizeType num_infos;
    rSerializer.load("size", num_infos);

    if (mrInterfaceInfos.size() != num_infos) {
        mrInterfaceInfos.resize(num_infos);
    }

    for (IndexType i=0; i<num_infos; ++i) {
        // first we create a new object, then we load its data
        // this is needed bcs of the polymorphic behavior of the InterfaceInfos
        // i.e. in order to create the correct type
        // => the vector contains baseclass-pointers!
        // Jordi I am quite sure that I could get around it by registering it, what do you think... TODO
        // I think doing it manually is more efficient, which I want so I would probably leave it ...
        // The serializer does some nasty casting when pointers are serialized...
        // I could do a benchmark at some point but I highly doubt that the serializer is faster ...
        mrInterfaceInfos[i] = mrpRefInterfaceInfo->Create();
        rSerializer.load("E", *(mrInterfaceInfos[i])); // NOT serializing the shared_ptr!
    }
}

} // namespace MapperUtilities
} // namespace Kratos.
