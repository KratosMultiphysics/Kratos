//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
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
#include "processes/skin_detection_process.h"
#include "utilities/variable_utils.h"
#include "includes/key_hash.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

namespace Kratos
{
template<SizeType TDim>
SkinDetectionProcess<TDim>::SkinDetectionProcess(
    ModelPart& rModelPart,
    Parameters ThisParameters
    ) : mrModelPart(rModelPart),
        mThisParameters(ThisParameters)
{
    mThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void SkinDetectionProcess<TDim>::Execute()
{
    KRATOS_TRY;

    // First assign MPI ids if needed
    std::unordered_set<IndexType> set_node_ids_interface;
    this->GenerateSetNodeIdsInterface(set_node_ids_interface);

    // Generate face maps
    HashMapVectorIntType inverse_face_map;
    HashMapVectorIntIdsType properties_face_map;
    this->GenerateFaceMaps(inverse_face_map, properties_face_map);

    // Filter local nodes
    this->FilterMPIInterfaceNodes(set_node_ids_interface, inverse_face_map);

    // Generate skin conditions
    ModelPart& r_work_model_part = this->SetUpAuxiliaryModelPart();
    this->FillAuxiliaryModelPart(r_work_model_part, inverse_face_map, properties_face_map);
    this->SetUpAdditionalSubModelParts(r_work_model_part);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void SkinDetectionProcess<TDim>::GenerateFaceMaps(
    HashMapVectorIntType& rInverseFaceMap,
    HashMapVectorIntIdsType& rPropertiesFaceMap
    ) const
{
    // Auxiliary values
    auto& r_elements_array = mrModelPart.Elements();
    const SizeType number_of_elements = r_elements_array.size();
    const auto it_elem_begin = r_elements_array.begin();

    for(IndexType i = 0; i < number_of_elements; ++i) {
        auto it_elem = it_elem_begin + i;

        // If the element is active
        if (it_elem->IsActive()) {
            GeometryType& r_geometry = it_elem->GetGeometry();

            const auto r_boundary_geometries = r_geometry.GenerateBoundariesEntities();
            const SizeType potential_number_neighbours = r_boundary_geometries.size();

            for (IndexType i_face = 0; i_face < potential_number_neighbours; ++i_face) {

                /* FACES/EDGES */
                const SizeType number_nodes = r_boundary_geometries[i_face].size();
                VectorIndexType vector_ids(number_nodes);
                VectorIndexType ordered_vector_ids(number_nodes);

                /* FACE/EDGE */
                for (IndexType i_node = 0; i_node < number_nodes; ++i_node) {
                    vector_ids[i_node] = r_boundary_geometries[i_face][i_node].Id();
                    ordered_vector_ids[i_node] = vector_ids[i_node];
                }

                /*** THE ARRAY OF IDS MUST BE ORDERED!!! ***/
                std::sort(vector_ids.begin(), vector_ids.end());
                // Check if the elements already exist in the HashMapVectorIntType
                HashMapVectorIntTypeIteratorType it_check = rInverseFaceMap.find(vector_ids);

                if(it_check == rInverseFaceMap.end() ) {
                    // If it doesn't exist it is added to the database
                    rInverseFaceMap.insert(std::pair<VectorIndexType, VectorIndexType>(vector_ids, ordered_vector_ids));
                    rPropertiesFaceMap.insert(std::pair<VectorIndexType, IndexType>(vector_ids, (it_elem->pGetProperties())->Id()));
                }
            }
        }
    }

    // Create the face_set
    HashSetVectorIntType face_set;

    for(IndexType i = 0; i < number_of_elements; ++i) {
        auto it_elem = it_elem_begin + i;

        // If the element is active
        if (it_elem->IsActive()) {
            GeometryType& r_geometry = it_elem->GetGeometry();

            const auto r_boundary_geometries = r_geometry.GenerateBoundariesEntities();
            const SizeType potential_number_neighbours = r_boundary_geometries.size();

            for (IndexType i_face = 0; i_face < potential_number_neighbours; ++i_face) {

                /* FACES/EDGES */
                const SizeType number_nodes = r_boundary_geometries[i_face].size();
                VectorIndexType vector_ids(number_nodes);

                /* FACE/EDGE */
                for (IndexType i_node = 0; i_node < number_nodes; ++i_node) {
                    vector_ids[i_node] = r_boundary_geometries[i_face][i_node].Id();
                }

                /*** THE ARRAY OF IDS MUST BE ORDERED!!! ***/
                std::sort(vector_ids.begin(), vector_ids.end());
                // Check if the elements already exist in the HashSetVectorIntType
                HashSetVectorIntTypeIteratorType it_check = face_set.find(vector_ids);

                if(it_check != face_set.end() ) {
                    // If it exists we remove from the inverse map
                    rInverseFaceMap.erase(vector_ids);
                    rPropertiesFaceMap.erase(vector_ids);
                } else {
                    // If it doesn't exist it is added to the database
                    face_set.insert(vector_ids);
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
ModelPart& SkinDetectionProcess<TDim>::SetUpAuxiliaryModelPart()
{
    // We create the auxiliary ModelPart
    const std::string& name_auxiliar_model_part = mThisParameters["name_auxiliar_model_part"].GetString();
    if (!(mrModelPart.HasSubModelPart(name_auxiliar_model_part))) {
        mrModelPart.CreateSubModelPart(name_auxiliar_model_part);
    } else {
        auto& r_conditions_array = mrModelPart.GetSubModelPart(name_auxiliar_model_part).Conditions();

        VariableUtils().SetFlag(TO_ERASE, true, r_conditions_array);

        mrModelPart.GetSubModelPart(name_auxiliar_model_part).RemoveConditionsFromAllLevels(TO_ERASE);

        mrModelPart.RemoveSubModelPart(name_auxiliar_model_part);
        mrModelPart.CreateSubModelPart(name_auxiliar_model_part);
    }
    return mrModelPart.GetSubModelPart(name_auxiliar_model_part);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void SkinDetectionProcess<TDim>::FillAuxiliaryModelPart(
    ModelPart& rAuxiliaryModelPart,
    HashMapVectorIntType& rInverseFaceMap,
    HashMapVectorIntIdsType& rPropertiesFaceMap
    )
{
    // The auxiliary name of the condition
    const std::string& r_name_condition = mThisParameters["name_auxiliar_condition"].GetString();
    std::string pre_name = "";
    if (r_name_condition == "Condition") {
        if constexpr (TDim == 3) {
            pre_name = "Surface";
        } else if constexpr (TDim == 2) {
            pre_name = "Line";
        }
    }
    const std::string base_name = pre_name + r_name_condition;

    // The number of conditions
    auto& r_condition_array = mrModelPart.GetRootModelPart().Conditions();
    const auto it_cond_begin = r_condition_array.begin();
    IndexPartition<IndexType>(r_condition_array.size()).for_each([&it_cond_begin](IndexType i) {
        (it_cond_begin + i)->SetId(i + 1);
    });

    // The indexes of the nodes of the skin
    std::unordered_set<IndexType> nodes_in_the_skin;

    this->CreateConditions(mrModelPart, rAuxiliaryModelPart, rInverseFaceMap, rPropertiesFaceMap, nodes_in_the_skin, base_name);

    // Adding to the auxiliary model part
    VectorIndexType indexes_skin;
    indexes_skin.insert(indexes_skin.end(), nodes_in_the_skin.begin(), nodes_in_the_skin.end());
    rAuxiliaryModelPart.AddNodes(indexes_skin);

    const SizeType echo_level = mThisParameters["echo_level"].GetInt();
    KRATOS_INFO_IF("SkinDetectionProcess", echo_level > 0) << rInverseFaceMap.size() << " have been created" << std::endl;

    // Verify that nodes are in the corresponding model part
    auto& r_nodes_array = rAuxiliaryModelPart.Nodes();
    const auto& r_data_communicator = mrModelPart.GetCommunicator().GetDataCommunicator();
    if (r_data_communicator.IsDistributed()) {
        // Wait until all the nodes are created
        r_data_communicator.Barrier();

        // First look into potential nodes to ensure
        const int rank = r_data_communicator.Rank();
        const int world_size = r_data_communicator.Size();
        std::vector<std::vector<std::size_t>> nodes_to_ensure(world_size);
        for (auto& r_node : r_nodes_array) {
            const int partition_index = r_node.FastGetSolutionStepValue(PARTITION_INDEX);
            if (rank != partition_index) {
                nodes_to_ensure[partition_index].push_back(r_node.Id());
            }
        }

        // Now we ensure the nodes
        int tag_send = 0;
        int tag_recv = 0;
        std::vector<std::size_t> nodes_to_bring;
        for (int i = 0; i < world_size; ++i) {
            if (i != rank) {
                // Compute the tag hashing origin rank and destination rank
                tag_send = 1000 * rank + i;
                r_data_communicator.Send(nodes_to_ensure[i], i, tag_send);
            } else {
                for (int j = 0; j < world_size; ++j) {
                    if (j != rank) {
                        // Compute the tag hashing origin rank and destination rank
                        tag_recv = 1000 * j + rank;
                        std::vector<std::size_t> rec_nodes_to_ensure;
                        r_data_communicator.Recv(rec_nodes_to_ensure, j, tag_recv);
                        for (auto node_id : rec_nodes_to_ensure) {
                            if (!rAuxiliaryModelPart.HasNode(node_id)) {
                                nodes_to_bring.push_back(node_id);
                            }
                        }
                    }
                }
            }
        }
        rAuxiliaryModelPart.AddNodes(nodes_to_bring);
    }

    // Now we set the flag on the nodes. The list of nodes of the auxiliary model part
    VariableUtils().SetFlag(INTERFACE, true, r_nodes_array);

    // In case we are in MPI we synchronize the INTERFACE flag
    mrModelPart.GetCommunicator().SynchronizeOrNodalFlags(INTERFACE);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void SkinDetectionProcess<TDim>::CreateConditions(
    ModelPart& rMainModelPart,
    ModelPart& rSkinModelPart,
    HashMapVectorIntType& rInverseFaceMap,
    HashMapVectorIntIdsType& rPropertiesFaceMap,
    std::unordered_set<IndexType>& rNodesInTheSkin,
    const std::string& rConditionName) const
{

    IndexType condition_id = rMainModelPart.GetRootModelPart().Conditions().size();
    const auto& r_process_info = rMainModelPart.GetProcessInfo();

    // Create the auxiliary conditions
    for (auto& r_map : rInverseFaceMap) {
        condition_id += 1;

        const VectorIndexType& r_nodes_face = r_map.second;

        Properties::Pointer p_prop = nullptr;
        const IndexType property_id = rPropertiesFaceMap[r_map.first];
         if (rMainModelPart.RecursivelyHasProperties(property_id)) {
             p_prop = rMainModelPart.pGetProperties(property_id);
         } else {
             p_prop = rMainModelPart.CreateNewProperties(property_id);
         }

        for (auto& r_index : r_nodes_face) {
            rNodesInTheSkin.insert(r_index);
        }

        const std::string complete_name = rConditionName + std::to_string(TDim) + "D" + std::to_string(r_nodes_face.size()) + "N"; // If the condition doesn't follow this structure...sorry, we then need to modify this...
        auto p_cond = rMainModelPart.CreateNewCondition(complete_name, condition_id, r_nodes_face, p_prop);
        rSkinModelPart.AddCondition(p_cond);
        p_cond->Set(INTERFACE, true);
        p_cond->Initialize(r_process_info);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void SkinDetectionProcess<TDim>::SetUpAdditionalSubModelParts(const ModelPart& rAuxiliaryModelPart)
{
    // We detect the conditions in the boundary model parts
    const SizeType n_model_parts = mThisParameters["list_model_parts_to_assign_conditions"].size();
    if (n_model_parts > 0) {

        // We build a database of indexes
        std::unordered_map<IndexType, std::unordered_set<IndexType>> conditions_nodes_ids_map;

        for (auto& r_cond : rAuxiliaryModelPart.Conditions()) {
            auto& r_geom = r_cond.GetGeometry();

            for (auto& r_node : r_geom) {
                auto it_set_found = conditions_nodes_ids_map.find(r_node.Id());
                if(it_set_found != conditions_nodes_ids_map.end()) {
                    conditions_nodes_ids_map[r_node.Id()].insert(r_cond.Id());
                } else {
                    std::unordered_set<IndexType> cond_index_ids ( {r_cond.Id()} );;
                    conditions_nodes_ids_map.insert({r_node.Id(), cond_index_ids});
                }
            }
        }

        ModelPart& r_root_model_part = mrModelPart.GetRootModelPart();
        for (IndexType i_mp = 0; i_mp < n_model_parts; ++i_mp) {
            const std::string& r_model_part_name = mThisParameters["list_model_parts_to_assign_conditions"].GetArrayItem(i_mp).GetString();
            ModelPart& r_sub_model_part = r_root_model_part.GetSubModelPart(r_model_part_name);

            // We iterate over the nodes of this model part
            auto& r_sub_nodes_array = r_sub_model_part.Nodes();
            const auto it_node_begin = r_sub_nodes_array.begin();
            std::vector<std::vector<IndexType>> conditions_ids_aux = IndexPartition<IndexType>(r_sub_nodes_array.size()).for_each<AccumReduction<std::vector<IndexType>>>([&](IndexType i) {
                auto it_node = it_node_begin + i;

                // Creating a buffer for parallel vector fill
                std::vector<IndexType> conditions_ids_buffer;

                auto it_set_found = conditions_nodes_ids_map.find(it_node->Id());
                if(it_set_found != conditions_nodes_ids_map.end()) {
                    for (auto& r_cond_id : conditions_nodes_ids_map[it_node->Id()]) {
                        auto& r_condition = mrModelPart.GetCondition(r_cond_id);
                        auto& r_geom = r_condition.GetGeometry();
                        bool has_nodes = true;
                        for (auto& r_node : r_geom) {
                            if (!r_sub_model_part.GetMesh().HasNode(r_node.Id())) {
                                has_nodes = false;
                                break;
                            }
                        }
                        // We append to the vector
                        if (has_nodes) conditions_ids_buffer.push_back(r_condition.Id());
                    }
                }

                return conditions_ids_buffer;
            });

            // We flatten the vector
            std::vector<IndexType> conditions_ids;

            // Calculate total size needed for conditions_ids to avoid reallocations
            const std::size_t total_size = block_for_each<SumReduction<std::size_t>>(conditions_ids_aux, [](const auto& r_vec) {
                return r_vec.size();
            });

            // Reserve memory for conditions_ids
            conditions_ids.reserve(total_size);

            // Flatten the vector of vectors
            for (const auto& r_vec : conditions_ids_aux) {
                conditions_ids.insert(conditions_ids.end(), r_vec.begin(), r_vec.end());
            }

            r_sub_model_part.AddConditions(conditions_ids);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void SkinDetectionProcess<TDim>::GenerateSetNodeIdsInterface(std::unordered_set<IndexType>& rSetNodeIdsInterface)
{
    // Only in case we are executing in MPI
    if (mrModelPart.IsDistributed()) {
        auto& r_communicator = mrModelPart.GetCommunicator();
        // const auto rank = r_communicator.GetDataCommunicator().Rank();
        const auto& r_nodes_interface = r_communicator.InterfaceMesh().Nodes();
        const std::size_t number_of_interface_nodes = r_nodes_interface.size();
        const auto it_interface_node_begin = r_nodes_interface.begin();
        std::vector<IndexType> node_ids_interface(number_of_interface_nodes);
        IndexPartition<std::size_t>(number_of_interface_nodes).for_each(
        [&node_ids_interface, it_interface_node_begin](std::size_t i) {
            auto it_interface_node = it_interface_node_begin + i;
            node_ids_interface[i] = it_interface_node->Id();
        });
        std::copy(node_ids_interface.begin(), node_ids_interface.end(), std::inserter(rSetNodeIdsInterface, rSetNodeIdsInterface.end()));
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void SkinDetectionProcess<TDim>::FilterMPIInterfaceNodes(
    const std::unordered_set<IndexType>& rSetNodeIdsInterface,
    HashMapVectorIntType& rInverseFaceMap
    )
{
    // First determine with the nodes in the MPI interface which faces are potentially removable
    std::vector<VectorIndexType> faces_to_remove;
    bool to_remove;
    for (auto& r_map : rInverseFaceMap) {
        to_remove = true;
        const VectorIndexType& r_vector_ids = r_map.first;
        const VectorIndexType& r_nodes_face = r_map.second;
        for (auto& r_index : r_nodes_face) {
            if (rSetNodeIdsInterface.find(r_index) == rSetNodeIdsInterface.end()) {
                to_remove = false;
                continue;
            }
        }
        if (to_remove) {
            faces_to_remove.push_back(r_vector_ids);
        }
    }

    /* Not all the faces are going to be removed, only the ones which are repeated in different processes. So we need to filter then. */

    // First we determine the rank and the size of the world
    const auto& r_communicator = mrModelPart.GetCommunicator();
    const auto& r_data_communicator = r_communicator.GetDataCommunicator();
    // const auto rank = r_data_communicator.Rank();
    // const auto world_size = r_data_communicator.Size();
    const auto& r_neighbour_indices = r_communicator.NeighbourIndices();
    std::vector<int> neighbour_indices;
    for (auto& r_index : r_neighbour_indices) {
        if (r_index >= 0) {
            neighbour_indices.push_back(r_index);
        }
    }

    // We define a scope, so everything will be removed at the end, except the clean up of the faces_to_remove
    {
        // Define the send tag
        const int tag_send = 1;

        // We generate the hash of the faces to use the communicator to send
        std::unordered_map<std::size_t, bool> faces_mpi_counter;
        std::unordered_map<std::size_t, VectorIndexType> faces_hash_map;
        std::vector<std::size_t> faces_to_remove_hash;
        VectorIndexHasher<std::vector<std::size_t>> vector_hasher;
        faces_to_remove_hash.reserve(faces_to_remove.size());
        for (auto& r_face_to_remove : faces_to_remove) {
            const std::size_t hash_face = vector_hasher(r_face_to_remove);
            faces_to_remove_hash.push_back(hash_face);
            faces_mpi_counter.insert({hash_face, false});
            faces_hash_map.insert({hash_face, r_face_to_remove});
        }

        // We send the list of the faces to be removed from the current ramk to the neighbour ranks
        for (auto& r_destination_rank : neighbour_indices) {
            // We send the faces_to_remove_hash to the other processes
            r_data_communicator.Send(faces_to_remove_hash, r_destination_rank, tag_send);
        }

        // Now we receive the faces_to_remove from the rest of the processes
        for (auto& r_origin_rank : neighbour_indices) {
            std::vector<std::size_t> rec_faces_to_remove_hash;
            r_data_communicator.Recv(rec_faces_to_remove_hash, r_origin_rank, tag_send);

            // Update the faces to be removed
            for (auto& r_hash_hash : rec_faces_to_remove_hash) {
                auto it_find_face = faces_mpi_counter.find(r_hash_hash);
                if (it_find_face != faces_mpi_counter.end()) {
                    it_find_face->second = true;
                }
            }
        }

        // Now we create the vector of ids to be removed
        std::vector<std::size_t> final_faces_to_remove;
        for (auto& r_face_pair : faces_mpi_counter) {
            if (r_face_pair.second) {
                final_faces_to_remove.push_back(r_face_pair.first);
            }
        }

        // Finally filter the faces to be removed
        faces_to_remove.clear();
        for (auto& r_face_to_remove : final_faces_to_remove) {
            auto it_find_face = faces_hash_map.find(r_face_to_remove);
            if (it_find_face != faces_hash_map.end()) {
                faces_to_remove.push_back(it_find_face->second);
            }
        }
    }

    // Finally we remove the MPi interface faces
    for (auto& r_face_to_remove : faces_to_remove) {
        rInverseFaceMap.erase(r_face_to_remove);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
const Parameters SkinDetectionProcess<TDim>::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "name_auxiliar_model_part"              : "SkinModelPart",
        "name_auxiliar_condition"               : "Condition",
        "list_model_parts_to_assign_conditions" : [],
        "echo_level"                            : 0
    })" );
    return default_parameters;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
ModelPart& SkinDetectionProcess<TDim>::GetModelPart() const
{
    return this->mrModelPart;
}

/***********************************************************************************/
/***********************************************************************************/
template<SizeType TDim>
Parameters SkinDetectionProcess<TDim>::GetSettings() const
{
    return this->mThisParameters;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
SkinDetectionProcess<TDim>::SkinDetectionProcess(
    ModelPart& rModelPart,
    Parameters Settings,
    Parameters DefaultSettings)
    : Process()
    , mrModelPart(rModelPart)
    , mThisParameters(Settings)
{
    mThisParameters.ValidateAndAssignDefaults(DefaultSettings);
}

/***********************************************************************************/
/***********************************************************************************/

template class SkinDetectionProcess<2>;
template class SkinDetectionProcess<3>;
// class SkinDetectionProcess

} // namespace Kratos
