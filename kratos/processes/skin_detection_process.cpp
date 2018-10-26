//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "processes/skin_detection_process.h"

namespace Kratos
{
template<SizeType TDim>
SkinDetectionProcess<TDim>::SkinDetectionProcess(
    ModelPart& rModelPart,
    Parameters ThisParameters
    ) : mrModelPart(rModelPart),
        mThisParameters(ThisParameters)
{
    Parameters default_parameters = Parameters(R"(
    {
        "name_auxiliar_model_part"              : "SkinModelPart",
        "name_auxiliar_condition"               : "Condition",
        "list_model_parts_to_assign_conditions" : [],
        "echo_level"                            : 0
    })" );

    mThisParameters.ValidateAndAssignDefaults(default_parameters);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void SkinDetectionProcess<TDim>::Execute()
{
    KRATOS_TRY;

    // Auxiliar values
    const SizeType number_of_elements = mrModelPart.Elements().size();
    const SizeType echo_level = mThisParameters["echo_level"].GetInt();

    /* NEIGHBOUR ELEMENTS */
    // Create the inverse_face_map
    HashMapVectorIntType inverse_face_map;
    HashMapVectorIntIdsType properties_face_map;

    for(IndexType i = 0; i < number_of_elements; ++i) {
        auto it_elem = mrModelPart.Elements().begin() + i;

        GeometryType& geom = it_elem->GetGeometry();

        const SizeType potential_number_neighbours = ComputePotentialNeighboursSize(it_elem);

        for (IndexType i_face = 0; i_face < potential_number_neighbours; ++i_face) {

            /* FACES/EDGES */
            const SizeType number_nodes = TDim == 2 ? geom.Edges()[i_face].size() : geom.Faces()[i_face].size();
            VectorIndexType vector_ids(number_nodes);
            VectorIndexType ordered_vector_ids(number_nodes);

            /* FACE/EDGE */
            if (TDim == 2) {
                for (IndexType i_node = 0; i_node < number_nodes; ++i_node) {
                    vector_ids[i_node] = geom.Edges()[i_face][i_node].Id();
                    ordered_vector_ids[i_node] = vector_ids[i_node];
                }
            } else {
                for (IndexType i_node = 0; i_node < number_nodes; ++i_node) {
                    vector_ids[i_node] = geom.Faces()[i_face][i_node].Id();
                    ordered_vector_ids[i_node] = vector_ids[i_node];
                }
            }

            /*** THE ARRAY OF IDS MUST BE ORDERED!!! ***/
            std::sort(vector_ids.begin(), vector_ids.end());
            // Check if the elements already exist in the HashMapVectorIntType
            HashMapVectorIntTypeIteratorType it_check = inverse_face_map.find(vector_ids);

            if(it_check == inverse_face_map.end() ) {
                // If it doesn't exist it is added to the database
                inverse_face_map.insert(std::pair<VectorIndexType, VectorIndexType>(vector_ids, ordered_vector_ids));
                properties_face_map.insert(std::pair<VectorIndexType, IndexType>(vector_ids, (it_elem->pGetProperties())->Id()));
            }
        }
    }

    // Create the face_set
    HashSetVectorIntType face_set;

    for(IndexType i = 0; i < number_of_elements; ++i) {
        auto it_elem = mrModelPart.Elements().begin() + i;

        GeometryType& geom = it_elem->GetGeometry();

        const SizeType potential_number_neighbours = ComputePotentialNeighboursSize(it_elem);

        for (IndexType i_face = 0; i_face < potential_number_neighbours; ++i_face) {

            /* FACES/EDGES */
            const SizeType number_nodes = TDim == 2 ? geom.Edges()[i_face].size() : geom.Faces()[i_face].size();
            VectorIndexType vector_ids(number_nodes);

            /* FACE/EDGE */
            if (TDim == 2) {
                for (IndexType i_node = 0; i_node < number_nodes; ++i_node) {
                    vector_ids[i_node] = geom.Edges()[i_face][i_node].Id();
                }
            } else {
                for (IndexType i_node = 0; i_node < number_nodes; ++i_node) {
                    vector_ids[i_node] = geom.Faces()[i_face][i_node].Id();
                }
            }

            /*** THE ARRAY OF IDS MUST BE ORDERED!!! ***/
            std::sort(vector_ids.begin(), vector_ids.end());
            // Check if the elements already exist in the HashSetVectorIntType
            HashSetVectorIntTypeIteratorType it_check = face_set.find(vector_ids);

            if(it_check != face_set.end() ) {
                // If it exists we remove from the inverse map
                inverse_face_map.erase(vector_ids);
                properties_face_map.erase(vector_ids);
            } else {
                // If it doesn't exist it is added to the database
                face_set.insert(vector_ids);
            }
        }
    }

    // We create the auxiliar ModelPart
    const std::string& name_auxiliar_model_part = mThisParameters["name_auxiliar_model_part"].GetString();
    if (!(mrModelPart.HasSubModelPart(name_auxiliar_model_part))) {
        mrModelPart.CreateSubModelPart(name_auxiliar_model_part);
    } else {
        auto& conditions_array = mrModelPart.GetSubModelPart(name_auxiliar_model_part).Conditions();

        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(conditions_array.size()); ++i)
            (conditions_array.begin() + i)->Set(TO_ERASE, true);

        mrModelPart.GetSubModelPart(name_auxiliar_model_part).RemoveConditionsFromAllLevels(TO_ERASE);

        mrModelPart.RemoveSubModelPart(name_auxiliar_model_part);
        mrModelPart.CreateSubModelPart(name_auxiliar_model_part);
    } 
    ModelPart& r_auxiliar_model_part = mrModelPart.GetSubModelPart(name_auxiliar_model_part);

    // The auxiliar name of the condition
    const std::string& name_condition = mThisParameters["name_auxiliar_condition"].GetString();
    std::string pre_name = "";
    if (TDim == 3 && name_condition == "Condition") 
        pre_name = "Surface";

    // The number of conditions
    IndexType condition_id = mrModelPart.GetRootModelPart().Conditions().size();

    // The indexes of the nodes of the skin
    std::unordered_set<IndexType> nodes_in_the_skin;

    // Create the auxiliar conditions
    for (auto& map : inverse_face_map) {
        condition_id += 1;

        const VectorIndexType& nodes_face = map.second;
        Properties::Pointer p_prop = mrModelPart.pGetProperties(properties_face_map[map.first]);

        for (auto& index : nodes_face)
            nodes_in_the_skin.insert(index);

        const std::string complete_name = pre_name + name_condition + std::to_string(TDim) + "D" + std::to_string(nodes_face.size()) + "N"; // If the condition doesn't follow this structure...sorry, we then need to modify this...
        auto p_cond = mrModelPart.CreateNewCondition(complete_name, condition_id, nodes_face, p_prop);
        r_auxiliar_model_part.AddCondition(p_cond);
        p_cond->Set(INTERFACE, true);
        p_cond->Initialize();
    }

    // Adding to the auxiliar model part
    VectorIndexType indexes_skin;
    indexes_skin.insert(indexes_skin.end(), nodes_in_the_skin.begin(), nodes_in_the_skin.end());
    r_auxiliar_model_part.AddNodes(indexes_skin);

    KRATOS_INFO_IF("SkinDetectionProcess", echo_level > 0) << inverse_face_map.size() << " have been created" << std::endl;

    // Now we set the falg on the nodes. The list of nodes of the auxiliar model part
    auto& nodes_array = r_auxiliar_model_part.Nodes();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
        auto it_node = nodes_array.begin() + i;
        it_node->Set(INTERFACE, true);
    }

    // We detect the conditions in the boundary model parts
    const SizeType n_model_parts = mThisParameters["list_model_parts_to_assign_conditions"].size();
    if (n_model_parts > 0) {

        // We build a database of indexes
        std::unordered_map<IndexType, std::unordered_set<IndexType>> conditions_nodes_ids_map;

        for (auto& cond : r_auxiliar_model_part.Conditions()) {
            auto& geom = cond.GetGeometry();

            for (auto& node : geom) {
                auto set = conditions_nodes_ids_map.find(node.Id());
                if(set != conditions_nodes_ids_map.end()) {
                    conditions_nodes_ids_map[node.Id()].insert(cond.Id());
                } else {
                    std::unordered_set<IndexType> cond_index_ids ( {cond.Id()} );;
                    conditions_nodes_ids_map.insert({node.Id(), cond_index_ids});
                }
            }
        }

        ModelPart& root_model_part = mrModelPart.GetRootModelPart();
        for (IndexType i_mp = 0; i_mp < n_model_parts; ++i_mp){
            const std::string& model_part_name = mThisParameters["list_model_parts_to_assign_conditions"].GetArrayItem(i_mp).GetString();
            ModelPart& sub_model_part = root_model_part.GetSubModelPart(model_part_name);

            std::vector<IndexType> conditions_ids;

            #pragma omp parallel
            {
                // Creating a buffer for parallel vector fill
                std::vector<IndexType> conditions_ids_buffer;

                // We iterate over the nodes of this model part
                auto& sub_nodes_array = sub_model_part.Nodes();
                #pragma omp for
                for(int i = 0; i < static_cast<int>(sub_nodes_array.size()); ++i) {
                    auto it_node = sub_nodes_array.begin() + i;

                    auto set = conditions_nodes_ids_map.find(it_node->Id());
                    if(set != conditions_nodes_ids_map.end()) {
                        for (auto& cond_id : conditions_nodes_ids_map[it_node->Id()]) {
                            auto& r_condition = mrModelPart.GetCondition(cond_id);
                            auto& geom = r_condition.GetGeometry();
                            bool has_nodes = true;
                            for (auto& node : geom) {
                                if (!sub_model_part.GetMesh().HasNode(node.Id())) {
                                    has_nodes = false;
                                    break;
                                }
                            }
                            // We append to the vector
                            if (has_nodes) conditions_ids_buffer.push_back(r_condition.Id());
                        }
                    }
                }

                // Combine buffers together
                #pragma omp critical
                {
                    std::move(conditions_ids_buffer.begin(),conditions_ids_buffer.end(),back_inserter(conditions_ids));
                }
            }

            sub_model_part.AddConditions(conditions_ids);
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
SizeType SkinDetectionProcess<2>::ComputePotentialNeighboursSize(ElementsIteratorType itElem)
{
    const auto& geometry = itElem->GetGeometry();
    return geometry.EdgesNumber();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
SizeType SkinDetectionProcess<3>::ComputePotentialNeighboursSize(ElementsIteratorType itElem)
{
    const auto& geometry = itElem->GetGeometry();
    return geometry.FacesNumber();
}

/***********************************************************************************/
/***********************************************************************************/

template class SkinDetectionProcess<2>;
template class SkinDetectionProcess<3>;
// class SkinDetectionProcess

} // namespace Kratos
