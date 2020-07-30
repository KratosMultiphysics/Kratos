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
    mThisParameters.ValidateAndAssignDefaults(this->GetDefaultSettings());
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void SkinDetectionProcess<TDim>::Execute()
{
    KRATOS_TRY;

    HashMapVectorIntType inverse_face_map;
    HashMapVectorIntIdsType properties_face_map;
    this->GenerateFaceMaps(inverse_face_map, properties_face_map);

    // Generate skin conditions
    ModelPart& r_work_model_part = this->SetUpAuxiliaryModelPart();
    this->FillAuxiliaryModelPart(r_work_model_part, inverse_face_map, properties_face_map);
    this->SetUpAdditionalSubModelParts(r_work_model_part);

    KRATOS_CATCH("");
}

template<SizeType TDim>
void SkinDetectionProcess<TDim>::GenerateFaceMaps(
    HashMapVectorIntType& rInverseFaceMap,
    HashMapVectorIntIdsType& rPropertiesFaceMap) const
{
    // Auxiliar values
    auto& r_elements_array = mrModelPart.Elements();
    const SizeType number_of_elements = r_elements_array.size();
    const auto it_elem_begin = r_elements_array.begin();

    for(IndexType i = 0; i < number_of_elements; ++i) {
        auto it_elem = it_elem_begin + i;

        // Detect if the element is active or not. If the user did not make any choice the element is active by default
        bool element_is_active = true;
        if (it_elem->IsDefined(ACTIVE))
            element_is_active = it_elem->Is(ACTIVE);

        if (element_is_active) {
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

        // Detect if the element is active or not. If the user did not make any choice the element is active by default
        bool element_is_active = true;
        if (it_elem->IsDefined(ACTIVE))
            element_is_active = it_elem->Is(ACTIVE);

        if (element_is_active) {
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

template<SizeType TDim>
ModelPart& SkinDetectionProcess<TDim>::SetUpAuxiliaryModelPart()
{
    // We create the auxiliar ModelPart
    const std::string& name_auxiliar_model_part = mThisParameters["name_auxiliar_model_part"].GetString();
    if (!(mrModelPart.HasSubModelPart(name_auxiliar_model_part))) {
        mrModelPart.CreateSubModelPart(name_auxiliar_model_part);
    } else {
        auto& r_conditions_array = mrModelPart.GetSubModelPart(name_auxiliar_model_part).Conditions();
        const SizeType number_of_conditions = r_conditions_array.size();
        const auto it_cond_begin = r_conditions_array.begin();

        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(number_of_conditions); ++i)
            (it_cond_begin + i)->Set(TO_ERASE, true);

        mrModelPart.GetSubModelPart(name_auxiliar_model_part).RemoveConditionsFromAllLevels(TO_ERASE);

        mrModelPart.RemoveSubModelPart(name_auxiliar_model_part);
        mrModelPart.CreateSubModelPart(name_auxiliar_model_part);
    }
    return mrModelPart.GetSubModelPart(name_auxiliar_model_part);
}

template<SizeType TDim>
void SkinDetectionProcess<TDim>::FillAuxiliaryModelPart(
    ModelPart& rAuxiliaryModelPart,
    HashMapVectorIntType& rInverseFaceMap,
    HashMapVectorIntIdsType& rPropertiesFaceMap)
{
    // The auxiliar name of the condition
    const std::string& name_condition = mThisParameters["name_auxiliar_condition"].GetString();
    std::string pre_name = "";
    if (TDim == 3 && name_condition == "Condition") {
        pre_name = "Surface";
    } else if (TDim == 2 && name_condition == "Condition") {
        pre_name = "Line";
    }
    const std::string base_name = pre_name + name_condition;

    // The number of conditions
    ConditionsArrayType& condition_array = mrModelPart.GetRootModelPart().Conditions();
    const auto& it_begin = condition_array.begin();
    for(IndexType i = 0; i < condition_array.size(); ++i)
        (it_begin + i)->SetId(i + 1);

    // The indexes of the nodes of the skin
    std::unordered_set<IndexType> nodes_in_the_skin;

    this->CreateConditions(mrModelPart, rAuxiliaryModelPart, rInverseFaceMap, rPropertiesFaceMap, nodes_in_the_skin, base_name);

    // Adding to the auxiliar model part
    VectorIndexType indexes_skin;
    indexes_skin.insert(indexes_skin.end(), nodes_in_the_skin.begin(), nodes_in_the_skin.end());
    rAuxiliaryModelPart.AddNodes(indexes_skin);

    const SizeType echo_level = mThisParameters["echo_level"].GetInt();
    KRATOS_INFO_IF("SkinDetectionProcess", echo_level > 0) << rInverseFaceMap.size() << " have been created" << std::endl;

    // Now we set the flag on the nodes. The list of nodes of the auxiliar model part
    auto& r_nodes_array = rAuxiliaryModelPart.Nodes();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
        auto it_node = r_nodes_array.begin() + i;
        it_node->Set(INTERFACE, true);
    }
}

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

    // Create the auxiliar conditions
    for (auto& map : rInverseFaceMap) {
        condition_id += 1;

        const VectorIndexType& nodes_face = map.second;

        Properties::Pointer p_prop;
        const IndexType property_id = rPropertiesFaceMap[map.first];
         if (rMainModelPart.RecursivelyHasProperties(property_id)) {
             p_prop = rMainModelPart.pGetProperties(property_id);
         } else {
             p_prop = rMainModelPart.CreateNewProperties(property_id);
         }

        for (auto& index : nodes_face)
            rNodesInTheSkin.insert(index);

        const std::string complete_name = rConditionName + std::to_string(TDim) + "D" + std::to_string(nodes_face.size()) + "N"; // If the condition doesn't follow this structure...sorry, we then need to modify this...
        auto p_cond = rMainModelPart.CreateNewCondition(complete_name, condition_id, nodes_face, p_prop);
        rSkinModelPart.AddCondition(p_cond);
        p_cond->Set(INTERFACE, true);
        p_cond->Initialize();
    }
}

template<SizeType TDim>
void SkinDetectionProcess<TDim>::SetUpAdditionalSubModelParts(const ModelPart& rAuxiliaryModelPart)
{
    // We detect the conditions in the boundary model parts
    const SizeType n_model_parts = mThisParameters["list_model_parts_to_assign_conditions"].size();
    if (n_model_parts > 0) {

        // We build a database of indexes
        std::unordered_map<IndexType, std::unordered_set<IndexType>> conditions_nodes_ids_map;

        for (auto& cond : rAuxiliaryModelPart.Conditions()) {
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
}


/***********************************************************************************/
/***********************************************************************************/
template<SizeType TDim>
Parameters SkinDetectionProcess<TDim>::GetDefaultSettings() const
{
    Parameters default_parameters = Parameters(R"(
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
