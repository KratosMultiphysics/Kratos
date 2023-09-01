//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes

// External includes

// Project includes
#include "includes/parallel_environment.h"
#include "utilities/entities_utilities.h"
#include "utilities/variable_utils.h"
#include "create_entities_from_geometries_modeler.h"

namespace Kratos
{

namespace {

/**
 * @brief Create entities from geometries.
 * @param EntityName the name of the entity
 * @param rModelPart the model part
 */
template <class TEntitiesContainerType, class TEntityIdentifier>
void CreateEntitiesFromGeometries(
    const TEntityIdentifier& rEntityIdentifier,
    ModelPart& rModelPart
    )
{
    // Create the entities container and allocate space
    TEntitiesContainerType entities_to_add;
    entities_to_add.reserve(rModelPart.NumberOfGeometries());

    // Get current max element id
    using SizeType = std::size_t;
    SizeType max_id;
    const auto& r_root_model_part = rModelPart.GetRootModelPart();
    if constexpr (std::is_same<typename TEntitiesContainerType::value_type, Element>::value) {
        max_id = block_for_each<MaxReduction<SizeType>>(r_root_model_part.Elements(), [](auto& rElement){
            return rElement.Id();
        });
    } else {
        max_id = block_for_each<MaxReduction<SizeType>>(r_root_model_part.Conditions(), [](auto& rCondition){
            return rCondition.Id();
        });
    }

    // Loop geometries to create the corresponding entities from them
    // Note that we retrieve the corresponding prototype entity from the entities idenfifier
    // This makes possible to loop and create entities from different type geometries
    for (auto& r_geom : rModelPart.Geometries()) {
        if (rEntityIdentifier.HasPrototypeEntity(r_geom)) {
            const auto& r_ref_entity = rEntityIdentifier.GetPrototypeEntity(r_geom);
            auto p_entity = r_ref_entity.Create(++max_id, r_geom, nullptr);
            entities_to_add.push_back(p_entity);
        }
    }

    // Add the created entities to current submodelpart
    if constexpr (std::is_same<typename TEntitiesContainerType::value_type, Element>::value) {
        rModelPart.AddElements(entities_to_add.begin(), entities_to_add.end());
    } else {
        rModelPart.AddConditions(entities_to_add.begin(), entities_to_add.end());
    }
}

}

template <>
void CreateEntitiesFromGeometriesModeler::RemoveModelPartEntities<Element>(ModelPart &rModelPart)
{
    const SizeType n_elements = rModelPart.NumberOfElements();
    KRATOS_WARNING_IF("CreateEntitiesFromGeometriesModeler", n_elements != 0)
        << "There are " << n_elements << " elements in '" << rModelPart.FullName() << "' model part. These are to be removed." << std::endl;
    VariableUtils().SetFlag(TO_ERASE, true, rModelPart.Elements());
    rModelPart.RemoveElementsFromAllLevels(TO_ERASE);
}

template <>
void CreateEntitiesFromGeometriesModeler::RemoveModelPartEntities<Condition>(ModelPart &rModelPart)
{
    const SizeType n_conditions = rModelPart.NumberOfConditions();
    KRATOS_WARNING_IF("CreateEntitiesFromGeometriesModeler", n_conditions != 0)
        << "There are " << n_conditions << " conditions in '" << rModelPart.FullName() << "' model part. These are to be removed." << std::endl;
    VariableUtils().SetFlag(TO_ERASE, true, rModelPart.Conditions());
    rModelPart.RemoveConditionsFromAllLevels(TO_ERASE);
}

template <>
void CreateEntitiesFromGeometriesModeler::LoopEntitiesList<Element>(Parameters EntitiesList)
{
    for (auto& r_data : EntitiesList) {
        // Get model part in which the entities are to be created
        const std::string& r_model_part_name = r_data["model_part_name"].GetString();
        auto& r_model_part = mpModel->GetModelPart(r_model_part_name);

        // Get entities substitution list
        // Note that the entities identifier allow to substitue different geometry types
        const std::string& r_element_name = r_data["element_name"].GetString();
        auto entity_identifier = EntitiesUtilities::EntitityIdentifier<Element>(r_element_name);

        // Wipe current model part entities
        RemoveModelPartEntities<Element>(r_model_part);

        // Create submodelpart elements from geometries
        CreateEntitiesFromGeometries<ModelPart::ElementsContainerType, EntitiesUtilities::EntitityIdentifier<Element>>(entity_identifier, r_model_part);
    }
}

template <>
void CreateEntitiesFromGeometriesModeler::LoopEntitiesList<Condition>(Parameters EntitiesList)
{
    for (auto& r_data : EntitiesList) {
        // Get model part in which the entities are to be created
        const std::string& r_model_part_name = r_data["model_part_name"].GetString();
        auto& r_model_part = mpModel->GetModelPart(r_model_part_name);

        // Get entities substitution list
        // Note that the entities identifier allow to substitue different geometry types
        const std::string& r_element_name = r_data["condition_name"].GetString();
        auto entity_identifier = EntitiesUtilities::EntitityIdentifier<Condition>(r_element_name);

        // Wipe current model part entities
        RemoveModelPartEntities<Condition>(r_model_part);

        // Create submodelpart conditions from geometries
        CreateEntitiesFromGeometries<ModelPart::ConditionsContainerType, EntitiesUtilities::EntitityIdentifier<Condition>>(entity_identifier, r_model_part);
    }
}

void CreateEntitiesFromGeometriesModeler::SetupModelPart()
{
    // Check parallelism
    //TODO: Make this MPI compatible by renumbering the local ids
    KRATOS_ERROR_IF(ParallelEnvironment::GetDefaultDataCommunicator().IsDistributed())
        << "CreateEntitiesFromGeometriesModeler is not MPI compatible yet." << std::endl;

    // Get the elements list from input settings
    const auto& r_elements_list = mParameters.GetValue("elements_list");
    const SizeType n_element_pairs = r_elements_list.size();
    KRATOS_INFO_IF("CreateEntitiesFromGeometriesModeler", mEchoLevel != 0 && n_element_pairs == 0)
        << "No elements found in element list." << std::endl;

    // Loop the element list to create the correspoding entities
    LoopEntitiesList<Element>(r_elements_list);

    // Get the conditions list from input settings
    const auto& r_conditions_list = mParameters.GetValue("conditions_list");
    const SizeType n_condition_pairs = r_conditions_list.size();
    KRATOS_INFO_IF("CreateEntitiesFromGeometriesModeler", mEchoLevel != 0 && n_condition_pairs == 0)
        << "No conditions found in condition list." << std::endl;

    // Loop the condition list to create the correspoding entities
    LoopEntitiesList<Condition>(r_conditions_list);
}

}  // namespace Kratos.


