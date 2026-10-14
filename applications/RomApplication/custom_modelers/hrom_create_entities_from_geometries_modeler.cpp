//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Marco Antonio Zuñiga Perez
//

// System includes

// External includes

// Project includes
#include "includes/parallel_environment.h"
#include "utilities/entities_utilities.h"
#include "utilities/variable_utils.h"
#include "custom_utilities/rom_auxiliary_utilities.h"
#include "hrom_create_entities_from_geometries_modeler.h"

namespace Kratos
{

namespace {

/**
 * @brief Create entities from geometries.
 * @param rEntityIdentifier the entity identifier (@see EntitiesUtilities)
 * @param rModelPart the model part
 * @param rGeomToEntityId map from geometry id to entity id (used when AssignIdsFromFile is true)
 * @param AssignIdsFromFile whether the entity ids are to be taken from the map or generated
 */
template <class TEntitiesContainerType, class TEntityIdentifier>
void CreateEntitiesFromGeometries(
    const TEntityIdentifier& rEntityIdentifier,
    ModelPart& rModelPart,
    const std::unordered_map<std::size_t, std::size_t>& rGeomToEntityId,
    const bool AssignIdsFromFile)
{
    // Create the entities container and allocate space
    std::vector<typename TEntitiesContainerType::pointer> entities_to_add;
    entities_to_add.reserve(rModelPart.NumberOfGeometries());

    // Get current max entity id
    std::size_t max_id = 0;
    const auto& r_root_model_part = rModelPart.GetRootModelPart();
    if constexpr (std::is_same<typename TEntitiesContainerType::value_type, Element>::value) {
        max_id = block_for_each<MaxReduction<std::size_t>>(r_root_model_part.Elements(), [](auto& rElement){
            return rElement.Id();
        });
    } else if constexpr (std::is_same<typename TEntitiesContainerType::value_type, Condition>::value) {
        max_id = block_for_each<MaxReduction<std::size_t>>(r_root_model_part.Conditions(), [](auto& rCondition){
            return rCondition.Id();
        });
    }

    // Loop geometries to create the corresponding entities from them
    // Note that we retrieve the corresponding prototype entity from the entities idenfifier
    // This makes possible to loop and create entities from different type geometries
    for (auto it = rModelPart.GeometriesBegin(); it != rModelPart.GeometriesEnd(); ++it) {
        auto p_geom = *(it.base());
        if (rEntityIdentifier.HasPrototypeEntity(*p_geom)) {
            const auto& r_ref_entity = rEntityIdentifier.GetPrototypeEntity(*p_geom);
            std::size_t new_id;
            if (AssignIdsFromFile) {
                const auto it_id = rGeomToEntityId.find(p_geom->Id());
                if (it_id != rGeomToEntityId.end()) {
                    new_id = it_id->second;
                } else {
                    new_id = ++max_id;
                }
            } else {
                new_id = ++max_id;
            }
            auto p_entity = r_ref_entity.Create(new_id, p_geom, nullptr);
            entities_to_add.push_back(p_entity);
        }
    }

    // Add the created entities to current submodelpart
    if constexpr (std::is_same<typename TEntitiesContainerType::value_type, Element>::value) {
        rModelPart.AddElements(entities_to_add.begin(), entities_to_add.end());
    } else if constexpr (std::is_same<typename TEntitiesContainerType::value_type, Condition>::value) {
        rModelPart.AddConditions(entities_to_add.begin(), entities_to_add.end());
    }

    // Entities created from geometries have no properties, which is required to write the mdpa format
    RomAuxiliaryUtilities::AssignPropertiesToEntitiesIfMissing(rModelPart);
}

}

template <class TEntitiesContainerType>
void HRomCreateEntitiesFromGeometriesModeler::RemoveModelPartEntities(ModelPart &rModelPart)
{
    using EntityType = typename TEntitiesContainerType::value_type;

    // Get number of entities (elements or conditions)
    SizeType n_entities;
    std::string entity_name;
    if constexpr (std::is_same<EntityType, Element>::value) {
        entity_name = "elements";
        n_entities = rModelPart.NumberOfElements();
    } else if constexpr (std::is_same<EntityType, Condition>::value) {
        entity_name = "conditions";
        n_entities = rModelPart.NumberOfConditions();
    }

    // Inform the user of the entities to be removed
    KRATOS_WARNING_IF("HRomCreateEntitiesFromGeometriesModeler", n_entities != 0)
        << "There are " << n_entities << " " << entity_name << " in '" << rModelPart.FullName() << "' model part. These are to be removed." << std::endl;

    // Wipe the corresponding container
    if (n_entities != 0) {
        if constexpr (std::is_same<EntityType, Element>::value) {
            VariableUtils().SetFlag(TO_ERASE, true, rModelPart.Elements());
            rModelPart.RemoveElementsFromAllLevels(TO_ERASE);
        } else if constexpr (std::is_same<EntityType, Condition>::value) {
            VariableUtils().SetFlag(TO_ERASE, true, rModelPart.Conditions());
            rModelPart.RemoveConditionsFromAllLevels(TO_ERASE);
        }
    }
}

template <class TEntitiesContainerType>
void HRomCreateEntitiesFromGeometriesModeler::LoopEntitiesList(Parameters EntitiesList)
{
    // Get entity type from container type
    using EntityType = typename TEntitiesContainerType::value_type;

    // Loop the entities to be created list
    for (auto& r_data : EntitiesList) {
        // Get model part in which the entities are to be created
        const std::string& r_model_part_name = r_data["model_part_name"].GetString();
        auto& r_model_part = mpModel->GetModelPart(r_model_part_name);

        // Get the name of the entity
        std::string entity_name_key;
        if constexpr (std::is_same<EntityType, Element>::value) {
            entity_name_key = "element_name";
        } else if constexpr (std::is_same<EntityType, Condition>::value) {
            entity_name_key = "condition_name";
        }

        // Get entities substitution list
        // Note that the entities identifier allow to substitue different geometry types
        const std::string& r_entity_name = r_data[entity_name_key].GetString();
        auto entity_identifier = EntitiesUtilities::EntitityIdentifier<EntityType>(r_entity_name);

        // Wipe current model part entities
        RemoveModelPartEntities<TEntitiesContainerType>(r_model_part);

        // Create submodelpart elements from geometries
        CreateEntitiesFromGeometries<TEntitiesContainerType, EntitiesUtilities::EntitityIdentifier<EntityType>>(
            entity_identifier,
            r_model_part,
            mGeomToEntityId,
            mParameters["assign_ids_from_file"].GetBool());
    }
}

void HRomCreateEntitiesFromGeometriesModeler::BuildMapToAssignEntityIdsFromFile()
{
    KRATOS_TRY

    const std::string ids_filename = mParameters["geometry_ids_file"].GetString();

    KRATOS_ERROR_IF(ids_filename == "")
        << "geometry_ids_file is empty but assign_ids_from_file is true." << std::endl;

    std::ifstream infile(ids_filename);
    KRATOS_ERROR_IF_NOT(infile.is_open())
        << "Could not open file: " << ids_filename << std::endl;

    mGeomToEntityId.clear();

    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty()) {
            continue;
        }

        std::stringstream ss(line);
        std::string geom_id_str, entity_id_str;
        std::getline(ss, geom_id_str, ',');
        std::getline(ss, entity_id_str, ',');

        const IndexType geom_id = std::stoul(geom_id_str);
        const IndexType entity_id = std::stoul(entity_id_str);

        KRATOS_ERROR_IF(mGeomToEntityId.find(geom_id) != mGeomToEntityId.end())
            << "Duplicate geometry id " << geom_id << " in geometry_ids_file." << std::endl;

        mGeomToEntityId[geom_id] = entity_id;
    }

    infile.close();

    KRATOS_CATCH("")
}

void HRomCreateEntitiesFromGeometriesModeler::CreateMissingSubModelParts(const Parameters& r_entities_list)
{
    for (auto& r_data : r_entities_list) {
        const std::string& r_model_part_name = r_data["model_part_name"].GetString();

        // Split the model part name into its path components
        std::vector<std::string> path;
        std::stringstream ss(r_model_part_name);
        std::string component;
        while (std::getline(ss, component, '.')) {
            path.push_back(component);
        }

        // Attain (and create if needed) the sub model parts along the path
        ModelPart* p_model_part = &(mpModel->GetModelPart(path[0]));
        for (std::size_t i = 1; i < path.size(); ++i) {
            if (p_model_part->HasSubModelPart(path[i])) {
                p_model_part = &(p_model_part->GetSubModelPart(path[i]));
            } else {
                p_model_part = &(p_model_part->CreateSubModelPart(path[i]));
            }
        }
    }
}

void HRomCreateEntitiesFromGeometriesModeler::SetupModelPart()
{
    if (mParameters["assign_ids_from_file"].GetBool()) {
        BuildMapToAssignEntityIdsFromFile();
    }

    // Check parallelism
    //TODO: Make this MPI compatible by renumbering the local ids
    KRATOS_ERROR_IF(ParallelEnvironment::GetDefaultDataCommunicator().IsDistributed())
        << "HRomCreateEntitiesFromGeometriesModeler is not MPI compatible yet." << std::endl;

    // Get the elements list from input settings
    const auto& r_elements_list = mParameters.GetValue("elements_list");
    const SizeType n_element_pairs = r_elements_list.size();
    KRATOS_INFO_IF("HRomCreateEntitiesFromGeometriesModeler", mEchoLevel != 0 && n_element_pairs == 0)
        << "No elements found in element list." << std::endl;

    // Create the referenced sub model parts if they are missing (e.g. empty ones in a reduced mesh)
    CreateMissingSubModelParts(r_elements_list);

    // Loop the element list to create the correspoding entities
    LoopEntitiesList<ElementsContainerType>(r_elements_list);

    // Get the conditions list from input settings
    const auto& r_conditions_list = mParameters.GetValue("conditions_list");
    const SizeType n_condition_pairs = r_conditions_list.size();
    KRATOS_INFO_IF("HRomCreateEntitiesFromGeometriesModeler", mEchoLevel != 0 && n_condition_pairs == 0)
        << "No conditions found in condition list." << std::endl;

    // Create the referenced sub model parts if they are missing (e.g. empty ones in a reduced mesh)
    CreateMissingSubModelParts(r_conditions_list);

    // Loop the condition list to create the correspoding entities
    LoopEntitiesList<ConditionsContainerType>(r_conditions_list);
}

}  // namespace Kratos.
