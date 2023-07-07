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
#include "utilities/variable_utils.h"
#include "create_entities_from_geometries_modeler.h"


namespace Kratos
{

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
            // Get current substitution settings
            const std::string entity_name = r_data["element_name"].GetString();
            const std::string model_part_name = r_data["model_part_name"].GetString();
            auto &r_model_part = mpModel->GetModelPart(model_part_name);

            // Wipe current model part entities
            RemoveModelPartEntities<Element>(r_model_part);

            // Create submodelpart elements from geometries
            CreateEntitiesFromGeometries<Element, ModelPart::ElementsContainerType>(entity_name, r_model_part);
        }
    }

    template <>
    void CreateEntitiesFromGeometriesModeler::LoopEntitiesList<Condition>(Parameters EntitiesList)
    {
        for (auto& r_data : EntitiesList) {
            // Get current substitution settings
            const std::string entity_name = r_data["condition_name"].GetString();
            const std::string model_part_name = r_data["model_part_name"].GetString();
            auto &r_model_part = mpModel->GetModelPart(model_part_name);

            // Wipe current model part entities
            RemoveModelPartEntities<Condition>(r_model_part);

            // Create submodelpart conditions from geometries
            CreateEntitiesFromGeometries<Condition, ModelPart::ConditionsContainerType>(entity_name, r_model_part);
        }
    }

    void CreateEntitiesFromGeometriesModeler::SetupModelPart()
    {
        // Check parallelism
        //TODO: Make this MPI compatible by renumbering the local ids
        KRATOS_ERROR_IF(ParallelEnvironment::GetDefaultDataCommunicator().IsDistributed())
            << "CreateEntitiesFromGeometriesModeler is not MPI compatible yet." << std::endl;

        // Get the elements list from input settings
        const auto &r_elements_list = mParameters.GetValue("elements_list");
        const SizeType n_element_pairs = r_elements_list.size();
        KRATOS_WARNING_IF("CreateEntitiesFromGeometriesModeler", mEchoLevel != 0 && n_element_pairs == 0)
            << "No elements found in element list." << std::endl;

        // Loop the element list to create the correspoding entities
        LoopEntitiesList<Element>(r_elements_list);

        // Get the conditions list from input settings
        const auto &r_conditions_list = mParameters.GetValue("conditions_list");
        const SizeType n_condition_pairs = r_conditions_list.size();
        KRATOS_WARNING_IF("CreateEntitiesFromGeometriesModeler", mEchoLevel != 0 && n_condition_pairs == 0)
            << "No conditions found in condition list." << std::endl;

        // Loop the condition list to create the correspoding entities
        LoopEntitiesList<Condition>(r_conditions_list);
    }

}  // namespace Kratos.


