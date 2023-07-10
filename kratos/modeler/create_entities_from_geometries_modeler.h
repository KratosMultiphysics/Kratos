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
//

#pragma once

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/define.h"
#include "includes/kratos_components.h"
#include "modeler/modeler.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/**
 * @class CreateEntitiesFromGeometriesModeler
 * @ingroup KratosCore
 * @brief Modeler to create element/condition entities from geometries
 * @details Given a list of pairs of the type submodelpart name and entity (element/condition) name,
 * this modeler creates the corresponding entities from the submodelpart geometries.
 * @author Ruben Zorrilla
 */
class KRATOS_API(KRATOS_CORE) CreateEntitiesFromGeometriesModeler 
    : public Modeler
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Modeler
    KRATOS_CLASS_POINTER_DEFINITION(CreateEntitiesFromGeometriesModeler);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Model
    CreateEntitiesFromGeometriesModeler(
        Model& rModel,
        Parameters Settings)
        : Modeler(rModel, Settings)
    {
        mpModel = &rModel;
        mParameters.ValidateAndAssignDefaults(GetDefaultParameters());
    }

    /// Destructor.
    ~CreateEntitiesFromGeometriesModeler() = default;

    /// Creates the CreateEntitiesFromGeometriesModeler Pointer
    Modeler::Pointer Create(
        Model& rModel,
        const Parameters Settings) const override
    {
        return Kratos::make_shared<CreateEntitiesFromGeometriesModeler>(rModel, Settings);
    }

    ///@}
    ///@name Operations
    ///@{

    /// Convert the geometry model or import analysis suitable models.
    void SetupModelPart() override;

    /// This method provides the defaults parameters to avoid conflicts between the different constructors
    const Parameters GetDefaultParameters() const override
    {
        const Parameters default_parameters = Parameters(R"({
            "echo_level" : 0,
            "elements_list" : [],
            "conditions_list" : []
        })");
        return default_parameters;
    }

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "CreateEntitiesFromGeometriesModeler";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
private:
    ///@name Private members
    ///@{

    Model* mpModel = nullptr; /// The model considered in the geometry replacement

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Loops through a list of entities.
     * @tparam TEntityType The type of entity.
     * @param EntitiesList The list of entities.
     */
    template <class TEntityType>
    void LoopEntitiesList(Parameters EntitiesList);

    /**
     * @brief Removes entities from a model part.
     * @tparam TEntityType The type of entity.
     * @param rModelPart The model part from which entities will be removed.
     */
    template <class TEntityType>
    void RemoveModelPartEntities(ModelPart& rModelPart);

    /**
     * @brief Create entities from geometries.
     * @param EntityName the name of the entity
     * @param rModelPart the model part
     */
    template <class TEntityType, class TEntitiesContainerType>
    void CreateEntitiesFromGeometries(
        const std::string EntityName,
        ModelPart& rModelPart
        )
    {
        // Get entity prototype from KratosComponents
        const auto& r_ref_entity = KratosComponents<TEntityType>::Get(EntityName);

        // Create the entities container and allocate space
        TEntitiesContainerType entities_to_add;
        entities_to_add.reserve(rModelPart.NumberOfGeometries());

        // Get current max element id
        SizeType max_id;
        const auto& r_root_model_part = rModelPart.GetRootModelPart();
        if constexpr (std::is_same<TEntityType, Element>::value) {
            max_id = block_for_each<MaxReduction<SizeType>>(r_root_model_part.Elements(), [](auto& rElement){
                return rElement.Id();
            });
        } else {
            max_id = block_for_each<MaxReduction<SizeType>>(r_root_model_part.Conditions(), [](auto& rCondition){
                return rCondition.Id();
            });
        }

        // Loop geometries to create the corresponding entities from them
        for (auto& r_geom : rModelPart.Geometries()) {
            auto p_entity = r_ref_entity.Create(++max_id, r_geom, nullptr);
            entities_to_add.push_back(p_entity);
        }

        // Add the created entities to current submodelpart
        if constexpr (std::is_same<TEntityType, Element>::value) {
            rModelPart.AddElements(entities_to_add.begin(), entities_to_add.end());
        } else {
            rModelPart.AddConditions(entities_to_add.begin(), entities_to_add.end());
        }
    }

    ///@}
}; // Class CreateEntitiesFromGeometriesModeler

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (
    std::istream& rIStream,
    CreateEntitiesFromGeometriesModeler& rThis);

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const CreateEntitiesFromGeometriesModeler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.
