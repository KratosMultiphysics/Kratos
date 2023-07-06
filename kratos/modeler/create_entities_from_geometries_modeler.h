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
#include "utilities/variable_utils.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/**
 * @brief Modeler to create element/condition entities from geometries
 * Given a list of pairs of the type submodelpart name and entity (element/condition) name,
 * this modeler creates the corresponding entities from the submodelpart geometries.
 */
class CreateEntitiesFromGeometriesModeler : public Modeler
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Modeler
    KRATOS_CLASS_POINTER_DEFINITION(CreateEntitiesFromGeometriesModeler);

    // typedef std::size_t SizeType;
    // typedef std::size_t IndexType;

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
    void SetupModelPart() override
    {
        // Get the elements list from input settings
        const auto& r_elements_list = mParameters.GetValue("elements_list");
        const SizeType n_element_pairs = r_elements_list.size();
        KRATOS_WARNING_IF("CreateEntitiesFromGeometriesModeler", mEchoLevel != 0 && n_element_pairs == 0)
            << "No elements found in element list." << std::endl;

        // Loop the element list to create the correspoding entities
        for (auto& r_data : r_elements_list) {
            // Get current substitution settings
            const std::string entity_name = r_data["element_name"].GetString();
            const std::string model_part_name = r_data["model_part_name"].GetString();
            auto& r_model_part = mpModel->GetModelPart(model_part_name);

            // Wipe current model part entities
            RemoveModelPartEntities(r_model_part);

            // Loop submodelpart geometries to create the corresponding entities from them
            const SizeType n_geom = r_model_part.NumberOfGeometries();
            const auto& r_ref_entity = KratosComponents<Element>::Get(entity_name);

            ModelPart::ElementsContainerType entities_to_add_ids;
            entities_to_add_ids.reserve(n_geom);

            // auto elements_to_add = IndexPartition<IndexType>(n_geom).for_each<AccumReduction<Element::Pointer, ModelPart::ElementsContainerType>>(
            //     IndexType iGeom, [](int &rValue) {
            //         auto it_geom = r_model_part.GeometriesBegin() + iGeom;
            //         r_ref_entity.Create(it_geom->Id(), it_geom->pGetGeometry(), nullptr);
            // });

            // IndexPartition<IndexType>(n_geom).for_each([&](IndexType iGeom){
            //     auto it_geom = r_model_part.GeometriesBegin() + iGeom;
            //     // r_model_part.CreateNewElement(entity_name, it_geom->Id(), it_geom->pGetGeometry(), nullptr);
            //     // r_model_part.CreateNewElement(entity_name, it_geom->Id(), it_geom->pGetGeometry(), nullptr);
            // });

        }
    }

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

    Model* mpModel = nullptr;

    ///@}
    ///@name Private Operations
    ///@{

    void RemoveModelPartEntities(ModelPart& rModelPart)
    {
        auto& r_root_model_part = rModelPart.GetRootModelPart();
        const SizeType n_elements = r_root_model_part.NumberOfElements();
        const SizeType n_conditions = r_root_model_part.NumberOfElements();
        KRATOS_WARNING_IF("CreateEntitiesFromGeometriesModeler", n_elements != 0)
            << "There are " << n_elements << " elements in '" << r_root_model_part.FullName() << "' model part. These are to be removed." << std::endl;
        KRATOS_WARNING_IF("CreateEntitiesFromGeometriesModeler", n_conditions != 0)
            << "There are " << n_conditions << " conditions in '" << r_root_model_part.FullName() << "' model part. These are to be removed." << std::endl;
        VariableUtils().SetFlag(TO_ERASE, true, r_root_model_part.Elements());
        VariableUtils().SetFlag(TO_ERASE, true, r_root_model_part.Conditions());
        r_root_model_part.RemoveElementsFromAllLevels(TO_ERASE);
        r_root_model_part.RemoveConditionsFromAllLevels(TO_ERASE);
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
