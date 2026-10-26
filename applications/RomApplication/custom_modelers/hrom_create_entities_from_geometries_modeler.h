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
//

#pragma once

// System includes
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/define.h"
#include "includes/define_registry.h"
#include "includes/kratos_components.h"
#include "modeler/create_entities_from_geometries_modeler.h"
#include "modeler/modeler.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/**
 * @class HRomCreateEntitiesFromGeometriesModeler
 * @ingroup RomApplication
 * @brief Modeler to create element/condition entities from geometries preserving a given numbering
 * @details Given a list of pairs of the type submodelpart name and entity (element/condition) name,
 * this modeler creates the corresponding entities from the submodelpart geometries.
 * This is an extension of the CreateEntitiesFromGeometriesModeler that adds the possibility of
 * assigning the entities ids from a file mapping the geometry ids to the entity ids.
 * This is required to recover the original entity numbering when the model part is imported from a
 * geometry based format (e.g. MED) in which the entities ids are not stored.
 * @author Marco Antonio Zuñiga Perez
 */
class KRATOS_API(ROM_APPLICATION) HRomCreateEntitiesFromGeometriesModeler
    : public CreateEntitiesFromGeometriesModeler
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of HRomCreateEntitiesFromGeometriesModeler
    KRATOS_CLASS_POINTER_DEFINITION(HRomCreateEntitiesFromGeometriesModeler);

    /// Size type definition
    using SizeType = std::size_t;

    /// Index type definition
    using IndexType = std::size_t;

    /// Elements container type
    using ElementsContainerType = ModelPart::ElementsContainerType;

    /// Conditions container type
    using ConditionsContainerType = ModelPart::ConditionsContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    HRomCreateEntitiesFromGeometriesModeler() : CreateEntitiesFromGeometriesModeler() {};

    /// Constructor with Model
    HRomCreateEntitiesFromGeometriesModeler(
        Model &rModel,
        Parameters Settings)
        : CreateEntitiesFromGeometriesModeler()
    {
        mpModel = &rModel;
        mParameters = Settings;
        mEchoLevel = Settings.Has("echo_level") ? Settings["echo_level"].GetInt() : 0;
        mParameters.ValidateAndAssignDefaults(GetDefaultParameters());
    }

    /// Destructor.
    ~HRomCreateEntitiesFromGeometriesModeler() = default;

    /// Creates the HRomCreateEntitiesFromGeometriesModeler Pointer
    Modeler::Pointer Create(
        Model& rModel,
        const Parameters Settings) const override
    {
        return Kratos::make_shared<HRomCreateEntitiesFromGeometriesModeler>(rModel, Settings);
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
            "conditions_list" : [],
            "assign_ids_from_file" : false,
            "geometry_ids_file" : ""
        })");
        return default_parameters;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "HRomCreateEntitiesFromGeometriesModeler";
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
    ///@name Static Member Variables
    ///@{

    KRATOS_REGISTRY_ADD_PROTOTYPE("Modelers.KratosMultiphysics.RomApplication", Modeler, HRomCreateEntitiesFromGeometriesModeler)
    KRATOS_REGISTRY_ADD_PROTOTYPE("Modelers.All", Modeler, HRomCreateEntitiesFromGeometriesModeler)

    ///@}
    ///@name Private members
    ///@{

    Model* mpModel = nullptr; /// The model considered in the geometry replacement

    std::unordered_map<IndexType, IndexType> mGeomToEntityId; /// Map from geometry id to entity id (used when assign_ids_from_file is true)

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Loops through a list of entities.
     * @tparam TEntitiesContainerType The entities container type.
     * @param EntitiesList The list of entities.
     */
    template <class TEntitiesContainerType>
    void LoopEntitiesList(Parameters EntitiesList);

    /**
     * @brief Removes entities from a model part.
     * @tparam TEntitiesContainerType The entities container type.
     * @param rModelPart The model part from which entities will be removed.
     */
    template <class TEntitiesContainerType>
    void RemoveModelPartEntities(ModelPart& rModelPart);

    /**
     * @brief Builds the map from geometry id to entity id from the file in geometry_ids_file.
     * @details The file is expected to contain one pair per line with the geometry id and the
     * entity id separated by a comma.
     */
    void BuildMapToAssignEntityIdsFromFile();

    /**
     * @brief Creates the sub model parts referenced in a list of entities if they do not exist.
     * @details In a reduced mesh some of the sub model parts referenced in the elements or
     * conditions list might be empty and, thus, not present in the imported file. This method
     * creates them empty before creating the entities so the base class logic can run.
     * @param r_entities_list The list of entities.
     */
    void CreateMissingSubModelParts(const Parameters& r_entities_list);

    ///@}
}; // Class HRomCreateEntitiesFromGeometriesModeler

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (
    std::istream& rIStream,
    HRomCreateEntitiesFromGeometriesModeler& rThis);

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const HRomCreateEntitiesFromGeometriesModeler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.
