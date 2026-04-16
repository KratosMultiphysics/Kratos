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
#include "includes/define_registry.h"
#include "includes/kratos_components.h"
#include "modeler/modeler.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/**
 * @class RomCreateEntitiesFromGeometriesModeler
 * @ingroup KratosCore
 * @brief Modeler to create element/condition entities from geometries
 * @details Given a list of pairs of the type submodelpart name and entity (element/condition) name,
 * this modeler creates the corresponding entities from the submodelpart geometries.
 * @author Ruben Zorrilla
 */
class KRATOS_API(KRATOS_CORE) RomCreateEntitiesFromGeometriesModeler
    : public Modeler
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Modeler
    KRATOS_CLASS_POINTER_DEFINITION(RomCreateEntitiesFromGeometriesModeler);

    /// Size type definition
    using SizeType = std::size_t;

    /// Elements container type
    using ElementsContainerType = ModelPart::ElementsContainerType;

    /// Conditions container type
    using ConditionsContainerType = ModelPart::ConditionsContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    RomCreateEntitiesFromGeometriesModeler() : Modeler() {};

    /// Constructor with Model
    RomCreateEntitiesFromGeometriesModeler(
        Model &rModel,
        Parameters Settings)
        : Modeler(rModel, Settings)
    {
        mpModel = &rModel;
        mParameters.ValidateAndAssignDefaults(GetDefaultParameters());
    }

    /// Destructor.
    ~RomCreateEntitiesFromGeometriesModeler() = default;

    /// Creates the RomCreateEntitiesFromGeometriesModeler Pointer
    Modeler::Pointer Create(
        Model& rModel,
        const Parameters Settings) const override
    {
        return Kratos::make_shared<RomCreateEntitiesFromGeometriesModeler>(rModel, Settings);
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
    ///@name Operators
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "RomCreateEntitiesFromGeometriesModeler";
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

    KRATOS_REGISTRY_ADD_PROTOTYPE("Modelers.KratosMultiphysics.RomApplication", Modeler, RomCreateEntitiesFromGeometriesModeler)
    KRATOS_REGISTRY_ADD_PROTOTYPE("Modelers.All", Modeler, RomCreateEntitiesFromGeometriesModeler)

    ///@name Private members
    ///@{

    Model* mpModel = nullptr; /// The model considered in the geometry replacement
    std::unordered_map<IndexType, IndexType> geom_to_entity_id;
    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Loops through a list of entities.
     * @tparam TEntitiesContainerType The entities container type.
     * @param EntitiesList The list of entities.
     */
    template <class TEntitytiesContainerType>
    void LoopEntitiesList(Parameters EntitiesList);

    /**
     * @brief Removes entities from a model part.
     * @tparam TEntitiesContainerType The entities container type.
     * @param rModelPart The model part from which entities will be removed.
     */
    template <class TEntitiesContainerType>
    void RemoveModelPartEntities(ModelPart& rModelPart);

    /**
     * @brief Build map to assign entities ids from a file.
     * @param rModelPart The model part.
     */
    void BuldMapToAssignEntityIdsFromFile();

    ///@}
}; // Class RomCreateEntitiesFromGeometriesModeler

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (
    std::istream& rIStream,
    RomCreateEntitiesFromGeometriesModeler& rThis);

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const RomCreateEntitiesFromGeometriesModeler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.
