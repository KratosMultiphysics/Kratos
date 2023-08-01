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
