//
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
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
#include "processes/assign_scalar_input_to_entities_process.h"

namespace Kratos
{

template<class TEntity>
AssignScalarInputToEntitiesProcess<TEntity>::AssignScalarInputToEntitiesProcess(
    ModelPart& rModelPart,
    Parameters rParameters
    ) : Process(Flags()) ,
        mrModelPart(rModelPart)
{
    KRATOS_TRY

    // Validate against defaults -- this ensures no type mismatch
    const Parameters default_parameters = GetDefaultParameters();
    rParameters.ValidateAndAssignDefaults(default_parameters);

    mpVariable = &KratosComponents<Variable<double>>::Get(rParameters["variable_name"].GetString());

//     // Read json string in materials file, create Parameters
//     const std::string& materials_filename = Params["Parameters"]["materials_filename"].GetString();
//     std::ifstream infile(materials_filename);
//     KRATOS_ERROR_IF_NOT(infile.good()) << "Materials file: " << materials_filename << " cannot be found" << std::endl;
//     std::stringstream buffer;
//     buffer << infile.rdbuf();
//     Parameters json_input(buffer.str());

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void AssignScalarInputToEntitiesProcess<TEntity>::Execute()
{
    KRATOS_TRY;

    InternalAssignValue<>(*mpVariable, 0.0);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
const Parameters AssignScalarInputToEntitiesProcess<TEntity>::GetDefaultParameters() const
{
    const Parameters default_parameters( R"(
    {
        "model_part_name" : "MODEL_PART_NAME",
        "mesh_id"         : 0,
        "variable_name"   : "VARIABLE_NAME",
        "file"            : ""
    }  )" );
    return default_parameters;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Node<3>, IndexedObject>& AssignScalarInputToEntitiesProcess<Node<3>>::GetEntitiesContainer()
{
    return mrModelPart.GetMesh().Nodes();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Condition, IndexedObject>& AssignScalarInputToEntitiesProcess<Condition>::GetEntitiesContainer()
{
    return mrModelPart.GetMesh().Conditions();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Element, IndexedObject>& AssignScalarInputToEntitiesProcess<Element>::GetEntitiesContainer()
{
    return mrModelPart.GetMesh().Elements();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<MasterSlaveConstraint, IndexedObject>& AssignScalarInputToEntitiesProcess<MasterSlaveConstraint>::GetEntitiesContainer()
{
    return mrModelPart.GetMesh().MasterSlaveConstraints();
}

/***********************************************************************************/
/***********************************************************************************/

template class AssignScalarInputToEntitiesProcess<Node<3>>;
template class AssignScalarInputToEntitiesProcess<Condition>;
template class AssignScalarInputToEntitiesProcess<Element>;
template class AssignScalarInputToEntitiesProcess<MasterSlaveConstraint>;

}  // namespace Kratos.
