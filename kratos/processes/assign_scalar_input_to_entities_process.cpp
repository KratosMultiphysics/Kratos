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
#include "containers/model.h"
#include "utilities/string_utilities.h"
#include "processes/assign_scalar_input_to_entities_process.h"

namespace Kratos
{

/// Local Flags
template<class TEntity>
const Kratos::Flags AssignScalarInputToEntitiesProcess<TEntity>::GEOMETRIC_DEFINITION(Kratos::Flags::Create(0));

/***********************************************************************************/
/***********************************************************************************/

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

    const std::string& r_variable_name = rParameters["variable_name"].GetString();
    KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Get(r_variable_name)) << "The variable " << r_variable_name << " does not exist" << std::endl;
    mpVariable = &KratosComponents<Variable<double>>::Get(r_variable_name);

    // Define the working dimension
    mDimension = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];

    // We have two options, or the values are defined in the entities or in a geometry, we don't know until we read the database
    mpDataModelPart = &mrModelPart.GetModel().CreateModelPart("AUXILIAR_MODEL_PART_INPUT_VARIABLE_" + r_variable_name);

    // Read the input file
    const std::string& r_filename = rParameters["file"].GetString();
    if (StringUtilities::SearchPartialString(r_filename, ".txt")) {
        ReadDataTXT(r_filename);
    } else if (StringUtilities::SearchPartialString(r_filename, ".json")) {
        ReadDataJSON(r_filename);
    } else {
        KRATOS_ERROR << "The process is only compatible with JSON and TXT" << std::endl;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void AssignScalarInputToEntitiesProcess<TEntity>::Execute()
{
    KRATOS_TRY;

    // TODO
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

template<class TEntity>
void AssignScalarInputToEntitiesProcess<TEntity>::ReadDataTXT(const std::string& rFileName)
{
    std::ifstream infile(rFileName);
    KRATOS_ERROR_IF_NOT(infile.good()) << "TXT file: " << rFileName << " cannot be found" << std::endl;
    std::stringstream buffer;
    buffer << infile.rdbuf();
    const std::string& r_string_file = buffer.str();

    // TODO
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void AssignScalarInputToEntitiesProcess<TEntity>::ReadDataJSON(const std::string& rFileName)
{
    std::ifstream infile(rFileName);
    KRATOS_ERROR_IF_NOT(infile.good()) << "JSON file: " << rFileName << " cannot be found" << std::endl;
    std::stringstream buffer;
    buffer << infile.rdbuf();
    Parameters json_input(buffer.str());

    // TODO
}

/***********************************************************************************/
/***********************************************************************************/

template class AssignScalarInputToEntitiesProcess<Node<3>>;
template class AssignScalarInputToEntitiesProcess<Condition>;
template class AssignScalarInputToEntitiesProcess<Element>;
template class AssignScalarInputToEntitiesProcess<MasterSlaveConstraint>;

}  // namespace Kratos.
