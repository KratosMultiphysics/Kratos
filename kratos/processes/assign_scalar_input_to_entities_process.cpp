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

    // Getting algorithm
    mAlgorithm = ConvertAlgorithmString(rParameters["transfer_algorithm"].GetString());

    // Get the geometry or entities
    const std::string& r_filename = rParameters["file"].GetString();
    if (StringUtilities::ContainsPartialString(r_filename, ".txt")) {
        IdentifyDataTXT(r_filename);
    } else if (StringUtilities::ContainsPartialString(r_filename, ".json")) {
        IdentifyDataJSON(r_filename);
    } else {
        KRATOS_ERROR << "The process is only compatible with JSON and TXT" << std::endl;
    }

    // Read the input file
    if (StringUtilities::ContainsPartialString(r_filename, ".txt")) {
        ReadDataTXT(r_filename);
    } else if (StringUtilities::ContainsPartialString(r_filename, ".json")) {
        ReadDataJSON(r_filename);
    } else {
        KRATOS_ERROR << "The process is only compatible with JSON and TXT" << std::endl;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void AssignScalarInputToEntitiesProcess<TEntity>::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY;

    // Get time
    const double time = mrModelPart.GetProcessInfo().GetValue(TIME);

    // Case of only one entity defined
    const SizeType number_of_entities = mCoordinates.size();
    const auto& r_var_database = mDatabase.GetVariableData(*mpVariable);
    if (number_of_entities == 1) {
        InternalAssignValue<>(*mpVariable, r_var_database.GetValue(0, time));
    } else {
        // TODO
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
const Parameters AssignScalarInputToEntitiesProcess<TEntity>::GetDefaultParameters() const
{
    const Parameters default_parameters( R"(
    {
        "model_part_name"    : "MODEL_PART_NAME",
        "mesh_id"            : 0,
        "variable_name"      : "VARIABLE_NAME",
        "file"               : "",
        "transfer_algorithm" : "nearest_neighbour"
    }  )" );
    return default_parameters;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
const std::string AssignScalarInputToEntitiesProcess<Node<3>>::GetEntitiesLabel()
{
    return "NODE_";
}

/***********************************************************************************/
/***********************************************************************************/

template<>
const std::string AssignScalarInputToEntitiesProcess<Condition>::GetEntitiesLabel()
{
    return "CONDITION_";
}

/***********************************************************************************/
/***********************************************************************************/

template<>
const std::string AssignScalarInputToEntitiesProcess<Element>::GetEntitiesLabel()
{
    return "ELEMENT_";
}

/***********************************************************************************/
/***********************************************************************************/

template<>
const std::string AssignScalarInputToEntitiesProcess<MasterSlaveConstraint>::GetEntitiesLabel()
{
    return "CONSTRAINT";
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
void AssignScalarInputToEntitiesProcess<TEntity>::IdentifyDataTXT(const std::string& rFileName)
{
    // Read txt
    std::ifstream infile(rFileName);
    KRATOS_ERROR_IF_NOT(infile.good()) << "TXT file: " << rFileName << " cannot be found" << std::endl;
    std::stringstream buffer;
    buffer << infile.rdbuf();

    // First line
    std::string line;
    std::getline(buffer, line);

    // Checking if geometric definition or entity identifier
    if (StringUtilities::ContainsPartialString(line, "(") && StringUtilities::ContainsPartialString(line, ")")) {
        this->Set(GEOMETRIC_DEFINITION, true);
    } else {
        this->Set(GEOMETRIC_DEFINITION, false);
    }

    std::istringstream iss(line);
    std::string token;
    SizeType counter = 0;
    if (this->Is(GEOMETRIC_DEFINITION)) {
        while(std::getline(iss, token, '\t')) {
            if (counter > 0) {
//                 std::stod(token, &sz);
            }
            ++counter;
        }
    } else {
        while(std::getline(iss, token, '\t')) {
            if (counter > 0) {
//                 std::stod(token, &sz);
            }
            ++counter;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void AssignScalarInputToEntitiesProcess<TEntity>::IdentifyDataJSON(const std::string& rFileName)
{

}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void AssignScalarInputToEntitiesProcess<TEntity>::ReadDataTXT(const std::string& rFileName)
{
    // Initialize the databases
    std::vector<IndexType> variables_ids(1);
    variables_ids[0] = mpVariable->Key();
    std::vector<IndexType> values_sizes(1, 1);
    const SizeType number_of_entities = mCoordinates.size();
    mDatabase.Initialize(variables_ids, values_sizes, number_of_entities);

    // Read txt
    std::ifstream infile(rFileName);
    KRATOS_ERROR_IF_NOT(infile.good()) << "TXT file: " << rFileName << " cannot be found" << std::endl;
    std::stringstream buffer;
    buffer << infile.rdbuf();

    // First line
    std::string line;
    std::getline(buffer, line);

    // The other lines
    SizeType number_time_steps = 0;
    while(std::getline(buffer, line)) {
        ++number_time_steps;
    }

    Vector time = ZeroVector(number_time_steps);
    std::vector<Vector> values(number_of_entities, time);

    // Reset buffer
    buffer.str("");
    buffer << infile.rdbuf();

    // First line
    std::getline(buffer, line);

    // The other lines
    SizeType counter = 0;
    SizeType sub_counter = 0;
    std::string::size_type sz;     // alias of size_t
    while(std::getline(buffer, line)) {
        std::istringstream iss(line);
        std::string token;
        sub_counter = 0;
        while(std::getline(iss, token, '\t')) {
            if (sub_counter == 0) {
                time[counter] = std::stod(token, &sz);
            } else {
                values[sub_counter - 1][counter] = std::stod(token, &sz);
            }
            ++sub_counter;
        }
        ++counter;
    }

    // Set the time table
    mDatabase.SetCommonColumn(time);

    // Set the entities values
    auto& r_var_database = mDatabase.GetVariableData(*mpVariable);
    for (IndexType i = 0; i < values.size(); ++i) {
        r_var_database.SetValues(time, values[i], i);
    }
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

    // Initialize the databases
    std::vector<IndexType> variables_ids(1);
    variables_ids[0] = mpVariable->Key();
    std::vector<IndexType> values_sizes(1, 1);
    const SizeType number_of_entities = mCoordinates.size();
    mDatabase.Initialize(variables_ids, values_sizes, number_of_entities);

    // Get the time vector
    const Vector& r_time = json_input["TIME"].GetVector();
    mDatabase.SetCommonColumn(r_time);

    // Fill database
    const std::string ent_label = GetEntitiesLabel();
    auto& r_var_database = mDatabase.GetVariableData(*mpVariable);
    const std::string& r_variable_name = mpVariable->Name();
//     if (this->Is(GEOMETRIC_DEFINITION)) {
//         // TODO: UPDATE
//         for (int i = 0; i < static_cast<int>(number_of_entities); ++i) {
//             auto it_ent = it_ent_begin + i;
//             const std::string identifier = ent_label + std::to_string(it_ent->Id());
//             const auto& r_vector = json_input[identifier][r_variable_name].GetVector();
//             r_var_database.SetValues(r_time, r_vector, i);
//         }
//     } else {
//         for (int i = 0; i < static_cast<int>(number_of_entities); ++i) {
//             auto it_ent = it_ent_begin + i;
//             const std::string identifier = ent_label + std::to_string(it_ent->Id());
//             const auto& r_vector = json_input[identifier][r_variable_name].GetVector();
//             r_var_database.SetValues(r_time, r_vector, i);
//         }
//     }
}

/***********************************************************************************/
/***********************************************************************************/

template class AssignScalarInputToEntitiesProcess<Node<3>>;
template class AssignScalarInputToEntitiesProcess<Condition>;
template class AssignScalarInputToEntitiesProcess<Element>;

}  // namespace Kratos.
