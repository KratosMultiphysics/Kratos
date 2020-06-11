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
#include "utilities/variable_utils.h"
#include "processes/assign_scalar_input_to_entities_process.h"

namespace Kratos
{

/// Local Flags
template<class TEntity, bool THistorical>
const Kratos::Flags AssignScalarInputToEntitiesProcess<TEntity, THistorical>::GEOMETRIC_DEFINITION(Kratos::Flags::Create(0));

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity, bool THistorical>
AssignScalarInputToEntitiesProcess<TEntity, THistorical>::AssignScalarInputToEntitiesProcess(
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

    // Compute the extrpolation weights
    ComputeExtrapolationWeight();

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity, bool THistorical>
void AssignScalarInputToEntitiesProcess<TEntity, THistorical>::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY;

    // Get time
    const double time = mrModelPart.GetProcessInfo().GetValue(TIME);

    // Case of only one entity defined
    const SizeType number_of_databases = mCoordinates.size();
    const auto& r_var_database = mDatabase.GetVariableData(*mpVariable);
    if (number_of_databases == 1) {
        InternalAssignValue(*mpVariable, r_var_database.GetValue(0, time));
    } else {
        // Getting entities array
        auto& r_entities_array = GetEntitiesContainer();
        const int number_of_entities = static_cast<int>(r_entities_array.size());

        // Initialize values
        ResetValues();

        if(number_of_entities != 0) {
            const auto it_begin = r_entities_array.begin();

            #pragma omp parallel for
            for(int i = 0; i < number_of_entities; i++) {
                auto it_entity = it_begin + i;

                const auto& r_weights = mWeightExtrapolation[i];
                double& r_value = GetValue(*it_entity, *mpVariable);
                for (auto& r_weight : r_weights) {
                    r_value += r_weight.second * r_var_database.GetValue(r_weight.first, time);
                }
            }
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity, bool THistorical>
const Parameters AssignScalarInputToEntitiesProcess<TEntity, THistorical>::GetDefaultParameters() const
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
array_1d<double, 3> AssignScalarInputToEntitiesProcess<Node<3>, AssignScalarInputToEntitiesProcessSettings::SaveAsNonHistoricalVariable>:: GetCoordinatesEntity(const IndexType Id)
{
    return mrModelPart.pGetNode(Id)->Coordinates();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
array_1d<double, 3> AssignScalarInputToEntitiesProcess<Node<3>, AssignScalarInputToEntitiesProcessSettings::SaveAsHistoricalVariable>:: GetCoordinatesEntity(const IndexType Id)
{
    return mrModelPart.pGetNode(Id)->Coordinates();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
array_1d<double, 3> AssignScalarInputToEntitiesProcess<Condition, AssignScalarInputToEntitiesProcessSettings::SaveAsNonHistoricalVariable>:: GetCoordinatesEntity(const IndexType Id)
{
    return mrModelPart.pGetCondition(Id)->GetGeometry().Center().Coordinates();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
array_1d<double, 3> AssignScalarInputToEntitiesProcess<Element, AssignScalarInputToEntitiesProcessSettings::SaveAsNonHistoricalVariable>:: GetCoordinatesEntity(const IndexType Id)
{
    return mrModelPart.pGetElement(Id)->GetGeometry().Center().Coordinates();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Node<3>, IndexedObject>& AssignScalarInputToEntitiesProcess<Node<3>, AssignScalarInputToEntitiesProcessSettings::SaveAsNonHistoricalVariable>::GetEntitiesContainer()
{
    return mrModelPart.GetMesh().Nodes();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Node<3>, IndexedObject>& AssignScalarInputToEntitiesProcess<Node<3>, AssignScalarInputToEntitiesProcessSettings::SaveAsHistoricalVariable>::GetEntitiesContainer()
{
    return mrModelPart.GetMesh().Nodes();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Condition, IndexedObject>& AssignScalarInputToEntitiesProcess<Condition, AssignScalarInputToEntitiesProcessSettings::SaveAsNonHistoricalVariable>::GetEntitiesContainer()
{
    return mrModelPart.GetMesh().Conditions();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Element, IndexedObject>& AssignScalarInputToEntitiesProcess<Element, AssignScalarInputToEntitiesProcessSettings::SaveAsNonHistoricalVariable>::GetEntitiesContainer()
{
    return mrModelPart.GetMesh().Elements();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarInputToEntitiesProcess<Node<3>, AssignScalarInputToEntitiesProcessSettings::SaveAsNonHistoricalVariable>::ResetValues()
{
    VariableUtils().SetNonHistoricalVariable(*mpVariable, 0.0, GetEntitiesContainer());
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarInputToEntitiesProcess<Node<3>, AssignScalarInputToEntitiesProcessSettings::SaveAsHistoricalVariable>::ResetValues()
{
    VariableUtils().SetVariable(*mpVariable, 0.0, GetEntitiesContainer());
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarInputToEntitiesProcess<Condition, AssignScalarInputToEntitiesProcessSettings::SaveAsNonHistoricalVariable>::ResetValues()
{
    VariableUtils().SetNonHistoricalVariable(*mpVariable, 0.0, GetEntitiesContainer());
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarInputToEntitiesProcess<Element, AssignScalarInputToEntitiesProcessSettings::SaveAsNonHistoricalVariable>::ResetValues()
{
    VariableUtils().SetNonHistoricalVariable(*mpVariable, 0.0, GetEntitiesContainer());
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarInputToEntitiesProcess<Node<3>, AssignScalarInputToEntitiesProcessSettings::SaveAsNonHistoricalVariable>::SetValue(
    Node<3>& rEntity,
    const Variable<double>& rVariable,
    const double Value
    )
{
    rEntity.SetValue(rVariable, Value);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarInputToEntitiesProcess<Node<3>, AssignScalarInputToEntitiesProcessSettings::SaveAsHistoricalVariable>::SetValue(
    Node<3>& rEntity,
    const Variable<double>& rVariable,
    const double Value
    )
{
    rEntity.FastGetSolutionStepValue(rVariable) = Value;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarInputToEntitiesProcess<Condition, AssignScalarInputToEntitiesProcessSettings::SaveAsNonHistoricalVariable>::SetValue(
    Condition& rEntity,
    const Variable<double>& rVariable,
    const double Value
    )
{
    rEntity.SetValue(rVariable, Value);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarInputToEntitiesProcess<Element, AssignScalarInputToEntitiesProcessSettings::SaveAsNonHistoricalVariable>::SetValue(
    Element& rEntity,
    const Variable<double>& rVariable,
    const double Value
    )
{
    rEntity.SetValue(rVariable, Value);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
double& AssignScalarInputToEntitiesProcess<Node<3>, AssignScalarInputToEntitiesProcessSettings::SaveAsNonHistoricalVariable>::GetValue(
    Node<3>& rEntity,
    const Variable<double>& rVariable
    )
{
    return rEntity.GetValue(rVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
double& AssignScalarInputToEntitiesProcess<Node<3>, AssignScalarInputToEntitiesProcessSettings::SaveAsHistoricalVariable>::GetValue(
    Node<3>& rEntity,
    const Variable<double>& rVariable
    )
{
    return rEntity.FastGetSolutionStepValue(rVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
double& AssignScalarInputToEntitiesProcess<Condition, AssignScalarInputToEntitiesProcessSettings::SaveAsNonHistoricalVariable>::GetValue(
    Condition& rEntity,
    const Variable<double>& rVariable
    )
{
    return rEntity.GetValue(rVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
double& AssignScalarInputToEntitiesProcess<Element, AssignScalarInputToEntitiesProcessSettings::SaveAsNonHistoricalVariable>::GetValue(
    Element& rEntity,
    const Variable<double>& rVariable
    )
{
    return rEntity.GetValue(rVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity, bool THistorical>
void AssignScalarInputToEntitiesProcess<TEntity, THistorical>::IdentifyDataTXT(const std::string& rFileName)
{
    KRATOS_TRY;

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
    std::string::size_type sz;     // alias of size_t
    if (this->Is(GEOMETRIC_DEFINITION)) {
        array_1d<double, 3> aux_array;
        while(std::getline(iss, token, '\t')) {
            if (counter > 0) {
                std::string aux_string = StringUtilities::ErasePartialString(token, "(");
                aux_string = StringUtilities::ErasePartialString(aux_string, ")");
                std::stringstream s_stream(aux_string); // Create string stream from the string
                SizeType sub_counter = 0;
                std::string substr;
                while(s_stream.good()) {
                    std::getline(s_stream, substr, ','); // Get first string delimited by comma
                    aux_array[sub_counter] = std::stod(substr, &sz);
                    ++sub_counter;
                }
                mCoordinates.push_back(aux_array);
            }
            ++counter;
        }
    } else {
        while(std::getline(iss, token, '\t')) {
            if (counter > 0) {
                const IndexType id = std::stod(token, &sz);
                mCoordinates.push_back(GetCoordinatesEntity(id));
            }
            ++counter;
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity, bool THistorical>
void AssignScalarInputToEntitiesProcess<TEntity, THistorical>::IdentifyDataJSON(const std::string& rFileName)
{
    KRATOS_TRY;

    // Reading json file
    std::ifstream infile(rFileName);
    KRATOS_ERROR_IF_NOT(infile.good()) << "JSON file: " << rFileName << " cannot be found" << std::endl;
    std::stringstream buffer;
    buffer << infile.rdbuf();
    Parameters json_input(buffer.str());

    // Getting number of definitions
    SizeType number_of_definitions = 0;
    for (auto& r_param : json_input) {
        if (!r_param.IsVector()) {  // Removing TIME
            ++number_of_definitions;
        }
    }

    // Check number of definitions
    KRATOS_ERROR_IF(number_of_definitions == 0) << "Definitions must be superior to 0" << std::endl;

    // Reserve
    if (mCoordinates.size() != number_of_definitions) {
        mCoordinates.resize(number_of_definitions);
    }

    KRATOS_ERROR_IF_NOT(json_input.Has("1")) << "Input not properly defined. Input must have values defined ordered" << std::endl;
    if (json_input["1"].Has("ID")) {
        this->Set(GEOMETRIC_DEFINITION, false);
    } else if (json_input["1"].Has("COORDINATES")) {
        this->Set(GEOMETRIC_DEFINITION, true);
    } else {
        KRATOS_ERROR << "ID or COORDINATES must be defined" << std::endl;
    }

    // Iterate over parameters
    for (IndexType i = 0; i < number_of_definitions; ++i) {
        const std::string identifier = std::to_string(i + 1);
        if (this->Is(GEOMETRIC_DEFINITION)) {
            mCoordinates[i] =  array_1d<double, 3>(json_input[identifier]["COORDINATES"].GetVector());
        } else {
            const IndexType id = static_cast<IndexType>(json_input[identifier]["ID"].GetInt());
            mCoordinates[i] = GetCoordinatesEntity(id);
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity, bool THistorical>
void AssignScalarInputToEntitiesProcess<TEntity, THistorical>::ReadDataTXT(const std::string& rFileName)
{
    KRATOS_TRY;

    // Initialize the databases
    std::vector<IndexType> variables_ids(1);
    variables_ids[0] = mpVariable->Key();
    std::vector<IndexType> values_sizes(1, 1);
    const SizeType number_of_definitions = mCoordinates.size();
    mDatabase.Initialize(variables_ids, values_sizes, number_of_definitions);

    // Define the number of time steps
    SizeType number_time_steps = 0;

    // Definition of auxiliar line
    std::string line;

    // Initial read
    {
        // Read txt
        std::ifstream infile(rFileName);
        KRATOS_ERROR_IF_NOT(infile.good()) << "TXT file: " << rFileName << " cannot be found" << std::endl;

        std::stringstream buffer;
        buffer << infile.rdbuf();

        // First line
        std::getline(buffer, line);

        // The other lines
        while(std::getline(buffer, line)) {
            ++number_time_steps;
        }
    }

    Vector time = ZeroVector(number_time_steps);
    std::vector<Vector> values(number_of_definitions, time);

    // Second read txt
    std::ifstream infile(rFileName);
    KRATOS_ERROR_IF_NOT(infile.good()) << "TXT file: " << rFileName << " cannot be found" << std::endl;
    std::stringstream buffer;
    buffer << infile.rdbuf();

    // First line
    std::getline(buffer, line);

    // The other lines
    SizeType counter = 0;
    std::string::size_type sz;     // alias of size_t
    while(std::getline(buffer, line)) {
        std::istringstream iss(line);
        std::string token;
        SizeType sub_counter = 0;

        while(std::getline(iss, token, '\t')) {
            const double value = std::stod(token, &sz);
            if (sub_counter == 0) {
                time[counter] = value;
            } else {
                values[sub_counter - 1][counter] = value;
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

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity, bool THistorical>
void AssignScalarInputToEntitiesProcess<TEntity, THistorical>::ReadDataJSON(const std::string& rFileName)
{
    KRATOS_TRY;

    // Reading json file
    std::ifstream infile(rFileName);
    KRATOS_ERROR_IF_NOT(infile.good()) << "JSON file: " << rFileName << " cannot be found" << std::endl;
    std::stringstream buffer;
    buffer << infile.rdbuf();
    Parameters json_input(buffer.str());

    // Initialize the databases
    std::vector<IndexType> variables_ids(1);
    variables_ids[0] = mpVariable->Key();
    std::vector<IndexType> values_sizes(1, 1);
    const SizeType number_of_definitions = mCoordinates.size();
    mDatabase.Initialize(variables_ids, values_sizes, number_of_definitions);

    // Get the time vector
    const Vector& r_time = json_input["TIME"].GetVector();
    mDatabase.SetCommonColumn(r_time);

    // Fill database
    auto& r_var_database = mDatabase.GetVariableData(*mpVariable);
    const std::string& r_variable_name = mpVariable->Name();
    for (IndexType i = 0; i < number_of_definitions; ++i) {
        const std::string identifier = std::to_string(i + 1);
        const auto& r_vector = json_input[identifier]["VALUES"][r_variable_name].GetVector();
        r_var_database.SetValues(r_time, r_vector, i);
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity, bool THistorical>
void AssignScalarInputToEntitiesProcess<TEntity, THistorical>::ComputeExtrapolationWeight()
{
    KRATOS_TRY;

    // Some definitions
    const auto& r_entities_array = GetEntitiesContainer();
    const auto it_ent_begin = r_entities_array.begin();
    const SizeType number_of_entities = r_entities_array.size();

    // Resize the weight extrapolation vector
    if (mWeightExtrapolation.size() != number_of_entities) {
        mWeightExtrapolation.resize(number_of_entities);
    }

    // Considering different algorithms to fill the weights
    const SizeType number_of_definitions = mCoordinates.size();
    if (mAlgorithm == Algorithm::NEAREST_NEIGHBOUR) {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(number_of_entities); ++i) {
            auto it_ent = it_ent_begin + i;
            const IndexType id = it_ent->Id();
            const array_1d<double, 3> coordinates = GetCoordinatesEntity(id);
            double distance = 1.0e24;
            IndexType index = 0;
            for (IndexType i = 0; i < number_of_definitions; ++i) {
                const double aux_distance = norm_2(coordinates - mCoordinates[i]);
                if (aux_distance < distance) {
                    distance = aux_distance;
                    index = i;
                }
            }
            std::unordered_map<IndexType, double> aux_map({{index, 1.0}});
            mWeightExtrapolation[i] = aux_map;
        }
    } else {
        KRATOS_ERROR << "Algorithm not defined" << std::endl;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity, bool THistorical>
void AssignScalarInputToEntitiesProcess<TEntity, THistorical>::InternalAssignValue(
    const Variable<double>& rVariable,
    const double Value
    )
{
    KRATOS_TRY;

    auto& r_entities_array = GetEntitiesContainer();
    const int number_of_entities = static_cast<int>(r_entities_array.size());

    if(number_of_entities != 0) {
        const auto it_begin = r_entities_array.begin();

        #pragma omp parallel for
        for(int i = 0; i<number_of_entities; i++) {
            auto it_entity = it_begin + i;
            SetValue(*it_entity, rVariable, Value);
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template class AssignScalarInputToEntitiesProcess<Node<3>, AssignScalarInputToEntitiesProcessSettings::SaveAsNonHistoricalVariable>;
template class AssignScalarInputToEntitiesProcess<Node<3>, AssignScalarInputToEntitiesProcessSettings::SaveAsHistoricalVariable>;
template class AssignScalarInputToEntitiesProcess<Condition>;
template class AssignScalarInputToEntitiesProcess<Element>;

}  // namespace Kratos.
