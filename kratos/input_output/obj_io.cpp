//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "input_output/obj_io.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

namespace Kratos
{

ObjIO::ObjIO(const std::filesystem::path& rFilename, Parameters ThisParameters)
    : mParameters(ThisParameters)
{
    mParameters.ValidateAndAssignDefaults(GetDefaultParameters());
    Kratos::shared_ptr<std::fstream> p_file = Kratos::make_shared<std::fstream>();

    // Check if the file extension is .obj, if not add it
    std::filesystem::path filePath = rFilename;
    if(filePath.extension() != ".obj") {
        filePath += ".obj";
    }

    // set default mode to read
    std::fstream::openmode open_mode;

    const std::string open_mode_str = mParameters["open_mode"].GetString();
    // handle other mode options
    if (open_mode_str == "read") {
        open_mode = std::fstream::in;
    } else if (open_mode_str == "append") {
        open_mode = std::fstream::in | std::fstream::app;
    } else if (open_mode_str == "write") {
        open_mode = std::fstream::out;
    } else {
        KRATOS_ERROR << "Unsupported open mode: " << open_mode_str << std::endl;
    }

    p_file->open(filePath.c_str(), open_mode);

    // Checking read/write status
    if (open_mode_str == "write") {
        KRATOS_ERROR_IF_NOT(p_file->is_open()) << "Could not create the output file  : " << filePath << std::endl;
    } else if (open_mode_str == "append") {
        KRATOS_ERROR_IF_NOT(p_file->is_open()) << "Could not open the output file  : " << filePath << std::endl;
    } else {
        KRATOS_ERROR_IF_NOT(p_file->is_open()) << "Could not open the input file  : " << filePath << std::endl;
    }

    // Store the pointer as a regular std::iostream
    mpInputStream = p_file;
}

/***********************************************************************************/
/***********************************************************************************/

ObjIO::ObjIO(
    Kratos::shared_ptr<std::iostream> pInputStream,
    Parameters ThisParameters
    )
    : IO(),
      mParameters(ThisParameters),
      mpInputStream(pInputStream)
{
    mParameters.ValidateAndAssignDefaults(GetDefaultParameters());
}

/***********************************************************************************/
/***********************************************************************************/

Parameters ObjIO::GetDefaultParameters()
{
    return Parameters(R"({
        "open_mode"            : "read",
        "new_entity_type"      : "element",
        "normal_as_historical" : false
    })" );
}

/***********************************************************************************/
/***********************************************************************************/

void ObjIO::ReadModelPart(ModelPart& rThisModelPart)
{
    // Check if the model part has properties
    if(!rThisModelPart.RecursivelyHasProperties(0)) {
        rThisModelPart.CreateNewProperties(0);
    }

    // Get the next IDs for nodes, elements, and conditions
    mNextNodeId = block_for_each<MaxReduction<std::size_t>>(
        rThisModelPart.GetRootModelPart().Nodes(),
        [](NodeType& rNode) { return rNode.Id();}) + 1;
    mFirstNodeId = mNextNodeId;

    mNextElementId = block_for_each<MaxReduction<std::size_t>>(
        rThisModelPart.GetRootModelPart().Elements(),
        [](Element& rElement) { return rElement.Id();}) + 1;

    mNextConditionId = block_for_each<MaxReduction<std::size_t>>(
        rThisModelPart.GetRootModelPart().Conditions(),
        [](Condition& rCondition) { return rCondition.Id();}) + 1;

    // TODO: Read header to be able to deduce number of vertices and faces

    // Read vertices, normals, and faces
    const bool normal_as_historical = mParameters["normal_as_historical"].GetBool();
    ReadVerticesAndFaces(rThisModelPart, normal_as_historical);
}

/***********************************************************************************/
/***********************************************************************************/

void ObjIO::WriteModelPart(const ModelPart& rThisModelPart)
{
    // To know if we are retrieving normals as historical variables
    const bool normal_as_historical = mParameters["normal_as_historical"].GetBool();
    KRATOS_ERROR << "WriteModelPart is not implemented for ObjIO." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void ObjIO::ReadVerticesAndFaces(
    ModelPart& rThisModelPart,
    const bool NormalAsHistoricalVariable
    )
{
    std::string line;
    while (std::getline(*mpInputStream, line)) {
        // Trim leading and trailing whitespace
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);

        // Ignore empty lines and comments
        if (line.empty() || line[0] == '#')
            continue;

        // Check line type
        if (line.substr(0, 2) == "v ") {
            ParseVertexLine(rThisModelPart, line);
        } else if (line.substr(0, 3) == "vn ") {
            ParseNormalLine(rThisModelPart, line, NormalAsHistoricalVariable);
        } else if (line.substr(0, 2) == "f ") {
            ParseFaceLine(rThisModelPart, line);
        }
        // Other types can be handled if needed
    }
}

/***********************************************************************************/
/***********************************************************************************/

void ObjIO::ParseVertexLine(
    ModelPart& rThisModelPart,
    const std::string& rLine
    )
{
    // The line is expected to be: v x y z [w]
    const std::vector<std::string> tokens = Tokenize(rLine);

    KRATOS_ERROR_IF(tokens.size() < 4) << "Invalid vertex line: " << rLine << std::endl;

    const double x = std::stod(tokens[1]);
    const double y = std::stod(tokens[2]);
    const double z = std::stod(tokens[3]);

    // Create new node
    rThisModelPart.CreateNewNode(mNextNodeId++, x, y, z);
}

/***********************************************************************************/
/***********************************************************************************/

void ObjIO::ParseNormalLine(
    ModelPart& rThisModelPart,
    const std::string& rLine,
    const bool NormalAsHistoricalVariable
    )
{
    // The line is expected to be: vn x y z
    const std::vector<std::string> tokens = Tokenize(rLine);

    KRATOS_ERROR_IF(tokens.size() < 4) << "Invalid vertex normal line: " << rLine << std::endl;

    // Parse the normal
    array_1d<double, 3> normal;
    normal[0] = std::stod(tokens[1]);
    normal[1] = std::stod(tokens[2]);
    normal[2] = std::stod(tokens[3]);

    // Set the normal
    const IndexType node_id = mNormalCounter + mFirstNodeId + 1;
    auto p_node = rThisModelPart.pGetNode(node_id); // TODO: Think about using iterators instead
    if (NormalAsHistoricalVariable) {
       noalias(p_node->FastGetSolutionStepValue(NORMAL)) = normal;
    } else {
        p_node->SetValue(NORMAL, normal);
    }
    mNormalCounter++;
}

/***********************************************************************************/
/***********************************************************************************/

void ObjIO::ParseFaceLine(
    ModelPart& rThisModelPart,
    const std::string& rLine
    )
{
    // The line is expected to be: f v1 v2 v3 ...
    const std::vector<std::string> tokens = Tokenize(rLine);

    // Create a list of nodes for the face
    const std::size_t number_of_tokens = tokens.size();
    KRATOS_ERROR_IF(number_of_tokens < 4) << "Invalid face line: " << rLine << std::endl;

    // Create a list of nodes for the face
    std::vector<IndexType> nodes_ids;
    nodes_ids.reserve(number_of_tokens - 1);

    // For each vertex index in the face
    for (std::size_t i = 1; i < number_of_tokens; ++i) {
        std::string vertex_token = tokens[i];
        // Obj index starts at 1, so we need to subtract 1
        const int vertex_index = std::stoi(vertex_token) + mFirstNodeId - 1;
        nodes_ids.push_back(vertex_index);
    }

    // Create element or condition or geometry
    const std::string new_entity_type = mParameters["new_entity_type"].GetString();
    if (new_entity_type == "geometry") {
        KRATOS_ERROR_IF(number_of_tokens > 5) << "Only support for triangles and quads in geometry creation. Number of nodes: " << number_of_tokens - 1 << std::endl;
        rThisModelPart.CreateNewGeometry(SupportedGeometries[number_of_tokens - 4], nodes_ids);
    } else if (new_entity_type == "element") {
        rThisModelPart.CreateNewElement("Element3D" + std::to_string(number_of_tokens - 1) + "N", mNextElementId++, nodes_ids, rThisModelPart.pGetProperties(0));
    } else if (new_entity_type == "condition") {
        rThisModelPart.CreateNewCondition("SurfaceCondition3D" + std::to_string(number_of_tokens - 1) + "N", mNextConditionId++, nodes_ids, rThisModelPart.pGetProperties(0));
    } else  {
        KRATOS_ERROR << "Invalid new entity type " << new_entity_type << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

std::vector<std::string> ObjIO::Tokenize(const std::string& rLine)
{
    std::vector<std::string> tokens;
    std::string::size_type start = 0;
    std::string::size_type end = 0;

    while ((end = rLine.find_first_of(" \t\r\n", start)) != std::string::npos) {
        if (end > start) {
            tokens.emplace_back(rLine.substr(start, end - start));
        }
        start = end + 1;
    }
    if (start < rLine.length()){
        tokens.emplace_back(rLine.substr(start));
    }
    return tokens;
}

/***********************************************************************************/
/***********************************************************************************/

std::string ObjIO::Info() const
{
    return "OBJ IO";
}

/***********************************************************************************/
/***********************************************************************************/

void ObjIO::PrintInfo(std::ostream& rOStream) const
{
    rOStream << Info();
}

/***********************************************************************************/
/***********************************************************************************/

void ObjIO::PrintData(std::ostream& rOStream) const
{

}

}  // namespace Kratos.