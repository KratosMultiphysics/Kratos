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
    Parameters ThisParameters) 
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
        "open_mode"       : "read",
        "new_entity_type" : "element"
    })" );
}

/***********************************************************************************/
/***********************************************************************************/

void ObjIO::ReadModelPart(ModelPart& rThisModelPart)
{
    if(!rThisModelPart.RecursivelyHasProperties(0)) {
        rThisModelPart.CreateNewProperties(0);
    }

    mNextNodeId = block_for_each<MaxReduction<std::size_t>>(
        rThisModelPart.GetRootModelPart().Nodes(),
        [](NodeType& rNode) { return rNode.Id();}) + 1;

    mNextElementId = block_for_each<MaxReduction<std::size_t>>(
        rThisModelPart.GetRootModelPart().Elements(),
        [](Element& rElement) { return rElement.Id();}) + 1;

    mNextConditionId = block_for_each<MaxReduction<std::size_t>>(
        rThisModelPart.GetRootModelPart().Conditions(),
        [](Condition& rCondition) { return rCondition.Id();}) + 1;

    // Clear the node index map and normals list
    mNodeIndexMap.clear();
    mNormalsList.clear();

    // We don't know how many vertices and normals there are, but we can reserve some space
    mNodeIndexMap.reserve(1000);
    mNormalsList.reserve(1000);

    ReadVerticesAndFaces(rThisModelPart);
}

/***********************************************************************************/
/***********************************************************************************/

void ObjIO::WriteModelPart(const ModelPart& rThisModelPart)
{
    KRATOS_ERROR << "WriteModelPart is not implemented for ObjIO." << std::endl;
}

void ObjIO::ReadVerticesAndFaces(ModelPart& rThisModelPart)
{
    std::string line;
    while (std::getline(*mpInputStream, line))
    {
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
            ParseNormalLine(line);
        } else if (line.substr(0, 2) == "f ") {
            ParseFaceLine(rThisModelPart, line);
        }
        // Other types can be handled if needed
    }
}

/***********************************************************************************/
/***********************************************************************************/

void ObjIO::ParseVertexLine(ModelPart& rThisModelPart, const std::string& line)
{
    // The line is expected to be: v x y z [w]
    std::vector<std::string> tokens = Tokenize(line);

    KRATOS_ERROR_IF(tokens.size() < 4) << "Invalid vertex line: " << line << std::endl;

    const double x = std::stod(tokens[1]);
    const double y = std::stod(tokens[2]);
    const double z = std::stod(tokens[3]);

    // Create new node
    NodeType::Pointer p_node = rThisModelPart.CreateNewNode(mNextNodeId++, x, y, z);

    // Store the node ID in the index map
    mNodeIndexMap.push_back(p_node->Id());
}

/***********************************************************************************/
/***********************************************************************************/

void ObjIO::ParseNormalLine(const std::string& line)
{
    // The line is expected to be: vn x y z
    std::vector<std::string> tokens = Tokenize(line);

    KRATOS_ERROR_IF(tokens.size() < 4) << "Invalid vertex normal line: " << line << std::endl;

    array_1d<double, 3> normal;
    normal[0] = std::stod(tokens[1]);
    normal[1] = std::stod(tokens[2]);
    normal[2] = std::stod(tokens[3]);

    // Store the normal in the normals list
    mNormalsList.push_back(normal);
}

/***********************************************************************************/
/***********************************************************************************/

void ObjIO::ParseFaceLine(ModelPart& rThisModelPart, const std::string& line)
{
    // The line is expected to be: f v1 v2 v3 ...
    std::vector<std::string> tokens = Tokenize(line);

    KRATOS_ERROR_IF(tokens.size() < 4) << "Invalid face line: " << line << std::endl;

    NodesArrayType temp_geom_nodes;

    // For each vertex index in the face
    for (std::size_t i = 1; i < tokens.size(); ++i) {
        std::string vertex_token = tokens[i];
        // OBJ indices can be of the form v, v/vt, v//vn, v/vt/vn
        // We need to extract the vertex index and the normal index
        std::string::size_type first_slash = vertex_token.find('/');
        std::string::size_type second_slash = vertex_token.find('/', first_slash + 1);

        std::string vertex_index_str;
        std::string normal_index_str;

        if (first_slash == std::string::npos) {
            // Format: v
            vertex_index_str = vertex_token;
        } else if (second_slash == first_slash + 1) {
            // Format: v//vn
            vertex_index_str = vertex_token.substr(0, first_slash);
            normal_index_str = vertex_token.substr(second_slash + 1);
        } else {
            // Format: v/vt or v/vt/vn
            vertex_index_str = vertex_token.substr(0, first_slash);
            if (second_slash != std::string::npos) {
                // v/vt/vn
                normal_index_str = vertex_token.substr(second_slash + 1);
            }
        }

        // OBJ indices start at 1, adjust to zero-based index
        const int vertex_index = std::stoi(vertex_index_str) - 1;
        KRATOS_ERROR_IF(vertex_index < 0 || static_cast<std::size_t>(vertex_index) >= mNodeIndexMap.size()) << "Invalid vertex index in face line: " << line << std::endl;

        IndexType node_id = mNodeIndexMap[vertex_index];
        NodeType::Pointer p_node = rThisModelPart.pGetNode(node_id);

        // Assign normal if available
        if (!normal_index_str.empty()) {
            int normal_index = std::stoi(normal_index_str) - 1;
            KRATOS_ERROR_IF(normal_index < 0 || static_cast<std::size_t>(normal_index) >= mNormalsList.size()) << "Invalid normal index in face line: " << line << std::endl;

            p_node->SetValue(NORMAL, mNormalsList[normal_index]);
        }

        temp_geom_nodes.push_back(p_node);
    }

    // Create element or condition or geometry (TODO: Support proper naming in geometries)
    const std::string new_entity_type = mParameters["new_entity_type"].GetString();
    // if (new_entity_type == "geometry") {
    //     rThisModelPart.CreateNewGeometry("", temp_geom_nodes);
    // } else 
    if (new_entity_type == "element") {
        rThisModelPart.CreateNewElement("Element3D3N", mNextElementId++, temp_geom_nodes, rThisModelPart.pGetProperties(0));
    } else if (new_entity_type == "condition") {
        rThisModelPart.CreateNewCondition("SurfaceCondition3D3N", mNextConditionId++, temp_geom_nodes, rThisModelPart.pGetProperties(0));
    } else  {
        KRATOS_ERROR << "Invalid new entity type " << new_entity_type << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

std::vector<std::string> ObjIO::Tokenize(const std::string& line)
{
    std::vector<std::string> tokens;
    std::string::size_type start = 0;
    std::string::size_type end = 0;

    while ((end = line.find_first_of(" \t\r\n", start)) != std::string::npos) {
        if (end > start) {
            tokens.emplace_back(line.substr(start, end - start));
        }
        start = end + 1;
    }
    if (start < line.length()){
        tokens.emplace_back(line.substr(start));
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