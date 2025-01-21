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
        "open_mode"                       : "read",
        "entity_type"                     : "element",
        "decompose_quads_into_triangles"  : false,
        "normal_as_historical"            : false
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
    mFirstElementId = mNextElementId;

    mNextConditionId = block_for_each<MaxReduction<std::size_t>>(
        rThisModelPart.GetRootModelPart().Conditions(),
        [](Condition& rCondition) { return rCondition.Id();}) + 1;
    mFirstConditionId = mNextConditionId;

    // Read vertices, normals, and faces
    const bool normal_as_historical = mParameters["normal_as_historical"].GetBool();
    const std::string entity_type = mParameters["entity_type"].GetString();
    const bool decompose_quads_into_triangles = mParameters["decompose_quads_into_triangles"].GetBool();
    ReadVerticesAndFaces(rThisModelPart, entity_type, normal_as_historical, decompose_quads_into_triangles);
}

/***********************************************************************************/
/***********************************************************************************/

void ObjIO::WriteModelPart(const ModelPart& rThisModelPart)
{
    // To know if we are retrieving normals as historical variables
    const bool normal_as_historical = mParameters["normal_as_historical"].GetBool();

    // Write the header
    //  #    Number of vertices: XXXXX
    //  #    Number of faces: YYYYYY
    *mpInputStream << "# Number of vertices: " << rThisModelPart.NumberOfNodes() << "\n";
    *mpInputStream << "# Number of faces: " << rThisModelPart.NumberOfElements() << "\n\n";

    // Write vertices
    *mpInputStream << "# Vertices" << std::endl;
    for (const auto& r_node : rThisModelPart.Nodes()) {
        *mpInputStream << "v " << r_node.X() << " " << r_node.Y() << " " << r_node.Z() << std::endl;
    }
    *mpInputStream << "\n";

    // Write faces
    // NOTE: We will assume that the nodes are ordered and start at 1
    *mpInputStream << "# Faces" << std::endl;
    std::string entity_type = mParameters["entity_type"].GetString();
    const std::size_t number_of_geometries = rThisModelPart.NumberOfGeometries();
    const std::size_t number_of_elements = rThisModelPart.NumberOfElements();
    const std::size_t number_of_conditions = rThisModelPart.NumberOfConditions();
    // Check if the entity type is valid
    if (entity_type == "geometry" && number_of_geometries == 0) {
        if (number_of_elements > 0) {
            entity_type = "element";
        } else if (number_of_conditions > 0) {
            entity_type = "condition";
        }
    } else if (entity_type == "element" && number_of_elements == 0) {
        if (number_of_geometries > 0) {
            entity_type = "geometry";
        } else if (number_of_conditions > 0) {
            entity_type = "condition";
        }
    } else if (entity_type == "condition" && number_of_conditions == 0) {
        if (number_of_geometries > 0) {
            entity_type = "geometry";
        } else if (number_of_elements > 0) {
            entity_type = "element";
        }
    }
    // Finally, write faces
    if (entity_type == "geometry") {
        for (const auto& r_geometry : rThisModelPart.Geometries()) {
            *mpInputStream << "f";
            for (const auto& r_node : r_geometry) {
                *mpInputStream << " " << r_node.Id();
            }
            *mpInputStream << "\n";
        }
    } else if (entity_type == "element") {
        for (const auto& r_element : rThisModelPart.Elements()) {
            *mpInputStream << "f";
            for (const auto& r_node : r_element.GetGeometry()) {
                *mpInputStream << " " << r_node.Id();
            }
            *mpInputStream << "\n";
        }
    } else if (entity_type == "condition") {
        for (const auto& r_condition : rThisModelPart.Conditions()) {
            *mpInputStream << "f";
            for (const auto& r_node : r_condition.GetGeometry()) {
                *mpInputStream << " " << r_node.Id();
            }
            *mpInputStream << "\n";
        }
    } else  {
        KRATOS_ERROR << "Invalid entity type " << entity_type << std::endl;
    }
    *mpInputStream << "\n";

    // Write normals
    *mpInputStream << "# Normals" << std::endl;
    if (normal_as_historical) {
        for (const auto& r_node : rThisModelPart.Nodes()) {
            const array_1d<double, 3>& r_normal = r_node.FastGetSolutionStepValue(NORMAL);
            *mpInputStream << "vn " << r_normal[0] << " " << r_normal[1] << " " << r_normal[2] << std::endl;
        }
    } else {
        for (const auto& r_node : rThisModelPart.Nodes()) {
            const array_1d<double, 3>& r_normal = r_node.GetValue(NORMAL);
            *mpInputStream << "vn " << r_normal[0] << " " << r_normal[1] << " " << r_normal[2] << std::endl;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void ObjIO::ReadVerticesAndFaces(
    ModelPart& rThisModelPart,
    const std::string& rEntityType,
    const bool NormalAsHistoricalVariable,
    const bool DecomposeQuadrilateral
    )
{
    std::string line;
    while (std::getline(*mpInputStream, line)) {
        // Trim leading and trailing whitespace
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);

        // Ignore empty lines and comments
        // But we check if: #    Number of vertices: XXXXX
        //                  #    Number of faces: YYYYYY
        //  are defined in the header
        if (line.empty()) {
            continue;
        } else if (line[0] == '#') {
            if (line.find("Number of vertices") != std::string::npos) {
                const std::size_t num_vertices = std::stoul(line.substr(line.find(":") + 1));
                rThisModelPart.Nodes().reserve(rThisModelPart.NumberOfNodes() + num_vertices);
            } else if (line.find("Number of faces") != std::string::npos) {
                const std::size_t num_faces = std::stoul(line.substr(line.find(":") + 1));
                if (rEntityType == "geometry") {
                    rThisModelPart.Geometries().reserve(rThisModelPart.NumberOfGeometries() + num_faces);
                } else if (rEntityType == "element") {
                    rThisModelPart.Elements().reserve(rThisModelPart.NumberOfElements() + num_faces);
                } else if (rEntityType == "condition") {
                    rThisModelPart.Conditions().reserve(rThisModelPart.NumberOfConditions() + num_faces);
                }
            }
            continue;
        }

        // Check line type
        if (line.substr(0, 2) == "v ") {
            ParseVertexLine(rThisModelPart, line);
        } else if (line.substr(0, 3) == "vn ") {
            ParseNormalLine(rThisModelPart, line, NormalAsHistoricalVariable);
        } else if (line.substr(0, 2) == "f ") {
            ParseFaceLine(rThisModelPart, line, rEntityType, DecomposeQuadrilateral);
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
    const IndexType node_id = mNormalCounter + mFirstNodeId;
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
    const std::string& rLine,
    const std::string& rEntityType,
    const bool DecomposeQuadrilateral
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
    if (rEntityType == "geometry") {
        KRATOS_ERROR_IF(number_of_tokens > 5) << "Only support for triangles and quads in geometry creation. Number of nodes: " << number_of_tokens - 1 << std::endl;
        if (DecomposeQuadrilateral && number_of_tokens == 5) {
            // Decompose quad into two triangles
            const std::vector<IndexType> tri1 = {nodes_ids[0], nodes_ids[1], nodes_ids[2]};
            const std::vector<IndexType> tri2 = {nodes_ids[0], nodes_ids[2], nodes_ids[3]};
            rThisModelPart.CreateNewGeometry("Triangle3D3", tri1);
            rThisModelPart.CreateNewGeometry("Triangle3D3", tri2);
        } else {
            rThisModelPart.CreateNewGeometry(SupportedGeometries[number_of_tokens - 4], nodes_ids);
        }
    } else if (rEntityType == "element") {
        auto p_prop = rThisModelPart.pGetProperties(0);
        if (DecomposeQuadrilateral && number_of_tokens == 5) {
            // Decompose quad into two triangles
            const std::vector<IndexType> tri1 = {nodes_ids[0], nodes_ids[1], nodes_ids[2]};
            const std::vector<IndexType> tri2 = {nodes_ids[0], nodes_ids[2], nodes_ids[3]};
            rThisModelPart.CreateNewElement("Element3D3N", mNextElementId++, tri1, p_prop);
            rThisModelPart.CreateNewElement("Element3D3N", mNextElementId++, tri2, p_prop);
        } else {
            rThisModelPart.CreateNewElement("Element3D" + std::to_string(number_of_tokens - 1) + "N", mNextElementId++, nodes_ids, p_prop);
        }
    } else if (rEntityType == "condition") {
        auto p_prop = rThisModelPart.pGetProperties(0);
        if (DecomposeQuadrilateral && number_of_tokens == 5) {
            // Decompose quad into two triangles
            const std::vector<IndexType> tri1 = {nodes_ids[0], nodes_ids[1], nodes_ids[2]};
            const std::vector<IndexType> tri2 = {nodes_ids[0], nodes_ids[2], nodes_ids[3]};
            rThisModelPart.CreateNewCondition("SurfaceCondition3D3N", mNextConditionId++, tri1, p_prop);
            rThisModelPart.CreateNewCondition("SurfaceCondition3D3N", mNextConditionId++, tri2, p_prop);
        } else {
            rThisModelPart.CreateNewCondition("SurfaceCondition3D" + std::to_string(number_of_tokens - 1) + "N", mNextConditionId++, nodes_ids, p_prop);
        }
    } else  {
        KRATOS_ERROR << "Invalid new entity type " << rEntityType << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

std::vector<std::string> ObjIO::Tokenize(const std::string& rLine)
{
    // Tokenize the line by spaces
    std::vector<std::string> tokens;
    std::string::size_type start = 0;
    std::string::size_type end = 0;

    // Find the position of the comment character '#'
    std::string line_no_comments = rLine;
    std::string::size_type comment_pos = rLine.find('#');
    if (comment_pos != std::string::npos) {
        line_no_comments = rLine.substr(0, comment_pos); // Remove the comment part
    }

    // Tokenize the line by spaces
    while ((end = line_no_comments.find_first_of(" \t\r\n", start)) != std::string::npos) {
        if (end > start) {
            tokens.emplace_back(line_no_comments.substr(start, end - start));
        }
        start = end + 1;
    }
    if (start < line_no_comments.length()){
        tokens.emplace_back(line_no_comments.substr(start));
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