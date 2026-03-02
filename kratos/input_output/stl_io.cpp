//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

// System includes

// External includes

// Project includes
#include "input_output/stl_io.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

namespace Kratos
{

StlIO::StlIO(const std::filesystem::path& rFilename, Parameters ThisParameters)
    : mParameters(ThisParameters)
{
    mParameters.ValidateAndAssignDefaults(GetDefaultParameters());
    Kratos::shared_ptr<std::fstream> p_file = Kratos::make_shared<std::fstream>();

    // Check if the file extension is .stl, if not add it
    std::filesystem::path filePath = rFilename;
    if(filePath.extension() != ".stl") {
        filePath += ".stl";
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

StlIO::StlIO(
    Kratos::shared_ptr<std::iostream> pInputStream,
    Parameters ThisParameters)
    : IO(),
      mParameters(ThisParameters),
      mpInputStream(pInputStream)
{
    mParameters.ValidateAndAssignDefaults(GetDefaultParameters());
}

Parameters StlIO::GetDefaultParameters()
{
    return Parameters(R"({
        "open_mode"       : "read",
        "new_entity_type" : "geometry"
    })" );
}

void StlIO::ReadModelPart(ModelPart& rThisModelPart)
{
    // Ensure the model part has properties
    if(!rThisModelPart.RecursivelyHasProperties(0)) {
        rThisModelPart.CreateNewProperties(0);
    }

    // Compute the next IDs for nodes, elements, and conditions
    mNextNodeId = block_for_each<MaxReduction<std::size_t>>(
        rThisModelPart.GetRootModelPart().Nodes(),
        [](NodeType& rNode) { return rNode.Id();}) + 1;

    mNextElementId = block_for_each<MaxReduction<std::size_t>>(
        rThisModelPart.GetRootModelPart().Elements(),
        [](Element& rElement) { return rElement.Id();}) + 1;

    mNextConditionId = block_for_each<MaxReduction<std::size_t>>(
        rThisModelPart.GetRootModelPart().Conditions(),
        [](Condition& rCondition) { return rCondition.Id();}) + 1;

    // Read the file until the end of the file is reached. The ReadSolid function will read one solid at a time, which corresponds to one submodelpart in the model part.
    while(!mpInputStream->eof()) {
        ReadSolid(rThisModelPart);
    }
}

void StlIO::WriteModelPart(const ModelPart& rThisModelPart)
{
    // Get data communicator
    const auto& r_data_communicator = rThisModelPart.GetCommunicator().GetDataCommunicator();

    // If rank is not zeo we remove the stream to avoid conflict
    if (r_data_communicator.Rank() != 0) {
        mpInputStream = nullptr;
    }
    r_data_communicator.Barrier();

    /* Write the solid block */
    // Write header of the file
    if (r_data_communicator.Rank() == 0) {
        (*mpInputStream) << "solid " << rThisModelPart.Name() << "\n";
    }

    // Depending on MPI/serial
    if (r_data_communicator.IsDistributed()) { // MPI
        // Write the elements blocks
        WriteEntityBlockMPI(rThisModelPart.Elements(), r_data_communicator);

        // Write the conditions blocks
        WriteEntityBlockMPI(rThisModelPart.Conditions(), r_data_communicator);

        // Write the geometries blocks
        WriteGeometryBlockMPI(rThisModelPart.Geometries(), r_data_communicator);
    } else { // Serial
        // Write the elements blocks
        WriteEntityBlock(rThisModelPart.Elements());

        // Write the conditions blocks
        WriteEntityBlock(rThisModelPart.Conditions());

        // Write the geometries blocks
        WriteGeometryBlock(rThisModelPart.Geometries());
    }

    // Write footer of the file
    if (r_data_communicator.Rank() == 0) {
        (*mpInputStream) << "endsolid\n";
    }
}

template<class TContainerType>
void StlIO::WriteEntityBlock(const TContainerType& rThisEntities)
{
    // Retrieve reference of the stream
    auto& r_stream = *mpInputStream;

    // Write facets
    std::size_t num_degenerate_geometries = 0;
    for (auto& r_entity : rThisEntities) {
        const auto& r_geometry = r_entity.GetGeometry();
        if (IsValidGeometry(r_geometry, num_degenerate_geometries)) {
            WriteFacet(r_geometry, r_stream);
        }
    }
    KRATOS_WARNING_IF("STL-IO", num_degenerate_geometries > 0)
        << "Model part contained " << num_degenerate_geometries
        << " geometries with area = 0.0, skipping these geometries." << std::endl;
}

void StlIO::WriteGeometryBlock(const GeometryContainerType& rThisGeometries)
{
    // Retrieve reference of the stream
    auto& r_stream = *mpInputStream;

    // Write facets
    std::size_t num_degenerate_geometries = 0;
    for (auto& r_geometry : rThisGeometries) {
        if (IsValidGeometry(r_geometry, num_degenerate_geometries)) {
            WriteFacet(r_geometry, r_stream);
        }
    }
    KRATOS_WARNING_IF("STL-IO", num_degenerate_geometries > 0)
        << "Model part contained " << num_degenerate_geometries
        << " geometries with area = 0.0, skipping these geometries." << std::endl;
}

template<class TContainerType>
void StlIO::WriteEntityBlockMPI(
    const TContainerType& rThisEntities,
    const DataCommunicator& rDataCommunicator
    )
{
    // Current partition stream
    std::stringstream ss;

    // Write facets
    std::size_t num_degenerate_geometries = 0;
    for (auto& r_entity : rThisEntities) {
        const auto& r_geometry = r_entity.GetGeometry();
        if (IsValidGeometry(r_geometry, num_degenerate_geometries)) {
            WriteFacet(r_geometry, ss);
        }
    }

    // Sum all partitions
    unsigned int converted_num_degenerate_geometries(num_degenerate_geometries);
    converted_num_degenerate_geometries = rDataCommunicator.SumAll(converted_num_degenerate_geometries);

    KRATOS_WARNING_IF("STL-IO", converted_num_degenerate_geometries > 0)
        << "Model part contained " << converted_num_degenerate_geometries
        << " geometries with area = 0.0, skipping these geometries." << std::endl;

    // Getting number of entities
    unsigned int number_of_entities = rThisEntities.size();
    number_of_entities = rDataCommunicator.SumAll(number_of_entities);

    // Retrieve rank and pass to rank 0
    if (number_of_entities > 0) {
        const int rank = rDataCommunicator.Rank();
        const int tag = 0;
        if (rank == 0) {
            // Retrieve reference of the stream
            auto& r_stream = *mpInputStream;

            // Add the corresponding 0 rank string
            r_stream << ss.str();

            // Now receive from other partitions
            for (int i = 1; i < rDataCommunicator.Size(); ++i) {
                std::string recv_string;
                rDataCommunicator.Recv(recv_string, i, tag);
                r_stream << recv_string;
            }
        } else {
            // Send the stream to rank 0
            rDataCommunicator.Send(ss.str(), 0, tag);
        }
    }
}

void StlIO::WriteGeometryBlockMPI(
    const GeometryContainerType& rThisGeometries,
    const DataCommunicator& rDataCommunicator
    )
{
    // Current partition stream
    std::stringstream ss;

    // Write facets
    std::size_t num_degenerate_geometries = 0;
    for (auto& r_geometry : rThisGeometries) {
        if (IsValidGeometry(r_geometry, num_degenerate_geometries)) {
            WriteFacet(r_geometry, ss);
        }
    }

    // Sum all partitions
    unsigned int converted_num_degenerate_geometries(num_degenerate_geometries); 
    converted_num_degenerate_geometries = rDataCommunicator.SumAll(converted_num_degenerate_geometries);

    KRATOS_WARNING_IF("STL-IO", converted_num_degenerate_geometries > 0) 
        << "Model part contained " << converted_num_degenerate_geometries
        << " geometries with area = 0.0, skipping these geometries." << std::endl;

    // Getting number of entities
    unsigned int number_of_geometries = rThisGeometries.size();
    number_of_geometries = rDataCommunicator.SumAll(number_of_geometries);

    // Retrieve rank and pass to rank 0
    if (number_of_geometries > 0) {
        const int rank = rDataCommunicator.Rank();
        const int tag = 0;
        if (rank == 0) {
            // Retrieve reference of the stream
            auto& r_stream = *mpInputStream;

            // Add the corresponding 0 rank string
            r_stream << ss.str();

            // Now receive from other partitions
            for (int i = 1; i < rDataCommunicator.Size(); ++i) {
                std::string recv_string;
                rDataCommunicator.Recv(recv_string, i, tag);
                r_stream << recv_string;
            }
        } else {
            // Send the stream to rank 0
            rDataCommunicator.Send(ss.str(), 0, tag);
        }
    }
}

template<class TStreamType>
void StlIO::WriteFacet(
    const GeometryType& rGeom,
    TStreamType& rStream
    )
{
    const auto& r_unit_normal = rGeom.UnitNormal(rGeom.Center());
    rStream << "    facet normal " << r_unit_normal[0] << " " << r_unit_normal[1] << " " << r_unit_normal[2] << "\n";
    rStream << "        outer loop\n";

    for (int i = 0; i < 3; i++) {
        const auto& r_node = rGeom[i];
        rStream << "           vertex " << r_node.X() << " " << r_node.Y() << " " << r_node.Z() << "\n";
    }

    rStream << "        endloop\n";
    rStream << "    endfacet\n";
}

/// Turn back information as a string.
std::string StlIO::Info() const
{
    return "STL IO";
}

/// Print information about this object.
void StlIO::PrintInfo(std::ostream& rOStream) const{
    rOStream << Info();
}

/// Print object's data.
void StlIO::PrintData(std::ostream& rOStream) const{

}

bool StlIO::IsValidGeometry(
    const Geometry<Node>& rGeometry,
    IndexType& rNumDegenerateGeos
    ) const
{
    // restrict to triangles only for now
    const bool is_triangle = (
        rGeometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D3 ||
        rGeometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D6);
    const bool area_greater_than_zero = rGeometry.Area() > std::numeric_limits<double>::epsilon();
    if (!area_greater_than_zero && is_triangle) {
        rNumDegenerateGeos++;
    }
    return (is_triangle && area_greater_than_zero);
}

void StlIO::ReadSolid(ModelPart& rThisModelPart)
{
    std::string word;

    *mpInputStream >> word; // Reading solid or eof
    if(mpInputStream->eof()) {
        return;
    }

    // Estimate number of lines of the current solid block, then restore stream position
    const std::streampos initial_position = mpInputStream->tellg();
    std::size_t num_lines = 0;
    if (initial_position != std::streampos(-1)) {
        std::string line;
        while (std::getline(*mpInputStream, line)) {
            ++num_lines;
            if (line.find("endsolid") != std::string::npos) {
                break;
            }
        }

        // Reset stream state/position after counting
        mpInputStream->clear();
        mpInputStream->seekg(initial_position);
        KRATOS_ERROR_IF_NOT(mpInputStream->good()) << "Failed to restore STL stream position after line counting." << std::endl;
    }
    const std::size_t estimated_number_of_facets = num_lines / 7; // Each facet block has 7 lines (facet normal, outer loop, 3 vertices, endloop, endfacet)
    const std::size_t estimated_nodes = estimated_number_of_facets * 3;
    std::vector<Node::Pointer> new_nodes;
    new_nodes.reserve(estimated_nodes);

    KRATOS_ERROR_IF(word != "solid") << "Invalid stl file. Solid block should begin with \"solid\" keyword but \"" << word << "\" was found" << std::endl;
    std::getline(*mpInputStream, word); // Reading solid name to be the model part name

    // Remove Comments
    for (const auto& symbol : {"COMMENT:", ";"}) {
        const auto position = word.find(symbol);
        if (position != word.npos) {
            word.erase(position);
        }
    }
    word.erase(std::remove_if(
        word.begin(),
        word.end(),
        [](auto character) -> bool {
            return character == '\r' || character == '\n' || character == ' '; //remove eol's and white spaces
        }
    ), word.end());

    // Empty solid name is valid in STL format
    if(word == "") {
        word = "main";
    }

    auto& r_sub_model_part = rThisModelPart.CreateSubModelPart(word);

    *mpInputStream >> word; // Reading facet or endsolid

    while(word == "facet") {
        ReadFacet(r_sub_model_part, new_nodes);

        *mpInputStream >> word; // Reading facet or endsolid
    }

    // Add nodes to the model part after reading all facets to avoid searching for existing nodes during the reading process
    const std::size_t num_new_nodes = new_nodes.size();
    auto& r_nodes = r_sub_model_part.Nodes();
    auto& r_nodes_data = r_nodes.GetContainer();
    r_nodes_data.insert(r_nodes_data.end(), new_nodes.begin(), new_nodes.end());
    r_nodes.SetSortedPartSize(r_nodes_data.size());
    auto& r_root_nodes = rThisModelPart.Nodes();
    auto& r_root_nodes_data = r_root_nodes.GetContainer();
    r_root_nodes_data.insert(r_root_nodes_data.end(), new_nodes.begin(), new_nodes.end());
    r_root_nodes.SetSortedPartSize(r_root_nodes_data.size());
    const std::string new_entity_type = mParameters["new_entity_type"].GetString();

    // Define TLS for parallel execution of entity creation
    struct TLS {
        NodesArrayType triangle_nodes;
        Properties::Pointer p_properties = nullptr;

        TLS(Properties::Pointer pProperties = nullptr) : p_properties(pProperties) {
            triangle_nodes.reserve(3); // reserving for 3 nodes as STL only supports triangles
            for (std::size_t i = 0; i < 3; ++i) {
                triangle_nodes.push_back(nullptr);
            }
        }
    };

    // Create entities in parallel using the created nodes. Each facet corresponds to 3 consecutive nodes in the new_nodes array
    if (new_entity_type == "geometry") {
        // const GeometryType& r_clone_geometry = KratosComponents<GeometryType>::Get("Triangle3D3N");

        // // Create geometries in parallel using the created nodes. Each facet corresponds to 3 consecutive nodes in the new_nodes array
        // const std::vector<GeometryType::Pointer> new_geometries = IndexPartition<IndexType>(num_new_nodes / 3).for_each<AccumReduction<GeometryType::Pointer>>(TLS(), [&, this](const IndexType i, TLS& rTLS) {
        //     const std::size_t node_index = i * 3;
        //     std::swap(rTLS.triangle_nodes(0), new_nodes[node_index]);
        //     std::swap(rTLS.triangle_nodes(1), new_nodes[node_index + 1]);
        //     std::swap(rTLS.triangle_nodes(2), new_nodes[node_index + 2]);
        //     return r_clone_geometry.Create(rTLS.triangle_nodes);
        // });

        // // Add geometries
        // r_sub_model_part.AddGeometries(new_geometries.begin(), new_geometries.end());

        // Memory issue with parallel creation of geometries, creating geometries sequentially for now with ModelPart's CreateNewGeometry
        TLS tls;
         for (IndexType i = 0; i < num_new_nodes; i += 3) {
            std::swap(tls.triangle_nodes(0), new_nodes[i]);
            std::swap(tls.triangle_nodes(1), new_nodes[i + 1]);
            std::swap(tls.triangle_nodes(2), new_nodes[i + 2]);
            r_sub_model_part.CreateNewGeometry("Triangle3D3", tls.triangle_nodes);
        }
    } else if (new_entity_type == "element") {
        const Element& r_clone_element = KratosComponents<Element>::Get("Element3D3N");

        // Create elements in parallel using the created nodes. Each facet corresponds to 3 consecutive nodes in the new_nodes array
        const std::vector<Element::Pointer> new_elements = IndexPartition<IndexType>(num_new_nodes / 3).for_each<AccumReduction<Element::Pointer>>(TLS(r_sub_model_part.pGetProperties(0)), [&, this](const IndexType i, TLS& rTLS) {
            const std::size_t node_index = i * 3;
            std::swap(rTLS.triangle_nodes(0), new_nodes[node_index]);
            std::swap(rTLS.triangle_nodes(1), new_nodes[node_index + 1]);
            std::swap(rTLS.triangle_nodes(2), new_nodes[node_index + 2]);
            return r_clone_element.Create(mNextElementId + i, rTLS.triangle_nodes, rTLS.p_properties);
        });
        mNextElementId += new_elements.size();

        // Add elements
        r_sub_model_part.AddElements(new_elements.begin(), new_elements.end());
    } else if (new_entity_type == "condition") {
        const Condition& r_clone_condition = KratosComponents<Condition>::Get("SurfaceCondition3D3N");

        // Create conditions in parallel using the created nodes. Each facet corresponds to 3 consecutive nodes in the new_nodes array
        const std::vector<Condition::Pointer> new_conditions = IndexPartition<IndexType>(num_new_nodes / 3).for_each<AccumReduction<Condition::Pointer>>(TLS(r_sub_model_part.pGetProperties(0)), [&, this](const IndexType i, TLS& rTLS) {
            const std::size_t node_index = i * 3;
            std::swap(rTLS.triangle_nodes(0), new_nodes[node_index]);
            std::swap(rTLS.triangle_nodes(1), new_nodes[node_index + 1]);
            std::swap(rTLS.triangle_nodes(2), new_nodes[node_index + 2]);
            return r_clone_condition.Create(mNextConditionId + i, rTLS.triangle_nodes, rTLS.p_properties);
        });
        mNextConditionId += new_conditions.size();

        // Add conditions
        r_sub_model_part.AddConditions(new_conditions.begin(), new_conditions.end());
    } else  {
        KRATOS_ERROR << "Invalid new entity type " << new_entity_type << std::endl;
    }

    KRATOS_ERROR_IF(word != "endsolid") << "Invalid stl file. Solid block should be closed with \"endsolid\" keyword but \"" << word << "\" was found" << std::endl;
    std::getline(*mpInputStream, word); // Reading solid name
}
void StlIO::ReadFacet(
    ModelPart& rThisModelPart,
    std::vector<Node::Pointer>& rNewNodes
    )
{
    std::string word;

    ReadKeyword("normal");

    std::getline(*mpInputStream, word); // Reading n_i n_j n_k

    *mpInputStream >> word; // Reading outer or endfacet

    while(word == "outer"){
        ReadLoop(rThisModelPart, rNewNodes);
        *mpInputStream >> word; // Reading outer or endfacet
    }

    KRATOS_ERROR_IF(word != "endfacet") << "Invalid stl file. facet block should be closed with \"endfacet\" keyword but \"" << word << "\" was found" << std::endl;
}

void StlIO::ReadLoop(
    ModelPart & rThisModelPart,
    std::vector<Node::Pointer>& rNewNodes
    )
{
    std::string word;

    ReadKeyword("loop");

    *mpInputStream >> word; // Reading vertex or endloop

    // Adding the nodes array of the facet directly to the provided array
    std::array<double, 3> coordinates;
    while(word == "vertex") {
        ReadPoint(coordinates);
        rNewNodes.push_back(Kratos::make_intrusive<Node>(mNextNodeId++, coordinates[0], coordinates[1], coordinates[2] ));
        *mpInputStream >> word; // Reading vertex or endloop
    }
    KRATOS_ERROR_IF(word != "endloop") << "Invalid stl file. loop block should be closed with \"endloop\" keyword but \"" << word << "\" was found" << std::endl;
}

void StlIO::ReadPoint(std::array<double, 3>& rCoordinates)
{
    for (double& r_coordinate : rCoordinates) {
        *mpInputStream >> r_coordinate;
    }
}

void StlIO::ReadKeyword(std::string const& Keyword)
{
    std::string word;

    *mpInputStream >> word; // Reading keyword
    KRATOS_ERROR_IF(word != Keyword) << "Invalid stl file. Looking for  \"" << Keyword << "\" keyword but \"" << word << "\" were found." << std::endl;
}

}  // namespace Kratos.