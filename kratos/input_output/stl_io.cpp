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
#include "includes/define.h"
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

    std::function<void(ModelPart&, NodesArrayType&)> create_entity_func;
    const std::string new_entity_type = mParameters["new_entity_type"].GetString();
    if (new_entity_type == "geometry") {
        create_entity_func = [](
            ModelPart& rThisModelPart,
            NodesArrayType& rIndexes) {
                rThisModelPart.CreateNewGeometry("Triangle3D3", rIndexes);
            };
    } else if (new_entity_type == "element") {
        create_entity_func = [this](
            ModelPart& rThisModelPart,
            NodesArrayType& rIndexes) {
                rThisModelPart.CreateNewElement("Element3D3N", this->mNextElementId++, rIndexes, rThisModelPart.pGetProperties(0));
            };
    } else if (new_entity_type == "condition") {
        create_entity_func = [this](
            ModelPart& rThisModelPart,
            NodesArrayType& rIndexes) {
                rThisModelPart.CreateNewCondition("SurfaceCondition3D3N", this->mNextConditionId++, rIndexes, rThisModelPart.pGetProperties(0));
            };
    } else  {
        KRATOS_ERROR << "Invalid new entity type " << new_entity_type << std::endl;
    }

    while(!mpInputStream->eof()) {
        ReadSolid(rThisModelPart,create_entity_func);
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

void StlIO::ReadSolid(
    ModelPart& rThisModelPart,
    const std::function<void(ModelPart&, NodesArrayType&)>& rCreateEntityFunctor
    )
{
    std::string word;

    *mpInputStream >> word; // Reading solid or eof
    if(mpInputStream->eof())
        return;

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

    if(word == "") // empty solid name is valid in STL format
        word = "main";

    auto& sub_model_part = rThisModelPart.CreateSubModelPart(word);

    *mpInputStream >> word; // Reading facet or endsolid

    while(word == "facet"){
        ReadFacet(sub_model_part, rCreateEntityFunctor);
        *mpInputStream >> word; // Reading facet or endsolid
    }

    KRATOS_ERROR_IF(word != "endsolid") << "Invalid stl file. Solid block should be closed with \"endsolid\" keyword but \"" << word << "\" was found" << std::endl;
    std::getline(*mpInputStream, word); // Reading solid name 
}

void StlIO::ReadFacet(
    ModelPart& rThisModelPart,
    const std::function<void(ModelPart&, NodesArrayType&)>& rCreateEntityFunctor
    )
{
    std::string word;

    ReadKeyword("normal");

    std::getline(*mpInputStream, word); // Reading n_i n_j n_k

    *mpInputStream >> word; // Reading outer or endfacet

    while(word == "outer"){
        ReadLoop(rThisModelPart, rCreateEntityFunctor);
        *mpInputStream >> word; // Reading outer or endfacet
    }

    KRATOS_ERROR_IF(word != "endfacet") << "Invalid stl file. facet block should be closed with \"endfacet\" keyword but \"" << word << "\" was found" << std::endl;
}

void StlIO::ReadLoop(
    ModelPart & rThisModelPart,
    const std::function<void(ModelPart&, NodesArrayType&)>& rCreateEntityFunctor
    )
{
    std::string word;

    ReadKeyword("loop");

    *mpInputStream >> word; // Reading vertex or endloop

    NodesArrayType temp_geom_nodes;
    while(word == "vertex"){
        Point coordinates = ReadPoint();
        temp_geom_nodes.push_back(rThisModelPart.CreateNewNode(mNextNodeId++, coordinates[0], coordinates[1], coordinates[2] ));
        *mpInputStream >> word; // Reading vertex or endloop
    }
    const std::string new_entity_type = mParameters["new_entity_type"].GetString();
    rCreateEntityFunctor(rThisModelPart,temp_geom_nodes);
    KRATOS_ERROR_IF(word != "endloop") << "Invalid stl file. loop block should be closed with \"endloop\" keyword but \"" << word << "\" was found" << std::endl;
}

Point StlIO::ReadPoint()
{
    Point result;
    std::string word;
    for(int i = 0 ; i < 3 ; i++){
        *mpInputStream >> word;
        result[i] = std::stod(word);
    }
    return result;
}

void StlIO::ReadKeyword(std::string const& Keyword)
{
    std::string word;

    *mpInputStream >> word; // Reading keyword
    KRATOS_ERROR_IF(word != Keyword) << "Invalid stl file. Looking for  \"" << Keyword << "\" keyword but \"" << word << "\" were found." << std::endl;
}

}  // namespace Kratos.