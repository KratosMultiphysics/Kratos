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

namespace Kratos
{

/// Constructor with  filenames.
StlIO::StlIO(std::filesystem::path const& Filename)
    : IO()
{
    std::fstream* pFile = new std::fstream();
    std::fstream::openmode OpenMode = std::fstream::in;

    pFile->open(Filename.c_str(), OpenMode);

    KRATOS_ERROR_IF_NOT(pFile->is_open()) << "Could not open the input file  : " << Filename << std::endl;

    // Store the pointer as a regular std::iostream
    mpInputStream = pFile;
}

StlIO::StlIO(std::filesystem::path const& Filename, const Flags Options)
    : mOptions(Options)
{
    std::fstream* pFile = new std::fstream();
    std::fstream::openmode OpenMode = std::fstream::in;

    // Set the mode
    if (mOptions.Is(IO::READ)) {
        OpenMode = std::fstream::in;
    } else if (mOptions.Is(IO::APPEND)) {
        OpenMode = std::fstream::in | std::fstream::app;
    } else if (mOptions.Is(IO::WRITE)) {
        OpenMode = std::fstream::out;
    } else {
        // If none of the READ, WRITE or APPEND are defined we will take READ as
        // default.
        OpenMode = std::fstream::in;
    }

    pFile->open(Filename.c_str(), OpenMode);

    KRATOS_ERROR_IF_NOT(pFile->is_open()) << "Could not open the input file  : " << Filename << std::endl;

    // Store the pointer as a regular std::iostream
    mpInputStream = pFile;
}

StlIO::StlIO(std::iostream* pInputStream) 
    : IO(), 
      mpInputStream(pInputStream)
{

}

void StlIO::ReadModelPart(ModelPart & rThisModelPart)
{
    if(!rThisModelPart.HasProperties(0))
        rThisModelPart.CreateNewProperties(0);
        
    while(!mpInputStream->eof())
        ReadSolid(rThisModelPart);
}

void StlIO::WriteModelPart(ModelPart & rThisModelPart)
{    
    (*mpInputStream) << "solid compensated_skin_part\n";


    for (auto & rElem : rThisModelPart.Elements()) {
        const auto & rGeom = rElem.GetGeometry();
        const auto & rUnitNormal = rGeom.UnitNormal(rGeom.Center());
        (*mpInputStream) << "    facet normal " << rUnitNormal[0] << " " << rUnitNormal[1] << " " << rUnitNormal[2] << "\n";
    
        (*mpInputStream) << "        outer loop\n";

        for(std::size_t i = 0; i < rGeom.size(); i++ )
        {
            auto & rNode = rGeom.GetPoint(i);
            (*mpInputStream) << "           vertex " << rNode.X() << " " << rNode.Y() << " " << rNode.Z() << "\n";
        }

        (*mpInputStream) << "        endloop\n";
        (*mpInputStream) << "    endfacet\n";
    }

    (*mpInputStream) << "endsolid\n";
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

void StlIO::ReadSolid(ModelPart & rThisModelPart)
{
    std::string word;

    *mpInputStream >> word; // Reading solid or eof
    if(mpInputStream->eof())
        return;

    KRATOS_ERROR_IF(word != "solid") << "Invalid stl file. Solid block should begin with \"solid\" keyword but \"" << word << "\" was found" << std::endl;
    std::getline(*mpInputStream, word); // Reading solid name to be the model part name

    word.erase(word.begin(), std::find_if(word.begin(), word.end(), [](int ch) {return !std::isspace(ch);})); // Triming the leading spaces

    if(word == "") // empty solid name is valid in STL format
        word = "main";

    auto& sub_model_part = rThisModelPart.CreateSubModelPart(word);

    *mpInputStream >> word; // Reading facet or endsolid
    
    KRATOS_WATCH(word);
    while(word == "facet"){
        ReadFacet(sub_model_part);
        *mpInputStream >> word; // Reading facet or endsolid
    }

    KRATOS_ERROR_IF(word != "endsolid") << "Invalid stl file. Solid block should be closed with \"endsolid\" keyword but \"" << word << "\" was found" << std::endl;
    std::getline(*mpInputStream, word); // Reading solid name 
}

void StlIO::ReadFacet(ModelPart & rThisModelPart)
{
    std::string word;

    ReadKeyword("normal");

    std::getline(*mpInputStream, word); // Reading n_i n_j n_k

    *mpInputStream >> word; // Reading outer or endfacet

    while(word == "outer"){
        ReadLoop(rThisModelPart);
        *mpInputStream >> word; // Reading outer or endfacet
    }

    KRATOS_ERROR_IF(word != "endfacet") << "Invalid stl file. facet block should be closed with \"endfacet\" keyword but \"" << word << "\" was found" << std::endl;
}

void StlIO::ReadLoop(ModelPart & rThisModelPart)
{
    std::string word;

    ReadKeyword("loop");

    *mpInputStream >> word; // Reading vertex or endloop
    std::size_t node_id = rThisModelPart.GetRootModelPart().NumberOfNodes() + 1;
    std::size_t element_id = rThisModelPart.GetRootModelPart().NumberOfElements() + 1;
    Element::NodesArrayType temp_element_nodes;
    while(word == "vertex"){
        Point coordinates = ReadPoint();
        temp_element_nodes.push_back(rThisModelPart.CreateNewNode(node_id++, coordinates[0], coordinates[1], coordinates[2] ));
        *mpInputStream >> word; // Reading vertex or endloop
    }
    rThisModelPart.CreateNewElement("Element3D3N", element_id, temp_element_nodes, rThisModelPart.pGetProperties(0));
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


