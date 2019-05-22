//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "input_output/stl_io.h"


namespace Kratos
{

    /// Constructor with  filenames.
    StlIO::StlIO(std::string const& Filename)
        : mpInputStream(new std::fstream(Filename, std::ios::in))  {
                KRATOS_ERROR_IF(mpInputStream->fail()) << "Could not open the input file : " << Filename << std::endl;
        }

    StlIO::StlIO(std::iostream* pInputStream) : IO(), mpInputStream(pInputStream){}

    void StlIO::ReadModelPart(ModelPart & rThisModelPart)
    {
        if(!rThisModelPart.HasProperties(0))
            rThisModelPart.CreateNewProperties(0);
            
        while(!mpInputStream->eof())
            ReadSolid(rThisModelPart);
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

        KRATOS_ERROR_IF(word != "solid") << "Invalid stl file. Solid block should begin with \"solid\" keyword but \"" << word << "\" was founded" << std::endl;
        std::getline(*mpInputStream, word); // Reading solid name to be the model part name

        word.erase(word.begin(), std::find_if(word.begin(), word.end(), [](int ch) {return !std::isspace(ch);})); // Triming the leading spaces
    	size_t comment_position = word.find("COMMENT");
 
	    if (comment_position != std::string::npos)
            word.erase(comment_position); // Removing the comment
        
        word.erase(word.find_last_not_of(" \t\n\r\f\v") + 1); // triming the trailing whitespaces 

	

        if(word == "") // empty solid name is valid in STL format
            word = "main";

        KRATOS_INFO("Input") << "Reading Solid '" << word << "'" << std::endl;

        auto& sub_model_part = rThisModelPart.CreateSubModelPart(word);

        *mpInputStream >> word; // Reading facet or endsolid
        
        while(word == "facet"){
            ReadFacet(sub_model_part);
            *mpInputStream >> word; // Reading facet or endsolid
        }

        KRATOS_ERROR_IF(word != "endsolid") << "Invalid stl file. Solid block should be closed with \"endsolid\" keyword but \"" << word << "\" was founded" << std::endl;
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

        KRATOS_ERROR_IF(word != "endfacet") << "Invalid stl file. facet block should be closed with \"endfacet\" keyword but \"" << word << "\" was founded" << std::endl;
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
        KRATOS_ERROR_IF(word != "endloop") << "Invalid stl file. loop block should be closed with \"endloop\" keyword but \"" << word << "\" was founded" << std::endl;
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

    void StlIO::ReadKeyword(std::string const& Keyword){
        std::string word;

        *mpInputStream >> word; // Reading keyword
        KRATOS_ERROR_IF(word != Keyword) << "Invalid stl file. Looking for  \"" << Keyword << "\" keyword but \"" << word << "\" were found." << std::endl;

    }

}  // namespace Kratos.


