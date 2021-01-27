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
//

// System includes
#include <fstream>

// External includes


// Project includes
#include "includes/define.h"
#include "input_output/xml_io.h"


namespace Kratos
{


    /// Constructor with  filenames.
    XmlIO::XmlIO(std::string const& Filename)
        : mpInputStream(new std::fstream(Filename, std::ios::in))  {
                KRATOS_ERROR_IF(mpInputStream->fail()) << "Could not open the input file : " << Filename << std::endl;
        }

    XmlIO::XmlIO(std::iostream* pInputStream) : mpInputStream(pInputStream){}

    void XmlIO::Read(){

        std::size_t starting_nested_level = mBlockInfoStack.size();
        while((!mpInputStream->eof()) && (starting_nested_level <= mBlockInfoStack.size())){
            auto tag = ReadTag();
            if(mpInputStream->eof())
                break;

            if(tag.IsStartTag()){
                std::cout << "Starting block \"" << tag.Name() << "\"" <<  std::endl;
                mBlockInfoStack.push_back(tag);
                auto& current_block_name = GetCurrentBlockInfo().Name();
                KRATOS_ERROR_IF(!HasBlockAction(current_block_name)) << " The block \"" << current_block_name << "\" is not supported" << std::endl;
                auto& action = mBlockActions[current_block_name];
                action(*this);
            }
            else {
                
                std::cout << "Closing block " << tag.Name() << std::endl;
                std::string current_block_end_tag = "/" + GetCurrentBlockInfo().Name();

                KRATOS_ERROR_IF(tag.Name() != current_block_end_tag) << "Closing the block \"" << GetCurrentBlockInfo().Name() << "\" with differnt end tag: <" << tag.Name() << ">" << std::endl;

                mBlockInfoStack.pop_back();

                if (mBlockInfoStack.empty()) // Closing the document
                    break;
            }

        }        
    }
 
    void XmlIO::ReadBlock(std::string const& BlockName, BlockActionFunctionType ActionFunction){
        ReadStartTag();
        auto& current_block_name = GetCurrentBlockInfo().Name();
        KRATOS_ERROR_IF(current_block_name != BlockName) << "Expected \"" << BlockName << "\" block but found \"" << current_block_name << "\"" << std::endl;
        ActionFunction(*this);

        ReadEndTag();
    }

    /// Turn back information as a string.
    std::string XmlIO::Info() const
    {
        return "XML IO";
    }

    /// Print information about this object.
    void XmlIO::PrintInfo(std::ostream& rOStream) const{
        rOStream << Info();
    }

    /// Print object's data.
    void XmlIO::PrintData(std::ostream& rOStream) const{

    }

    void XmlIO::ReadStartTag(){
        std::string tag;
        std::getline(*mpInputStream, tag, '<');
        std::getline(*mpInputStream, tag, '>');

        std::stringstream tag_stream(tag);

        BlockInfo current_block_info(tag_stream);
        mBlockInfoStack.push_back(current_block_info);
    }

    void XmlIO::ReadEndTag(){
        std::string tag;
        std::getline(*mpInputStream, tag, '<');
        std::getline(*mpInputStream, tag, '>');
        
        std::string current_block_end_tag = "/" + GetCurrentBlockInfo().Name();

        KRATOS_ERROR_IF(tag != current_block_end_tag) << "Closing the block \"" << GetCurrentBlockInfo().Name() << "\" with differnt end tag: <" << tag << ">" << std::endl;

        mBlockInfoStack.pop_back();
    }

    XmlIO::BlockInfo XmlIO::ReadTag(){
        std::string tag;
        std::getline(*mpInputStream, tag, '<');
        std::getline(*mpInputStream, tag, '>');

        std::stringstream tag_stream(tag);

        BlockInfo current_block_info(tag_stream);
        return current_block_info;
    }


}  // namespace Kratos.


