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

#if !defined(KRATOS_XML_IO_H_INCLUDED )
#define  KRATOS_XML_IO_H_INCLUDED


// System includes
#include <string>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <functional>
#include <algorithm>


// External includes


// Project includes
#include "includes/define.h"



namespace Kratos
{
///@addtogroup ApplicationNameApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class XmlIO
{
public:
    ///@name Type Definitions
    ///@{

    using BlockActionFunctionType=std::function<void(XmlIO&)>;

    class BlockInfo{
        std::string mBlockName;
        std::unordered_map<std::string, std::string> mAttributes;
    public:
        BlockInfo(std::stringstream& TheStartTag){
            std::cout << std::endl << "\"" << TheStartTag.str() << "\"" << std::endl;
            TheStartTag >> std::skipws >> mBlockName;
            //std::cout << "Reading block \'" << mBlockName << "\' with attributes: ";
            while (!TheStartTag.eof())
            {
                std::string key;
                std::string value;
                std::getline(TheStartTag, key, '=');
                key.erase(key.begin(), std::find_if(key.begin(), key.end(), [](unsigned char ch) {return !std::isspace(ch);}));
                if(key.empty())
                    continue; // This should be the end of the tag
                std::getline(TheStartTag, value, '"');
                std::getline(TheStartTag, value, '"');

                //std::cout << "[" << key << ":" << value << "]";
                mAttributes[key]=value;
            }
            std::cout << std::endl;            
        }
        std::string const& Name() const {
            return mBlockName;
        }

        bool HasAttribute(std::string const& AttributeName) const {
            return (mAttributes.find(AttributeName) != mAttributes.end());
        }

        bool IsStartTag() const {
            return (mBlockName[0] != '/');
        }

        std::string const& GetAttribute(std::string const& AttributeName) const {
            auto i_attribute = mAttributes.find(AttributeName);
            KRATOS_DEBUG_ERROR_IF(i_attribute == mAttributes.end()) << "The attribute \"" << AttributeName << "\" was not found in this block" << std::endl;

            return i_attribute->second;
        }

    };

    /// Pointer definition of XmlIO
    KRATOS_CLASS_POINTER_DEFINITION(XmlIO);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    XmlIO() = delete;

    /// Constructor with filename.
    XmlIO(std::string const& Filename);

    /// Constructor with stream.
    XmlIO(std::iostream* pInputStream);

    /// Destructor.
    virtual ~XmlIO(){
        delete mpInputStream;
    }

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void SetBlockAction(std::string const& BlockName, BlockActionFunctionType ActionFunction){
        mBlockActions[BlockName] = ActionFunction;
    }

    void Read();
    
    void ReadBlock(std::string const& BlockName, BlockActionFunctionType ActionFunction);

    template<typename TDataType>
    void ReadBlockContent(std::vector<TDataType>& Results){
        for(auto& i_result : Results){
            *mpInputStream >> i_result;
        }
    }

    template<typename TDataType>
    void ReadBlockContent(TDataType& Results){
        *mpInputStream >> Results;
    }

    ///@}
    ///@name Access
    ///@{

    BlockInfo const& GetCurrentBlockInfo(){
        KRATOS_DEBUG_ERROR_IF(mBlockInfoStack.empty()) << "Asking info before reading the block header" << std::endl;
        return mBlockInfoStack.back();
    }


    ///@}
    ///@name Inquiry
    ///@{

    bool HasBlockAction(std::string const& BlockName){
        return (mBlockActions.find(BlockName) != mBlockActions.end());
    }


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const;


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    std::iostream* mpInputStream;
    std::unordered_map<std::string, BlockActionFunctionType> mBlockActions;
    std::vector<BlockInfo> mBlockInfoStack;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void ReadStartTag();
    void ReadEndTag();
    BlockInfo ReadTag();


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    XmlIO& operator=(XmlIO const& rOther);

    /// Copy constructor.
    XmlIO(XmlIO const& rOther);


    ///@}

}; // Class XmlIO

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                XmlIO& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const XmlIO& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_XML_IO_H_INCLUDED  defined
