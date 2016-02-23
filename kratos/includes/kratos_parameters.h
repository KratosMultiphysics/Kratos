//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//


#if !defined(KRATOS_KRATOS_PARAMETERS_H_INCLUDED )
#define  KRATOS_KRATOS_PARAMETERS_H_INCLUDED



// System includes

#include <string>
#include <iostream>
#include <sstream>

// External includes


// Project includes
#include "includes/define.h"
#include "rapidjson/document.h"
#include "rapidjson/error/en.h"
#include "rapidjson/writer.h"
#include "rapidjson/prettywriter.h"
#include "rapidjson/stringbuffer.h"

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
class ParameterValue
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ParameterValue);
        
    //ATTENTION: please DO NOT use this constructor. It assumes rapidjson and hence it should be considered as an implementation detail
    ParameterValue(rapidjson::Value& rvalue): mrvalue(rvalue)
    {
    }
    
    const std::string WriteJsonString() const
    {
        rapidjson::StringBuffer buffer;
        rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
        mrvalue.Accept(writer);
        return buffer.GetString();
    }


    const  std::string PrettyPrintJsonString() const   
    {
        rapidjson::StringBuffer buffer;
        rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(buffer);
        mrvalue.Accept(writer);
        return buffer.GetString();
    }

    
    //*******************************************************************************************************
    ParameterValue GetValue(const std::string entry)
    {
        return ParameterValue(mrvalue[entry.c_str()]);
    }
    ParameterValue operator[](const std::string entry)
    {
        return ParameterValue(mrvalue[entry.c_str()]);
    }
    void SetValue(const std::string entry, const ParameterValue& other_value)
    {
        KRATOS_THROW_ERROR(std::logic_error,"it is not allowed to set directly a value ","")
//         rapidjson::ParseResult ok = mrvalue.Parse<0>(other_value.WriteJsonString().c_str());
    }
    //*******************************************************************************************************
    bool Has(const std::string entry)
    {
       return mrvalue.HasMember(entry.c_str());
    }
    
    
    bool IsNumber()
    {
        return mrvalue.IsNumber();
    }
    bool IsDouble()
    {
        return mrvalue.IsDouble();
    }
    bool IsInt()
    {
        return mrvalue.IsInt();
    }
    bool IsBool()
    {
        return mrvalue.IsBool();
    }
    bool IsString()
    {
        return mrvalue.IsString();
    }
    bool IsArray()
    {
        return mrvalue.IsArray();
    }
    bool IsSubParameter()
    {
        return mrvalue.IsObject();
    }
    double GetDouble()
    {
        if(mrvalue.IsNumber() == false) KRATOS_THROW_ERROR(std::invalid_argument,"argument must be a number","");
        return mrvalue.GetDouble();
    }
    int GetInt()
    {
        if(mrvalue.IsNumber() == false) KRATOS_THROW_ERROR(std::invalid_argument,"argument must be a number","");
        return mrvalue.GetInt();
    }
    bool GetBool()
    {
        if(mrvalue.IsBool() == false) KRATOS_THROW_ERROR(std::invalid_argument,"argument must be a bool","");
        return mrvalue.GetBool();
    }
    std::string GetString()
    {
        if(mrvalue.IsString() == false) KRATOS_THROW_ERROR(std::invalid_argument,"argument must be a string","");
        return mrvalue.GetString();
    }

    void SetDouble(const double value)
    {
        mrvalue.SetDouble(value);
    }
    void SetInt(const int value)
    {
        mrvalue.SetInt(value);
    }
    void SetBool(const bool value)
    {
        mrvalue.SetBool(value);
    }
    void SetString(const std::string value)
    {
         mrvalue.SetString(rapidjson::StringRef(value.c_str()));
    }
   
    
    //*******************************************************************************************************
    //methods for array
    unsigned int size()
    {
        if(mrvalue.IsArray() == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"size can only be queried if the value if of Array type","");
        return mrvalue.Size();
    }
    
    ParameterValue GetArrayItem(unsigned int index)
    { 
        if(mrvalue.IsArray() == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"GetArrayItem only makes sense if the value if of Array type","")
        else
        {
            if(index >= mrvalue.Size())
                 KRATOS_THROW_ERROR(std::invalid_argument,"index exceeds array size. Index value is : ",index)
            return ParameterValue(mrvalue[index]);
        }
    }
     
    void SetArrayItem(unsigned int index, const ParameterValue& other_array_item)
    {
        KRATOS_THROW_ERROR(std::logic_error,"not allowed to set directly an array value","")
    }
       
    ParameterValue operator[](unsigned int index)
    {
        return this->GetArrayItem(index);
    }
private:
    rapidjson::Value& mrvalue;
};

/// Short class definition.
/** This class provides a wrapper to Json
 * the idea is to abstract the kratos from the underlying implementation of the json reader
 *
 * initially using rapidjson lib, but may well vary in the future
*/
class KratosParameters
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosParameters
    KRATOS_CLASS_POINTER_DEFINITION(KratosParameters);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosParameters(std::string json_string)
    {
        rapidjson::ParseResult ok = md.Parse<0>(json_string.c_str());
        if( !ok )
        {
            std::stringstream msg;
            msg << rapidjson::GetParseError_En(ok.Code()) << " offset of the error from the beginning of the string = " << ok.Offset() << std::endl;
            msg << "a much more explicative error message can be obtained by analysing the input string with an online analyzer such for example json lint" << std::endl;
            msg << "the value of the string that was attempted to parse is :" << std::endl << std::endl;
            msg << json_string;
            KRATOS_THROW_ERROR(std::invalid_argument, "error found in parsing the json_string, the value of the json string was: \n", msg.str());
        }
    }

    /// Assignment operator.
    KratosParameters& operator=(KratosParameters const& rOther)
    {   
        rapidjson::ParseResult ok = md.Parse<0>(rOther.WriteJsonString().c_str());
        if( !ok )
        {
            std::stringstream msg;
            msg << rapidjson::GetParseError_En(ok.Code()) << " offset of the error from the beginning of the string = " << ok.Offset() << std::endl;
            msg << "a much more explicative error message can be obtained by analysing the input string " << std::endl;
            msg << "with an online analyzer such for example json lint" << std::endl;
            msg << "the value of the string that was attempted to parse is :" << std::endl << std::endl;
            msg << rOther.WriteJsonString();
            KRATOS_THROW_ERROR(std::invalid_argument, "error found in parsing the json_string, the value of the json string was: \n", msg.str());
        } 
        return *this;
    }
    /// Copy constructor.
    KratosParameters(KratosParameters const& rOther)
    {
        rapidjson::ParseResult ok = md.Parse<0>(rOther.WriteJsonString().c_str());
        if( !ok )
        {
            std::stringstream msg;
            msg << rapidjson::GetParseError_En(ok.Code()) << " offset of the error from the beginning of the string = " << ok.Offset() << std::endl;
            msg << "a much more explicative error message can be obtained by analysing the input string " << std::endl;
            msg << "with an online analyzer such for example json lint" << std::endl;
            msg << "the value of the string that was attempted to parse is :" << std::endl << std::endl;
            msg << rOther.WriteJsonString();
            KRATOS_THROW_ERROR(std::invalid_argument, "error found in parsing the json_string, the value of the json string was: \n", msg.str());
        }
    }

    /// Destructor.
    virtual ~KratosParameters() {}

    
    ///@}
    ///@name Operators
    ///@{
    ParameterValue GetValue(const std::string entry)
    {
        return ParameterValue(md[entry.c_str()]);        
    }

    bool Has(const std::string entry)
    {
       return md.HasMember(entry.c_str());
    }



    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{
    const std::string WriteJsonString() const
    {
        rapidjson::StringBuffer buffer;
        rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
        md.Accept(writer);
        return buffer.GetString();
    }


    const  std::string PrettyPrintJsonString() const   
    {
        rapidjson::StringBuffer buffer;
        rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(buffer);
        md.Accept(writer);
        return buffer.GetString();
    }


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const {return this->PrettyPrintJsonString();}

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "KratosParameters Object " << Info();}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {};


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
    rapidjson::Document md;


    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{




    ///@}

}; // Class KratosParameters

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  KratosParameters& rThis) {return rIStream;}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const KratosParameters& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block





}  // namespace Kratos.

#endif // KRATOS_KRATOS_PARAMETERS_H_INCLUDED  defined 


