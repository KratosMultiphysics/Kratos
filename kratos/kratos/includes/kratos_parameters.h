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
class Parameters
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(Parameters);

    //ATTENTION: please DO NOT use this constructor. It assumes rapidjson and hence it should be considered as an implementation detail
    Parameters(rapidjson::Value* pvalue, boost::shared_ptr<rapidjson::Document> pdoc): mpvalue(pvalue),mpdoc(pdoc)
    {
    }

    Parameters(std::string json_string)
    {

        mpdoc =  boost::shared_ptr<rapidjson::Document>(new rapidjson::Document() );
        rapidjson::ParseResult ok = mpdoc->Parse<0>(json_string.c_str());

        if( !ok )
        {
            std::stringstream msg;
            msg << rapidjson::GetParseError_En(ok.Code()) << " offset of the error from the beginning of the string = " << ok.Offset() << std::endl;
            msg << "a much more explicative error message can be obtained by analysing the input string with an online analyzer such for example json lint" << std::endl;
            msg << "the value of the string that was attempted to parse is :" << std::endl << std::endl;
            msg << json_string;
            KRATOS_THROW_ERROR(std::invalid_argument, "error found in parsing the json_string, the value of the json string was: \n", msg.str());
        }

        mpvalue = (mpdoc.get());

    }

    /// Assignment operator.
    Parameters& operator=(Parameters const& rOther)
    {
        mpvalue->CopyFrom(*(rOther.GetUnderlyingStorage()), mpdoc->GetAllocator());

        return *this;
    }
    /// Copy constructor.
    Parameters(Parameters const& rOther)
    {
        mpdoc =  rOther.mpdoc;
        mpvalue = rOther.mpvalue;
    }

    //generates a clone of the current document
    Parameters Clone()
    {
        boost::shared_ptr<rapidjson::Document> pnew_cloned_doc =  boost::shared_ptr<rapidjson::Document>(new rapidjson::Document() );
        rapidjson::ParseResult ok = pnew_cloned_doc->Parse<0>(this->WriteJsonString().c_str());
        if( !ok )
        {
            std::stringstream msg;
            msg << rapidjson::GetParseError_En(ok.Code()) << " offset of the error from the beginning of the string = " << ok.Offset() << std::endl;
            msg << "a much more explicative error message can be obtained by analysing the input string " << std::endl;
            msg << "with an online analyzer such for example json lint" << std::endl;
            msg << "the value of the string that was attempted to parse is :" << std::endl << std::endl;
            msg << this->WriteJsonString();
            KRATOS_THROW_ERROR(std::invalid_argument, "error found in parsing the json_string, the value of the json string was: \n", msg.str());
        }
        return Parameters(pnew_cloned_doc.get(),pnew_cloned_doc);
    }

    /// Destructor.
    virtual ~Parameters() {}

    const std::string WriteJsonString() const
    {
        rapidjson::StringBuffer buffer;
        rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
        mpvalue->Accept(writer);
        return buffer.GetString();
    }


    const  std::string PrettyPrintJsonString() const
    {
        rapidjson::StringBuffer buffer;
        rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(buffer);
        mpvalue->Accept(writer);
        return buffer.GetString();
    }


    //*******************************************************************************************************
    Parameters GetValue(const std::string entry)
    {
        rapidjson::Value* pvalue = &((*mpvalue)[entry.c_str()]);

        return Parameters(pvalue, mpdoc);
    }
    Parameters operator[](const std::string entry)
    {
        return Parameters(&(*mpvalue)[entry.c_str()],mpdoc);
    }
    void SetValue(const std::string entry, const Parameters& other_value)
    {
        Parameters tmp(&(*mpvalue)[entry.c_str()],mpdoc);
        tmp.InternalSetValue(other_value);
    }



    //*******************************************************************************************************
    bool Has(const std::string entry)
    {
        return mpvalue->HasMember(entry.c_str());
    }


    bool IsNumber()
    {
        return mpvalue->IsNumber();
    }
    bool IsDouble()
    {
        return mpvalue->IsDouble();
    }
    bool IsInt()
    {
        return mpvalue->IsInt();
    }
    bool IsBool()
    {
        return mpvalue->IsBool();
    }
    bool IsString()
    {
        return mpvalue->IsString();
    }
    bool IsArray()
    {
        return mpvalue->IsArray();
    }
    bool IsSubParameter()
    {
        return mpvalue->IsObject();
    }
    double GetDouble()
    {
        if(mpvalue->IsNumber() == false) KRATOS_THROW_ERROR(std::invalid_argument,"argument must be a number","");
        return mpvalue->GetDouble();
    }
    int GetInt()
    {
        if(mpvalue->IsNumber() == false) KRATOS_THROW_ERROR(std::invalid_argument,"argument must be a number","");
        return mpvalue->GetInt();
    }
    bool GetBool()
    {
        if(mpvalue->IsBool() == false) KRATOS_THROW_ERROR(std::invalid_argument,"argument must be a bool","");
        return mpvalue->GetBool();
    }
    std::string GetString()
    {
        if(mpvalue->IsString() == false) KRATOS_THROW_ERROR(std::invalid_argument,"argument must be a string","");
        return mpvalue->GetString();
    }

    void SetDouble(const double value)
    {
        mpvalue->SetDouble(value);
    }
    void SetInt(const int value)
    {
        mpvalue->SetInt(value);
    }
    void SetBool(const bool value)
    {
        mpvalue->SetBool(value);
    }
    void SetString(const std::string value)
    {
        mpvalue->SetString(rapidjson::StringRef(value.c_str()));
    }


    //*******************************************************************************************************
    //methods for array
    unsigned int size()
    {
        if(mpvalue->IsArray() == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"size can only be queried if the value if of Array type","");
        return mpvalue->Size();
    }

    Parameters GetArrayItem(unsigned int index)
    {
        if(mpvalue->IsArray() == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"GetArrayItem only makes sense if the value if of Array type","")
            else
            {
                if(index >= mpvalue->Size())
                    KRATOS_THROW_ERROR(std::invalid_argument,"index exceeds array size. Index value is : ",index)
                    return Parameters(&(*mpvalue)[index],mpdoc);
            }
    }

    void SetArrayItem(unsigned int index, const Parameters& other_array_item)
    {
        if(mpvalue->IsArray() == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"SetArrayItem only makes sense if the value if of Array type","")
            else
            {
                if(index >= mpdoc->Size())
                    KRATOS_THROW_ERROR(std::invalid_argument,"index exceeds array size. Index value is : ",index)
#if RAPIDJSON_HAS_CXX11_RVALUE_REFS
                    (*mpvalue)[index] = rapidjson::Value(*other_array_item.GetUnderlyingStorage(), mpdoc->GetAllocator());
#else
                    (*mpvalue)[index].CopyFrom(*other_array_item.GetUnderlyingStorage(), mpdoc->GetAllocator());
#endif
            }
    }
    Parameters operator[](unsigned int index)
    {
        return this->GetArrayItem(index);
    }

    //ATTENTION: please DO NOT use this method. It is a low level accessor, and may change in the future
    rapidjson::Value* GetUnderlyingStorage() const
    {
        return mpvalue;
    }

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return this->PrettyPrintJsonString();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Parameters Object " << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {};

private:
    rapidjson::Value* mpvalue;
    boost::shared_ptr<rapidjson::Document> mpdoc;

    void InternalSetValue(const Parameters& other_value)
    {
        mpvalue->CopyFrom(*(other_value.GetUnderlyingStorage()), mpdoc->GetAllocator());
    }
};

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  Parameters& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Parameters& rThis)
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
