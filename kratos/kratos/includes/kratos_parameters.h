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
        if(this->Has(entry) == false) KRATOS_THROW_ERROR(std::invalid_argument,"getting a value that does not exist. entry string : ",entry);
        rapidjson::Value* pvalue = &((*mpvalue)[entry.c_str()]);

        return Parameters(pvalue, mpdoc);
    }
    Parameters operator[](const std::string entry)
    {
        return Parameters(&(*mpvalue)[entry.c_str()],mpdoc);
    }
    void SetValue(const std::string entry, const Parameters& other_value)
    {
        if(this->Has(entry) == false) KRATOS_THROW_ERROR(std::invalid_argument,"value must exist to be set. Use AddValue instead","");
        Parameters tmp(&(*mpvalue)[entry.c_str()],mpdoc);
        tmp.InternalSetValue(other_value);
    }

    void AddValue(const std::string entry, const Parameters& other_value)
    {
        if(this->Has(entry) == false)
        {
            rapidjson::Value tmp;
            tmp.CopyFrom(*(other_value.GetUnderlyingStorage()), mpdoc->GetAllocator()); //this will be moved away
            rapidjson::Value name(entry.c_str(), mpdoc->GetAllocator()); //rhis will be moved away
            this->mpvalue->AddMember(name, tmp, mpdoc->GetAllocator());
        }
    }
    Parameters AddEmptyValue(const std::string entry)
    {
        if(this->Has(entry) == false)
        {
            rapidjson::Value tmp;
            rapidjson::Value name(entry.c_str(), mpdoc->GetAllocator()); //rhis will be moved away
            this->mpvalue->AddMember(name, tmp, mpdoc->GetAllocator());
        }
        return this->GetValue(entry);
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
        rapidjson::Value tmp(value.c_str(), mpdoc->GetAllocator());
        *mpvalue = tmp;
//         mpvalue->SetString(rapidjson::StringRef(value.c_str()));
//        mpvalue->SetString(value.c_str(), value.length());
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

    /**This function is designed to verify that the parameters under testing match the
     * form prescribed by the defaults.
     * If the parameters contain values that do not appear in the defaults, an error is thrown, 
     * whereas if a parameter is found in the defaults but not in the Parameters been tested,
     * it is copied to the parameters.
     * 
     * this version of the function only walks one level, without descending in the branches
     */
    void ValidateAndAssignDefaults(Parameters& defaults)
    {
        KRATOS_TRY


        //first verifies that all the enries in the current parameters
        //have a correspondance in the defaults.
        //if it is not the case throw an error
        for (rapidjson::Value::ConstMemberIterator itr = this->mpvalue->MemberBegin(); itr != this->mpvalue->MemberEnd(); ++itr)
        {
            std::string item_name = itr->name.GetString();

            if(!defaults.Has(item_name) )
            {
                std::stringstream msg;
                msg << "******************************************************************************************************" << std::endl;
                msg << "the item with name \"" << item_name << "\" is present in this Parameters but NOT in the default values" << std::endl;
                msg << "******************************************************************************************************" << std::endl;
                msg << "hence Validation fails" << std::endl;
                msg << "parameters being validated are : " << std::endl;
                msg << this->PrettyPrintJsonString() << std::endl;
                msg << "defaults against which the current parameters are validated are :" << std::endl;
                msg << defaults.PrettyPrintJsonString() << std::endl;
                KRATOS_THROW_ERROR(std::invalid_argument,"",msg.str());
            }

            bool type_coincides = false;
            rapidjson::Value* value_defaults = (defaults[item_name.c_str()]).GetUnderlyingStorage();
            if(itr->value.IsInt() && value_defaults->IsNumber()) type_coincides = true;
            if(itr->value.IsBool() && value_defaults->IsBool()) type_coincides = true;
            if(itr->value.IsDouble() && value_defaults->IsDouble()) type_coincides = true;
            if(itr->value.IsArray() && value_defaults->IsArray()) type_coincides = true;
            if(itr->value.IsString() && value_defaults->IsString()) type_coincides = true;
            if(itr->value.IsObject() && value_defaults->IsObject()) type_coincides = true;

            if(type_coincides == false)
            {
                std::stringstream msg;
                msg << "******************************************************************************************************" << std::endl;
                msg << "the item with name :\"" << item_name << "\" does not have the same type as the corresponding one in the default values" << std::endl;
                msg << "******************************************************************************************************" << std::endl;
                msg << "parameters being validated are : " << std::endl;
                msg << this->PrettyPrintJsonString() << std::endl;
                msg << "defaults against which the current parameters are validated are :" << std::endl;
                msg << defaults.PrettyPrintJsonString() << std::endl;
                KRATOS_THROW_ERROR(std::invalid_argument,"",msg.str());
            }
        }



        //now iterate over all the defaults. In the case a default value is not assigned in the current Parameters
        //add an item copying its value
        if(defaults.IsSubParameter())
        {
            for (rapidjson::Value::MemberIterator itr = defaults.mpvalue->MemberBegin(); itr != defaults.mpvalue->MemberEnd(); ++itr)
            {
                std::string item_name = itr->name.GetString();
                if(!this->Has(item_name))
                {
                    rapidjson::Value* pvalue = &itr->value;

                    this->AddValue(item_name, Parameters(pvalue, defaults.mpdoc));
                }
            }
        }


        KRATOS_CATCH("")
    }

    /**This function is designed to verify that the parameters under testing match the
     * form prescribed by the defaults.
     * If the parameters contain values that do not appear in the defaults, an error is thrown, 
     * whereas if a parameter is found in the defaults but not in the Parameters been tested,
     * it is copied to the parameters.
     * 
     * this version walks and validates the entire json tree below 
     * the point at which the function is called
    */
    void RecursivelyValidateAndAssignDefaults(Parameters& defaults)
    {
        KRATOS_TRY


        //first verifies that all the enries in the current parameters
        //have a correspondance in the defaults.
        //if it is not the case throw an error
        for (rapidjson::Value::ConstMemberIterator itr = this->mpvalue->MemberBegin(); itr != this->mpvalue->MemberEnd(); ++itr)
        {
            std::string item_name = itr->name.GetString();

            if(!defaults.Has(item_name) )
            {
                std::stringstream msg;
                msg << "the item with name \"" << item_name << "\" is present in this Parameters but NOT in the default values" << std::endl;
                msg << "hence Validation fails" << std::endl;
                msg << "parameters being validated are : " << std::endl;
                msg << this->PrettyPrintJsonString() << std::endl;
                msg << "defaults against which the current parameters are validated are :" << std::endl;
                msg << defaults.PrettyPrintJsonString() << std::endl;
                KRATOS_THROW_ERROR(std::invalid_argument,"",msg.str());
            }

            bool type_coincides = false;
            rapidjson::Value* value_defaults = (defaults[item_name.c_str()]).GetUnderlyingStorage();
            if(itr->value.IsInt() && value_defaults->IsInt()) type_coincides = true;
            if(itr->value.IsBool() && value_defaults->IsBool()) type_coincides = true;
            if(itr->value.IsDouble() && value_defaults->IsDouble()) type_coincides = true;
            if(itr->value.IsArray() && value_defaults->IsArray()) type_coincides = true;
            if(itr->value.IsString() && value_defaults->IsString()) type_coincides = true;
            if(itr->value.IsObject() && value_defaults->IsObject()) type_coincides = true;

            if(type_coincides == false)
            {
                std::stringstream msg;
                msg << "the item with name :\"" << item_name << "\" does not have the same type as the corresponding one in the default values" << std::endl;
                msg << "parameters being validated are : " << std::endl;
                msg << this->PrettyPrintJsonString() << std::endl;
                msg << "defaults against which the current parameters are validated are :" << std::endl;
                msg << defaults.PrettyPrintJsonString() << std::endl;
                KRATOS_THROW_ERROR(std::invalid_argument,"",msg.str());
            }
            //now walk the tree recursively
            if(itr->value.IsObject())
            {
                Parameters subobject = (*this)[item_name];
                Parameters defaults_subobject = defaults[item_name];
                subobject.ValidateAndAssignDefaults(defaults_subobject);
            }
        }



        //now iterate over all the defaults. In the case a default value is not assigned in the current Parameters
        //add an item copying its value
        if(defaults.IsSubParameter())
        {
            for (rapidjson::Value::MemberIterator itr = defaults.mpvalue->MemberBegin(); itr != defaults.mpvalue->MemberEnd(); ++itr)
            {
                std::string item_name = itr->name.GetString();
                if(!this->Has(item_name))
                {
                    rapidjson::Value* pvalue = &itr->value;

                    this->AddValue(item_name, Parameters(pvalue, defaults.mpdoc));
                }

                //now walk the tree recursively
                if(itr->value.IsObject())
                {
                    Parameters subobject = (*this)[item_name];
                    Parameters defaults_subobject = defaults[item_name];
                    subobject.ValidateAndAssignDefaults(defaults_subobject);
                }
            }
        }


        KRATOS_CATCH("")
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
