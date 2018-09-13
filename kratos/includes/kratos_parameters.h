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
#include "includes/ublas_interface.h"
#include "json/json.hpp"                      //import nlohmann json library

namespace Kratos
{
///@addtogroup KratosCore
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
/**
 * @class Parameters
 * @ingroup KratosCore
 * @brief This class provides to Kratos a data structure for I/O based on the standard of JSON
 * @details In computing, JavaScript Object Notation or JSON is an open-standard file format that uses human-readable text to transmit data objects consisting of attribute–value pairs and array data types (or any other serializable value). It is a very common data format used for asynchronous browser–server communication, including as a replacement for XML in some AJAX-style systems. More info: https://json.org/
 * This class uses nlohmann JSON header only library
 * @author Riccardo Rossi
 */
class Parameters
{
private:
    using json = nlohmann::json;
    ///@name Nested clases
    ///@{
    class iterator_adaptor : public std::iterator<std::forward_iterator_tag, Parameters>
    {
        using value_iterator = json::iterator;
        value_iterator mValueIterator;
        std::unique_ptr<Parameters> mpParameters;
    public:
        iterator_adaptor(value_iterator it,  Kratos::shared_ptr<json> proot) :mValueIterator(it), mpParameters(new Parameters(&(*it), proot)) {}
        iterator_adaptor(const iterator_adaptor& it) : mValueIterator(it.mValueIterator),  mpParameters(new Parameters(*(it.mpParameters))) {}
        iterator_adaptor& operator++()
        {
            mValueIterator++;
            return *this;
        }
        iterator_adaptor operator++(int)
        {
            iterator_adaptor tmp(*this);
            operator++();
            return tmp;
        }
        bool operator==(const iterator_adaptor& rhs) const
        {
            return mValueIterator == rhs.mValueIterator;
        }
        bool operator!=(const iterator_adaptor& rhs) const
        {
            return mValueIterator != rhs.mValueIterator;
        }
        Parameters& operator*() const
        {
            mpParameters->mpValue = &(*mValueIterator);
            return *mpParameters;
        }
        Parameters* operator->() const
        {
            mpParameters->mpValue = &(*mValueIterator);
            return mpParameters.get();
        }
        value_iterator& base()
        {
            return mValueIterator;
        }
        value_iterator const& base() const
        {
            return mValueIterator;
        }
        const std::string name()
        {
            return mValueIterator.key();
        }
    };

    class const_iterator_adaptor : public std::iterator<std::forward_iterator_tag, Parameters>
    {
        using value_iterator = json::const_iterator;
        value_iterator mValueIterator;
        std::unique_ptr<Parameters> mpParameters;
    public:
        const_iterator_adaptor(value_iterator it,  Kratos::shared_ptr<json> proot) :mValueIterator(it), mpParameters(new Parameters(const_cast<json*>(&(*it)), proot)) {}
        //TODO: use copy constructor in the following method
        const_iterator_adaptor(const const_iterator_adaptor& it) : mValueIterator(it.mValueIterator), mpParameters(new Parameters(*(it.mpParameters))) {}
        const_iterator_adaptor& operator++()
        {
            mValueIterator++;
            mpParameters->mpValue = const_cast<json*>( &(*mValueIterator) );
            return *this;
        }
        const_iterator_adaptor operator++(int)
        {
            const_iterator_adaptor tmp(*this);
            operator++();
            return tmp;
        }
        bool operator==(const const_iterator_adaptor& rhs) const
        {
            return mValueIterator == rhs.mValueIterator;
        }
        bool operator!=(const const_iterator_adaptor& rhs) const
        {
            return mValueIterator != rhs.mValueIterator;
        }
        const Parameters& operator*() const
        {
            //mpParameters->mpValue = &(*mValueIterator);
            return *mpParameters;

        }
        const Parameters* operator->() const
        {
            return mpParameters.get();
        }
        value_iterator& base()
        {
            return mValueIterator;
        }
        value_iterator const& base() const
        {
            return mValueIterator;
        }
        const std::string name()
        {
            return mValueIterator.key();
        }
    };

    ///@}

public:
    ///@name Type Definitions
    ///@{

    /// Index definition
    typedef std::size_t IndexType;

    /// Size definition
    typedef std::size_t SizeType;

    /// Pointer definition of MmgProcess
    KRATOS_CLASS_POINTER_DEFINITION(Parameters);

    /// Definiton of the iterators
    using iterator = iterator_adaptor;
    using const_iterator = const_iterator_adaptor;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    Parameters()
    {
    }

    /**
     * @brief String constructor. It takes a string as input, which parses into a nlohmann::json class
     * @param rJsonString The string to be parsed into a nlohmann::json class
     */
    Parameters(const std::string& json_string)
    {
        mpRoot = Kratos::shared_ptr<json>(new json( json::parse( json_string )));
        mpValue = mpRoot.get();
    }

    /// Copy constructor.
    Parameters(Parameters const& rOther)
    {
        //TODO: verify if mpValue is not null and eventually destruct correctly the data
        mpRoot = rOther.mpRoot;
        mpValue = rOther.mpValue;
    }

    /// Destructor.
    virtual ~Parameters()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    Parameters& operator=(Parameters const& rOther)
    {
        if(mpRoot.get() ==  mpValue || mpRoot == nullptr) {
            mpRoot = Kratos::shared_ptr<json>(new json( json::parse( rOther.WriteJsonString() )));
            mpValue = mpRoot.get();
        } else {
            *mpValue = json( json::parse( rOther.WriteJsonString() ) );
            // note that mpRoot is unchanged
        }

        return *this;
    }

    Parameters operator[](const std::string& entry)
    {
        return this->GetValue(entry);
    }

    Parameters operator[](unsigned int index)
    {
        return this->GetArrayItem(index);
    }

    ///@}
    ///@name Operations
    ///@{

    //generates a clone of the current document
    Parameters Clone()
    {
        //TODO: make a clone
        //TODO: find a better way to make the copy
        return Parameters(mpValue->dump());                     //new json(*mpValue));
    }

    const std::string WriteJsonString() const
    {
        return mpValue->dump();
    }

    const  std::string PrettyPrintJsonString() const
    {
        return mpValue->dump(4);
    }

    //*******************************************************************************************************
    Parameters GetValue(const std::string& entry)
    {
        auto j = mpValue->find(entry);
        KRATOS_ERROR_IF(j == mpValue->end()) << "Getting a value that does not exist. entry string : " << entry << std::endl;
        return Parameters(&(*j), mpRoot);
    }

    void SetValue(const std::string& entry, const Parameters& other_value)
    {
        KRATOS_ERROR_IF(mpValue->find(entry) == mpValue->end()) << "Value must exist to be set. Use AddValue instead" << std::endl;
        (*mpValue)[entry] = *(other_value.mpValue);
    }

    void AddValue(const std::string& entry, const Parameters& other_value)
    {
        if(mpValue->find(entry) == mpValue->end()) {
            (*mpValue)[entry] = *(other_value.mpValue);
        }
    }

    Parameters AddEmptyValue(const std::string& entry)
    {
        if(this->Has(entry) == false) {
            return Parameters(&(*mpValue)[entry],  mpRoot);
        }
        return this->GetValue(entry);
    }

    //*******************************************************************************************************

    bool RemoveValue(const std::string& entry)
    {
        // TODO: Implement this!!
        return false;
    }

    bool Has(const std::string& entry) const
    {
        return mpValue->find(entry) != mpValue->end();
    }

    bool IsNull() const
    {
        return mpValue->is_null();
    }
    bool IsNumber() const
    {
        return mpValue->is_number();
    }
    bool IsDouble() const
    {
        return mpValue->is_number_float();
    }
    bool IsInt() const
    {
        return mpValue->is_number_integer();
    }
    bool IsBool() const
    {
        return mpValue->is_boolean();
    }
    bool IsString() const
    {
        return mpValue->is_string();
    }
    bool IsArray() const
    {
        return mpValue->is_array();
    }
    bool IsVector() const
    {
        return false; //WIP
//         return mpValue->is_array();
    }
    bool IsMatrix() const
    {
        return false; //WIP
//         return mpValue->is_array();
    }
    bool IsSubParameter() const
    {
        return mpValue->is_object();
    }

    double GetDouble() const
    {
        return mpValue->get<double>();
    }
    int GetInt() const
    {
        KRATOS_ERROR_IF_NOT(mpValue->is_number()) << "Argument must be a number" << std::endl;
        return mpValue->get<int>();
    }
    bool GetBool() const
    {
        if (mpValue->is_boolean() == false) {
            //RecursivelyFindValue(*mpdoc, *mpValue);
            KRATOS_ERROR << "Argument must be a bool" << std::endl;
        }
        return mpValue->get<bool>();
    }
    std::string GetString() const
    {
        KRATOS_ERROR_IF_NOT(mpValue->is_string()) << "Argument must be a string" << std::endl;
        return mpValue->get<std::string>();
    }
    Vector GetVector() const
    {
//         if(mpValue->is_string() == false) KRATOS_THROW_ERROR(std::invalid_argument,"argument must be a string","");
//         return mpValue->get<std::string>();
        return Vector(); // WIP
    }
    Matrix GetMatrix() const
    {
//         if(mpValue->is_string() == false) KRATOS_THROW_ERROR(std::invalid_argument,"argument must be a string","");
//         return mpValue->get<std::string>();
        return Matrix(); // WIP
    }

    void SetDouble(const double value)
    {
        *mpValue=value;
    }
    void SetInt(const int value)
    {
        *mpValue=value;
    }
    void SetBool(const bool value)
    {
        *mpValue=value;
    }
    void SetString(const std::string& value)
    {
        *mpValue=value;
    }
    void SetVector(const Vector& value)
    {
        // TODO: Finish this
//         *mpValue=value;
    }
    void SetMatrix(const Matrix& value)
    {
        // TODO: Finish this
//         *mpValue=value;
    }

    iterator begin()
    {
        return iterator(mpValue->begin(),  mpRoot);
    }

    iterator end()
    {
        return iterator(mpValue->end(),  mpRoot);
    }

    const_iterator begin() const
    {
        return const_iterator(mpValue->cbegin(),  mpRoot);
    }

    const_iterator end() const
    {
        return const_iterator(mpValue->cend(),  mpRoot);
    }

    //*******************************************************************************************************
    //methods for array
    SizeType size() const
    {
        KRATOS_ERROR_IF_NOT(mpValue->is_array())  << "Size can only be queried if the value if of Array type" << std::endl;
        return mpValue->size();
    }

    void swap(Parameters& rOther) noexcept
    {
        std::swap(mpValue, rOther.mpValue);
        std::swap(mpRoot, rOther.mpRoot);
    }

    void Reset() noexcept
    {
        Parameters p;
        swap(p);
    }

    Parameters GetArrayItem(unsigned int index)
    {
        if(mpValue->is_array() == false)
            KRATOS_ERROR << "GetArrayItem only makes sense if the value if of Array type" << std::endl;
        else {
            KRATOS_ERROR_IF(index >= mpValue->size()) << "Index exceeds array size. Index value is : " << index << std::endl;
            return Parameters(&((*mpValue)[index]),  mpRoot);
        }
    }

    void SetArrayItem(unsigned int index, const Parameters& other_array_item)
    {
        if(mpValue->is_array() == false) {
            KRATOS_ERROR << "SetArrayItem only makes sense if the value if of Array type" << std::endl;
        } else {
            KRATOS_ERROR_IF(index >= mpValue->size()) << "Index exceeds array size. Index value is : " << index <<std::endl;
            (*mpValue)[index] = *other_array_item.mpValue;
        }
    }

    void AddEmptyArray(const std::string& rEntry)
    {
        // TODO: Implement this
    }

    void Append(const double Value)
    {
        // TODO: Implement this
    }
    void Append(const int Value)
    {
        // TODO: Implement this
    }

    void Append(const bool Value)
    {
        // TODO: Implement this
    }

    void Append(const std::string& rValue)
    {
        // TODO: Implement this
    }

    void Append(const Vector& rValue)
    {
        // TODO: Implement this
    }

    void Append(const Matrix& rValue)
    {
        // TODO: Implement this
    }

    void Append(const Parameters& rValue)
    {
        // TODO: Implement this
    }

    /**
     * @brief Checks if the names and values are the same, no importance to the order.
     * @details Lists have to be ordered, though! Take into account that in Kratos some physical vectors are represented with a list.
     * @param rParameters The parameters to be checked
     * @return True if it has, false othersise
     */
    bool IsEquivalentTo(Parameters& rParameters)
    {
        // TODO: Implement this
        return false;
    }

    /**
     * @brief  Checks if the names and the type of values are the same, no importance to the order.
     * @details Lists have to be ordered, though! Take into account that in Kratos some physical vectors are represented with a list.
     * @param rParameters The parameters to be checked
     * @return True if it has, false othersise
     */
    bool HasSameKeysAndTypeOfValuesAs(Parameters& rParameters)
    {
        // TODO: Implement this
        return false;
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
        for (auto itr = this->mpValue->begin(); itr != this->mpValue->end(); ++itr) {
            std::string item_name = itr.key();
            if(!defaults.Has(item_name) ) {
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
            auto value_defaults = (defaults[item_name]).GetUnderlyingStorage();
            if(itr->is_number() && value_defaults->is_number()) type_coincides = true;

//             if(itr->is_number_float() && value_defaults->is_number_float()) type_coincides = true;
            if(itr->is_array() && value_defaults->is_array()) type_coincides = true;
            if(itr->is_string() && value_defaults->is_string()) type_coincides = true;
            if(itr->is_object() && value_defaults->is_object()) type_coincides = true;

            //
            //both must be bool to be acceptable
            if(itr->is_boolean() && value_defaults->is_boolean()) type_coincides = true;

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

            //now iterate over all the defaults. In the case a default value is not assigned in the current Parameters
            //add an item copying its value
            if(defaults.IsSubParameter())
            {
                for (json::iterator itr = defaults.mpValue->begin(); itr != defaults.mpValue->end(); ++itr)
                {
                    std::string item_name = itr.key();
                    if(!this->Has(item_name))
                    {
                        (*mpValue)[item_name] = itr.value();
                    }
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
        for (auto itr = this->mpValue->cbegin(); itr != this->mpValue->cend(); ++itr)
        {
            std::string item_name = itr.key();

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
            auto value_defaults = (defaults[item_name]).GetUnderlyingStorage();
            if(itr->is_number_integer() && value_defaults->is_number_integer()) type_coincides = true;
            if(itr->is_boolean() && value_defaults->is_boolean()) type_coincides = true;
            if(itr->is_number_float() && value_defaults->is_number_float()) type_coincides = true;
            if(itr->is_array() && value_defaults->is_array()) type_coincides = true;
            if(itr->is_string() && value_defaults->is_string()) type_coincides = true;
            if(itr->is_object() && value_defaults->is_object()) type_coincides = true;

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
            if(itr->is_object())
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
            for (auto itr = defaults.mpValue->begin(); itr != defaults.mpValue->end(); ++itr)
            {
                std::string item_name = itr.key();
                if(mpValue->find(item_name) ==  mpValue->end())
                {
                    (*mpValue)[item_name] = itr.value();
                }

                //now walk the tree recursively
                if(itr->is_object())
                {
                    Parameters subobject = (*this)[item_name];
                    Parameters defaults_subobject = defaults[item_name];
                    subobject.ValidateAndAssignDefaults(defaults_subobject);
                }
            }
        }


        KRATOS_CATCH("")
    }

//     void RecursivelyFindValue(
//         const rapidjson::Value& rbase_value,
//         const rapidjson::Value& rvalue_to_find) const
//     {
//         for (rapidjson::Value::ConstMemberIterator itr = rbase_value.MemberBegin(); itr != rbase_value.MemberEnd(); ++itr)
//         {
//             if (&(itr->value) == &rvalue_to_find)
//             {
//                 rapidjson::StringBuffer buffer;
//                 rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(buffer);
//                 mpValue->Accept(writer);
//                 std::cout << "base = " << buffer.GetString() << std::endl;
//                 std::cout << "problematic var name " << itr->name.GetString() << " value " << itr->value.GetString() << std::endl;
//             }
//             else
//             {
//                 if (itr->value.IsObject()) RecursivelyFindValue(itr->value, rvalue_to_find);
//                 //TODO: it could be an array
//             }
//         }
//     }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

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

    json* mpValue;                   // This is where the json is actually stored
    Kratos::shared_ptr<json> mpRoot; // This is a shared pointer to the root structure (this is what allows us to acces in a tree structure to the JSON database)

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Direct constructor. It takes as parameters the "member" variables of the Parameters class
     * @param pValue The nlohmann::json class raw pointer
     * @param pRoot A shared pointer to a nlohmann::json class
     * @warning Please DO NOT use this constructor. It assumes nlohmann::json and hence it should be considered as an implementation detail
     */
    Parameters(json* pvalue, Kratos::shared_ptr<json> proot): mpValue(pvalue), mpRoot(proot)
    {}

    //ATTENTION: please DO NOT use this constructor. It assumes rapidjson and hence it should be considered as an implementation detail
//     Parameters(const json::iterator& it): mpValue(*it)
//     {
//         mis_owner = false;
//     }
//     //ATTENTION: please DO NOT use this constructor. It assumes rapidjson and hence it should be considered as an implementation detail
//     Parameters(const json::const_iterator& it): mpValue(*it)
//     {
//         mis_owner = false;
//     }

    /**
     * @brief This method is created in order to access from the iterators to the databae
     * @return mpValue The database storage
     * @warning Please DO NOT use this method. It is a low level accessor, and may change in the future
     */
    json* GetUnderlyingStorage()
    {
        return mpValue;
    }

    /**
     * @brief This method is created in order to access from the iterators to the databae
     * @return mpValue The database storage
     * @warning Please DO NOT use this method. It is a low level accessor, and may change in the future
     */
    Kratos::shared_ptr<json> GetUnderlyingRootStorage()
    {
        return mpRoot;
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
