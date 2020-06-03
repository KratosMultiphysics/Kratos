//     ______     _____ _           ________
//    / ____/___ / ___/(_)___ ___  /  _/ __ |
//   / /   / __ \\__ \/ / __ `__ \ / // / / /
//  / /___/ /_/ /__/ / / / / / / // // /_/ /
//  \____/\____/____/_/_/ /_/ /_/___/\____/
//  Kratos CoSimulationApplication
//
//  License:         BSD License, see license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

#ifndef CO_SIM_IO_INFO_INCLUDED
#define CO_SIM_IO_INFO_INCLUDED

// System includes
#include <string>
#include <map>
#include <memory>
#include <iostream>
#include <type_traits>

// Project includes
#include "macros.hpp"

namespace CoSimIO {
namespace Internals {

inline std::string Name(int dummy)         {return "int";}
inline std::string Name(double dummy)      {return "double";}
inline std::string Name(bool dummy)        {return "bool";}
inline std::string Name(std::string dummy) {return "string";}

class InfoDataBase
{
public:
    virtual ~InfoDataBase() = default;
    virtual const void* GetData() const = 0;
    virtual std::string GetDataTypeName() const = 0;
    virtual std::shared_ptr<InfoDataBase> Clone() const = 0;
    // virtual void Print(const void* pSource, std::ostream& rOStream) const;
    virtual void Save(std::ostream& O_OutStream) const = 0;
    virtual void Load(std::istream& I_InStream) = 0;
    virtual void Print(std::ostream& rOStream) const = 0;
};

/// output stream function
inline std::ostream & operator <<(
    std::ostream& rOStream,
    const InfoDataBase& rThis)
{
    rThis.Print(rOStream);

    return rOStream;
}

template<class TDataType>
class InfoData : public InfoDataBase
{
public:
    explicit InfoData(const TDataType Source) : mData(Source) {}

    std::string GetDataTypeName() const override {return Internals::Name(TDataType());}

    const void* GetData() const override
    {
        return &mData;
    }

    std::shared_ptr<InfoDataBase> Clone() const override
    {
        return std::make_shared<InfoData<TDataType>>(TDataType());
    }

    void Save(std::ostream& O_OutStream) const override
    {
        O_OutStream << mData;
    }
    void Load(std::istream& I_InStream) override
    {
        I_InStream >> mData;
    }

    void Print(std::ostream& rOStream) const override
    {
        rOStream << "value: " << mData << " | type: " << GetDataTypeName();
    }


private:
    TDataType mData;
};

} // namespace Internals


class Info
{
public:
    Info() : mOptions(){}

    Info(const Info& Other) : mOptions(Other.mOptions){}

    Info& operator=(const Info& Other){
        mOptions = Other.mOptions;
        return *this;
    }

    virtual ~Info() = default;

    template<typename TDataType>
    const TDataType& Get(const std::string& I_Key) const
    {
        CO_SIM_IO_ERROR_IF_NOT(Has(I_Key)) << "Trying to get \"" << I_Key << "\" which does not exist!" << std::endl;
        const auto& r_val = mOptions.at(I_Key);
        CO_SIM_IO_ERROR_IF(r_val->GetDataTypeName() != Internals::Name(TDataType())) << "Wrong DataType! Trying to get \"" << I_Key << "\" which is of type \"" << r_val->GetDataTypeName() << "\" with \"" << Internals::Name(TDataType()) << "\"!" << std::endl;
        return *static_cast<const TDataType*>(r_val->GetData());
    }

    bool Has(const std::string& I_Key) const
    {
        return mOptions.count(I_Key)>0;
    }

    template<typename TDataType>
    void Set(const std::string& I_Key, TDataType I_Value)
    {
        static_assert(
            std::is_same<TDataType, double>::value ||
            std::is_same<TDataType, int>::value    ||
            std::is_same<TDataType, bool>::value   ||
            // std::is_same<TDataType, Info>::value   || // makes it recursive
            std::is_same<TDataType, std::string>::value,
                "Only allowed types are double, int, bool, string");

        mOptions[I_Key] = std::make_shared<Internals::InfoData<TDataType>>(I_Value);
    }

    void Set(const std::string& I_Key, const char * I_Value)
    {
       Set(I_Key, std::string(I_Value));
    }

    void Erase(const std::string& I_Key)
    {
        mOptions.erase(I_Key);
    }

    void Clear()
    {
        mOptions.clear();
    }

    std::size_t Size() const
    {
        return mOptions.size();
    }

    void Save(std::ostream& O_OutStream) const
    {
        static std::map<std::string, std::string> s_registered_object_names {
            {typeid(Internals::InfoData<int>).name(),         "InfoData_int"},
            {typeid(Internals::InfoData<double>).name(),      "InfoData_double"},
            {typeid(Internals::InfoData<bool>).name(),        "InfoData_bool"},
            {typeid(Internals::InfoData<std::string>).name(), "InfoData_string"}
        };

        O_OutStream << Size() << "\n";
        for (const auto& r_pair: mOptions) {
            const auto& r_val = *(r_pair.second);
            auto it_obj = s_registered_object_names.find(typeid(r_val).name());
            CO_SIM_IO_ERROR_IF(it_obj == s_registered_object_names.end()) << "No name registered" << std::endl;
            O_OutStream << r_pair.first << "\n";
            O_OutStream << it_obj->second << "\n";
            r_pair.second->Save(O_OutStream);
            O_OutStream << "\n";
        }
    }
    void Load(std::istream& I_InStream)
    {
        static std::map<std::string, std::shared_ptr<Internals::InfoDataBase>> s_registered_object_prototypes {
            {"InfoData_int"    , std::make_shared<Internals::InfoData<int>>(1)},
            {"InfoData_double" , std::make_shared<Internals::InfoData<double>>(1)},
            {"InfoData_bool"   , std::make_shared<Internals::InfoData<bool>>(1)},
            {"InfoData_string" , std::make_shared<Internals::InfoData<std::string>>("")}
        };

        std::string key, registered_name;

        int size;
        I_InStream >> size;

        for (int i=0; i<size; ++i) {
            I_InStream >> key;
            I_InStream >> registered_name;
            auto it_prototype = s_registered_object_prototypes.find(registered_name);
            CO_SIM_IO_ERROR_IF(it_prototype == s_registered_object_prototypes.end()) << "No prototype registered for " << registered_name << std::endl;

            auto p_clone = it_prototype->second->Clone();
            p_clone->Load(I_InStream);
            mOptions[key] = p_clone;
        }
    }

    virtual void Print(std::ostream& rOStream) const
    {
        rOStream << "CoSimIO-Info; containing " << Size() << " entries\n";

        for (const auto& r_pair: mOptions) {
            rOStream << "  name: " << r_pair.first << " | " << *(r_pair.second) << std::endl;
        }
    }

private:
    std::map<std::string, std::shared_ptr<Internals::InfoDataBase>> mOptions;
};

/// output stream function
inline std::ostream & operator <<(
    std::ostream& rOStream,
    const Info& rThis)
{
    rThis.Print(rOStream);
    return rOStream;
}

} // namespace CoSimIO

#endif // CO_SIM_IO_INFO_INCLUDED
