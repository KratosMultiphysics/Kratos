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

// Project includes
#include "define.hpp"
#include "serializer.hpp"

namespace CoSimIO {
namespace Internals {

template<typename T>
std::string Name();

class InfoDataBase
{
public:
    virtual ~InfoDataBase() = default;
    virtual const void* GetData() const {CO_SIM_IO_ERROR << "This is the baseclass!" << std::endl;};
    virtual std::string GetDataTypeName() const {CO_SIM_IO_ERROR << "This is the baseclass!" << std::endl;};
    virtual void Print(
        std::ostream& rOStream,
        const std::string& rPrefixString) const {CO_SIM_IO_ERROR << "This is the baseclass!" << std::endl;};

protected:
    InfoDataBase() = default; // needed for Serializer
private:
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const {};
    virtual void load(Serializer& rSerializer) {};
};


template<class TDataType>
class InfoData : public InfoDataBase
{
public:
    explicit InfoData(const TDataType Source) : mData(Source) {}

    std::string GetDataTypeName() const override {return Internals::Name<TDataType>();}

    const void* GetData() const override
    {
        return &mData;
    }

    void Print(
        std::ostream& rOStream,
        const std::string& rPrefixString) const override
    {
        rOStream << "value: " << mData << " | type: " << GetDataTypeName() << "\n";
    }

private:
    TDataType mData;

    InfoData() = default;

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        CO_SIM_IO_SERIALIZE_SAVE_BASE_CLASS(rSerializer, InfoDataBase)
        rSerializer.save("mData", mData);

    }
    void load(Serializer& rSerializer) override
    {
        CO_SIM_IO_SERIALIZE_LOAD_BASE_CLASS(rSerializer, InfoDataBase)
        rSerializer.load("mData", mData);
    }
};

} // namespace Internals


class CO_SIM_IO_API Info
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
        ValidateType<TDataType>();

        CO_SIM_IO_ERROR_IF_NOT(Has(I_Key)) << "Trying to get \"" << I_Key << "\" which does not exist!\nCurrently available:\n" << *this << std::endl;
        return GetExistingKey<TDataType>(I_Key);
    }

    template<typename TDataType>
    const TDataType& Get(const std::string& I_Key, const TDataType& I_Default) const
    {
        ValidateType<TDataType>();

        if (Has(I_Key)) {
            return GetExistingKey<TDataType>(I_Key);
        } else {
            // this does NOT insert the value! (same behavior as in python)
            return I_Default;
        }
    }

    bool Has(const std::string& I_Key) const
    {
        return mOptions.count(I_Key)>0;
    }

    template<typename TDataType>
    void Set(const std::string& I_Key, TDataType I_Value)
    {
        ValidateType<TDataType>();

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

    void Print(
        std::ostream& rOStream,
        const std::string& rPrefixString="") const;

private:
    std::map<std::string, std::shared_ptr<Internals::InfoDataBase>> mOptions;

    static bool* mpSerializerTypesRegistered;

    template<typename TDataType>
    const TDataType& GetExistingKey(const std::string& I_Key) const
    {
        ValidateType<TDataType>();

        const auto& r_val = mOptions.at(I_Key);
        CO_SIM_IO_ERROR_IF(r_val->GetDataTypeName() != Internals::Name<TDataType>()) << "Wrong DataType! Trying to get \"" << I_Key << "\" which is of type \"" << r_val->GetDataTypeName() << "\" with \"" << Internals::Name<TDataType>() << "\"!" << std::endl;
        return *static_cast<const TDataType*>(r_val->GetData());
    }

    template<typename TDataType>
    static void ValidateType()
    {
        static_assert(
            std::is_same<TDataType, double>::value      ||
            std::is_same<TDataType, int>::value         ||
            std::is_same<TDataType, std::size_t>::value ||
            std::is_same<TDataType, bool>::value        ||
            std::is_same<TDataType, Info>::value        || // makes it recursive
            std::is_same<TDataType, std::string>::value,
                "Only allowed types are double, int, size_t, bool, string, Info");
    }

    friend class CoSimIO::Internals::Serializer; // needs "CoSimIO::Internals::" because it is in different namespace

    void save(CoSimIO::Internals::Serializer& rSerializer) const;

    void load(CoSimIO::Internals::Serializer& rSerializer);

    static void RegisterTypesInSerializer();
};

/// output stream function
inline std::ostream & operator <<(
    std::ostream& rOStream,
    const Info& rThis)
{
    rThis.Print(rOStream);
    return rOStream;
}

namespace Internals {

template<> inline std::string Name<int>()         {return "int";}
template<> inline std::string Name<std::size_t>() {return "size_t";}
template<> inline std::string Name<double>()      {return "double";}
template<> inline std::string Name<bool>()        {return "bool";}
template<> inline std::string Name<std::string>() {return "string";}
template<> inline std::string Name<Info>()        {return "info";}

template<>
inline void InfoData<bool>::Print(
    std::ostream& rOStream,
    const std::string& rPrefixString) const
{
    rOStream << "value: " << std::boolalpha << mData << std::noboolalpha << " | type: " << GetDataTypeName() << "\n";
}

template<>
inline void InfoData<Info>::Print(
    std::ostream& rOStream,
    const std::string& rPrefixString) const
{
    rOStream << "type: ";
    mData.Print(rOStream, rPrefixString);
}

} // namespace Internals
} // namespace CoSimIO

#endif // CO_SIM_IO_INFO_INCLUDED
