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

#ifndef CO_SIM_IO_INTERNALS_H_INCLUDED
#define CO_SIM_IO_INTERNALS_H_INCLUDED

// System includes
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>

// Project includes
#include "define.hpp"
#include "macros.hpp"

namespace CoSimIO {
namespace Internals {

template<typename TDataType>
class DataContainer
{
public:
    virtual ~DataContainer() = default;

    virtual std::size_t size() const = 0;
    virtual void resize(const std::size_t NewSize) = 0;

    virtual TDataType* data() = 0;
    virtual const TDataType* data() const = 0;
    const TDataType* data_const() const
    {
        const auto& r_const_this = *this;
        return this->data();
    }

    const TDataType& operator[](const std::size_t Index) const
    {
        return this->data()[Index];
    }

    TDataType& operator[](const std::size_t Index)
    {
        return this->data()[Index];
    }
};


/// output stream function
template<typename TDataType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const DataContainer<TDataType>& rThis)
{
    const std::size_t size = rThis.size();

    rOStream << "[";
    if(size>0) rOStream << rThis[0];
    if(size>1) {
        for(std::size_t i=1; i<size; ++i)
            rOStream<<", "<<rThis[i];
    }
    rOStream << "]";

    return rOStream;
}

template<typename TDataType>
class DataContainerStdVector : public DataContainer<TDataType>
{
public:
    explicit DataContainerStdVector(std::vector<TDataType>& rVector)
        : mrVector(rVector) {}

    std::size_t size() const override {return mrVector.size();}
    void resize(const std::size_t NewSize) override {mrVector.resize(NewSize);} // resize does not change the capacity if resized to a smaller size
    const TDataType* data() const override {return mrVector.data();}
    TDataType* data() override {return mrVector.data();}

private:
    std::vector<TDataType>& mrVector;
};

template<typename TDataType>
class DataContainerStdVectorReadOnly : public DataContainer<TDataType>
{
public:
    explicit DataContainerStdVectorReadOnly(const std::vector<TDataType>& rVector)
        : mrVector(rVector) {}

    std::size_t size() const override {return mrVector.size();}
    void resize(const std::size_t NewSize) override {CO_SIM_IO_ERROR << "Using non-const member of readonly object!" << std::endl;}
    const TDataType* data() const override {return mrVector.data();}
    TDataType* data() override {CO_SIM_IO_ERROR << "Using non-const member of readonly object!" << std::endl; return nullptr;}

private:
    const std::vector<TDataType>& mrVector;
};

template<typename TDataType>
class DataContainerRawMemory : public DataContainer<TDataType>
{
public:
    explicit DataContainerRawMemory(TDataType** ppData, const std::size_t Size)
        : mppData(ppData), mSize(Size) {}

    std::size_t size() const override {return mSize;};
    void resize(const std::size_t NewSize) override
    {
        if (NewSize > mSize) { // only increase the capacity if too small => same behavior as std::vector
            if(mSize != 0) // Maybe is not null neigther allocated: double *data. Pooyan.
                free(*mppData); // this is ok according to the standard, no matter if it is null or allocated //also check if using "std::"

            *mppData = (TDataType *)malloc((NewSize)*sizeof(TDataType)); // TODO maybe use realloc? //also check if using "std::"
            CO_SIM_IO_ERROR_IF_NOT(*mppData) << "Memory reallocation failed";
        }
        mSize = NewSize;

    };
    const TDataType* data() const override {return *mppData;}
    TDataType* data() override {return *mppData;}

private:
    TDataType** mppData;
    std::size_t mSize;
};

template<typename TDataType>
class DataContainerRawMemoryReadOnly : public DataContainer<TDataType>
{
public:
    explicit DataContainerRawMemoryReadOnly(const TDataType* pData, const std::size_t Size)
        : mpData(pData), mSize(Size) {}

    std::size_t size() const override {return mSize;};
    void resize(const std::size_t NewSize) override {CO_SIM_IO_ERROR << "Using non-const member of readonly object!" << std::endl;};
    const TDataType* data() const override {return mpData;}
    TDataType* data() override {CO_SIM_IO_ERROR << "Using non-const member of readonly object!" << std::endl; return nullptr;}

private:
    const TDataType* mpData;
    const std::size_t mSize;
};

inline void AddMissingSettings(const SettingsType& rDefaultSettings, SettingsType& rSettings)
{
    for (const auto& r_setting : rDefaultSettings) {
        if (rSettings.count(r_setting.first) == 0) {
            rSettings[r_setting.first] = r_setting.second;
        }
    }
}

inline SettingsType ReadSettingsFile(const std::string& rSettingsFileName)
{
    std::ifstream settings_file(rSettingsFileName);

    if (!settings_file.good()) {
        CO_SIM_IO_INFO("CoSimIO") << "Input file \"" << rSettingsFileName << "\" could not be read, using default configuration" << std::endl;
        return SettingsType();
    }

    SettingsType settings;
    std::string current_line;
    while (std::getline(settings_file, current_line)) {
        // TODO implement this
    }
    return settings;
}

} // namespace Internals
} // namespace CoSimIO

#endif // CO_SIM_IO_INTERNALS_H_INCLUDED
