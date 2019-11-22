// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//                   license: CoSimulationApplication/license.txt
//
//  Main authors:    Philipp Bucher
//

#ifndef KRATOS_CO_SIM_IO_INTERNALS_H_INCLUDED
#define KRATOS_CO_SIM_IO_INTERNALS_H_INCLUDED

// System includes
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>

// Project includes
#include "co_sim_io_define.h"

namespace CoSimIO {
namespace Internals {

template<typename TDataType>
class DataContainer
{
public:
    virtual ~DataContainer() = default;

    virtual std::size_t size() const = 0;
    virtual void resize(const std::size_t NewSize) = 0;

    virtual const TDataType* data() const  = 0;
    TDataType* data()
    {
        return const_cast<TDataType*>(const_cast<const DataContainer*>(this)->data());
    }

    const TDataType& operator[](const std::size_t Index) const
    {
        return this->data()[Index];
    }

    TDataType& operator[](const std::size_t Index)
    {
        return const_cast<TDataType&>(const_cast<const DataContainer*>(this)->operator[](Index));
    }

    void resize_if_smaller(const std::size_t MinRequiredSize)
    {
        if (MinRequiredSize > this->size()) {
            this-resize(MinRequiredSize);
        }
    }
};


/// output stream function
template<typename TDataType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const DataContainer<TDataType>& rThis)
{
    // rThis.PrintInfo(rOStream);
    // rOStream << std::endl;
    // rThis.PrintData(rOStream);

    return rOStream;
}

template<typename TDataType>
class DataContainerStdVector : public DataContainer<TDataType>
{
public:
    DataContainerStdVector(std::vector<TDataType>& rVector) : mrVector(rVector) {}

    std::size_t size() const override {return mrVector.size();}
    void resize(const std::size_t NewSize) override {mrVector.resize(NewSize);}
    const TDataType* data() const override {return mrVector.data();}

private:
    std::vector<TDataType>& mrVector;
};

template<typename TDataType>
class DataContainerRawMemory : public DataContainer<TDataType>
{
public:
    DataContainerRawMemory(TDataType** ppData, const std::size_t Size) : mppData(ppData), mSize(Size) {}

    std::size_t size() const override {return mSize;};
    void resize(const std::size_t NewSize) override
    {
        // TODO use C-functions (malloc & free, or realloc)
        std::cout << "Before Resize" << std::endl;
        mSize = NewSize;
        delete [] mppData[0]; // delete the old memory and allocate new with different size
        std::cout << "Before New allocation" << std::endl;
        *mppData = new TDataType[mSize]; // TODO is this correct? shouldn't be "*mppData[0] = new TDataType[mSize]"?
        std::cout << "After New allocation" << std::endl;
        KRATOS_CO_SIM_ERROR_IF_NOT(*mppData) << "Memory reallocation failed";

        std::cout << "Exiting ..." << std::endl;

    };
    const TDataType* data() const override {return *mppData;}

private:
    TDataType** mppData;
    std::size_t mSize;
};

enum class ControlSignal
{
    Dummy,
    BreakSolutionLoop,
    ConvergenceAchieved,

    AdvanceInTime,
    InitializeSolutionStep,
    SolveSolutionStep,
    FinalizeSolutionStep,

    ImportGeometry,
    ExportGeometry,
    ImportMesh,
    ExportMesh,
    ImportData,
    ExportData,
};

typedef std::unordered_map<std::string, std::string> SettingsType;

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
        std::cout << "Input file \"" << rSettingsFileName << "\" could not be read, using default configuration" << std::endl;
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

#endif /* KRATOS_CO_SIM_IO_INTERNALS_H_INCLUDED */
