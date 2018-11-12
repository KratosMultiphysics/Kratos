//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Jordi Cotela
//

#include "includes/parallel_environment.h"
#include "includes/kratos_components.h"
#include "input_output/logger.h"

namespace Kratos {

ParallelEnvironment::ParallelEnvironment()
{
    RegisterDataCommunicator("Serial", DataCommunicator(), MakeDefault);
}

ParallelEnvironment& ParallelEnvironment::GetInstance()
{
    static ParallelEnvironment parallel_environment;

    return parallel_environment;
}

void ParallelEnvironment::RegisterDataCommunicator(
    const std::string Name,
    const DataCommunicator& rPrototype,
    bool Default)
{
    auto found = mDataCommunicators.find(Name);
    if (found == mDataCommunicators.end())
    {
        auto result = mDataCommunicators.emplace(Name, rPrototype.Clone());
        // result.first returns the created pair, pair_iterator->second the cloned prototype (which is a UniquePointer)
        auto pair_iterator = result.first;
        KratosComponents<DataCommunicator>::Add(Name, *(pair_iterator->second));

        if (Default == MakeDefault)
        {
            mDefaultCommunicator = pair_iterator;
        }
    }
    else {
        KRATOS_WARNING("ParallelEnvironment")
        << "Trying to register a new DataCommunicator with name " << Name
        << " but a DataCommunicator with the same name already exists: "
        << *(found->second)
        << " The second object has not been added." << std::endl;
    }
}

DataCommunicator& ParallelEnvironment::GetDataCommunicator(const std::string Name) const
{
    auto found = mDataCommunicators.find(Name);
    KRATOS_ERROR_IF(found == mDataCommunicators.end()) << "Requesting unknown DataCommunicator " << Name << "." << std::endl;
    return *(found->second);
}

DataCommunicator& ParallelEnvironment::GetDefaultDataCommunicator() const
{
    return *(mDefaultCommunicator->second);
}

void ParallelEnvironment::SetDefaultDataCommunicator(const std::string Name)
{
    auto found = mDataCommunicators.find(Name);
    KRATOS_ERROR_IF(found == mDataCommunicators.end())
    << "Trying to set \"" << Name << "\" as the default DataCommunicator,"
    << " but no registered DataCommunicator with that name has been found." << std::endl;

    mDefaultCommunicator = found;
}

bool ParallelEnvironment::HasDataCommunicator(const std::string Name) const
{
    return (mDataCommunicators.find(Name) != mDataCommunicators.end());
}


std::string ParallelEnvironment::Info() const
{
    std::stringstream buffer;
    PrintInfo(buffer);
    return buffer.str();
}

void ParallelEnvironment::PrintInfo(std::ostream &rOStream) const
{
    rOStream << "ParallelEnvironment";
}

void ParallelEnvironment::PrintData(std::ostream &rOStream) const
{
    rOStream << "Number of DataCommunicators: " << mDataCommunicators.size() << std::endl;
    for (auto it_prototype = mDataCommunicators.begin(); it_prototype != mDataCommunicators.end(); ++it_prototype)
    {
        rOStream << "  \"" <<  it_prototype->first << "\": " << *(it_prototype->second);
    }
    rOStream << "Default communicator: \"" << mDefaultCommunicator->first << "\": " << *(mDefaultCommunicator->second);
}

}