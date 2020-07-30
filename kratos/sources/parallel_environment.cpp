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

#ifdef _OPENMP
#include "omp.h"
#endif

#include "includes/parallel_environment.h"
#include "includes/kratos_components.h"
#include "input_output/logger.h"

namespace Kratos {

// Public interface of ParallelEnvironment ////////////////////////////////////

DataCommunicator& ParallelEnvironment::GetDataCommunicator(const std::string& rName)
{
    const ParallelEnvironment& env = GetInstance();
    return env.GetDataCommunicatorDetail(rName);
}

DataCommunicator& ParallelEnvironment::GetDefaultDataCommunicator()
{
    const ParallelEnvironment& env = GetInstance();
    return env.GetDefaultDataCommunicatorDetail();
}

void ParallelEnvironment::SetDefaultDataCommunicator(const std::string& rName)
{
    ParallelEnvironment& env = GetInstance();
    env.SetDefaultDataCommunicatorDetail(rName);
}

int ParallelEnvironment::GetDefaultRank()
{
    const ParallelEnvironment& env = GetInstance();
    return env.mDefaultRank;
}

int ParallelEnvironment::GetDefaultSize()
{
    const ParallelEnvironment& env = GetInstance();
    return env.mDefaultSize;
}

void ParallelEnvironment::SetUpMPIEnvironment(EnvironmentManager::Pointer pEnvironmentManager)
{
    ParallelEnvironment& env = GetInstance();
    env.SetUpMPIEnvironmentDetail(std::move(pEnvironmentManager));
}

void ParallelEnvironment::RegisterDataCommunicator(
    const std::string& Name,
    DataCommunicator::UniquePointer pPrototype,
    const bool Default)
{
    ParallelEnvironment& env = GetInstance();
    env.RegisterDataCommunicatorDetail(Name, std::move(pPrototype), Default);
}

void ParallelEnvironment::UnregisterDataCommunicator(const std::string& Name)
{
    ParallelEnvironment& env = GetInstance();
    env.UnregisterDataCommunicatorDetail(Name);
}

bool ParallelEnvironment::HasDataCommunicator(const std::string& rName)
{
    const ParallelEnvironment& env = GetInstance();
    return env.HasDataCommunicatorDetail(rName);
}

std::string ParallelEnvironment::GetDefaultDataCommunicatorName()
{
    const ParallelEnvironment& env = GetInstance();
    return env.mDefaultCommunicator->first;
}

bool ParallelEnvironment::MPIIsInitialized()
{
    const ParallelEnvironment& env = GetInstance();
    return env.MPIIsInitializedDetail();
}

bool ParallelEnvironment::MPIIsFinalized()
{
    const ParallelEnvironment& env = GetInstance();
    return env.MPIIsFinalizedDetail();
}

std::string ParallelEnvironment::Info()
{
    const ParallelEnvironment& env = GetInstance();
    return env.InfoDetail();
}

void ParallelEnvironment::PrintInfo(std::ostream &rOStream)
{
    const ParallelEnvironment& env = GetInstance();
    return env.PrintInfoDetail(rOStream);
}

void ParallelEnvironment::PrintData(std::ostream &rOStream)
{
    const ParallelEnvironment& env = GetInstance();
    return env.PrintDataDetail(rOStream);
}

// Implementation details /////////////////////////////////////////////////////

ParallelEnvironment::ParallelEnvironment()
{
    RegisterDataCommunicatorDetail("Serial", DataCommunicator::Create(), MakeDefault);
}

ParallelEnvironment::~ParallelEnvironment()
{
    // First release the registered DataCommunicators
    mDataCommunicators.clear();

    // Then finalize MPI if necessary by freeing the manager instance
    if (mpEnvironmentManager)
    {
        mpEnvironmentManager.reset();
    }    

    mDestroyed = true;
    mpInstance = nullptr;
}

ParallelEnvironment& ParallelEnvironment::GetInstance()
{
    // Using double-checked locking to ensure thread safety in the first creation of the singleton.
    if (mpInstance == nullptr)
    {
        #ifdef _OPENMP
        #pragma omp critical
        if (mpInstance == nullptr)
        {
        #endif
            KRATOS_ERROR_IF(mDestroyed) << "Accessing ParallelEnvironment after its destruction" << std::endl;
            Create();
        #ifdef _OPENMP
        }
        #endif
    }

    return *mpInstance;
}

void ParallelEnvironment::SetAsDefault(
    std::unordered_map<std::string, DataCommunicator::UniquePointer>::iterator& rThisCommunicator)
{
    mDefaultCommunicator = rThisCommunicator;
    const auto& r_comm = *(rThisCommunicator->second);
    mDefaultRank = r_comm.Rank();
    mDefaultSize = r_comm.Size();
}

void ParallelEnvironment::Create()
{
    static ParallelEnvironment parallel_environment;
    mpInstance = &parallel_environment;
}

void ParallelEnvironment::SetUpMPIEnvironmentDetail(EnvironmentManager::Pointer pEnvironmentManager)
{
    KRATOS_ERROR_IF(MPIIsInitialized() || MPIIsFinalized())
    << "Trying to configure run for MPI twice. This should not be happening!" << std::endl;

    mpEnvironmentManager = std::move(pEnvironmentManager);
}

void ParallelEnvironment::RegisterDataCommunicatorDetail(
    const std::string& Name,
    DataCommunicator::UniquePointer pPrototype,
    const bool Default)
{
    auto found = mDataCommunicators.find(Name);
    if (found == mDataCommunicators.end())
    {
        auto result = mDataCommunicators.emplace(Name, std::move(pPrototype));
        // result.first returns the created pair, pair_iterator->second the cloned prototype (which is a UniquePointer)
        auto pair_iterator = result.first;
        KratosComponents<DataCommunicator>::Add(Name, *(pair_iterator->second));

        if (Default == MakeDefault)
        {
            SetAsDefault(pair_iterator);
        }
    }
    else {
        KRATOS_WARNING("ParallelEnvironment")
        << "Trying to register a new DataCommunicator with name " << Name
        << " but a DataCommunicator with the same name already exists: "
        << *(found->second)
        << " The provided DataCommunicator has not been added." << std::endl;
    }
}

void ParallelEnvironment::UnregisterDataCommunicatorDetail(const std::string& Name)
{
    KRATOS_ERROR_IF(Name == mDefaultCommunicator->first)
    << "Trying to unregister the default DataCommunicator \"" << Name
    << "\". Please define a new default before unregistering the current one."
    << std::endl;
    int num_erased = mDataCommunicators.erase(Name);
    KRATOS_WARNING_IF("ParallelEnvironment", num_erased == 0)
    << "Trying to unregister a DataCommunicator with name " << Name
    << " but no DataCommunicator of that name exsits."
    << " No changes were made." << std::endl;
    if (num_erased == 1)
    {
        KratosComponents<DataCommunicator>::Remove(Name);
    }
}

DataCommunicator& ParallelEnvironment::GetDataCommunicatorDetail(const std::string& rName) const
{
    auto found = mDataCommunicators.find(rName);
    KRATOS_ERROR_IF(found == mDataCommunicators.end())
    << "Requesting unknown DataCommunicator " << rName << "." << std::endl;
    return *(found->second);
}

DataCommunicator& ParallelEnvironment::GetDefaultDataCommunicatorDetail() const
{
    return *(mDefaultCommunicator->second);
}

void ParallelEnvironment::SetDefaultDataCommunicatorDetail(const std::string& rName)
{
    auto found = mDataCommunicators.find(rName);
    KRATOS_ERROR_IF(found == mDataCommunicators.end())
    << "Trying to set \"" << rName << "\" as the default DataCommunicator,"
    << " but no registered DataCommunicator with that name has been found." << std::endl;

    SetAsDefault(found);
}

bool ParallelEnvironment::HasDataCommunicatorDetail(const std::string& rName) const
{
    return (mDataCommunicators.find(rName) != mDataCommunicators.end());
}

bool ParallelEnvironment::MPIIsInitializedDetail() const
{
    return (mpEnvironmentManager == nullptr) ? false : mpEnvironmentManager->IsInitialized();
}

bool ParallelEnvironment::MPIIsFinalizedDetail() const
{
    return (mpEnvironmentManager == nullptr) ? false : mpEnvironmentManager->IsFinalized();
}

std::string ParallelEnvironment::InfoDetail() const
{
    std::stringstream buffer;
    PrintInfo(buffer);
    return buffer.str();
}

void ParallelEnvironment::PrintInfoDetail(std::ostream &rOStream) const
{
    rOStream << "ParallelEnvironment";
}

void ParallelEnvironment::PrintDataDetail(std::ostream &rOStream) const
{
    rOStream << "Number of DataCommunicators: " << mDataCommunicators.size() << std::endl;
    for (auto it_prototype = mDataCommunicators.begin(); it_prototype != mDataCommunicators.end(); ++it_prototype)
    {
        rOStream << "  \"" <<  it_prototype->first << "\": " << *(it_prototype->second);
    }
    rOStream << "Default communicator: \"" << mDefaultCommunicator->first << "\": " << *(mDefaultCommunicator->second);
}

ParallelEnvironment* ParallelEnvironment::mpInstance = nullptr;
bool ParallelEnvironment::mDestroyed = false;

}