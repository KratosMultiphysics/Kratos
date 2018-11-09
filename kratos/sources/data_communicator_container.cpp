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

#include "containers/data_communicator_container.h"
#include "includes/kratos_components.h"
#include "input_output/logger.h"

namespace Kratos {

DataCommunicatorContainer::DataCommunicatorContainer()
{}

DataCommunicatorContainer::~DataCommunicatorContainer()
{}

void DataCommunicatorContainer::Register(const std::string Name, DataCommunicator& rPrototype)
{
    auto found = mpDataCommunicators.find(Name);
    if (found == mpDataCommunicators.end())
    {
        auto result = mpDataCommunicators.emplace(Name, rPrototype.Clone());
        // result.first returns the created pair, pair_iterator->second the cloned prototype (which is a UniquePointer)
        auto pair_iterator = result.first;
        KratosComponents<DataCommunicator>::Add(Name, *(pair_iterator->second));
    }
    else {
        KRATOS_WARNING("DataCommunicatorContainer")
        << "Trying to register a new DataCommunicator with name " << Name
        << " but a DataCommunicator with the same name already exists: "
        << *(found->second)
        << " The second object has not been added." << std::endl;
    }
}

DataCommunicator& DataCommunicatorContainer::Get(const std::string Name) const
{
    auto found = mpDataCommunicators.find(Name);
    KRATOS_ERROR_IF(found == mpDataCommunicators.end()) << "Requesting unknown DataCommunicator " << Name << "." << std::endl;
    return *(found->second);
}

bool DataCommunicatorContainer::Has(const std::string Name) const
{
    return (mpDataCommunicators.find(Name) != mpDataCommunicators.end());
}


std::string DataCommunicatorContainer::Info() const
{
    std::stringstream buffer;
    PrintInfo(buffer);
    return buffer.str();
}

void DataCommunicatorContainer::PrintInfo(std::ostream &rOStream) const
{
    rOStream << "DataCommunicatorContainer";
}

void DataCommunicatorContainer::PrintData(std::ostream &rOStream) const
{
    rOStream << "Number of DataCommunicators: " << mpDataCommunicators.size() << std::endl;
    for (auto it_prototype = mpDataCommunicators.begin(); it_prototype != mpDataCommunicators.end(); ++it_prototype)
    {
        rOStream << "  \"" <<  it_prototype->first << "\": " << *(it_prototype->second);
    }
}

}