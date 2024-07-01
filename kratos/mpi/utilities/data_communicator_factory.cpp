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

#include "data_communicator_factory.h"

#include "mpi.h"

#include "includes/parallel_environment.h"
#include "mpi/includes/mpi_data_communicator.h"

namespace Kratos
{
namespace DataCommunicatorFactory
{

const DataCommunicator& DuplicateAndRegister(
    const DataCommunicator& rOriginalCommunicator,
    const std::string& rNewCommunicatorName)
{
    MPI_Comm origin_mpi_comm = MPIDataCommunicator::GetMPICommunicator(rOriginalCommunicator);
    MPI_Comm duplicate_comm;
    MPI_Comm_dup(origin_mpi_comm, &duplicate_comm);

    ParallelEnvironment::RegisterDataCommunicator(
        rNewCommunicatorName, MPIDataCommunicator::Create(duplicate_comm), ParallelEnvironment::DoNotMakeDefault);
    return ParallelEnvironment::GetDataCommunicator(rNewCommunicatorName);
}

const DataCommunicator& SplitAndRegister(
    const DataCommunicator& rOriginalCommunicator,
    int Color, int Key,
    const std::string& rNewCommunicatorName)
{
    MPI_Comm origin_mpi_comm = MPIDataCommunicator::GetMPICommunicator(rOriginalCommunicator);
    MPI_Comm split_mpi_comm;
    MPI_Comm_split(origin_mpi_comm, Color, Key, &split_mpi_comm);

    ParallelEnvironment::RegisterDataCommunicator(
        rNewCommunicatorName, MPIDataCommunicator::Create(split_mpi_comm), ParallelEnvironment::DoNotMakeDefault);
    return ParallelEnvironment::GetDataCommunicator(rNewCommunicatorName);
}

const DataCommunicator& CreateFromRanksAndRegister(
        const DataCommunicator& rOriginalCommunicator,
        const std::vector<int>& rRanks,
        const std::string& rNewCommunicatorName)
{
    MPI_Comm origin_mpi_comm = MPIDataCommunicator::GetMPICommunicator(rOriginalCommunicator);
    MPI_Group all_ranks, selected_ranks;

    MPI_Comm_group(origin_mpi_comm, &all_ranks);
    MPI_Group_incl(all_ranks, rRanks.size(), rRanks.data(), &selected_ranks);

    MPI_Comm comm_from_ranks;
    int tag = 0;
    MPI_Comm_create_group(origin_mpi_comm, selected_ranks, tag, &comm_from_ranks);

    MPI_Group_free(&all_ranks);
    MPI_Group_free(&selected_ranks);
    ParallelEnvironment::RegisterDataCommunicator(
        rNewCommunicatorName, MPIDataCommunicator::Create(comm_from_ranks), ParallelEnvironment::DoNotMakeDefault);
    return ParallelEnvironment::GetDataCommunicator(rNewCommunicatorName);
}

const DataCommunicator& CreateUnionAndRegister(
    const DataCommunicator& rFirstDataCommunicator,
    const DataCommunicator& rSecondDataCommunicator,
    const DataCommunicator& rParentDataCommunicator,
    const std::string& rNewCommunicatorName)
{
    MPI_Comm parent_comm = MPIDataCommunicator::GetMPICommunicator(rParentDataCommunicator);
    MPI_Comm first_mpi_comm = MPIDataCommunicator::GetMPICommunicator(rFirstDataCommunicator);
    MPI_Comm second_mpi_comm = MPIDataCommunicator::GetMPICommunicator(rSecondDataCommunicator);

    MPI_Comm combined_mpi_comm;
    constexpr int key = 0; // With key = 0 we are preserving the rank ordering on the new communicator

    if ( (first_mpi_comm != MPI_COMM_NULL) || (second_mpi_comm != MPI_COMM_NULL) )
    {
        constexpr int color = 0;
        MPI_Comm_split(parent_comm, color, key, &combined_mpi_comm);
    }
    else
    {
        MPI_Comm_split(parent_comm, MPI_UNDEFINED, key, &combined_mpi_comm);
    }

    ParallelEnvironment::RegisterDataCommunicator(
        rNewCommunicatorName, MPIDataCommunicator::Create(combined_mpi_comm), ParallelEnvironment::DoNotMakeDefault);
    return ParallelEnvironment::GetDataCommunicator(rNewCommunicatorName);
}

const DataCommunicator& CreateIntersectionAndRegister(
    const DataCommunicator& rFirstDataCommunicator,
    const DataCommunicator& rSecondDataCommunicator,
    const DataCommunicator& rParentDataCommunicator,
    const std::string& rNewCommunicatorName)
{
    MPI_Comm parent_comm = MPIDataCommunicator::GetMPICommunicator(rParentDataCommunicator);
    MPI_Comm first_mpi_comm = MPIDataCommunicator::GetMPICommunicator(rFirstDataCommunicator);
    MPI_Comm second_mpi_comm = MPIDataCommunicator::GetMPICommunicator(rSecondDataCommunicator);

    MPI_Comm combined_mpi_comm;
    constexpr int key = 0; // With key = 0 we are preserving the rank ordering on the new communicator

    if ( (first_mpi_comm != MPI_COMM_NULL) && (second_mpi_comm != MPI_COMM_NULL) )
    {
        constexpr int color = 0;
        MPI_Comm_split(parent_comm, color, key, &combined_mpi_comm);
    }
    else
    {
        MPI_Comm_split(parent_comm, MPI_UNDEFINED, key, &combined_mpi_comm);
    }

    ParallelEnvironment::RegisterDataCommunicator(
        rNewCommunicatorName, MPIDataCommunicator::Create(combined_mpi_comm), ParallelEnvironment::DoNotMakeDefault);
    return ParallelEnvironment::GetDataCommunicator(rNewCommunicatorName);
}

}
}