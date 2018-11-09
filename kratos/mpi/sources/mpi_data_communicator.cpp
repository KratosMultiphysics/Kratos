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

#include "mpi/includes/mpi_data_communicator.h"

namespace Kratos {

namespace Internals {

// MPI_Datatype wrappers

template<> inline MPI_Datatype MPIDatatype<int>(const int&)
{
    return MPI_INT;
}

template<> inline MPI_Datatype MPIDatatype<double>(const double&)
{
    return MPI_DOUBLE;
}

template<> inline MPI_Datatype MPIDatatype<array_1d<double,3>>(const array_1d<double,3>&)
{
    return MPI_DOUBLE;
}

}

// MPIDataCommunicator implementation

// Reduce operations

int MPIDataCommunicator::Sum(const int rLocalValue, const int Root) const
{
    int global_value(rLocalValue);
    MPI_Reduce(&rLocalValue, &global_value, 1, Internals::MPIDatatype(rLocalValue), MPI_SUM, Root, mComm);
    return global_value;
}

double MPIDataCommunicator::Sum(const double rLocalValue, const int Root) const
{
    double global_value(rLocalValue);
    MPI_Reduce(&rLocalValue, &global_value, 1, Internals::MPIDatatype(rLocalValue), MPI_SUM, Root, mComm);
    return global_value;
}

array_1d<double,3> MPIDataCommunicator::Sum(const array_1d<double,3>& rLocalValue, const int Root) const
{
    array_1d<double,3> global_value(rLocalValue);
    MPI_Reduce(
        rLocalValue.data().data(), global_value.data().data(), 3,
        Internals::MPIDatatype(rLocalValue), MPI_SUM, Root, mComm);
    return global_value;
}

int MPIDataCommunicator::Min(const int rLocalValue, const int Root) const
{
    int global_value(rLocalValue);
    MPI_Reduce(&rLocalValue, &global_value, 1, Internals::MPIDatatype(rLocalValue), MPI_MIN, Root, mComm);
    return global_value;
}

double MPIDataCommunicator::Min(const double rLocalValue, const int Root) const
{
    double global_value(rLocalValue);
    MPI_Reduce(&rLocalValue, &global_value, 1, Internals::MPIDatatype(rLocalValue), MPI_MIN, Root, mComm);
    return global_value;
}

array_1d<double,3> MPIDataCommunicator::Min(const array_1d<double,3>& rLocalValue, const int Root) const
{
    array_1d<double,3> global_value(rLocalValue);
    MPI_Reduce(
        rLocalValue.data().data(), global_value.data().data(), 3,
        Internals::MPIDatatype(rLocalValue), MPI_MIN, Root, mComm);
    return global_value;
}

int MPIDataCommunicator::Max(const int rLocalValue, const int Root) const
{
    int global_value(rLocalValue);
    MPI_Reduce(&rLocalValue, &global_value, 1, Internals::MPIDatatype(rLocalValue), MPI_MAX, Root, mComm);
    return global_value;
}

double MPIDataCommunicator::Max(const double rLocalValue, const int Root) const
{
    double global_value(rLocalValue);
    MPI_Reduce(&rLocalValue, &global_value, 1, Internals::MPIDatatype(rLocalValue), MPI_MAX, Root, mComm);
    return global_value;
}

array_1d<double,3> MPIDataCommunicator::Max(const array_1d<double,3>& rLocalValue, const int Root) const
{
    array_1d<double,3> global_value(rLocalValue);
    MPI_Reduce(
        rLocalValue.data().data(), global_value.data().data(), 3,
        Internals::MPIDatatype(rLocalValue), MPI_MAX, Root, mComm);
    return global_value;
}

// Allreduce operations

int MPIDataCommunicator::SumAll(const int rLocalValue) const
{
    int global_value(rLocalValue);
    MPI_Allreduce(&rLocalValue, &global_value, 1, Internals::MPIDatatype(rLocalValue), MPI_SUM, mComm);
    return global_value;
}

double MPIDataCommunicator::SumAll(const double rLocalValue) const
{
    double global_value(rLocalValue);
    MPI_Allreduce(&rLocalValue, &global_value, 1, Internals::MPIDatatype(rLocalValue), MPI_SUM, mComm);
    return global_value;
}

array_1d<double,3> MPIDataCommunicator::SumAll(const array_1d<double,3>& rLocalValue) const
{
    array_1d<double,3> global_value(rLocalValue);
    MPI_Allreduce(
        rLocalValue.data().data(), global_value.data().data(), 3,
        Internals::MPIDatatype(rLocalValue), MPI_SUM, mComm);
    return global_value;
}

int MPIDataCommunicator::MinAll(const int rLocalValue) const
{
    int global_value(rLocalValue);
    MPI_Allreduce(&rLocalValue, &global_value, 1, Internals::MPIDatatype(rLocalValue), MPI_MIN, mComm);
    return global_value;
}

double MPIDataCommunicator::MinAll(const double rLocalValue) const
{
    double global_value(rLocalValue);
    MPI_Allreduce(&rLocalValue, &global_value, 1, Internals::MPIDatatype(rLocalValue), MPI_MIN, mComm);
    return global_value;
}

array_1d<double,3> MPIDataCommunicator::MinAll(const array_1d<double,3>& rLocalValue) const
{
    array_1d<double,3> global_value(rLocalValue);
    MPI_Allreduce(
        rLocalValue.data().data(), global_value.data().data(), 3,
        Internals::MPIDatatype(rLocalValue), MPI_MIN, mComm);
    return global_value;
}

int MPIDataCommunicator::MaxAll(const int rLocalValue) const
{
    int global_value(rLocalValue);
    MPI_Allreduce(&rLocalValue, &global_value, 1, Internals::MPIDatatype(rLocalValue), MPI_MAX, mComm);
    return global_value;
}

double MPIDataCommunicator::MaxAll(const double rLocalValue) const
{
    double global_value(rLocalValue);
    MPI_Allreduce(&rLocalValue, &global_value, 1, Internals::MPIDatatype(rLocalValue), MPI_MAX, mComm);
    return global_value;
}

array_1d<double,3> MPIDataCommunicator::MaxAll(const array_1d<double,3>& rLocalValue) const
{
    array_1d<double,3> global_value(rLocalValue);
    MPI_Allreduce(
        rLocalValue.data().data(), global_value.data().data(), 3,
        Internals::MPIDatatype(rLocalValue), MPI_MAX, mComm);
    return global_value;
}

// Scan operations

int MPIDataCommunicator::ScanSum(const int rLocalValue) const
{
    int partial_total;
    MPI_Scan(&rLocalValue, &partial_total, 1, Internals::MPIDatatype(rLocalValue), MPI_SUM, mComm);
    return partial_total;
}

double MPIDataCommunicator::ScanSum(const double rLocalValue) const
{
    double partial_total;
    MPI_Scan(&rLocalValue, &partial_total, 1, Internals::MPIDatatype(rLocalValue), MPI_SUM, mComm);
    return partial_total;
}


// Sendrecv operations

void MPIDataCommunicator::SendRecv(
    const std::vector<int>& rSendValues, const unsigned int SendDestination,
    std::vector<int>& rRecvValues, const unsigned int RecvSource) const
{
    const int send_tag = 0;
    const int recv_tag = 0;
    MPI_Sendrecv(
        rSendValues.data(), rSendValues.size(), MPI_INT, SendDestination, send_tag,
        rRecvValues.data(), rRecvValues.size(), MPI_INT, RecvSource, recv_tag,
        mComm, MPI_STATUS_IGNORE);
}

void MPIDataCommunicator::SendRecv(
    const std::vector<double>& rSendValues, const unsigned int SendDestination,
    std::vector<double>& rRecvValues, const unsigned int RecvSource) const
{
    const int send_tag = 0;
    const int recv_tag = 0;
    MPI_Sendrecv(
        rSendValues.data(), rSendValues.size(), MPI_DOUBLE, SendDestination, send_tag,
        rRecvValues.data(), rRecvValues.size(), MPI_DOUBLE, RecvSource, recv_tag,
        mComm, MPI_STATUS_IGNORE);
}

// Broadcast

void MPIDataCommunicator::Broadcast(
    std::vector<int>& rBuffer,
    const unsigned int SourceRank) const
{
    MPI_Bcast(rBuffer.data(), rBuffer.size(), MPI_INT, SourceRank, mComm);
}

void MPIDataCommunicator::Broadcast(
    std::vector<double>& rBuffer,
    const unsigned int SourceRank) const
{
    MPI_Bcast(rBuffer.data(), rBuffer.size(), MPI_DOUBLE, SourceRank, mComm);
}

// Scatter operations

void MPIDataCommunicator::Scatter(
    const std::vector<int>& rSendValues,
    std::vector<int>& rRecvValues,
    const unsigned int SourceRank) const
{
    const int sends_per_rank = rRecvValues.size();
    MPI_Scatter(
        rSendValues.data(), sends_per_rank, MPI_INT,
        rRecvValues.data(), sends_per_rank, MPI_INT,
        SourceRank, mComm);
}

void MPIDataCommunicator::Scatter(
    const std::vector<double>& rSendValues,
    std::vector<double>& rRecvValues,
    const unsigned int SourceRank) const
{
    const int sends_per_rank = rRecvValues.size();
    MPI_Scatter(
        rSendValues.data(), sends_per_rank, MPI_DOUBLE,
        rRecvValues.data(), sends_per_rank, MPI_DOUBLE,
        SourceRank, mComm);
}

void MPIDataCommunicator::Scatterv(
    const std::vector<int>& rSendValues,
    const std::vector<int>& rSendCounts,
    const std::vector<int>& rSendOffsets,
    std::vector<int>& rRecvValues,
    const unsigned int SourceRank) const
{
    MPI_Scatterv(
        rSendValues.data(), rSendCounts.data(), rSendOffsets.data(), MPI_INT,
        rRecvValues.data(), rRecvValues.size(), MPI_INT,
        SourceRank, mComm);
}

void MPIDataCommunicator::Scatterv(
    const std::vector<double>& rSendValues,
    const std::vector<int>& rSendCounts,
    const std::vector<int>& rSendOffsets,
    std::vector<double>& rRecvValues,
    const unsigned int SourceRank) const
{
    MPI_Scatterv(
        rSendValues.data(), rSendCounts.data(), rSendOffsets.data(), MPI_DOUBLE,
        rRecvValues.data(), rRecvValues.size(), MPI_DOUBLE,
        SourceRank, mComm);
}

// Gather operations

void MPIDataCommunicator::Gather(
    const std::vector<int>& rSendValues,
    std::vector<int>& rRecvValues,
    const unsigned int DestinationRank) const
{
    const int sends_per_rank = rSendValues.size();
    MPI_Gather(
        rSendValues.data(), sends_per_rank, MPI_INT,
        rRecvValues.data(), sends_per_rank, MPI_INT,
        DestinationRank, mComm);
}

void MPIDataCommunicator::Gather(
    const std::vector<double>& rSendValues,
    std::vector<double>& rRecvValues,
    const unsigned int DestinationRank) const
{
    const int sends_per_rank = rSendValues.size();
    MPI_Gather(
        rSendValues.data(), sends_per_rank, MPI_DOUBLE,
        rRecvValues.data(), sends_per_rank, MPI_DOUBLE,
        DestinationRank, mComm);
}

}