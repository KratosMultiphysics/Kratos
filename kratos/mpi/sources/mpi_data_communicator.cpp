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

void MPIDataCommunicator::Sum(int& rValue, const int Root) const
{
    int local_value(rValue);
    MPI_Reduce(&local_value, &rValue, 1, Internals::MPIDatatype(rValue), MPI_SUM, Root, mComm);
}

void MPIDataCommunicator::Sum(double& rValue, const int Root) const
{
    double local_value(rValue);
    MPI_Reduce(&local_value, &rValue, 1, Internals::MPIDatatype(rValue), MPI_SUM, Root, mComm);
}

void MPIDataCommunicator::Sum(array_1d<double,3>& rValue, const int Root) const
{
    array_1d<double,3> local_value(rValue);
    MPI_Reduce(
        local_value.data().data(), rValue.data().data(), 3,
        Internals::MPIDatatype(rValue), MPI_SUM, Root, mComm);
}

void MPIDataCommunicator::Min(int& rValue, const int Root) const
{
    int local_value(rValue);
    MPI_Reduce(&local_value, &rValue, 1, Internals::MPIDatatype(rValue), MPI_MIN, Root, mComm);
}

void MPIDataCommunicator::Min(double& rValue, const int Root) const
{
    double local_value(rValue);
    MPI_Reduce(&local_value, &rValue, 1, Internals::MPIDatatype(rValue), MPI_MIN, Root, mComm);
}

void MPIDataCommunicator::Min(array_1d<double,3>& rValue, const int Root) const
{
    array_1d<double,3> local_value(rValue);
    MPI_Reduce(
        local_value.data().data(), rValue.data().data(), 3,
        Internals::MPIDatatype(rValue), MPI_MIN, Root, mComm);
}

void MPIDataCommunicator::Max(int& rValue, const int Root) const
{
    int local_value(rValue);
    MPI_Reduce(&local_value, &rValue, 1, Internals::MPIDatatype(rValue), MPI_MAX, Root, mComm);
}

void MPIDataCommunicator::Max(double& rValue, const int Root) const
{
    double local_value(rValue);
    MPI_Reduce(&local_value, &rValue, 1, Internals::MPIDatatype(rValue), MPI_MAX, Root, mComm);
}

void MPIDataCommunicator::Max(array_1d<double,3>& rValue, const int Root) const
{
    array_1d<double,3> local_value(rValue);
    MPI_Reduce(
        local_value.data().data(), rValue.data().data(), 3,
        Internals::MPIDatatype(rValue), MPI_MAX, Root, mComm);
}

// Allreduce operations

void MPIDataCommunicator::SumAll(int& rValue) const
{
    int local_value(rValue);
    MPI_Allreduce(&local_value, &rValue, 1, Internals::MPIDatatype(rValue), MPI_SUM, mComm);
}

void MPIDataCommunicator::SumAll(double& rValue) const
{
    double local_value(rValue);
    MPI_Allreduce(&local_value, &rValue, 1, Internals::MPIDatatype(rValue), MPI_SUM, mComm);
}

void MPIDataCommunicator::SumAll(array_1d<double,3>& rValue) const
{
    array_1d<double,3> local_value(rValue);
    MPI_Allreduce(
        local_value.data().data(), rValue.data().data(), 3,
        Internals::MPIDatatype(rValue), MPI_SUM, mComm);
}

void MPIDataCommunicator::MinAll(int& rValue) const
{
    int local_value(rValue);
    MPI_Allreduce(&local_value, &rValue, 1, Internals::MPIDatatype(rValue), MPI_MIN, mComm);
}

void MPIDataCommunicator::MinAll(double& rValue) const
{
    double local_value(rValue);
    MPI_Allreduce(&local_value, &rValue, 1, Internals::MPIDatatype(rValue), MPI_MIN, mComm);
}

void MPIDataCommunicator::MinAll(array_1d<double,3>& rValue) const
{
    array_1d<double,3> local_value(rValue);
    MPI_Allreduce(
        local_value.data().data(), rValue.data().data(), 3,
        Internals::MPIDatatype(rValue), MPI_MIN, mComm);
}

void MPIDataCommunicator::MaxAll(int& rValue) const
{
    int local_value(rValue);
    MPI_Allreduce(&local_value, &rValue, 1, Internals::MPIDatatype(rValue), MPI_MAX, mComm);
}

void MPIDataCommunicator::MaxAll(double& rValue) const
{
    double local_value(rValue);
    MPI_Allreduce(&local_value, &rValue, 1, Internals::MPIDatatype(rValue), MPI_MAX, mComm);
}

void MPIDataCommunicator::MaxAll(array_1d<double,3>& rValue) const
{
    array_1d<double,3> local_value(rValue);
    MPI_Allreduce(
        local_value.data().data(), rValue.data().data(), 3,
        Internals::MPIDatatype(rValue), MPI_MAX, mComm);
}

// Scan operations

int MPIDataCommunicator::ScanSum(const int& rLocalValue) const
{
    int partial_total;
    MPI_Scan(&rLocalValue, &partial_total, 1, Internals::MPIDatatype(rLocalValue), MPI_SUM, mComm);
    return partial_total;
}

double MPIDataCommunicator::ScanSum(const double& rLocalValue) const
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

}