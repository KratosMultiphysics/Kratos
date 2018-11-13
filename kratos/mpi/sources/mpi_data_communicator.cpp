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

// Life cycle

MPIDataCommunicator::MPIDataCommunicator(MPI_Comm MPIComm):
    DataCommunicator(),
    mComm(MPIComm)
{}

MPIDataCommunicator::~MPIDataCommunicator()
{}

DataCommunicator::UniquePointer MPIDataCommunicator::Clone() const
{
    return Kratos::make_unique<MPIDataCommunicator>(mComm);
}

// Barrier wrapper

void MPIDataCommunicator::Barrier() const
{
    MPI_Barrier(mComm);
}

// Reduce operations

int MPIDataCommunicator::Sum(const int rLocalValue, const int Root) const
{
    int global_value(rLocalValue);
    int ierr = MPI_Reduce(&rLocalValue, &global_value, 1, Internals::MPIDatatype(rLocalValue), MPI_SUM, Root, mComm);
    CheckMPIErrorCode(ierr, "MPI_Reduce");
    return global_value;
}

double MPIDataCommunicator::Sum(const double rLocalValue, const int Root) const
{
    double global_value(rLocalValue);
    int ierr = MPI_Reduce(&rLocalValue, &global_value, 1, Internals::MPIDatatype(rLocalValue), MPI_SUM, Root, mComm);
    CheckMPIErrorCode(ierr, "MPI_Reduce");
    return global_value;
}

array_1d<double,3> MPIDataCommunicator::Sum(const array_1d<double,3>& rLocalValue, const int Root) const
{
    array_1d<double,3> global_value(rLocalValue);
    int ierr = MPI_Reduce(
        rLocalValue.data().data(), global_value.data().data(), 3,
        Internals::MPIDatatype(rLocalValue), MPI_SUM, Root, mComm);
    CheckMPIErrorCode(ierr, "MPI_Reduce");
    return global_value;
}

void MPIDataCommunicator::Sum(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rGlobalValues,
        const int Root) const
{
    const int message_size = rLocalValues.size();
    int ierr = MPI_Reduce(
        rLocalValues.data(), rGlobalValues.data(), message_size, MPI_INT,
        MPI_SUM, Root, mComm);
    CheckMPIErrorCode(ierr, "MPI_Reduce");
}

void MPIDataCommunicator::Sum(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rGlobalValues,
        const int Root) const
{
    const int message_size = rLocalValues.size();
    int ierr = MPI_Reduce(
        rLocalValues.data(), rGlobalValues.data(), message_size, MPI_DOUBLE,
        MPI_SUM, Root, mComm);
    CheckMPIErrorCode(ierr, "MPI_Reduce");
}

int MPIDataCommunicator::Min(const int rLocalValue, const int Root) const
{
    int global_value(rLocalValue);
    int ierr = MPI_Reduce(&rLocalValue, &global_value, 1, Internals::MPIDatatype(rLocalValue), MPI_MIN, Root, mComm);
    CheckMPIErrorCode(ierr, "MPI_Reduce");
    return global_value;
}

double MPIDataCommunicator::Min(const double rLocalValue, const int Root) const
{
    double global_value(rLocalValue);
    int ierr = MPI_Reduce(&rLocalValue, &global_value, 1, Internals::MPIDatatype(rLocalValue), MPI_MIN, Root, mComm);
    CheckMPIErrorCode(ierr, "MPI_Reduce");
    return global_value;
}

array_1d<double,3> MPIDataCommunicator::Min(const array_1d<double,3>& rLocalValue, const int Root) const
{
    array_1d<double,3> global_value(rLocalValue);
    int ierr = MPI_Reduce(
        rLocalValue.data().data(), global_value.data().data(), 3,
        Internals::MPIDatatype(rLocalValue), MPI_MIN, Root, mComm);
    CheckMPIErrorCode(ierr, "MPI_Reduce");
    return global_value;
}

void MPIDataCommunicator::Min(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rGlobalValues,
        const int Root) const
{
    const int message_size = rLocalValues.size();
    int ierr = MPI_Reduce(
        rLocalValues.data(), rGlobalValues.data(), message_size, MPI_INT,
        MPI_MIN, Root, mComm);
    CheckMPIErrorCode(ierr, "MPI_Reduce");
}

void MPIDataCommunicator::Min(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rGlobalValues,
        const int Root) const
{
    const int message_size = rLocalValues.size();
    int ierr = MPI_Reduce(
        rLocalValues.data(), rGlobalValues.data(), message_size, MPI_DOUBLE,
        MPI_MIN, Root, mComm);
    CheckMPIErrorCode(ierr, "MPI_Reduce");
}

int MPIDataCommunicator::Max(const int rLocalValue, const int Root) const
{
    int global_value(rLocalValue);
    int ierr = MPI_Reduce(&rLocalValue, &global_value, 1, Internals::MPIDatatype(rLocalValue), MPI_MAX, Root, mComm);
    CheckMPIErrorCode(ierr, "MPI_Reduce");
    return global_value;
}

double MPIDataCommunicator::Max(const double rLocalValue, const int Root) const
{
    double global_value(rLocalValue);
    int ierr = MPI_Reduce(&rLocalValue, &global_value, 1, Internals::MPIDatatype(rLocalValue), MPI_MAX, Root, mComm);
    CheckMPIErrorCode(ierr, "MPI_Reduce");
    return global_value;
}

array_1d<double,3> MPIDataCommunicator::Max(const array_1d<double,3>& rLocalValue, const int Root) const
{
    array_1d<double,3> global_value(rLocalValue);
    int ierr = MPI_Reduce(
        rLocalValue.data().data(), global_value.data().data(), 3,
        Internals::MPIDatatype(rLocalValue), MPI_MAX, Root, mComm);
    CheckMPIErrorCode(ierr, "MPI_Reduce");
    return global_value;
}

void MPIDataCommunicator::Max(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rGlobalValues,
        const int Root) const
{
    const int message_size = rLocalValues.size();
    int ierr = MPI_Reduce(
        rLocalValues.data(), rGlobalValues.data(), message_size, MPI_INT,
        MPI_MAX, Root, mComm);
    CheckMPIErrorCode(ierr, "MPI_Reduce");
}

void MPIDataCommunicator::Max(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rGlobalValues,
        const int Root) const
{
    const int message_size = rLocalValues.size();
    int ierr = MPI_Reduce(
        rLocalValues.data(), rGlobalValues.data(), message_size, MPI_DOUBLE,
        MPI_MAX, Root, mComm);
    CheckMPIErrorCode(ierr, "MPI_Reduce");
}

// Allreduce operations

int MPIDataCommunicator::SumAll(const int rLocalValue) const
{
    int global_value(rLocalValue);
    int ierr = MPI_Allreduce(&rLocalValue, &global_value, 1, Internals::MPIDatatype(rLocalValue), MPI_SUM, mComm);
    CheckMPIErrorCode(ierr, "MPI_Allreduce");
    return global_value;
}

double MPIDataCommunicator::SumAll(const double rLocalValue) const
{
    double global_value(rLocalValue);
    int ierr = MPI_Allreduce(&rLocalValue, &global_value, 1, Internals::MPIDatatype(rLocalValue), MPI_SUM, mComm);
    CheckMPIErrorCode(ierr, "MPI_Allreduce");
    return global_value;
}

array_1d<double,3> MPIDataCommunicator::SumAll(const array_1d<double,3>& rLocalValue) const
{
    array_1d<double,3> global_value(rLocalValue);
    int ierr = MPI_Allreduce(
        rLocalValue.data().data(), global_value.data().data(), 3,
        Internals::MPIDatatype(rLocalValue), MPI_SUM, mComm);
    CheckMPIErrorCode(ierr, "MPI_Allreduce");
    return global_value;
}

void MPIDataCommunicator::SumAll(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rGlobalValues) const
{
    const int message_size = rLocalValues.size();
    int ierr = MPI_Allreduce(
        rLocalValues.data(), rGlobalValues.data(), message_size, MPI_INT,
        MPI_SUM, mComm);
    CheckMPIErrorCode(ierr, "MPI_Allreduce");
}

void MPIDataCommunicator::SumAll(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rGlobalValues) const
{
    const int message_size = rLocalValues.size();
    int ierr = MPI_Allreduce(
        rLocalValues.data(), rGlobalValues.data(), message_size, MPI_DOUBLE,
        MPI_SUM, mComm);
    CheckMPIErrorCode(ierr, "MPI_Allreduce");
}

int MPIDataCommunicator::MinAll(const int rLocalValue) const
{
    int global_value(rLocalValue);
    int ierr = MPI_Allreduce(&rLocalValue, &global_value, 1, Internals::MPIDatatype(rLocalValue), MPI_MIN, mComm);
    CheckMPIErrorCode(ierr, "MPI_Allreduce");
    return global_value;
}

double MPIDataCommunicator::MinAll(const double rLocalValue) const
{
    double global_value(rLocalValue);
    int ierr = MPI_Allreduce(&rLocalValue, &global_value, 1, Internals::MPIDatatype(rLocalValue), MPI_MIN, mComm);
    CheckMPIErrorCode(ierr, "MPI_Allreduce");
    return global_value;
}

array_1d<double,3> MPIDataCommunicator::MinAll(const array_1d<double,3>& rLocalValue) const
{
    array_1d<double,3> global_value(rLocalValue);
    int ierr = MPI_Allreduce(
        rLocalValue.data().data(), global_value.data().data(), 3,
        Internals::MPIDatatype(rLocalValue), MPI_MIN, mComm);
    CheckMPIErrorCode(ierr, "MPI_Allreduce");
    return global_value;
}

void MPIDataCommunicator::MinAll(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rGlobalValues) const
{
    const int message_size = rLocalValues.size();
    int ierr = MPI_Allreduce(
        rLocalValues.data(), rGlobalValues.data(), message_size, MPI_INT,
        MPI_MIN, mComm);
    CheckMPIErrorCode(ierr, "MPI_Allreduce");
}

void MPIDataCommunicator::MinAll(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rGlobalValues) const
{
    const int message_size = rLocalValues.size();
    int ierr = MPI_Allreduce(
        rLocalValues.data(), rGlobalValues.data(), message_size, MPI_DOUBLE,
        MPI_MIN, mComm);
    CheckMPIErrorCode(ierr, "MPI_Allreduce");
}

int MPIDataCommunicator::MaxAll(const int rLocalValue) const
{
    int global_value(rLocalValue);
    int ierr = MPI_Allreduce(&rLocalValue, &global_value, 1, Internals::MPIDatatype(rLocalValue), MPI_MAX, mComm);
    CheckMPIErrorCode(ierr, "MPI_Allreduce");
    return global_value;
}

double MPIDataCommunicator::MaxAll(const double rLocalValue) const
{
    double global_value(rLocalValue);
    int ierr = MPI_Allreduce(&rLocalValue, &global_value, 1, Internals::MPIDatatype(rLocalValue), MPI_MAX, mComm);
    CheckMPIErrorCode(ierr, "MPI_Allreduce");
    return global_value;
}

array_1d<double,3> MPIDataCommunicator::MaxAll(const array_1d<double,3>& rLocalValue) const
{
    array_1d<double,3> global_value(rLocalValue);
    int ierr = MPI_Allreduce(
        rLocalValue.data().data(), global_value.data().data(), 3,
        Internals::MPIDatatype(rLocalValue), MPI_MAX, mComm);
    CheckMPIErrorCode(ierr, "MPI_Allreduce");
    return global_value;
}

void MPIDataCommunicator::MaxAll(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rGlobalValues) const
{
    const int message_size = rLocalValues.size();
    int ierr = MPI_Allreduce(
        rLocalValues.data(), rGlobalValues.data(), message_size, MPI_INT,
        MPI_MAX, mComm);
    CheckMPIErrorCode(ierr, "MPI_Allreduce");
}

void MPIDataCommunicator::MaxAll(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rGlobalValues) const
{
    const int message_size = rLocalValues.size();
    int ierr = MPI_Allreduce(
        rLocalValues.data(), rGlobalValues.data(), message_size, MPI_DOUBLE,
        MPI_MAX, mComm);
    CheckMPIErrorCode(ierr, "MPI_Allreduce");
}

// Scan operations

int MPIDataCommunicator::ScanSum(const int rLocalValue) const
{
    int partial_total;
    int ierr = MPI_Scan(&rLocalValue, &partial_total, 1, Internals::MPIDatatype(rLocalValue), MPI_SUM, mComm);
    CheckMPIErrorCode(ierr, "MPI_Scan");
    return partial_total;
}

double MPIDataCommunicator::ScanSum(const double rLocalValue) const
{
    double partial_total;
    int ierr = MPI_Scan(&rLocalValue, &partial_total, 1, Internals::MPIDatatype(rLocalValue), MPI_SUM, mComm);
    CheckMPIErrorCode(ierr, "MPI_Scan");
    return partial_total;
}


// Sendrecv operations

void MPIDataCommunicator::SendRecv(
    const std::vector<int>& rSendValues, const unsigned int SendDestination,
    std::vector<int>& rRecvValues, const unsigned int RecvSource) const
{
    const int send_tag = 0;
    const int recv_tag = 0;
    int ierr = MPI_Sendrecv(
        rSendValues.data(), rSendValues.size(), MPI_INT, SendDestination, send_tag,
        rRecvValues.data(), rRecvValues.size(), MPI_INT, RecvSource, recv_tag,
        mComm, MPI_STATUS_IGNORE);
    CheckMPIErrorCode(ierr, "MPI_Sendrecv");
}

void MPIDataCommunicator::SendRecv(
    const std::vector<double>& rSendValues, const unsigned int SendDestination,
    std::vector<double>& rRecvValues, const unsigned int RecvSource) const
{
    const int send_tag = 0;
    const int recv_tag = 0;
    int ierr = MPI_Sendrecv(
        rSendValues.data(), rSendValues.size(), MPI_DOUBLE, SendDestination, send_tag,
        rRecvValues.data(), rRecvValues.size(), MPI_DOUBLE, RecvSource, recv_tag,
        mComm, MPI_STATUS_IGNORE);
    CheckMPIErrorCode(ierr, "MPI_Sendrecv");
}

// Broadcast

void MPIDataCommunicator::Broadcast(
    std::vector<int>& rBuffer,
    const unsigned int SourceRank) const
{
    int ierr = MPI_Bcast(rBuffer.data(), rBuffer.size(), MPI_INT, SourceRank, mComm);
    CheckMPIErrorCode(ierr, "MPI_Bcast");
}

void MPIDataCommunicator::Broadcast(
    std::vector<double>& rBuffer,
    const unsigned int SourceRank) const
{
    int ierr = MPI_Bcast(rBuffer.data(), rBuffer.size(), MPI_DOUBLE, SourceRank, mComm);
    CheckMPIErrorCode(ierr, "MPI_Bcast");
}

// Scatter operations

void MPIDataCommunicator::Scatter(
    const std::vector<int>& rSendValues,
    std::vector<int>& rRecvValues,
    const unsigned int SourceRank) const
{
    const int sends_per_rank = rRecvValues.size();
    int ierr = MPI_Scatter(
        rSendValues.data(), sends_per_rank, MPI_INT,
        rRecvValues.data(), sends_per_rank, MPI_INT,
        SourceRank, mComm);
    CheckMPIErrorCode(ierr, "MPI_Scatter");
}

void MPIDataCommunicator::Scatter(
    const std::vector<double>& rSendValues,
    std::vector<double>& rRecvValues,
    const unsigned int SourceRank) const
{
    const int sends_per_rank = rRecvValues.size();
    int ierr = MPI_Scatter(
        rSendValues.data(), sends_per_rank, MPI_DOUBLE,
        rRecvValues.data(), sends_per_rank, MPI_DOUBLE,
        SourceRank, mComm);
    CheckMPIErrorCode(ierr, "MPI_Scatter");
}

void MPIDataCommunicator::Scatterv(
    const std::vector<int>& rSendValues,
    const std::vector<int>& rSendCounts,
    const std::vector<int>& rSendOffsets,
    std::vector<int>& rRecvValues,
    const unsigned int SourceRank) const
{
    int ierr = MPI_Scatterv(
        rSendValues.data(), rSendCounts.data(), rSendOffsets.data(), MPI_INT,
        rRecvValues.data(), rRecvValues.size(), MPI_INT,
        SourceRank, mComm);
    CheckMPIErrorCode(ierr, "MPI_Scatterv");
}

void MPIDataCommunicator::Scatterv(
    const std::vector<double>& rSendValues,
    const std::vector<int>& rSendCounts,
    const std::vector<int>& rSendOffsets,
    std::vector<double>& rRecvValues,
    const unsigned int SourceRank) const
{
    int ierr = MPI_Scatterv(
        rSendValues.data(), rSendCounts.data(), rSendOffsets.data(), MPI_DOUBLE,
        rRecvValues.data(), rRecvValues.size(), MPI_DOUBLE,
        SourceRank, mComm);
    CheckMPIErrorCode(ierr, "MPI_Scatterv");
}

// Gather operations

void MPIDataCommunicator::Gather(
    const std::vector<int>& rSendValues,
    std::vector<int>& rRecvValues,
    const unsigned int DestinationRank) const
{
    const int sends_per_rank = rSendValues.size();
    int ierr = MPI_Gather(
        rSendValues.data(), sends_per_rank, MPI_INT,
        rRecvValues.data(), sends_per_rank, MPI_INT,
        DestinationRank, mComm);
    CheckMPIErrorCode(ierr, "MPI_Gather");
}

void MPIDataCommunicator::Gather(
    const std::vector<double>& rSendValues,
    std::vector<double>& rRecvValues,
    const unsigned int DestinationRank) const
{
    const int sends_per_rank = rSendValues.size();
    int ierr = MPI_Gather(
        rSendValues.data(), sends_per_rank, MPI_DOUBLE,
        rRecvValues.data(), sends_per_rank, MPI_DOUBLE,
        DestinationRank, mComm);
    CheckMPIErrorCode(ierr, "MPI_Gather");
}

void MPIDataCommunicator::Gatherv(
    const std::vector<int>& rSendValues,
    std::vector<int>& rRecvValues,
    const std::vector<int>& rRecvCounts,
    const std::vector<int>& rRecvOffsets,
    const unsigned int DestinationRank) const
{
    int ierr = MPI_Gatherv(
        rSendValues.data(), rSendValues.size(), MPI_INT,
        rRecvValues.data(), rRecvCounts.data(), rRecvOffsets.data(), MPI_INT,
        DestinationRank, mComm);
    CheckMPIErrorCode(ierr, "MPI_Gatherv");
}

void MPIDataCommunicator::Gatherv(
    const std::vector<double>& rSendValues,
    std::vector<double>& rRecvValues,
    const std::vector<int>& rRecvCounts,
    const std::vector<int>& rRecvOffsets,
    const unsigned int DestinationRank) const
{
    int ierr = MPI_Gatherv(
        rSendValues.data(), rSendValues.size(), MPI_DOUBLE,
        rRecvValues.data(), rRecvCounts.data(), rRecvOffsets.data(), MPI_DOUBLE,
        DestinationRank, mComm);
    CheckMPIErrorCode(ierr, "MPI_Gatherv");
}

void MPIDataCommunicator::AllGather(
    const std::vector<int>& rSendValues,
    std::vector<int>& rRecvValues) const
{
    const int sends_per_rank = rSendValues.size();
    int ierr = MPI_Allgather(
        rSendValues.data(), sends_per_rank, MPI_INT,
        rRecvValues.data(), sends_per_rank, MPI_INT,
        mComm);
    CheckMPIErrorCode(ierr, "MPI_Allgather");
}

void MPIDataCommunicator::AllGather(
    const std::vector<double>& rSendValues,
    std::vector<double>& rRecvValues) const
{
    const int sends_per_rank = rSendValues.size();
    int ierr = MPI_Allgather(
        rSendValues.data(), sends_per_rank, MPI_DOUBLE,
        rRecvValues.data(), sends_per_rank, MPI_DOUBLE,
        mComm);
    CheckMPIErrorCode(ierr, "MPI_Allgather");
}

// Access

MPI_Comm MPIDataCommunicator::GetMPICommunicator() const
{
    return mComm;
}

// Inquiry

int MPIDataCommunicator::Rank() const
{
    int rank;
    int ierr = MPI_Comm_rank(mComm, &rank);
    CheckMPIErrorCode(ierr, "MPI_Comm_rank");
    return rank;
}

int MPIDataCommunicator::Size() const
{
    int size;
    int ierr = MPI_Comm_size(mComm, &size);
    CheckMPIErrorCode(ierr, "MPI_Comm_size");
    return size;
}

bool MPIDataCommunicator::IsDistributed() const
{
    return true;
}

// IO

std::string MPIDataCommunicator::Info() const
{
    std::stringstream buffer;
    PrintInfo(buffer);
    return buffer.str();
}

void MPIDataCommunicator::PrintInfo(std::ostream &rOStream) const
{
    rOStream << "MPIDataCommunicator";
}

void MPIDataCommunicator::PrintData(std::ostream &rOStream) const
{}

// Error checking

void MPIDataCommunicator::CheckMPIErrorCode(const int ierr, const std::string MPICallName) const
{
    KRATOS_ERROR_IF_NOT(ierr == MPI_SUCCESS) << MPICallName << " failed with error code " << ierr << "." << std::endl;
}

}