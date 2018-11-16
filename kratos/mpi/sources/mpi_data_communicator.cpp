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

template<> inline MPI_Datatype MPIDatatype<std::vector<double>>(const std::vector<double>&)
{
    return MPI_DOUBLE;
}

template<> inline MPI_Datatype MPIDatatype<double>(const double&)
{
    return MPI_DOUBLE;
}

template<> inline MPI_Datatype MPIDatatype<array_1d<double,3>>(const array_1d<double,3>&)
{
    return MPI_DOUBLE;
}

template<> inline MPI_Datatype MPIDatatype<std::vector<int>>(const std::vector<int>&)
{
    return MPI_INT;
}

// Buffer argument deduction

template<> inline void* GetData(int& rValues)
{
    return &rValues;
}

template<> inline const void* GetData(const int& rValues)
{
    return &rValues;
}

template<> inline void* GetData(double& rValues)
{
    return &rValues;
}

template<> inline const void* GetData(const double& rValues)
{
    return &rValues;
}

template<> inline void* GetData(array_1d<double,3>& rValues)
{
    return rValues.data().data();
}

template<> inline const void* GetData(const array_1d<double,3>& rValues)
{
    return rValues.data().data();
}

template<> inline void* GetData(std::vector<int>& rValues)
{
    return rValues.data();
}

template<> inline const void* GetData(const std::vector<int>& rValues)
{
    return rValues.data();
}

template<> inline void* GetData(std::vector<double>& rValues)
{
    return rValues.data();
}

template<> inline const void* GetData(const std::vector<double>& rValues)
{
    return rValues.data();
}

// MPI message size deduction

template<> inline int MessageSize(const int& rValues)
{
    return 1;
}

template<> inline int MessageSize(const double& rValues)
{
    return 1;
}

template<> inline int MessageSize(const array_1d<double,3>& rValues)
{
    return 3;
}

template<> inline int MessageSize(const std::vector<int>& rValues)
{
    return rValues.size();
}

template<> inline int MessageSize(const std::vector<double>& rValues)
{
    return rValues.size();
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
    int ierr = MPI_Barrier(mComm);
    CheckMPIErrorCode(ierr,"MPI_Barrier");
}

// Reduce operations

int MPIDataCommunicator::Sum(const int rLocalValue, const int Root) const
{
    int global_value(rLocalValue);
    ReduceDetail(rLocalValue, global_value, MPI_SUM, Root);
    return global_value;
}

double MPIDataCommunicator::Sum(const double rLocalValue, const int Root) const
{
    double global_value(rLocalValue);
    ReduceDetail(rLocalValue, global_value, MPI_SUM, Root);
    return global_value;
}

array_1d<double,3> MPIDataCommunicator::Sum(const array_1d<double,3>& rLocalValue, const int Root) const
{
    array_1d<double,3> global_value(rLocalValue);
    ReduceDetail(rLocalValue, global_value, MPI_SUM, Root);
    return global_value;
}

void MPIDataCommunicator::Sum(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rGlobalValues,
        const int Root) const
{
    ReduceDetail(rLocalValues, rGlobalValues, MPI_SUM, Root);
}

void MPIDataCommunicator::Sum(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rGlobalValues,
        const int Root) const
{
    ReduceDetail(rLocalValues, rGlobalValues, MPI_SUM, Root);
}

int MPIDataCommunicator::Min(const int rLocalValue, const int Root) const
{
    int global_value(rLocalValue);
    ReduceDetail(rLocalValue, global_value, MPI_MIN, Root);
    return global_value;
}

double MPIDataCommunicator::Min(const double rLocalValue, const int Root) const
{
    double global_value(rLocalValue);
    ReduceDetail(rLocalValue, global_value, MPI_MIN, Root);
    return global_value;
}

array_1d<double,3> MPIDataCommunicator::Min(const array_1d<double,3>& rLocalValue, const int Root) const
{
    array_1d<double,3> global_value(rLocalValue);
    ReduceDetail(rLocalValue, global_value, MPI_MIN, Root);
    return global_value;
}

void MPIDataCommunicator::Min(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rGlobalValues,
        const int Root) const
{
    ReduceDetail(rLocalValues,rGlobalValues,MPI_MIN,Root);
}

void MPIDataCommunicator::Min(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rGlobalValues,
        const int Root) const
{
    ReduceDetail(rLocalValues,rGlobalValues,MPI_MIN,Root);
}

int MPIDataCommunicator::Max(const int rLocalValue, const int Root) const
{
    int global_value(rLocalValue);
    ReduceDetail(rLocalValue,global_value,MPI_MAX,Root);
    return global_value;
}

double MPIDataCommunicator::Max(const double rLocalValue, const int Root) const
{
    double global_value(rLocalValue);
    ReduceDetail(rLocalValue,global_value,MPI_MAX,Root);
    return global_value;
}

array_1d<double,3> MPIDataCommunicator::Max(const array_1d<double,3>& rLocalValue, const int Root) const
{
    array_1d<double,3> global_value(rLocalValue);
    ReduceDetail(rLocalValue,global_value,MPI_MAX,Root);
    return global_value;
}

void MPIDataCommunicator::Max(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rGlobalValues,
        const int Root) const
{
    ReduceDetail(rLocalValues,rGlobalValues,MPI_MAX,Root);
}

void MPIDataCommunicator::Max(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rGlobalValues,
        const int Root) const
{
    ReduceDetail(rLocalValues,rGlobalValues,MPI_MAX,Root);
}

// Allreduce operations

int MPIDataCommunicator::SumAll(const int rLocalValue) const
{
    int global_value(rLocalValue);
    AllReduceDetail(rLocalValue,global_value,MPI_SUM);
    return global_value;
}

double MPIDataCommunicator::SumAll(const double rLocalValue) const
{
    double global_value(rLocalValue);
    AllReduceDetail(rLocalValue,global_value,MPI_SUM);
    return global_value;
}

array_1d<double,3> MPIDataCommunicator::SumAll(const array_1d<double,3>& rLocalValue) const
{
    array_1d<double,3> global_value(rLocalValue);
    AllReduceDetail(rLocalValue,global_value,MPI_SUM);
    return global_value;
}

void MPIDataCommunicator::SumAll(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rGlobalValues) const
{
    AllReduceDetail(rLocalValues,rGlobalValues,MPI_SUM);
}

void MPIDataCommunicator::SumAll(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rGlobalValues) const
{
    AllReduceDetail(rLocalValues,rGlobalValues,MPI_SUM);
}

int MPIDataCommunicator::MinAll(const int rLocalValue) const
{
    int global_value(rLocalValue);
    AllReduceDetail(rLocalValue,global_value,MPI_MIN);
    return global_value;
}

double MPIDataCommunicator::MinAll(const double rLocalValue) const
{
    double global_value(rLocalValue);
    AllReduceDetail(rLocalValue,global_value,MPI_MIN);
    return global_value;
}

array_1d<double,3> MPIDataCommunicator::MinAll(const array_1d<double,3>& rLocalValue) const
{
    array_1d<double,3> global_value(rLocalValue);
    AllReduceDetail(rLocalValue,global_value,MPI_MIN);
    return global_value;
}

void MPIDataCommunicator::MinAll(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rGlobalValues) const
{
    AllReduceDetail(rLocalValues,rGlobalValues,MPI_MIN);
}

void MPIDataCommunicator::MinAll(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rGlobalValues) const
{
    AllReduceDetail(rLocalValues,rGlobalValues,MPI_MIN);
}

int MPIDataCommunicator::MaxAll(const int rLocalValue) const
{
    int global_value(rLocalValue);
    AllReduceDetail(rLocalValue,global_value,MPI_MAX);
    return global_value;
}

double MPIDataCommunicator::MaxAll(const double rLocalValue) const
{
    double global_value(rLocalValue);
    AllReduceDetail(rLocalValue,global_value,MPI_MAX);
    return global_value;
}

array_1d<double,3> MPIDataCommunicator::MaxAll(const array_1d<double,3>& rLocalValue) const
{
    array_1d<double,3> global_value(rLocalValue);
    AllReduceDetail(rLocalValue,global_value,MPI_MAX);
    return global_value;
}

void MPIDataCommunicator::MaxAll(
        const std::vector<int>& rLocalValues,
        std::vector<int>& rGlobalValues) const
{
    AllReduceDetail(rLocalValues,rGlobalValues,MPI_MAX);
}

void MPIDataCommunicator::MaxAll(
        const std::vector<double>& rLocalValues,
        std::vector<double>& rGlobalValues) const
{
    AllReduceDetail(rLocalValues,rGlobalValues,MPI_MAX);
}

// Scan operations

int MPIDataCommunicator::ScanSum(const int rLocalValue) const
{
    int partial_total;
    ScanDetail(rLocalValue,partial_total,MPI_SUM);
    return partial_total;
}

double MPIDataCommunicator::ScanSum(const double rLocalValue) const
{
    double partial_total;
    ScanDetail(rLocalValue,partial_total,MPI_SUM);
    return partial_total;
}

void MPIDataCommunicator::ScanSum(
    const std::vector<int>& rLocalValues, std::vector<int>& rPartialSums) const
{
    ScanDetail(rLocalValues,rPartialSums,MPI_SUM);
}

void MPIDataCommunicator::ScanSum(
    const std::vector<double>& rLocalValues, std::vector<double>& rPartialSums) const
{
    ScanDetail(rLocalValues,rPartialSums,MPI_SUM);
}

// Sendrecv operations

void MPIDataCommunicator::SendRecv(
    const std::vector<int>& rSendValues, const int SendDestination,
    std::vector<int>& rRecvValues, const int RecvSource) const
{
    SendRecvDetail(rSendValues,SendDestination,rRecvValues,RecvSource);
}

void MPIDataCommunicator::SendRecv(
    const std::vector<double>& rSendValues, const int SendDestination,
    std::vector<double>& rRecvValues, const int RecvSource) const
{
    SendRecvDetail(rSendValues,SendDestination,rRecvValues,RecvSource);
}

// Broadcast

void MPIDataCommunicator::Broadcast(
    int& rBuffer,
    const int SourceRank) const
{
    BroadcastDetail(rBuffer,SourceRank);
}

void MPIDataCommunicator::Broadcast(
        double& rBuffer,
        const int SourceRank) const
{
    BroadcastDetail(rBuffer,SourceRank);
}

void MPIDataCommunicator::Broadcast(
    std::vector<int>& rBuffer,
    const int SourceRank) const
{
    BroadcastDetail(rBuffer,SourceRank);
}

void MPIDataCommunicator::Broadcast(
    std::vector<double>& rBuffer,
    const int SourceRank) const
{
    BroadcastDetail(rBuffer,SourceRank);
}

// Scatter operations

void MPIDataCommunicator::Scatter(
    const std::vector<int>& rSendValues,
    std::vector<int>& rRecvValues,
    const int SourceRank) const
{
    ScatterDetail(rSendValues,rRecvValues,SourceRank);
}

void MPIDataCommunicator::Scatter(
    const std::vector<double>& rSendValues,
    std::vector<double>& rRecvValues,
    const int SourceRank) const
{
    ScatterDetail(rSendValues,rRecvValues,SourceRank);
}

void MPIDataCommunicator::Scatterv(
    const std::vector<int>& rSendValues,
    const std::vector<int>& rSendCounts,
    const std::vector<int>& rSendOffsets,
    std::vector<int>& rRecvValues,
    const int SourceRank) const
{
    ScattervDetail(rSendValues,rSendCounts,rSendOffsets,rRecvValues,SourceRank);
}

void MPIDataCommunicator::Scatterv(
    const std::vector<double>& rSendValues,
    const std::vector<int>& rSendCounts,
    const std::vector<int>& rSendOffsets,
    std::vector<double>& rRecvValues,
    const int SourceRank) const
{
    ScattervDetail(rSendValues,rSendCounts,rSendOffsets,rRecvValues,SourceRank);
}

// Gather operations

void MPIDataCommunicator::Gather(
    const std::vector<int>& rSendValues,
    std::vector<int>& rRecvValues,
    const int DestinationRank) const
{
    GatherDetail(rSendValues, rRecvValues, DestinationRank);
}

void MPIDataCommunicator::Gather(
    const std::vector<double>& rSendValues,
    std::vector<double>& rRecvValues,
    const int DestinationRank) const
{
    GatherDetail(rSendValues, rRecvValues, DestinationRank);
}

void MPIDataCommunicator::Gatherv(
    const std::vector<int>& rSendValues,
    std::vector<int>& rRecvValues,
    const std::vector<int>& rRecvCounts,
    const std::vector<int>& rRecvOffsets,
    const int DestinationRank) const
{
    GathervDetail(rSendValues,rRecvValues,rRecvCounts,rRecvOffsets,DestinationRank);
}

void MPIDataCommunicator::Gatherv(
    const std::vector<double>& rSendValues,
    std::vector<double>& rRecvValues,
    const std::vector<int>& rRecvCounts,
    const std::vector<int>& rRecvOffsets,
    const int DestinationRank) const
{
    GathervDetail(rSendValues,rRecvValues,rRecvCounts,rRecvOffsets,DestinationRank);
}

void MPIDataCommunicator::AllGather(
    const std::vector<int>& rSendValues,
    std::vector<int>& rRecvValues) const
{
    AllGatherDetail(rSendValues,rRecvValues);
}

void MPIDataCommunicator::AllGather(
    const std::vector<double>& rSendValues,
    std::vector<double>& rRecvValues) const
{
    AllGatherDetail(rSendValues,rRecvValues);
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

// Implementation details of MPI calls

template<class TDataType> void MPIDataCommunicator::ReduceDetail(
    const TDataType& rLocalValues, TDataType& rReducedValues,
    MPI_Op Operation, const int Root) const
{
    int ierr = MPI_Reduce(
        Internals::GetData(rLocalValues), Internals::GetData(rReducedValues),
        Internals::MessageSize(rLocalValues), Internals::MPIDatatype(rLocalValues),
        Operation, Root, mComm);
    CheckMPIErrorCode(ierr, "MPI_Reduce");
}

template<class TDataType> void MPIDataCommunicator::AllReduceDetail(
    const TDataType& rLocalValues, TDataType& rReducedValues,
    MPI_Op Operation) const
{
    int ierr = MPI_Allreduce(
        Internals::GetData(rLocalValues), Internals::GetData(rReducedValues),
        Internals::MessageSize(rLocalValues), Internals::MPIDatatype(rLocalValues),
        Operation, mComm);
    CheckMPIErrorCode(ierr, "MPI_Allreduce");
}

template<class TDataType> void MPIDataCommunicator::ScanDetail(
    const TDataType& rLocalValues, TDataType& rReducedValues,
    MPI_Op Operation) const
{
    int ierr = MPI_Scan(
        Internals::GetData(rLocalValues), Internals::GetData(rReducedValues),
        Internals::MessageSize(rLocalValues), Internals::MPIDatatype(rLocalValues),
        Operation, mComm);
    CheckMPIErrorCode(ierr, "MPI_Scan");
}

template<class TDataType> void MPIDataCommunicator::SendRecvDetail(
    const TDataType& rSendMessage, const int SendRank,
    TDataType& rRecvMessage, const int RecvRank) const
{
    const int send_tag = 0;
    const int recv_tag = 0;
    int ierr = MPI_Sendrecv(
        Internals::GetData(rSendMessage), Internals::MessageSize(rSendMessage),
        Internals::MPIDatatype(rSendMessage), SendRank, send_tag,
        Internals::GetData(rRecvMessage), Internals::MessageSize(rRecvMessage),
        Internals::MPIDatatype(rRecvMessage), RecvRank, recv_tag,
        mComm, MPI_STATUS_IGNORE);
    CheckMPIErrorCode(ierr, "MPI_Sendrecv");
}

template<class TDataType> void MPIDataCommunicator::BroadcastDetail(
    TDataType& rBuffer, const int SourceRank) const
{
    int ierr = MPI_Bcast(
        Internals::GetData(rBuffer), Internals::MessageSize(rBuffer),
        Internals::MPIDatatype(rBuffer), SourceRank, mComm);
    CheckMPIErrorCode(ierr, "MPI_Bcast");
}

template<class TDataType> void MPIDataCommunicator::ScatterDetail(
    const TDataType& rSendValues, TDataType& rRecvValues, const int SourceRank) const
{
    const int sends_per_rank = Internals::MessageSize(rRecvValues);
    int ierr = MPI_Scatter(
        Internals::GetData(rSendValues), sends_per_rank, Internals::MPIDatatype(rSendValues),
        Internals::GetData(rRecvValues), sends_per_rank, Internals::MPIDatatype(rRecvValues),
        SourceRank, mComm);
    CheckMPIErrorCode(ierr, "MPI_Scatter");
}

template<class TDataType> void MPIDataCommunicator::ScattervDetail(
        const TDataType& rSendValues, const std::vector<int>& rSendCounts, const std::vector<int>& rSendOffsets,
        TDataType& rRecvValues, const int SourceRank) const
{
    int ierr = MPI_Scatterv(
        Internals::GetData(rSendValues), rSendCounts.data(), rSendOffsets.data(), Internals::MPIDatatype(rSendValues),
        Internals::GetData(rRecvValues), Internals::MessageSize(rRecvValues), Internals::MPIDatatype(rRecvValues),
        SourceRank, mComm);
    CheckMPIErrorCode(ierr, "MPI_Scatterv");
}

template<class TDataType> void MPIDataCommunicator::GatherDetail(
    const TDataType& rSendValues, TDataType& rRecvValues, const int RecvRank) const
{
    const int sends_per_rank = Internals::MessageSize(rSendValues);
    int ierr = MPI_Gather(
        Internals::GetData(rSendValues), sends_per_rank, Internals::MPIDatatype(rSendValues),
        Internals::GetData(rRecvValues), sends_per_rank, Internals::MPIDatatype(rRecvValues),
        RecvRank, mComm);
    CheckMPIErrorCode(ierr, "MPI_Gather");
}


template<class TDataType> void MPIDataCommunicator::GathervDetail(
    const TDataType& rSendValues, TDataType& rRecvValues,
    const std::vector<int>& rRecvCounts, const std::vector<int>& rRecvOffsets,
    const int RecvRank) const
{
    int ierr = MPI_Gatherv(
        Internals::GetData(rSendValues), Internals::MessageSize(rSendValues), Internals::MPIDatatype(rSendValues),
        Internals::GetData(rRecvValues), rRecvCounts.data(), rRecvOffsets.data(), Internals::MPIDatatype(rRecvValues),
        RecvRank, mComm);
    CheckMPIErrorCode(ierr, "MPI_Gatherv");
}

template<class TDataType> void MPIDataCommunicator::AllGatherDetail(
    const TDataType& rSendValues, TDataType& rRecvValues) const
{
    const int sends_per_rank = Internals::MessageSize(rSendValues);
    int ierr = MPI_Allgather(
        Internals::GetData(rSendValues), sends_per_rank, Internals::MPIDatatype(rSendValues),
        Internals::GetData(rRecvValues), sends_per_rank, Internals::MPIDatatype(rRecvValues),
        mComm);
    CheckMPIErrorCode(ierr, "MPI_Allgather");
}

}