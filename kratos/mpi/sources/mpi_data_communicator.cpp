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

#include "mpi/includes/mpi_data_communicator.h"
#include "mpi/includes/mpi_manager.h"
#include "mpi/includes/mpi_message.h"

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_REDUCE_INTERFACE_FOR_TYPE
#define KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_REDUCE_INTERFACE_FOR_TYPE(type)                                 \
type MPIDataCommunicator::Sum(const type rLocalValue, const int Root) const {                               \
    return ReduceDetail(rLocalValue, MPI_SUM, Root);                                                        \
}                                                                                                           \
std::vector<type> MPIDataCommunicator::Sum(const std::vector<type>& rLocalValues, const int Root) const {   \
    return ReduceDetailVector(rLocalValues, MPI_SUM, Root);                                                 \
}                                                                                                           \
void MPIDataCommunicator::Sum(                                                                              \
    const std::vector<type>& rLocalValues, std::vector<type>& rGlobalValues, const int Root) const {        \
    ReduceDetail(rLocalValues, rGlobalValues, MPI_SUM, Root);                                               \
}                                                                                                           \
type MPIDataCommunicator::Min(const type rLocalValue, const int Root) const {                               \
    return ReduceDetail(rLocalValue, MPI_MIN, Root);                                                        \
}                                                                                                           \
std::vector<type> MPIDataCommunicator::Min(const std::vector<type>& rLocalValues, const int Root) const {   \
    return ReduceDetailVector(rLocalValues, MPI_MIN, Root);                                                 \
}                                                                                                           \
void MPIDataCommunicator::Min(                                                                              \
    const std::vector<type>& rLocalValues, std::vector<type>& rGlobalValues, const int Root) const {        \
    ReduceDetail(rLocalValues, rGlobalValues, MPI_MIN, Root);                                               \
}                                                                                                           \
type MPIDataCommunicator::Max(const type rLocalValue, const int Root) const {                               \
    return ReduceDetail(rLocalValue, MPI_MAX, Root);                                                        \
}                                                                                                           \
std::vector<type> MPIDataCommunicator::Max(const std::vector<type>& rLocalValues, const int Root) const {   \
    return ReduceDetailVector(rLocalValues, MPI_MAX, Root);                                                 \
}                                                                                                           \
void MPIDataCommunicator::Max(                                                                              \
    const std::vector<type>& rLocalValues, std::vector<type>& rGlobalValues, const int Root) const {        \
    ReduceDetail(rLocalValues, rGlobalValues, MPI_MAX, Root);                                               \
}                                                                                                           \

#endif

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_ALLREDUCE_INTERFACE_FOR_TYPE
#define KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_ALLREDUCE_INTERFACE_FOR_TYPE(type)               \
type MPIDataCommunicator::SumAll(const type rLocalValue) const {                             \
    return AllReduceDetail(rLocalValue, MPI_SUM);                                            \
}                                                                                            \
std::vector<type> MPIDataCommunicator::SumAll(const std::vector<type>& rLocalValues) const { \
    return AllReduceDetailVector(rLocalValues, MPI_SUM);                                     \
}                                                                                            \
void MPIDataCommunicator::SumAll(                                                            \
    const std::vector<type>& rLocalValues, std::vector<type>& rGlobalValues) const {         \
    AllReduceDetail(rLocalValues, rGlobalValues, MPI_SUM);                                   \
}                                                                                            \
type MPIDataCommunicator::MinAll(const type rLocalValue) const {                             \
    return AllReduceDetail(rLocalValue, MPI_MIN);                                            \
}                                                                                            \
std::vector<type> MPIDataCommunicator::MinAll(const std::vector<type>& rLocalValues) const { \
    return AllReduceDetailVector(rLocalValues, MPI_MIN);                                     \
}                                                                                            \
void MPIDataCommunicator::MinAll(                                                            \
    const std::vector<type>& rLocalValues, std::vector<type>& rGlobalValues) const {         \
    AllReduceDetail(rLocalValues, rGlobalValues, MPI_MIN);                                   \
}                                                                                            \
type MPIDataCommunicator::MaxAll(const type rLocalValue) const {                             \
    return AllReduceDetail(rLocalValue, MPI_MAX);                                            \
}                                                                                            \
std::vector<type> MPIDataCommunicator::MaxAll(const std::vector<type>& rLocalValues) const { \
    return AllReduceDetailVector(rLocalValues, MPI_MAX);                                     \
}                                                                                            \
void MPIDataCommunicator::MaxAll(                                                            \
    const std::vector<type>& rLocalValues, std::vector<type>& rGlobalValues) const {         \
    AllReduceDetail(rLocalValues, rGlobalValues, MPI_MAX);                                   \
}                                                                                            \

#endif

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_SCANSUM_INTERFACE_FOR_TYPE
#define KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_SCANSUM_INTERFACE_FOR_TYPE(type)                  \
type MPIDataCommunicator::ScanSum(const type LocalValue) const {                              \
    return ScanDetail(LocalValue, MPI_SUM);                                                   \
}                                                                                             \
std::vector<type> MPIDataCommunicator::ScanSum(const std::vector<type>& rLocalValues) const { \
    return ScanDetail(rLocalValues, MPI_SUM);                                                 \
}                                                                                             \
void MPIDataCommunicator::ScanSum(                                                            \
    const std::vector<type>& rLocalValues, std::vector<type>& rPartialSums) const {           \
    ScanDetail(rLocalValues,rPartialSums,MPI_SUM);                                            \
}                                                                                             \

#endif

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_SENDRECV_INTERFACE_FOR_TYPE
#define KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_SENDRECV_INTERFACE_FOR_TYPE(type)           \
type MPIDataCommunicator::SendRecvImpl(                                                 \
    const type SendValue,                                                               \
    const int SendDestination, const int SendTag,                                       \
    const int RecvSource, const int RecvTag) const {                                    \
    return SendRecvDetail(SendValue, SendDestination, SendTag, RecvSource, RecvTag);    \
}                                                                                       \
std::vector<type> MPIDataCommunicator::SendRecvImpl(                                    \
    const std::vector<type>& rSendValues,                                               \
    const int SendDestination, const int SendTag,                                       \
    const int RecvSource, const int RecvTag) const {                                    \
    return SendRecvDetail(rSendValues, SendDestination, SendTag, RecvSource, RecvTag);  \
}                                                                                       \
void MPIDataCommunicator::SendRecvImpl(                                                 \
    const std::vector<type>& rSendValues, const int SendDestination, const int SendTag, \
    std::vector<type>& rRecvValues, const int RecvSource, const int RecvTag) const {    \
    SendRecvDetail(rSendValues,SendDestination,SendTag,rRecvValues,RecvSource,RecvTag); \
}                                                                                       \
void MPIDataCommunicator::SendRecvImpl(                                                 \
    const type SendValue, const int SendDestination, const int SendTag,                 \
    type& rRecvValue, const int RecvSource, const int RecvTag) const {                  \
    SendRecvDetail(SendValue,SendDestination,SendTag,rRecvValue,RecvSource,RecvTag);    \
}                                                                                       \
void MPIDataCommunicator::SendImpl(const std::vector<type>& rSendValues,                \
    const int SendDestination, const int SendTag) const {                               \
    SendDetail(rSendValues, SendDestination, SendTag);                                  \
}                                                                                       \
void MPIDataCommunicator::RecvImpl(std::vector<type>& rRecvValues,                      \
    const int RecvSource, const int RecvTag) const {                                    \
    RecvDetail(rRecvValues, RecvSource, RecvTag);                                       \
}                                                                                       \

#endif

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_BROADCAST_INTERFACE_FOR_TYPE
#define KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_BROADCAST_INTERFACE_FOR_TYPE(type)                      \
void MPIDataCommunicator::BroadcastImpl(type& rBuffer, const int SourceRank) const {                \
    BroadcastDetail(rBuffer,SourceRank);                                                            \
}                                                                                                   \
void MPIDataCommunicator::BroadcastImpl(std::vector<type>& rBuffer, const int SourceRank) const {   \
    BroadcastDetail(rBuffer,SourceRank);                                                            \
}                                                                                                   \

#endif

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_SCATTER_INTERFACE_FOR_TYPE
#define KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_SCATTER_INTERFACE_FOR_TYPE(type)            \
std::vector<type> MPIDataCommunicator::Scatter(                                         \
    const std::vector<type>& rSendValues, const int SourceRank) const {                 \
    return ScatterDetail(rSendValues, SourceRank);                                      \
}                                                                                       \
void MPIDataCommunicator::Scatter(                                                      \
    const std::vector<type>& rSendValues, std::vector<type>& rRecvValues,               \
    const int SourceRank) const {                                                       \
    ScatterDetail(rSendValues,rRecvValues,SourceRank);                                  \
}                                                                                       \
std::vector<type> MPIDataCommunicator::Scatterv(                                        \
    const std::vector<std::vector<type>>& rSendValues, const int SourceRank) const {    \
    return ScattervDetail(rSendValues, SourceRank);                                     \
}                                                                                       \
void MPIDataCommunicator::Scatterv(                                                     \
    const std::vector<type>& rSendValues, const std::vector<int>& rSendCounts,          \
    const std::vector<int>& rSendOffsets, std::vector<type>& rRecvValues,               \
    const int SourceRank) const {                                                       \
    ScattervDetail(rSendValues,rSendCounts,rSendOffsets,rRecvValues,SourceRank);        \
}                                                                                       \

#endif

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_GATHER_INTERFACE_FOR_TYPE
#define KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_GATHER_INTERFACE_FOR_TYPE(type)         \
std::vector<type> MPIDataCommunicator::Gather(                                      \
    const std::vector<type>& rSendValues, const int DestinationRank) const {        \
    return GatherDetail(rSendValues, DestinationRank);                              \
}                                                                                   \
void MPIDataCommunicator::Gather(                                                   \
    const std::vector<type>& rSendValues, std::vector<type>& rRecvValues,           \
    const int DestinationRank) const {                                              \
    GatherDetail(rSendValues, rRecvValues, DestinationRank);                        \
}                                                                                   \
std::vector<std::vector<type>> MPIDataCommunicator::Gatherv(                        \
    const std::vector<type>& rSendValues, const int DestinationRank) const {        \
    return GathervDetail(rSendValues, DestinationRank);                             \
}                                                                                   \
void MPIDataCommunicator::Gatherv(                                                  \
    const std::vector<type>& rSendValues, std::vector<type>& rRecvValues,           \
    const std::vector<int>& rRecvCounts, const std::vector<int>& rRecvOffsets,      \
    const int DestinationRank) const {                                              \
    GathervDetail(rSendValues,rRecvValues,rRecvCounts,rRecvOffsets,DestinationRank);\
}                                                                                   \
std::vector<type> MPIDataCommunicator::AllGather(                                   \
    const std::vector<type>& rSendValues) const {                                   \
    return AllGatherDetail(rSendValues);                                            \
}                                                                                   \
void MPIDataCommunicator::AllGather(                                                \
    const std::vector<type>& rSendValues, std::vector<type>& rRecvValues) const {   \
    AllGatherDetail(rSendValues,rRecvValues);                                       \
}                                                                                   \

#endif

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_PUBLIC_INTERFACE_FOR_TYPE
#define KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_PUBLIC_INTERFACE_FOR_TYPE(type)   \
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_REDUCE_INTERFACE_FOR_TYPE(type)    \
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_ALLREDUCE_INTERFACE_FOR_TYPE(type) \
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_SCANSUM_INTERFACE_FOR_TYPE(type)   \
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_SCATTER_INTERFACE_FOR_TYPE(type)   \
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_GATHER_INTERFACE_FOR_TYPE(type)    \

#endif

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_IMPLEMENTATION_FOR_TYPE
#define KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_IMPLEMENTATION_FOR_TYPE(type)   \
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_SENDRECV_INTERFACE_FOR_TYPE(type)  \
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_BROADCAST_INTERFACE_FOR_TYPE(type) \

#endif

namespace Kratos {
// MPIDataCommunicator implementation

// Life cycle

MPIDataCommunicator::MPIDataCommunicator(MPI_Comm MPIComm):
    DataCommunicator(),
    mComm(MPIComm)
{
    if (!ParallelEnvironment::MPIIsInitialized())
    {
        ParallelEnvironment::SetUpMPIEnvironment(MPIManager::Create());
    }
}

MPIDataCommunicator::~MPIDataCommunicator()
{
    // If the MPI_Comm object is not one of the standard ones, it is our responsibility to manage its lifetime.
    if(mComm != MPI_COMM_WORLD && mComm != MPI_COMM_SELF && mComm != MPI_COMM_NULL)
    {
        MPI_Comm_free(&mComm);
    }
}

// New communicator creation

MPIDataCommunicator::UniquePointer MPIDataCommunicator::Create(MPI_Comm Comm)
{
    return Kratos::make_unique<MPIDataCommunicator>(Comm);
}

// Barrier wrapper

void MPIDataCommunicator::Barrier() const
{
    int ierr = MPI_Barrier(mComm);
    CheckMPIErrorCode(ierr,"MPI_Barrier");
}

// Complete interface for basic types

KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_PUBLIC_INTERFACE_FOR_TYPE(int)
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_PUBLIC_INTERFACE_FOR_TYPE(unsigned int)
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_PUBLIC_INTERFACE_FOR_TYPE(long unsigned int)
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_PUBLIC_INTERFACE_FOR_TYPE(double)

// Reduce operations

array_1d<double,3> MPIDataCommunicator::Sum(const array_1d<double,3>& rLocalValue, const int Root) const
{
    return ReduceDetail(rLocalValue, MPI_SUM, Root);
}

array_1d<double,3> MPIDataCommunicator::Min(const array_1d<double,3>& rLocalValue, const int Root) const
{
    array_1d<double,3> global_value(rLocalValue);
    ReduceDetail(rLocalValue, global_value, MPI_MIN, Root);
    return global_value;
}

array_1d<double,3> MPIDataCommunicator::Max(const array_1d<double,3>& rLocalValue, const int Root) const
{
    array_1d<double,3> global_value(rLocalValue);
    ReduceDetail(rLocalValue,global_value,MPI_MAX,Root);
    return global_value;
}

bool MPIDataCommunicator::AndReduce(const bool Value, const int Root) const
{
    return ReduceDetail(Value, MPI_LAND, Root);
}

Kratos::Flags MPIDataCommunicator::AndReduce(const Kratos::Flags Values, const Kratos::Flags Mask, const int Root) const
{
    Flags::BlockType local_active_flags = Values.GetDefined() & Mask.GetDefined();
    Flags::BlockType active_flags = local_active_flags;
    ReduceDetail(local_active_flags, active_flags, MPI_BOR, Root);

    Flags::BlockType flags = Values.GetFlags();
    Flags::BlockType reduced_flags = flags;
    ReduceDetail(flags, reduced_flags, MPI_BAND, Root);

    Flags out;
    out.SetDefined(active_flags | Values.GetDefined());
    out.SetFlags( (reduced_flags & active_flags) | (Values.GetFlags() & ~active_flags) );
    return out;
}

bool MPIDataCommunicator::OrReduce(const bool Value, const int Root) const
{
    return ReduceDetail(Value, MPI_LOR, Root);
}

Kratos::Flags MPIDataCommunicator::OrReduce(const Kratos::Flags Values, const Kratos::Flags Mask, const int Root) const
{
    Flags::BlockType local_active_flags = Values.GetDefined() & Mask.GetDefined();
    Flags::BlockType active_flags = local_active_flags;
    ReduceDetail(local_active_flags, active_flags, MPI_BOR, Root);

    Flags::BlockType flags = Values.GetFlags();
    Flags::BlockType reduced_flags = flags;
    ReduceDetail(flags, reduced_flags, MPI_BOR, Root);

    Flags out;
    out.SetDefined(active_flags | Values.GetDefined());
    out.SetFlags( (reduced_flags & active_flags) | (Values.GetFlags() & ~active_flags) );
    return out;
}

// Allreduce operations

array_1d<double,3> MPIDataCommunicator::SumAll(const array_1d<double,3>& rLocalValue) const
{
    array_1d<double,3> global_value(rLocalValue);
    AllReduceDetail(rLocalValue,global_value,MPI_SUM);
    return global_value;
}

array_1d<double,3> MPIDataCommunicator::MinAll(const array_1d<double,3>& rLocalValue) const
{
    array_1d<double,3> global_value(rLocalValue);
    AllReduceDetail(rLocalValue,global_value,MPI_MIN);
    return global_value;
}

array_1d<double,3> MPIDataCommunicator::MaxAll(const array_1d<double,3>& rLocalValue) const
{
    array_1d<double,3> global_value(rLocalValue);
    AllReduceDetail(rLocalValue,global_value,MPI_MAX);
    return global_value;
}

bool MPIDataCommunicator::AndReduceAll(const bool Value) const
{
    return AllReduceDetail(Value, MPI_LAND);
}

Kratos::Flags MPIDataCommunicator::AndReduceAll(const Kratos::Flags Values, const Kratos::Flags Mask) const
{
    Flags::BlockType local_active_flags = Values.GetDefined() & Mask.GetDefined();
    Flags::BlockType active_flags;
    AllReduceDetail(local_active_flags, active_flags, MPI_BOR);

    Flags::BlockType flags = Values.GetFlags();
    Flags::BlockType reduced_flags;
    AllReduceDetail(flags, reduced_flags, MPI_BAND);

    Flags out;
    out.SetDefined(active_flags | Values.GetDefined());
    out.SetFlags( (reduced_flags & active_flags) | (Values.GetFlags() & ~active_flags) );
    return out;
}

bool MPIDataCommunicator::OrReduceAll(const bool Value) const
{
    return AllReduceDetail(Value, MPI_LOR);
}

Kratos::Flags MPIDataCommunicator::OrReduceAll(const Kratos::Flags Values, const Kratos::Flags Mask) const
{
    Flags::BlockType local_active_flags = Values.GetDefined() & Mask.GetDefined();
    Flags::BlockType active_flags = local_active_flags;
    AllReduceDetail(local_active_flags, active_flags, MPI_BOR);

    Flags::BlockType flags = Values.GetFlags();
    Flags::BlockType reduced_flags = flags;
    AllReduceDetail(flags, reduced_flags, MPI_BOR);

    Flags out;
    out.SetDefined(active_flags | Values.GetDefined());
    out.SetFlags( (reduced_flags & active_flags) | (Values.GetFlags() & ~active_flags) );
    return out;
}

// Access

MPI_Comm MPIDataCommunicator::GetMPICommunicator(const DataCommunicator& rDataCommunicator)
{
    if (rDataCommunicator.IsDistributed())
    {
        const MPIDataCommunicator& r_mpi_data_comm = static_cast<const MPIDataCommunicator&>(rDataCommunicator);
        return r_mpi_data_comm.mComm;
    }
    else
    {
        return MPI_COMM_SELF;
    }
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

bool MPIDataCommunicator::IsDefinedOnThisRank() const
{
    return !this->IsNullOnThisRank();
}

bool MPIDataCommunicator::IsNullOnThisRank() const
{
    return mComm == MPI_COMM_NULL;
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
{
    rOStream << "This is rank " << Rank() << " of " << Size() << "." << std::endl;
}

// Error checking

void MPIDataCommunicator::CheckMPIErrorCode(const int ierr, const std::string& MPICallName) const
{
    KRATOS_ERROR_IF_NOT(ierr == MPI_SUCCESS) << MPICallName << " failed with error code " << ierr << "." << std::endl;
}

// Protected interface reimplementing DataCommunicator calls

KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_IMPLEMENTATION_FOR_TYPE(int)
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_IMPLEMENTATION_FOR_TYPE(unsigned int)
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_IMPLEMENTATION_FOR_TYPE(long unsigned int)
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_IMPLEMENTATION_FOR_TYPE(double)

// Broadcast operations

void MPIDataCommunicator::BroadcastImpl(std::string& rBroadcastValues, const int SourceRank) const
{
    BroadcastDetail(rBroadcastValues, SourceRank);
}

// Sendrecv operations

std::string MPIDataCommunicator::SendRecvImpl(
    const std::string& rSendValues, const int SendDestination, const int SendTag,
    const int RecvSource, const int RecvTag) const
{
    int send_size = rSendValues.size();
    int recv_size;
    SendRecvDetail(send_size, SendDestination, SendTag, recv_size, RecvSource, RecvTag);

    std::string recv_values;
    recv_values.resize(recv_size);
    SendRecvDetail(rSendValues,SendDestination, SendTag, recv_values, RecvSource, RecvTag);
    return recv_values;
}

void MPIDataCommunicator::SendRecvImpl(
        const std::string& rSendValues, const int SendDestination, const int SendTag,
        std::string& rRecvValues, const int RecvSource, const int RecvTag) const
{
    SendRecvDetail(rSendValues,SendDestination,SendTag,rRecvValues,RecvSource,RecvTag);
}

void MPIDataCommunicator::SendImpl(const std::string& rSendValues, const int SendDestination, const int SendTag) const
{
    SendDetail(rSendValues, SendDestination, SendTag);
}

void MPIDataCommunicator::RecvImpl(std::string& rRecvValues, const int RecvSource, const int RecvTag) const
{
    RecvDetail(rRecvValues, RecvSource, RecvTag);
}

// Implementation details of MPI calls

template<class TDataType> void MPIDataCommunicator::ReduceDetail(
    const TDataType& rLocalValues, TDataType& rReducedValues,
    MPI_Op Operation, const int Root) const
{
    #ifdef KRATOS_DEBUG
    KRATOS_ERROR_IF_NOT(ErrorIfFalseOnAnyRank(IsValidRank(Root)))
    << "In call to MPI_Reduce: " << Root << " is not a valid rank." << std::endl;
    const int local_size = MPIMessageSize(rLocalValues);
    const int reduced_size = MPIMessageSize(rReducedValues);
    KRATOS_ERROR_IF_NOT(IsEqualOnAllRanks(local_size))
    << "Input error in call to MPI_Reduce: "
    << "There should be the same amount of local values to send from each rank." << std::endl;
    KRATOS_ERROR_IF(BroadcastErrorIfTrue(local_size != reduced_size,Root))
    << "Input error in call to MPI_Reduce for rank " << Root << ": "
    << "Sending " << local_size << " values " << "but receiving " << reduced_size << " values." << std::endl;
    #endif // KRATOS_DEBUG

    int ierr = MPI_Reduce(
        MPIBuffer(rLocalValues), MPIBuffer(rReducedValues),
        MPIMessageSize(rLocalValues), MPIDatatype(rLocalValues),
        Operation, Root, mComm);
    CheckMPIErrorCode(ierr, "MPI_Reduce");
}

template<class TDataType> TDataType MPIDataCommunicator::ReduceDetail(
    const TDataType& rLocalValues, MPI_Op Operation, const int Root) const
{
    TDataType global_values(rLocalValues);
    ReduceDetail(rLocalValues, global_values, Operation, Root);
    return global_values;
}

template<class TDataType>
std::vector<TDataType> MPIDataCommunicator::ReduceDetailVector(
    const std::vector<TDataType>& rLocalValues,
    MPI_Op Operation,
    const int Root) const
{
    std::vector<TDataType> reduced_values;
    if (Rank() == Root)
    {
        reduced_values.resize(rLocalValues.size());
    }
    ReduceDetail(rLocalValues, reduced_values, Operation, Root);
    return reduced_values;
}


template<class TDataType> void MPIDataCommunicator::AllReduceDetail(
    const TDataType& rLocalValues, TDataType& rReducedValues,
    MPI_Op Operation) const
{
    #ifdef KRATOS_DEBUG
    const int local_size = MPIMessageSize(rLocalValues);
    const int reduced_size = MPIMessageSize(rReducedValues);
    KRATOS_ERROR_IF_NOT(IsEqualOnAllRanks(local_size))
    << "Input error in call to MPI_Allreduce: "
    << "There should be the same amount of local values to send from each rank." << std::endl;
    KRATOS_ERROR_IF(ErrorIfTrueOnAnyRank(local_size != reduced_size))
    << "Input error in call to MPI_Allreduce for rank " << Rank() << ": "
    << "Sending " << local_size << " values " << "but receiving " << reduced_size << " values." << std::endl;
    #endif // KRATOS_DEBUG

    int ierr = MPI_Allreduce(
        MPIBuffer(rLocalValues), MPIBuffer(rReducedValues),
        MPIMessageSize(rLocalValues), MPIDatatype(rLocalValues),
        Operation, mComm);
    CheckMPIErrorCode(ierr, "MPI_Allreduce");
}

template<class TDataType> TDataType MPIDataCommunicator::AllReduceDetail(
    const TDataType& rLocalValues, MPI_Op Operation) const
{
    TDataType global_values(rLocalValues);
    AllReduceDetail(rLocalValues, global_values, Operation);
    return global_values;
}

template<class TDataType>
std::vector<TDataType> MPIDataCommunicator::AllReduceDetailVector(
    const std::vector<TDataType>& rLocalValues,
    MPI_Op Operation) const
{
    std::vector<TDataType> reduced_values(rLocalValues.size());
    AllReduceDetail(rLocalValues, reduced_values, Operation);
    return reduced_values;
}

template<class TDataType> void MPIDataCommunicator::ScanDetail(
    const TDataType& rLocalValues, TDataType& rReducedValues,
    MPI_Op Operation) const
{
    #ifdef KRATOS_DEBUG
    const int local_size = MPIMessageSize(rLocalValues);
    const int reduced_size = MPIMessageSize(rReducedValues);
    KRATOS_ERROR_IF_NOT(IsEqualOnAllRanks(local_size))
    << "Input error in call to MPI_Scan: "
    << "There should be the same amount of local values to send from each rank." << std::endl;
    KRATOS_ERROR_IF(ErrorIfTrueOnAnyRank(local_size != reduced_size))
    << "Input error in call to MPI_Scan for rank " << Rank() << ": "
    << "Sending " << local_size << " values " << "but receiving " << reduced_size << " values." << std::endl;
    #endif // KRATOS_DEBUG

    int ierr = MPI_Scan(
        MPIBuffer(rLocalValues), MPIBuffer(rReducedValues),
        MPIMessageSize(rLocalValues), MPIDatatype(rLocalValues),
        Operation, mComm);
    CheckMPIErrorCode(ierr, "MPI_Scan");
}

template<class TDataType> TDataType MPIDataCommunicator::ScanDetail(
    const TDataType LocalValues, MPI_Op Operation) const
{
    TDataType global_value;
    ScanDetail(LocalValues, global_value, MPI_SUM);
    return global_value;
}

template<class TDataType> std::vector<TDataType> MPIDataCommunicator::ScanDetail(
    const std::vector<TDataType>& rLocalValues, MPI_Op Operation) const
{
    std::vector<TDataType> global_values(rLocalValues.size());
    ScanDetail(rLocalValues, global_values, MPI_SUM);
    return global_values;
}

template<class TDataType> void MPIDataCommunicator::SendRecvDetail(
    const TDataType& rSendMessage, const int SendDestination, const int SendTag,
    TDataType& rRecvMessage, const int RecvSource, const int RecvTag) const
{
    int ierr = MPI_Sendrecv(
        MPIBuffer(rSendMessage), MPIMessageSize(rSendMessage),
        MPIDatatype(rSendMessage), SendDestination, SendTag,
        MPIBuffer(rRecvMessage), MPIMessageSize(rRecvMessage),
        MPIDatatype(rRecvMessage), RecvSource, RecvTag,
        mComm, MPI_STATUS_IGNORE);
    CheckMPIErrorCode(ierr, "MPI_Sendrecv");
}

template<class TDataType> TDataType MPIDataCommunicator::SendRecvDetail(
    const TDataType& rSendMessage,
    const int SendDestination, const int SendTag,
    const int RecvSource, const int RecvTag) const
{
    TDataType recv_values;
    SendRecvDetail(rSendMessage,SendDestination, SendTag ,recv_values,RecvSource, RecvTag);
    return recv_values;
}

template<class TDataType> std::vector<TDataType> MPIDataCommunicator::SendRecvDetail(
    const std::vector<TDataType>& rSendMessage,
    const int SendDestination, const int SendTag,
    const int RecvSource, const int RecvTag) const
{
    int send_size = rSendMessage.size();
    int recv_size;
    SendRecvDetail(send_size, SendDestination, SendTag, recv_size, RecvSource, RecvTag);

    std::vector<TDataType> recv_values(recv_size);
    SendRecvDetail(rSendMessage,SendDestination, SendTag ,recv_values,RecvSource, RecvTag);
    return recv_values;
}

template<class TDataType> void MPIDataCommunicator::SendDetail(
    const TDataType& rSendValues, const int SendDestination, const int SendTag) const
{
    int ierr = MPI_Send(MPIBuffer(rSendValues), MPIMessageSize(rSendValues),
        MPIDatatype(rSendValues), SendDestination, SendTag, mComm);
    CheckMPIErrorCode(ierr, "MPI_Send");
}

template<class TDataType> void MPIDataCommunicator::RecvDetail(
    TDataType& rRecvValues, const int RecvSource, const int RecvTag) const
{
    MPI_Status status;
    int ierr = MPI_Probe(RecvSource,RecvTag,mComm,&status);
    CheckMPIErrorCode(ierr, "MPI_Probe");

    int recv_size;
    ierr = MPI_Get_count(&status, MPIDatatype(rRecvValues), &recv_size);

    if (rRecvValues.size() != (unsigned int)recv_size) rRecvValues.resize(recv_size);

    ierr = MPI_Recv(MPIBuffer(rRecvValues), recv_size, MPIDatatype(rRecvValues),
        RecvSource, RecvTag, mComm, MPI_STATUS_IGNORE);
    CheckMPIErrorCode(ierr, "MPI_Recv");
}

template<class TDataType> void MPIDataCommunicator::BroadcastDetail(
    TDataType& rBuffer, const int SourceRank) const
{
    #ifdef KRATOS_DEBUG
    KRATOS_ERROR_IF_NOT(ErrorIfFalseOnAnyRank(IsValidRank(SourceRank)))
    << "In call to MPI_Bcast: " << SourceRank << " is not a valid rank." << std::endl;
    KRATOS_ERROR_IF_NOT(IsEqualOnAllRanks(MPIMessageSize(rBuffer)))
    << "Input error in call to MPI_Bcast: "
    << "The buffer does not have the same size on all ranks." << std::endl;
    #endif

    int ierr = MPI_Bcast(
        MPIBuffer(rBuffer), MPIMessageSize(rBuffer),
        MPIDatatype(rBuffer), SourceRank, mComm);
    CheckMPIErrorCode(ierr, "MPI_Bcast");
}

template<class TSendDataType, class TRecvDataType> void MPIDataCommunicator::ScatterDetail(
    const TSendDataType& rSendValues, TRecvDataType& rRecvValues, const int SourceRank) const
{
    #ifdef KRATOS_DEBUG
    KRATOS_ERROR_IF_NOT(ErrorIfFalseOnAnyRank(IsValidRank(SourceRank)))
    << "In call to MPI_Scatter: " << SourceRank << " is not a valid rank." << std::endl;
    const int send_size = MPIMessageSize(rSendValues);
    const int recv_size = MPIMessageSize(rRecvValues);
    KRATOS_ERROR_IF_NOT(IsEqualOnAllRanks(recv_size))
    << "Input error in call to MPI_Scatter: "
    << "The destination buffer does not have the same size on all ranks." << std::endl;
    KRATOS_ERROR_IF(BroadcastErrorIfTrue(send_size != recv_size*Size(),SourceRank))
    << "Input error in call to MPI_Scatter for rank " << SourceRank << ": "
    << "Sending " << send_size << " values " << "but receiving " << recv_size << " values ("
    << recv_size * Size() << " values to send expected)." << std::endl;
    #endif // KRATOS_DEBUG

    const int sends_per_rank = MPIMessageSize(rRecvValues);
    int ierr = MPI_Scatter(
        MPIBuffer(rSendValues), sends_per_rank, MPIDatatype(rSendValues),
        MPIBuffer(rRecvValues), sends_per_rank, MPIDatatype(rRecvValues),
        SourceRank, mComm);
    CheckMPIErrorCode(ierr, "MPI_Scatter");
}

template<class TDataType> std::vector<TDataType> MPIDataCommunicator::ScatterDetail(
        const std::vector<TDataType>& rSendValues, const int SourceRank) const
{
    const int send_size = rSendValues.size();
    const int world_size = Size();
    KRATOS_ERROR_IF_NOT( send_size % world_size == 0 )
    << "In call to MPI_Scatter: A message of size " << send_size
    << " cannot be evenly distributed amongst " << world_size << " ranks." << std::endl;
    int message_size = send_size / world_size;

    Broadcast(message_size, SourceRank);

    std::vector<TDataType> message(message_size);
    ScatterDetail(rSendValues, message, SourceRank);
    return message;
}

template<class TDataType> void MPIDataCommunicator::ScattervDetail(
        const TDataType& rSendValues, const std::vector<int>& rSendCounts, const std::vector<int>& rSendOffsets,
        TDataType& rRecvValues, const int SourceRank) const
{
    #ifdef KRATOS_DEBUG
    ValidateScattervInput(rSendValues, rSendCounts, rSendOffsets, rRecvValues, SourceRank);
    #endif

    int ierr = MPI_Scatterv(
        MPIBuffer(rSendValues), rSendCounts.data(), rSendOffsets.data(), MPIDatatype(rSendValues),
        MPIBuffer(rRecvValues), MPIMessageSize(rRecvValues), MPIDatatype(rRecvValues),
        SourceRank, mComm);
    CheckMPIErrorCode(ierr, "MPI_Scatterv");
}

template<class TDataType> std::vector<TDataType> MPIDataCommunicator::ScattervDetail(
    const std::vector<std::vector<TDataType>>& rSendValues,const int SourceRank) const
{
    std::vector<TDataType> message;
    std::vector<int> message_lengths;
    std::vector<int> message_offsets;
    std::vector<TDataType> result;
    PrepareScattervBuffers(
        rSendValues, message, message_lengths, message_offsets, result, SourceRank);

    ScattervDetail(message, message_lengths, message_offsets, result, SourceRank);
    return result;
}

template<class TSendDataType, class TRecvDataType> void MPIDataCommunicator::GatherDetail(
    const TSendDataType& rSendValues, TRecvDataType& rRecvValues, const int RecvRank) const
{
    #ifdef KRATOS_DEBUG
    KRATOS_ERROR_IF_NOT(ErrorIfFalseOnAnyRank(IsValidRank(RecvRank)))
    << "In call to MPI_Gather: " << RecvRank << " is not a valid rank." << std::endl;
    const int send_size = MPIMessageSize(rSendValues);
    const int recv_size = MPIMessageSize(rRecvValues);
    KRATOS_ERROR_IF_NOT(IsEqualOnAllRanks(send_size))
    << "Input error in call to MPI_Gather: "
    << "There should be the same amount of local values to send from each rank." << std::endl;
    KRATOS_ERROR_IF(BroadcastErrorIfTrue(send_size*Size() != recv_size,RecvRank))
    << "Input error in call to MPI_Gather for rank " << RecvRank << ": "
    << "Sending " << send_size << " values " << "but receiving " << recv_size << " values ("
    << send_size * Size() << " values to receive expected)." << std::endl;
    #endif // KRATOS_DEBUG

    const int sends_per_rank = MPIMessageSize(rSendValues);
    int ierr = MPI_Gather(
        MPIBuffer(rSendValues), sends_per_rank, MPIDatatype(rSendValues),
        MPIBuffer(rRecvValues), sends_per_rank, MPIDatatype(rRecvValues),
        RecvRank, mComm);
    CheckMPIErrorCode(ierr, "MPI_Gather");
}

template<class TDataType> std::vector<TDataType> MPIDataCommunicator::GatherDetail(
    const std::vector<TDataType>& rSendValues, const int DestinationRank) const
{
    int message_size = rSendValues.size();
    std::vector<TDataType> gathered_values;
    if (Rank() == DestinationRank)
    {
        gathered_values.resize(message_size*Size());
    }
    GatherDetail(rSendValues, gathered_values, DestinationRank);
    return gathered_values;
}

template<class TDataType> void MPIDataCommunicator::GathervDetail(
    const TDataType& rSendValues, TDataType& rRecvValues,
    const std::vector<int>& rRecvCounts, const std::vector<int>& rRecvOffsets,
    const int RecvRank) const
{
    #ifdef KRATOS_DEBUG
    ValidateGathervInput(rSendValues, rRecvValues, rRecvCounts, rRecvOffsets, RecvRank);
    #endif

    int ierr = MPI_Gatherv(
        MPIBuffer(rSendValues), MPIMessageSize(rSendValues), MPIDatatype(rSendValues),
        MPIBuffer(rRecvValues), rRecvCounts.data(), rRecvOffsets.data(), MPIDatatype(rRecvValues),
        RecvRank, mComm);
    CheckMPIErrorCode(ierr, "MPI_Gatherv");
}

template<class TDataType> std::vector<std::vector<TDataType>>
MPIDataCommunicator::GathervDetail(const std::vector<TDataType>& rSendValues, const int DestinationRank) const
{
    std::vector<TDataType> message;
    std::vector<int> message_lengths;
    std::vector<int> message_offsets;
    PrepareGathervBuffers(rSendValues, message, message_lengths, message_offsets, DestinationRank);

    Gatherv(rSendValues, message, message_lengths, message_offsets, DestinationRank);

    std::vector<std::vector<TDataType>> output_message;
    PrepareGathervReturn(message, message_lengths, message_offsets, output_message, DestinationRank);
    return output_message;
}

template<class TDataType> void MPIDataCommunicator::AllGatherDetail(
    const TDataType& rSendValues, TDataType& rRecvValues) const
{
    #ifdef KRATOS_DEBUG
    const int send_size = MPIMessageSize(rSendValues);
    const int recv_size = MPIMessageSize(rRecvValues);
    KRATOS_ERROR_IF_NOT(IsEqualOnAllRanks(send_size))
    << "Input error in call to MPI_Allgather: "
    << "There should be the same amount of local values to send from each rank." << std::endl;
    KRATOS_ERROR_IF(ErrorIfTrueOnAnyRank(send_size*Size() != recv_size))
    << "Input error in call to MPI_Allgather for rank " << Rank() << ": "
    << "Sending " << send_size << " values " << "but receiving " << recv_size << " values ("
    << send_size * Size() << " values to receive expected)." << std::endl;
    #endif // KRATOS_DEBUG

    const int sends_per_rank = MPIMessageSize(rSendValues);
    int ierr = MPI_Allgather(
        MPIBuffer(rSendValues), sends_per_rank, MPIDatatype(rSendValues),
        MPIBuffer(rRecvValues), sends_per_rank, MPIDatatype(rRecvValues),
        mComm);
    CheckMPIErrorCode(ierr, "MPI_Allgather");
}

template<class TDataType> std::vector<TDataType> MPIDataCommunicator::AllGatherDetail(
    const std::vector<TDataType>& rSendValues) const
{
    std::vector<TDataType> output(rSendValues.size()*Size());
    AllGatherDetail(rSendValues, output);
    return output;
}

bool MPIDataCommunicator::BroadcastErrorIfTrue(bool Condition, const int SourceRank) const
{
    int ierr = MPI_Bcast(&Condition,1,MPI_C_BOOL,SourceRank,mComm);
    CheckMPIErrorCode(ierr, "MPI_Bcast");
    const int rank = Rank();
    KRATOS_ERROR_IF(Condition && (rank != SourceRank) )
    << "Rank " << rank << ": Stopping because of error in rank " << SourceRank << "." << std::endl;
    return Condition;
}

bool MPIDataCommunicator::BroadcastErrorIfFalse(bool Condition, const int SourceRank) const
{
    int ierr = MPI_Bcast(&Condition,1,MPI_C_BOOL,SourceRank,mComm);
    CheckMPIErrorCode(ierr, "MPI_Bcast");
    const int rank = Rank();
    KRATOS_ERROR_IF(!Condition && (rank != SourceRank) )
    << "Rank " << rank << ": Stopping because of error in rank " << SourceRank << "." << std::endl;
    return Condition;
}

bool MPIDataCommunicator::ErrorIfTrueOnAnyRank(bool Condition) const
{
    // Note: this function cannot use the helper function AllReduceDetail
    // even if it implements the same funtionality. AllReduceDetail calls
    // ErrorIfTrueOnAnyRank in debug mode for consistency checking
    // and that would result on a circular call.
    bool or_condition;
    int ierr = MPI_Allreduce(&Condition, &or_condition, 1, MPI_C_BOOL, MPI_LOR, mComm);
    CheckMPIErrorCode(ierr, "MPI_Allreduce");
    KRATOS_ERROR_IF(or_condition && !Condition)
    << "Rank " << Rank() << ": Stopping because an error was detected on a different rank." << std::endl;
    return or_condition;
}

bool MPIDataCommunicator::ErrorIfFalseOnAnyRank(bool Condition) const
{
    bool and_condition;
    int ierr = MPI_Allreduce(&Condition, &and_condition, 1, MPI_C_BOOL, MPI_LAND, mComm);
    CheckMPIErrorCode(ierr, "MPI_Allreduce");
    KRATOS_ERROR_IF(!and_condition && Condition)
    << "Rank " << Rank() << ": Stopping because an error was detected on a different rank." << std::endl;
    return and_condition;
}

bool MPIDataCommunicator::IsEqualOnAllRanks(const int LocalValue) const
{
    int local_buffer[2]{LocalValue, -LocalValue};
    int min_buffer[2]{0, 0};
    int ierr = MPI_Allreduce(&local_buffer,&min_buffer,2,MPI_INT,MPI_MIN,mComm);
    CheckMPIErrorCode(ierr, "MPI_Allreduce");
    int min_value = min_buffer[0];
    int max_value = -min_buffer[1];
    return min_value == max_value;
}

bool MPIDataCommunicator::IsValidRank(const int Rank) const
{
    return (Rank >= 0) && (Rank < Size());
}

template<class TDataType> void MPIDataCommunicator::ValidateScattervInput(
    const TDataType& rSendValues,
    const std::vector<int>& rSendCounts, const std::vector<int>& rSendOffsets,
    TDataType& rRecvValues, const int SourceRank) const
{
    KRATOS_ERROR_IF_NOT(ErrorIfFalseOnAnyRank(IsValidRank(SourceRank)))
    << "In call to MPI_Scatterv: " << SourceRank << " is not a valid rank." << std::endl;

    // All ranks expect a message of the correct size
    int expected_size = 0;
    const int available_recv_size = MPIMessageSize(rRecvValues);
    int ierr = MPI_Scatter(rSendCounts.data(), 1, MPI_INT, &expected_size, 1, MPI_INT, SourceRank, mComm);
    CheckMPIErrorCode(ierr, "MPI_Scatter");
    KRATOS_ERROR_IF(ErrorIfTrueOnAnyRank(expected_size != available_recv_size))
    << "Input error in call to MPI_Scatterv for rank " << Rank() << ": "
    << "This rank will receive " << expected_size << " values but the receive buffer has size "
    << available_recv_size << "." << std::endl;

    // Message size is not smaller than total expected size (can only check for too small, since the source message may be padded).
    int total_size = 0;
    const int message_size = MPIMessageSize(rSendValues);
    ierr = MPI_Reduce(&available_recv_size, &total_size, 1, MPI_INT, MPI_SUM, SourceRank, mComm);
    CheckMPIErrorCode(ierr, "MPI_Reduce");
    KRATOS_ERROR_IF(BroadcastErrorIfTrue(total_size > message_size, SourceRank))
    << "Input error in call to MPI_Scatterv for rank " << SourceRank << ": "
    << "The sent message contains " << message_size << " values, but " << available_recv_size
    << " values are expected in total across all ranks." << std::endl;

    // No overflow in sent buffer.
    std::stringstream message;
    bool failed = false;
    if (Rank() == SourceRank)
    {
        for (int i = 0; i < Size(); i++)
        {
            if(rSendOffsets[i]+rSendCounts[i] > message_size) {
                message
                << "Input error in call to MPI_Scatterv for rank " << SourceRank << ": "
                << "Reading past sent message end when sending message for rank " << i << "." << std::endl;
                failed = true;
                break;
            }
        }
    }
    KRATOS_ERROR_IF(BroadcastErrorIfTrue(failed, SourceRank)) << message.str();
}

template<class TDataType> void MPIDataCommunicator::ValidateGathervInput(
    const TDataType& rSendValues, TDataType& rRecvValues,
    const std::vector<int>& rRecvCounts, const std::vector<int>& rRecvOffsets,
    const int RecvRank) const
{
    KRATOS_ERROR_IF_NOT(ErrorIfFalseOnAnyRank(IsValidRank(RecvRank)))
    << "In call to MPI_Gatherv: " << RecvRank << " is not a valid rank." << std::endl;

    // All ranks send a message of the correct size
    int expected_recv_size = 0;
    const int send_size = MPIMessageSize(rSendValues);
    int ierr = MPI_Scatter(rRecvCounts.data(), 1, MPI_INT, &expected_recv_size, 1, MPI_INT, RecvRank, mComm);
    CheckMPIErrorCode(ierr, "MPI_Scatter");
    KRATOS_ERROR_IF(ErrorIfTrueOnAnyRank(send_size != expected_recv_size))
    << "Input error in call to MPI_Gatherv for rank " << Rank() << ": "
    << "This rank will send " << send_size << " values but " << RecvRank << " expects "
    << expected_recv_size << " values from it." << std::endl;

    // Message size is not larger than total expected size (can only check for too large, since the recv message may be padded).
    int total_size = 0;
    const int message_size = MPIMessageSize(rSendValues);
    const int expected_message_size = MPIMessageSize(rRecvValues);
    ierr = MPI_Reduce(&message_size, &total_size, 1, MPI_INT, MPI_SUM, RecvRank, mComm);
    CheckMPIErrorCode(ierr, "MPI_Reduce");
    KRATOS_ERROR_IF(BroadcastErrorIfTrue(total_size > expected_message_size, RecvRank))
    << "Input error in call to MPI_Gatherv for rank " << RecvRank << ": "
    << "The sent messages contain " << total_size << " values in total, but only "
    << expected_message_size << " values are expected in rank " << RecvRank << "." << std::endl;

    // No overflow in recv buffer.
    std::stringstream message;
    bool failed = false;
    if (Rank() == RecvRank)
    {
        for (int i = 0; i < Size(); i++)
        {
            if (rRecvOffsets[i]+rRecvCounts[i] > expected_message_size) {
                message
                << "Input error in call to MPI_Gatherv for rank " << RecvRank << ": "
                << "Writing past buffer end when sending message for rank " << i << "." << std::endl;
                failed = true;
                break;
            }
        }
    }
    KRATOS_ERROR_IF(BroadcastErrorIfTrue(failed, RecvRank)) << message.str();
}

template <class TDataType> void MPIDataCommunicator::PrepareScattervBuffers(
    const std::vector<std::vector<TDataType>>& rInputMessage,
    std::vector<TDataType>& rScattervMessage,
    std::vector<int>& rMessageLengths,
    std::vector<int>& rMessageDistances,
    std::vector<TDataType>& rResult,
    const int SourceRank) const
{
    if (Rank() == SourceRank)
    {
        unsigned int size = Size();
        KRATOS_ERROR_IF_NOT(rInputMessage.size() == size)
        << "Input error in call to MPI_Scatterv: Expected " << size << " vectors as input, got " << rInputMessage.size() << "." << std::endl;

        rMessageLengths.resize(size);
        rMessageDistances.resize(size);
        unsigned int message_size = 0;
        for (unsigned int i = 0; i < rInputMessage.size(); i++)
        {
            rMessageDistances[i] = message_size;
            unsigned int rank_size = rInputMessage[i].size();
            rMessageLengths[i] = rank_size;
            message_size += rank_size;
        }

        rScattervMessage.resize(message_size);
        for (unsigned int i = 0, counter = 0; i < rInputMessage.size(); i++)
        {
            for (unsigned int j = 0; j < rInputMessage[i].size(); j++, counter++)
            {
                rScattervMessage[counter] = rInputMessage[i][j];
            }
        }
    }

    int result_size;
    ScatterDetail(rMessageLengths, result_size, SourceRank);
    rResult.resize(result_size);
}


template<class TDataType> void MPIDataCommunicator::PrepareGathervBuffers(
    const std::vector<TDataType>& rGathervInput,
    std::vector<TDataType>& rGathervMessage,
    std::vector<int>& rMessageLengths,
    std::vector<int>& rMessageDistances,
    const int DestinationRank) const
{
    int message_size_send = rGathervInput.size();
    const int rank = Rank();
    const int size = Size();
    if (rank == DestinationRank)
    {
        rMessageLengths.resize(size);
    }
    GatherDetail(message_size_send, rMessageLengths, DestinationRank);

    if (rank == DestinationRank)
    {
        rMessageDistances.resize(size);
        int message_size = 0;
        for (int i = 0; i < size; i++)
        {
            rMessageDistances[i] = message_size;
            message_size += rMessageLengths[i];
        }
        rGathervMessage.resize(message_size);
    }
}

template<class TDataType> void MPIDataCommunicator::PrepareGathervReturn(
    const std::vector<TDataType>& rGathervMessage,
    const std::vector<int>& rMessageLengths,
    const std::vector<int>& rMessageDistances,
    std::vector<std::vector<TDataType>>& rOutputMessage,
    const int DestinationRank) const
{
    const int size = Size();
    rOutputMessage.resize(size);
    if (Rank() == DestinationRank)
    {
        for (int i = 0, counter = 0; i < size; i++)
        {
            rOutputMessage[i].resize(rMessageLengths[i]);
            for (int j = 0; j < rMessageLengths[i]; j++, counter++)
            {
                rOutputMessage[i][j] = rGathervMessage[counter];
            }
        }
    }
}

// MPI_Datatype wrapper
template<class TValue> inline MPI_Datatype MPIDataCommunicator::MPIDatatype(const TValue&) const
{
    return MPIMessage<TValue>::DataType();
}

// Buffer argument deduction
template<class TContainer> inline void* MPIDataCommunicator::MPIBuffer(TContainer& rValues) const
{
    return MPIMessage<TContainer>::Buffer(rValues);
}

template<class TContainer> inline const void* MPIDataCommunicator::MPIBuffer(const TContainer& rValues) const
{
    return MPIMessage<TContainer>::Buffer(rValues);
}

// MPI message size deduction
template<class TContainer> inline int MPIDataCommunicator::MPIMessageSize(const TContainer& rValues) const
{
    return MPIMessage<TContainer>::Size(rValues);
}

}

#undef KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_REDUCE_INTERFACE_FOR_TYPE
#undef KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_ALLREDUCE_INTERFACE_FOR_TYPE
#undef KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_SCANSUM_INTERFACE_FOR_TYPE
#undef KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_SENDRECV_INTERFACE_FOR_TYPE
#undef KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_BROADCAST_INTERFACE_FOR_TYPE
#undef KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_SCATTER_INTERFACE_FOR_TYPE
#undef KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_GATHER_INTERFACE_FOR_TYPE
#undef KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_PUBLIC_INTERFACE_FOR_TYPE
#undef KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_IMPLEMENTATION_FOR_TYPE
