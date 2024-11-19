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

// System includes
#include <algorithm>

// External includes

// Project includes
#include "includes/parallel_environment.h"
#include "mpi/includes/mpi_data_communicator.h"
#include "mpi/includes/mpi_manager.h"
#include "mpi/includes/mpi_message.h"
#include "mpi/utilities/data_communicator_factory.h"

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_SYNC_SHAPE_INTERFACE_FOR_TYPE
#define KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_SYNC_SHAPE_INTERFACE_FOR_TYPE(...)       \
    bool MPIDataCommunicator::SynchronizeShape(__VA_ARGS__& rValue) const {          \
        return SynchronizeShapeDetail(rValue);                                       \
    }                                                                                \
    bool MPIDataCommunicator::SynchronizeShape(                                      \
        const __VA_ARGS__& rSendValue, const int SendDestination, const int SendTag, \
        __VA_ARGS__& rRecvValue, const int RecvSource, const int RecvTag) const {    \
        return SynchronizeShapeDetail(rSendValue, SendDestination, SendTag,          \
                                      rRecvValue, RecvSource, RecvTag);              \
    }

#endif

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_REDUCE_INTERFACE_FOR_TYPE
#define KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_REDUCE_INTERFACE_FOR_TYPE(...)                                              \
__VA_ARGS__ MPIDataCommunicator::Sum(const __VA_ARGS__& rLocalValue, const int Root) const {                            \
    __VA_ARGS__ global_values(rLocalValue);                                                                             \
    ReduceDetail(rLocalValue, global_values, MPI_SUM, Root);                                                            \
    return global_values;                                                                                               \
}                                                                                                                       \
std::vector<__VA_ARGS__> MPIDataCommunicator::Sum(const std::vector<__VA_ARGS__>& rLocalValues, const int Root) const { \
    return ReduceDetailVector(rLocalValues, MPI_SUM, Root);                                                             \
}                                                                                                                       \
void MPIDataCommunicator::Sum(                                                                                          \
    const std::vector<__VA_ARGS__>& rLocalValues, std::vector<__VA_ARGS__>& rGlobalValues, const int Root) const {      \
    ReduceDetail(rLocalValues, rGlobalValues, MPI_SUM, Root);                                                           \
}                                                                                                                       \
__VA_ARGS__ MPIDataCommunicator::Min(const __VA_ARGS__& rLocalValue, const int Root) const {                            \
    __VA_ARGS__ global_values(rLocalValue);                                                                             \
    ReduceDetail(rLocalValue, global_values, MPI_MIN, Root);                                                            \
    return global_values;                                                                                               \
}                                                                                                                       \
std::vector<__VA_ARGS__> MPIDataCommunicator::Min(const std::vector<__VA_ARGS__>& rLocalValues, const int Root) const { \
    return ReduceDetailVector(rLocalValues, MPI_MIN, Root);                                                             \
}                                                                                                                       \
void MPIDataCommunicator::Min(                                                                                          \
    const std::vector<__VA_ARGS__>& rLocalValues, std::vector<__VA_ARGS__>& rGlobalValues, const int Root) const {      \
    ReduceDetail(rLocalValues, rGlobalValues, MPI_MIN, Root);                                                           \
}                                                                                                                       \
__VA_ARGS__ MPIDataCommunicator::Max(const __VA_ARGS__& rLocalValue, const int Root) const {                            \
    __VA_ARGS__ global_values(rLocalValue);                                                                             \
    ReduceDetail(rLocalValue, global_values, MPI_MAX, Root);                                                            \
    return global_values;                                                                                               \
}                                                                                                                       \
std::vector<__VA_ARGS__> MPIDataCommunicator::Max(const std::vector<__VA_ARGS__>& rLocalValues, const int Root) const { \
    return ReduceDetailVector(rLocalValues, MPI_MAX, Root);                                                             \
}                                                                                                                       \
void MPIDataCommunicator::Max(                                                                                          \
    const std::vector<__VA_ARGS__>& rLocalValues, std::vector<__VA_ARGS__>& rGlobalValues, const int Root) const {      \
    ReduceDetail(rLocalValues, rGlobalValues, MPI_MAX, Root);                                                           \
}                                                                                                                       \

#endif

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_ALLREDUCE_INTERFACE_FOR_TYPE
#define KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_ALLREDUCE_INTERFACE_FOR_TYPE(...)                               \
__VA_ARGS__ MPIDataCommunicator::SumAll(const __VA_ARGS__& rLocalValue) const {                             \
    __VA_ARGS__ global_values(rLocalValue);                                                                 \
    AllReduceDetail(rLocalValue, global_values, MPI_SUM);                                                   \
    return global_values;                                                                                   \
}                                                                                                           \
std::vector<__VA_ARGS__> MPIDataCommunicator::SumAll(const std::vector<__VA_ARGS__>& rLocalValues) const {  \
    return AllReduceDetailVector(rLocalValues, MPI_SUM);                                                    \
}                                                                                                           \
void MPIDataCommunicator::SumAll(                                                                           \
    const std::vector<__VA_ARGS__>& rLocalValues, std::vector<__VA_ARGS__>& rGlobalValues) const {          \
    AllReduceDetail(rLocalValues, rGlobalValues, MPI_SUM);                                                  \
}                                                                                                           \
__VA_ARGS__ MPIDataCommunicator::MinAll(const __VA_ARGS__& rLocalValue) const {                             \
    __VA_ARGS__ global_values(rLocalValue);                                                                 \
    AllReduceDetail(rLocalValue, global_values, MPI_MIN);                                                   \
    return global_values;                                                                                   \
}                                                                                                           \
std::vector<__VA_ARGS__> MPIDataCommunicator::MinAll(const std::vector<__VA_ARGS__>& rLocalValues) const {  \
    return AllReduceDetailVector(rLocalValues, MPI_MIN);                                                    \
}                                                                                                           \
void MPIDataCommunicator::MinAll(                                                                           \
    const std::vector<__VA_ARGS__>& rLocalValues, std::vector<__VA_ARGS__>& rGlobalValues) const {          \
    AllReduceDetail(rLocalValues, rGlobalValues, MPI_MIN);                                                  \
}                                                                                                           \
__VA_ARGS__ MPIDataCommunicator::MaxAll(const __VA_ARGS__& rLocalValue) const {                             \
    __VA_ARGS__ global_values(rLocalValue);                                                                 \
    AllReduceDetail(rLocalValue, global_values, MPI_MAX);                                                   \
    return global_values;                                                                                   \
}                                                                                                           \
std::vector<__VA_ARGS__> MPIDataCommunicator::MaxAll(const std::vector<__VA_ARGS__>& rLocalValues) const {  \
    return AllReduceDetailVector(rLocalValues, MPI_MAX);                                                    \
}                                                                                                           \
void MPIDataCommunicator::MaxAll(                                                                           \
    const std::vector<__VA_ARGS__>& rLocalValues, std::vector<__VA_ARGS__>& rGlobalValues) const {          \
    AllReduceDetail(rLocalValues, rGlobalValues, MPI_MAX);                                                  \
}
#endif

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_ALLREDUCE_LOC_IMPLEMENTATION_FOR_TYPE
#define KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_ALLREDUCE_LOC_IMPLEMENTATION_FOR_TYPE(...)                      \
std::pair<__VA_ARGS__, int> MPIDataCommunicator::MinLocAll(const __VA_ARGS__& rLocalValue) const {          \
    std::pair<__VA_ARGS__, int> local_values({rLocalValue, Rank()});                                        \
    return AllReduceDetailWithLocation(local_values, MPI_MINLOC);                                           \
}                                                                                                           \
std::pair<__VA_ARGS__, int> MPIDataCommunicator::MaxLocAll(const __VA_ARGS__& rLocalValue) const {          \
    std::pair<__VA_ARGS__, int> local_values({rLocalValue, Rank()});                                        \
    return AllReduceDetailWithLocation(local_values, MPI_MAXLOC);                                           \
}
#endif

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_SCANSUM_INTERFACE_FOR_TYPE
#define KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_SCANSUM_INTERFACE_FOR_TYPE(...)                                 \
__VA_ARGS__ MPIDataCommunicator::ScanSum(const __VA_ARGS__& LocalValue) const {                             \
    __VA_ARGS__ partial_sums(LocalValue);                                                                   \
    ScanDetail(LocalValue, partial_sums, MPI_SUM);                                                          \
    return partial_sums;                                                                                    \
}                                                                                                           \
std::vector<__VA_ARGS__> MPIDataCommunicator::ScanSum(const std::vector<__VA_ARGS__>& rLocalValues) const { \
    return ScanDetail(rLocalValues, MPI_SUM);                                                               \
}                                                                                                           \
void MPIDataCommunicator::ScanSum(                                                                          \
    const std::vector<__VA_ARGS__>& rLocalValues, std::vector<__VA_ARGS__>& rPartialSums) const {           \
    ScanDetail(rLocalValues,rPartialSums,MPI_SUM);                                                          \
}                                                                                                           \

#endif

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_SENDRECV_INTERFACE_FOR_TYPE
#define KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_SENDRECV_INTERFACE_FOR_TYPE(...)                    \
__VA_ARGS__ MPIDataCommunicator::SendRecvImpl(                                                  \
    const __VA_ARGS__& SendValue,                                                               \
    const int SendDestination, const int SendTag,                                               \
    const int RecvSource, const int RecvTag) const {                                            \
    return SendRecvDetail(SendValue, SendDestination, SendTag, RecvSource, RecvTag);            \
}                                                                                               \
std::vector<__VA_ARGS__> MPIDataCommunicator::SendRecvImpl(                                     \
    const std::vector<__VA_ARGS__>& rSendValues,                                                \
    const int SendDestination, const int SendTag,                                               \
    const int RecvSource, const int RecvTag) const {                                            \
    return SendRecvDetail(rSendValues, SendDestination, SendTag, RecvSource, RecvTag);          \
}                                                                                               \
void MPIDataCommunicator::SendRecvImpl(                                                         \
    const std::vector<__VA_ARGS__>& rSendValues, const int SendDestination, const int SendTag,  \
    std::vector<__VA_ARGS__>& rRecvValues, const int RecvSource, const int RecvTag) const {     \
    SendRecvDetail(rSendValues,SendDestination,SendTag,rRecvValues,RecvSource,RecvTag);         \
}                                                                                               \
void MPIDataCommunicator::SendRecvImpl(                                                         \
    const __VA_ARGS__& SendValue, const int SendDestination, const int SendTag,                 \
    __VA_ARGS__& rRecvValue, const int RecvSource, const int RecvTag) const {                   \
    SendRecvDetail(SendValue,SendDestination,SendTag,rRecvValue,RecvSource,RecvTag);            \
}                                                                                               \
void MPIDataCommunicator::SendImpl(const __VA_ARGS__& rSendValues,                              \
    const int SendDestination, const int SendTag) const {                                       \
    std::vector<__VA_ARGS__> send{rSendValues};                                                 \
    SendDetail(send, SendDestination, SendTag);                                                 \
}                                                                                               \
void MPIDataCommunicator::SendImpl(const std::vector<__VA_ARGS__>& rSendValues,                 \
    const int SendDestination, const int SendTag) const {                                       \
    SendDetail(rSendValues, SendDestination, SendTag);                                          \
}                                                                                               \
void MPIDataCommunicator::RecvImpl(__VA_ARGS__& rRecvValues,                                    \
    const int RecvSource, const int RecvTag) const {                                            \
    std::vector<__VA_ARGS__> recv(1);                                                           \
    RecvDetail(recv, RecvSource, RecvTag);                                                      \
    rRecvValues = recv[0];                                                                      \
}                                                                                               \
void MPIDataCommunicator::RecvImpl(std::vector<__VA_ARGS__>& rRecvValues,                       \
    const int RecvSource, const int RecvTag) const {                                            \
    RecvDetail(rRecvValues, RecvSource, RecvTag);                                               \
}                                                                                               \

#endif

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_BROADCAST_INTERFACE_FOR_TYPE
#define KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_BROADCAST_INTERFACE_FOR_TYPE(...)                               \
void MPIDataCommunicator::BroadcastImpl(__VA_ARGS__& rBuffer, const int SourceRank) const {                 \
    BroadcastDetail(rBuffer,SourceRank);                                                                    \
}                                                                                                           \
void MPIDataCommunicator::BroadcastImpl(std::vector<__VA_ARGS__>& rBuffer, const int SourceRank) const {    \
    BroadcastDetail(rBuffer,SourceRank);                                                                    \
}                                                                                                           \

#endif

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_SCATTER_INTERFACE_FOR_TYPE
#define KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_SCATTER_INTERFACE_FOR_TYPE(...)                 \
std::vector<__VA_ARGS__> MPIDataCommunicator::Scatter(                                      \
    const std::vector<__VA_ARGS__>& rSendValues, const int SourceRank) const {              \
    return ScatterDetail(rSendValues, SourceRank);                                          \
}                                                                                           \
void MPIDataCommunicator::Scatter(                                                          \
    const std::vector<__VA_ARGS__>& rSendValues, std::vector<__VA_ARGS__>& rRecvValues,     \
    const int SourceRank) const {                                                           \
    ScatterDetail(rSendValues,rRecvValues,SourceRank);                                      \
}                                                                                           \
std::vector<__VA_ARGS__> MPIDataCommunicator::Scatterv(                                     \
    const std::vector<std::vector<__VA_ARGS__>>& rSendValues, const int SourceRank) const { \
    return ScattervDetail(rSendValues, SourceRank);                                         \
}                                                                                           \
void MPIDataCommunicator::Scatterv(                                                         \
    const std::vector<__VA_ARGS__>& rSendValues, const std::vector<int>& rSendCounts,       \
    const std::vector<int>& rSendOffsets, std::vector<__VA_ARGS__>& rRecvValues,            \
    const int SourceRank) const {                                                           \
    ScattervDetail(rSendValues,rSendCounts,rSendOffsets,rRecvValues,SourceRank);            \
}                                                                                           \

#endif

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_GATHER_INTERFACE_FOR_TYPE
#define KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_GATHER_INTERFACE_FOR_TYPE(...)                      \
std::vector<__VA_ARGS__> MPIDataCommunicator::Gather(                                           \
    const std::vector<__VA_ARGS__>& rSendValues, const int DestinationRank) const {             \
    return GatherDetail(rSendValues, DestinationRank);                                          \
}                                                                                               \
void MPIDataCommunicator::Gather(                                                               \
    const std::vector<__VA_ARGS__>& rSendValues, std::vector<__VA_ARGS__>& rRecvValues,         \
    const int DestinationRank) const {                                                          \
    GatherDetail(rSendValues, rRecvValues, DestinationRank);                                    \
}                                                                                               \
std::vector<std::vector<__VA_ARGS__>> MPIDataCommunicator::Gatherv(                             \
    const std::vector<__VA_ARGS__>& rSendValues, const int DestinationRank) const {             \
    return GathervDetail(rSendValues, DestinationRank);                                         \
}                                                                                               \
void MPIDataCommunicator::Gatherv(                                                              \
    const std::vector<__VA_ARGS__>& rSendValues, std::vector<__VA_ARGS__>& rRecvValues,         \
    const std::vector<int>& rRecvCounts, const std::vector<int>& rRecvOffsets,                  \
    const int DestinationRank) const {                                                          \
    GathervDetail(rSendValues,rRecvValues,rRecvCounts,rRecvOffsets,DestinationRank);            \
}                                                                                               \
std::vector<__VA_ARGS__> MPIDataCommunicator::AllGather(                                        \
    const std::vector<__VA_ARGS__>& rSendValues) const {                                        \
    return AllGatherDetail(rSendValues);                                                        \
}                                                                                               \
void MPIDataCommunicator::AllGather(                                                            \
    const std::vector<__VA_ARGS__>& rSendValues, std::vector<__VA_ARGS__>& rRecvValues) const { \
    AllGatherDetail(rSendValues,rRecvValues);                                                   \
}                                                                                               \
std::vector<std::vector<__VA_ARGS__>> MPIDataCommunicator::AllGatherv(                          \
    const std::vector<__VA_ARGS__>& rSendValues) const {                                        \
    return AllGathervDetail(rSendValues);                                                       \
}                                                                                               \
void MPIDataCommunicator::AllGatherv(const std::vector<__VA_ARGS__>& rSendValues,               \
    std::vector<__VA_ARGS__>& rRecvValues, const std::vector<int>& rRecvCounts,                 \
    const std::vector<int>& rRecvOffsets) const {                                               \
    AllGathervDetail(rSendValues,rRecvValues,rRecvCounts,rRecvOffsets);                         \
}
#endif

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_PUBLIC_INTERFACE_FOR_TYPE
#define KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_PUBLIC_INTERFACE_FOR_TYPE(...)      \
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_REDUCE_INTERFACE_FOR_TYPE(__VA_ARGS__)      \
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_ALLREDUCE_INTERFACE_FOR_TYPE(__VA_ARGS__)   \
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_SCANSUM_INTERFACE_FOR_TYPE(__VA_ARGS__)     \
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_SCATTER_INTERFACE_FOR_TYPE(__VA_ARGS__)     \
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_GATHER_INTERFACE_FOR_TYPE(__VA_ARGS__)      \
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_SYNC_SHAPE_INTERFACE_FOR_TYPE(__VA_ARGS__)  \

#endif

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_IMPLEMENTATION_FOR_TYPE
#define KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_IMPLEMENTATION_FOR_TYPE(...)        \
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_SENDRECV_INTERFACE_FOR_TYPE(__VA_ARGS__)    \
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_BROADCAST_INTERFACE_FOR_TYPE(__VA_ARGS__)   \

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
    const int ierr = MPI_Barrier(mComm);
    CheckMPIErrorCode(ierr,"MPI_Barrier");
}

// Complete interface for basic types

KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_PUBLIC_INTERFACE_FOR_TYPE(char)
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_PUBLIC_INTERFACE_FOR_TYPE(int)
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_PUBLIC_INTERFACE_FOR_TYPE(unsigned int)
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_PUBLIC_INTERFACE_FOR_TYPE(long unsigned int)
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_PUBLIC_INTERFACE_FOR_TYPE(double)
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_PUBLIC_INTERFACE_FOR_TYPE(array_1d<double, 3>)
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_PUBLIC_INTERFACE_FOR_TYPE(array_1d<double, 4>)
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_PUBLIC_INTERFACE_FOR_TYPE(array_1d<double, 6>)
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_PUBLIC_INTERFACE_FOR_TYPE(array_1d<double, 9>)
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_PUBLIC_INTERFACE_FOR_TYPE(Vector)
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_PUBLIC_INTERFACE_FOR_TYPE(Matrix)

// Reduce operations

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
    const int ierr = MPI_Comm_rank(mComm, &rank);
    CheckMPIErrorCode(ierr, "MPI_Comm_rank");
    return rank;
}

int MPIDataCommunicator::Size() const
{
    int size;
    const int ierr = MPI_Comm_size(mComm, &size);
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

const DataCommunicator& MPIDataCommunicator::GetSubDataCommunicator(
    const std::vector<int>& rRanks,
    const std::string& rNewCommunicatorName
    ) const
{
    // Initial check
    const int rank = Rank();
    const int total_size = Size();
    KRATOS_ERROR_IF_NOT(rRanks.size() <= static_cast<std::size_t>(total_size)) << "Inconsistency between the communicator total world size: " << total_size << " and the number of ranks required: " << rRanks.size() << std::endl;

    // Retrieve data communicator
    const DataCommunicator& r_data_communicator = ParallelEnvironment::HasDataCommunicator(rNewCommunicatorName) ? ParallelEnvironment::GetDataCommunicator(rNewCommunicatorName) : DataCommunicatorFactory::CreateFromRanksAndRegister(*this, rRanks, rNewCommunicatorName);
    const auto it_find = std::find(rRanks.begin(), rRanks.end(), rank);
    if (it_find != rRanks.end()) {
        KRATOS_ERROR_IF_NOT(r_data_communicator.IsDefinedOnThisRank()) << "The rank " << rank << " does not participate in the existing data communicator " << rNewCommunicatorName  << " despite being in the provided rank list" << std::endl;
        const std::size_t world_size = static_cast<std::size_t >(r_data_communicator.Size());
        KRATOS_ERROR_IF_NOT(rRanks.size() == world_size) << "Inconsistency between the communicator world size: " << world_size << " and the number of ranks required: " << rRanks.size() << std::endl;
        std::size_t number_of_active_ranks = 1;
        number_of_active_ranks = r_data_communicator.SumAll(number_of_active_ranks);
        KRATOS_ERROR_IF_NOT(number_of_active_ranks == rRanks.size()) << "Inconsistency between the number of active ranks: " << number_of_active_ranks << " and the number of ranks required: " << rRanks.size() << std::endl;
    } else {
        KRATOS_ERROR_IF_NOT(r_data_communicator.IsNullOnThisRank()) << "The rank " << rank << " participates in the existing data communicator " << rNewCommunicatorName << " despite not being in the provided rank list" << std::endl;
    }
    return r_data_communicator;
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

KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_IMPLEMENTATION_FOR_TYPE(char)
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_IMPLEMENTATION_FOR_TYPE(int)
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_IMPLEMENTATION_FOR_TYPE(unsigned int)
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_IMPLEMENTATION_FOR_TYPE(long unsigned int)
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_IMPLEMENTATION_FOR_TYPE(double)
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_IMPLEMENTATION_FOR_TYPE(array_1d<double, 3>)
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_IMPLEMENTATION_FOR_TYPE(array_1d<double, 4>)
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_IMPLEMENTATION_FOR_TYPE(array_1d<double, 6>)
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_IMPLEMENTATION_FOR_TYPE(array_1d<double, 9>)
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_IMPLEMENTATION_FOR_TYPE(Vector)
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_IMPLEMENTATION_FOR_TYPE(Matrix)

// MinLoc and MaxLoc AllReduce operations
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_ALLREDUCE_LOC_IMPLEMENTATION_FOR_TYPE(char)
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_ALLREDUCE_LOC_IMPLEMENTATION_FOR_TYPE(int)
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_ALLREDUCE_LOC_IMPLEMENTATION_FOR_TYPE(unsigned int)
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_ALLREDUCE_LOC_IMPLEMENTATION_FOR_TYPE(long unsigned int)
KRATOS_MPI_DATA_COMMUNICATOR_DEFINE_ALLREDUCE_LOC_IMPLEMENTATION_FOR_TYPE(double)

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

template<class TDataType> bool MPIDataCommunicator::SynchronizeShapeDetail(TDataType& rValue) const
{
    if constexpr(MPIMessage<TDataType>::HasDynamicMemoryAllocation) {
        MPIMessage<TDataType> mpi_message;
        const auto& shape = mpi_message.Shape(rValue);
        const auto& reduced_shape = MaxAll(shape);
        return mpi_message.Resize(rValue, reduced_shape);
    } else {
        return false;
    }
}

template<class TDataType> bool MPIDataCommunicator::SynchronizeShapeDetail(
    const TDataType& rSendValue,
    const int SendDestination,
    const int SendTag,
    TDataType& rRecvValue,
    const int RecvSource,
    const int RecvTag) const
{
    if constexpr(MPIMessage<TDataType>::HasDynamicMemoryAllocation) {
        // first shapes needs to be communicated.
        const std::vector<unsigned int>& send_shape = MPIMessage<TDataType>().Shape(rSendValue);

        // first communicate the number of dimensions.
        const unsigned int send_dims = send_shape.size();
        unsigned int recv_dims = 0;
        int ierr = MPI_Sendrecv(&send_dims, 1, MPI_UNSIGNED, SendDestination, SendTag, &recv_dims, 1, MPI_UNSIGNED, RecvSource, RecvTag, mComm, MPI_STATUS_IGNORE);
        CheckMPIErrorCode(ierr, "MPI_Sendrecv");

        // // now communicate the shapes
        std::vector<unsigned int> recv_shape(recv_dims);
        MPI_Sendrecv(send_shape.data(), send_shape.size(), MPI_UNSIGNED, SendDestination, SendTag, recv_shape.data(), recv_dims, MPI_UNSIGNED, RecvSource, RecvTag, mComm, MPI_STATUS_IGNORE);

        return MPIMessage<TDataType>().Resize(rRecvValue, recv_shape);
    } else {
        return false;
    }
}

template<class TDataType> void MPIDataCommunicator::ReduceDetail(
    const TDataType& rLocalValues, TDataType& rReducedValues,
    MPI_Op Operation, const int Root) const
{
    MPIMessage<TDataType> mpi_send_msg, mpi_recv_msg;

#ifdef KRATOS_DEBUG
    KRATOS_ERROR_IF_NOT(ErrorIfFalseOnAnyRank(IsValidRank(Root)))
    << "In call to MPI_Reduce: " << Root << " is not a valid rank." << std::endl;
    const int local_size = mpi_send_msg.Size(rLocalValues);
    const int reduced_size = mpi_recv_msg.Size(rReducedValues);
    KRATOS_ERROR_IF_NOT(IsEqualOnAllRanks(local_size))
    << "Input error in call to MPI_Reduce: "
    << "There should be the same amount of local values to send from each rank." << std::endl;
    KRATOS_ERROR_IF(BroadcastErrorIfTrue(local_size != reduced_size,Root))
    << "Input error in call to MPI_Reduce for rank " << Root << ": "
    << "Sending " << local_size << " values " << "but receiving " << reduced_size << " values." << std::endl;
#endif // KRATOS_DEBUG

    const int ierr = MPI_Reduce(
        mpi_send_msg.Buffer(rLocalValues), mpi_recv_msg.Buffer(rReducedValues),
        mpi_send_msg.Size(rLocalValues), mpi_send_msg.DataType(),
        Operation, Root, mComm);

    CheckMPIErrorCode(ierr, "MPI_Reduce");

    if (Rank() == Root) {
        mpi_recv_msg.Update(rReducedValues);
    }
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

    TDataType temp{};
    if (rLocalValues.size() > 0) {
        temp = rLocalValues.front();
    }

    SynchronizeShape(temp);

    if (Rank() == Root) {
        reduced_values.resize(rLocalValues.size(), temp);
    }

    ReduceDetail(rLocalValues, reduced_values, Operation, Root);
    return reduced_values;
}


template<class TDataType> void MPIDataCommunicator::AllReduceDetail(
    const TDataType& rLocalValues, TDataType& rReducedValues,
    MPI_Op Operation) const
{
    MPIMessage<TDataType> mpi_send_msg, mpi_recv_msg;

#ifdef KRATOS_DEBUG
    const int local_size = mpi_send_msg.Size(rLocalValues);
    const int reduced_size = mpi_recv_msg.Size(rReducedValues);
    KRATOS_ERROR_IF_NOT(IsEqualOnAllRanks(local_size))
    << "Input error in call to MPI_Allreduce: "
    << "There should be the same amount of local values to send from each rank." << std::endl;
    KRATOS_ERROR_IF(ErrorIfTrueOnAnyRank(local_size != reduced_size))
    << "Input error in call to MPI_Allreduce for rank " << Rank() << ": "
    << "Sending " << local_size << " values " << "but receiving " << reduced_size << " values." << std::endl;
#endif // KRATOS_DEBUG

    const int ierr = MPI_Allreduce(
        mpi_send_msg.Buffer(rLocalValues), mpi_recv_msg.Buffer(rReducedValues),
        mpi_send_msg.Size(rLocalValues), mpi_send_msg.DataType(),
        Operation, mComm);
    CheckMPIErrorCode(ierr, "MPI_Allreduce");

    mpi_recv_msg.Update(rReducedValues);
}

template<class TDataType>
TDataType MPIDataCommunicator::AllReduceDetail(
    const TDataType& rLocalValues,
    MPI_Op Operation
    ) const
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
    TDataType temp{};
    if (rLocalValues.size() > 0) {
        temp = rLocalValues.front();
    }

    SynchronizeShape(temp);

    std::vector<TDataType> reduced_values(rLocalValues.size(), temp);
    AllReduceDetail(rLocalValues, reduced_values, Operation);
    return reduced_values;
}

template<class TDataType>
std::pair<TDataType, int> MPIDataCommunicator::AllReduceDetailWithLocation(
    const std::pair<TDataType, int>& rLocalValues,
    MPI_Op Operation
    ) const
{
    struct {
        std::conditional_t<std::is_same<TDataType,char>::value, int, TDataType> value;
        int rank;
    } local_reduce, global_reduce;
    local_reduce.value = rLocalValues.first;
    local_reduce.rank = rLocalValues.second;
    MPI_Datatype data_type = MPIMessage<std::pair<TDataType, int>>().DataType();
    MPI_Allreduce(&local_reduce, &global_reduce, 1, data_type, Operation, mComm);
    std::pair<TDataType, int> global_values({static_cast<TDataType>(global_reduce.value), global_reduce.rank});
    return global_values;
}

template<class TDataType> void MPIDataCommunicator::ScanDetail(
    const TDataType& rLocalValues, TDataType& rReducedValues,
    MPI_Op Operation) const
{
    MPIMessage<TDataType> mpi_send_msg, mpi_recv_msg;

#ifdef KRATOS_DEBUG
    const int local_size = mpi_send_msg.Size(rLocalValues);
    const int reduced_size = mpi_recv_msg.Size(rReducedValues);
    KRATOS_ERROR_IF_NOT(IsEqualOnAllRanks(local_size))
    << "Input error in call to MPI_Scan: "
    << "There should be the same amount of local values to send from each rank." << std::endl;
    KRATOS_ERROR_IF(ErrorIfTrueOnAnyRank(local_size != reduced_size))
    << "Input error in call to MPI_Scan for rank " << Rank() << ": "
    << "Sending " << local_size << " values " << "but receiving " << reduced_size << " values." << std::endl;
#endif // KRATOS_DEBUG

    const int ierr = MPI_Scan(
        mpi_send_msg.Buffer(rLocalValues), mpi_recv_msg.Buffer(rReducedValues),
        mpi_send_msg.Size(rLocalValues), mpi_send_msg.DataType(),
        Operation, mComm);
    CheckMPIErrorCode(ierr, "MPI_Scan");

    mpi_recv_msg.Update(rReducedValues);
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
    TDataType temp{};
    if (rLocalValues.size() > 0) {
        temp = rLocalValues.front();
    }

    SynchronizeShape(temp);

    std::vector<TDataType> global_values(rLocalValues.size(), temp);
    ScanDetail(rLocalValues, global_values, MPI_SUM);
    return global_values;
}

template<class TDataType> void MPIDataCommunicator::SendRecvDetail(
    const TDataType& rSendMessage, const int SendDestination, const int SendTag,
    TDataType& rRecvMessage, const int RecvSource, const int RecvTag) const
{
    MPIMessage<TDataType> mpi_send_msg, mpi_recv_msg;

    const int ierr = MPI_Sendrecv(
        mpi_send_msg.Buffer(rSendMessage), mpi_send_msg.Size(rSendMessage), mpi_send_msg.DataType(), SendDestination, SendTag,
        mpi_recv_msg.Buffer(rRecvMessage), mpi_recv_msg.Size(rRecvMessage), mpi_recv_msg.DataType(), RecvSource, RecvTag,
        mComm, MPI_STATUS_IGNORE);
    CheckMPIErrorCode(ierr, "MPI_Sendrecv");

    mpi_recv_msg.Update(rRecvMessage);
}

template<class TDataType> TDataType MPIDataCommunicator::SendRecvDetail(
    const TDataType& rSendMessage,
    const int SendDestination, const int SendTag,
    const int RecvSource, const int RecvTag) const
{
    TDataType recv_message;

    if constexpr(MPIMessage<TDataType>::HasDynamicMemoryAllocation) {
        SynchronizeShape(rSendMessage, SendDestination, SendTag, recv_message, RecvSource,RecvTag);
    }

    SendRecvDetail(rSendMessage,SendDestination, SendTag, recv_message, RecvSource, RecvTag);
    return recv_message;
}

template<class TDataType> std::vector<TDataType> MPIDataCommunicator::SendRecvDetail(
    const std::vector<TDataType>& rSendMessage,
    const int SendDestination, const int SendTag,
    const int RecvSource, const int RecvTag) const
{
    int send_size = rSendMessage.size();
    int recv_size;
    SendRecvDetail(send_size, SendDestination, SendTag, recv_size, RecvSource, RecvTag);

    TDataType recv_temp{};

    if constexpr(MPIMessage<TDataType>::HasDynamicMemoryAllocation) {
        TDataType send_temp{};
        if (rSendMessage.size() > 0) {
            send_temp = rSendMessage.front();
        }
        SynchronizeShape(send_temp, SendDestination, SendTag, recv_temp, RecvSource,RecvTag);
    }

    std::vector<TDataType> recv_values(recv_size, recv_temp);
    SendRecvDetail(rSendMessage,SendDestination, SendTag, recv_values, RecvSource, RecvTag);
    return recv_values;
}

template<class TDataType> void MPIDataCommunicator::SendDetail(
    const TDataType& rSendValues, const int SendDestination, const int SendTag) const
{
    using sub_data_type = typename MPIMessage<TDataType>::SubDataType;

    MPIMessage<TDataType> mpi_send_message;

    if constexpr(MPIMessage<sub_data_type>::HasDynamicMemoryAllocation) {
        // first send the shape
        std::vector<unsigned int> shape;
        if (rSendValues.size() > 0) {
            shape = MPIMessage<sub_data_type>().Shape(rSendValues.front());
        } else {
            shape = MPIMessage<sub_data_type>().Shape(sub_data_type{});
        }

        // increases the SendTag by one to ensure the communicated data
        // is received by the shape recieve.
        const int ierr = MPI_Send(shape.data(), shape.size(), MPI_UNSIGNED, SendDestination, SendTag + 1, mComm);
        CheckMPIErrorCode(ierr, "MPI_Send");
    }

    const int ierr = MPI_Send(mpi_send_message.Buffer(rSendValues), mpi_send_message.Size(rSendValues),
        mpi_send_message.DataType(), SendDestination, SendTag, mComm);
    CheckMPIErrorCode(ierr, "MPI_Send");
}

template<class TDataType> void MPIDataCommunicator::RecvDetail(
    TDataType& rRecvValues, const int RecvSource, const int RecvTag) const
{
    using sub_data_type = typename MPIMessage<TDataType>::SubDataType;

    MPIMessage<TDataType> mpi_recv_message;
    MPI_Status status;

    sub_data_type temp{};

    if constexpr(MPIMessage<sub_data_type>::HasDynamicMemoryAllocation) {
        // first receive the shape dims of sub data type
        // recieves with RecvTag + 1, becaues shape information is sent
        // with this mdofied tag.
        int ierr = MPI_Probe(RecvSource, RecvTag + 1, mComm,&status);
        CheckMPIErrorCode(ierr, "MPI_Probe");
        int recv_dims;
        ierr = MPI_Get_count(&status, MPI_UNSIGNED, &recv_dims);
        CheckMPIErrorCode(ierr, "MPI_Get_count");

        // now recieve the shape of the sub data type
        std::vector<unsigned int> recv_shape(recv_dims);
        ierr = MPI_Recv(recv_shape.data(), recv_dims, MPI_UNSIGNED, RecvSource, RecvTag + 1, mComm, MPI_STATUS_IGNORE);
        CheckMPIErrorCode(ierr, "MPI_Recv");
        MPIMessage<sub_data_type>().Resize(temp, recv_shape);
    }

    int ierr = MPI_Probe(RecvSource, RecvTag, mComm, &status);
    CheckMPIErrorCode(ierr, "MPI_Probe");

    int recv_size;
    ierr = MPI_Get_count(&status, mpi_recv_message.DataType(), &recv_size);
    CheckMPIErrorCode(ierr, "MPI_Get_count");

    const unsigned int sub_data_type_size = MPIMessage<sub_data_type>().Size(temp);
    recv_size /= (sub_data_type_size > 0 ? sub_data_type_size : 1);
    if (rRecvValues.size() != (unsigned int)recv_size) {
        rRecvValues.resize(recv_size, temp);
    } else {
        for (auto& r_sub_item : rRecvValues) {
            MPIMessage<sub_data_type>().Resize(r_sub_item, MPIMessage<sub_data_type>().Shape(temp));
        }
    }

    ierr = MPI_Recv(mpi_recv_message.Buffer(rRecvValues), mpi_recv_message.Size(rRecvValues), mpi_recv_message.DataType(),
        RecvSource, RecvTag, mComm, MPI_STATUS_IGNORE);
    CheckMPIErrorCode(ierr, "MPI_Recv");

    mpi_recv_message.Update(rRecvValues);
}

template<class TDataType> void MPIDataCommunicator::BroadcastDetail(
    TDataType& rBuffer, const int SourceRank) const
{
    MPIMessage<TDataType> mpi_message;
#ifdef KRATOS_DEBUG
    KRATOS_ERROR_IF_NOT(ErrorIfFalseOnAnyRank(IsValidRank(SourceRank)))
    << "In call to MPI_Bcast: " << SourceRank << " is not a valid rank." << std::endl;
    KRATOS_ERROR_IF_NOT(IsEqualOnAllRanks(mpi_message.Size(rBuffer)))
    << "Input error in call to MPI_Bcast: "
    << "The buffer does not have the same size on all ranks." << std::endl;
#endif

    const int ierr = MPI_Bcast(
        mpi_message.Buffer(rBuffer), mpi_message.Size(rBuffer),
        mpi_message.DataType(), SourceRank, mComm);
    CheckMPIErrorCode(ierr, "MPI_Bcast");

    if (Rank() != SourceRank) {
        mpi_message.Update(rBuffer);
    }
}

template<class TSendDataType, class TRecvDataType> void MPIDataCommunicator::ScatterDetail(
    const TSendDataType& rSendValues, TRecvDataType& rRecvValues, const int SourceRank) const
{
    MPIMessage<TSendDataType> mpi_send_msg;
    MPIMessage<TRecvDataType> mpi_recv_msg;

#ifdef KRATOS_DEBUG
    KRATOS_ERROR_IF_NOT(ErrorIfFalseOnAnyRank(IsValidRank(SourceRank)))
    << "In call to MPI_Scatter: " << SourceRank << " is not a valid rank." << std::endl;
    const int send_size = mpi_send_msg.Size(rSendValues);
    const int recv_size = mpi_recv_msg.Size(rRecvValues);
    KRATOS_ERROR_IF_NOT(IsEqualOnAllRanks(recv_size))
    << "Input error in call to MPI_Scatter: "
    << "The destination buffer does not have the same size on all ranks." << std::endl;
    KRATOS_ERROR_IF(BroadcastErrorIfTrue(send_size != recv_size*Size(),SourceRank))
    << "Input error in call to MPI_Scatter for rank " << SourceRank << ": "
    << "Sending " << send_size << " values " << "but receiving " << recv_size << " values ("
    << recv_size * Size() << " values to send expected)." << std::endl;
#endif // KRATOS_DEBUG

    const int sends_per_rank = mpi_recv_msg.Size(rRecvValues);
    const int ierr = MPI_Scatter(
        mpi_send_msg.Buffer(rSendValues), sends_per_rank, mpi_send_msg.DataType(),
        mpi_recv_msg.Buffer(rRecvValues), sends_per_rank, mpi_recv_msg.DataType(),
        SourceRank, mComm);
    CheckMPIErrorCode(ierr, "MPI_Scatter");

    mpi_recv_msg.Update(rRecvValues);
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

    std::vector<TDataType> message;
    if (message_size > 0) {
        TDataType temp{};
        if (Rank() == SourceRank) {
            temp = rSendValues.front();
        }
        SynchronizeShape(temp);
        message.resize(message_size, temp);
        ScatterDetail(rSendValues, message, SourceRank);
    }

    return message;
}

template<class TDataType> void MPIDataCommunicator::ScattervDetail(
        const TDataType& rSendValues, const std::vector<int>& rSendCounts, const std::vector<int>& rSendOffsets,
        TDataType& rRecvValues, const int SourceRank) const
{
#ifdef KRATOS_DEBUG
    ValidateScattervInput(rSendValues, rSendCounts, rSendOffsets, rRecvValues, SourceRank);
#endif

    MPIMessage<TDataType> mpi_send_msg, mpi_recv_msg;

    if constexpr(MPIMessage<TDataType>::HasContiguousPrimitiveData) {
        const int ierr = MPI_Scatterv(
            mpi_send_msg.Buffer(rSendValues), rSendCounts.data(), rSendOffsets.data(), mpi_send_msg.DataType(),
            mpi_recv_msg.Buffer(rRecvValues), mpi_recv_msg.Size(rRecvValues), mpi_recv_msg.DataType(),
            SourceRank, mComm);
        CheckMPIErrorCode(ierr, "MPI_Scatterv");
    } else {
        // now we have update the rSendCounts and rSendOffsets properly for primitive data type sizes
        // because the TDataType is not contiguous
        int sub_data_type_length = mpi_send_msg.SubDataTypeSize(rSendValues);

        std::vector<int> primitive_send_counts(rSendCounts.size());
        std::vector<int> primitive_send_offsets(rSendOffsets.size());

        std::transform(rSendCounts.begin(), rSendCounts.end(), primitive_send_counts.begin(), [sub_data_type_length](const auto SendCount) { return SendCount * sub_data_type_length; });
        std::transform(rSendOffsets.begin(), rSendOffsets.end(), primitive_send_offsets.begin(), [sub_data_type_length](const auto SendOffset) { return SendOffset * sub_data_type_length; });

        const int ierr = MPI_Scatterv(
            mpi_send_msg.Buffer(rSendValues), primitive_send_counts.data(), primitive_send_offsets.data(), mpi_send_msg.DataType(),
            mpi_recv_msg.Buffer(rRecvValues), mpi_recv_msg.Size(rRecvValues), mpi_recv_msg.DataType(),
            SourceRank, mComm);
        CheckMPIErrorCode(ierr, "MPI_Scatterv");
    }

    mpi_recv_msg.Update(rRecvValues);
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
    MPIMessage<TSendDataType> mpi_send_msg;
    MPIMessage<TRecvDataType> mpi_recv_msg;

#ifdef KRATOS_DEBUG
    KRATOS_ERROR_IF_NOT(ErrorIfFalseOnAnyRank(IsValidRank(RecvRank)))
    << "In call to MPI_Gather: " << RecvRank << " is not a valid rank." << std::endl;
    const int send_size = mpi_send_msg.Size(rSendValues);
    const int recv_size = mpi_recv_msg.Size(rRecvValues);
    KRATOS_ERROR_IF_NOT(IsEqualOnAllRanks(send_size))
    << "Input error in call to MPI_Gather: "
    << "There should be the same amount of local values to send from each rank." << std::endl;
    KRATOS_ERROR_IF(BroadcastErrorIfTrue(send_size*Size() != recv_size,RecvRank))
    << "Input error in call to MPI_Gather for rank " << RecvRank << ": "
    << "Sending " << send_size << " values " << "but receiving " << recv_size << " values ("
    << send_size * Size() << " values to receive expected)." << std::endl;
#endif // KRATOS_DEBUG

    const int sends_per_rank = mpi_send_msg.Size(rSendValues);
    const int ierr = MPI_Gather(
        mpi_send_msg.Buffer(rSendValues), sends_per_rank, mpi_send_msg.DataType(),
        mpi_recv_msg.Buffer(rRecvValues), sends_per_rank, mpi_recv_msg.DataType(),
        RecvRank, mComm);
    CheckMPIErrorCode(ierr, "MPI_Gather");

    if (Rank() == RecvRank) {
        mpi_recv_msg.Update(rRecvValues);
    }
}

template<class TDataType> std::vector<TDataType> MPIDataCommunicator::GatherDetail(
    const std::vector<TDataType>& rSendValues, const int DestinationRank) const
{
    int message_size = rSendValues.size();

    TDataType temp{};
    if (rSendValues.size() > 0) {
        temp = rSendValues.front();
    }
    SynchronizeShape(temp);

    std::vector<TDataType> gathered_values;
    if (Rank() == DestinationRank) {
        gathered_values.resize(message_size*Size(), temp);
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

    MPIMessage<TDataType> mpi_send_msg, mpi_recv_msg;

    if constexpr(MPIMessage<TDataType>::HasContiguousPrimitiveData) {
        const int ierr = MPI_Gatherv(
            mpi_send_msg.Buffer(rSendValues), mpi_send_msg.Size(rSendValues), mpi_send_msg.DataType(),
            mpi_recv_msg.Buffer(rRecvValues), rRecvCounts.data(), rRecvOffsets.data(), mpi_recv_msg.DataType(),
            RecvRank, mComm);
        CheckMPIErrorCode(ierr, "MPI_Gatherv");
    } else {
        // now we have update the rRecvCounts and rRecvOffsets properly for primitive data type sizes
        // because the TDataType is not contiguous

        int sub_data_type_length = mpi_recv_msg.SubDataTypeSize(rRecvValues);

        std::vector<int> primitive_recv_counts(rRecvCounts.size());
        std::vector<int> primitive_recv_offsets(rRecvOffsets.size());

        std::transform(rRecvCounts.begin(), rRecvCounts.end(), primitive_recv_counts.begin(), [sub_data_type_length](const auto RecvCount) { return RecvCount * sub_data_type_length; });
        std::transform(rRecvOffsets.begin(), rRecvOffsets.end(), primitive_recv_offsets.begin(), [sub_data_type_length](const auto RecvOffset) { return RecvOffset * sub_data_type_length; });

        const int ierr = MPI_Gatherv(
            mpi_send_msg.Buffer(rSendValues), mpi_send_msg.Size(rSendValues), mpi_send_msg.DataType(),
            mpi_recv_msg.Buffer(rRecvValues), primitive_recv_counts.data(), primitive_recv_offsets.data(), mpi_recv_msg.DataType(),
            RecvRank, mComm);
        CheckMPIErrorCode(ierr, "MPI_Scatterv");
    }

    if (Rank() == RecvRank) {
        mpi_recv_msg.Update(rRecvValues);
    }
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
    MPIMessage<TDataType> mpi_send_msg, mpi_recv_msg;

#ifdef KRATOS_DEBUG
    const int send_size = mpi_send_msg.Size(rSendValues);
    const int recv_size = mpi_recv_msg.Size(rRecvValues);
    KRATOS_ERROR_IF_NOT(IsEqualOnAllRanks(send_size))
    << "Input error in call to MPI_Allgather: "
    << "There should be the same amount of local values to send from each rank." << std::endl;
    KRATOS_ERROR_IF(ErrorIfTrueOnAnyRank(send_size*Size() != recv_size))
    << "Input error in call to MPI_Allgather for rank " << Rank() << ": "
    << "Sending " << send_size << " values " << "but receiving " << recv_size << " values ("
    << send_size * Size() << " values to receive expected)." << std::endl;
#endif // KRATOS_DEBUG

    const int sends_per_rank = mpi_send_msg.Size(rSendValues);
    const int ierr = MPI_Allgather(
        mpi_send_msg.Buffer(rSendValues), sends_per_rank, mpi_send_msg.DataType(),
        mpi_recv_msg.Buffer(rRecvValues), sends_per_rank, mpi_recv_msg.DataType(),
        mComm);
    CheckMPIErrorCode(ierr, "MPI_Allgather");

    mpi_recv_msg.Update(rRecvValues);
}

template<class TDataType> std::vector<TDataType> MPIDataCommunicator::AllGatherDetail(
    const std::vector<TDataType>& rSendValues) const
{
    TDataType temp{};
    if (rSendValues.size() > 0) {
        temp = rSendValues.front();
    }
    SynchronizeShape(temp);
    std::vector<TDataType> output(rSendValues.size()*Size(), temp);
    AllGatherDetail(rSendValues, output);
    return output;
}

template<class TDataType>
void MPIDataCommunicator::AllGathervDetail(
    const TDataType& rSendValues, TDataType& rRecvValues,
    const std::vector<int>& rRecvCounts, const std::vector<int>& rRecvOffsets) const
{
#ifdef KRATOS_DEBUG
    ValidateAllGathervInput(rSendValues, rRecvValues, rRecvCounts, rRecvOffsets);
#endif // KRATOS_DEBUG

    MPIMessage<TDataType> mpi_send_msg, mpi_recv_msg;

    if constexpr(MPIMessage<TDataType>::HasContiguousPrimitiveData) {
        const int ierr = MPI_Allgatherv(
            mpi_send_msg.Buffer(rSendValues), mpi_send_msg.Size(rSendValues), mpi_send_msg.DataType(),
            mpi_recv_msg.Buffer(rRecvValues), rRecvCounts.data(), rRecvOffsets.data(), mpi_recv_msg.DataType(),
            mComm);
        CheckMPIErrorCode(ierr, "MPI_Allgatherv");
    } else {
        // now we have update the rRecvCounts and rRecvOffsets properly for primitive data type sizes
        // because the TDataType is not contiguous

        int sub_data_type_length = mpi_recv_msg.SubDataTypeSize(rRecvValues);

        std::vector<int> primitive_recv_counts(rRecvCounts.size());
        std::vector<int> primitive_recv_offsets(rRecvOffsets.size());

        std::transform(rRecvCounts.begin(), rRecvCounts.end(), primitive_recv_counts.begin(), [sub_data_type_length](const auto RecvCount) { return RecvCount * sub_data_type_length; });
        std::transform(rRecvOffsets.begin(), rRecvOffsets.end(), primitive_recv_offsets.begin(), [sub_data_type_length](const auto RecvOffset) { return RecvOffset * sub_data_type_length; });

        const int ierr = MPI_Allgatherv(
            mpi_send_msg.Buffer(rSendValues), mpi_send_msg.Size(rSendValues), mpi_send_msg.DataType(),
            mpi_recv_msg.Buffer(rRecvValues), primitive_recv_counts.data(), primitive_recv_offsets.data(), mpi_recv_msg.DataType(),
            mComm);
        CheckMPIErrorCode(ierr, "MPI_Allgatherv");
    }

    mpi_recv_msg.Update(rRecvValues);
}

template<class TDataType>
std::vector<std::vector<TDataType>> MPIDataCommunicator::AllGathervDetail(
    const std::vector<TDataType>& rSendValues) const
{
    std::vector<TDataType> message;
    std::vector<int> message_lengths;
    std::vector<int> message_offsets;
    PrepareAllGathervBuffers(rSendValues, message, message_lengths, message_offsets);
    AllGatherv(rSendValues, message, message_lengths, message_offsets);
    std::vector<std::vector<TDataType>> output_message;
    PrepareAllGathervReturn(message, message_lengths, message_offsets, output_message);
    return output_message;
}

bool MPIDataCommunicator::BroadcastErrorIfTrue(bool Condition, const int SourceRank) const
{
    const int ierr = MPI_Bcast(&Condition,1,MPI_C_BOOL,SourceRank,mComm);
    CheckMPIErrorCode(ierr, "MPI_Bcast");
    const int rank = Rank();
    KRATOS_ERROR_IF(Condition && (rank != SourceRank) )
    << "Rank " << rank << ": Stopping because of error in rank " << SourceRank << "." << std::endl;
    return Condition;
}

bool MPIDataCommunicator::BroadcastErrorIfFalse(bool Condition, const int SourceRank) const
{
    const int ierr = MPI_Bcast(&Condition,1,MPI_C_BOOL,SourceRank,mComm);
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
    const int ierr = MPI_Allreduce(&Condition, &or_condition, 1, MPI_C_BOOL, MPI_LOR, mComm);
    CheckMPIErrorCode(ierr, "MPI_Allreduce");
    KRATOS_ERROR_IF(or_condition && !Condition)
    << "Rank " << Rank() << ": Stopping because an error was detected on a different rank." << std::endl;
    return or_condition;
}

bool MPIDataCommunicator::ErrorIfFalseOnAnyRank(bool Condition) const
{
    bool and_condition;
    const int ierr = MPI_Allreduce(&Condition, &and_condition, 1, MPI_C_BOOL, MPI_LAND, mComm);
    CheckMPIErrorCode(ierr, "MPI_Allreduce");
    KRATOS_ERROR_IF(!and_condition && Condition)
    << "Rank " << Rank() << ": Stopping because an error was detected on a different rank." << std::endl;
    return and_condition;
}

bool MPIDataCommunicator::IsEqualOnAllRanks(const int LocalValue) const
{
    int local_buffer[2]{LocalValue, -LocalValue};
    int min_buffer[2]{0, 0};
    const int ierr = MPI_Allreduce(&local_buffer,&min_buffer,2,MPI_INT,MPI_MIN,mComm);
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
    const int available_recv_size = rRecvValues.size();
    int ierr = MPI_Scatter(rSendCounts.data(), 1, MPI_INT, &expected_size, 1, MPI_INT, SourceRank, mComm);
    CheckMPIErrorCode(ierr, "MPI_Scatter");
    KRATOS_ERROR_IF(ErrorIfTrueOnAnyRank(expected_size != available_recv_size))
    << "Input error in call to MPI_Scatterv for rank " << Rank() << ": "
    << "This rank will receive " << expected_size << " values but the receive buffer has size "
    << available_recv_size << "." << std::endl;

    // Message size is not smaller than total expected size (can only check for too small, since the source message may be padded).
    int total_size = 0;
    const int message_size = rSendValues.size();
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
    const int send_size = rSendValues.size();
    int ierr = MPI_Scatter(rRecvCounts.data(), 1, MPI_INT, &expected_recv_size, 1, MPI_INT, RecvRank, mComm);
    CheckMPIErrorCode(ierr, "MPI_Scatter");
    KRATOS_ERROR_IF(ErrorIfTrueOnAnyRank(send_size != expected_recv_size))
    << "Input error in call to MPI_Gatherv for rank " << Rank() << ": "
    << "This rank will send " << send_size << " values but " << RecvRank << " expects "
    << expected_recv_size << " values from it." << std::endl;

    // Message size is not larger than total expected size (can only check for too large, since the recv message may be padded).
    int total_size = 0;
    const int message_size = rSendValues.size();
    const int expected_message_size = rRecvValues.size();
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

template<class TDataType>
void MPIDataCommunicator::ValidateAllGathervInput(
    const TDataType& rSendValues, TDataType& rRecvValues,
    const std::vector<int>& rRecvCounts, const std::vector<int>& rRecvOffsets) const
{
    // All ranks send a message of the expected size
    const int send_size = rSendValues.size();
    const int comm_size = Size();
    std::vector<int> effective_recv_sizes(comm_size);
    const int ierr = MPI_Allgather(&send_size, 1, MPI_INT, effective_recv_sizes.data(), 1, MPI_INT, mComm);
    CheckMPIErrorCode(ierr, "MPI_Allgather");
    for (int i = 0; i < comm_size; i++) {
        KRATOS_ERROR_IF(ErrorIfTrueOnAnyRank(rRecvCounts[i] != effective_recv_sizes[i]))
        << "Input error in call to MPI_Allgatherv for rank " << Rank() << ": "
        << "This rank will receive " << effective_recv_sizes[i] << " values from rank " << i
        << " values but expects " << rRecvCounts[i] << " values from it." << std::endl;
    }

    // Message size is not larger than total expected size (can only check for too large, since the recv message may be padded).
    const int message_size = rSendValues.size();
    const int expected_message_size = rRecvValues.size();
    const int total_size = this->SumAll(message_size);
    KRATOS_ERROR_IF(total_size > expected_message_size)
    << "Input error in call to MPI_Allgatherv for rank " << Rank() << ": "
    << "The sent messages contain " << total_size << " values in total, but only "
    << expected_message_size << " values are expected in all ranks." << std::endl;

    // No overflow in recv buffer.
    std::stringstream message;
    bool failed = false;
    for (int i = 0; i < Size(); i++) {
        if (rRecvOffsets[i]+rRecvCounts[i] > expected_message_size) {
            message
            << "Input error in call to MPI_Allgatherv for rank " << Rank() << ": "
            << "Writing past buffer end when sending message for rank " << i << "." << std::endl;
            failed = true;
            break;
        }
    }
    KRATOS_ERROR_IF(failed) << message.str();
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

    TDataType temp{};
    if (rScattervMessage.size() > 0) {
        temp = rScattervMessage.front();
    }
    SynchronizeShape(temp);

    int result_size;
    ScatterDetail(rMessageLengths, result_size, SourceRank);
    rResult.resize(result_size, temp);
}


template<class TDataType> void MPIDataCommunicator::PrepareGathervBuffers(
    const std::vector<TDataType>& rGathervInput,
    std::vector<TDataType>& rGathervMessage,
    std::vector<int>& rMessageLengths,
    std::vector<int>& rMessageDistances,
    const int DestinationRank) const
{
    const int message_size_send = rGathervInput.size();
    const int rank = Rank();
    const int size = Size();
    if (rank == DestinationRank)
    {
        rMessageLengths.resize(size);
    }
    GatherDetail(message_size_send, rMessageLengths, DestinationRank);

    TDataType temp{};
    if (rGathervInput.size() > 0) {
        temp = rGathervInput.front();
    }

    SynchronizeShape(temp);

    if (rank == DestinationRank)
    {
        rMessageDistances.resize(size);
        int message_size = 0;
        for (int i = 0; i < size; i++)
        {
            rMessageDistances[i] = message_size;
            message_size += rMessageLengths[i];
        }
        rGathervMessage.resize(message_size, temp);
    }
}

template<class TDataType>
void MPIDataCommunicator::PrepareAllGathervBuffers(
    const std::vector<TDataType>& rGathervInput,
    std::vector<TDataType>& rGathervMessage,
    std::vector<int>& rMessageLengths,
    std::vector<int>& rMessageDistances) const
{
    const int size = Size();
    std::vector<int> message_size_send(1, rGathervInput.size());
    rMessageLengths.resize(size);
    AllGatherDetail(message_size_send, rMessageLengths);
    rMessageDistances.resize(size);
    int message_size = 0;
    for (int i = 0; i < size; i++) {
        rMessageDistances[i] = message_size;
        message_size += rMessageLengths[i];
    }

    TDataType temp{};
    if (rGathervInput.size() > 0) {
        temp = rGathervInput.front();
    }

    SynchronizeShape(temp);
    rGathervMessage.resize(message_size, temp);
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

template<class TDataType>
void MPIDataCommunicator::PrepareAllGathervReturn(
    const std::vector<TDataType>& rGathervMessage,
    const std::vector<int>& rMessageLengths,
    const std::vector<int>& rMessageDistances,
    std::vector<std::vector<TDataType>>& rOutputMessage) const
{
    const int size = Size();
    rOutputMessage.resize(size);
    for (int i = 0, counter = 0; i < size; i++) {
        rOutputMessage[i].resize(rMessageLengths[i]);
        for (int j = 0; j < rMessageLengths[i]; j++, counter++) {
            rOutputMessage[i][j] = rGathervMessage[counter];
        }
    }
}


// MPI_Datatype wrapper
template<class TValue> inline MPI_Datatype MPIDataCommunicator::MPIDatatype(const TValue&) const
{
    return MPIMessage<TValue>().DataType();
}

// Buffer argument deduction
template<class TContainer> inline void* MPIDataCommunicator::MPIBuffer(TContainer& rValues) const
{
    return MPIMessage<TContainer>().Buffer(rValues);
}

template<class TContainer> inline const void* MPIDataCommunicator::MPIBuffer(const TContainer& rValues) const
{
    return MPIMessage<TContainer>().Buffer(rValues);
}

// MPI message size deduction
template<class TContainer> inline int MPIDataCommunicator::MPIMessageSize(const TContainer& rValues) const
{
    return MPIMessage<TContainer>().Size(rValues);
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
