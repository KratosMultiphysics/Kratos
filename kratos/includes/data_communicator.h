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

#ifndef KRATOS_DATA_COMMUNICATOR_H_INCLUDED
#define KRATOS_DATA_COMMUNICATOR_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "containers/array_1d.h"
#include "containers/flags.h"
#include "includes/define.h"

// Using a macro instead of a function to get the correct line in the error message.
#ifndef KRATOS_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK
#define KRATOS_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(Size1, Size2, CheckedFunction) \
    KRATOS_DEBUG_ERROR_IF(Size1 != Size2) \
    << "Input error in call to DataCommunicator::" << CheckedFunction \
    << ": The sizes of the local and distributed buffers do not match." << std::endl;
#endif

// Methods based on MPI_Reduce, supporting sum, max or min operations.
/* Variants for each method are provided, either returning the reduced value or filling a provided vector buffer.
 * The returned value is only meaningful on the Root rank.
 */
#ifndef KRATOS_BASE_DATA_COMMUNICATOR_DECLARE_REDUCE_INTERFACE_FOR_TYPE
#define KRATOS_BASE_DATA_COMMUNICATOR_DECLARE_REDUCE_INTERFACE_FOR_TYPE(type)                                       \
virtual type Sum(const type rLocalValue, const int Root) const { return rLocalValue; }                              \
virtual std::vector<type> Sum(const std::vector<type>& rLocalValues, const int Root) const {                        \
    return rLocalValues;                                                                                            \
}                                                                                                                   \
virtual void Sum(const std::vector<type>& rLocalValues, std::vector<type>& rGlobalValues, const int Root) const {   \
    KRATOS_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rLocalValues.size(), rGlobalValues.size(), "Sum");                    \
    rGlobalValues = Sum(rLocalValues, Root);                                                                        \
}                                                                                                                   \
virtual type Min(const type rLocalValue, const int Root) const { return rLocalValue; }                              \
virtual std::vector<type> Min(const std::vector<type>& rLocalValues, const int Root) const {                        \
    return rLocalValues;                                                                                            \
}                                                                                                                   \
virtual void Min(const std::vector<type>& rLocalValues, std::vector<type>& rGlobalValues, const int Root) const {   \
    KRATOS_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rLocalValues.size(), rGlobalValues.size(), "Min");                    \
    rGlobalValues = Min(rLocalValues, Root);                                                                        \
}                                                                                                                   \
virtual type Max(const type rLocalValue, const int Root) const { return rLocalValue; }                              \
virtual std::vector<type> Max(const std::vector<type>& rLocalValues, const int Root) const {                        \
    return rLocalValues;                                                                                            \
}                                                                                                                   \
virtual void Max(const std::vector<type>& rLocalValues, std::vector<type>& rGlobalValues, const int Root) const {   \
    KRATOS_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rLocalValues.size(), rGlobalValues.size(), "Max");                    \
    rGlobalValues = Max(rLocalValues, Root);                                                                        \
}                                                                                                                   \

#endif

// Methods based on MPI_Allreduce, supporting sum, max or min operations.
/* Variants for each method are provided, either returning the reduced value or filling a provided vector buffer.
 * The returned value is defined on all ranks.
 */
#ifndef KRATOS_BASE_DATA_COMMUNICATOR_DECLARE_ALLREDUCE_INTERFACE_FOR_TYPE
#define KRATOS_BASE_DATA_COMMUNICATOR_DECLARE_ALLREDUCE_INTERFACE_FOR_TYPE(type)                        \
virtual type SumAll(const type rLocalValue) const { return rLocalValue; }                               \
virtual std::vector<type> SumAll(const std::vector<type>& rLocalValues) const {                         \
    return rLocalValues;                                                                                \
}                                                                                                       \
virtual void SumAll(const std::vector<type>& rLocalValues, std::vector<type>& rGlobalValues) const {    \
    KRATOS_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rLocalValues.size(), rGlobalValues.size(), "SumAll");     \
    rGlobalValues = SumAll(rLocalValues);                                                               \
}                                                                                                       \
virtual type MinAll(const type rLocalValue) const { return rLocalValue; }                               \
virtual std::vector<type> MinAll(const std::vector<type>& rLocalValues) const {                         \
    return rLocalValues;                                                                                \
}                                                                                                       \
virtual void MinAll(const std::vector<type>& rLocalValues, std::vector<type>& rGlobalValues) const {    \
    KRATOS_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rLocalValues.size(), rGlobalValues.size(), "MinAll");     \
    rGlobalValues = MinAll(rLocalValues);                                                               \
}                                                                                                       \
virtual type MaxAll(const type rLocalValue) const { return rLocalValue; }                               \
virtual std::vector<type> MaxAll(const std::vector<type>& rLocalValues) const {                         \
    return rLocalValues;                                                                                \
}                                                                                                       \
virtual void MaxAll(const std::vector<type>& rLocalValues, std::vector<type>& rGlobalValues) const {    \
    KRATOS_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rLocalValues.size(), rGlobalValues.size(), "MaxAll");     \
    rGlobalValues = MaxAll(rLocalValues);                                                               \
}                                                                                                       \

#endif

// Compute the partial sum of the given quantity from rank 0 to the current rank (included).
/* This is a wrapper to MPI_Scan.
 * Variants for each method are provided, either returning the reduced value or filling a provided vector buffer.
 * The returned value is defined on all ranks.
 */
#ifndef KRATOS_BASE_DATA_COMMUNICATOR_DECLARE_SCANSUM_INTERFACE_FOR_TYPE
#define KRATOS_BASE_DATA_COMMUNICATOR_DECLARE_SCANSUM_INTERFACE_FOR_TYPE(type)                      \
virtual type ScanSum(const type rLocalValue) const { return rLocalValue; }                          \
virtual std::vector<type> ScanSum(const std::vector<type>& rLocalValues) const {                    \
    return rLocalValues;                                                                            \
}                                                                                                   \
virtual void ScanSum(const std::vector<type>& rLocalValues, std::vector<type>& rPartialSums) const {\
    KRATOS_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rLocalValues.size(), rPartialSums.size(), "ScanSum"); \
    rPartialSums = ScanSum(rLocalValues);                                                           \
}                                                                                                   \

#endif

// Exchange data with other ranks. This is a wrapper for MPI_Sendrecv.
/* Versions which outputting the result as a return argument or by filling an output buffer argument are provided.
 * The return version has a performance overhead, since the dimensions of the receiving buffer have to be
 * communicated. If the dimensions of the receiving buffer are known at the destination rank, the output buffer
 * variant should be preferred.
 */
#ifndef KRATOS_BASE_DATA_COMMUNICATOR_DECLARE_SENDRECV_INTERFACE_FOR_TYPE
#define KRATOS_BASE_DATA_COMMUNICATOR_DECLARE_SENDRECV_INTERFACE_FOR_TYPE(type)                                 \
virtual std::vector<type> SendRecv(                                                                             \
    const std::vector<type>& rSendValues, const int SendDestination, const int RecvSource) const {              \
    KRATOS_ERROR_IF( (Rank() != SendDestination) || (Rank() != RecvSource))                                     \
    << "Communication between different ranks is not possible with a serial DataCommunicator." << std::endl;    \
    return rSendValues;                                                                                         \
}                                                                                                               \
virtual std::vector<type> SendRecv(                                                                             \
    const std::vector<type>& rSendValues, const int SendDestination, const int SendTag,                         \
    const int RecvSource, const int RecvTag) const {                                                            \
    KRATOS_ERROR_IF( (Rank() != SendDestination) || (Rank() != RecvSource))                                     \
    << "Communication between different ranks is not possible with a serial DataCommunicator." << std::endl;    \
    return rSendValues;                                                                                         \
}                                                                                                               \
virtual void SendRecv(                                                                                          \
    const std::vector<type>& rSendValues, const int SendDestination, const int SendTag,                         \
    std::vector<type>& rRecvValues, const int RecvSource, const int RecvTag) const {                            \
    KRATOS_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rSendValues.size(), rRecvValues.size(), "SendRecv");              \
    rRecvValues = SendRecv(rSendValues, SendDestination, RecvSource);                                           \
}                                                                                                               \

#endif

// Synchronize a buffer to the value held by the broadcasting rank.
/* This is a wrapper for MPI_Bcast.
 *  @param[in/out] The broadcast value (input on SourceRank, output on all other ranks).
 *  @param[in] SourceRank The rank transmitting the value.
 */
#ifndef KRATOS_BASE_DATA_COMMUNICATOR_DECLARE_BROADCAST_INTERFACE_FOR_TYPE
#define KRATOS_BASE_DATA_COMMUNICATOR_DECLARE_BROADCAST_INTERFACE_FOR_TYPE(type)    \
virtual void Broadcast(type& rBuffer, const int SourceRank) const {}                \
virtual void Broadcast(std::vector<type>& rBuffer, const int SourceRank) const {}   \

#endif

/// Wrappers for MPI_Scatter and MPI_Scatterv calls.
/* Versions which outputting the result as a return argument or by filling an output buffer argument are provided.
 * The return version has a performance overhead, since the dimensions of the receiving buffer have to be
 * communicated. If the dimensions of the receiving buffers are known at the destination rank, the output buffer
 * variant should be preferred.
 */
#ifndef KRATOS_BASE_DATA_COMMUNICATOR_DECLARE_SCATTER_INTERFACE_FOR_TYPE
#define KRATOS_BASE_DATA_COMMUNICATOR_DECLARE_SCATTER_INTERFACE_FOR_TYPE(type)                                                              \
virtual std::vector<type> Scatter(const std::vector<type>& rSendValues, const int SourceRank) const {                                       \
     KRATOS_ERROR_IF( Rank() != SourceRank )                                                                                                \
    << "Communication between different ranks is not possible with a serial DataCommunicator." << std::endl;                                \
    return rSendValues;                                                                                                                     \
}                                                                                                                                           \
virtual void Scatter(                                                                                                                       \
    const std::vector<type>& rSendValues, std::vector<type>& rRecvValues, const int SourceRank) const {                                     \
    KRATOS_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rSendValues.size(),rRecvValues.size(),"Scatter");                                             \
    rRecvValues = Scatter(rSendValues, SourceRank);                                                                                         \
}                                                                                                                                           \
virtual std::vector<type> Scatterv(const std::vector<std::vector<type>>& rSendValues, const int SourceRank) const {                         \
    KRATOS_ERROR_IF( Rank() != SourceRank )                                                                                                 \
    << "Communication between different ranks is not possible with a serial DataCommunicator." << std::endl;                                \
    KRATOS_ERROR_IF( static_cast<unsigned int>(Size()) != rSendValues.size() )                                                              \
    << "Unexpected number of sends in DataCommuncatior::Scatterv (serial DataCommunicator always assumes a single process)." << std::endl;  \
    return rSendValues[0];                                                                                                                  \
}                                                                                                                                           \
virtual void Scatterv(                                                                                                                      \
    const std::vector<type>& rSendValues, const std::vector<int>& rSendCounts, const std::vector<int>& rSendOffsets,                        \
    std::vector<type>& rRecvValues, const int SourceRank) const {                                                                           \
    KRATOS_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rRecvValues.size(), rSendValues.size(), "Scatterv (values check)");                           \
    KRATOS_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rSendCounts.size(), 1, "Scatterv (counts check)");                                            \
    KRATOS_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rSendOffsets.size(), 1, "Scatterv (offsets check)");                                          \
    KRATOS_ERROR_IF( Rank() != SourceRank )                                                                                                 \
    << "Communication between different ranks is not possible with a serial DataCommunicator." << std::endl;                                \
    rRecvValues = rSendValues;                                                                                                              \
}                                                                                                                                           \

#endif

/// Wrappers for MPI_Gather, MPI_Gatherv and MPI_Allgather calls.
/* Versions which outputting the result as a return argument or by filling an output buffer argument are provided.
 * The return version has a performance overhead, since the dimensions of the receiving buffer have to be
 * communicated. If the dimensions of the receiving buffers are known at the destination rank, the output buffer
 * variant should be preferred.
 */
#ifndef KRATOS_BASE_DATA_COMMUNICATOR_DECLARE_GATHER_INTERFACE_FOR_TYPE
#define KRATOS_BASE_DATA_COMMUNICATOR_DECLARE_GATHER_INTERFACE_FOR_TYPE(type)                                       \
virtual std::vector<type> Gather(const std::vector<type>& rSendValues, const int DestinationRank) const {           \
    KRATOS_ERROR_IF( Rank() != DestinationRank )                                                                    \
    << "Communication between different ranks is not possible with a serial DataCommunicator." << std::endl;        \
    return rSendValues;                                                                                             \
}                                                                                                                   \
virtual void Gather(                                                                                                \
    const std::vector<type>& rSendValues, std::vector<type>& rRecvValues, const int DestinationRank) const {        \
    KRATOS_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rSendValues.size(),rRecvValues.size(),"Gather");                      \
    rRecvValues = Gather(rSendValues, DestinationRank);                                                             \
}                                                                                                                   \
virtual std::vector<std::vector<type>> Gatherv(                                                                     \
    const std::vector<type>& rSendValues, const int DestinationRank) const {                                        \
    KRATOS_ERROR_IF( Rank() != DestinationRank )                                                                    \
    << "Communication between different ranks is not possible with a serial DataCommunicator." << std::endl;        \
    return std::vector<std::vector<type>>{rSendValues};                                                             \
}                                                                                                                   \
virtual void Gatherv(                                                                                               \
    const std::vector<type>& rSendValues, std::vector<type>& rRecvValues,                                           \
    const std::vector<int>& rRecvCounts, const std::vector<int>& rRecvOffsets, const int DestinationRank) const {   \
    KRATOS_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rRecvValues.size(), rSendValues.size(), "Gatherv (values check)");    \
    KRATOS_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rRecvCounts.size(), 1, "Gatherv (counts check)");                     \
    KRATOS_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rRecvOffsets.size(), 1, "Gatherv (offset check)");                    \
    KRATOS_ERROR_IF( Rank() != DestinationRank )                                                                    \
    << "Communication between different ranks is not possible with a serial DataCommunicator." << std::endl;        \
    rRecvValues = rSendValues;                                                                                      \
}                                                                                                                   \
virtual std::vector<type> AllGather(const std::vector<type>& rSendValues) const { return rSendValues; }             \
virtual void AllGather(const std::vector<type>& rSendValues, std::vector<type>& rRecvValues) const {                \
    KRATOS_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rSendValues.size(),rRecvValues.size(),"AllGather");                   \
    rRecvValues = AllGather(rSendValues);                                                                           \
}                                                                                                                   \

#endif

#ifndef KRATOS_BASE_DATA_COMMUNICATOR_DECLARE_INTERFACE_FOR_TYPE
#define KRATOS_BASE_DATA_COMMUNICATOR_DECLARE_INTERFACE_FOR_TYPE(type)   \
KRATOS_BASE_DATA_COMMUNICATOR_DECLARE_REDUCE_INTERFACE_FOR_TYPE(type)    \
KRATOS_BASE_DATA_COMMUNICATOR_DECLARE_ALLREDUCE_INTERFACE_FOR_TYPE(type) \
KRATOS_BASE_DATA_COMMUNICATOR_DECLARE_SCANSUM_INTERFACE_FOR_TYPE(type)   \
KRATOS_BASE_DATA_COMMUNICATOR_DECLARE_SENDRECV_INTERFACE_FOR_TYPE(type)  \
KRATOS_BASE_DATA_COMMUNICATOR_DECLARE_BROADCAST_INTERFACE_FOR_TYPE(type) \
KRATOS_BASE_DATA_COMMUNICATOR_DECLARE_SCATTER_INTERFACE_FOR_TYPE(type)   \
KRATOS_BASE_DATA_COMMUNICATOR_DECLARE_GATHER_INTERFACE_FOR_TYPE(type)    \

#endif

namespace Kratos
{
///@addtogroup Kratos Core
///@{

///@name Kratos Classes
///@{

/// Serial (do-nothing) version of a wrapper class for MPI communication.
/** @see MPIDataCommunicator for a working distributed memory implementation.
  */
class KRATOS_API(KRATOS_CORE) DataCommunicator
{
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DataCommunicator
    KRATOS_CLASS_POINTER_DEFINITION(DataCommunicator);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DataCommunicator() {}

    /// Destructor.
    virtual ~DataCommunicator() {}

    ///@}
    ///@name Operations
    ///@{

    /// Create a new DataCommunicator as a copy of this one.
    /** This method is used in ParallelEnvironment to register DataCommunicators
     *  @see ParallelEnvironment.
     *  @return a unique pointer to the new DataCommunicator.
     */
    virtual DataCommunicator::UniquePointer Clone() const
    {
        return Kratos::make_unique<DataCommunicator>();
    }

    /// Pause program exectution until all threads reach this call.
    /** Wrapper for MPI_Barrier. */
    virtual void Barrier() const {}

    // Complete interface for basic types

    KRATOS_BASE_DATA_COMMUNICATOR_DECLARE_INTERFACE_FOR_TYPE(int)
    KRATOS_BASE_DATA_COMMUNICATOR_DECLARE_INTERFACE_FOR_TYPE(unsigned int)
    KRATOS_BASE_DATA_COMMUNICATOR_DECLARE_INTERFACE_FOR_TYPE(double)

    // Reduce operations

    /// Sum rLocalValue across all ranks in the Communicator (array_1d<double,3> version).
    /** This is a wrapper to MPI_Reduce.
     *  @param[in] rLocalValue Local contribution to the sum.
     *  @param[in] Root The rank where the result will be computed.
     *  @return The summed quantity (meaningful only in Root).
     */
    virtual array_1d<double,3> Sum(const array_1d<double,3>& rLocalValue, const int Root) const
    {
        return rLocalValue;
    }


    /// Obtain the minimum of rLocalValue across all ranks in the Communicator (array_1d<double,3> version).
    /** This is a wrapper to MPI_Reduce.
     *  @param[in] rLocalValue Local value to consider in computing the minimum.
     *  @param[in] Root The rank where the result will be computed.
     *  @return The minimum value (meaningful only in Root).
     */
    virtual array_1d<double,3> Min(const array_1d<double,3>& rLocalValue, const int Root) const
    {
        return rLocalValue;
    }

    /// Obtain the maximum of rLocalValue across all ranks in the Communicator (array_1d<double,3> version).
    /** This is a wrapper to MPI_Reduce.
     *  @param[in] rLocalValue Local value to consider in computing the maximum.
     *  @param[in] Root The rank where the result will be computed.
     *  @return The maximum value (meaningful only in Root).
     */
    virtual array_1d<double,3> Max(const array_1d<double,3>& rLocalValue, const int Root) const
    {
        return rLocalValue;
    }

    virtual Kratos::Flags AndReduce(
        const Kratos::Flags Values,
        const Kratos::Flags Mask,
        const int Root) const
    {
        return Values;
    }

    virtual Kratos::Flags OrReduce(
        const Kratos::Flags Values,
        const Kratos::Flags Mask,
        const int Root) const
    {
        return Values;
    }

    // Allreduce operations

    /// Sum rLocalValue across all ranks in the Communicator (array_1d<double,3> version).
    /** This is a wrapper to MPI_Alleduce.
     *  @param[in] rLocalValue Local contribution to the sum.
     *  @return The summed quantity.
     */
    virtual array_1d<double,3> SumAll(const array_1d<double,3>& rLocalValue) const
    {
        return rLocalValue;
    }

    /// Obtain the minimum of rLocalValue across all ranks in the Communicator (array_1d<double,3> version).
    /** This is a wrapper to MPI_Allreduce.
     *  @param[in] rLocalValue Local value to consider in computing the minimum.
     *  @return The minimum value.
     */
    virtual array_1d<double,3> MinAll(const array_1d<double,3>& rLocalValue) const
    {
        return rLocalValue;
    }

    /// Obtain the maximum of rLocalValue across all ranks in the Communicator (array_1d<double,3> version).
    /** This is a wrapper to MPI_Allreduce.
     *  @param[in] rLocalValue Local value to consider in computing the maximum.
     *  @return The maximum value.
     */
    virtual array_1d<double,3> MaxAll(const array_1d<double,3>& rLocalValue) const
    {
        return rLocalValue;
    }


    virtual Kratos::Flags AndReduceAll(const Kratos::Flags Values, const Kratos::Flags Mask) const
    {
        return Values;
    }

    virtual Kratos::Flags OrReduceAll(const Kratos::Flags Values, const Kratos::Flags Mask) const
    {
        return Values;
    }

    // Sendrecv operations

    /// Exchange data with other ranks (string version).
    /** This is a wrapper for MPI_Sendrecv.
     *  @param[in] rSendValues String to send to rank SendDestination.
     *  @param[in] SendDestination Rank the string will be sent to.
     *  @param[in] RecvSource Rank the string is expected from.
     *  @return Received string from rank RecvSource.
     *  @note This version has a performance penalty compared to the variant
     *  taking both input and output buffers, since the dimensions of the
     *  receiving buffer have to be communicated. If the dimensions of the
     *  receiving buffer are known at the destination rank, the other variant
     *  should be preferred.
     */
    virtual std::string SendRecv(
        const std::string& rSendValues,
        const int SendDestination,
        const int RecvSource) const
    {
        KRATOS_ERROR_IF( (Rank() != SendDestination) || (Rank() != RecvSource))
        << "Communication between different ranks is not possible with a serial DataCommunicator." << std::endl;

        return rSendValues;
    }

    /// Exchange data with other ranks (string version).
    /** This is a wrapper for MPI_Sendrecv.
     *  @param[in] rSendValues String to send to rank SendDestination.
     *  @param[in] SendDestination Rank the string will be sent to.
     *  @param[in] SendTag Message tag for sent values.
     *  @param[out] rRecvValues Received string from rank RecvSource.
     *  @param[in] RecvSource Rank the string is expected from.
     *  @param[in] RecvTag Message tag for received values.
     */
    virtual void SendRecv(
        const std::string& rSendValues, const int SendDestination, const int SendTag,
        std::string& rRecvValues, const int RecvSource, const int RecvTag) const
    {
        KRATOS_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rSendValues.size(), rRecvValues.size(), "SendRecv");
        rRecvValues = SendRecv(rSendValues, SendDestination, RecvSource);
    }

    ///@}
    ///@name Inquiry
    ///@{

    /// Retrun the parallel rank for this DataCommunicator.
    /** This is a wrapper for calls to MPI_Comm_rank. */
    virtual int Rank() const
    {
        return 0;
    }

    /// Retrun the parallel size for this DataCommunicator.
    /** This is a wrapper for calls to MPI_Comm_size. */
    virtual int Size() const
    {
        return 1;
    }

    /// Check whether this DataCommunicator is aware of parallelism.
    virtual bool IsDistributed() const
    {
        return false;
    }

    ///@}
    ///@name Access
    ///@{

    /// Convenience function to retireve the current default DataCommunicator.
    /** @return A reference to the DataCommunicator instance registered as default in ParallelEnvironment.
     */
    static DataCommunicator& GetDefault();

    ///@}
    ///@name Helper functions for error checking in MPI
    ///@{

    /// This function throws an error on ranks != Sourcerank if Condition evaluates to true.
    /** This method is intended as a helper function to force ranks to stop after an error
     *  is detected on a given rank. A typical use case would be to completely stop the simulation
     *  if an error is detected on the root process.
     *  The intended usage is something like:
     *
     *  KRATOS_ERROR_IF( data_communicator_instance.BroadcastErrorIfTrue(Condition, Root) )
     *  << "Detailed error message in Root rank";
     *
     *  If an error is detected, ranks other than Root will fail with a generic error message.
     *  Failing on the Root rank is left to the caller, so that a detailed error message can be
     *  produced.
     *
     *  @note This method should be called from all ranks, it will deadlock if called within
     *  an if(rank == some_rank) statement.
     *  @see MPIDataCommunicator.
     *  @param Condition The condition to check.
     *  @param SourceRank The rank where the condition is meaningful.
     *  @return The result of evaluating Condition.
     */
    virtual bool BroadcastErrorIfTrue(bool Condition, const int SourceRank) const
    {
        return Condition;
    }

    /// This function throws an error on ranks != Sourcerank if Condition evaluates to false.
    /** This method is intended as a helper function to force ranks to stop after an error
     *  is detected on a given rank. A typical use case would be to completely stop the simulation
     *  if an error is detected on the root process.
     *  The intended usage is something like:
     *
     *  KRATOS_ERROR_IF_NOT( data_communicator_instance.BroadcastErrorIfFalse(Condition, Root) )
     *  << "Detailed error message in Root rank";
     *
     *  If an error is detected, ranks other than Root will fail with a generic error message.
     *  Failing on the Root rank is left to the caller, so that a detailed error message can be
     *  produced.
     *
     *  @note This method should be called from all ranks, it will deadlock if called within
     *  an if(rank == some_rank) statement.
     *  @see MPIDataCommunicator.
     *  @param Condition The condition to check.
     *  @param SourceRank The rank where the condition is meaningful.
     *  @return The result of evaluating Condition.
     */
    virtual bool BroadcastErrorIfFalse(bool Condition, const int SourceRank) const
    {
        return Condition;
    }

    /// This function throws an error on ranks where Condition evaluates to false, if it evaluated to true on a different rank.
    /** This method is intended as a helper function to force ranks to stop after an error
     *  is detected on one or more ranks.
     *  The intended usage is something like:
     *
     *  KRATOS_ERROR_IF( data_communicator_instance.ErrorIfTrueOnAnyRank(Condition) )
     *  << "Detailed error message in ranks where Condition == true.";
     *
     *  If an error is detected, ranks other than those where it was detected will fail with
     *  a generic error message.
     *  Failing on the ranks where the condition is true is left to the caller,
     *  so that a detailed error message can be produced.
     *
     *  @note This method should be called from all ranks, it will deadlock if called within
     *  an if(rank == some_rank) statement.
     *  @see MPIDataCommunicator.
     *  @param Condition The condition to check.
     *  @return The result of evaluating Condition.
     */
    virtual bool ErrorIfTrueOnAnyRank(bool Condition) const
    {
        return Condition;
    }

    /// This function throws an error on ranks where Condition evaluates to true, if it evaluated to false on a different rank.
    /** This method is intended as a helper function to force ranks to stop after an error
     *  is detected on one or more ranks.
     *  The intended usage is something like:
     *
     *  KRATOS_ERROR_IF_NOT( data_communicator_instance.ErrorIfFalseOnAnyRank(Condition) )
     *  << "Detailed error message in ranks where Condition == false.";
     *
     *  If an error is detected, ranks other than those where it was detected will fail with
     *  a generic error message.
     *  Failing on the ranks where the condition is false is left to the caller,
     *  so that a detailed error message can be produced.
     *
     *  @note This method should be called from all ranks, it will deadlock if called within
     *  an if(rank == some_rank) statement.
     *  @see MPIDataCommunicator.
     *  @param Condition The condition to check.
     *  @return The result of evaluating Condition.
     */
    virtual bool ErrorIfFalseOnAnyRank(bool Condition) const
    {
        return Condition;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        PrintInfo(buffer);
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << "DataCommunicator";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const
    {
        rOStream
        << "Serial do-nothing version of the Kratos wrapper for MPI communication.\n"
        << "Rank 0 of 1 assumed." << std::endl;
    }

    ///@}

  private:

    ///@name Un accessible methods
    ///@{

    /// Copy constructor.
    DataCommunicator(DataCommunicator const &rOther) = delete;

    /// Assignment operator.
    DataCommunicator &operator=(DataCommunicator const &rOther) = delete;

    ///@}

}; // Class DataCommunicator

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                DataCommunicator &rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const DataCommunicator &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

} // namespace Kratos.

#undef KRATOS_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK

#undef KRATOS_BASE_DATA_COMMUNICATOR_DECLARE_REDUCE_INTERFACE_FOR_TYPE
#undef KRATOS_BASE_DATA_COMMUNICATOR_DECLARE_ALLREDUCE_INTERFACE_FOR_TYPE
#undef KRATOS_BASE_DATA_COMMUNICATOR_DECLARE_SCANSUM_INTERFACE_FOR_TYPE
#undef KRATOS_BASE_DATA_COMMUNICATOR_DECLARE_SENDRECV_INTERFACE_FOR_TYPE
#undef KRATOS_BASE_DATA_COMMUNICATOR_DECLARE_BROADCAST_INTERFACE_FOR_TYPE
#undef KRATOS_BASE_DATA_COMMUNICATOR_DECLARE_SCATTER_INTERFACE_FOR_TYPE
#undef KRATOS_BASE_DATA_COMMUNICATOR_DECLARE_GATHER_INTERFACE_FOR_TYPE
#undef KRATOS_BASE_DATA_COMMUNICATOR_DECLARE_INTERFACE_FOR_TYPE

#endif // KRATOS_DATA_COMMUNICATOR_H_INCLUDED  defined
