//     ______     _____ _           ________
//    / ____/___ / ___/(_)___ ___  /  _/ __ |
//   / /   / __ \\__ \/ / __ `__ \ / // / / /
//  / /___/ /_/ /__/ / / / / / / // // /_/ /
//  \____/\____/____/_/_/ /_/ /_/___/\____/
//  Kratos CoSimulationApplication
//
//  License:         BSD License, see license.txt
//
//  Main authors:    Jordi Cotela
//                   Philipp Bucher (https://github.com/philbucher)
//
// Ported from "kratos/includes/data_communicator.h"
//

#ifndef CO_SIM_IO_DATA_COMMUNICATOR_INCLUDED
#define CO_SIM_IO_DATA_COMMUNICATOR_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <type_traits>

// External includes

// Project includes
#include "define.hpp"
#include "stream_serializer.hpp"

// Using a macro instead of a function to get the correct line in the error message.
#ifndef CO_SIM_IO_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK
#define CO_SIM_IO_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(Size1, Size2, CheckedFunction) \
    CO_SIM_IO_DEBUG_ERROR_IF(Size1 != Size2) \
    << "Input error in call to DataCommunicator::" << CheckedFunction \
    << ": The sizes of the local and distributed buffers do not match." << std::endl;
#endif

// Methods based on MPI_Reduce, supporting sum, max or min operations.
/* Variants for each method are provided, either returning the reduced value or filling a provided vector buffer.
 * The returned value is only meaningful on the Root rank.
 */
#ifndef CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_REDUCE_INTERFACE_FOR_TYPE
#define CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_REDUCE_INTERFACE_FOR_TYPE(type)                                       \
virtual type Sum(const type rLocalValue, const int Root) const { return rLocalValue; }                              \
virtual std::vector<type> Sum(const std::vector<type>& rLocalValues, const int Root) const {                        \
    return rLocalValues;                                                                                            \
}                                                                                                                   \
virtual void Sum(const std::vector<type>& rLocalValues, std::vector<type>& rGlobalValues, const int Root) const {   \
    CO_SIM_IO_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rLocalValues.size(), rGlobalValues.size(), "Sum");                    \
    rGlobalValues = Sum(rLocalValues, Root);                                                                        \
}                                                                                                                   \
virtual type Min(const type rLocalValue, const int Root) const { return rLocalValue; }                              \
virtual std::vector<type> Min(const std::vector<type>& rLocalValues, const int Root) const {                        \
    return rLocalValues;                                                                                            \
}                                                                                                                   \
virtual void Min(const std::vector<type>& rLocalValues, std::vector<type>& rGlobalValues, const int Root) const {   \
    CO_SIM_IO_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rLocalValues.size(), rGlobalValues.size(), "Min");                    \
    rGlobalValues = Min(rLocalValues, Root);                                                                        \
}                                                                                                                   \
virtual type Max(const type rLocalValue, const int Root) const { return rLocalValue; }                              \
virtual std::vector<type> Max(const std::vector<type>& rLocalValues, const int Root) const {                        \
    return rLocalValues;                                                                                            \
}                                                                                                                   \
virtual void Max(const std::vector<type>& rLocalValues, std::vector<type>& rGlobalValues, const int Root) const {   \
    CO_SIM_IO_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rLocalValues.size(), rGlobalValues.size(), "Max");                    \
    rGlobalValues = Max(rLocalValues, Root);                                                                        \
}                                                                                                                   \

#endif

// Methods based on MPI_Allreduce, supporting sum, max or min operations.
/* Variants for each method are provided, either returning the reduced value or filling a provided vector buffer.
 * The returned value is defined on all ranks.
 */
#ifndef CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_ALLREDUCE_INTERFACE_FOR_TYPE
#define CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_ALLREDUCE_INTERFACE_FOR_TYPE(type)                        \
virtual type SumAll(const type rLocalValue) const { return rLocalValue; }                               \
virtual std::vector<type> SumAll(const std::vector<type>& rLocalValues) const {                         \
    return rLocalValues;                                                                                \
}                                                                                                       \
virtual void SumAll(const std::vector<type>& rLocalValues, std::vector<type>& rGlobalValues) const {    \
    CO_SIM_IO_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rLocalValues.size(), rGlobalValues.size(), "SumAll");     \
    rGlobalValues = SumAll(rLocalValues);                                                               \
}                                                                                                       \
virtual type MinAll(const type rLocalValue) const { return rLocalValue; }                               \
virtual std::vector<type> MinAll(const std::vector<type>& rLocalValues) const {                         \
    return rLocalValues;                                                                                \
}                                                                                                       \
virtual void MinAll(const std::vector<type>& rLocalValues, std::vector<type>& rGlobalValues) const {    \
    CO_SIM_IO_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rLocalValues.size(), rGlobalValues.size(), "MinAll");     \
    rGlobalValues = MinAll(rLocalValues);                                                               \
}                                                                                                       \
virtual type MaxAll(const type rLocalValue) const { return rLocalValue; }                               \
virtual std::vector<type> MaxAll(const std::vector<type>& rLocalValues) const {                         \
    return rLocalValues;                                                                                \
}                                                                                                       \
virtual void MaxAll(const std::vector<type>& rLocalValues, std::vector<type>& rGlobalValues) const {    \
    CO_SIM_IO_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rLocalValues.size(), rGlobalValues.size(), "MaxAll");     \
    rGlobalValues = MaxAll(rLocalValues);                                                               \
}                                                                                                       \

#endif

// Compute the partial sum of the given quantity from rank 0 to the current rank (included).
/* This is a wrapper to MPI_Scan.
 * Variants for each method are provided, either returning the reduced value or filling a provided vector buffer.
 * The returned value is defined on all ranks.
 */
#ifndef CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_SCANSUM_INTERFACE_FOR_TYPE
#define CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_SCANSUM_INTERFACE_FOR_TYPE(type)                      \
virtual type ScanSum(const type rLocalValue) const { return rLocalValue; }                          \
virtual std::vector<type> ScanSum(const std::vector<type>& rLocalValues) const {                    \
    return rLocalValues;                                                                            \
}                                                                                                   \
virtual void ScanSum(const std::vector<type>& rLocalValues, std::vector<type>& rPartialSums) const {\
    CO_SIM_IO_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rLocalValues.size(), rPartialSums.size(), "ScanSum"); \
    rPartialSums = ScanSum(rLocalValues);                                                           \
}                                                                                                   \

#endif

// Exchange data with other ranks. This is a wrapper for MPI_Sendrecv, MPI_Send and MPI_Recv.
/* Versions which outputting the result as a return argument or by filling an output buffer argument are provided.
 * The return version has a performance overhead, since the dimensions of the receiving buffer have to be
 * communicated. If the dimensions of the receiving buffer are known at the destination rank, the output buffer
 * variant should be preferred.
 */
#ifndef CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_SENDRECV_INTERFACE_FOR_TYPE
#define CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_SENDRECV_INTERFACE_FOR_TYPE(type)                                 \
virtual type SendRecvImpl(                                                                                      \
    const type rSendValues, const int SendDestination, const int SendTag,                                       \
    const int RecvSource, const int RecvTag) const {                                                            \
    CO_SIM_IO_ERROR_IF( (Rank() != SendDestination) || (Rank() != RecvSource))                                     \
    << "Communication between different ranks is not possible with a serial DataCommunicator." << std::endl;    \
    return rSendValues;                                                                                         \
}                                                                                                               \
virtual std::vector<type> SendRecvImpl(                                                                         \
    const std::vector<type>& rSendValues, const int SendDestination, const int SendTag,                         \
    const int RecvSource, const int RecvTag) const {                                                            \
    CO_SIM_IO_ERROR_IF( (Rank() != SendDestination) || (Rank() != RecvSource))                                     \
    << "Communication between different ranks is not possible with a serial DataCommunicator." << std::endl;    \
    return rSendValues;                                                                                         \
}                                                                                                               \
virtual void SendRecvImpl(                                                                                      \
    const type rSendValues, const int SendDestination, const int SendTag,                                       \
    type& rRecvValues, const int RecvSource, const int RecvTag) const {                                         \
    rRecvValues = SendRecvImpl(rSendValues, SendDestination, SendTag, RecvSource, RecvTag);                     \
}                                                                                                               \
virtual void SendRecvImpl(                                                                                      \
    const std::vector<type>& rSendValues, const int SendDestination, const int SendTag,                         \
    std::vector<type>& rRecvValues, const int RecvSource, const int RecvTag) const {                            \
    CO_SIM_IO_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rSendValues.size(), rRecvValues.size(), "SendRecv");              \
    rRecvValues = SendRecvImpl(rSendValues, SendDestination, SendTag, RecvSource, RecvTag);                     \
}                                                                                                               \
virtual void SendImpl(                                                                                          \
    const std::vector<type>& rSendValues, const int SendDestination, const int SendTag = 0) const {             \
    CO_SIM_IO_ERROR_IF(Rank() != SendDestination)                                                                  \
    << "Communication between different ranks is not possible with a serial DataCommunicator." << std::endl;    \
}                                                                                                               \
virtual void RecvImpl(std::vector<type>& rRecvValues, const int RecvSource, const int RecvTag = 0) const {      \
    CO_SIM_IO_ERROR << "Calling serial DataCommunicator::Recv, which has no meaningful return." << std::endl;      \
}                                                                                                               \

#endif

// Synchronize a buffer to the value held by the broadcasting rank.
/* This is a wrapper for MPI_Bcast.
 *  @param[in/out] The broadcast value (input on SourceRank, output on all other ranks).
 *  @param[in] SourceRank The rank transmitting the value.
 */
#ifndef CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_BROADCAST_INTERFACE_FOR_TYPE
#define CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_BROADCAST_INTERFACE_FOR_TYPE(type)        \
virtual void BroadcastImpl(type& rBuffer, const int SourceRank) const {}                \
virtual void BroadcastImpl(std::vector<type>& rBuffer, const int SourceRank) const {}   \

#endif

/// Wrappers for MPI_Scatter and MPI_Scatterv calls.
/* Versions which outputting the result as a return argument or by filling an output buffer argument are provided.
 * The return version has a performance overhead, since the dimensions of the receiving buffer have to be
 * communicated. If the dimensions of the receiving buffers are known at the destination rank, the output buffer
 * variant should be preferred.
 */
#ifndef CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_SCATTER_INTERFACE_FOR_TYPE
#define CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_SCATTER_INTERFACE_FOR_TYPE(type)                                                              \
virtual std::vector<type> Scatter(const std::vector<type>& rSendValues, const int SourceRank) const {                                       \
     CO_SIM_IO_ERROR_IF( Rank() != SourceRank )                                                                                                \
    << "Communication between different ranks is not possible with a serial DataCommunicator." << std::endl;                                \
    return rSendValues;                                                                                                                     \
}                                                                                                                                           \
virtual void Scatter(                                                                                                                       \
    const std::vector<type>& rSendValues, std::vector<type>& rRecvValues, const int SourceRank) const {                                     \
    CO_SIM_IO_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rSendValues.size(),rRecvValues.size(),"Scatter");                                             \
    rRecvValues = Scatter(rSendValues, SourceRank);                                                                                         \
}                                                                                                                                           \
virtual std::vector<type> Scatterv(const std::vector<std::vector<type>>& rSendValues, const int SourceRank) const {                         \
    CO_SIM_IO_ERROR_IF( Rank() != SourceRank )                                                                                                 \
    << "Communication between different ranks is not possible with a serial DataCommunicator." << std::endl;                                \
    CO_SIM_IO_ERROR_IF( static_cast<unsigned int>(Size()) != rSendValues.size() )                                                              \
    << "Unexpected number of sends in DataCommuncatior::Scatterv (serial DataCommunicator always assumes a single process)." << std::endl;  \
    return rSendValues[0];                                                                                                                  \
}                                                                                                                                           \
virtual void Scatterv(                                                                                                                      \
    const std::vector<type>& rSendValues, const std::vector<int>& rSendCounts, const std::vector<int>& rSendOffsets,                        \
    std::vector<type>& rRecvValues, const int SourceRank) const {                                                                           \
    CO_SIM_IO_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rRecvValues.size(), rSendValues.size(), "Scatterv (values check)");                           \
    CO_SIM_IO_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rSendCounts.size(), 1, "Scatterv (counts check)");                                            \
    CO_SIM_IO_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rSendOffsets.size(), 1, "Scatterv (offsets check)");                                          \
    CO_SIM_IO_ERROR_IF( Rank() != SourceRank )                                                                                                 \
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
#ifndef CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_GATHER_INTERFACE_FOR_TYPE
#define CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_GATHER_INTERFACE_FOR_TYPE(type)                                       \
virtual std::vector<type> Gather(const std::vector<type>& rSendValues, const int DestinationRank) const {           \
    CO_SIM_IO_ERROR_IF( Rank() != DestinationRank )                                                                    \
    << "Communication between different ranks is not possible with a serial DataCommunicator." << std::endl;        \
    return rSendValues;                                                                                             \
}                                                                                                                   \
virtual void Gather(                                                                                                \
    const std::vector<type>& rSendValues, std::vector<type>& rRecvValues, const int DestinationRank) const {        \
    CO_SIM_IO_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rSendValues.size(),rRecvValues.size(),"Gather");                      \
    rRecvValues = Gather(rSendValues, DestinationRank);                                                             \
}                                                                                                                   \
virtual std::vector<std::vector<type>> Gatherv(                                                                     \
    const std::vector<type>& rSendValues, const int DestinationRank) const {                                        \
    CO_SIM_IO_ERROR_IF( Rank() != DestinationRank )                                                                    \
    << "Communication between different ranks is not possible with a serial DataCommunicator." << std::endl;        \
    return std::vector<std::vector<type>>{rSendValues};                                                             \
}                                                                                                                   \
virtual void Gatherv(                                                                                               \
    const std::vector<type>& rSendValues, std::vector<type>& rRecvValues,                                           \
    const std::vector<int>& rRecvCounts, const std::vector<int>& rRecvOffsets, const int DestinationRank) const {   \
    CO_SIM_IO_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rRecvValues.size(), rSendValues.size(), "Gatherv (values check)");    \
    CO_SIM_IO_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rRecvCounts.size(), 1, "Gatherv (counts check)");                     \
    CO_SIM_IO_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rRecvOffsets.size(), 1, "Gatherv (offset check)");                    \
    CO_SIM_IO_ERROR_IF( Rank() != DestinationRank )                                                                    \
    << "Communication between different ranks is not possible with a serial DataCommunicator." << std::endl;        \
    rRecvValues = rSendValues;                                                                                      \
}                                                                                                                   \
virtual std::vector<type> AllGather(const std::vector<type>& rSendValues) const { return rSendValues; }             \
virtual void AllGather(const std::vector<type>& rSendValues, std::vector<type>& rRecvValues) const {                \
    CO_SIM_IO_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rSendValues.size(),rRecvValues.size(),"AllGather");                   \
    rRecvValues = AllGather(rSendValues);                                                                           \
}                                                                                                                   \

#endif

#ifndef CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_PUBLIC_INTERFACE_FOR_TYPE
#define CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_PUBLIC_INTERFACE_FOR_TYPE(type)   \
CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_REDUCE_INTERFACE_FOR_TYPE(type)    \
CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_ALLREDUCE_INTERFACE_FOR_TYPE(type) \
CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_SCANSUM_INTERFACE_FOR_TYPE(type)   \
CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_SCATTER_INTERFACE_FOR_TYPE(type)   \
CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_GATHER_INTERFACE_FOR_TYPE(type)    \

#endif

#ifndef CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_IMPLEMENTATION_FOR_TYPE
#define CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_IMPLEMENTATION_FOR_TYPE(type)   \
CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_SENDRECV_INTERFACE_FOR_TYPE(type)  \
CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_BROADCAST_INTERFACE_FOR_TYPE(type) \

#endif


namespace CoSimIO {
namespace Internals {

///@addtogroup CoSimIO Core
///@{

///@name CoSimIO Classes
///@{

/// Serial (do-nothing) version of a wrapper class for MPI communication.
/** @see MPIDataCommunicator for a working distributed memory implementation.
  */
class DataCommunicator
{
  private:

    template<typename T> class serialization_is_required {
    private:

        template<typename U> struct serialization_traits {
            constexpr static bool is_std_vector = false;
            constexpr static bool value_type_is_compound = false;
        };

        template<typename U> struct serialization_traits<std::vector<U>> {
            constexpr static bool is_std_vector = true;
            constexpr static bool value_type_is_compound = std::is_compound<U>::value;
        };

        constexpr static bool is_vector_of_simple_types = serialization_traits<T>::is_std_vector && !serialization_traits<T>::value_type_is_compound;

    public:
        constexpr static bool value = std::is_compound<T>::value && !is_vector_of_simple_types;
    };

    template<bool value> struct TypeFromBool {};

    template<typename T> void CheckSerializationForSimpleType(const T& rSerializedType, TypeFromBool<true>) const {}
    template<typename T>
    CO_SIM_IO_DEPRECATED//("Calling serialization-based communication for a simple type. Please implement direct communication support for this type.")
    void CheckSerializationForSimpleType(const T& rSerializedType, TypeFromBool<false>) const {}

  public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DataCommunicator() = default;

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
    static std::unique_ptr<DataCommunicator> Create()
    {
        return CoSimIO::make_unique<DataCommunicator>();
    }

    /// Pause program exectution until all threads reach this call.
    /** Wrapper for MPI_Barrier. */
    virtual void Barrier() const {}

    // Complete interface for basic types

    CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_PUBLIC_INTERFACE_FOR_TYPE(int)
    CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_PUBLIC_INTERFACE_FOR_TYPE(unsigned int)
    CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_PUBLIC_INTERFACE_FOR_TYPE(long unsigned int)
    CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_PUBLIC_INTERFACE_FOR_TYPE(double)

    // Reduce operations

    virtual bool AndReduce(
        const bool Value,
        const int Root) const
    {
        return Value;
    }

    virtual bool OrReduce(
        const bool Value,
        const int Root) const
    {
        return Value;
    }

    // Allreduce operations

    virtual bool AndReduceAll(const bool Value) const
    {
        return Value;
    }

    virtual bool OrReduceAll(const bool Value) const
    {
        return Value;
    }

    // Broadcast operations

    /// Synchronize a buffer to the value held by the broadcasting rank.
    /** This is a wrapper for MPI_Bcast.
     *  @param[in/out] The broadcast value (input on SourceRank, output on all other ranks).
     *  @param[in] SourceRank The rank transmitting the value.
     *  This function will transfer basic types (or std::vectors of basic types) directly.
     *  For complex classes, serialization will be used to package the object(s) before communication.
     */
    template<typename TObject>
    void Broadcast(TObject& rBroadcastObject, const int SourceRank) const
    {
        this->BroadcastImpl(rBroadcastObject, SourceRank);
    }

    // Sendrecv operations

    /// Exchange data with other ranks.
    /** This is a wrapper for MPI_Sendrecv.
     *  @param[in] rSendValues Values to send to rank SendDestination.
     *  @param[in] SendDestination Rank the values will be sent to.
     *  @param[in] SendTag Message tag for sent values.
     *  @param[in] RecvSource Rank values are expected from.
     *  @param[in] RecvTag Message tag for received values.
     *  @return Received values from rank RecvSource.
     *  This function will transfer basic types (or std::vectors of basic types) directly.
     *  For complex classes, serialization will be used to package the object(s) before communication.
     */
    template<typename TObject>
    TObject SendRecv(
        const TObject& rSendObject, const int SendDestination, const int SendTag,
        const int RecvSource, const int RecvTag) const
    {
        return this->SendRecvImpl(rSendObject, SendDestination, SendTag, RecvSource, RecvTag);
    }

    /// Exchange data with other ranks.
    /** This is a wrapper for MPI_Sendrecv.
     *  @param[in] rSendValues Values to send to rank SendDestination.
     *  @param[in] RecvSource Rank values are expected from.
     *  @param[in] SendDestination Rank the values will be sent to.
     *  @return Received values from rank RecvSource.
     *  This function will transfer basic types (or std::vectors of basic types) directly.
     *  For complex classes, serialization will be used to package the object(s) before communication.
     */
    template<class TObject>
    TObject SendRecv(
        const TObject& rSendObject, const int SendDestination, const int RecvSource) const
    {
        return this->SendRecvImpl(rSendObject, SendDestination, 0, RecvSource, 0);
    }

    /// Exchange data with other ranks.
    /** This is a wrapper for MPI_Sendrecv.
     *  @param[in] rSendValues Values to send to rank SendDestination.
     *  @param[in] SendDestination Rank the values will be sent to.
     *  @param[in] SendTag Message tag for sent values.
     *  @param[out] rRecvValues Received values from rank RecvSource.
     *  @param[in] RecvSource Rank values are expected from.
     *  @param[in] RecvTag Message tag for received values.
     *  This function will transfer basic types (or std::vectors of basic types) directly.
     *  For complex classes, serialization will be used to package the object(s) before communication.
     */
    template<class TObject>
    void SendRecv(
        const TObject& rSendObject, const int SendDestination, const int SendTag,
        TObject& rRecvObject, const int RecvSource, const int RecvTag) const
    {
        this->SendRecvImpl(rSendObject, SendDestination, SendTag, rRecvObject, RecvSource, RecvTag);
    }

    /// Exchange data with other ranks.
    /** This is a wrapper for MPI_Sendrecv.
     *  @param[in] rSendValues Values to send to rank SendDestination.
     *  @param[in] SendDestination Rank the values will be sent to.
     *  @param[out] rRecvValues Received values from rank RecvSource.
     *  @param[in] RecvSource Rank values are expected from.
     *  This function will transfer basic types (or std::vectors of basic types) directly.
     *  For complex classes, serialization will be used to package the object(s) before communication.
     */
    template<class TObject>
    void SendRecv(
        const TObject& rSendObject, const int SendDestination, TObject& rRecvObject, const int RecvSource) const
    {
        this->SendRecvImpl(rSendObject, SendDestination, 0, rRecvObject, RecvSource, 0);
    }

    /// Exchange data with other ranks.
    /** This is a wrapper for MPI_Send.
     *  @param[in] rSendValues Objects to send to rank SendDestination.
     *  @param[in] SendDestination Rank the data will be sent to.
     *  @param[in] SendTag Message tag for sent values.
     *  This function will transfer basic types (or std::vectors of basic types) directly.
     *  For complex classes, serialization will be used to package the object(s) before communication.
     */
    template<typename TObject>
    void Send(const TObject& rSendValues, const int SendDestination, const int SendTag = 0) const
    {
        this->SendImpl(rSendValues, SendDestination, SendTag);
    }

    /// Exchange data with other ranks.
    /** This is a wrapper for MPI_Recv.
     *  @param[out] rRecvObject Objects to receive from rank RecvSource.
     *  @param[in] RecvSource Rank the data will be received from.
     *  @param[in] RecvTag Message tag for received values.
     *  This function will transfer basic types (or std::vectors of basic types) directly.
     *  For complex classes, serialization will be used to package the object(s) before communication.
     */
    template<typename TObject>
    void Recv(TObject& rRecvObject, const int RecvSource, const int RecvTag = 0) const
    {
        this->RecvImpl(rRecvObject, RecvSource, RecvTag);
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

    /// Check whether this DataCommunicator involves the current rank.
    /** In MPI, if the rank is not involved in communication, the communicator
     *  is MPI_COMM_NULL and is not a valid argument for most MPI calls.
     */
    virtual bool IsDefinedOnThisRank() const
    {
        return true;
    }

    /// Check whether this DataCommunicator is MPI_COMM_NULL.
    /** In MPI, if the rank is not involved in communication, the communicator
     *  is MPI_COMM_NULL and is not a valid argument for most MPI calls.
     */
    virtual bool IsNullOnThisRank() const
    {
        return false;
    }

    ///@}
    ///@name Helper functions for error checking in MPI
    ///@{

    /// This function throws an error on ranks != Sourcerank if Condition evaluates to true.
    /** This method is intended as a helper function to force ranks to stop after an error
     *  is detected on a given rank. A typical use case would be to completely stop the simulation
     *  if an error is detected on the root process.
     *  The intended usage is something like:
     *
     *  CO_SIM_IO_ERROR_IF( data_communicator_instance.BroadcastErrorIfTrue(Condition, Root) )
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
     *  CO_SIM_IO_ERROR_IF_NOT( data_communicator_instance.BroadcastErrorIfFalse(Condition, Root) )
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
     *  CO_SIM_IO_ERROR_IF( data_communicator_instance.ErrorIfTrueOnAnyRank(Condition) )
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
     *  CO_SIM_IO_ERROR_IF_NOT( data_communicator_instance.ErrorIfFalseOnAnyRank(Condition) )
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
        << "Serial do-nothing version of the CoSimIO wrapper for MPI communication.\n"
        << "Rank 0 of 1 assumed." << std::endl;
    }

    ///@}

  protected:

    ///@name Protected operations
    ///@{

    CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_IMPLEMENTATION_FOR_TYPE(int)
    CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_IMPLEMENTATION_FOR_TYPE(unsigned int)
    CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_IMPLEMENTATION_FOR_TYPE(long unsigned int)
    CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_IMPLEMENTATION_FOR_TYPE(double)

    /// Synchronize a buffer to the value held by the broadcasting rank (string version).
    /** This is a wrapper for MPI_Bcast.
     *  @param[in/out] The broadcast value (input on SourceRank, output on all other ranks).
     *  @param[in] SourceRank The rank transmitting the value.
     */
    virtual void BroadcastImpl(std::string& rBuffer, const int SourceRank) const {};

    /// Synchronize a buffer to the value held by the broadcasting rank (generic version).
    /** This is a wrapper for MPI_Bcast, using serialization to package complex objects.
     *  @param[in/out] The broadcast value (input on SourceRank, output on all other ranks).
     *  @param[in] SourceRank The rank transmitting the value.
     */
    template<class TObject>
    void BroadcastImpl(TObject& rBroadcastObject, const int SourceRank) const
    {
        CheckSerializationForSimpleType(rBroadcastObject, TypeFromBool<serialization_is_required<TObject>::value>());
        if (this->IsDistributed())
        {
            unsigned int message_size;
            std::string broadcast_message;
            int rank = this->Rank();
            if (rank == SourceRank)
            {
                StreamSerializer send_serializer;
                send_serializer.save("data", rBroadcastObject);
                broadcast_message = send_serializer.GetStringRepresentation();

                message_size = static_cast<unsigned int>(broadcast_message.size());
            }

            this->Broadcast(message_size, SourceRank);

            if (rank != SourceRank)
            {
                broadcast_message.resize(message_size);
            }

            this->Broadcast(broadcast_message, SourceRank);

            if (rank != SourceRank)
            {
                StreamSerializer recv_serializer(broadcast_message);
                recv_serializer.load("data", rBroadcastObject);
            }
        }
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
    virtual void SendRecvImpl(
        const std::string& rSendValues, const int SendDestination, const int SendTag,
        std::string& rRecvValues, const int RecvSource, const int RecvTag) const
    {
        CO_SIM_IO_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK(rSendValues.size(), rRecvValues.size(), "SendRecv");
        rRecvValues = SendRecvImpl(rSendValues, SendDestination, SendTag, RecvSource, RecvTag);
    }

    /// Exchange data with other ranks (string version).
    /** This is a wrapper for MPI_Sendrecv.
     *  @param[in] rSendValues String to send to rank SendDestination.
     *  @param[in] SendDestination Rank the string will be sent to.
     *  @param[in] SendTag Message tag for sent values.
     *  @param[in] RecvSource Rank the string is expected from.
     *  @param[in] RecvTag Message tag for received values.
     *  @return Received string from rank RecvSource.
     */
    virtual std::string SendRecvImpl(
        const std::string& rSendValues, const int SendDestination, const int SendTag,
        const int RecvSource, const int RecvTag) const
    {
        CO_SIM_IO_ERROR_IF( (Rank() != SendDestination) || (Rank() != RecvSource))
        << "Communication between different ranks is not possible with a serial DataCommunicator." << std::endl;
        return rSendValues;
    }

    /// Exchange data with other ranks (generic version).
    /** This is a wrapper for MPI_Sendrecv that uses serialization to tranfer arbitrary objects.
     *  The objects are expected to be serializable and come in an stl-like container supporting size() and resize()
     *  @param[in] rSendValues Objects to send to rank SendDestination.
     *  @param[in] SendDestination Rank the data will be sent to.
     *  @param[in] SendTag Message tag for sent values.
     *  @param[in] RecvSource Rank the data is expected from.
     *  @param[in] RecvTag Message tag for received values.
     *  @return Received data from rank RecvSource.
     */
    template<class TObject> TObject SendRecvImpl(
        const TObject& rSendObject,
        const int SendDestination, const int SendTag,
        const int RecvSource, const int RecvTag) const
    {
        CheckSerializationForSimpleType(rSendObject, TypeFromBool<serialization_is_required<TObject>::value>());
        if (this->IsDistributed())
        {
            StreamSerializer send_serializer;
            send_serializer.save("data", rSendObject);
            std::string send_message = send_serializer.GetStringRepresentation();

            std::string recv_message = this->SendRecv(send_message, SendDestination, RecvSource);

            StreamSerializer recv_serializer(recv_message);
            TObject recv_object;
            recv_serializer.load("data", recv_object);
            return recv_object;
        }
        else
        {
            CO_SIM_IO_ERROR_IF( (Rank() != SendDestination) || (Rank() != RecvSource))
            << "Communication between different ranks is not possible with a serial DataCommunicator." << std::endl;

            return rSendObject;
        }
    }

    /// Send data to other ranks (string version).
    /** This is a wrapper for MPI_Send.
     *  @param[in] rSendValues String to send to rank SendDestination.
     *  @param[in] SendDestination Rank the string will be sent to.
     *  @param[in] SendTag Message tag for sent values.
     */
    virtual void SendImpl(const std::string& rSendValues, const int SendDestination, const int SendTag) const
    {
        CO_SIM_IO_ERROR_IF(Rank() != SendDestination)
        << "Communication between different ranks is not possible with a serial DataCommunicator." << std::endl;
    }

    /// Exchange data with other ranks (generic version).
    /** This is a wrapper for MPI_Send that uses serialization to tranfer arbitrary objects.
     *  The objects are expected to be serializable and come in an stl-like container supporting size() and resize()
     *  @param[in] rSendValues Objects to send to rank SendDestination.
     *  @param[in] SendDestination Rank the data will be sent to.
     *  @param[in] SendTag Message tag for sent values.
     */
    template<class TObject> void SendImpl(
        const TObject& rSendObject, const int SendDestination, const int SendTag) const
    {
        CheckSerializationForSimpleType(rSendObject, TypeFromBool<serialization_is_required<TObject>::value>());
        if (this->IsDistributed())
        {
            StreamSerializer send_serializer;
            send_serializer.save("data", rSendObject);
            std::string send_message = send_serializer.GetStringRepresentation();

            this->SendImpl(send_message, SendDestination, SendTag);
        }
        else
        {
            CO_SIM_IO_ERROR_IF(Rank() != SendDestination)
            << "Communication between different ranks is not possible with a serial DataCommunicator." << std::endl;
        }
    }

    /// Receive data from other ranks (string version).
    /** This is a wrapper for MPI_Recv.
     *  @param[out] rRecvValues Received string from rank RecvSource.
     *  @param[in] RecvSource Rank the string is expected from.
     *  @param[in] RecvTag Message tag for received values.
     */
    virtual void RecvImpl(std::string& rRecvValues, const int RecvSource, const int RecvTag = 0) const
    {
        CO_SIM_IO_ERROR << "Calling serial DataCommunicator::Recv, which has no meaningful return." << std::endl;
    }

    /// Exchange data with other ranks (generic version).
    /** This is a wrapper for MPI_Recv that uses serialization to tranfer arbitrary objects.
     *  The objects are expected to be serializable and come in an stl-like container supporting size() and resize()
     *  @param[out] rRecvObject Objects to receive from rank RecvSource.
     *  @param[in] RecvSource Rank the data will be received from.
     *  @param[in] RecvTag Message tag for received values.
     */
    template<class TObject> void RecvImpl(
        TObject& rRecvObject, const int RecvSource, const int RecvTag = 0) const
    {
        CheckSerializationForSimpleType(rRecvObject, TypeFromBool<serialization_is_required<TObject>::value>());
        if (this->IsDistributed())
        {
            std::string recv_message;

            this->Recv(recv_message, RecvSource, RecvTag);

            StreamSerializer recv_serializer(recv_message);
            recv_serializer.load("data", rRecvObject);
        }
        else
        {
            CO_SIM_IO_ERROR_IF(Rank() != RecvSource)
            << "Communication between different ranks is not possible with a serial DataCommunicator." << std::endl;
        }
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

} // namespace Internals
} // namespace CoSimIO

#undef CO_SIM_IO_DATA_COMMUNICATOR_DEBUG_SIZE_CHECK

#undef CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_REDUCE_INTERFACE_FOR_TYPE
#undef CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_ALLREDUCE_INTERFACE_FOR_TYPE
#undef CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_SCANSUM_INTERFACE_FOR_TYPE
#undef CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_SENDRECV_INTERFACE_FOR_TYPE
#undef CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_BROADCAST_INTERFACE_FOR_TYPE
#undef CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_SCATTER_INTERFACE_FOR_TYPE
#undef CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_GATHER_INTERFACE_FOR_TYPE
#undef CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_PUBLIC_INTERFACE_FOR_TYPE
#undef CO_SIM_IO_BASE_DATA_COMMUNICATOR_DECLARE_IMPLEMENTATION_FOR_TYPE

#endif // CO_SIM_IO_DATA_COMMUNICATOR_INCLUDED
